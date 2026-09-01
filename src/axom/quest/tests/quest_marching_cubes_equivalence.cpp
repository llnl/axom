// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file quest_marching_cubes_equivalence.cpp
 *
 * @brief Compares the legacy and bump MarchingCubes backends on structured meshes.
 *
 * The outputs do not match one for one.
 *   - The legacy backend emits unwelded facets. Each facet has its own DIM nodes,
       and adjacent facets duplicate shared vertices.
 *   - The bump backend welds the coordset. In 3D it can produce polygons,
 *     which the adaptor fan-triangulates when filling the legacy triangle output.
 *
 * Counts can differ, but underlying geometry of the extracted mesh cannot.
 *
 * Checks:
 *
 * E1. CROSSING-CELL SET.  The set of parent cell ids that produce at least one
 *     facet must match. This is a relatively cheap test.
 *     Supported structured layouts use the same i-fastest cell numbering in both backends,
 *     and permuted field layouts are rejected.
 *
 * E2. VERTEX SET.  An edge of a cell receives a contour vertex if and only if
 *     its two endpoints lie on opposite sides of the isovalue.
 *     That is a property of the sign pattern not of the case table.
 *
 *     The comparison is a two-sided Hausdorff check with a spatial hash.
 *     Any binning scheme has boundary cases where two near-identical points land in different bins.
 *     The tolerance covers that and bump's single-precision edge interpolation.
 *
 * E3. TOTAL AREA (3D) / LENGTH (2D).  On an ambiguous cell the two case tables may
 *     triangulate the same vertex set differently, and triangulating a non-planar bump polygon
 *     gives an area that depends on the fan origin.
 */

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "quest_marching_cubes_equivalence.cpp requires conduit"
#endif
#ifndef AXOM_USE_BUMP
  #error "quest_marching_cubes_equivalence.cpp requires bump"
#endif

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/quest/MarchingCubes.hpp"

#include "conduit_blueprint.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <cstdint>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
using RuntimePolicy = axom::runtime_policy::Policy;

int hostAllocatorID() { return axom::execution_space<axom::SEQ_EXEC>::allocatorID(); }

void copyBlueprintToPolicy(conduit::Node& dst,
                           const conduit::Node& src,
                           RuntimePolicy policy,
                           int allocatorID)
{
  namespace bputils = axom::bump::utilities;

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(policy == RuntimePolicy::cuda)
  {
    bputils::copy<axom::CUDA_EXEC<256>>(dst, src, allocatorID);
    return;
  }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(policy == RuntimePolicy::hip)
  {
    bputils::copy<axom::HIP_EXEC<256>>(dst, src, allocatorID);
    return;
  }
#endif

  AXOM_UNUSED_VAR(policy);
  AXOM_UNUSED_VAR(allocatorID);
  dst.set(src);
}

void copyBlueprintToHost(conduit::Node& dst, const conduit::Node& src)
{
  axom::bump::utilities::copy<axom::SEQ_EXEC>(dst, src, hostAllocatorID());
}

//---------------------------------------------------------------------------
// Analytic fields
//---------------------------------------------------------------------------

//! @brief Signed distance to a plane. On an axis-aligned grid it produces no
//!   ambiguous cells, so it is the case where E3 is safe to assert.
struct PlanarField
{
  double nx, ny, nz, d;
  double operator()(double x, double y, double z) const { return nx * x + ny * y + nz * z - d; }
};

//! @brief Signed distance to a sphere (3D) / circle (2D).
struct RoundField
{
  double cx, cy, cz, r;
  double operator()(double x, double y, double z) const
  {
    const double dx = x - cx, dy = y - cy, dz = z - cz;
    return std::sqrt(dx * dx + dy * dy + dz * dz) - r;
  }
};

/*!
 * @brief A gyroid.
 *
 * Its curvature produces non-planar cut polygons. That makes it the harshest
 * test of the fan-triangulation tolerance in E3.
 */
struct GyroidField
{
  double scale;
  double operator()(double x, double y, double z) const
  {
    const double sx = scale * x, sy = scale * y, sz = scale * z;
    return std::sin(sx) * std::cos(sy) + std::sin(sy) * std::cos(sz) + std::sin(sz) * std::cos(sx);
  }
};

//---------------------------------------------------------------------------
// Mesh construction. Structured plus explicit meshes that both backends accept.
//---------------------------------------------------------------------------

//! @brief Build a single-domain structured explicit mesh on [0,1]^DIM with
//!   @a n cells per side, sampling @a f at nodes in i-fastest order.
template <int DIM, typename Field>
void buildStructured(conduit::Node& mesh, int n, const Field& f, const std::string& fieldName)
{
  static_assert(DIM == 2 || DIM == 3, "DIM must be 2 or 3");
  const int nn = n + 1;
  conduit::index_t N = static_cast<conduit::index_t>(nn) * nn;
  if(DIM == 3)
  {
    N *= nn;
  }

  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "explicit";
  cs["values/x"].set(conduit::DataType::float64(N));
  cs["values/y"].set(conduit::DataType::float64(N));
  auto* x = cs["values/x"].as_float64_ptr();
  auto* y = cs["values/y"].as_float64_ptr();
  double* z = nullptr;
  if(DIM == 3)
  {
    cs["values/z"].set(conduit::DataType::float64(N));
    z = cs["values/z"].as_float64_ptr();
  }

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "structured";
  topo["coordset"] = "coords";
  topo["elements/dims/i"] = n;
  topo["elements/dims/j"] = n;
  if(DIM == 3)
  {
    topo["elements/dims/k"] = n;
  }

  conduit::Node& fld = mesh["fields/" + fieldName];
  fld["topology"] = "mesh";
  fld["association"] = "vertex";
  fld["values"].set(conduit::DataType::float64(N));
  auto* fv = fld["values"].as_float64_ptr();

  const int nk = (DIM == 3) ? nn : 1;
  conduit::index_t idx = 0;
  for(int k = 0; k < nk; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i, ++idx)
      {
        const double px = double(i) / n;
        const double py = double(j) / n;
        const double pz = (DIM == 3) ? double(k) / n : 0.0;
        x[idx] = px;
        y[idx] = py;
        if(z != nullptr)
        {
          z[idx] = pz;
        }
        fv[idx] = f(px, py, pz);
      }
    }
  }
}

/*!
 * @brief Build the same box as buildStructured<3>, but as a uniform coordset + uniform topology.
 *
 * @note n must be a power of two.  The explicit builder writes node coordinates
 *   as double(i)/n while a uniform coordset is evaluated as origin + i*spacing.
 *   Those agree bit-for-bit only when 1/n is exactly representable, and the test
 *   below compares vertex sets, so a non-power-of-two n would fail for a reason
 *   that has nothing to do with the code under test.
 */
template <typename Field>
void buildUniform3D(conduit::Node& mesh, int n, const Field& f, const std::string& fieldName)
{
  const int nn = n + 1;
  const conduit::index_t N = static_cast<conduit::index_t>(nn) * nn * nn;
  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "uniform";
  cs["dims/i"] = nn;
  cs["dims/j"] = nn;
  cs["dims/k"] = nn;
  cs["origin/x"] = 0.0;
  cs["origin/y"] = 0.0;
  cs["origin/z"] = 0.0;
  cs["spacing/dx"] = 1.0 / n;
  cs["spacing/dy"] = 1.0 / n;
  cs["spacing/dz"] = 1.0 / n;

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "uniform";
  topo["coordset"] = "coords";

  conduit::Node& fld = mesh["fields/" + fieldName];
  fld["topology"] = "mesh";
  fld["association"] = "vertex";
  fld["values"].set(conduit::DataType::float64(N));
  auto* fv = fld["values"].as_float64_ptr();
  conduit::index_t idx = 0;
  for(int k = 0; k < nn; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i, ++idx)
      {
        fv[idx] = f(double(i) / n, double(j) / n, double(k) / n);
      }
    }
  }
}

//! @brief The same box again, as a rectilinear coordset + rectilinear topology.
template <typename Field>
void buildRectilinear3D(conduit::Node& mesh, int n, const Field& f, const std::string& fieldName)
{
  const int nn = n + 1;
  const conduit::index_t N = static_cast<conduit::index_t>(nn) * nn * nn;
  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "rectilinear";
  for(const char* comp : {"x", "y", "z"})
  {
    cs[std::string("values/") + comp].set(conduit::DataType::float64(nn));
    auto* v = cs[std::string("values/") + comp].as_float64_ptr();
    for(int i = 0; i < nn; ++i)
    {
      v[i] = double(i) / n;
    }
  }

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "rectilinear";
  topo["coordset"] = "coords";

  conduit::Node& fld = mesh["fields/" + fieldName];
  fld["topology"] = "mesh";
  fld["association"] = "vertex";
  fld["values"].set(conduit::DataType::float64(N));
  auto* fv = fld["values"].as_float64_ptr();
  conduit::index_t idx = 0;
  for(int k = 0; k < nn; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i, ++idx)
      {
        fv[idx] = f(double(i) / n, double(j) / n, double(k) / n);
      }
    }
  }
}

/*!
 * @brief Build a strided-structured (ghost-padded) version of the same box.
 *
 * The real zone extent is n^3, but the coordset and field arrays are allocated
 * over a padded (n+2*g)^3 window, with a topology that has elements/dims/{offsets,strides}.
 * This is the layout quest_marching_cubes_example produces with --strided.
 * The test below exercises that layout directly without relying on the example driver.
 */
template <typename Field>
void buildStridedStructured3D(conduit::Node& mesh,
                              int n,
                              int g,
                              const Field& f,
                              const std::string& fieldName)
{
  const int nnReal = n + 1;          // real points per axis
  const int nnPad = nnReal + 2 * g;  // padded points per axis
  const conduit::index_t N = static_cast<conduit::index_t>(nnPad) * nnPad * nnPad;
  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "explicit";
  for(const char* comp : {"x", "y", "z"})
  {
    cs[std::string("values/") + comp].set(conduit::DataType::float64(N));
  }
  auto* x = cs["values/x"].as_float64_ptr();
  auto* y = cs["values/y"].as_float64_ptr();
  auto* z = cs["values/z"].as_float64_ptr();

  conduit::Node& fld = mesh["fields/" + fieldName];
  fld["topology"] = "mesh";
  fld["association"] = "vertex";
  fld["values"].set(conduit::DataType::float64(N));
  auto* fv = fld["values"].as_float64_ptr();
  fld["offsets"].set(std::vector<conduit::int32> {g, g, g});
  fld["strides"].set(std::vector<conduit::int32> {1, nnPad, nnPad * nnPad});

  // Fill the whole padded window; ghosts get values continuing the same field so
  // a ghost leak shows up as extra facets rather than as garbage.
  conduit::index_t idx = 0;
  for(int k = 0; k < nnPad; ++k)
  {
    for(int j = 0; j < nnPad; ++j)
    {
      for(int i = 0; i < nnPad; ++i, ++idx)
      {
        const double px = double(i - g) / n, py = double(j - g) / n, pz = double(k - g) / n;
        x[idx] = px;
        y[idx] = py;
        z[idx] = pz;
        fv[idx] = f(px, py, pz);
      }
    }
  }

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "structured";
  topo["coordset"] = "coords";
  topo["elements/dims/i"] = n;
  topo["elements/dims/j"] = n;
  topo["elements/dims/k"] = n;
  topo["elements/dims/offsets"].set(std::vector<conduit::int32> {g, g, g});
  topo["elements/dims/strides"].set(std::vector<conduit::int32> {1, nnPad, nnPad * nnPad});
}

//---------------------------------------------------------------------------
// Extracted result, normalized so the two backends are comparable
//---------------------------------------------------------------------------

struct BackendResult
{
  std::set<axom::IndexType> crossingCells;               //!< E1
  std::vector<axom::primal::Point<double, 3>> vertices;  //!< E2 (z==0 in 2D)
  double measure {0.0};                                  //!< E3: area in 3D, length in 2D
  axom::IndexType facetCount {0};
  axom::IndexType nodeCount {0};
};

template <int DIM>
BackendResult runBackend(const conduit::Node& mesh,
                         const std::string& fieldName,
                         double contourVal,
                         RuntimePolicy policy,
                         bool useBump,
                         conduit::Node* bumpBlueprint = nullptr)
{
  namespace quest = axom::quest;

  const int allocatorID = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(useBump);

  conduit::Node execMesh;
  copyBlueprintToPolicy(execMesh, mesh, policy, allocatorID);
  mc.setMesh(execMesh, "mesh");
  mc.setFunctionField(fieldName);
  mc.computeIsocontour(contourVal);

  if(useBump && bumpBlueprint != nullptr)
  {
    conduit::Node bumpBlueprintExec;
    mc.populateContourMeshBlueprint(bumpBlueprintExec);
    copyBlueprintToHost(*bumpBlueprint, bumpBlueprintExec);
  }

  BackendResult r;
  r.facetCount = mc.getContourCellCount();
  r.nodeCount = mc.getContourNodeCount();

  const axom::Array<double, 2> coordsHost(mc.getContourNodeCoords(), hostAllocatorID());
  const axom::Array<axom::IndexType, 2> cornersHost(mc.getContourFacetCorners(), hostAllocatorID());
  const axom::Array<axom::IndexType> parentsHost(mc.getContourFacetParents(), hostAllocatorID());
  const auto coords = coordsHost.view();
  const auto corners = cornersHost.view();
  const auto parents = parentsHost.view();

  for(axom::IndexType f = 0; f < r.facetCount; ++f)
  {
    r.crossingCells.insert(parents[f]);
  }

  for(axom::IndexType v = 0; v < r.nodeCount; ++v)
  {
    axom::primal::Point<double, 3> p {};
    p[0] = coords(v, 0);
    p[1] = coords(v, 1);
    p[2] = (DIM == 3) ? coords(v, 2) : 0.0;
    r.vertices.push_back(p);
  }

  for(axom::IndexType f = 0; f < r.facetCount; ++f)
  {
    if(DIM == 3)
    {
      axom::primal::Point<double, 3> a {}, b {}, c {};
      for(int d = 0; d < 3; ++d)
      {
        a[d] = coords(corners(f, 0), d);
        b[d] = coords(corners(f, 1), d);
        c[d] = coords(corners(f, 2), d);
      }
      const auto u = axom::primal::Vector<double, 3>(a, b);
      const auto w = axom::primal::Vector<double, 3>(a, c);
      r.measure += 0.5 * axom::primal::Vector<double, 3>::cross_product(u, w).norm();
    }
    else
    {
      axom::primal::Point<double, 2> a {}, b {};
      for(int d = 0; d < 2; ++d)
      {
        a[d] = coords(corners(f, 0), d);
        b[d] = coords(corners(f, 1), d);
      }
      r.measure += axom::primal::Vector<double, 2>(a, b).norm();
    }
  }

  return r;
}

//---------------------------------------------------------------------------
// E2 support: two-sided Hausdorff check via a spatial hash.
//
// Deliberately NOT quantize-and-compare: binning to a grid has a boundary
// problem where two points 1e-12 apart straddle a bin edge and compare
// unequal.  Hashing into cells of side `tol` and probing the 3^3 neighborhood
// finds any partner within `tol` regardless of where the bin edges fall.
//---------------------------------------------------------------------------

using CellKey = std::int64_t;

CellKey cellKey(std::int64_t i, std::int64_t j, std::int64_t k)
{
  // Small mixing hash; the coordinate range here is [0,1] so the cell indices
  // are bounded by 1/tol and collisions are handled by the bucket vector.
  const std::int64_t h = (i * 73856093) ^ (j * 19349663) ^ (k * 83492791);
  return h;
}

class PointLocator
{
public:
  PointLocator(const std::vector<axom::primal::Point<double, 3>>& pts, double tol)
    : m_pts(pts)
    , m_tol(tol)
  {
    for(std::size_t n = 0; n < pts.size(); ++n)
    {
      m_buckets[keyOf(pts[n])].push_back(n);
    }
  }

  //! @brief Distance from @a q to the nearest stored point, or infinity if none
  //!   lies within the search neighborhood.
  double nearestDistance(const axom::primal::Point<double, 3>& q) const
  {
    const auto ci = cellIndex(q);
    double best = std::numeric_limits<double>::infinity();
    for(std::int64_t di = -1; di <= 1; ++di)
    {
      for(std::int64_t dj = -1; dj <= 1; ++dj)
      {
        for(std::int64_t dk = -1; dk <= 1; ++dk)
        {
          const auto it = m_buckets.find(cellKey(ci[0] + di, ci[1] + dj, ci[2] + dk));
          if(it == m_buckets.end())
          {
            continue;
          }
          for(const auto n : it->second)
          {
            const double d = axom::primal::Vector<double, 3>(q, m_pts[n]).norm();
            best = std::min(best, d);
          }
        }
      }
    }
    return best;
  }

private:
  axom::StackArray<std::int64_t, 3> cellIndex(const axom::primal::Point<double, 3>& p) const
  {
    axom::StackArray<std::int64_t, 3> c;
    for(int d = 0; d < 3; ++d)
    {
      c[d] = static_cast<std::int64_t>(std::floor(p[d] / m_tol));
    }
    return c;
  }

  CellKey keyOf(const axom::primal::Point<double, 3>& p) const
  {
    const auto c = cellIndex(p);
    return cellKey(c[0], c[1], c[2]);
  }

  const std::vector<axom::primal::Point<double, 3>>& m_pts;
  double m_tol;
  std::unordered_map<CellKey, std::vector<std::size_t>> m_buckets;
};

//! @brief One-sided Hausdorff: max over @a from of the distance to the nearest
//!   point of @a to.  Returns the max and the index attaining it.
double oneSidedHausdorff(const std::vector<axom::primal::Point<double, 3>>& from,
                         const PointLocator& to,
                         std::size_t& argMax)
{
  double worst = 0.0;
  argMax = 0;
  for(std::size_t n = 0; n < from.size(); ++n)
  {
    const double d = to.nearestDistance(from[n]);
    if(d > worst)
    {
      worst = d;
      argMax = n;
    }
  }
  return worst;
}

//---------------------------------------------------------------------------
// E3 support 1/2. Ambiguous cells.
//
// Two reasons the case tables can legitimately triangulate the same vertex set
// differently:
//
//   Face ambiguity. A face has the checkerboard sign pattern. Its diagonal
//   pairs match internally and disagree with each other.
//
//   Body-diagonal ambiguity. The minority sign class is exactly a body-diagonal
//   pair. This is classic case 4. It has no ambiguous face, so a face-only
//   detector misses it.
//
// Expressed in corner SIGNS only, so it is independent of either backend's
// table indexing convention.  Corner n is (i,j,k) with i fastest: n = i+2j+4k.
//
// Exhaustively self-tested below (ambiguity_detector_selftest). A detector
// without a negative control is not much use.
//---------------------------------------------------------------------------

//! Cyclic corner order of each of the 6 faces.
constexpr int kFaces[6][4] =
  {{0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 2, 6, 4}, {1, 3, 7, 5}};
//! The 4 body diagonals.
constexpr int kBodyDiagonals[4][2] = {{0, 7}, {1, 6}, {2, 5}, {3, 4}};

bool cellHasFaceAmbiguity3D(const bool s[8])
{
  for(const auto& f : kFaces)
  {
    if(s[f[0]] == s[f[2]] && s[f[1]] == s[f[3]] && s[f[0]] != s[f[1]])
    {
      return true;
    }
  }
  return false;
}

bool cellHasBodyDiagonalAmbiguity3D(const bool s[8])
{
  int nPos = 0;
  for(int i = 0; i < 8; ++i)
  {
    nPos += s[i] ? 1 : 0;
  }
  if(nPos != 2 && nPos != 6)
  {
    return false;
  }
  const bool minority = (nPos == 2);
  int a = -1, b = -1;
  for(int i = 0; i < 8; ++i)
  {
    if(s[i] == minority)
    {
      (a < 0 ? a : b) = i;
    }
  }
  for(const auto& d : kBodyDiagonals)
  {
    if((a == d[0] && b == d[1]) || (a == d[1] && b == d[0]))
    {
      return true;
    }
  }
  return false;
}

bool cellIsAmbiguous3D(const bool s[8])
{
  return cellHasFaceAmbiguity3D(s) || cellHasBodyDiagonalAmbiguity3D(s);
}

//! 2D: corners cyclic 0,1,3,2 (same n = i+2j convention).  Only the two
//! checkerboard patterns are ambiguous.
bool cellIsAmbiguous2D(const bool s[4])
{
  return (s[0] == s[3]) && (s[1] == s[2]) && (s[0] != s[1]);
}

/*!
 * @brief Count ambiguous cells.
 *
 * @note The corner test uses `>=`, matching MarchingCubesImpl::computeCrossingCase
 *   and (after the isoValueForBump nudge) the bump backend.  Using `>` here was a
 *   bug in the first version: it disagreed with the code under test at exactly
 *   the nodes where the tie convention matters, so the detector could classify a
 *   different cell set than either backend actually cut.
 */
template <int DIM, typename Field>
axom::IndexType countAmbiguousCells(int n, const Field& f, double contourVal)
{
  auto sign = [&](int i, int j, int k) {
    const double px = double(i) / n, py = double(j) / n, pz = (DIM == 3) ? double(k) / n : 0.0;
    return f(px, py, pz) >= contourVal;
  };

  axom::IndexType count = 0;
  const int nk = (DIM == 3) ? n : 1;
  for(int k = 0; k < nk; ++k)
  {
    for(int j = 0; j < n; ++j)
    {
      for(int i = 0; i < n; ++i)
      {
        if(DIM == 2)
        {
          const bool s[4] = {sign(i, j, 0),
                             sign(i + 1, j, 0),
                             sign(i, j + 1, 0),
                             sign(i + 1, j + 1, 0)};
          count += cellIsAmbiguous2D(s) ? 1 : 0;
        }
        else
        {
          const bool s[8] = {sign(i, j, k),
                             sign(i + 1, j, k),
                             sign(i, j + 1, k),
                             sign(i + 1, j + 1, k),
                             sign(i, j, k + 1),
                             sign(i + 1, j, k + 1),
                             sign(i, j + 1, k + 1),
                             sign(i + 1, j + 1, k + 1)};
          count += cellIsAmbiguous3D(s) ? 1 : 0;
        }
      }
    }
  }
  return count;
}

//---------------------------------------------------------------------------
// E3 support 2/2. Fan-triangulation sensitivity.
//
// The quest adaptor fan-triangulates bump's polygonal cut faces from local
// corner 0 (MarchingCubesBumpAdaptor.hpp, adaptCutFieldOutputViews).  For a
// planar polygon every fan gives the same area. For a non-planar polygon the
// area depends on which corner the fan starts from, so the reported surface
// area is partly an artifact of an arbitrary choice. On a high-curvature field
// the cut polygons can be markedly non-planar. That, not table ambiguity, is
// what makes a total-area comparison against a differently triangulated backend
// questionable.
//
// Measure it instead of guessing. Re-fan each polygon from corner 1 and report
// the relative spread. Zero spread means every polygon is planar, or already a
// triangle. In that case E3 is a fair comparison.
//---------------------------------------------------------------------------

struct FanSensitivity
{
  double areaFan0 {0.0};
  double areaFan1 {0.0};
  double relSpread {0.0};         //!< Aggregate |A(fan@0) - A(fan@1)| / A(fan@0)
  double maxPolyRelSpread {0.0};  //!< Worst single-polygon spread. The aggregate can cancel.
  axom::IndexType polygonCount {0};
  axom::IndexType maxCorners {0};
};

double triArea(const axom::primal::Point<double, 3>& a,
               const axom::primal::Point<double, 3>& b,
               const axom::primal::Point<double, 3>& c)
{
  const auto u = axom::primal::Vector<double, 3>(a, b);
  const auto w = axom::primal::Vector<double, 3>(a, c);
  return 0.5 * axom::primal::Vector<double, 3>::cross_product(u, w).norm();
}

//! @brief Area of a polygon fan-triangulated starting at local corner @a origin.
double polygonFanArea(const std::vector<axom::primal::Point<double, 3>>& v, int origin)
{
  const int N = static_cast<int>(v.size());
  if(N < 3)
  {
    return 0.0;
  }
  double area = 0.0;
  for(int t = 0; t < N - 2; ++t)
  {
    area += triArea(v[origin], v[(origin + 1 + t) % N], v[(origin + 2 + t) % N]);
  }
  return area;
}

//! @brief Measure how much bump's polygonal output's area depends on the fan origin.
FanSensitivity measureFanSensitivity(const conduit::Node& contourDom)
{
  FanSensitivity fs;

  const conduit::Node& topo = contourDom.fetch_existing("topologies").child(0);
  const conduit::Node& elems = topo.fetch_existing("elements");
  if(!elems.has_child("sizes"))
  {
    return fs;
  }
  const auto sizes = elems.fetch_existing("sizes").as_index_t_accessor();
  const auto offsets = elems.fetch_existing("offsets").as_index_t_accessor();
  const auto conn = elems.fetch_existing("connectivity").as_index_t_accessor();

  const std::string csName = topo.fetch_existing("coordset").as_string();
  const conduit::Node& vals = contourDom.fetch_existing("coordsets/" + csName + "/values");
  const bool has3 = vals.has_child("z");
  const auto xs = vals.fetch_existing("x").as_double_accessor();
  const auto ys = vals.fetch_existing("y").as_double_accessor();

  for(conduit::index_t z = 0; z < sizes.number_of_elements(); ++z)
  {
    const auto nc = sizes[z];
    fs.maxCorners = std::max(fs.maxCorners, static_cast<axom::IndexType>(nc));
    if(nc < 3)
    {
      continue;
    }
    ++fs.polygonCount;
    std::vector<axom::primal::Point<double, 3>> v;
    for(conduit::index_t c = 0; c < nc; ++c)
    {
      const auto id = conn[offsets[z] + c];
      axom::primal::Point<double, 3> p {};
      p[0] = xs[id];
      p[1] = ys[id];
      p[2] = has3 ? vals.fetch_existing("z").as_double_accessor()[id] : 0.0;
      v.push_back(p);
    }
    const double a0 = polygonFanArea(v, 0);
    const double a1 = polygonFanArea(v, 1);
    fs.areaFan0 += a0;
    fs.areaFan1 += a1;
    fs.maxPolyRelSpread = std::max(fs.maxPolyRelSpread, std::abs(a0 - a1) / std::max(a0, 1.0e-300));
  }

  fs.relSpread = std::abs(fs.areaFan0 - fs.areaFan1) / std::max(fs.areaFan0, 1.0e-300);
  return fs;
}

//---------------------------------------------------------------------------
// The comparison itself
//---------------------------------------------------------------------------

template <int DIM, typename Field>
void compareBackends(int n,
                     const Field& f,
                     double contourVal,
                     RuntimePolicy policy,
                     const std::string& label,
                     double vertexTol = 1.0e-5,
                     // bump's FieldIntersector interpolates edge crossings in
                     // float (FieldType == float), so contour vertices carry
                     // ~1e-7 relative error and the measure inherits it.  An
                     // initial 1e-9 here was a TEST bug, not a code defect: it
                     // is below what single-precision interpolation can deliver.
                     // Measured relDiff on the passing cases is 5e-9 to 3e-6.
                     double measureRelTol = 1.0e-5)
{
  const std::string fieldName = "fcn";
  conduit::Node mesh;
  buildStructured<DIM>(mesh, n, f, fieldName);

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  const auto legacy = runBackend<DIM>(mesh, fieldName, contourVal, policy, /*useBump=*/false);
  conduit::Node bumpBp;
  const auto bump = runBackend<DIM>(mesh, fieldName, contourVal, policy, /*useBump=*/true, &bumpBp);

  const auto ambiguous = countAmbiguousCells<DIM>(n, f, contourVal);
  FanSensitivity fan;
  if(DIM == 3 && bumpBp.number_of_children() > 0)
  {
    fan = measureFanSensitivity(bumpBp.child(0));
  }

  SLIC_INFO(axom::fmt::format(
    "[{}] legacy: {} facets / {} nodes / {} cells; bump: {} facets / {} nodes / {} cells; "
    "ambiguous cells: {}",
    label,
    legacy.facetCount,
    legacy.nodeCount,
    legacy.crossingCells.size(),
    bump.facetCount,
    bump.nodeCount,
    bump.crossingCells.size(),
    ambiguous));

  ASSERT_GT(legacy.facetCount, 0) << "[" << label << "] legacy produced an empty contour";
  ASSERT_GT(bump.facetCount, 0) << "[" << label << "] bump produced an empty contour";

  // E1. crossing-cell set
  {
    std::vector<axom::IndexType> onlyLegacy, onlyBump;
    std::set_difference(legacy.crossingCells.begin(),
                        legacy.crossingCells.end(),
                        bump.crossingCells.begin(),
                        bump.crossingCells.end(),
                        std::back_inserter(onlyLegacy));
    std::set_difference(bump.crossingCells.begin(),
                        bump.crossingCells.end(),
                        legacy.crossingCells.begin(),
                        legacy.crossingCells.end(),
                        std::back_inserter(onlyBump));

    EXPECT_TRUE(onlyLegacy.empty())
      << "E1 [" << label << "]: " << onlyLegacy.size()
      << " cells produce facets in legacy but not bump (first: " << onlyLegacy.front()
      << "). A cell dropped by the bump crossing pre-filter is the likely cause.";
    EXPECT_TRUE(onlyBump.empty()) << "E1 [" << label << "]: " << onlyBump.size()
                                  << " cells produce facets in bump but not legacy (first: "
                                  << (onlyBump.empty() ? -1 : onlyBump.front()) << ").";
  }

  // E2. vertex set (table-independent; see file header)
  {
    const PointLocator legacyLoc(legacy.vertices, vertexTol);
    const PointLocator bumpLoc(bump.vertices, vertexTol);

    std::size_t argMax = 0;
    const double bumpToLegacy = oneSidedHausdorff(bump.vertices, legacyLoc, argMax);
    EXPECT_LT(bumpToLegacy, vertexTol)
      << "E2 [" << label << "]: a bump contour vertex has no legacy counterpart within tolerance"
      << " (worst distance " << bumpToLegacy << " at bump vertex " << argMax << ").";

    const double legacyToBump = oneSidedHausdorff(legacy.vertices, bumpLoc, argMax);
    EXPECT_LT(legacyToBump, vertexTol)
      << "E2 [" << label << "]: a legacy contour vertex has no bump counterpart within tolerance"
      << " (worst distance " << legacyToBump << " at legacy vertex " << argMax << ").";
  }

  // E2b. welding happened
  // Legacy stores DIM nodes per facet with no sharing; bump welds.  On a
  // surface with shared edges the welded count must be strictly smaller.
  EXPECT_LT(bump.nodeCount, legacy.nodeCount)
    << "[" << label << "]: bump node count is not smaller than legacy's. Welding regressed.";

  // E3. measure (assert only when unambiguous)
  const double relDiff = std::abs(bump.measure - legacy.measure) / std::max(legacy.measure, 1.0e-300);
  SLIC_INFO(axom::fmt::format("[{}] measure legacy={:.12g} bump={:.12g} relDiff={:.3e}",
                              label,
                              legacy.measure,
                              bump.measure,
                              relDiff));
  SLIC_INFO(axom::fmt::format(
    "[{}] fan sensitivity: {} polygons, max {} corners, area(fan@0)={:.12g} area(fan@1)={:.12g} "
    "relSpread={:.3e} maxPolyRelSpread={:.3e}",
    label,
    fan.polygonCount,
    fan.maxCorners,
    fan.areaFan0,
    fan.areaFan1,
    fan.relSpread,
    fan.maxPolyRelSpread));

  // E3 tolerance.
  //
  // Do not gate the check on fan spread. That creates cliffs where a tiny change
  // in resolution flips the assertion on and off. Instead fold the measured fan
  // spread into the tolerance. bump's reported area is only stable within that
  // spread. Differences below it do not say much. Differences above it do.
  const double e3Tol = std::max(measureRelTol, fan.relSpread);
  EXPECT_LT(fan.relSpread, 0.1)
    << "[" << label << "]: fan-origin spread is so large that E3 is nearly vacuous; bump's "
    << "polygons are extremely non-planar at this resolution.";
  if(ambiguous == 0)
  {
    EXPECT_LT(relDiff, e3Tol)
      << "E3 [" << label << "]: contour measure differs by more than the fan-origin ambiguity ("
      << fan.relSpread << ") can explain, and there are no ambiguous cells, so the two backends "
      << "disagree on the triangulation of an identical vertex set.";
  }
  else
  {
    SLIC_INFO(axom::fmt::format(
      "[{}] E3 not asserted: {} ambiguous cells, where the two case tables may legitimately "
      "triangulate the same vertex set differently and no principled tolerance exists.",
      label,
      ambiguous));
  }
}

/*!
 * @brief An isovalue outside the data range is valid input, not an error.
 *
 * Before the fix, runExtraction() allocated m_output before dispatching and left
 * it non-null-but-EMPTY on the no-crossing path.  hasContourMeshBlueprint() then
 * reported true, populateContourMeshBlueprint passed its guard, and
 * triangulateBlueprintMesh reached fetch_existing("topologies") on an empty node.
 */
void test_empty_contour(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  conduit::Node mesh;
  RoundField f {0.5, 0.5, 0.5, 0.25};
  buildStructured<3>(mesh, 6, f, "fcn");

  const int allocatorID = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(true);
  mc.setMesh(mesh, "mesh");
  mc.setFunctionField("fcn");
  mc.computeIsocontour(1000.0);  // far outside the range of the signed distance

  EXPECT_EQ(mc.getContourCellCount(), 0);
  EXPECT_EQ(mc.getContourNodeCount(), 0);

  // Neither of these may throw or abort.
  conduit::Node bp;
  mc.populateContourMeshBlueprint(bp);
  conduit::Node bpTri;
  mc.populateContourMeshBlueprint(bpTri, /*triangulate=*/true);

  conduit::Node relinquished;
  mc.relinquishContourDataBlueprint(relinquished);
  SUCCEED();
}

/*!
 * @brief Uniform and rectilinear input must work, and must agree with the
 *   explicit-structured mesh describing the same geometry.
 *
 * setDomain() accepts "uniform" and "rectilinear" and the sphinx/RELEASE-NOTES
 * advertise them, but every m_isStructured path constructed MeshViewUtil, which
 * requires a "structured" topology AND an "explicit" coordset and SLIC_ERRORs
 * otherwise. That used to hard-error inside the crossing pre-filter.
 *
 * The legacy backend cannot read uniform or rectilinear input, so there is no
 * direct legacy reference. Instead compare bump-on-uniform and
 * bump-on-rectilinear against bump-on-structured-explicit. The tests above
 * already pin bump-on-structured-explicit to the legacy backend. All three
 * describe the same box, so the results must be identical.
 */
void test_uniform_and_rectilinear(RuntimePolicy policy)
{
  const int n = 8;  // power of two: see buildUniform3D's note on coordinate agreement
  RoundField f {0.5, 0.5, 0.5, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node structured, uniform, rectilinear;
  buildStructured<3>(structured, n, f, fieldName);
  buildUniform3D(uniform, n, f, fieldName);
  buildRectilinear3D(rectilinear, n, f, fieldName);

  for(const auto& m : {&uniform, &rectilinear})
  {
    conduit::Node info;
    ASSERT_TRUE(conduit::blueprint::mesh::verify(*m, info)) << info.to_yaml();
  }

  const auto refRun = runBackend<3>(structured, fieldName, 0.0, policy, /*useBump=*/true);
  ASSERT_GT(refRun.facetCount, 0);

  const auto uniformRun = runBackend<3>(uniform, fieldName, 0.0, policy, /*useBump=*/true);
  const auto rectRun = runBackend<3>(rectilinear, fieldName, 0.0, policy, /*useBump=*/true);

  SLIC_INFO(
    axom::fmt::format("[uniform/rectilinear] structured-explicit: {} facets; uniform: {} facets; "
                      "rectilinear: {} facets",
                      refRun.facetCount,
                      uniformRun.facetCount,
                      rectRun.facetCount));

  for(const auto& kv :
      {std::make_pair("uniform", &uniformRun), std::make_pair("rectilinear", &rectRun)})
  {
    const std::string what = kv.first;
    const BackendResult& r = *kv.second;
    EXPECT_EQ(r.crossingCells, refRun.crossingCells)
      << what << " cut a different cell set than the equivalent structured-explicit mesh";
    EXPECT_EQ(r.facetCount, refRun.facetCount) << what << " produced a different facet count";
    EXPECT_EQ(r.nodeCount, refRun.nodeCount) << what << " produced a different node count";
    EXPECT_NEAR(r.measure, refRun.measure, 1.0e-12 * std::max(refRun.measure, 1.0))
      << what << " produced a different contour measure";
  }
}

/*!
 * @brief A float32 function field must be rejected, not silently misread.
 *
 * The structured pre-filter reads the field via
 * MeshViewUtil::getConstFieldView<double>(), which assumes the values are
 * double and does not check. A float32 field was reinterpreted as float64, so
 * the pre-filter selected a garbage cell set while bump's extractor read the
 * field correctly. That is a wrong answer with no error. The minimum bar here
 * is a loud rejection that names the dtype.
 */
void test_float32_field_rejected(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  const int n = 6;
  RoundField f {0.5, 0.5, 0.5, 0.25};
  conduit::Node mesh;
  buildStructured<3>(mesh, n, f, "fcn");

  // Re-write the function field as float32, keeping the same values.
  {
    const conduit::Node& n_old = mesh.fetch_existing("fields/fcn/values");
    const auto acc = n_old.as_double_accessor();
    const conduit::index_t N = n_old.dtype().number_of_elements();
    std::vector<float> tmp(static_cast<std::size_t>(N));
    for(conduit::index_t i = 0; i < N; ++i)
    {
      tmp[static_cast<std::size_t>(i)] = static_cast<float>(acc[i]);
    }
    mesh["fields/fcn/values"].set(tmp.data(), static_cast<conduit::index_t>(tmp.size()));
  }
  ASSERT_TRUE(mesh.fetch_existing("fields/fcn/values").dtype().is_float32());

  const int allocatorID = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(true);

  // SLIC's default handler aborts, so this is a death test.  SimpleLogger writes
  // to stdout while gtest's death test captures only stderr, so the message must
  // be routed to stderr INSIDE the forked child or the regex has nothing to
  // match (an empty "Actual msg" is the symptom).
  EXPECT_DEATH_IF_SUPPORTED(
    {
      axom::slic::addStreamToAllMsgLevels(
        new axom::slic::GenericOutputStream(&std::cerr, "[<LEVEL>] <MESSAGE>\n"));
      mc.setMesh(mesh, "mesh");
      mc.setFunctionField("fcn");
      mc.computeIsocontour(0.0);
    },
    "float64");
}

//! @brief Run the bump backend in a death-test child and require a field-layout error.
void expectBumpFieldLayoutRejected(const conduit::Node& mesh,
                                   RuntimePolicy policy,
                                   const char* expectedMessage)
{
  EXPECT_DEATH_IF_SUPPORTED(
    {
      axom::slic::addStreamToAllMsgLevels(
        new axom::slic::GenericOutputStream(&std::cerr, "[<LEVEL>] <MESSAGE>\n"));
      const int allocatorID = axom::policyToDefaultAllocatorID(policy);
      axom::quest::MarchingCubes mc(policy,
                                    allocatorID,
                                    axom::quest::MarchingCubesDataParallelism::byPolicy);
      mc.setUseBumpBackend(true);
      mc.setMesh(mesh, "mesh");
      mc.setFunctionField("fcn");
      mc.computeIsocontour(0.0);
    },
    expectedMessage);
}

/*!
 * @brief Unsupported strided function-field layouts must be rejected rather
 *   than silently indexed with the topology's compact zone numbering.
 *
 * A permuted field is not i-fastest.  A separately padded field can remain
 * i-fastest, but its offsets and strides differ from those of the topology.
 * bump's flat field view cannot represent either case correctly.
 */
void test_invalid_field_layouts_rejected(RuntimePolicy policy)
{
  constexpr int n = 6;
  constexpr int pad = 2;
  constexpr int nnPad = n + 1 + 2 * pad;
  RoundField f {0.5, 0.5, 0.5, 0.25};

  conduit::Node permuted;
  buildStridedStructured3D(permuted, n, pad, f, "fcn");
  permuted["fields/fcn/strides"].set(std::vector<conduit::int32> {nnPad * nnPad, nnPad, 1});
  expectBumpFieldLayoutRejected(permuted, policy, "i-fastest");

  conduit::Node mismatched;
  buildStridedStructured3D(mismatched, n, pad, f, "fcn");
  mismatched["fields/fcn/offsets"].set(std::vector<conduit::int32> {pad + 1, pad, pad});
  expectBumpFieldLayoutRejected(mismatched, policy, "same offsets and strides");
}

/*!
 * @brief An input mesh already carrying a field named "originalElements" must
 *   not hijack the reported parent cell ids.
 *
 * bump's TableBasedExtractor::makeOriginalElements branches on whether the INPUT
 * mesh has a field of the configured name and, if so, maps those values forward
 * instead of writing zone indices.  Any mesh produced by a prior bump operation
 * carries exactly that field, and the empty "fields" option does not suppress
 * the branch.  MarchingCubes therefore requests a private name and renames the
 * result back before anyone sees it.
 */
void test_original_elements_collision(RuntimePolicy policy)
{
  const int n = 8;
  RoundField f {0.5, 0.5, 0.5, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node clean, poisoned;
  buildStructured<3>(clean, n, f, fieldName);
  poisoned.set(clean);

  // Plant a decoy: an element field of the colliding name whose values are
  // deliberately nothing like zone indices.
  {
    const conduit::index_t nCells = static_cast<conduit::index_t>(n) * n * n;
    conduit::Node& fld = poisoned["fields/originalElements"];
    fld["topology"] = "mesh";
    fld["association"] = "element";
    fld["values"].set(conduit::DataType::int64(nCells));
    auto* v = fld["values"].as_int64_ptr();
    for(conduit::index_t i = 0; i < nCells; ++i)
    {
      v[i] = -7;  // if these leak through as parent ids, the check below fails
    }
  }
  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(poisoned, info)) << info.to_yaml();

  const auto cleanRun = runBackend<3>(clean, fieldName, 0.0, policy, /*useBump=*/true);
  const auto poisonedRun = runBackend<3>(poisoned, fieldName, 0.0, policy, /*useBump=*/true);

  ASSERT_GT(cleanRun.facetCount, 0);
  EXPECT_EQ(poisonedRun.crossingCells, cleanRun.crossingCells)
    << "a pre-existing 'originalElements' field on the input changed the reported parent ids";
  for(const auto id : poisonedRun.crossingCells)
  {
    EXPECT_GE(id, 0) << "decoy 'originalElements' values leaked through as parent cell ids";
  }
}

/*!
 * @brief Strided-structured input works with the bump backend.
 *
 * This directly covers the native strided path. Earlier example-driver
 * compaction hid it by densifying coordinates and stripping metadata before
 * MarchingCubes saw the mesh.
 *
 * Static reading says it should work. dispatch_only_structured_topology routes
 * elements/dims/{offsets,strides} to make_strided_structured_topology, and
 * StridedStructuredIndexing::indexToLogicalIndex uses the LOCAL zone dims, so
 * bump's zone index space is compact and i-fastest. That matches quest's
 * MDMapping(cellShape, COLUMN). This test settles it empirically.
 *
 * Oracle. Use the same geometry as a dense structured mesh, which the tests
 * above already pin to the legacy backend. Ghost values are poisoned, so an
 * implementation that ignores the offsets cannot accidentally agree.
 */
void test_strided_structured(RuntimePolicy policy)
{
  const int n = 8;
  const int pad = 2;  // ghost layers
  RoundField f {0.5, 0.5, 0.5, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node compact, strided;
  buildStructured<3>(compact, n, f, fieldName);
  buildStridedStructured3D(strided, n, pad, f, fieldName);

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(strided, info)) << info.to_yaml();

  /*
    Three-way agreement on ghost-padded input.

    This faulted before the dispatch-path fix in dispatch_structured_topology.hpp:
    the strided predicate probed topo.has_path("offsets") instead of
    "elements/dims/offsets", so a padded mesh was read as compact and makeTopology
    walked off the end of the coordset.

    Legacy-on-compact is the reference (pinned to bump-on-compact by the tests
    above).  Legacy-on-strided shows the padded fixture is well formed, so a
    bump-only failure cannot be blamed on the fixture.  Bump-on-strided is the
    claim.  The parent-id range check catches a ghost leak, which is the failure
    mode a facet count alone would miss.
  */
  const auto legacyCompact = runBackend<3>(compact, fieldName, 0.0, policy, /*useBump=*/false);
  const auto legacyStrided = runBackend<3>(strided, fieldName, 0.0, policy, /*useBump=*/false);
  const auto bumpStrided = runBackend<3>(strided, fieldName, 0.0, policy, /*useBump=*/true);

  SLIC_INFO(axom::fmt::format(
    "[strided] legacy/compact={} facets, legacy/strided={} facets, bump/strided={} facets",
    legacyCompact.facetCount,
    legacyStrided.facetCount,
    bumpStrided.facetCount));

  ASSERT_GT(legacyCompact.facetCount, 0);
  ASSERT_EQ(legacyStrided.crossingCells, legacyCompact.crossingCells)
    << "the legacy backend disagrees between strided and compact input, so the padded fixture "
    << "is malformed. Fix the fixture before reading anything into the bump result";

  EXPECT_EQ(bumpStrided.crossingCells, legacyCompact.crossingCells)
    << "bump on strided-structured input cut a different cell set than the equivalent compact "
    << "mesh: a ghost leak or a zone-numbering mismatch";

  const axom::IndexType nCells = static_cast<axom::IndexType>(n) * n * n;
  for(const auto id : bumpStrided.crossingCells)
  {
    EXPECT_GE(id, 0) << "parent id below the real zone range (ghost leak)";
    EXPECT_LT(id, nCells) << "parent id above the real zone range (ghost leak)";
  }
}

//---------------------------------------------------------------------------
// Test bodies
//---------------------------------------------------------------------------

void test_planar_3d(RuntimePolicy policy)
{
  // Plane z = 0.5. An axis-aligned planar field has no saddle cells, so E3 is
  // safe to assert here. This is the strictest case.
  PlanarField f {0.0, 0.0, 1.0, 0.5};
  compareBackends<3>(8, f, 0.0, policy, "planar3d");
}

void test_oblique_planar_3d(RuntimePolicy policy)
{
  // Oblique plane: still no saddles, but every case-table entry gets exercised
  // rather than just the axis-aligned ones.
  const double s = 1.0 / std::sqrt(1.0 + 0.16 + 1.44);
  PlanarField f {1.0 * s, 0.4 * s, 1.2 * s, 1.3 * s};
  compareBackends<3>(8, f, 0.0, policy, "oblique_planar3d");
}

void test_round_3d(RuntimePolicy policy)
{
  RoundField f {0.5, 0.5, 0.5, 0.25};
  compareBackends<3>(12, f, 0.0, policy, "round3d");
}

void test_gyroid_3d(RuntimePolicy policy)
{
  // Chosen to produce strongly non-planar polygons. E3 remains asserted, with
  // its tolerance widened by the independently measured fan-origin spread.
  GyroidField f {3.0 * M_PI};
  compareBackends<3>(10, f, 0.0, policy, "gyroid3d");
}

void test_planar_2d(RuntimePolicy policy)
{
  PlanarField f {0.0, 1.0, 0.0, 0.5};
  compareBackends<2>(8, f, 0.0, policy, "planar2d");
}

void test_round_2d(RuntimePolicy policy)
{
  RoundField f {0.5, 0.5, 0.0, 0.25};
  compareBackends<2>(12, f, 0.0, policy, "round2d");
}

/*!
 * @brief Falsification control for E1.
 *
 * The bump crossing pre-filter classifies corners in double
 * (`fcnView(...) > m_contourVal`) while bump's FieldIntersector classifies in
 * float (`FieldIntersector::FieldType == float`).  Since float() is monotone,
 * bump's positive set is always a subset of the pre-filter's, and the dangerous
 * direction is a cell whose corners are ALL strictly above the isovalue in
 * double but where at least one rounds to exactly float(isovalue): the
 * pre-filter excludes the cell as non-crossing, while bump would have emitted a
 * fragment.
 *
 * This test constructs exactly that cell.  It is expected to FAIL on the branch
 * as-is, and to pass once the pre-filter compares in the intersector's type.
 * A regression test that has never been observed to fail is not evidence; this
 * one demonstrates that E1 can detect the specific defect.
 */
void test_float_ulp_band_falsification(RuntimePolicy policy)
{
  const int n = 4;
  const std::string fieldName = "fcn";
  const double contourVal = 1.0;

  conduit::Node mesh;
  PlanarField f {0.0, 0.0, 1.0, 0.5};
  buildStructured<3>(mesh, n, f, fieldName);

  // Overwrite the field: put every node strictly above the isovalue in double,
  // then pull one cell's corners into the band [contourVal, contourVal+ulp)
  // where float rounds them down onto float(contourVal).
  auto* fv = mesh["fields/" + fieldName + "/values"].as_float64_ptr();
  const conduit::index_t N = mesh["fields/" + fieldName + "/values"].dtype().number_of_elements();
  const int nn = n + 1;
  auto nodeAt = [&](int i, int j, int k) { return i + j * nn + k * nn * nn; };
  for(conduit::index_t i = 0; i < N; ++i)
  {
    fv[i] = contourVal + 1.0;
  }
  // Put a genuine contour in the upper part of the mesh.  Without it the entire
  // field sits above the isovalue, both paths report zero facets, and the
  // comparison below passes as 0 == 0. That would not catch a regression that
  // breaks both paths together.
  for(int k = 3; k <= n; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i)
      {
        fv[nodeAt(i, j, k)] = contourVal - 1.0;
      }
    }
  }
  // Cell (0,0,0): seven corners well above, one corner just barely above in
  // double but equal to contourVal after rounding to float.
  const double tiny = std::nextafter(contourVal, 2.0) - contourVal;  // one double ULP
  fv[nodeAt(0, 0, 0)] = contourVal + tiny;
  fv[nodeAt(1, 0, 0)] = contourVal + tiny;
  fv[nodeAt(0, 1, 0)] = contourVal + tiny;
  fv[nodeAt(1, 1, 0)] = contourVal + tiny;

  // Sanity: in double every corner is strictly above; in float the perturbed
  // ones are not.  If this fails, the platform's float rounding differs and the
  // test premise is void.
  ASSERT_GT(fv[nodeAt(0, 0, 0)], contourVal);
  ASSERT_FALSE(static_cast<float>(fv[nodeAt(0, 0, 0)]) > static_cast<float>(contourVal))
    << "premise void: the perturbed value does not collapse onto float(contourVal)";

  // Build the SAME geometry as an unstructured hex topology.  This is the
  // control that makes the test discriminating: the unstructured path's
  // pre-filter calls intersectorView.determineTableCase() (i.e. bump's own
  // float classification), while the structured path re-implements the
  // classification in double. Same nodes, same field, same cells. Any
  // difference isolates the defect to the structured pre-filter rather than
  // resting on an unverified claim about what bump "would" emit.
  conduit::Node unstructuredMesh;
  unstructuredMesh.set(mesh);
  {
    conduit::Node& topo = unstructuredMesh["topologies/mesh"];
    topo.reset();
    topo["type"] = "unstructured";
    topo["coordset"] = "coords";
    topo["elements/shape"] = "hex";
    const conduit::index_t nCells = static_cast<conduit::index_t>(n) * n * n;
    topo["elements/connectivity"].set(conduit::DataType::int64(nCells * 8));
    auto* c = topo["elements/connectivity"].as_int64_ptr();
    conduit::index_t at = 0;
    for(int k = 0; k < n; ++k)
    {
      for(int j = 0; j < n; ++j)
      {
        for(int i = 0; i < n; ++i)
        {
          c[at++] = nodeAt(i, j, k);
          c[at++] = nodeAt(i + 1, j, k);
          c[at++] = nodeAt(i + 1, j + 1, k);
          c[at++] = nodeAt(i, j + 1, k);
          c[at++] = nodeAt(i, j, k + 1);
          c[at++] = nodeAt(i + 1, j, k + 1);
          c[at++] = nodeAt(i + 1, j + 1, k + 1);
          c[at++] = nodeAt(i, j + 1, k + 1);
        }
      }
    }
  }
  conduit::Node uinfo;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(unstructuredMesh, uinfo)) << uinfo.to_yaml();

  const auto structuredRun = runBackend<3>(mesh, fieldName, contourVal, policy, /*useBump=*/true);
  const auto unstructuredRun =
    runBackend<3>(unstructuredMesh, fieldName, contourVal, policy, /*useBump=*/true);

  SLIC_INFO(axom::fmt::format(
    "[ulp_band] bump/structured: {} facets; bump/unstructured (same geometry): {} facets",
    structuredRun.facetCount,
    unstructuredRun.facetCount));

  ASSERT_GT(structuredRun.facetCount, 0)
    << "the comparison below would be vacuous: this mesh must carry a real contour";
  EXPECT_EQ(structuredRun.crossingCells, unstructuredRun.crossingCells)
    << "the two paths cut different cell sets on identical geometry";
  EXPECT_EQ(structuredRun.facetCount, unstructuredRun.facetCount)
    << "E1 falsification: the bump backend gives different answers for the same geometry "
    << "depending on whether the topology is structured or unstructured.  The structured "
    << "crossing pre-filter classifies corners in double (fcnView(...) > m_contourVal) while "
    << "the unstructured path delegates to intersectorView.determineTableCase(), which "
    << "classifies in FieldIntersector::FieldType (float).  Fix: make the structured "
    << "pre-filter compare in the intersector's type.";
}

}  // namespace

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Self-test for the ambiguity detector.
//
// A detector used to gate an assertion needs its own negative controls, or a
// silently-always-false detector would make E3 look permanently trustworthy.
// Expected counts were derived independently (by enumeration outside this
// file) before being written here:
//   - 120 of 256 sign patterns have an ambiguous face;
//   -   8 of 256 are the body-diagonal (case 4) pattern and its complement;
//   - the two classes are disjoint. Case 4 has no ambiguous face, which is
//     exactly why a face-only detector misses it;
//   - 128 of 256 total, i.e. exactly half, and complement-symmetric.
//---------------------------------------------------------------------------
TEST(quest_marching_cubes_equivalence, ambiguity_detector_selftest)
{
  auto pattern = [](int mask, bool s[8]) {
    for(int i = 0; i < 8; ++i)
    {
      s[i] = ((mask >> i) & 1) != 0;
    }
  };

  int nFace = 0, nBody = 0, nBoth = 0, nAny = 0;
  for(int mask = 0; mask < 256; ++mask)
  {
    bool s[8];
    pattern(mask, s);
    const bool fa = cellHasFaceAmbiguity3D(s);
    const bool ba = cellHasBodyDiagonalAmbiguity3D(s);
    nFace += fa ? 1 : 0;
    nBody += ba ? 1 : 0;
    nBoth += (fa && ba) ? 1 : 0;
    nAny += cellIsAmbiguous3D(s) ? 1 : 0;

    // Complement symmetry: flipping every sign cannot change whether the
    // configuration is ambiguous.
    bool c[8];
    for(int i = 0; i < 8; ++i)
    {
      c[i] = !s[i];
    }
    EXPECT_EQ(cellIsAmbiguous3D(s), cellIsAmbiguous3D(c))
      << "complement symmetry violated for sign mask " << mask;
  }

  EXPECT_EQ(nFace, 120);
  EXPECT_EQ(nBody, 8);
  EXPECT_EQ(nBoth, 0) << "face and body-diagonal ambiguity should be disjoint classes";
  EXPECT_EQ(nAny, 128);

  // Named configurations.  Corner n is (i,j,k) with i fastest: n = i + 2j + 4k.
  auto mk = [&](std::initializer_list<int> on, bool s[8]) {
    for(int i = 0; i < 8; ++i)
    {
      s[i] = false;
    }
    for(int i : on)
    {
      s[i] = true;
    }
  };
  bool s[8];

  mk({}, s);
  EXPECT_FALSE(cellIsAmbiguous3D(s)) << "all-outside must not be ambiguous (negative control)";
  mk({0, 1, 2, 3, 4, 5, 6, 7}, s);
  EXPECT_FALSE(cellIsAmbiguous3D(s)) << "all-inside must not be ambiguous (negative control)";
  mk({0}, s);
  EXPECT_FALSE(cellIsAmbiguous3D(s)) << "case 1 (single corner)";
  mk({0, 1}, s);
  EXPECT_FALSE(cellIsAmbiguous3D(s)) << "case 2 (edge pair)";
  mk({0, 1, 2, 3}, s);
  EXPECT_FALSE(cellIsAmbiguous3D(s)) << "case 8 (whole face)";
  mk({0, 3}, s);
  EXPECT_TRUE(cellHasFaceAmbiguity3D(s)) << "case 3 (face diagonal) is face-ambiguous";
  mk({0, 7}, s);
  EXPECT_FALSE(cellHasFaceAmbiguity3D(s))
    << "case 4 (body diagonal) has no ambiguous face. This is why a face-only "
       "detector misses it, and why cellHasBodyDiagonalAmbiguity3D exists";
  EXPECT_TRUE(cellHasBodyDiagonalAmbiguity3D(s)) << "case 4 (body diagonal)";
  EXPECT_TRUE(cellIsAmbiguous3D(s));

  // 2D: exactly the two checkerboards of the 16 patterns.
  int n2 = 0;
  for(int mask = 0; mask < 16; ++mask)
  {
    bool q[4];
    for(int i = 0; i < 4; ++i)
    {
      q[i] = ((mask >> i) & 1) != 0;
    }
    n2 += cellIsAmbiguous2D(q) ? 1 : 0;
  }
  EXPECT_EQ(n2, 2);
}

TEST(quest_marching_cubes_equivalence, planar_3d_seq) { test_planar_3d(RuntimePolicy::seq); }
TEST(quest_marching_cubes_equivalence, oblique_planar_3d_seq)
{
  test_oblique_planar_3d(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, round_3d_seq) { test_round_3d(RuntimePolicy::seq); }
TEST(quest_marching_cubes_equivalence, gyroid_3d_seq) { test_gyroid_3d(RuntimePolicy::seq); }
TEST(quest_marching_cubes_equivalence, planar_2d_seq) { test_planar_2d(RuntimePolicy::seq); }
TEST(quest_marching_cubes_equivalence, round_2d_seq) { test_round_2d(RuntimePolicy::seq); }
TEST(quest_marching_cubes_equivalence, uniform_and_rectilinear_seq)
{
  test_uniform_and_rectilinear(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, strided_structured_seq)
{
  test_strided_structured(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, float32_field_rejected_seq)
{
  test_float32_field_rejected(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, invalid_field_layouts_rejected_seq)
{
  test_invalid_field_layouts_rejected(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, original_elements_collision_seq)
{
  test_original_elements_collision(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, empty_contour_seq)
{
  test_empty_contour(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, float_ulp_band_falsification_seq)
{
  test_float_ulp_band_falsification(RuntimePolicy::seq);
}

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(quest_marching_cubes_equivalence, planar_3d_omp) { test_planar_3d(RuntimePolicy::omp); }
TEST(quest_marching_cubes_equivalence, round_3d_omp) { test_round_3d(RuntimePolicy::omp); }
TEST(quest_marching_cubes_equivalence, gyroid_3d_omp) { test_gyroid_3d(RuntimePolicy::omp); }
TEST(quest_marching_cubes_equivalence, round_2d_omp) { test_round_2d(RuntimePolicy::omp); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(quest_marching_cubes_equivalence, round_3d_cuda) { test_round_3d(RuntimePolicy::cuda); }
TEST(quest_marching_cubes_equivalence, gyroid_3d_cuda) { test_gyroid_3d(RuntimePolicy::cuda); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(quest_marching_cubes_equivalence, round_3d_hip) { test_round_3d(RuntimePolicy::hip); }
TEST(quest_marching_cubes_equivalence, gyroid_3d_hip) { test_gyroid_3d(RuntimePolicy::hip); }
#endif

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
