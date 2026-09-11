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
 * Legacy output duplicates vertices by facet.
 * Bump output welds vertices and can contain polygons that the adaptor triangulates.
 * The checks compare:
 *   - E1: parent cell sets
 *   - E2: vertex sets
 *   - E3: area or length when ambiguity cannot change the triangulation
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

#include "axom/quest/tests/quest_marching_cubes_testing_helpers.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <iostream>
#include <cstdint>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
namespace mctest = axom::quest::testing::marching_cubes;

using mctest::copyBlueprintToHost;
using mctest::copyBlueprintToPolicy;
using mctest::GyroidField;
using mctest::hostAllocatorID;
using mctest::PlanarField;
using mctest::RoundField;

using RuntimePolicy = axom::runtime_policy::Policy;

/*!
 * @brief Build the same box as mctest::buildStructured<3> with a uniform coordset and topology.
 *
 * @note @a n must be a power of two so the uniform and explicit coordinates agree bit for bit.
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

//! Build the same box with a rectilinear coordset and topology
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
 * The coordset and field arrays cover a padded (n+2*g)^3 window.
 * Topology offsets and strides select the n^3 real zones.
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

  // Continue the field into ghost nodes so a ghost leak produces extra facets
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
// Comparable backend results
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
BackendResult extractBackendResult(const axom::quest::MarchingCubes& mc,
                                   axom::IndexType facet_begin = 0,
                                   axom::IndexType node_begin = 0)
{
  BackendResult r;
  const axom::IndexType facet_end = mc.getContourCellCount();
  const axom::IndexType node_end = mc.getContourNodeCount();
  r.facetCount = facet_end - facet_begin;
  r.nodeCount = node_end - node_begin;

  const axom::Array<double, 2> coordsHost(mc.getContourNodeCoords(), hostAllocatorID());
  const axom::Array<axom::IndexType, 2> cornersHost(mc.getContourFacetCorners(), hostAllocatorID());
  const axom::Array<axom::IndexType> parentsHost(mc.getContourFacetParents(), hostAllocatorID());
  const auto coords = coordsHost.view();
  const auto corners = cornersHost.view();
  const auto parents = parentsHost.view();

  for(axom::IndexType f = facet_begin; f < facet_end; ++f)
  {
    r.crossingCells.insert(parents[f]);
  }

  for(axom::IndexType v = node_begin; v < node_end; ++v)
  {
    axom::primal::Point<double, 3> p {};
    p[0] = coords(v, 0);
    p[1] = coords(v, 1);
    p[2] = (DIM == 3) ? coords(v, 2) : 0.0;
    r.vertices.push_back(p);
  }

  for(axom::IndexType f = facet_begin; f < facet_end; ++f)
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

template <int DIM>
BackendResult runBackend(const conduit::Node& mesh,
                         const std::string& field_name,
                         double contour_value,
                         RuntimePolicy policy,
                         bool use_bump,
                         conduit::Node* bump_blueprint = nullptr)
{
  namespace quest = axom::quest;

  const int allocator_id = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocator_id, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(use_bump);

  conduit::Node exec_mesh;
  copyBlueprintToPolicy(exec_mesh, mesh, policy, allocator_id);
  mc.setMesh(exec_mesh, "mesh");
  mc.setFunctionField(field_name);
  mc.computeIsocontour(contour_value);

  if(use_bump && bump_blueprint != nullptr)
  {
    conduit::Node bump_blueprint_exec;
    mc.populateContourMeshBlueprint(bump_blueprint_exec);
    copyBlueprintToHost(*bump_blueprint, bump_blueprint_exec);
  }

  return extractBackendResult<DIM>(mc);
}

template <int DIM>
std::array<BackendResult, 3> runAccumulatedBackend(const conduit::Node& mesh,
                                                   const std::array<std::string, 3>& field_names,
                                                   RuntimePolicy policy,
                                                   bool use_bump,
                                                   std::array<conduit::Node, 3>* bump_blueprints = nullptr)
{
  namespace quest = axom::quest;

  const int allocator_id = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocator_id, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(use_bump);

  conduit::Node exec_mesh;
  copyBlueprintToPolicy(exec_mesh, mesh, policy, allocator_id);
  mc.setMesh(exec_mesh, "mesh");

  std::array<BackendResult, 3> results;
  axom::IndexType facet_begin = 0;
  axom::IndexType node_begin = 0;
  for(int field = 0; field < 3; ++field)
  {
    mc.setFunctionField(field_names[field]);
    mc.computeIsocontour(0.0);
    results[field] = extractBackendResult<DIM>(mc, facet_begin, node_begin);
    if(use_bump && bump_blueprints != nullptr)
    {
      conduit::Node bump_blueprint_exec;
      mc.populateContourMeshBlueprint(bump_blueprint_exec);
      copyBlueprintToHost((*bump_blueprints)[field], bump_blueprint_exec);
    }
    facet_begin = mc.getContourFacetCount();
    node_begin = mc.getContourNodeCount();
  }
  return results;
}

//---------------------------------------------------------------------------
// Two-sided Hausdorff comparison. Search adjacent hash cells so nearby points
// still match when they straddle a cell boundary.
//---------------------------------------------------------------------------

using CellKey = std::int64_t;

CellKey cellKey(std::int64_t i, std::int64_t j, std::int64_t k)
{
  // Hash collisions remain in the bucket and are resolved by distance checks
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

//! @brief Return the one-sided Hausdorff distance and its source point index
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
// Ambiguous cells. Case tables can triangulate the same vertices differently
// for a checkerboard face or a body-diagonal minority pair. The latter is case
// 4 and has no ambiguous face. Corner n is (i,j,k), with n = i + 2j + 4k.
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

//! 2D corners are cyclic 0,1,3,2. Only the two
//! checkerboard patterns are ambiguous.
bool cellIsAmbiguous2D(const bool s[4])
{
  return (s[0] == s[3]) && (s[1] == s[2]) && (s[0] != s[1]);
}

/*!
 * @brief Count ambiguous cells.
 *
 * @note The corner test uses `>=` to match MarchingCubesImpl::computeCrossingCase and the adjusted Bump isovalue
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
// Fan-triangulation sensitivity. The adaptor fans each Bump polygon from corner 0.
// For a non-planar polygon, area depends on the fan origin. Compare fans from corners 0 and 1
// and use their relative spread as the measure tolerance.
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

//! Measure how much the Bump polygon area depends on the fan origin.
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
// Backend comparison
//---------------------------------------------------------------------------

template <int DIM>
void compareBackendResults(const BackendResult& legacy,
                           const BackendResult& bump,
                           const std::string& label,
                           axom::IndexType ambiguous,
                           const FanSensitivity& fan = {},
                           double vertex_tolerance = 1.0e-5,
                           double measure_relative_tolerance = 1.0e-5)
{
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

  std::vector<axom::IndexType> only_legacy;
  std::vector<axom::IndexType> only_bump;
  std::set_difference(legacy.crossingCells.begin(),
                      legacy.crossingCells.end(),
                      bump.crossingCells.begin(),
                      bump.crossingCells.end(),
                      std::back_inserter(only_legacy));
  std::set_difference(bump.crossingCells.begin(),
                      bump.crossingCells.end(),
                      legacy.crossingCells.begin(),
                      legacy.crossingCells.end(),
                      std::back_inserter(only_bump));

  EXPECT_TRUE(only_legacy.empty()) << "E1 [" << label << "]: " << only_legacy.size()
                                   << " cells produce facets in legacy but not bump (first: "
                                   << (only_legacy.empty() ? -1 : only_legacy.front()) << ").";
  EXPECT_TRUE(only_bump.empty()) << "E1 [" << label << "]: " << only_bump.size()
                                 << " cells produce facets in bump but not legacy (first: "
                                 << (only_bump.empty() ? -1 : only_bump.front()) << ").";

  const PointLocator legacy_locator(legacy.vertices, vertex_tolerance);
  const PointLocator bump_locator(bump.vertices, vertex_tolerance);
  std::size_t arg_max = 0;
  const double bump_to_legacy = oneSidedHausdorff(bump.vertices, legacy_locator, arg_max);
  EXPECT_LT(bump_to_legacy, vertex_tolerance)
    << "E2 [" << label << "]: a bump contour vertex has no legacy counterpart within tolerance"
    << " (worst distance " << bump_to_legacy << " at bump vertex " << arg_max << ").";

  const double legacy_to_bump = oneSidedHausdorff(legacy.vertices, bump_locator, arg_max);
  EXPECT_LT(legacy_to_bump, vertex_tolerance)
    << "E2 [" << label << "]: a legacy contour vertex has no bump counterpart within tolerance"
    << " (worst distance " << legacy_to_bump << " at legacy vertex " << arg_max << ").";

  EXPECT_LT(bump.nodeCount, legacy.nodeCount)
    << "[" << label << "]: bump node count is not smaller than legacy's. Welding regressed.";

  const double relative_difference =
    std::abs(bump.measure - legacy.measure) / std::max(legacy.measure, 1.0e-300);
  const double measure_tolerance = std::max(measure_relative_tolerance, fan.relSpread);
  EXPECT_LT(fan.relSpread, 0.1)
    << "[" << label << "]: fan-origin spread is too large for a useful measure comparison.";
  if(ambiguous == 0)
  {
    EXPECT_LT(relative_difference, measure_tolerance)
      << "E3 [" << label << "]: contour measure differs without an ambiguous cell to explain it";
  }
}

template <int DIM, typename Field>
void compareBackends(int n,
                     const Field& f,
                     double contourVal,
                     RuntimePolicy policy,
                     const std::string& label,
                     double vertexTol = 1.0e-5,
                     // Bump interpolates edge crossings in float.
                     // Passing cases have measured relative errors from 5e-9 to 3e-6.
                     double measureRelTol = 1.0e-5)
{
  const std::string fieldName = "fcn";
  conduit::Node mesh;
  mctest::buildStructured<DIM>(mesh, n, f, fieldName);

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

  compareBackendResults<DIM>(legacy, bump, label, ambiguous, fan, vertexTol, measureRelTol);
}

/*!
 * @brief An isovalue outside the data range is valid input, not an error.
 *
 * This covers Blueprint population and relinquishment after an extraction with no crossing cells.
 */
void test_empty_contour(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  conduit::Node mesh;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  mctest::buildStructured<3>(mesh, 6, f, "fcn");

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
 * @brief Compare uniform and rectilinear input with equivalent explicit input.
 *
 * The legacy backend cannot read uniform or rectilinear meshes, so the explicit
 * Bump result is the reference for all three representations of the same box.
 */
void test_uniform_and_rectilinear(RuntimePolicy policy)
{
  const int n = 8;  // power of two: see buildUniform3D's note on coordinate agreement
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node structured, uniform, rectilinear;
  mctest::buildStructured<3>(structured, n, f, fieldName);
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
 * @brief Check that MarchingCubes rejects a float32 function field.
 *
 * The structured pre-filter requires float64 values.
 * Reinterpreting float32 values as float64 produces an invalid crossing-cell set.
 */
void test_float32_field_rejected(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  const int n = 6;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  conduit::Node mesh;
  mctest::buildStructured<3>(mesh, n, f, "fcn");

  // Rewrite the function field as float32, keeping the same values
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

  // Route SLIC output to stderr in the child so gtest can match the diagnostic
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
 * @brief Check that MarchingCubes rejects unsupported strided field layouts.
 *
 * Bump's flat field view cannot represent a permuted field or one whose offsets
 * and strides differ from the topology.
 */
void test_invalid_field_layouts_rejected(RuntimePolicy policy)
{
  constexpr int n = 6;
  constexpr int pad = 2;
  constexpr int nnPad = n + 1 + 2 * pad;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};

  conduit::Node permuted;
  buildStridedStructured3D(permuted, n, pad, f, "fcn");
  permuted["fields/fcn/strides"].set(std::vector<conduit::int32> {nnPad * nnPad, nnPad, 1});
  // The diagnostic reports both the field and topology layouts.
  expectBumpFieldLayoutRejected(permuted, policy, "but its topology has strides");

  conduit::Node mismatched;
  buildStridedStructured3D(mismatched, n, pad, f, "fcn");
  mismatched["fields/fcn/offsets"].set(std::vector<conduit::int32> {pad + 1, pad, pad});
  expectBumpFieldLayoutRejected(mismatched, policy, "but its topology has offsets");
}

/*!
 * @brief Preserve parent ids when the input has an "originalElements" field.
 *
 * TableBasedExtractor otherwise propagates the input field instead of writing
 * source zone indices. MarchingCubes uses a private field name to avoid this collision.
 */
void test_original_elements_collision(RuntimePolicy policy)
{
  const int n = 8;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node clean, poisoned;
  mctest::buildStructured<3>(clean, n, f, fieldName);
  poisoned.set(clean);

  // Use negative decoy values that cannot be valid zone indices
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
 * @brief Compare strided and compact representations of the same mesh.
 *
 * Ghost values continue the field outside the real zone range,
 * so ignoring the topology offsets produces extra facets.
 */
void test_strided_structured(RuntimePolicy policy)
{
  const int n = 8;
  const int pad = 2;  // ghost layers
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node compact, strided;
  mctest::buildStructured<3>(compact, n, f, fieldName);
  buildStridedStructured3D(strided, n, pad, f, fieldName);

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(strided, info)) << info.to_yaml();

  // Legacy results establish that the compact and strided fixtures agree.
  // Parent-id bounds catch ghost cells even when facet counts happen to match.
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
    << "legacy results differ between strided and compact input";

  EXPECT_EQ(bumpStrided.crossingCells, legacyCompact.crossingCells)
    << "Bump strided input cut a different cell set than compact input";

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
  // An axis-aligned plane has no ambiguous cells, so compare its area
  PlanarField f {{0.0, 0.0, 1.0}, 0.5};
  compareBackends<3>(8, f, 0.0, policy, "planar3d");
}

void test_oblique_planar_3d(RuntimePolicy policy)
{
  // An oblique plane exercises more case-table entries without ambiguity
  const double s = 1.0 / std::sqrt(1.0 + 0.16 + 1.44);
  PlanarField f {{1.0 * s, 0.4 * s, 1.2 * s}, 1.3 * s};
  compareBackends<3>(8, f, 0.0, policy, "oblique_planar3d");
}

void test_round_3d(RuntimePolicy policy)
{
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  compareBackends<3>(12, f, 0.0, policy, "round3d");
}

void test_gyroid_3d(RuntimePolicy policy)
{
  // Use the measured fan-origin spread as the area tolerance
  GyroidField f {3.0 * M_PI};
  compareBackends<3>(10, f, 0.0, policy, "gyroid3d");
}

void test_planar_2d(RuntimePolicy policy)
{
  PlanarField f {{0.0, 1.0, 0.0}, 0.5};
  compareBackends<2>(8, f, 0.0, policy, "planar2d");
}

void test_round_2d(RuntimePolicy policy)
{
  RoundField f {{0.5, 0.5, 0.0}, 0.25};
  compareBackends<2>(12, f, 0.0, policy, "round2d");
}

template <int DIM>
void test_accumulated_fields(RuntimePolicy policy)
{
  constexpr int n = 12;
  PlanarField plane {{0.47, 0.43, 0.39}, {1.0, 0.4, 1.2}};
  RoundField round {{0.5, 0.5, DIM == 3 ? 0.5 : 0.0}, 0.27};
  GyroidField gyroid {{3.0, 3.0, DIM == 3 ? 1.5 : 0.0}};
  const std::array<std::string, 3> field_names {{"plane", "round", "gyroid"}};

  conduit::Node mesh;
  mctest::buildStructured<DIM>(mesh, n, plane, field_names[0]);
  mctest::addVertexField<DIM>(mesh, round, field_names[1]);
  mctest::addVertexField<DIM>(mesh, gyroid, field_names[2]);

  const auto legacy = runAccumulatedBackend<DIM>(mesh, field_names, policy, false);
  std::array<conduit::Node, 3> bump_blueprints;
  const auto bump = runAccumulatedBackend<DIM>(mesh, field_names, policy, true, &bump_blueprints);
  const std::array<axom::IndexType, 3> ambiguous {{countAmbiguousCells<DIM>(n, plane, 0.0),
                                                   countAmbiguousCells<DIM>(n, round, 0.0),
                                                   countAmbiguousCells<DIM>(n, gyroid, 0.0)}};

  for(int field = 0; field < 3; ++field)
  {
    FanSensitivity fan;
    if constexpr(DIM == 3)
    {
      ASSERT_EQ(bump_blueprints[field].number_of_children(), 1);
      fan = measureFanSensitivity(bump_blueprints[field].child(0));
    }
    compareBackendResults<DIM>(legacy[field],
                               bump[field],
                               "accumulated_" + field_names[field],
                               ambiguous[field],
                               fan);
  }
}

/*!
 * @brief Compare float and double classification near the isovalue.
 *
 * Values just above the isovalue in double can round to the isovalue in float.
 * Structured and unstructured forms of the same mesh must still select the same crossing cells.
 */
void test_float_ulp_band_falsification(RuntimePolicy policy)
{
  const int n = 4;
  const std::string fieldName = "fcn";
  const double contourVal = 1.0;

  conduit::Node mesh;
  PlanarField f {{0.0, 0.0, 1.0}, 0.5};
  mctest::buildStructured<3>(mesh, n, f, fieldName);

  // Put one cell's corners just above the isovalue in double but equal to it after conversion to float
  auto* fv = mesh["fields/" + fieldName + "/values"].as_float64_ptr();
  const conduit::index_t N = mesh["fields/" + fieldName + "/values"].dtype().number_of_elements();
  const int nn = n + 1;
  auto nodeAt = [&](int i, int j, int k) { return i + j * nn + k * nn * nn; };
  for(conduit::index_t i = 0; i < N; ++i)
  {
    fv[i] = contourVal + 1.0;
  }
  // Add a separate contour so equal empty results cannot pass the test
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
  // Four corners of cell (0,0,0) lie in the double-to-float rounding gap
  const double tiny = std::nextafter(contourVal, 2.0) - contourVal;  // one double ULP
  fv[nodeAt(0, 0, 0)] = contourVal + tiny;
  fv[nodeAt(1, 0, 0)] = contourVal + tiny;
  fv[nodeAt(0, 1, 0)] = contourVal + tiny;
  fv[nodeAt(1, 1, 0)] = contourVal + tiny;

  // Confirm the rounding behavior required by the test
  ASSERT_GT(fv[nodeAt(0, 0, 0)], contourVal);
  ASSERT_FALSE(static_cast<float>(fv[nodeAt(0, 0, 0)]) > static_cast<float>(contourVal))
    << "perturbed value did not round to float(contourVal)";

  // Reuse the coordinates and field with an unstructured hex topology.
  // Its pre-filter uses Bump's float classification.
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

  ASSERT_GT(structuredRun.facetCount, 0) << "test mesh must contain a contour";
  EXPECT_EQ(structuredRun.crossingCells, unstructuredRun.crossingCells)
    << "the two paths cut different cell sets on identical geometry";
  EXPECT_EQ(structuredRun.facetCount, unstructuredRun.facetCount)
    << "structured and unstructured paths disagree near the float isovalue";
}

}  // namespace

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Self-test for the ambiguity detector.
//
// Of the 256 sign patterns, 120 have face ambiguity and 8 have body-diagonal ambiguity.
// The classes are disjoint and invariant under sign reversal.
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

  // Named configurations. Corner n is (i,j,k), with n = i + 2j + 4k.
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
  EXPECT_FALSE(cellHasFaceAmbiguity3D(s)) << "case 4 has no ambiguous face";
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
TEST(quest_marching_cubes_equivalence, accumulated_fields_2d_seq)
{
  test_accumulated_fields<2>(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, accumulated_fields_3d_seq)
{
  test_accumulated_fields<3>(RuntimePolicy::seq);
}
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
TEST(quest_marching_cubes_equivalence, accumulated_fields_2d_omp)
{
  test_accumulated_fields<2>(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_equivalence, accumulated_fields_3d_omp)
{
  test_accumulated_fields<3>(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_equivalence, strided_structured_omp)
{
  test_strided_structured(RuntimePolicy::omp);
}
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(quest_marching_cubes_equivalence, round_3d_cuda) { test_round_3d(RuntimePolicy::cuda); }
TEST(quest_marching_cubes_equivalence, gyroid_3d_cuda) { test_gyroid_3d(RuntimePolicy::cuda); }
TEST(quest_marching_cubes_equivalence, accumulated_fields_2d_cuda)
{
  test_accumulated_fields<2>(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_equivalence, accumulated_fields_3d_cuda)
{
  test_accumulated_fields<3>(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_equivalence, strided_structured_cuda)
{
  test_strided_structured(RuntimePolicy::cuda);
}
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(quest_marching_cubes_equivalence, round_3d_hip) { test_round_3d(RuntimePolicy::hip); }
TEST(quest_marching_cubes_equivalence, gyroid_3d_hip) { test_gyroid_3d(RuntimePolicy::hip); }
TEST(quest_marching_cubes_equivalence, accumulated_fields_2d_hip)
{
  test_accumulated_fields<2>(RuntimePolicy::hip);
}
TEST(quest_marching_cubes_equivalence, accumulated_fields_3d_hip)
{
  test_accumulated_fields<3>(RuntimePolicy::hip);
}
TEST(quest_marching_cubes_equivalence, strided_structured_hip)
{
  test_strided_structured(RuntimePolicy::hip);
}
#endif

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
