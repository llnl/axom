// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file quest_marching_cubes_bump.cpp
 *
 * @brief Validation tests for the bump-based backend of quest::MarchingCubes,
 * covering structured AND unstructured single-shape (quad/hex) input
 * on every available execution space.
 *
 * This exercises:
 *  - The bump backend (MarchingCubes::setUseBumpBackend(true)).
 *  - Unstructured single-shape quad (2D) and hex (3D) topologies.
 *  - An explicit edge-manifoldness check on the extracted contour.
 *    This is relevant to the saddle-ambiguity behavior of the original MC tables.
 *    it should pass on the smooth analytic fields used here regardless of backend.
 *
 * Verification oracles (all backend-agnostic; applied to bump output):
 *   O1. On-surface value: every output facet node, evaluated in the analytic field,
 *       equals the contour value within tolerance (exact for planar).
 *   O2. Parent containment: each facet's parent cell id is in range,
 *       and the t centroid lies within the parent cell's axis-aligned bounds.
 *   O3. Edge manifoldness on bump's welded Blueprint output:
 *       every contour edge (3D) is shared by exactly 1 or 2 facets and no edge is used 3+ times.
 *       A closed smooth surface interior to the domain should have all-interior edges shared exactly twice;
 *       boundary-clipped edges may be shared once.
 */

#include "axom/config.hpp"

#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)

  #include "axom/core.hpp"
  #include "axom/bump/utilities/conduit_memory.hpp"
  #include "axom/primal.hpp"
  #include "axom/quest/MarchingCubes.hpp"
  #include "axom/quest/util/mesh_helpers.hpp"
  #include "axom/sidre.hpp"
  #include "axom/spin/MortonIndex.hpp"
  #include "axom/mint/mesh/UnstructuredMesh.hpp"

  #include "conduit_blueprint.hpp"

  #include "gtest/gtest.h"

  #include <cmath>
  #include <cstdint>
  #include <map>
  #include <unordered_map>
  #include <utility>

namespace
{
using RuntimePolicy = axom::runtime_policy::Policy;
using Point3D = axom::primal::Point<double, 3>;
using BoundingBox3D = axom::primal::BoundingBox<double, 3>;
using QuantizedPoint3D = axom::primal::Point<std::int64_t, 3>;

int hostAllocatorID() { return axom::execution_space<axom::SEQ_EXEC>::allocatorID(); }

void copyBlueprintToPolicy(conduit::Node& dst,
                           const conduit::Node& src,
                           RuntimePolicy policy,
                           int allocatorID)
{
  namespace bputils = axom::bump::utilities;

  // The bump backend reads Blueprint arrays in the execution space associated
  // with the runtime policy. Keep mesh construction host-side, then copy the
  // finished Blueprint tree into policy-compatible memory for MarchingCubes.
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
// Analytic fields (mirror the example's PlanarTestStrategy / RoundTestStrategy)
//---------------------------------------------------------------------------

//! @brief Signed distance to a plane through @a origin with unit normal @a n.
struct PlanarField
{
  using PointType = axom::primal::Point<double, 3>;
  using VectorType = axom::primal::Vector<double, 3>;
  using PlaneType = axom::primal::Plane<double, 3>;

  PlanarField(const PointType& origin, const VectorType& normal) : plane(normal, origin) { }

  double operator()(double x, double y, double z) const
  {
    return plane.signedDistance(PointType {x, y, z});
  }

  PlaneType plane;
};

//! @brief Signed distance to a sphere/circle.
struct RoundField
{
  using PointType = axom::primal::Point<double, 3>;
  using SphereType = axom::primal::Sphere<double, 3>;

  RoundField(const PointType& center, double radius) : sphere(center, radius) { }

  double operator()(double x, double y, double z) const
  {
    return sphere.computeSignedDistance(PointType {x, y, z});
  }

  SphereType sphere;
};

//---------------------------------------------------------------------------
// O3 helper: edge-manifoldness over a welded triangle/segment soup.
//---------------------------------------------------------------------------

/*!
 * @brief Count, for a 3D triangle soup, how many facets use each undirected edge
 * after welding coincident vertices (quantized hash).
 * Returns the max multiplicity and the count of edges used 3+ times.
 */
struct EdgeManifoldResult
{
  int maxMultiplicity = 0;
  axom::IndexType edgesUsed3PlusTimes = 0;
  axom::IndexType boundaryEdges = 0;  // used exactly once
  axom::IndexType interiorEdges = 0;  // used exactly twice
};

EdgeManifoldResult checkEdgeManifold3D(const axom::ArrayView<const double, 2>& nodeCoords,
                                       const axom::ArrayView<const axom::IndexType, 2>& facetCorners,
                                       double weldTol)
{
  const double inv = 1.0 / weldTol;
  auto quantize = [inv](double v) { return static_cast<std::int64_t>(std::llround(v * inv)); };

  // The legacy MarchingCubes arrays duplicate coordinates per facet.
  // Quantize coordinates here to recover welded vertex ids for the
  // helper self-test without changing the output contract.
  std::unordered_map<QuantizedPoint3D, axom::IndexType, axom::spin::PointHash<std::int64_t>> vmap;
  const axom::IndexType nFacets = facetCorners.shape()[0];

  auto weldedId = [&](axom::IndexType row) {
    QuantizedPoint3D key {quantize(nodeCoords(row, 0)),
                          quantize(nodeCoords(row, 1)),
                          quantize(nodeCoords(row, 2))};
    auto it = vmap.find(key);
    if(it != vmap.end())
    {
      return it->second;
    }
    const axom::IndexType id = static_cast<axom::IndexType>(vmap.size());
    vmap.emplace(key, id);
    return id;
  };

  std::map<std::pair<axom::IndexType, axom::IndexType>, int> edgeUse;
  for(axom::IndexType f = 0; f < nFacets; ++f)
  {
    axom::IndexType v[3];
    for(int c = 0; c < 3; ++c)
    {
      v[c] = weldedId(facetCorners(f, c));
    }
    for(int e = 0; e < 3; ++e)
    {
      axom::IndexType a = v[e], b = v[(e + 1) % 3];
      if(a == b)
      {
        continue;  // degenerate edge; ignore
      }
      if(a > b)
      {
        std::swap(a, b);
      }
      edgeUse[{a, b}]++;
    }
  }

  EdgeManifoldResult res;
  for(const auto& kv : edgeUse)
  {
    res.maxMultiplicity = std::max(res.maxMultiplicity, kv.second);
    if(kv.second == 1)
    {
      res.boundaryEdges++;
    }
    else if(kv.second == 2)
    {
      res.interiorEdges++;
    }
    else if(kv.second >= 3)
    {
      res.edgesUsed3PlusTimes++;
    }
  }
  return res;
}

EdgeManifoldResult checkBlueprintEdgeManifold3D(const conduit::Node& contourDom)
{
  // The Blueprint accessor exposes bump's native welded polygonal contour,
  // so this path can count edges directly from connectivity without coordinate
  // re-welding or fan-triangulating polygons.
  const conduit::Node& n_topo = contourDom.fetch_existing("topologies").child(0);
  const conduit::Node& n_elems = n_topo.fetch_existing("elements");
  const auto sizes = n_elems.fetch_existing("sizes").as_index_t_accessor();
  const auto offsets = n_elems.fetch_existing("offsets").as_index_t_accessor();
  const auto conn = n_elems.fetch_existing("connectivity").as_index_t_accessor();

  std::map<std::pair<axom::IndexType, axom::IndexType>, int> edgeUse;
  const conduit::index_t nZones = sizes.number_of_elements();
  for(conduit::index_t z = 0; z < nZones; ++z)
  {
    const auto nCorners = static_cast<axom::IndexType>(sizes[z]);
    const auto offset = static_cast<axom::IndexType>(offsets[z]);
    for(axom::IndexType e = 0; e < nCorners; ++e)
    {
      axom::IndexType a = static_cast<axom::IndexType>(conn[offset + e]);
      axom::IndexType b = static_cast<axom::IndexType>(conn[offset + ((e + 1) % nCorners)]);
      if(a == b)
      {
        continue;  // degenerate edge; ignore
      }
      if(a > b)
      {
        std::swap(a, b);
      }
      edgeUse[{a, b}]++;
    }
  }

  EdgeManifoldResult res;
  for(const auto& kv : edgeUse)
  {
    res.maxMultiplicity = std::max(res.maxMultiplicity, kv.second);
    if(kv.second == 1)
    {
      res.boundaryEdges++;
    }
    else if(kv.second == 2)
    {
      res.interiorEdges++;
    }
    else if(kv.second >= 3)
    {
      res.edgesUsed3PlusTimes++;
    }
  }
  return res;
}

//---------------------------------------------------------------------------
// Mesh builders
//---------------------------------------------------------------------------

//! @brief Identity warp (no displacement); the default coordinate map.
struct NoWarp
{
  void operator()(double& /*x*/, double& /*y*/, double& /*z*/) const { }
};

// A smooth sinusoidal shear that makes the hexes genuinely non-axis-aligned (curvilinear)
// while keeping them valid (small displacement, positive Jacobian).
// Boundaries z=0 and z=1 are pinned so the [0,1]^3 box is preserved enough
// for the interior sphere to stay strictly inside.
struct SinusoidalWarp
{
  double amp;
  void operator()(double& x, double& y, double& z) const
  {
    const double sx = std::sin(M_PI * x), sy = std::sin(M_PI * y), sz = std::sin(M_PI * z);
    // Displace interior nodes; boundary planes have a zero sin factor,
    // so the outer box edges/corners stay put.
    x += amp * sy * sz;
    y += amp * sx * sz;
    z += amp * sx * sy;
  }
};

//! @brief Build a single-domain structured explicit 3D mesh [0,1]^3
//!   with @a n cells per side and the analytic field @a f sampled at nodes.
template <typename Field, typename Warp = NoWarp>
void buildStructured3D(conduit::Node& mesh,
                       int n,
                       const Field& f,
                       const std::string& fieldName,
                       const Warp& warp = Warp {})
{
  const int nn = n + 1;
  const conduit::index_t N = static_cast<conduit::index_t>(nn) * nn * nn;

  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "explicit";
  cs["values/x"].set(conduit::DataType::float64(N));
  cs["values/y"].set(conduit::DataType::float64(N));
  cs["values/z"].set(conduit::DataType::float64(N));
  auto* x = cs["values/x"].as_float64_ptr();
  auto* y = cs["values/y"].as_float64_ptr();
  auto* z = cs["values/z"].as_float64_ptr();

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "structured";
  topo["coordset"] = "coords";
  topo["elements/dims/i"] = n;
  topo["elements/dims/j"] = n;
  topo["elements/dims/k"] = n;

  conduit::Node& fld = mesh["fields/" + fieldName];
  fld["topology"] = "mesh";
  fld["association"] = "vertex";
  fld["values"].set(conduit::DataType::float64(N));
  auto* fv = fld["values"].as_float64_ptr();

  // i-fastest node ordering, matching bump's StructuredIndexing.
  conduit::index_t idx = 0;
  for(int k = 0; k < nn; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i, ++idx)
      {
        double px = double(i) / n, py = double(j) / n, pz = double(k) / n;
        // Apply the (smooth) coordinate warp, then sample the analytic field at the resulting physical location
        // so the on-surface oracle (which evaluates f at output node coords) stays consistent for warped hexes.
        warp(px, py, pz);
        x[idx] = px;
        y[idx] = py;
        z[idx] = pz;
        fv[idx] = f(px, py, pz);
      }
    }
  }
}

void addStructuredMask3D(conduit::Node& mesh,
                         int n,
                         const std::string& maskFieldName,
                         int selectedValue,
                         int rejectedValue)
{
  const conduit::index_t nCells = static_cast<conduit::index_t>(n) * n * n;

  conduit::Node& mask = mesh["fields/" + maskFieldName];
  mask["topology"] = "mesh";
  mask["association"] = "element";
  mask["values"].set(conduit::DataType::int32(nCells));
  auto* values = mask["values"].as_int32_ptr();

  // Build an element-associated mask that selects the lower half of the
  // structured mesh in k. The masked test below then verifies that bump's
  // selectedZones path emits contour facets only from cells with this value.
  conduit::index_t idx = 0;
  for(int k = 0; k < n; ++k)
  {
    for(int j = 0; j < n; ++j)
    {
      for(int i = 0; i < n; ++i, ++idx)
      {
        AXOM_UNUSED_VAR(i);
        AXOM_UNUSED_VAR(j);
        values[idx] = (k < n / 2) ? selectedValue : rejectedValue;
      }
    }
  }
}

Point3D contourFacetCentroid3D(const axom::ArrayView<const double, 2>& nodeCoords,
                               const axom::ArrayView<const axom::IndexType, 2>& facetCorners,
                               axom::IndexType facetIndex)
{
  Point3D centroid {};
  for(int c = 0; c < 3; ++c)
  {
    const axom::IndexType nodeIndex = facetCorners(facetIndex, c);
    for(int d = 0; d < 3; ++d)
    {
      centroid[d] += nodeCoords(nodeIndex, d);
    }
  }
  for(int d = 0; d < 3; ++d)
  {
    centroid[d] /= 3.0;
  }
  return centroid;
}

void addCoordsetPointToBounds(const conduit::Node& n_values,
                              axom::IndexType nodeIndex,
                              BoundingBox3D& bounds)
{
  const auto x = n_values.fetch_existing("x").as_float64_accessor();
  const auto y = n_values.fetch_existing("y").as_float64_accessor();
  const auto z = n_values.fetch_existing("z").as_float64_accessor();
  bounds.addPoint(Point3D {static_cast<double>(x[nodeIndex]),
                           static_cast<double>(y[nodeIndex]),
                           static_cast<double>(z[nodeIndex])});
}

bool structuredCellBounds3D(const conduit::Node& mesh, axom::IndexType cellIndex, BoundingBox3D& bounds)
{
  const conduit::Node& topo = mesh.fetch_existing("topologies/mesh");
  const axom::IndexType ni =
    static_cast<axom::IndexType>(topo.fetch_existing("elements/dims/i").to_value());
  const axom::IndexType nj =
    static_cast<axom::IndexType>(topo.fetch_existing("elements/dims/j").to_value());
  const axom::IndexType nk =
    static_cast<axom::IndexType>(topo.fetch_existing("elements/dims/k").to_value());
  const axom::IndexType nCells = ni * nj * nk;
  if(cellIndex < 0 || cellIndex >= nCells)
  {
    return false;
  }

  const axom::IndexType i = cellIndex % ni;
  const axom::IndexType j = (cellIndex / ni) % nj;
  const axom::IndexType k = cellIndex / (ni * nj);

  const std::string coordsetName = topo.fetch_existing("coordset").as_string();
  const conduit::Node& n_values = mesh.fetch_existing("coordsets/" + coordsetName + "/values");

  const axom::IndexType nni = ni + 1;
  const axom::IndexType nnj = nj + 1;
  // Structured test meshes use i-fastest node ordering, matching the builder above
  // and bump's structured indexing. Enumerate the eight logical corners of the
  // reported parent cell to form a physical-space AABB.
  auto nodeIndex = [=](axom::IndexType ii, axom::IndexType jj, axom::IndexType kk) {
    return ii + jj * nni + kk * nni * nnj;
  };

  for(axom::IndexType dk = 0; dk <= 1; ++dk)
  {
    for(axom::IndexType dj = 0; dj <= 1; ++dj)
    {
      for(axom::IndexType di = 0; di <= 1; ++di)
      {
        addCoordsetPointToBounds(n_values, nodeIndex(i + di, j + dj, k + dk), bounds);
      }
    }
  }
  return bounds.isValid();
}

bool unstructuredHexCellBounds3D(const conduit::Node& mesh,
                                 axom::IndexType cellIndex,
                                 BoundingBox3D& bounds)
{
  const conduit::Node& topo = mesh.fetch_existing("topologies/mesh");
  if(topo.fetch_existing("elements/shape").as_string() != std::string("hex"))
  {
    return false;
  }

  constexpr axom::IndexType HEX_NODES = 8;
  const auto conn = topo.fetch_existing("elements/connectivity").as_index_t_accessor();
  const axom::IndexType firstConn = cellIndex * HEX_NODES;
  if(cellIndex < 0 || firstConn + HEX_NODES > conn.number_of_elements())
  {
    return false;
  }

  const std::string coordsetName = topo.fetch_existing("coordset").as_string();
  const conduit::Node& n_values = mesh.fetch_existing("coordsets/" + coordsetName + "/values");

  for(axom::IndexType c = 0; c < HEX_NODES; ++c)
  {
    const axom::IndexType nodeIndex = static_cast<axom::IndexType>(conn[firstConn + c]);
    addCoordsetPointToBounds(n_values, nodeIndex, bounds);
  }
  return bounds.isValid();
}

bool parentCellBounds3D(const conduit::Node& mesh, axom::IndexType cellIndex, BoundingBox3D& bounds)
{
  const std::string topoType = mesh.fetch_existing("topologies/mesh/type").as_string();
  if(topoType == "structured")
  {
    return structuredCellBounds3D(mesh, cellIndex, bounds);
  }
  if(topoType == "unstructured")
  {
    return unstructuredHexCellBounds3D(mesh, cellIndex, bounds);
  }
  return false;
}

//---------------------------------------------------------------------------
// Core check: run bump-backed MarchingCubes and apply tests O1, O2 and O3
//---------------------------------------------------------------------------

template <typename Field>
void runAndVerify3D(conduit::Node& mesh,
                    const Field& f,
                    double contourVal,
                    RuntimePolicy policy,
                    const std::string& fieldName,
                    bool expectClosedInterior,
                    double analyticSurfaceTol,
                    const std::string& maskFieldName = {},
                    int maskVal = 1,
                    axom::quest::MarchingCubesRobustnessPolicy robustness =
                      axom::quest::MarchingCubesRobustnessPolicy::standard)
{
  namespace quest = axom::quest;

  const int allocatorID = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(true);
  mc.setRobustnessPolicy(robustness);

  // MarchingCubes' public input contract is multi-domain.
  // Keep the wrapped node alive through computeIsocontour(),
  // since the single-domain objects cache pointers into it.
  conduit::Node mdMesh;
  mdMesh.append().set(mesh);
  conduit::Node execMdMesh;
  copyBlueprintToPolicy(execMdMesh, mdMesh, policy, allocatorID);
  mc.setMesh(execMdMesh, "mesh", maskFieldName);
  if(!maskFieldName.empty())
  {
    mc.setMaskValue(maskVal);
  }
  mc.setFunctionField(fieldName);
  mc.computeIsocontour(contourVal);

  conduit::Node contourBpExec;
  mc.populateContourMeshBlueprint(contourBpExec);
  conduit::Node contourBp;
  copyBlueprintToHost(contourBp, contourBpExec);
  // Validate both output surfaces: the richer welded Blueprint mesh for topology/manifoldness,
  // and the legacy fixed-stride arrays for compatibility with existing MarchingCubes callers.
  ASSERT_TRUE(conduit::blueprint::mesh::is_multi_domain(contourBp));
  ASSERT_EQ(conduit::blueprint::mesh::number_of_domains(contourBp), 1);
  const conduit::Node& contourDom = contourBp.child(0);
  ASSERT_TRUE(contourDom.has_path("state/domain_id"));
  EXPECT_EQ(contourDom["state/domain_id"].to_int32(), 0);
  ASSERT_TRUE(contourDom.has_path("topologies"));
  ASSERT_EQ(contourDom["topologies"].number_of_children(), 1);
  const conduit::Node& contourTopo = contourDom["topologies"].child(0);
  ASSERT_TRUE(contourTopo.has_path("elements/connectivity"));
  ASSERT_TRUE(contourTopo.has_path("elements/sizes"));
  ASSERT_TRUE(contourTopo.has_path("elements/offsets"));
  ASSERT_TRUE(contourDom.has_path("fields/originalElements/values"));

  const axom::Array<double, 2> coordsHost(mc.getContourNodeCoords(), hostAllocatorID());
  const axom::Array<axom::IndexType, 2> cornersHost(mc.getContourFacetCorners(), hostAllocatorID());
  const axom::Array<axom::IndexType> parentsHost(mc.getContourFacetParents(), hostAllocatorID());
  const auto coords = coordsHost.view();
  const auto corners = cornersHost.view();
  const auto parents = parentsHost.view();
  const axom::IndexType nFacets = mc.getContourCellCount();

  ASSERT_GT(nFacets, 0) << "Expected a non-empty contour.";
  // Legacy invariant: node count == facetCount * DIM.
  EXPECT_EQ(mc.getContourNodeCount(), nFacets * 3);

  // O1: on-surface value.
  double maxValErr = 0.0;
  for(axom::IndexType r = 0; r < mc.getContourNodeCount(); ++r)
  {
    const double v = f(coords(r, 0), coords(r, 1), coords(r, 2));
    maxValErr = std::max(maxValErr, std::abs(v - contourVal));
  }
  EXPECT_LT(maxValErr, analyticSurfaceTol) << "O1: a facet node is off the isosurface.";

  // O2: parent id range + facet centroid within parent cell bounds.
  const conduit::index_t nCells =
    conduit::blueprint::mesh::topology::length(mesh["topologies/mesh"]);
  const conduit::Node* maskValues = nullptr;
  if(!maskFieldName.empty())
  {
    ASSERT_TRUE(mesh.has_path("fields/" + maskFieldName + "/values"));
    maskValues = &mesh.fetch_existing("fields/" + maskFieldName + "/values");
  }
  for(axom::IndexType ff = 0; ff < nFacets; ++ff)
  {
    EXPECT_GE(parents[ff], 0);
    EXPECT_LT(parents[ff], static_cast<axom::IndexType>(nCells)) << "O2: parent id out of range.";
    if(!maskFieldName.empty() && parents[ff] >= 0 && parents[ff] < nCells)
    {
      EXPECT_EQ(maskValues->as_int32_accessor()[parents[ff]], maskVal)
        << "masked extraction emitted a facet from an unselected parent zone.";
    }
    if(parents[ff] >= 0 && parents[ff] < nCells)
    {
      BoundingBox3D parentBounds;
      ASSERT_TRUE(parentCellBounds3D(mesh, parents[ff], parentBounds))
        << "O2: could not compute parent cell bounds.";
      const Point3D centroid = contourFacetCentroid3D(coords, corners, ff);
      parentBounds.expand(1.0e-8);
      EXPECT_TRUE(parentBounds.contains(centroid))
        << "O2: facet centroid is outside the reported parent cell bounds.";
    }
  }

  // O3: edge manifoldness.
  const auto em = checkBlueprintEdgeManifold3D(contourDom);
  EXPECT_EQ(em.edgesUsed3PlusTimes, 0)
    << "O3: a contour edge is shared by 3+ facets (non-manifold).";
  EXPECT_LE(em.maxMultiplicity, 2) << "O3: max edge multiplicity exceeds 2.";
  if(expectClosedInterior)
  {
    // A closed surface strictly interior to the domain has no boundary edges.
    EXPECT_EQ(em.boundaryEdges, 0)
      << "O3: closed interior surface unexpectedly has boundary (once-used) edges.";
  }

  conduit::Node relinquishedBpExec;
  mc.relinquishContourDataBlueprint(relinquishedBpExec);
  conduit::Node relinquishedBp;
  copyBlueprintToHost(relinquishedBp, relinquishedBpExec);
  ASSERT_TRUE(conduit::blueprint::mesh::is_multi_domain(relinquishedBp));
  ASSERT_EQ(conduit::blueprint::mesh::number_of_domains(relinquishedBp), 1);
  EXPECT_EQ(mc.getContourCellCount(), 0);
}

//---------------------------------------------------------------------------
// Tests: structured and unstructured, planar and round, per policy.
//---------------------------------------------------------------------------

void test_structured_planar(RuntimePolicy policy)
{
  conduit::Node mesh;
  PlanarField f {{0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}};  // horizontal plane z=0.5
  buildStructured3D(mesh, 8, f, "fcn");
  // Planar contour clips the domain -> open surface (boundary edges expected).
  runAndVerify3D(mesh, f, 0.0, policy, "fcn", /*expectClosedInterior=*/false, 1.0e-6);
}

void test_structured_round(RuntimePolicy policy)
{
  conduit::Node mesh;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};  // sphere fully inside [0,1]^3
  buildStructured3D(mesh, 16, f, "fcn");
  // Sphere interior to the domain -> closed surface (no boundary edges).
  // The contour is exact for the linearly interpolated nodal field, so the
  // analytic signed-distance residual is O(h^2), not roundoff.
  runAndVerify3D(mesh, f, 0.0, policy, "fcn", /*expectClosedInterior=*/true, 5.0e-3);
}

void test_structured_planar_mask(RuntimePolicy policy)
{
  conduit::Node mesh;
  PlanarField f {{0.5, 0.5, 0.30}, {0.0, 0.0, 1.0}};  // horizontal plane z=0.30
  buildStructured3D(mesh, 8, f, "fcn");

  // Select only the lower k-slab. Since z=0.30 lies in that selected half,
  // the contour should be non-empty, and runAndVerify3D checks every reported
  // parent cell has the selected mask value.
  addStructuredMask3D(mesh, 8, "mask", /*selectedValue=*/7, /*rejectedValue=*/3);
  runAndVerify3D(mesh,
                 f,
                 0.0,
                 policy,
                 "fcn",
                 /*expectClosedInterior=*/false,
                 1.0e-6,
                 "mask",
                 7);
}

// Unstructured hex: build structured, then convert in-place via the Axom mesh helper,
// then run the same verification.  Uses a sidre Group because the converter operates on sidre.
void test_unstructured_hex_round(RuntimePolicy policy)
{
  axom::sidre::DataStore ds;
  axom::sidre::Group* meshGrp = ds.getRoot()->createGroup("mesh");

  // Build a structured mesh into a conduit node, import into sidre.
  conduit::Node structured;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  buildStructured3D(structured, 16, f, "fcn");
  meshGrp->importConduitTree(structured);

  // Convert structured -> unstructured single-shape hex in place.
  // Keep the test fixture mesh host-readable for the CPU-side oracle below.
  // runAndVerify3D copies the finished Blueprint mesh into policy-compatible
  // memory before invoking MarchingCubes.
  axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_3d(meshGrp,
                                                                              "mesh",
                                                                              RuntimePolicy::seq);

  // Re-export to a conduit node for MarchingCubes::setMesh.
  conduit::Node unstructured;
  meshGrp->createNativeLayout(unstructured);

  ASSERT_EQ(unstructured["topologies/mesh/type"].as_string(), std::string("unstructured"));
  runAndVerify3D(unstructured, f, 0.0, policy, "fcn", /*expectClosedInterior=*/true, 5.0e-3);
}

// Warped/curvilinear unstructured hex (the case Phase 2.3 called out but the
// initial test set omitted).  Exercises the non-box-hex geometric path: the
// extractor still classifies topology from corner signs, and edge crossings are
// interpolated along physical edges.  The interior sphere remains closed.
void test_unstructured_hex_round_warped(RuntimePolicy policy)
{
  axom::sidre::DataStore ds;
  axom::sidre::Group* meshGrp = ds.getRoot()->createGroup("mesh");

  conduit::Node structured;
  RoundField f {{0.5, 0.5, 0.5}, 0.22};  // slightly smaller r: stays inside after warp
  buildStructured3D(structured, 16, f, "fcn", SinusoidalWarp {0.015});
  meshGrp->importConduitTree(structured);

  axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_3d(meshGrp,
                                                                              "mesh",
                                                                              RuntimePolicy::seq);

  conduit::Node unstructured;
  meshGrp->createNativeLayout(unstructured);

  ASSERT_EQ(unstructured["topologies/mesh/type"].as_string(), std::string("unstructured"));
  // Looser analytic tolerance: on a warped mesh the piecewise-linear contour's
  // signed-distance residual grows with cell distortion, so this checks
  // validity/closedness, not high-accuracy reconstruction.
  runAndVerify3D(unstructured, f, 0.0, policy, "fcn", /*expectClosedInterior=*/true, 2.0e-2);
}

// Phase 6 seam: selecting the `robust` policy must, today, produce a valid
// contour identical to `standard` (robust currently aliases standard).
// This will need to be adjusted once the robust MC backend is added.
void test_robustness_seam(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};

  // Robust currently aliases standard.
  // A future robust intersector can update this expectation.
  auto facetCountFor = [&](quest::MarchingCubesRobustnessPolicy rp) {
    conduit::Node mesh;
    buildStructured3D(mesh, 16, f, "fcn");
    conduit::Node mdMesh;
    mdMesh.append().set(mesh);
    const int allocatorID = axom::policyToDefaultAllocatorID(policy);
    conduit::Node execMdMesh;
    copyBlueprintToPolicy(execMdMesh, mdMesh, policy, allocatorID);
    quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
    mc.setUseBumpBackend(true);
    mc.setRobustnessPolicy(rp);
    mc.setMesh(execMdMesh, "mesh");
    mc.setFunctionField("fcn");
    mc.computeIsocontour(0.0);
    return mc.getContourCellCount();
  };

  const auto stdCount = facetCountFor(quest::MarchingCubesRobustnessPolicy::standard);
  const auto robustCount = facetCountFor(quest::MarchingCubesRobustnessPolicy::robust);
  EXPECT_GT(stdCount, 0);
  EXPECT_EQ(stdCount, robustCount)
    << "robust policy is expected to alias standard until a robust intersector exists.";
}

//---------------------------------------------------------------------------
// GTest registration (sequential always; others when compiled in).
//---------------------------------------------------------------------------

TEST(quest_marching_cubes_bump, structured_planar_seq)
{
  test_structured_planar(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, structured_round_seq) { test_structured_round(RuntimePolicy::seq); }
TEST(quest_marching_cubes_bump, structured_planar_mask_seq)
{
  test_structured_planar_mask(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_seq)
{
  test_unstructured_hex_round(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_warped_seq)
{
  test_unstructured_hex_round_warped(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, robustness_seam_nfc_seq)
{
  test_robustness_seam(RuntimePolicy::seq);
}

  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP) && !defined(_WIN32)
TEST(quest_marching_cubes_bump, structured_round_omp) { test_structured_round(RuntimePolicy::omp); }
TEST(quest_marching_cubes_bump, structured_planar_mask_omp)
{
  test_structured_planar_mask(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_omp)
{
  test_unstructured_hex_round(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_warped_omp)
{
  test_unstructured_hex_round_warped(RuntimePolicy::omp);
}
  #endif

  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(quest_marching_cubes_bump, structured_round_cuda)
{
  test_structured_round(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_bump, structured_planar_mask_cuda)
{
  test_structured_planar_mask(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_cuda)
{
  test_unstructured_hex_round(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_warped_cuda)
{
  test_unstructured_hex_round_warped(RuntimePolicy::cuda);
}
  #endif

  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(quest_marching_cubes_bump, structured_round_hip) { test_structured_round(RuntimePolicy::hip); }
TEST(quest_marching_cubes_bump, structured_planar_mask_hip)
{
  test_structured_planar_mask(RuntimePolicy::hip);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_hip)
{
  test_unstructured_hex_round(RuntimePolicy::hip);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_warped_hip)
{
  test_unstructured_hex_round_warped(RuntimePolicy::hip);
}
  #endif

// Self-test of the O3 edge-manifold helper (independent of MarchingCubes).
TEST(quest_marching_cubes_bump, edge_manifold_helper_selftest)
{
  // Two triangles sharing edge (0,0,0)-(1,0,0): a manifold pair.
  axom::Array<double, 2> coords(axom::ArrayOptions::Uninitialized(), 6, 3);
  // tri 0: (0,0,0),(1,0,0),(0,1,0)
  coords(0, 0) = 0;
  coords(0, 1) = 0;
  coords(0, 2) = 0;
  coords(1, 0) = 1;
  coords(1, 1) = 0;
  coords(1, 2) = 0;
  coords(2, 0) = 0;
  coords(2, 1) = 1;
  coords(2, 2) = 0;
  // tri 1: (0,0,0),(1,0,0),(0,-1,0)  -> shares edge (0,0,0)-(1,0,0)
  coords(3, 0) = 0;
  coords(3, 1) = 0;
  coords(3, 2) = 0;
  coords(4, 0) = 1;
  coords(4, 1) = 0;
  coords(4, 2) = 0;
  coords(5, 0) = 0;
  coords(5, 1) = -1;
  coords(5, 2) = 0;

  axom::Array<axom::IndexType, 2> corners(axom::ArrayOptions::Uninitialized(), 2, 3);
  corners(0, 0) = 0;
  corners(0, 1) = 1;
  corners(0, 2) = 2;
  corners(1, 0) = 3;
  corners(1, 1) = 4;
  corners(1, 2) = 5;

  const auto em = checkEdgeManifold3D(coords.view(), corners.view(), 1.0e-9);
  EXPECT_EQ(em.maxMultiplicity, 2);  // shared edge used twice
  EXPECT_EQ(em.interiorEdges, 1);    // exactly one shared edge
  EXPECT_EQ(em.boundaryEdges, 4);    // the other four edges used once
  EXPECT_EQ(em.edgesUsed3PlusTimes, 0);
}

}  // namespace

#endif  // AXOM_USE_CONDUIT && AXOM_USE_BUMP

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  return result;
}
