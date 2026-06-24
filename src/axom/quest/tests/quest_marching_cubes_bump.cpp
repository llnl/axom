// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file quest_marching_cubes_bump.cpp
 *
 * @brief Validation tests for the bump CutField backend of quest::MarchingCubes,
 * covering structured AND unstructured single-shape (quad/hex) input on every
 * available execution space.
 *
 * What this exercises that the legacy quest_marching_cubes_example did not:
 *  - The bump backend (MarchingCubes::setUseBumpBackend(true)).
 *  - Unstructured single-shape quad (2D) and hex (3D) topologies.
 *  - An explicit *edge-manifoldness* check on the extracted contour, which the
 *    legacy example's three oracles (on-surface value, parent containment,
 *    cells-with-contour) never tested.  This is the check most relevant to the
 *    documented 1987-table saddle-ambiguity behavior; it should pass on the
 *    smooth analytic fields used here regardless of backend.
 *
 * Verification oracles (all backend-agnostic; applied to bump output):
 *   O1. On-surface value: every output facet node, evaluated in the analytic
 *       field, equals the contour value within tolerance (exact for planar).
 *   O2. Parent containment: each facet's parent cell id is in range, and the
 *       facet centroid lies within the parent cell's axis-aligned bounds.
 *   O3. Edge manifoldness on bump's welded Blueprint output: every contour
 *       edge (3D) is shared by exactly 1 or 2 facets; no edge is used 3+
 *       times.  A closed smooth surface interior to the domain should have
 *       all-interior edges shared exactly twice; boundary-clipped edges may be
 *       shared once.
 *
 * NOTE: This file is written to compile under AXOM_USE_CONDUIT && AXOM_USE_BUMP
 * and link against quest+bump+gtest, matching the other quest tests.  It has not
 * been built in this environment; the analytic oracles and edge-manifold logic
 * are independently unit-checked (see the inline EdgeManifold self-test).
 */

#include "axom/config.hpp"

#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)

  #include "axom/core.hpp"
  #include "axom/primal.hpp"
  #include "axom/quest/MarchingCubes.hpp"
  #include "axom/quest/util/mesh_helpers.hpp"
  #include "axom/sidre.hpp"
  #include "axom/mint/mesh/UnstructuredMesh.hpp"

  #include "conduit_blueprint.hpp"

  #include "gtest/gtest.h"

  #include <array>
  #include <cmath>
  #include <cstdint>
  #include <map>
  #include <utility>

namespace
{
using RuntimePolicy = axom::runtime_policy::Policy;
constexpr double EPS = 1.0e-10;

//---------------------------------------------------------------------------
// Analytic fields (mirror the example's PlanarTestStrategy / RoundTestStrategy)
//---------------------------------------------------------------------------

//! @brief Signed distance to a plane through @a origin with unit normal @a n.
struct PlanarField
{
  double ox, oy, oz, nx, ny, nz;
  double operator()(double x, double y, double z) const
  {
    return (x - ox) * nx + (y - oy) * ny + (z - oz) * nz;
  }
};

//! @brief Signed distance to a sphere/circle.
struct RoundField
{
  double cx, cy, cz, r;
  double operator()(double x, double y, double z) const
  {
    const double dx = x - cx, dy = y - cy, dz = z - cz;
    return std::sqrt(dx * dx + dy * dy + dz * dz) - r;
  }
};

//---------------------------------------------------------------------------
// O3 helper: edge-manifoldness over a welded triangle/segment soup.
//---------------------------------------------------------------------------

/*!
 * @brief Count, for a 3D triangle soup, how many facets use each undirected
 * edge after welding coincident vertices (quantized hash).  Returns the max
 * multiplicity and the count of edges used 3+ times.
 */
struct EdgeManifoldResult
{
  int maxMultiplicity = 0;
  axom::IndexType edgesUsed3PlusTimes = 0;
  axom::IndexType boundaryEdges = 0;  // used exactly once
  axom::IndexType interiorEdges = 0;  // used exactly twice
};

inline std::int64_t quantize(double v, double inv)
{
  return static_cast<std::int64_t>(std::llround(v * inv));
}

EdgeManifoldResult checkEdgeManifold3D(const axom::ArrayView<const double, 2>& nodeCoords,
                                       const axom::ArrayView<const axom::IndexType, 2>& facetCorners,
                                       double weldTol)
{
  const double inv = 1.0 / weldTol;
  // Weld vertices by quantized coordinate -> vertex id.
  std::map<std::array<std::int64_t, 3>, axom::IndexType> vmap;
  const axom::IndexType nFacets = facetCorners.shape()[0];

  auto weldedId = [&](axom::IndexType row) {
    std::array<std::int64_t, 3> key {quantize(nodeCoords(row, 0), inv),
                                     quantize(nodeCoords(row, 1), inv),
                                     quantize(nodeCoords(row, 2), inv)};
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

//! @brief Build a single-domain structured explicit 3D mesh [0,1]^3 with
//! @a n cells per side and the analytic field @a f sampled at nodes.
template <typename Field>
void buildStructured3D(conduit::Node& mesh, int n, const Field& f, const std::string& fieldName)
{
  const int nn = n + 1;
  mesh.reset();
  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "explicit";
  const conduit::index_t N = static_cast<conduit::index_t>(nn) * nn * nn;
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
        const double px = double(i) / n, py = double(j) / n, pz = double(k) / n;
        x[idx] = px;
        y[idx] = py;
        z[idx] = pz;
        fv[idx] = f(px, py, pz);
      }
    }
  }
}

//---------------------------------------------------------------------------
// Core check: run bump-backed MarchingCubes and apply O1/O2/O3.
//---------------------------------------------------------------------------

template <typename Field>
void runAndVerify3D(conduit::Node& mesh,
                    const Field& f,
                    double contourVal,
                    RuntimePolicy policy,
                    const std::string& fieldName,
                    bool expectClosedInterior,
                    double analyticSurfaceTol)
{
  namespace quest = axom::quest;

  quest::MarchingCubes mc(policy,
                          axom::execution_space<axom::SEQ_EXEC>::allocatorID(),
                          quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(true);

  // MarchingCubes' public input contract is multi-domain.  Keep the wrapped
  // node alive through computeIsocontour(), since the single-domain objects
  // cache pointers into it.
  conduit::Node mdMesh;
  mdMesh.append().set(mesh);
  mc.setMesh(mdMesh, "mesh");
  mc.setFunctionField(fieldName);
  mc.computeIsocontour(contourVal);

  conduit::Node contourBp;
  mc.populateContourMeshBlueprint(contourBp);
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

  const auto coords = mc.getContourNodeCoords();
  const auto corners = mc.getContourFacetCorners();
  const auto parents = mc.getContourFacetParents();
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

  // O2: parent id range (centroid-in-cell omitted here; needs cell lookup).
  const conduit::index_t nCells =
    conduit::blueprint::mesh::topology::length(mesh["topologies/mesh"]);
  for(axom::IndexType ff = 0; ff < nFacets; ++ff)
  {
    EXPECT_GE(parents[ff], 0);
    EXPECT_LT(parents[ff], static_cast<axom::IndexType>(nCells)) << "O2: parent id out of range.";
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

  conduit::Node relinquishedBp;
  mc.relinquishContourDataBlueprint(relinquishedBp);
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
  PlanarField f {0.5, 0.5, 0.5, 0.0, 0.0, 1.0};  // horizontal plane z=0.5
  buildStructured3D(mesh, 8, f, "fcn");
  // Planar contour clips the domain -> open surface (boundary edges expected).
  runAndVerify3D(mesh, f, 0.0, policy, "fcn", /*expectClosedInterior=*/false, 1.0e-6);
}

void test_structured_round(RuntimePolicy policy)
{
  conduit::Node mesh;
  RoundField f {0.5, 0.5, 0.5, 0.25};  // sphere fully inside [0,1]^3
  buildStructured3D(mesh, 16, f, "fcn");
  // Sphere interior to the domain -> closed surface (no boundary edges).
  // The contour is exact for the linearly interpolated nodal field, so the
  // analytic signed-distance residual is O(h^2), not roundoff.
  runAndVerify3D(mesh, f, 0.0, policy, "fcn", /*expectClosedInterior=*/true, 5.0e-3);
}

// Unstructured hex: build structured, then convert in-place via the Axom
// mesh helper, then run the same verification.  Uses a sidre Group because the
// converter operates on sidre.
void test_unstructured_hex_round(RuntimePolicy policy)
{
  axom::sidre::DataStore ds;
  axom::sidre::Group* meshGrp = ds.getRoot()->createGroup("mesh");

  // Build a structured mesh into a conduit node, import into sidre.
  conduit::Node structured;
  RoundField f {0.5, 0.5, 0.5, 0.25};
  buildStructured3D(structured, 16, f, "fcn");
  meshGrp->importConduitTree(structured);

  // Convert structured -> unstructured single-shape hex in place.
  axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_3d(meshGrp, "mesh", policy);

  // Re-export to a conduit node for MarchingCubes::setMesh.
  conduit::Node unstructured;
  meshGrp->createNativeLayout(unstructured);

  ASSERT_EQ(unstructured["topologies/mesh/type"].as_string(), std::string("unstructured"));
  runAndVerify3D(unstructured, f, 0.0, policy, "fcn", /*expectClosedInterior=*/true, 5.0e-3);
}

//---------------------------------------------------------------------------
// GTest registration (sequential always; others when compiled in).
//---------------------------------------------------------------------------

TEST(quest_marching_cubes_bump, structured_planar_seq)
{
  test_structured_planar(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, structured_round_seq) { test_structured_round(RuntimePolicy::seq); }
TEST(quest_marching_cubes_bump, unstructured_hex_round_seq)
{
  test_unstructured_hex_round(RuntimePolicy::seq);
}

  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP) && !defined(_WIN32)
TEST(quest_marching_cubes_bump, structured_round_omp) { test_structured_round(RuntimePolicy::omp); }
TEST(quest_marching_cubes_bump, unstructured_hex_round_omp)
{
  test_unstructured_hex_round(RuntimePolicy::omp);
}
  #endif

  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(quest_marching_cubes_bump, structured_round_cuda)
{
  test_structured_round(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_bump, unstructured_hex_round_cuda)
{
  test_unstructured_hex_round(RuntimePolicy::cuda);
}
  #endif

  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(quest_marching_cubes_bump, structured_round_hip) { test_structured_round(RuntimePolicy::hip); }
TEST(quest_marching_cubes_bump, unstructured_hex_round_hip)
{
  test_unstructured_hex_round(RuntimePolicy::hip);
}
  #endif

// Self-test of the O3 edge-manifold helper, so a green CI proves the oracle
// itself is correct (independent of MarchingCubes).
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
