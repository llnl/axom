// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file quest_marching_cubes_bump.cpp
 *
 * @brief Tests the Bump backend of quest::MarchingCubes.
 *
 * The tests cover structured meshes and single-shape unstructured hex meshes
 * in each enabled execution space. They check analytic residuals, parent-cell IDs,
 * crossing cells, and edge incidence in the welded output.
 */

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "quest_marching_cubes_bump.cpp requires conduit"
#endif
#ifndef AXOM_USE_BUMP
  #error "quest_marching_cubes_bump.cpp requires bump"
#endif
#ifndef AXOM_USE_SIDRE
  #error "quest_marching_cubes_bump.cpp requires sidre"
#endif

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/spin/MortonIndex.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/util/mesh_helpers.hpp"

#include "conduit_blueprint.hpp"

#include "axom/quest/tests/quest_marching_cubes_testing_helpers.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>

namespace
{
namespace mctest = axom::quest::testing::marching_cubes;

using mctest::copyBlueprintToHost;
using mctest::copyBlueprintToPolicy;
using mctest::GyroidField;
using mctest::hostAllocatorID;
using mctest::PlanarField;
using mctest::RoundField;
using mctest::SinusoidalWarp;

using RuntimePolicy = axom::runtime_policy::Policy;
using QuantizedPoint3D = axom::primal::Point<std::int64_t, 3>;

//---------------------------------------------------------------------------
// Edge-manifold checks
//---------------------------------------------------------------------------

/// Summary of edge incidence in a welded surface mesh.
struct EdgeManifoldResult
{
  int maxMultiplicity = 0;
  axom::IndexType edgesUsed3PlusTimes = 0;
  axom::IndexType boundaryEdges = 0;  // used exactly once
  axom::IndexType interiorEdges = 0;  // used exactly twice
};

/// Count triangle incidence after welding coincident coordinates with a quantized hash.
EdgeManifoldResult checkEdgeManifold3D(const axom::ArrayView<const double, 2>& nodeCoords,
                                       const axom::ArrayView<const axom::IndexType, 2>& facetCorners,
                                       double weldTol)
{
  const double inv = 1.0 / weldTol;
  auto quantize = [inv](double v) { return static_cast<std::int64_t>(std::llround(v * inv)); };

  // The legacy MarchingCubes arrays duplicate coordinates per facet.
  // Recover welded vertex ids for the helper self-test.
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

/// Count edge incidence directly from a welded Blueprint contour.
EdgeManifoldResult checkBlueprintEdgeManifold3D(const conduit::Node& contourDom)
{
  // Blueprint output is already welded, so count its connectivity directly.
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
  // structured mesh in k. The masked test below verifies that Bump's
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

template <int DIM>
using Point = axom::primal::Point<double, DIM>;

template <int DIM>
using BoundingBox = axom::primal::BoundingBox<double, DIM>;

template <int DIM>
bool parentCellNodeIds(const conduit::Node& mesh,
                       axom::IndexType cell_index,
                       std::array<axom::IndexType, 1 << DIM>& node_ids)
{
  const conduit::Node& topo = mesh.fetch_existing("topologies/mesh");
  const std::string topo_type = topo.fetch_existing("type").as_string();

  if(topo_type == "structured")
  {
    const axom::IndexType ni = topo.fetch_existing("elements/dims/i").to_index_t();
    const axom::IndexType nj = topo.fetch_existing("elements/dims/j").to_index_t();
    const axom::IndexType nk = DIM == 3 ? topo.fetch_existing("elements/dims/k").to_index_t() : 1;
    if(cell_index < 0 || cell_index >= ni * nj * nk)
    {
      return false;
    }

    const axom::IndexType i = cell_index % ni;
    const axom::IndexType j = (cell_index / ni) % nj;
    const axom::IndexType k = DIM == 3 ? cell_index / (ni * nj) : 0;
    const axom::IndexType nni = ni + 1;
    const axom::IndexType nnj = nj + 1;

    for(int corner = 0; corner < (1 << DIM); ++corner)
    {
      const axom::IndexType ii = i + ((corner & 1) != 0);
      const axom::IndexType jj = j + ((corner & 2) != 0);
      const axom::IndexType kk = k + ((corner & 4) != 0);
      node_ids[corner] = ii + jj * nni + kk * nni * nnj;
    }
    return true;
  }

  if(topo_type == "unstructured")
  {
    const std::string expected_shape = DIM == 3 ? "hex" : "quad";
    if(topo.fetch_existing("elements/shape").as_string() != expected_shape)
    {
      return false;
    }

    const auto conn = topo.fetch_existing("elements/connectivity").as_index_t_accessor();
    constexpr axom::IndexType nodes_per_cell = 1 << DIM;
    const axom::IndexType first = cell_index * nodes_per_cell;
    if(cell_index < 0 || first + nodes_per_cell > conn.number_of_elements())
    {
      return false;
    }
    for(int corner = 0; corner < nodes_per_cell; ++corner)
    {
      node_ids[corner] = static_cast<axom::IndexType>(conn[first + corner]);
    }
    return true;
  }

  return false;
}

template <int DIM>
bool parentCellBounds(const conduit::Node& mesh, axom::IndexType cell_index, BoundingBox<DIM>& bounds)
{
  std::array<axom::IndexType, 1 << DIM> node_ids;
  if(!parentCellNodeIds<DIM>(mesh, cell_index, node_ids))
  {
    return false;
  }

  const conduit::Node& topo = mesh.fetch_existing("topologies/mesh");
  const std::string coordset_name = topo.fetch_existing("coordset").as_string();
  const conduit::Node& values = mesh.fetch_existing("coordsets/" + coordset_name + "/values");
  const auto x = values.fetch_existing("x").as_double_accessor();
  const auto y = values.fetch_existing("y").as_double_accessor();
  const auto z =
    DIM == 3 ? values.fetch_existing("z").as_double_accessor() : conduit::double_accessor {};

  for(const auto node_id : node_ids)
  {
    Point<DIM> point;
    point[0] = x[node_id];
    point[1] = y[node_id];
    if(DIM == 3)
    {
      point[2] = z[node_id];
    }
    bounds.addPoint(point);
  }
  return bounds.isValid();
}

template <int DIM>
void verifyExpectedCrossingCells(const conduit::Node& mesh,
                                 const std::string& field_name,
                                 double contour_value,
                                 const std::set<axom::IndexType>& actual_cells,
                                 const std::string& mask_field_name = {},
                                 int mask_value = 1)
{
  const auto field = mesh.fetch_existing("fields/" + field_name + "/values").as_double_accessor();
  const conduit::index_t num_cells =
    conduit::blueprint::mesh::topology::length(mesh.fetch_existing("topologies/mesh"));
  const conduit::Node* mask = mask_field_name.empty()
    ? nullptr
    : &mesh.fetch_existing("fields/" + mask_field_name + "/values");

  for(axom::IndexType cell = 0; cell < num_cells; ++cell)
  {
    std::array<axom::IndexType, 1 << DIM> node_ids;
    ASSERT_TRUE(parentCellNodeIds<DIM>(mesh, cell, node_ids));

    bool below = false;
    bool at_or_above = false;
    for(const auto node_id : node_ids)
    {
      const float value = static_cast<float>(field[node_id]);
      below = below || value < static_cast<float>(contour_value);
      at_or_above = at_or_above || value >= static_cast<float>(contour_value);
    }

    const bool selected = mask == nullptr || mask->as_int32_accessor()[cell] == mask_value;
    const bool expected = selected && below && at_or_above;
    EXPECT_EQ(actual_cells.count(cell) != 0, expected)
      << "parent cell " << cell << " disagrees with its corner classification";
  }
}

template <int DIM, typename Field>
void verifyAnalyticFacet(const conduit::Node& mesh,
                         const Field& field,
                         double contour_value,
                         double surface_tolerance,
                         double geometry_tolerance,
                         const axom::ArrayView<const double, 2>& node_coords,
                         const axom::ArrayView<const axom::IndexType, 2>& facet_corners,
                         const axom::ArrayView<const axom::IndexType>& facet_parents,
                         axom::IndexType facet)
{
  const axom::IndexType num_cells = static_cast<axom::IndexType>(
    conduit::blueprint::mesh::topology::length(mesh.fetch_existing("topologies/mesh")));
  const axom::IndexType parent = facet_parents[facet];
  ASSERT_GE(parent, 0);
  ASSERT_LT(parent, num_cells);

  BoundingBox<DIM> bounds;
  ASSERT_TRUE(parentCellBounds<DIM>(mesh, parent, bounds));
  BoundingBox<DIM> expanded(bounds);
  expanded.expand(geometry_tolerance);
  const bool structured = mesh.fetch_existing("topologies/mesh/type").as_string() == "structured";

  for(int corner = 0; corner < DIM; ++corner)
  {
    const axom::IndexType node = facet_corners(facet, corner);
    ASSERT_GE(node, 0);
    ASSERT_LT(node, node_coords.shape()[0]);

    Point<DIM> point;
    point[0] = node_coords(node, 0);
    point[1] = node_coords(node, 1);
    if(DIM == 3)
    {
      point[2] = node_coords(node, 2);
    }

    const double z = DIM == 3 ? point[2] : 0.0;
    EXPECT_NEAR(field(point[0], point[1], z), contour_value, surface_tolerance)
      << "facet " << facet << ", corner " << corner << " is off the analytic contour";
    EXPECT_TRUE(expanded.contains(point))
      << "facet " << facet << ", corner " << corner << " lies outside parent cell " << parent;

    if(structured)
    {
      double distance_to_face = std::numeric_limits<double>::max();
      for(int dim = 0; dim < DIM; ++dim)
      {
        distance_to_face = std::min(distance_to_face, std::abs(point[dim] - bounds.getMin()[dim]));
        distance_to_face = std::min(distance_to_face, std::abs(point[dim] - bounds.getMax()[dim]));
      }
      EXPECT_LE(distance_to_face, geometry_tolerance)
        << "facet " << facet << ", corner " << corner
        << " lies in the interior of its structured parent cell";
    }
  }
}

template <int DIM, typename Field>
void verifyAnalyticFacets(const conduit::Node& mesh,
                          const Field& field,
                          const std::string& field_name,
                          double contour_value,
                          double surface_tolerance,
                          double geometry_tolerance,
                          const axom::ArrayView<const double, 2>& node_coords,
                          const axom::ArrayView<const axom::IndexType, 2>& facet_corners,
                          const axom::ArrayView<const axom::IndexType>& facet_parents,
                          axom::IndexType facet_begin,
                          axom::IndexType facet_end,
                          const std::string& mask_field_name = {},
                          int mask_value = 1)
{
  std::set<axom::IndexType> actual_cells;
  for(axom::IndexType facet = facet_begin; facet < facet_end; ++facet)
  {
    verifyAnalyticFacet<DIM>(mesh,
                             field,
                             contour_value,
                             surface_tolerance,
                             geometry_tolerance,
                             node_coords,
                             facet_corners,
                             facet_parents,
                             facet);
    actual_cells.insert(facet_parents[facet]);
  }

  verifyExpectedCrossingCells<DIM>(mesh,
                                   field_name,
                                   contour_value,
                                   actual_cells,
                                   mask_field_name,
                                   mask_value);
}

//---------------------------------------------------------------------------
// Run MarchingCubes and validate its Bump output.
//---------------------------------------------------------------------------

template <int DIM, typename Field>
void runAndVerify(conduit::Node& mesh,
                  const Field& field,
                  double contour_value,
                  RuntimePolicy policy,
                  const std::string& field_name,
                  bool expect_closed_interior,
                  double surface_tolerance,
                  const std::string& mask_field_name = {},
                  int mask_value = 1,
                  axom::quest::MarchingCubesRobustnessPolicy robustness =
                    axom::quest::MarchingCubesRobustnessPolicy::standard)
{
  namespace quest = axom::quest;

  const int allocatorID = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(true);
  mc.setRobustnessPolicy(robustness);

  conduit::Node execMesh;
  copyBlueprintToPolicy(execMesh, mesh, policy, allocatorID);
  mc.setMesh(execMesh, "mesh", mask_field_name);
  if(!mask_field_name.empty())
  {
    mc.setMaskValue(mask_value);
  }
  mc.setFunctionField(field_name);
  mc.computeIsocontour(contour_value);

  conduit::Node contourBpExec;
  mc.populateContourMeshBlueprint(contourBpExec);
  conduit::Node contourBp;
  copyBlueprintToHost(contourBp, contourBpExec);
  // Check the welded Blueprint output and the fixed-stride compatibility arrays.
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

  conduit::Node triContourBpExec;
  mc.populateContourMeshBlueprint(triContourBpExec, true);
  conduit::Node triContourBp;
  copyBlueprintToHost(triContourBp, triContourBpExec);
  ASSERT_TRUE(conduit::blueprint::mesh::is_multi_domain(triContourBp));
  ASSERT_EQ(conduit::blueprint::mesh::number_of_domains(triContourBp), 1);
  const conduit::Node& triContourDom = triContourBp.child(0);
  const conduit::Node& triContourTopo = triContourDom["topologies"].child(0);
  const conduit::Node& triElems = triContourTopo.fetch_existing("elements");
  const auto triSizes = triElems.fetch_existing("sizes").as_index_t_accessor();
  const auto triConn = triElems.fetch_existing("connectivity").as_index_t_accessor();
  ASSERT_EQ(triConn.number_of_elements(), triSizes.number_of_elements() * DIM);
  for(conduit::index_t z = 0; z < triSizes.number_of_elements(); ++z)
  {
    EXPECT_EQ(triSizes[z], DIM);
  }
  EXPECT_EQ(triContourDom["fields/originalElements/values"].dtype().number_of_elements(),
            triSizes.number_of_elements());
  const std::string contourCoordsetName = contourTopo.fetch_existing("coordset").as_string();
  const std::string triCoordsetName = triContourTopo.fetch_existing("coordset").as_string();
  EXPECT_EQ(contourDom.fetch_existing("coordsets/" + contourCoordsetName + "/values/x")
              .dtype()
              .number_of_elements(),
            triContourDom.fetch_existing("coordsets/" + triCoordsetName + "/values/x")
              .dtype()
              .number_of_elements());

  const axom::Array<double, 2> coordsHost(mc.getContourNodeCoords(), hostAllocatorID());
  const axom::Array<axom::IndexType, 2> cornersHost(mc.getContourFacetCorners(), hostAllocatorID());
  const axom::Array<axom::IndexType> parentsHost(mc.getContourFacetParents(), hostAllocatorID());
  const auto coords = coordsHost.view();
  const auto corners = cornersHost.view();
  const auto parents = parentsHost.view();
  const axom::IndexType nFacets = mc.getContourCellCount();

  ASSERT_GT(nFacets, 0) << "Expected a non-empty contour.";
  EXPECT_LE(mc.getContourNodeCount(), nFacets * DIM)
    << "Bump output did not reuse contour vertices.";
  for(axom::IndexType fIdx = 0; fIdx < nFacets; ++fIdx)
  {
    for(int c = 0; c < DIM; ++c)
    {
      EXPECT_GE(corners(fIdx, c), 0);
      EXPECT_LT(corners(fIdx, c), mc.getContourNodeCount());
    }
  }

  verifyAnalyticFacets<DIM>(mesh,
                            field,
                            field_name,
                            contour_value,
                            surface_tolerance,
                            1.0e-5,
                            coords,
                            corners,
                            parents,
                            0,
                            nFacets,
                            mask_field_name,
                            mask_value);

  // Check edge incidence on the welded Blueprint output.
  if(DIM == 3)
  {
    const auto em = checkBlueprintEdgeManifold3D(contourDom);
    EXPECT_EQ(em.edgesUsed3PlusTimes, 0) << "A contour edge is shared by at least three facets.";
    EXPECT_LE(em.maxMultiplicity, 2) << "Contour edge multiplicity exceeds two.";
    if(expect_closed_interior)
    {
      EXPECT_EQ(em.boundaryEdges, 0) << "Closed contour has boundary edges.";
    }
  }

  conduit::Node relinquishedBpExec;
  mc.relinquishContourDataBlueprint(relinquishedBpExec);
  conduit::Node relinquishedBp;
  copyBlueprintToHost(relinquishedBp, relinquishedBpExec);
  ASSERT_TRUE(conduit::blueprint::mesh::is_multi_domain(relinquishedBp));
  ASSERT_EQ(conduit::blueprint::mesh::number_of_domains(relinquishedBp), 1);
  EXPECT_EQ(mc.getContourCellCount(), 0);
  EXPECT_EQ(mc.getContourNodeCount(), 0);
}

//---------------------------------------------------------------------------
// Structured and unstructured tests
//---------------------------------------------------------------------------

void test_structured_planar(RuntimePolicy policy)
{
  conduit::Node mesh;
  PlanarField f {{0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}};  // horizontal plane z=0.5
  mctest::buildStructured<3>(mesh, 8, f, "fcn");
  // The plane clips the domain, so the contour is open.
  runAndVerify<3>(mesh, f, 0.0, policy, "fcn", /*expect_closed_interior=*/false, 1.0e-6);
}

void test_structured_round(RuntimePolicy policy)
{
  conduit::Node mesh;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};  // sphere fully inside [0,1]^3
  mctest::buildStructured<3>(mesh, 16, f, "fcn");
  // The sphere is inside the domain, so the contour is closed.
  // The contour is exact for the linearly interpolated nodal field, so the
  // analytic signed-distance residual is O(h^2), not roundoff.
  runAndVerify<3>(mesh, f, 0.0, policy, "fcn", /*expect_closed_interior=*/true, 5.0e-3);
}

void test_structured_planar_mask(RuntimePolicy policy)
{
  conduit::Node mesh;
  PlanarField f {{0.5, 0.5, 0.30}, {0.0, 0.0, 1.0}};  // horizontal plane z=0.30
  mctest::buildStructured<3>(mesh, 8, f, "fcn");

  // Select only the lower k-slab. Since z=0.30 lies in that selected half,
  // the contour should be non-empty, and runAndVerify checks every reported
  // parent cell has the selected mask value.
  addStructuredMask3D(mesh, 8, "mask", /*selectedValue=*/7, /*rejectedValue=*/3);
  runAndVerify<3>(mesh,
                  f,
                  0.0,
                  policy,
                  "fcn",
                  /*expect_closed_interior=*/false,
                  1.0e-6,
                  "mask",
                  7);
}

// The conversion helper operates on a Sidre group.
void test_unstructured_hex_round(RuntimePolicy policy)
{
  axom::sidre::DataStore ds;
  axom::sidre::Group* meshGrp = ds.getRoot()->createGroup("mesh");

  conduit::Node structured;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  mctest::buildStructured<3>(structured, 16, f, "fcn");
  meshGrp->importConduitTree(structured);

  // Keep the mesh in host memory for the analytic checks in runAndVerify().
  axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_3d(meshGrp,
                                                                              "mesh",
                                                                              RuntimePolicy::seq);

  conduit::Node unstructured;
  meshGrp->createNativeLayout(unstructured);

  ASSERT_EQ(unstructured["topologies/mesh/type"].as_string(), std::string("unstructured"));
  runAndVerify<3>(unstructured,
                  f,
                  0.0,
                  policy,
                  "fcn",
                  /*expect_closed_interior=*/true,
                  5.0e-3);
}

// Exercise physical-edge interpolation on a curvilinear hex mesh.
void test_unstructured_hex_round_warped(RuntimePolicy policy)
{
  axom::sidre::DataStore ds;
  axom::sidre::Group* meshGrp = ds.getRoot()->createGroup("mesh");

  conduit::Node structured;
  RoundField f {{0.5, 0.5, 0.5}, 0.22};  // stays inside the warped mesh
  mctest::buildStructured<3>(structured, 16, f, "fcn", SinusoidalWarp {0.015});
  meshGrp->importConduitTree(structured);

  axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_3d(meshGrp,
                                                                              "mesh",
                                                                              RuntimePolicy::seq);

  conduit::Node unstructured;
  meshGrp->createNativeLayout(unstructured);

  ASSERT_EQ(unstructured["topologies/mesh/type"].as_string(), std::string("unstructured"));
  // Cell distortion increases the piecewise-linear contour's residual.
  runAndVerify<3>(unstructured,
                  f,
                  0.0,
                  policy,
                  "fcn",
                  /*expect_closed_interior=*/true,
                  2.0e-2);
}

// The robust policy currently aliases the standard policy.
void test_robustness_policy(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};

  auto facetCountFor = [&](quest::MarchingCubesRobustnessPolicy rp) {
    conduit::Node mesh;
    mctest::buildStructured<3>(mesh, 16, f, "fcn");
    const int allocatorID = axom::policyToDefaultAllocatorID(policy);
    conduit::Node execMesh;
    copyBlueprintToPolicy(execMesh, mesh, policy, allocatorID);
    quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
    mc.setUseBumpBackend(true);
    mc.setRobustnessPolicy(rp);
    mc.setMesh(execMesh, "mesh");
    mc.setFunctionField("fcn");
    mc.computeIsocontour(0.0);
    return mc.getContourCellCount();
  };

  const auto stdCount = facetCountFor(quest::MarchingCubesRobustnessPolicy::standard);
  const auto robustCount = facetCountFor(quest::MarchingCubesRobustnessPolicy::robust);
  EXPECT_GT(stdCount, 0);
  EXPECT_EQ(stdCount, robustCount) << "Robust and standard policies should currently match.";
}

template <int DIM>
void test_accumulated_analytic_fields(RuntimePolicy policy)
{
  namespace quest = axom::quest;

  constexpr int n = 12;
  PlanarField plane {{0.47, 0.43, 0.39}, {1.0, 0.4, 1.2}};
  RoundField round {{0.5, 0.5, DIM == 3 ? 0.5 : 0.0}, 0.27};
  GyroidField gyroid {{3.0, 3.0, DIM == 3 ? 1.5 : 0.0}};
  const double mesh_spacing = 1.0 / n;
  const double geometry_tolerance = 1.0e-5 * mesh_spacing;
  const double round_tolerance = 0.1 * mesh_spacing;
  const double gyroid_scale_norm =
    std::sqrt(gyroid.scale[0] * gyroid.scale[0] + gyroid.scale[1] * gyroid.scale[1] +
              gyroid.scale[2] * gyroid.scale[2]);
  const double gyroid_tolerance = 0.1 * gyroid_scale_norm * mesh_spacing;

  conduit::Node mesh;
  mctest::buildStructured<DIM>(mesh, n, plane, "plane");
  mctest::addVertexField<DIM>(mesh, round, "round");
  mctest::addVertexField<DIM>(mesh, gyroid, "gyroid");

  const int allocator_id = axom::policyToDefaultAllocatorID(policy);
  conduit::Node exec_mesh;
  copyBlueprintToPolicy(exec_mesh, mesh, policy, allocator_id);

  quest::MarchingCubes mc(policy, allocator_id, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(true);
  mc.setMesh(exec_mesh, "mesh");

  mc.setFunctionField("plane");
  mc.computeIsocontour(0.0);
  const axom::IndexType plane_end = mc.getContourFacetCount();

  mc.setFunctionField("round");
  mc.computeIsocontour(0.0);
  const axom::IndexType round_end = mc.getContourFacetCount();

  mc.setFunctionField("gyroid");
  mc.computeIsocontour(0.0);
  const axom::IndexType gyroid_end = mc.getContourFacetCount();

  ASSERT_GT(plane_end, 0);
  ASSERT_GT(round_end, plane_end);
  ASSERT_GT(gyroid_end, round_end);

  const axom::Array<double, 2> coords_host(mc.getContourNodeCoords(), hostAllocatorID());
  const axom::Array<axom::IndexType, 2> corners_host(mc.getContourFacetCorners(), hostAllocatorID());
  const axom::Array<axom::IndexType> parents_host(mc.getContourFacetParents(), hostAllocatorID());
  const auto coords = coords_host.view();
  const auto corners = corners_host.view();
  const auto parents = parents_host.view();

  verifyAnalyticFacets<DIM>(mesh,
                            plane,
                            "plane",
                            0.0,
                            geometry_tolerance,
                            geometry_tolerance,
                            coords,
                            corners,
                            parents,
                            0,
                            plane_end);
  verifyAnalyticFacets<DIM>(mesh,
                            round,
                            "round",
                            0.0,
                            round_tolerance,
                            geometry_tolerance,
                            coords,
                            corners,
                            parents,
                            plane_end,
                            round_end);
  verifyAnalyticFacets<DIM>(mesh,
                            gyroid,
                            "gyroid",
                            0.0,
                            gyroid_tolerance,
                            geometry_tolerance,
                            coords,
                            corners,
                            parents,
                            round_end,
                            gyroid_end);
}

//---------------------------------------------------------------------------
// Multi-domain accumulation.
// Separate runs provide counts without output offsets.
// A combined run exercises offsets in the shared output arrays.
//---------------------------------------------------------------------------

/// Shift every x-coordinate in \a mesh by \a offset.
void translateMeshX(conduit::Node& mesh, double offset)
{
  conduit::Node& n_x = mesh.fetch_existing("coordsets/coords/values/x");
  auto* x = n_x.as_float64_ptr();
  for(conduit::index_t i = 0; i < n_x.dtype().number_of_elements(); ++i)
  {
    x[i] += offset;
  }
}

template <int DIM>
void test_multidomain_planar(RuntimePolicy policy)
{
  namespace quest = axom::quest;

  // Different plane heights reveal facets stored at the wrong domain offset.
  PlanarField f0 {{0.5, 0.5, DIM == 3 ? 0.5 : 0.0},
                  {0.0, DIM == 2 ? 1.0 : 0.0, DIM == 3 ? 1.0 : 0.0}};
  PlanarField f1 {{0.5, 0.3, DIM == 3 ? 0.3 : 0.0},
                  {0.0, DIM == 2 ? 1.0 : 0.0, DIM == 3 ? 1.0 : 0.0}};

  conduit::Node dom0, dom1;
  mctest::buildStructured<DIM>(dom0, 10, f0, "fcn");
  mctest::buildStructured<DIM>(dom1, 10, f1, "fcn");
  translateMeshX(dom1, 2.0);
  dom0["state/domain_id"] = 7;
  dom1["state/domain_id"] = 19;

  const int allocator_id = axom::policyToDefaultAllocatorID(policy);
  const std::array<const conduit::Node*, 2> domains {{&dom0, &dom1}};
  const std::array<const PlanarField*, 2> fields {{&f0, &f1}};
  const std::array<axom::IndexType, 2> domain_ids {{7, 19}};

  // Get reference counts without multi-domain offsets.
  axom::IndexType separate_facet_count = 0;
  axom::IndexType separate_node_count = 0;
  for(const conduit::Node* domain : domains)
  {
    conduit::Node exec_domain;
    copyBlueprintToPolicy(exec_domain, *domain, policy, allocator_id);
    quest::MarchingCubes mc(policy, allocator_id, quest::MarchingCubesDataParallelism::byPolicy);
    mc.setUseBumpBackend(true);
    mc.setMesh(exec_domain, "mesh");
    mc.setFunctionField("fcn");
    mc.computeIsocontour(0.0);
    separate_facet_count += mc.getContourCellCount();
    separate_node_count += mc.getContourNodeCount();
  }
  ASSERT_GT(separate_facet_count, 0);

  // Run both domains through one MarchingCubes instance.
  conduit::Node mdMesh;
  mdMesh.append().set_external(dom0);
  mdMesh.append().set_external(dom1);
  ASSERT_TRUE(conduit::blueprint::mesh::is_multi_domain(mdMesh));

  conduit::Node exec_mesh;
  copyBlueprintToPolicy(exec_mesh, mdMesh, policy, allocator_id);
  quest::MarchingCubes mc(policy, allocator_id, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setUseBumpBackend(true);
  mc.setMesh(exec_mesh, "mesh");
  mc.setFunctionField("fcn");
  mc.computeIsocontour(0.0);

  EXPECT_EQ(mc.getContourCellCount(), separate_facet_count)
    << "Multi-domain facet count differs from the separate runs.";
  EXPECT_EQ(mc.getContourNodeCount(), separate_node_count)
    << "Multi-domain node count differs from the separate runs.";

  axom::Array<axom::IndexType, 2> facetNodeIds;
  axom::Array<double, 2> facetNodeCoords;
  axom::Array<axom::IndexType, 1> facetParentIds;
  axom::Array<axom::IndexType> facetDomainIds;
  mc.relinquishContourData(facetNodeIds, facetNodeCoords, facetParentIds, facetDomainIds);
  EXPECT_EQ(mc.getContourCellCount(), 0);
  EXPECT_EQ(mc.getContourNodeCount(), 0);

  // Copy policy-allocated output to the host for validation.
  const axom::Array<axom::IndexType, 2> ids(facetNodeIds, hostAllocatorID());
  const axom::Array<double, 2> coords(facetNodeCoords, hostAllocatorID());
  const axom::Array<axom::IndexType> parents(facetParentIds, hostAllocatorID());
  const axom::Array<axom::IndexType> output_domain_ids(facetDomainIds, hostAllocatorID());

  const axom::IndexType num_facets = output_domain_ids.size();
  ASSERT_EQ(ids.shape()[0], num_facets);
  ASSERT_EQ(parents.size(), num_facets);

  std::array<std::set<axom::IndexType>, 2> actual_cells;
  std::array<bool, 2> saw_domain {{false, false}};
  for(axom::IndexType facet = 0; facet < num_facets; ++facet)
  {
    const auto domain_id = output_domain_ids[facet];
    const int domain_index = domain_id == domain_ids[0] ? 0 : domain_id == domain_ids[1] ? 1 : -1;
    ASSERT_GE(domain_index, 0) << "Unexpected parent domain id " << domain_id;
    saw_domain[domain_index] = true;
    actual_cells[domain_index].insert(parents[facet]);
    verifyAnalyticFacet<DIM>(*domains[domain_index],
                             *fields[domain_index],
                             0.0,
                             1.0e-6,
                             1.0e-5,
                             coords.view(),
                             ids.view(),
                             parents.view(),
                             facet);
  }
  EXPECT_TRUE(saw_domain[0] && saw_domain[1]) << "Both domains should contribute facets.";
  for(int domain = 0; domain < 2; ++domain)
  {
    verifyExpectedCrossingCells<DIM>(*domains[domain], "fcn", 0.0, actual_cells[domain]);
  }
}

//---------------------------------------------------------------------------
// gtest registration
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
TEST(quest_marching_cubes_bump, robust_matches_standard_seq)
{
  test_robustness_policy(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, multidomain_planar_2d_seq)
{
  test_multidomain_planar<2>(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, multidomain_planar_3d_seq)
{
  test_multidomain_planar<3>(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_2d_seq)
{
  test_accumulated_analytic_fields<2>(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_3d_seq)
{
  test_accumulated_analytic_fields<3>(RuntimePolicy::seq);
}

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP) && !defined(_WIN32)
TEST(quest_marching_cubes_bump, structured_round_omp) { test_structured_round(RuntimePolicy::omp); }
TEST(quest_marching_cubes_bump, structured_planar_mask_omp)
{
  test_structured_planar_mask(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_bump, multidomain_planar_2d_omp)
{
  test_multidomain_planar<2>(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_bump, multidomain_planar_3d_omp)
{
  test_multidomain_planar<3>(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_2d_omp)
{
  test_accumulated_analytic_fields<2>(RuntimePolicy::omp);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_3d_omp)
{
  test_accumulated_analytic_fields<3>(RuntimePolicy::omp);
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
TEST(quest_marching_cubes_bump, multidomain_planar_2d_cuda)
{
  test_multidomain_planar<2>(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_bump, multidomain_planar_3d_cuda)
{
  test_multidomain_planar<3>(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_2d_cuda)
{
  test_accumulated_analytic_fields<2>(RuntimePolicy::cuda);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_3d_cuda)
{
  test_accumulated_analytic_fields<3>(RuntimePolicy::cuda);
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
TEST(quest_marching_cubes_bump, multidomain_planar_2d_hip)
{
  test_multidomain_planar<2>(RuntimePolicy::hip);
}
TEST(quest_marching_cubes_bump, multidomain_planar_3d_hip)
{
  test_multidomain_planar<3>(RuntimePolicy::hip);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_2d_hip)
{
  test_accumulated_analytic_fields<2>(RuntimePolicy::hip);
}
TEST(quest_marching_cubes_bump, accumulated_analytic_fields_3d_hip)
{
  test_accumulated_analytic_fields<3>(RuntimePolicy::hip);
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

// Test the edge-manifold helper without MarchingCubes.
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

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();
  return result;
}
