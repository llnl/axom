// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_MFEM)

  #include "gtest/gtest.h"

  #include "axom/core.hpp"
  #include "axom/quest/SamplingShaper.hpp"
  #include "axom/quest/detail/shaping/shaping_helpers.hpp"
  #include "axom/bump/utilities/conduit_memory.hpp"
  #include "axom/bump/views/dispatch_coordset.hpp"
  #include "axom/sidre.hpp"

  #include "conduit.hpp"
  #include "conduit_blueprint.hpp"

namespace
{

template <typename T, typename U>
bool compareArrayView(axom::ArrayView<T> lhs, axom::ArrayView<U> rhs)
{
  if(lhs.size() != rhs.size())
  {
    return false;
  }

  for(axom::IndexType i = 0; i < lhs.size(); ++i)
  {
    if(lhs[i] != rhs[i])
    {
      return false;
    }
  }
  return true;
}

void setNodeValues(conduit::Node& node, axom::ArrayView<const double> values)
{
  node.set(conduit::DataType::float64(values.size()));
  auto* data = node.as_float64_ptr();
  for(axom::IndexType i = 0; i < values.size(); ++i)
  {
    data[i] = values[i];
  }
}

void setNodeValues(conduit::Node& node, axom::ArrayView<const conduit::index_t> values)
{
  node.set(conduit::DataType::index_t(values.size()));
  auto* data = node.as_index_t_ptr();
  for(axom::IndexType i = 0; i < values.size(); ++i)
  {
    data[i] = values[i];
  }
}

conduit::Node makeQuadMesh()
{
  conduit::Node mesh;

  mesh["coordsets/coords/type"] = "explicit";
  const axom::Array<double> x {{0., 1., 0., 1.}};
  const axom::Array<double> y {{0., 0., 1., 1.}};
  setNodeValues(mesh["coordsets/coords/values/x"], x.view());
  setNodeValues(mesh["coordsets/coords/values/y"], y.view());

  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "quad";
  const axom::Array<conduit::index_t> connectivity {{0, 1, 3, 2}};
  setNodeValues(mesh["topologies/mesh/elements/connectivity"], connectivity.view());

  return mesh;
}

conduit::Node makeHexMesh()
{
  conduit::Node mesh;

  mesh["coordsets/coords/type"] = "explicit";
  const axom::Array<double> x {{0., 1., 0., 1., 0., 1., 0., 1.}};
  const axom::Array<double> y {{0., 0., 1., 1., 0., 0., 1., 1.}};
  const axom::Array<double> z {{0., 0., 0., 0., 1., 1., 1., 1.}};
  setNodeValues(mesh["coordsets/coords/values/x"], x.view());
  setNodeValues(mesh["coordsets/coords/values/y"], y.view());
  setNodeValues(mesh["coordsets/coords/values/z"], z.view());

  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "hex";
  const axom::Array<conduit::index_t> connectivity {{0, 1, 3, 2, 4, 5, 7, 6}};
  setNodeValues(mesh["topologies/mesh/elements/connectivity"], connectivity.view());

  return mesh;
}

}  // namespace

TEST(quest_blueprint_quadrature_mesh, generate_closed_uniform_quad_mesh)
{
  conduit::Node mesh = makeQuadMesh();

  int sampleResolution[3] = {2, 3, 1};
  axom::quest::shaping::generateQuadraturePointMesh(mesh,
                                                    "mesh",
                                                    axom::execution_space<axom::SEQ_EXEC>::allocatorID(),
                                                    sampleResolution,
                                                    axom::numerics::QuadratureType::ClosedUniform);

  conduit::Node info;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  const conduit::Node& quadTopo = mesh["topologies/quadrature_points"];
  EXPECT_EQ(quadTopo["type"].as_string(), "unstructured");
  EXPECT_EQ(quadTopo["coordset"].as_string(), "quadrature_points");
  EXPECT_EQ(quadTopo["elements/shape"].as_string(), "point");

  namespace utils = axom::bump::utilities;
  const auto connView = utils::make_array_view<conduit::index_t>(
    mesh["topologies/quadrature_points/elements/connectivity"]);
  const auto sizesView = utils::make_array_view<conduit::index_t>(
    mesh["topologies/quadrature_points/elements/sizes"]);
  const auto offsetsView = utils::make_array_view<conduit::index_t>(
    mesh["topologies/quadrature_points/elements/offsets"]);
  const auto originalElementsView =
    utils::make_array_view<conduit::index_t>(mesh["fields/originalElements/values"]);
  const auto quadratureWeightsView =
    utils::make_array_view<double>(mesh["fields/quadratureWeights/values"]);

  const axom::Array<double> expectedX {{0., 1., 0., 1., 0., 1.}};
  const axom::Array<double> expectedY {{0., 0., 0.5, 0.5, 1., 1.}};
  const axom::Array<conduit::index_t> expectedConn {{0, 1, 2, 3, 4, 5}};
  const axom::Array<conduit::index_t> expectedSizes {{1, 1, 1, 1, 1, 1}};
  const axom::Array<conduit::index_t> expectedOffsets {{0, 1, 2, 3, 4, 5}};
  const axom::Array<conduit::index_t> expectedOriginalElements {{0, 0, 0, 0, 0, 0}};
  const axom::Array<double> expectedWeights {{1. / 12., 1. / 12., 1. / 3., 1. / 3., 1. / 12., 1. / 12.}};

  axom::bump::views::dispatch_explicit_coordset(
    mesh["coordsets/quadrature_points"], [&](auto coordsetView) {
      for(axom::IndexType i = 0; i < expectedX.size(); ++i)
      {
        EXPECT_NEAR(coordsetView[i][0], expectedX[i], 1e-12);
        EXPECT_NEAR(coordsetView[i][1], expectedY[i], 1e-12);
      }
    });
  EXPECT_TRUE(compareArrayView(expectedConn.view(), connView));
  EXPECT_TRUE(compareArrayView(expectedSizes.view(), sizesView));
  EXPECT_TRUE(compareArrayView(expectedOffsets.view(), offsetsView));
  EXPECT_TRUE(compareArrayView(expectedOriginalElements.view(), originalElementsView));
  for(axom::IndexType i = 0; i < expectedWeights.size(); ++i)
  {
    EXPECT_NEAR(expectedWeights[i], quadratureWeightsView[i], 1e-12);
  }
}

TEST(quest_blueprint_quadrature_mesh, generate_open_uniform_hex_mesh)
{
  conduit::Node mesh = makeHexMesh();

  int sampleResolution[3] = {2, 1, 2};
  axom::quest::shaping::generateQuadraturePointMesh(mesh,
                                                    "mesh",
                                                    axom::execution_space<axom::SEQ_EXEC>::allocatorID(),
                                                    sampleResolution,
                                                    axom::numerics::QuadratureType::OpenUniform);

  conduit::Node info;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/type")) << mesh.to_yaml();
  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/values/x")) << mesh.to_yaml();
  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/values/y")) << mesh.to_yaml();
  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/values/z")) << mesh.to_yaml();

  namespace utils = axom::bump::utilities;
  const auto originalElementsView =
    utils::make_array_view<conduit::index_t>(mesh["fields/originalElements/values"]);

  const axom::Array<double> expectedX {{1. / 3., 2. / 3., 1. / 3., 2. / 3.}};
  const axom::Array<double> expectedY {{0.5, 0.5, 0.5, 0.5}};
  const axom::Array<double> expectedZ {{1. / 3., 1. / 3., 2. / 3., 2. / 3.}};
  const axom::Array<conduit::index_t> expectedOriginalElements {{0, 0, 0, 0}};

  axom::bump::views::dispatch_explicit_coordset(
    mesh["coordsets/quadrature_points"], [&](auto coordsetView) {
      for(axom::IndexType i = 0; i < expectedX.size(); ++i)
      {
        EXPECT_NEAR(coordsetView[i][0], expectedX[i], 1e-6);
        EXPECT_NEAR(coordsetView[i][1], expectedY[i], 1e-6);
        EXPECT_NEAR(coordsetView[i][2], expectedZ[i], 1e-6);
      }
    });
  EXPECT_TRUE(compareArrayView(expectedOriginalElements.view(), originalElementsView));
}

TEST(quest_blueprint_quadrature_mesh, state_wrapper_generation_is_idempotent)
{
  conduit::Node mesh = makeQuadMesh();

  axom::quest::shaping::BlueprintState bpState;
  bpState.m_allocator_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  bpState.m_topology_name = "mesh";
  bpState.m_internal_node = mesh;

  int sampleResolution[3] = {2, 2, 1};
  axom::quest::shaping::generateSamplingPositions(
    bpState,
    sampleResolution,
    axom::numerics::QuadratureType::ClosedUniform);

  ASSERT_TRUE(bpState.m_internal_node.has_path("fields/originalElements/values"));
  conduit::Node savedOriginalElements;
  savedOriginalElements.set_external(
    bpState.m_internal_node["fields/originalElements/values"]);

  axom::quest::shaping::generateSamplingPositions(
    bpState,
    sampleResolution,
    axom::numerics::QuadratureType::OpenUniform);

  EXPECT_TRUE(bpState.m_internal_node.has_path("topologies/quadrature_points"));

  namespace utils = axom::bump::utilities;
  const auto originalElementsView = utils::make_array_view<conduit::index_t>(
    bpState.m_internal_node["fields/originalElements/values"]);
  const auto savedOriginalElementsView =
    utils::make_array_view<conduit::index_t>(savedOriginalElements);

  EXPECT_TRUE(compareArrayView(savedOriginalElementsView, originalElementsView));
}

TEST(quest_blueprint_quadrature_mesh, blueprint_state_field_helpers_support_replacement_ops)
{
  conduit::Node mesh = makeQuadMesh();

  axom::quest::shaping::BlueprintState bpState;
  bpState.m_allocator_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  bpState.m_topology_name = "mesh";
  bpState.m_internal_node = mesh;

  int sampleResolution[3] = {2, 2, 1};
  axom::quest::shaping::generateSamplingPositions(
    bpState,
    sampleResolution,
    axom::numerics::QuadratureType::ClosedUniform);

  conduit::Node& shapeField = bpState.m_internal_node["fields/inout_shape"];
  shapeField["association"] = "element";
  shapeField["topology"] = "quadrature_points";
  const axom::Array<double> shapeValues {{1., 0., 1., 0.}};
  setNodeValues(shapeField["values"], shapeValues.view());

  conduit::Node* materialField = bpState.createMaterialFunction("mat_inout_void");
  ASSERT_NE(materialField, nullptr);
  const axom::Array<double> materialValues {{1., 1., 0., 0.}};
  setNodeValues((*materialField)["values"], materialValues.view());

  conduit::Node* shapeCopy =
    axom::quest::shaping::cloneInOutFunction(bpState.getShapeFunction("inout_shape"));
  ASSERT_NE(shapeCopy, nullptr);

  axom::quest::shaping::replaceMaterial(shapeCopy, materialField, true);

  namespace utils = axom::bump::utilities;
  const auto replacedView = utils::make_array_view<double>((*materialField)["values"]);
  const axom::Array<double> expectedReplaced {{0., 1., 0., 0.}};
  EXPECT_TRUE(compareArrayView(expectedReplaced.view(), replacedView));

  conduit::Node* createdField = bpState.createMaterialFunction("mat_inout_created");
  ASSERT_NE(createdField, nullptr);
  axom::quest::shaping::copyShapeIntoMaterial(shapeCopy, createdField, false);

  const auto copiedView = utils::make_array_view<double>((*createdField)["values"]);
  EXPECT_TRUE(compareArrayView(shapeValues.view(), copiedView));

  delete shapeCopy;
}

TEST(quest_blueprint_quadrature_mesh, compute_volume_fractions_for_material_from_quadrature_weights)
{
  conduit::Node mesh = makeQuadMesh();

  axom::quest::shaping::BlueprintState bpState;
  bpState.m_allocator_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  bpState.m_topology_name = "mesh";
  bpState.m_internal_node = mesh;

  int sampleResolution[3] = {2, 2, 1};
  axom::quest::shaping::generateSamplingPositions(
    bpState,
    sampleResolution,
    axom::numerics::QuadratureType::ClosedUniform);

  conduit::Node* materialField = bpState.createMaterialFunction("mat_inout_test");
  ASSERT_NE(materialField, nullptr);

  const axom::Array<double> materialValues {{1., 0., 1., 0.}};
  setNodeValues((*materialField)["values"], materialValues.view());

  axom::quest::shaping::computeVolumeFractionsForMaterial(bpState, "mat_inout_test");

  ASSERT_TRUE(bpState.m_internal_node.has_path("fields/vol_frac_test/values"));
  namespace utils = axom::bump::utilities;
  const auto volFracValues =
    utils::make_array_view<double>(bpState.m_internal_node["fields/vol_frac_test/values"]);

  ASSERT_EQ(volFracValues.size(), 1);
  EXPECT_NEAR(volFracValues[0], 0.5, 1e-12);
}

TEST(quest_blueprint_quadrature_mesh, sampling_shaper_constructs_from_blueprint_node_and_group)
{
  conduit::Node mesh = makeQuadMesh();
  axom::klee::ShapeSet shapeSet;
  const auto policy = axom::runtime_policy::Policy::seq;
  const int allocatorId = axom::policyToDefaultAllocatorID(policy);

  axom::quest::SamplingShaper nodeShaper(policy, allocatorId, shapeSet, mesh, "mesh");

  axom::sidre::DataStore ds;
  auto* meshGroup = ds.getRoot()->createGroup("mesh");
  ASSERT_TRUE(meshGroup->importConduitTree(mesh));

  axom::quest::SamplingShaper groupShaper(policy, allocatorId, shapeSet, meshGroup, "mesh");
}

#endif
