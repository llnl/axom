// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/bump/views/StridedStructuredIndexing.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit_blueprint_mesh_examples.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

#include <string>
#include <vector>

namespace views = axom::bump::views;
namespace utils = axom::bump::utilities;

//----------------------------------------------------------------------
TEST(bump_views_indexing, strided_structured_indexing_2d)
{
  /*

   x---x---x---x---x---x---x
   |   |   |   |   |   |   |
   x---x---*---*---*---*---x   *=real node, x=ignored node, O=origin node 16
   |   |   |   |   |   |   |
   x---x---*---*---*---*---x
   |   |   |   |   |   |   |
   x---x---O---*---*---*---x
   |   |   |   |   |   |   |
   x---x---x---x---x---x---x
   |   |   |   |   |   |   |
   x---x---x---x---x---x---x

   */
  using Indexing = axom::bump::views::StridedStructuredIndexing<int, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  LogicalIndex dims {4, 3};  // window size in 4*3 elements in 7,6 overall
  LogicalIndex origin {2, 2};
  LogicalIndex stride {1, 7};
  Indexing indexing(dims, origin, stride);

  EXPECT_EQ(indexing.dimension(), 2);
  EXPECT_EQ(indexing.size(), dims[0] * dims[1]);

  // Iterate over local
  for(int j = 0; j < dims[1]; j++)
  {
    for(int i = 0; i < dims[0]; i++)
    {
      LogicalIndex logical {i, j};
      const auto flat = indexing.logicalIndexToIndex(logical);
      const auto logical2 = indexing.indexToLogicalIndex(flat);
      EXPECT_EQ(logical, logical2);

      EXPECT_EQ(logical, indexing.globalToLocal(indexing.localToGlobal(logical)));
      EXPECT_EQ(flat, indexing.globalToLocal(indexing.localToGlobal(flat)));
    }
  }

  // Iterate over global
  int index = 0;
  for(int j = 0; j < dims[1] + origin[1]; j++)
  {
    for(int i = 0; i < stride[1]; i++, index++)
    {
      LogicalIndex logical {i, j};
      const auto flat = indexing.globalToGlobal(logical);

      // flat should start at 0 and increase
      EXPECT_EQ(flat, index);

      // Global flat back to logical.
      LogicalIndex logical2 = indexing.globalToGlobal(flat);
      EXPECT_EQ(logical, logical2);

      // If we're in a valid region for the local window, try some other things.
      if(i >= origin[0] && i < origin[0] + dims[0] && j >= origin[1] && j < origin[1] + dims[1])
      {
        const auto logicalLocal = indexing.globalToLocal(logical);
        const auto flatLocal = indexing.globalToLocal(flat);

        // Flat local back to flat global
        EXPECT_EQ(flat, indexing.localToGlobal(flatLocal));

        // Logical local back to logical global
        EXPECT_EQ(logical, indexing.localToGlobal(logicalLocal));
      }
    }
  }

  // Check whether these local points exist.
  EXPECT_TRUE(indexing.contains(LogicalIndex {0, 0}));
  EXPECT_TRUE(indexing.contains(LogicalIndex {dims[0] - 1, dims[1] - 1}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 0}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 3}));
}

//----------------------------------------------------------------------
TEST(bump_views_indexing, strided_structured_indexing_3d)
{
  using Indexing3D = axom::bump::views::StridedStructuredIndexing<int, 3>;
  using LogicalIndex = typename axom::bump::views::StridedStructuredIndexing<int, 3>::LogicalIndex;
  LogicalIndex dims {4, 3, 3};  // window size in 4*3*3 elements in 6*5*5 overall
  LogicalIndex origin {2, 2, 2};
  LogicalIndex stride {1, 6, 30};
  Indexing3D indexing(dims, origin, stride);

  EXPECT_EQ(indexing.dimension(), 3);
  EXPECT_EQ(indexing.size(), dims[0] * dims[1] * dims[2]);

  const LogicalIndex logical0_0_0 {0, 0, 0};
  const auto index0_0_0 = indexing.logicalIndexToIndex(logical0_0_0);
  EXPECT_EQ(index0_0_0, 0);

  const LogicalIndex logical2_2_2 {2, 2, 2};
  const auto index2_2_2 = indexing.logicalIndexToIndex(logical2_2_2);
  EXPECT_EQ(index2_2_2, 2 + 2 * dims[0] + 2 * dims[0] * dims[1]);

  LogicalIndex logical = indexing.indexToLogicalIndex(index2_2_2);
  EXPECT_TRUE(logical == logical2_2_2);

  for(int k = 0; k < dims[2]; k++)
  {
    for(int j = 0; j < dims[1]; j++)
    {
      for(int i = 0; i < dims[0]; i++)
      {
        LogicalIndex logical {i, j, k};
        const auto flat = indexing.logicalIndexToIndex(logical);
        const auto logical2 = indexing.indexToLogicalIndex(flat);
        EXPECT_EQ(logical, logical2);
      }
    }
  }

  EXPECT_TRUE(indexing.contains(logical0_0_0));
  EXPECT_TRUE(indexing.contains(LogicalIndex {dims[0] - 1, dims[1] - 1, dims[2] - 1}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 0, 0}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 3, 0}));
}

//----------------------------------------------------------------------
// ElementFieldIndexing maps compact zone indices to element field indices.
//
// Conduit's strided_structured example fills each real element with its 1-based
// zone number and padding with 0. Therefore, values[indexing[z]] equals z + 1.
//----------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
void test_element_field_indexing_strided()
{
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::strided_structured<NDIMS>(hostMesh);

  // The field and topology have different strides.
  // The field mapping must use the metadata on ele_vals, not the topology's node strides.
  const auto topoStrides = hostMesh["topologies/mesh/elements/dims/strides"].as_int_accessor();
  const auto fieldStrides = hostMesh["fields/ele_vals/strides"].as_index_t_accessor();
  ASSERT_NE(topoStrides[1] - 1, fieldStrides[1]);

  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  const conduit::Node& n_field = deviceMesh["fields/ele_vals"];
  auto topoView =
    views::make_strided_structured_topology<NDIMS>::view(deviceMesh["topologies/mesh"]);
  const auto indexing = axom::bump::makeElementFieldIndexing(topoView, n_field);

  // Capture the mapping in the selected execution space
  // and gather the padded field at each compact zone index.
  const axom::IndexType nzones = topoView.numberOfZones();
  ASSERT_GT(nzones, 0);
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<double> gathered(nzones, nzones, allocatorID);
  axom::Array<axom::IndexType> indices(nzones, nzones, allocatorID);
  auto gatheredView = gathered.view();
  auto indicesView = indices.view();
  auto valuesView = utils::make_array_view<double>(n_field["values"]);
  axom::for_all<ExecSpace>(nzones, [=] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
    indicesView[zoneIndex] = indexing[zoneIndex];
    gatheredView[zoneIndex] = valuesView[indexing[zoneIndex]];
  });

  axom::Array<double> hostGathered(gathered, axom::execution_space<axom::SEQ_EXEC>::allocatorID());
  axom::Array<axom::IndexType> hostIndices(indices,
                                           axom::execution_space<axom::SEQ_EXEC>::allocatorID());

  // TableBasedExtractor uses this policy to slice element fields.
  using Indexing = typename decltype(topoView)::IndexingPolicy;
  axom::bump::SSElementFieldIndexing<Indexing> reference;
  reference.m_indexing = topoView.indexing();
  reference.update(hostMesh["fields/ele_vals"]);

  // The Conduit strided_structured example gives each real element the 1-based zone number.
  for(axom::IndexType z = 0; z < nzones; z++)
  {
    EXPECT_EQ(hostGathered[z], static_cast<double>(z + 1)) << "zone " << z;
    EXPECT_EQ(hostIndices[z], reference[z]) << "zone " << z;
  }
}

TEST(bump_views_indexing, element_field_indexing_strided_seq)
{
  test_element_field_indexing_strided<seq_exec, 2>();
  test_element_field_indexing_strided<seq_exec, 3>();
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_views_indexing, element_field_indexing_strided_omp)
{
  test_element_field_indexing_strided<omp_exec, 2>();
  test_element_field_indexing_strided<omp_exec, 3>();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_views_indexing, element_field_indexing_strided_cuda)
{
  test_element_field_indexing_strided<cuda_exec, 2>();
  test_element_field_indexing_strided<cuda_exec, 3>();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_views_indexing, element_field_indexing_strided_hip)
{
  test_element_field_indexing_strided<hip_exec, 2>();
  test_element_field_indexing_strided<hip_exec, 3>();
}
#endif

TEST(bump_views_indexing, element_field_indexing_strided_partial_metadata)
{
  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);
  auto topoView = views::make_strided_structured_topology<2>::view(mesh["topologies/mesh"]);
  const axom::IndexType nzones = topoView.numberOfZones();
  const auto dims = topoView.indexing().logicalDimensions();

  // No layout metadata: the field is stored in zone order.
  {
    conduit::Node n_field;
    n_field["association"] = "element";
    n_field["topology"] = "mesh";
    n_field["values"].set(std::vector<double>(nzones, 0.));
    const auto indexing = axom::bump::makeElementFieldIndexing(topoView, n_field);
    for(axom::IndexType z = 0; z < nzones; z++)
    {
      EXPECT_EQ(indexing[z], z);
    }
  }

  // With strides only, offsets default to 0.
  {
    conduit::Node n_field;
    n_field["association"] = "element";
    n_field["topology"] = "mesh";
    n_field["strides"].set(std::vector<conduit::index_t> {1, 5});
    n_field["values"].set(std::vector<double>(5 * dims[1], 0.));
    const auto indexing = axom::bump::makeElementFieldIndexing(topoView, n_field);
    for(axom::IndexType z = 0; z < nzones; z++)
    {
      EXPECT_EQ(indexing[z], (z % dims[0]) + (z / dims[0]) * 5);
    }
  }

  // With offsets only, strides default to the compact zone layout.
  {
    conduit::Node n_field;
    n_field["association"] = "element";
    n_field["topology"] = "mesh";
    n_field["offsets"].set(std::vector<conduit::index_t> {0, 1});
    n_field["values"].set(std::vector<double>(dims[0] * (dims[1] + 1), 0.));
    const auto indexing = axom::bump::makeElementFieldIndexing(topoView, n_field);
    for(axom::IndexType z = 0; z < nzones; z++)
    {
      EXPECT_EQ(indexing[z], z + dims[0]);
    }
  }
}

// Views whose element fields are stored in zone order return the identity.
TEST(bump_views_indexing, element_field_indexing_zone_order_views)
{
  const std::vector<std::pair<std::string, std::vector<int>>> cases {{"uniform", {4, 5, 3}},
                                                                     {"rectilinear", {4, 5}},
                                                                     {"structured", {4, 5, 3}},
                                                                     {"quads", {4, 5}},
                                                                     {"hexs", {3, 3, 3}},
                                                                     {"mixed_2d", {4, 5}}};
  for(const auto& c : cases)
  {
    SCOPED_TRACE(c.first);
    conduit::Node mesh;
    conduit::blueprint::mesh::examples::braid(c.first,
                                              c.second[0],
                                              c.second[1],
                                              c.second.size() > 2 ? c.second[2] : 0,
                                              mesh);
    if(!mesh.has_path("fields/radial"))
    {
      // The mixed-shape Braid examples have no element field.
      // Create one with the Blueprint element count.
      const auto nelem = conduit::blueprint::mesh::utils::topology::length(mesh["topologies/mesh"]);
      mesh["fields/radial/association"] = "element";
      mesh["fields/radial/topology"] = "mesh";
      mesh["fields/radial/values"].set(std::vector<double>(nelem, 0.));
    }
    const conduit::Node& n_field = mesh["fields/radial"];
    ASSERT_EQ(n_field["association"].as_string(), std::string("element"));

    int calls = 0;
    views::dispatch_topology(mesh["topologies/mesh"], [&](const std::string&, auto topoView) {
      calls++;
      const auto indexing = axom::bump::makeElementFieldIndexing(topoView, n_field);
      const axom::IndexType nzones = topoView.numberOfZones();
      EXPECT_EQ(nzones, n_field["values"].dtype().number_of_elements());
      for(axom::IndexType z = 0; z < nzones; z++)
      {
        EXPECT_EQ(indexing[z], z);
      }
    });
    EXPECT_EQ(calls, 1);
  }
}

TEST(bump_views_indexing, element_field_indexing_rejects_bad_fields)
{
  axom::slic::ScopedAbortToThrow abort_guard;

  conduit::Node strided;
  axom::blueprint::testing::data::strided_structured<2>(strided);
  auto stridedView = views::make_strided_structured_topology<2>::view(strided["topologies/mesh"]);

  // A vertex field is not an element field.
  EXPECT_THROW(axom::bump::makeElementFieldIndexing(stridedView, strided["fields/vert_vals"]),
               axom::slic::SlicAbortException);

  // A layout that maps the last zone past the end of the values.
  {
    conduit::Node n_field;
    n_field.set(strided["fields/ele_vals"]);
    conduit::index_t* offsets = n_field["offsets"].value();
    offsets[1] += 3;
    EXPECT_THROW(axom::bump::makeElementFieldIndexing(stridedView, n_field),
                 axom::slic::SlicAbortException);
  }

  // Nonpositive strides are not supported.
  {
    conduit::Node n_field;
    n_field.set(strided["fields/ele_vals"]);
    conduit::index_t* strides = n_field["strides"].value();
    strides[1] = 0;
    EXPECT_THROW(axom::bump::makeElementFieldIndexing(stridedView, n_field),
                 axom::slic::SlicAbortException);
  }

  // Layout metadata must have one entry per topology dimension.
  {
    conduit::Node n_field;
    n_field.set(strided["fields/ele_vals"]);
    n_field.remove("offsets");
    n_field["offsets"].set(std::vector<conduit::index_t> {0});
    EXPECT_THROW(axom::bump::makeElementFieldIndexing(stridedView, n_field),
                 axom::slic::SlicAbortException);
  }
  {
    conduit::Node n_field;
    n_field.set(strided["fields/ele_vals"]);
    n_field["strides"].set(std::vector<conduit::index_t> {1, 5, 25});
    EXPECT_THROW(axom::bump::makeElementFieldIndexing(stridedView, n_field),
                 axom::slic::SlicAbortException);
  }

  // Layout metadata on a view that cannot apply it.
  {
    conduit::Node mesh;
    conduit::blueprint::mesh::examples::braid("uniform", 4, 5, 0, mesh);
    conduit::Node& n_field = mesh["fields/radial"];
    n_field["offsets"].set(std::vector<conduit::index_t> {1, 1});
    n_field["strides"].set(std::vector<conduit::index_t> {1, 5});
    auto uniformView = views::make_uniform_topology<2>::view(mesh["topologies/mesh"]);
    EXPECT_THROW(axom::bump::makeElementFieldIndexing(uniformView, n_field),
                 axom::slic::SlicAbortException);
  }
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
