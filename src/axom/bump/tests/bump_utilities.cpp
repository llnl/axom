// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit_blueprint_mesh_examples.hpp>

#include <iostream>
#include <algorithm>
#include <vector>

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_conduit_allocate
{
  static void test()
  {
    const auto conduitAllocatorId = axom::sidre::ConduitMemory::axomAllocIdToConduit(
      axom::execution_space<ExecSpace>::allocatorID());

    EXPECT_TRUE(conduitAllocatorId > 0);

    constexpr int nValues = 100;
    conduit::Node n;
    n.set_allocator(conduitAllocatorId);
    n.set(conduit::DataType::int32(nValues));

    // Make sure we can store some values into the data that were allocated.
    auto nview = utils::make_array_view<int>(n);
    axom::for_all<ExecSpace>(nValues, AXOM_LAMBDA(axom::IndexType index) { nview[index] = index; });

    EXPECT_EQ(n.dtype().number_of_elements(), nValues);

    // Get the values back to the host.
    std::vector<int> hostValues(nValues);
    axom::copy(hostValues.data(), n.data_ptr(), sizeof(int) * nValues);

    // Check that the values were set.
    for(int i = 0; i < nValues; i++)
    {
      EXPECT_EQ(hostValues[i], i);
    }

    // Check zero allocation.
    n.reset();
    n.set_allocator(conduitAllocatorId);
    n.set(conduit::DataType::int32(0));
    EXPECT_EQ(n.dtype().number_of_elements(), 0);
  }
};

TEST(bump_utilities, allocate_seq) { test_conduit_allocate<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, allocate_omp) { test_conduit_allocate<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, allocate_cuda) { test_conduit_allocate<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, allocate_hip) { test_conduit_allocate<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(bump_utilities, make_array_view_interleaved_seq)
{
  constexpr conduit::index_t n = 4;
  axom::Array<double> interleaved {{-1., 10., -2., 20., -3., 30., -4., 40., -5.}};
  conduit::Node n_data;
  n_data.set_external(conduit::DataType(conduit::DataType::FLOAT64_ID,
                                        n,
                                        sizeof(double),
                                        2 * sizeof(double),
                                        sizeof(double),
                                        conduit::Endianness::DEFAULT_ID),
                      interleaved.data());

  auto view = utils::make_array_view<double>(n_data);
  EXPECT_EQ(view.size(), n);

  axom::for_all<seq_exec>(
    n,
    AXOM_LAMBDA(axom::IndexType index) { view[index] = static_cast<double>((index + 1) * 100); });

  EXPECT_EQ(interleaved[0], -1.);
  EXPECT_EQ(interleaved[1], 100.);
  EXPECT_EQ(interleaved[2], -2.);
  EXPECT_EQ(interleaved[3], 200.);
  EXPECT_EQ(interleaved[4], -3.);
  EXPECT_EQ(interleaved[5], 300.);
  EXPECT_EQ(interleaved[6], -4.);
  EXPECT_EQ(interleaved[7], 400.);
  EXPECT_EQ(interleaved[8], -5.);
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_copy_braid
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Run some minmax operations on device (proves that the data was in the right place) and check the results.

    constexpr double eps = 1.e-7;

    auto x = bump::MinMax<ExecSpace, double>::execute(deviceMesh["coordsets/coords/values/x"]);
    //std::cout << std::setw(16) << "x={" << x.first << ", " << x.second << "}\n";
    EXPECT_NEAR(x.first, -10., eps);
    EXPECT_NEAR(x.second, 10., eps);

    auto y = bump::MinMax<ExecSpace, double>::execute(deviceMesh["coordsets/coords/values/y"]);
    //std::cout << std::setw(16) << "y={" << y.first << ", " << y.second << "}\n";
    EXPECT_NEAR(y.first, -10., eps);
    EXPECT_NEAR(y.second, 10., eps);

    auto c =
      bump::MinMax<ExecSpace, double>::execute(deviceMesh["topologies/mesh/elements/connectivity"]);
    //std::cout << std::setw(16) << "conn={" << c.first << ", " << c.second << "}\n";
    EXPECT_NEAR(c.first, 0., eps);
    EXPECT_NEAR(c.second, 999., eps);

    auto r = bump::MinMax<ExecSpace, double>::execute(deviceMesh["fields/radial/values"]);
    //std::cout << std::setw(16) << "radial={" << r.first << ", " << r.second << "}\n";
    EXPECT_NEAR(r.first, 19.2450089729875, eps);
    EXPECT_NEAR(r.second, 173.205080756888, eps);

    // Test copying an empty Node from device to host.
    conduit::Node emptyDeviceMesh, emptyHostMesh;
    utils::copy<ExecSpace>(emptyHostMesh, emptyDeviceMesh);
    EXPECT_TRUE(emptyHostMesh.dtype().is_empty());
  }

  static void create(conduit::Node& mesh)
  {
    const int d[3] = {10, 10, 10};
    conduit::blueprint::mesh::examples::braid("hexs", d[0], d[1], d[2], mesh);
  }
};

TEST(bump_utilities, copy_seq) { test_copy_braid<seq_exec>::test(); }

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, copy_omp) { test_copy_braid<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, copy_cuda) { test_copy_braid<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, copy_hip) { test_copy_braid<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
// Check the indexing contract between a topology and its vertex fields.
//
// Bump indexes a flat field with the node ids from TopologyView::zone().
// For a strided structured topology, any field-specific offsets and strides
// must match the topology.
//------------------------------------------------------------------------------
TEST(bump_utilities, validate_vertex_field_indexing_compact)
{
  // This compact mesh has no offset or stride metadata.
  conduit::Node mesh;
  conduit::blueprint::mesh::examples::braid("hexs", 4, 4, 4, mesh);

  EXPECT_NO_FATAL_FAILURE(
    utils::validateVertexFieldIndexing(mesh["topologies/mesh"], mesh["fields/braid"], "braid"));
}

TEST(bump_utilities, validate_vertex_field_indexing_strided_matching)
{
  // Conduit's strided_structured example gives the field and topology the same layout.
  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);

  ASSERT_TRUE(mesh.has_path("topologies/mesh/elements/dims/offsets"));
  ASSERT_TRUE(mesh.has_path("fields/vert_vals/offsets"));

  EXPECT_NO_FATAL_FAILURE(utils::validateVertexFieldIndexing(mesh["topologies/mesh"],
                                                             mesh["fields/vert_vals"],
                                                             "vert_vals"));
}

TEST(bump_utilities, validate_vertex_field_indexing_rejects_mismatch)
{
  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);

  // Change one field offset. Without this check, Bump would read shifted values
  // and place the contour incorrectly.
  {
    conduit::Node perturbed;
    perturbed.set(mesh);
    conduit::int_accessor off = perturbed["fields/vert_vals/offsets"].as_int_accessor();
    std::vector<int> shifted;
    for(conduit::index_t i = 0; i < off.number_of_elements(); i++)
    {
      shifted.push_back(off[i] + (i == 0 ? 1 : 0));
    }
    perturbed["fields/vert_vals/offsets"].set(shifted);

    // SLIC's default handler aborts.
    // SimpleLogger writes to stdout, while GTest death tests capture stderr,
    // so add a stream inside the child.
    EXPECT_DEATH_IF_SUPPORTED(
      {
        axom::slic::addStreamToAllMsgLevels(
          new axom::slic::GenericOutputStream(&std::cerr, "[<LEVEL>] <MESSAGE>\n"));
        utils::validateVertexFieldIndexing(perturbed["topologies/mesh"],
                                           perturbed["fields/vert_vals"],
                                           "vert_vals");
      },
      "but its topology has offsets");
  }

  // Reject field layout metadata on a compact topology.
  {
    conduit::Node compact;
    conduit::blueprint::mesh::examples::braid("hexs", 4, 4, 4, compact);
    compact["fields/braid/offsets"].set(std::vector<int> {1, 0, 0});
    compact["fields/braid/strides"].set(std::vector<int> {1, 4, 16});

    EXPECT_DEATH_IF_SUPPORTED(
      {
        axom::slic::addStreamToAllMsgLevels(
          new axom::slic::GenericOutputStream(&std::cerr, "[<LEVEL>] <MESSAGE>\n"));
        utils::validateVertexFieldIndexing(compact["topologies/mesh"],
                                           compact["fields/braid"],
                                           "braid");
      },
      "its topology does not");
  }
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();
  return result;
}
