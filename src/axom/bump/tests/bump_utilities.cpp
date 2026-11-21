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
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const auto conduitAllocatorId = axom::sidre::ConduitMemory::axomAllocIdToConduit(
      allocatorID);

    EXPECT_TRUE(conduitAllocatorId > 0);

    constexpr int nValues = 100;
    conduit::Node n;
    n.set_allocator(conduitAllocatorId);
    n.set(conduit::DataType::int32(nValues));

    // Make sure we can store some values into the data that were allocated.
    auto nview = utils::make_array_view<int>(n);
    axom::for_all<ExecSpace>(nValues,
                             [=] AXOM_HOST_DEVICE(axom::IndexType index) { nview[index] = index; });

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

    // Test that nested Nodes inherit their parent's allocator.
    conduit::Node tree;
    tree.set_allocator(conduitAllocatorId);

    auto expectAllocator = [allocatorID](const conduit::Node& node) {
      EXPECT_EQ(axom::sidre::ConduitMemory::conduitAllocIdToAxom(node.allocator()),
                allocatorID);
    };

    conduit::Node& nested = tree["nested"];
    conduit::Node& childValues = tree["nested/values"];

    expectAllocator(tree);
    expectAllocator(nested);
    expectAllocator(childValues);

    childValues.set(conduit::DataType::int32(nValues));
    EXPECT_EQ(axom::getAllocatorIDFromPointer(childValues.data_ptr()), allocatorID);

    auto childView = utils::make_array_view<int>(childValues);
    axom::for_all<ExecSpace>(nValues, [=] AXOM_HOST_DEVICE(axom::IndexType index) {
      childView[index] = 2 * index;
    });

    std::vector<int> hostChildValues(nValues);
    axom::copy(hostChildValues.data(), childValues.data_ptr(), sizeof(int) * nValues);
    for(int i = 0; i < nValues; i++)
    {
      EXPECT_EQ(hostChildValues[i], 2 * i);
    }
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
template <typename ExecSpace>
struct test_conduit_compact_to
{
  static void test()
  {
    constexpr conduit::index_t n = 4;
    axom::Array<double> interleaved {{-1., 1.25, -2., 2.5, -3., 3.75, -4., 5., -5.}};

    const int sourceAllocatorID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    int destAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_UMPIRE)
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      const int hostAllocatorID = axom::detail::getAllocatorID<axom::MemorySpace::Host>();
      if(hostAllocatorID != sourceAllocatorID)
      {
        destAllocatorID = hostAllocatorID;
      }
    }
#endif
    const auto sourceConduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(sourceAllocatorID);

    conduit::Node host;
    host.set_allocator(sourceConduitAllocatorId);
    host["array"].set_external(conduit::DataType(conduit::DataType::FLOAT64_ID,
                                                 n,
                                                 sizeof(double),
                                                 2 * sizeof(double),
                                                 sizeof(double),
                                                 conduit::Endianness::DEFAULT_ID),
                               interleaved.data());
    host["scalar"] = 42;
    host["label"] = "compact_to";
    host["nested/ids"].set(std::vector<int> {7, 9});

    EXPECT_FALSE(host.is_compact());

    const auto conduitAllocatorId = axom::sidre::ConduitMemory::axomAllocIdToConduit(destAllocatorID);
    conduit::Node device;
    device.set_allocator(conduitAllocatorId);
    host.compact_to(device);

    EXPECT_TRUE(device.is_compact());
    EXPECT_TRUE(device.is_contiguous());
    EXPECT_EQ(device.total_bytes_allocated(), device.total_bytes_compact());
    EXPECT_EQ(device.total_bytes_compact(), host.total_bytes_compact());
    EXPECT_EQ(axom::sidre::ConduitMemory::conduitAllocIdToAxom(device.allocator()), destAllocatorID);
    EXPECT_EQ(axom::getAllocatorIDFromPointer(device.data_ptr()), destAllocatorID);
    EXPECT_EQ(axom::getAllocatorIDFromPointer(device["array"].data_ptr()), destAllocatorID);
    EXPECT_EQ(axom::getAllocatorIDFromPointer(device["scalar"].data_ptr()), destAllocatorID);
    EXPECT_EQ(axom::getAllocatorIDFromPointer(device["label"].data_ptr()), destAllocatorID);
    EXPECT_EQ(axom::getAllocatorIDFromPointer(device["nested/ids"].data_ptr()), destAllocatorID);

    EXPECT_EQ(device["array"].dtype().stride(), device["array"].dtype().element_bytes());
    EXPECT_EQ(device["nested/ids"].dtype().stride(), device["nested/ids"].dtype().element_bytes());

    auto arrayView = utils::make_array_view<double>(device["array"]);
    auto scalarView = utils::make_array_view<int>(device["scalar"]);
    auto idsView = utils::make_array_view<int>(device["nested/ids"]);

    axom::ReduceSum<ExecSpace, double> arraySumReduce(0.);
    axom::for_all<ExecSpace>(n, [=] AXOM_HOST_DEVICE(axom::IndexType index) {
      arraySumReduce += arrayView[index];
    });
    EXPECT_DOUBLE_EQ(arraySumReduce.get(), 12.5);

    axom::ReduceSum<ExecSpace, int> scalarReduce(0);
    axom::for_all<ExecSpace>(1, [=] AXOM_HOST_DEVICE(axom::IndexType) {
      scalarReduce += scalarView[0];
    });
    EXPECT_EQ(scalarReduce.get(), 42);

    axom::ReduceSum<ExecSpace, int> idsReduce(0);
    axom::for_all<ExecSpace>(idsView.size(), [=] AXOM_HOST_DEVICE(axom::IndexType index) {
      idsReduce += idsView[index];
    });
    EXPECT_EQ(idsReduce.get(), 16);

    conduit::Node newHost;
    device.compact_to(newHost);

    EXPECT_TRUE(newHost.is_compact());
    EXPECT_TRUE(newHost.is_contiguous());
    EXPECT_EQ(newHost["array"].dtype().stride(), newHost["array"].dtype().element_bytes());

    auto newHostArray = newHost["array"].as_double_accessor();
    ASSERT_EQ(newHostArray.number_of_elements(), n);
    EXPECT_DOUBLE_EQ(newHostArray[0], 1.25);
    EXPECT_DOUBLE_EQ(newHostArray[1], 2.5);
    EXPECT_DOUBLE_EQ(newHostArray[2], 3.75);
    EXPECT_DOUBLE_EQ(newHostArray[3], 5.);
    EXPECT_EQ(newHost["scalar"].to_int(), 42);
    EXPECT_EQ(newHost["label"].as_string(), std::string("compact_to"));

    auto newHostIds = newHost["nested/ids"].as_int_accessor();
    ASSERT_EQ(newHostIds.number_of_elements(), 2);
    EXPECT_EQ(newHostIds[0], 7);
    EXPECT_EQ(newHostIds[1], 9);
  }
};

TEST(bump_utilities, compact_to_seq) { test_conduit_compact_to<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, compact_to_omp) { test_conduit_compact_to<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, compact_to_cuda) { test_conduit_compact_to<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, compact_to_hip) { test_conduit_compact_to<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_conduit_host_read
{
  static void test(bool checkDeviceAllocated = false)
  {
    conduit::Node host;
    host["int_scalar"] = 42;
    host["double_scalar"] = 3.25;
    host["label"] = "host_read";
    host["int_array"].set(std::vector<int> {5, 10, 15, 20});
    host["double_array"].set(std::vector<double> {1.25, 2.5, 3.75});

    EXPECT_EQ(utils::hostRead<int>(host["int_scalar"]), 42);
    EXPECT_EQ(utils::hostReadInt(host["int_scalar"]), 42);
    EXPECT_DOUBLE_EQ(utils::hostRead<double>(host["double_scalar"]), 3.25);
    EXPECT_DOUBLE_EQ(utils::hostReadDouble(host["double_scalar"]), 3.25);
    EXPECT_EQ(utils::hostRead<std::string>(host["label"]), std::string("host_read"));
    EXPECT_EQ(utils::hostReadString(host["label"]), std::string("host_read"));

    conduit::Node device;
    utils::copy<ExecSpace>(device, host);

    EXPECT_EQ(utils::hostReadInt(device["int_scalar"]), 42);
    EXPECT_DOUBLE_EQ(utils::hostReadDouble(device["double_scalar"]), 3.25);
    EXPECT_EQ(utils::hostReadString(device["label"]), std::string("host_read"));

    const auto intValues = utils::hostReadArray<int>(device["int_array"]);
    ASSERT_EQ(intValues.size(), 4);
    EXPECT_EQ(intValues[0], 5);
    EXPECT_EQ(intValues[1], 10);
    EXPECT_EQ(intValues[2], 15);
    EXPECT_EQ(intValues[3], 20);

    const auto doubleValues = utils::hostReadArray<double>(device["double_array"]);
    ASSERT_EQ(doubleValues.size(), 3);
    EXPECT_DOUBLE_EQ(doubleValues[0], 1.25);
    EXPECT_DOUBLE_EQ(doubleValues[1], 2.5);
    EXPECT_DOUBLE_EQ(doubleValues[2], 3.75);

    if(checkDeviceAllocated)
    {
      EXPECT_TRUE(utils::isDeviceAllocated(device["int_scalar"]));
      EXPECT_TRUE(utils::isDeviceAllocated(device["double_scalar"]));
      EXPECT_TRUE(utils::isDeviceAllocated(device["label"]));
      EXPECT_TRUE(utils::isDeviceAllocated(device["int_array"]));
      EXPECT_TRUE(utils::isDeviceAllocated(device["double_array"]));
    }
  }
};

TEST(bump_utilities, host_read_seq) { test_conduit_host_read<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, host_read_omp) { test_conduit_host_read<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, host_read_cuda) { test_conduit_host_read<cuda_exec>::test(true); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, host_read_hip) { test_conduit_host_read<hip_exec>::test(true); }
#endif

TEST(bump_utilities, host_read_rejects_non_leaf_nodes)
{
  axom::slic::ScopedAbortToThrow abort_guard;

  conduit::Node empty, object, list, array;
  object["child"] = 5;
  list.append() = 7;
  array.set(std::vector<int> {1, 2});

  EXPECT_THROW(utils::hostRead<int>(empty), axom::slic::SlicAbortException);
  EXPECT_THROW(utils::hostRead<int>(object), axom::slic::SlicAbortException);
  EXPECT_THROW(utils::hostRead<int>(list), axom::slic::SlicAbortException);
  EXPECT_THROW(utils::hostRead<int>(array), axom::slic::SlicAbortException);
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_conduit_host_read_nested_copy
{
  static void test()
  {
    conduit::Node host;
    host["metadata/count"] = 42;
    host["metadata/scale"] = 3.25;
    host["metadata/label"] = "nested_host_read";
    host["fields/ids"].set(std::vector<int> {2, 4, 6, 8});
    host["fields/weights"].set(std::vector<double> {1.5, 2.5, 3.5});

    conduit::Node copied;
    utils::copy<ExecSpace>(copied, host);

    EXPECT_EQ(utils::hostReadInt(copied["metadata/count"]), 42);
    EXPECT_DOUBLE_EQ(utils::hostReadDouble(copied["metadata/scale"]), 3.25);
    EXPECT_EQ(utils::hostReadString(copied["metadata/label"]),
              std::string("nested_host_read"));

    const auto ids = utils::hostReadArray<int>(copied["fields/ids"]);
    ASSERT_EQ(ids.size(), 4);
    EXPECT_EQ(ids[0], 2);
    EXPECT_EQ(ids[1], 4);
    EXPECT_EQ(ids[2], 6);
    EXPECT_EQ(ids[3], 8);

    const auto weights = utils::hostReadArray<double>(copied["fields/weights"]);
    ASSERT_EQ(weights.size(), 3);
    EXPECT_DOUBLE_EQ(weights[0], 1.5);
    EXPECT_DOUBLE_EQ(weights[1], 2.5);
    EXPECT_DOUBLE_EQ(weights[2], 3.5);
  }
};

TEST(bump_utilities, host_read_nested_copy_seq)
{
  test_conduit_host_read_nested_copy<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, host_read_nested_copy_omp)
{
  test_conduit_host_read_nested_copy<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, host_read_nested_copy_cuda)
{
  test_conduit_host_read_nested_copy<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, host_read_nested_copy_hip)
{
  test_conduit_host_read_nested_copy<hip_exec>::test();
}
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

  axom::for_all<seq_exec>(n, [=] AXOM_HOST_DEVICE(axom::IndexType index) {
    view[index] = static_cast<double>((index + 1) * 100);
  });

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

TEST(bump_utilities, fill_from_node_lengths)
{
  axom::slic::ScopedAbortToThrow abort_guard;

  conduit::Node n;
  n["exact"].set(std::vector<int> {2, 3});
  n["long"].set(std::vector<int> {4, 5, 6});
  n["short"].set(std::vector<int> {7});
  n["double_long"].set(std::vector<double> {1.5, 2.5, 3.5});

  axom::StackArray<int, 2> arr {{-1, -1}};
  EXPECT_FALSE(utils::fillFromNode(n, "missing", arr));
  EXPECT_EQ(arr[0], -1);

  EXPECT_TRUE(utils::fillFromNode(n, "exact", arr));
  EXPECT_EQ(arr[0], 2);
  EXPECT_EQ(arr[1], 3);

  // Extra values are ignored.
  EXPECT_TRUE(utils::fillFromNode(n, "long", arr));
  EXPECT_EQ(arr[0], 4);
  EXPECT_EQ(arr[1], 5);

  axom::StackArray<double, 2> darr {{-1., -1.}};
  EXPECT_TRUE(utils::fillFromNode(n, "double_long", darr));
  EXPECT_DOUBLE_EQ(darr[0], 1.5);
  EXPECT_DOUBLE_EQ(darr[1], 2.5);

  // Too few values would read past the end of the node's data.
  EXPECT_THROW(utils::fillFromNode(n, "short", arr), axom::slic::SlicAbortException);
}

template <typename ExecSpace>
struct test_fill_from_node_after_copy
{
  static void test()
  {
    conduit::Node host;
    host["ints"].set(std::vector<int> {8, 9, 10});
    host["weights"].set(std::vector<double> {0.25, 0.5});

    conduit::Node copied;
    utils::copy<ExecSpace>(copied, host);

    axom::StackArray<int, 2> ints {{-1, -1}};
    EXPECT_TRUE(utils::fillFromNode(copied, "ints", ints));
    EXPECT_EQ(ints[0], 8);
    EXPECT_EQ(ints[1], 9);

    axom::StackArray<double, 2> weights {{-1., -1.}};
    EXPECT_TRUE(utils::fillFromNode(copied, "weights", weights));
    EXPECT_DOUBLE_EQ(weights[0], 0.25);
    EXPECT_DOUBLE_EQ(weights[1], 0.5);
  }
};

TEST(bump_utilities, fill_from_node_after_copy_seq)
{
  test_fill_from_node_after_copy<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, fill_from_node_after_copy_omp)
{
  test_fill_from_node_after_copy<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, fill_from_node_after_copy_cuda)
{
  test_fill_from_node_after_copy<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, fill_from_node_after_copy_hip)
{
  test_fill_from_node_after_copy<hip_exec>::test();
}
#endif

TEST(bump_utilities, conduit_single_value_vector_matches_scalar_representation)
{
  conduit::Node scalar, vector_int, vector_float;
  scalar.set(7);
  vector_int.set(std::vector<int> {7});
  vector_float.set(std::vector<float> {3.5f});

  EXPECT_EQ(scalar.dtype().id(), vector_int.dtype().id());
  EXPECT_EQ(scalar.dtype().number_of_elements(), 1);
  EXPECT_EQ(vector_int.dtype().number_of_elements(), 1);
  EXPECT_EQ(scalar.dtype().element_bytes(), vector_int.dtype().element_bytes());
  EXPECT_EQ(scalar.total_bytes_allocated(), vector_int.total_bytes_allocated());
  EXPECT_FALSE(scalar.is_data_external());
  EXPECT_FALSE(vector_int.is_data_external());

  EXPECT_EQ(vector_int.to_int(), 7);
  EXPECT_EQ(vector_float.dtype().number_of_elements(), 1);
  EXPECT_FLOAT_EQ(vector_float.to_float(), 3.5f);
}

//------------------------------------------------------------------------------
// Bump uses topology node ids to access flat vertex fields.
// A strided field's offsets and strides must match the topology.
//------------------------------------------------------------------------------
TEST(bump_utilities, validate_vertex_field_indexing_compact)
{
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
  axom::slic::ScopedAbortToThrow abort_guard;

  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);

  // Change one field offset.
  // Without this check, Bump would read shifted values and place the contour incorrectly.
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

    EXPECT_THROW(utils::validateVertexFieldIndexing(perturbed["topologies/mesh"],
                                                    perturbed["fields/vert_vals"],
                                                    "vert_vals"),
                 axom::slic::SlicAbortException);
  }

  // Strided fields require strides. Offsets are optional.
  {
    conduit::Node offsets_only;
    offsets_only.set(mesh);
    offsets_only["fields/vert_vals"].remove("strides");
    EXPECT_THROW(utils::validateVertexFieldIndexing(offsets_only["topologies/mesh"],
                                                    offsets_only["fields/vert_vals"],
                                                    "vert_vals"),
                 axom::slic::SlicAbortException);
  }

  // Reject field layout metadata on a compact topology.
  {
    conduit::Node compact;
    conduit::blueprint::mesh::examples::braid("hexs", 4, 4, 4, compact);
    compact["fields/braid/offsets"].set(std::vector<int> {1, 0, 0});
    compact["fields/braid/strides"].set(std::vector<int> {1, 4, 16});

    EXPECT_THROW(utils::validateVertexFieldIndexing(compact["topologies/mesh"],
                                                    compact["fields/braid"],
                                                    "braid"),
                 axom::slic::SlicAbortException);
  }
}

TEST(bump_utilities, validate_vertex_field_indexing_resolves_defaults)
{
  axom::slic::ScopedAbortToThrow abort_guard;

  // The topology omits offsets, so they default to 0. The field explicitly uses [0, 0].
  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);
  mesh["topologies/mesh/elements/dims"].remove("offsets");
  conduit::Node& n_field = mesh["fields/vert_vals"];
  conduit::index_t* offsets = n_field["offsets"].value();
  offsets[0] = 0;
  offsets[1] = 0;
  EXPECT_NO_THROW(utils::validateVertexFieldIndexing(mesh["topologies/mesh"], n_field, "vert_vals"));

  // An explicit compact field layout matches a compact structured topology.
  conduit::Node compact;
  conduit::blueprint::mesh::examples::braid("structured", 4, 5, 0, compact);
  conduit::Node& n_compactField = compact["fields/braid"];
  n_compactField["offsets"].set(std::vector<conduit::index_t> {0, 0});
  n_compactField["strides"].set(std::vector<conduit::index_t> {1, 4});
  EXPECT_NO_THROW(
    utils::validateVertexFieldIndexing(compact["topologies/mesh"], n_compactField, "braid"));

  // Strided topology metadata requires strides, like strided field metadata.
  conduit::Node offsets_only_topology;
  offsets_only_topology.set(mesh);
  offsets_only_topology["topologies/mesh/elements/dims/offsets"].set(
    std::vector<conduit::index_t> {0, 0});
  offsets_only_topology["topologies/mesh/elements/dims"].remove("strides");
  EXPECT_THROW(utils::validateVertexFieldIndexing(offsets_only_topology["topologies/mesh"],
                                                  offsets_only_topology["fields/vert_vals"],
                                                  "vert_vals"),
               axom::slic::SlicAbortException);

  // A different j-stride on the same compact topology is rejected.
  conduit::index_t* strides = n_compactField["strides"].value();
  strides[1] = 5;
  EXPECT_THROW(
    utils::validateVertexFieldIndexing(compact["topologies/mesh"], n_compactField, "braid"),
    axom::slic::SlicAbortException);
}

TEST(bump_utilities, validate_vertex_field_indexing_rejects_stride_mismatch)
{
  axom::slic::ScopedAbortToThrow abort_guard;

  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);

  // Same offsets, different j-stride.
  conduit::Node& n_field = mesh["fields/vert_vals"];
  conduit::index_t* strides = n_field["strides"].value();
  strides[1] += 1;

  EXPECT_THROW(utils::validateVertexFieldIndexing(mesh["topologies/mesh"], n_field, "vert_vals"),
               axom::slic::SlicAbortException);
}

TEST(bump_utilities, validate_vertex_field_indexing_rejects_topology_i_stride)
{
  axom::slic::ScopedAbortToThrow abort_guard;

  // StructuredTopologyView::zone() assumes adjacent corners differ by 1 in i.
  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);
  conduit::int_accessor topoStrides = mesh["topologies/mesh/elements/dims/strides"].as_int_accessor();
  const std::vector<int> wide {2, topoStrides[1] * 2};
  mesh["topologies/mesh/elements/dims/strides"].set(wide);
  mesh["fields/vert_vals"].remove("offsets");
  mesh["fields/vert_vals"].remove("strides");

  EXPECT_THROW(utils::validateVertexFieldIndexing(mesh["topologies/mesh"],
                                                  mesh["fields/vert_vals"],
                                                  "vert_vals"),
               axom::slic::SlicAbortException);
}

TEST(bump_utilities, validate_vertex_field_indexing_strided_topology_plain_field)
{
  // Without field layout metadata, topology node ids index the values directly,
  // which is correct when the field is sized like the coordset.
  conduit::Node mesh;
  axom::blueprint::testing::data::strided_structured<2>(mesh);
  conduit::Node& n_field = mesh["fields/vert_vals"];
  n_field.remove("offsets");
  n_field.remove("strides");

  EXPECT_NO_FATAL_FAILURE(
    utils::validateVertexFieldIndexing(mesh["topologies/mesh"], n_field, "vert_vals"));
}

template <typename ExecSpace>
void test_validate_vertex_field_indexing_memory()
{
  axom::slic::ScopedAbortToThrow abort_guard;

  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::strided_structured<3>(hostMesh);
  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  EXPECT_NO_THROW(utils::validateVertexFieldIndexing(deviceMesh["topologies/mesh"],
                                                     deviceMesh["fields/vert_vals"],
                                                     "vert_vals"));

  // Shift the field's k-offset, then copy the bad metadata to ExecSpace memory.
  conduit::index_t* offsets = hostMesh["fields/vert_vals/offsets"].value();
  offsets[2] += 1;
  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  EXPECT_THROW(utils::validateVertexFieldIndexing(deviceMesh["topologies/mesh"],
                                                  deviceMesh["fields/vert_vals"],
                                                  "vert_vals"),
               axom::slic::SlicAbortException);
}

TEST(bump_utilities, validate_vertex_field_indexing_memory_seq)
{
  test_validate_vertex_field_indexing_memory<seq_exec>();
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_utilities, validate_vertex_field_indexing_memory_omp)
{
  test_validate_vertex_field_indexing_memory<omp_exec>();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_utilities, validate_vertex_field_indexing_memory_cuda)
{
  test_validate_vertex_field_indexing_memory<cuda_exec>();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_utilities, validate_vertex_field_indexing_memory_hip)
{
  test_validate_vertex_field_indexing_memory<hip_exec>();
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();
  return result;
}
