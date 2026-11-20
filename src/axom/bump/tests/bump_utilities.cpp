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

#include <cstdint>

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

  axom::StackArray<int, 2> arr {{-1, -1}};
  EXPECT_FALSE(utils::fillFromNode(n, "missing", arr));
  EXPECT_EQ(arr[0], -1);

  EXPECT_TRUE(utils::fillFromNode(n, "exact", arr));
  EXPECT_EQ(arr[0], 2);
  EXPECT_EQ(arr[1], 3);

  // Extra values are ignored.
  EXPECT_TRUE(utils::fillFromNode(n, "long", arr, true));
  EXPECT_EQ(arr[0], 4);
  EXPECT_EQ(arr[1], 5);

  // Too few values would read past the end of the node's data.
  EXPECT_THROW(utils::fillFromNode(n, "short", arr), axom::slic::SlicAbortException);
  EXPECT_THROW(utils::fillFromNode(n, "short", arr, true), axom::slic::SlicAbortException);
}

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

namespace
{
void expectSameLeafType(const conduit::Node& src, const conduit::Node& copied)
{
  EXPECT_EQ(src.dtype().id(), copied.dtype().id());
  EXPECT_EQ(src.dtype().number_of_elements(), copied.dtype().number_of_elements());
  EXPECT_EQ(src.dtype().element_bytes(), copied.dtype().element_bytes());
}

[[maybe_unused]] void expectPointerInBucket(const conduit::Node& leaf, const conduit::Node& bucket)
{
  const auto leaf_begin = reinterpret_cast<std::uintptr_t>(leaf.data_ptr());
  const auto leaf_end = leaf_begin + static_cast<std::uintptr_t>(leaf.dtype().bytes_compact());
  const auto bucket_begin = reinterpret_cast<std::uintptr_t>(bucket.data_ptr());
  const auto bucket_end = bucket_begin + static_cast<std::uintptr_t>(bucket.dtype().bytes_compact());

  EXPECT_GE(leaf_begin, bucket_begin);
  EXPECT_LE(leaf_end, bucket_end);
}
}  // namespace

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_copy_single_value_blueprint_leaves
{
  static void create(conduit::Node& n)
  {
    n["coordsets/coords/type"] = "explicit";
    n["coordsets/coords/values/x"].set(std::vector<float> {1.5f});

    n["topologies/mesh/type"] = "unstructured";
    n["topologies/mesh/coordset"] = "coords";
    n["topologies/mesh/elements/shape"] = "point";
    n["topologies/mesh/elements/connectivity"].set(std::vector<int> {0});
    n["topologies/mesh/elements/sizes"].set(std::vector<int> {1});
    n["topologies/mesh/elements/offsets"].set(std::vector<int> {0});

    n["fields/one/topology"] = "mesh";
    n["fields/one/association"] = "element";
    n["fields/one/values"].set(std::vector<double> {2.5});

    n["fields/mixed/topology"] = "mesh";
    n["fields/mixed/association"] = "element";
    n["fields/mixed/matset"] = "mat";
    n["fields/mixed/matset_values"].set(std::vector<float> {3.5f});

    n["fields/mixed_nested/topology"] = "mesh";
    n["fields/mixed_nested/association"] = "element";
    n["fields/mixed_nested/matset"] = "mat_multi";
    n["fields/mixed_nested/matset_values/A/values"].set(std::vector<float> {4.5f});

    n["fields/flag/topology"] = "mesh";
    n["fields/flag/association"] = "element";
    n["fields/flag/values"].set(std::vector<conduit::uint8> {9});

    n["fields/small/topology"] = "mesh";
    n["fields/small/association"] = "element";
    n["fields/small/values"].set(std::vector<conduit::uint16> {12});

    n["matsets/mat/topology"] = "mesh";
    n["matsets/mat/material_ids"].set(std::vector<int> {22});
    n["matsets/mat/volume_fractions"].set(std::vector<float> {1.f});
    n["matsets/mat/sizes"].set(std::vector<int> {1});
    n["matsets/mat/offsets"].set(std::vector<int> {0});
    n["matsets/mat/indices"].set(std::vector<int> {0});

    n["matsets/mat_multi/topology"] = "mesh";
    n["matsets/mat_multi/volume_fractions/A/values"].set(std::vector<float> {0.75f});
    n["matsets/mat_multi/volume_fractions/A/indices"].set(std::vector<int> {0});

    n["metadata/scalar_int"] = 7;
  }

  static void test()
  {
    conduit::Node host, copied;
    create(host);
    utils::copy<ExecSpace>(copied, host);

    if constexpr(axom::execution_space<ExecSpace>::onDevice())
    {
      ASSERT_TRUE(copied.has_child("_axom_single_elements"));
      const conduit::Node& storage = copied["_axom_single_elements"];
      ASSERT_TRUE(storage.has_child("bytes_1"));
      ASSERT_TRUE(storage.has_child("bytes_2"));
      ASSERT_TRUE(storage.has_child("bytes_4"));
      ASSERT_TRUE(storage.has_child("bytes_8"));

      EXPECT_TRUE(utils::isDeviceAllocated(storage["bytes_1"]));
      EXPECT_TRUE(utils::isDeviceAllocated(storage["bytes_2"]));
      EXPECT_TRUE(utils::isDeviceAllocated(storage["bytes_4"]));
      EXPECT_TRUE(utils::isDeviceAllocated(storage["bytes_8"]));

      EXPECT_TRUE(utils::isDeviceAllocated(copied["coordsets/coords/values/x"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["topologies/mesh/elements/connectivity"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["topologies/mesh/elements/sizes"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["topologies/mesh/elements/offsets"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["fields/one/values"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["fields/mixed/matset_values"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["fields/mixed_nested/matset_values/A/values"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["fields/flag/values"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["fields/small/values"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["matsets/mat/material_ids"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["matsets/mat/volume_fractions"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["matsets/mat/sizes"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["matsets/mat/offsets"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["matsets/mat/indices"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["matsets/mat_multi/volume_fractions/A/values"]));
      EXPECT_TRUE(utils::isDeviceAllocated(copied["matsets/mat_multi/volume_fractions/A/indices"]));
      EXPECT_FALSE(utils::isDeviceAllocated(copied["metadata/scalar_int"]));

      EXPECT_TRUE(copied["coordsets/coords/values/x"].is_data_external());
      EXPECT_TRUE(copied["topologies/mesh/elements/connectivity"].is_data_external());
      EXPECT_TRUE(copied["fields/one/values"].is_data_external());
      EXPECT_TRUE(copied["fields/mixed/matset_values"].is_data_external());
      EXPECT_TRUE(copied["fields/mixed_nested/matset_values/A/values"].is_data_external());
      EXPECT_TRUE(copied["fields/flag/values"].is_data_external());
      EXPECT_TRUE(copied["fields/small/values"].is_data_external());
      EXPECT_TRUE(copied["matsets/mat/material_ids"].is_data_external());
      EXPECT_TRUE(copied["matsets/mat/volume_fractions"].is_data_external());
      EXPECT_TRUE(copied["matsets/mat/sizes"].is_data_external());
      EXPECT_TRUE(copied["matsets/mat/offsets"].is_data_external());
      EXPECT_TRUE(copied["matsets/mat/indices"].is_data_external());
      EXPECT_TRUE(copied["matsets/mat_multi/volume_fractions/A/values"].is_data_external());
      EXPECT_TRUE(copied["matsets/mat_multi/volume_fractions/A/indices"].is_data_external());
      EXPECT_FALSE(copied["metadata/scalar_int"].is_data_external());

      expectPointerInBucket(copied["fields/flag/values"], storage["bytes_1"]);
      expectPointerInBucket(copied["fields/small/values"], storage["bytes_2"]);
      expectPointerInBucket(copied["coordsets/coords/values/x"], storage["bytes_4"]);
      expectPointerInBucket(copied["fields/one/values"], storage["bytes_8"]);
      expectPointerInBucket(copied["fields/mixed_nested/matset_values/A/values"], storage["bytes_4"]);
      expectPointerInBucket(copied["matsets/mat_multi/volume_fractions/A/values"],
                            storage["bytes_4"]);
      expectPointerInBucket(copied["matsets/mat_multi/volume_fractions/A/indices"],
                            storage["bytes_4"]);
    }
    else
    {
      EXPECT_FALSE(copied.has_child("_axom_single_elements"));
      EXPECT_FALSE(copied["coordsets/coords/values/x"].is_data_external());
      EXPECT_FALSE(copied["fields/one/values"].is_data_external());
      EXPECT_FALSE(copied["fields/flag/values"].is_data_external());
    }

    expectSameLeafType(host["coordsets/coords/values/x"], copied["coordsets/coords/values/x"]);
    expectSameLeafType(host["topologies/mesh/elements/connectivity"],
                       copied["topologies/mesh/elements/connectivity"]);
    expectSameLeafType(host["fields/one/values"], copied["fields/one/values"]);
    expectSameLeafType(host["fields/mixed/matset_values"], copied["fields/mixed/matset_values"]);
    expectSameLeafType(host["fields/mixed_nested/matset_values/A/values"],
                       copied["fields/mixed_nested/matset_values/A/values"]);
    expectSameLeafType(host["fields/flag/values"], copied["fields/flag/values"]);
    expectSameLeafType(host["fields/small/values"], copied["fields/small/values"]);
    expectSameLeafType(host["matsets/mat/material_ids"], copied["matsets/mat/material_ids"]);
    expectSameLeafType(host["matsets/mat/volume_fractions"], copied["matsets/mat/volume_fractions"]);
    expectSameLeafType(host["matsets/mat/sizes"], copied["matsets/mat/sizes"]);
    expectSameLeafType(host["matsets/mat/offsets"], copied["matsets/mat/offsets"]);
    expectSameLeafType(host["matsets/mat/indices"], copied["matsets/mat/indices"]);
    expectSameLeafType(host["matsets/mat_multi/volume_fractions/A/values"],
                       copied["matsets/mat_multi/volume_fractions/A/values"]);
    expectSameLeafType(host["matsets/mat_multi/volume_fractions/A/indices"],
                       copied["matsets/mat_multi/volume_fractions/A/indices"]);
    expectSameLeafType(host["metadata/scalar_int"], copied["metadata/scalar_int"]);

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<double> results(16, 16, allocatorID);
    auto resultsView = results.view();
    auto xView = utils::make_array_view<float>(copied["coordsets/coords/values/x"]);
    auto connView = utils::make_array_view<int>(copied["topologies/mesh/elements/connectivity"]);
    auto sizesView = utils::make_array_view<int>(copied["topologies/mesh/elements/sizes"]);
    auto offsetsView = utils::make_array_view<int>(copied["topologies/mesh/elements/offsets"]);
    auto fieldView = utils::make_array_view<double>(copied["fields/one/values"]);
    auto mixedFieldView = utils::make_array_view<float>(copied["fields/mixed/matset_values"]);
    auto nestedMixedFieldView =
      utils::make_array_view<float>(copied["fields/mixed_nested/matset_values/A/values"]);
    auto flagFieldView = utils::make_array_view<conduit::uint8>(copied["fields/flag/values"]);
    auto smallFieldView = utils::make_array_view<conduit::uint16>(copied["fields/small/values"]);
    auto materialIDsView = utils::make_array_view<int>(copied["matsets/mat/material_ids"]);
    auto volumeFractionsView =
      utils::make_array_view<float>(copied["matsets/mat/volume_fractions"]);
    auto matSizesView = utils::make_array_view<int>(copied["matsets/mat/sizes"]);
    auto matOffsetsView = utils::make_array_view<int>(copied["matsets/mat/offsets"]);
    auto matIndicesView = utils::make_array_view<int>(copied["matsets/mat/indices"]);
    auto matMultiValuesView =
      utils::make_array_view<float>(copied["matsets/mat_multi/volume_fractions/A/values"]);
    auto matMultiIndicesView =
      utils::make_array_view<int>(copied["matsets/mat_multi/volume_fractions/A/indices"]);

    axom::for_all<ExecSpace>(1, [=] AXOM_HOST_DEVICE(axom::IndexType) {
      resultsView[0] = xView[0];
      resultsView[1] = fieldView[0];
      resultsView[2] = mixedFieldView[0];
      resultsView[3] = nestedMixedFieldView[0];
      resultsView[4] = flagFieldView[0];
      resultsView[5] = smallFieldView[0];
      resultsView[6] = connView[0];
      resultsView[7] = sizesView[0];
      resultsView[8] = offsetsView[0];
      resultsView[9] = materialIDsView[0];
      resultsView[10] = volumeFractionsView[0];
      resultsView[11] = matSizesView[0];
      resultsView[12] = matOffsetsView[0];
      resultsView[13] = matIndicesView[0];
      resultsView[14] = matMultiValuesView[0];
      resultsView[15] = matMultiIndicesView[0];
    });

    double resultsHost[16] = {};
    axom::copy(resultsHost, results.data(), sizeof(resultsHost));
    EXPECT_DOUBLE_EQ(resultsHost[0], 1.5);
    EXPECT_DOUBLE_EQ(resultsHost[1], 2.5);
    EXPECT_DOUBLE_EQ(resultsHost[2], 3.5);
    EXPECT_DOUBLE_EQ(resultsHost[3], 4.5);
    EXPECT_DOUBLE_EQ(resultsHost[4], 9.);
    EXPECT_DOUBLE_EQ(resultsHost[5], 12.);
    EXPECT_DOUBLE_EQ(resultsHost[6], 0.);
    EXPECT_DOUBLE_EQ(resultsHost[7], 1.);
    EXPECT_DOUBLE_EQ(resultsHost[8], 0.);
    EXPECT_DOUBLE_EQ(resultsHost[9], 22.);
    EXPECT_DOUBLE_EQ(resultsHost[10], 1.);
    EXPECT_DOUBLE_EQ(resultsHost[11], 1.);
    EXPECT_DOUBLE_EQ(resultsHost[12], 0.);
    EXPECT_DOUBLE_EQ(resultsHost[13], 0.);
    EXPECT_DOUBLE_EQ(resultsHost[14], 0.75);
    EXPECT_DOUBLE_EQ(resultsHost[15], 0.);
  }
};

TEST(bump_utilities, copy_single_value_blueprint_leaves_seq)
{
  test_copy_single_value_blueprint_leaves<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, copy_single_value_blueprint_leaves_omp)
{
  test_copy_single_value_blueprint_leaves<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, copy_single_value_blueprint_leaves_cuda)
{
  test_copy_single_value_blueprint_leaves<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, copy_single_value_blueprint_leaves_hip)
{
  test_copy_single_value_blueprint_leaves<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_copy_single_value_blueprint_subtrees
{
  static void test()
  {
    conduit::Node hostCoordset;
    hostCoordset["type"] = "explicit";
    hostCoordset["values/x"].set(std::vector<float> {6.5f});

    conduit::Node hostTopology;
    hostTopology["type"] = "unstructured";
    hostTopology["coordset"] = "coords";
    hostTopology["elements/shape"] = "point";
    hostTopology["elements/connectivity"].set(std::vector<int> {0});
    hostTopology["elements/sizes"].set(std::vector<int> {1});
    hostTopology["elements/offsets"].set(std::vector<int> {0});

    conduit::Node hostField;
    hostField["topology"] = "mesh";
    hostField["association"] = "element";
    hostField["values"].set(std::vector<double> {7.5});

    conduit::Node hostMixedField;
    hostMixedField["topology"] = "mesh";
    hostMixedField["association"] = "element";
    hostMixedField["matset"] = "mat";
    hostMixedField["matset_values/A/values"].set(std::vector<float> {8.5f});

    conduit::Node hostMatset;
    hostMatset["topology"] = "mesh";
    hostMatset["volume_fractions/A/values"].set(std::vector<float> {1.f});
    hostMatset["volume_fractions/A/indices"].set(std::vector<int> {0});

    conduit::Node copiedCoordset, copiedTopology, copiedField, copiedMixedField, copiedMatset;
    utils::copy<ExecSpace>(copiedCoordset, hostCoordset);
    utils::copy<ExecSpace>(copiedTopology, hostTopology);
    utils::copy<ExecSpace>(copiedField, hostField);
    utils::copy<ExecSpace>(copiedMixedField, hostMixedField);
    utils::copy<ExecSpace>(copiedMatset, hostMatset);

    if constexpr(axom::execution_space<ExecSpace>::onDevice())
    {
      EXPECT_TRUE(copiedCoordset.has_child("_axom_single_elements"));
      EXPECT_TRUE(copiedTopology.has_child("_axom_single_elements"));
      EXPECT_TRUE(copiedField.has_child("_axom_single_elements"));
      EXPECT_TRUE(copiedMixedField.has_child("_axom_single_elements"));
      EXPECT_TRUE(copiedMatset.has_child("_axom_single_elements"));

      EXPECT_TRUE(copiedCoordset["values/x"].is_data_external());
      EXPECT_TRUE(copiedTopology["elements/connectivity"].is_data_external());
      EXPECT_TRUE(copiedTopology["elements/sizes"].is_data_external());
      EXPECT_TRUE(copiedTopology["elements/offsets"].is_data_external());
      EXPECT_TRUE(copiedField["values"].is_data_external());
      EXPECT_TRUE(copiedMixedField["matset_values/A/values"].is_data_external());
      EXPECT_TRUE(copiedMatset["volume_fractions/A/values"].is_data_external());
      EXPECT_TRUE(copiedMatset["volume_fractions/A/indices"].is_data_external());
    }
    else
    {
      EXPECT_FALSE(copiedCoordset.has_child("_axom_single_elements"));
      EXPECT_FALSE(copiedTopology.has_child("_axom_single_elements"));
      EXPECT_FALSE(copiedField.has_child("_axom_single_elements"));
      EXPECT_FALSE(copiedMixedField.has_child("_axom_single_elements"));
      EXPECT_FALSE(copiedMatset.has_child("_axom_single_elements"));
    }

    expectSameLeafType(hostCoordset["values/x"], copiedCoordset["values/x"]);
    expectSameLeafType(hostTopology["elements/connectivity"],
                       copiedTopology["elements/connectivity"]);
    expectSameLeafType(hostTopology["elements/sizes"], copiedTopology["elements/sizes"]);
    expectSameLeafType(hostTopology["elements/offsets"], copiedTopology["elements/offsets"]);
    expectSameLeafType(hostField["values"], copiedField["values"]);
    expectSameLeafType(hostMixedField["matset_values/A/values"],
                       copiedMixedField["matset_values/A/values"]);
    expectSameLeafType(hostMatset["volume_fractions/A/values"],
                       copiedMatset["volume_fractions/A/values"]);
    expectSameLeafType(hostMatset["volume_fractions/A/indices"],
                       copiedMatset["volume_fractions/A/indices"]);

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<double> results(8, 8, allocatorID);
    auto resultsView = results.view();
    auto xView = utils::make_array_view<float>(copiedCoordset["values/x"]);
    auto connView = utils::make_array_view<int>(copiedTopology["elements/connectivity"]);
    auto sizesView = utils::make_array_view<int>(copiedTopology["elements/sizes"]);
    auto offsetsView = utils::make_array_view<int>(copiedTopology["elements/offsets"]);
    auto fieldView = utils::make_array_view<double>(copiedField["values"]);
    auto mixedFieldView = utils::make_array_view<float>(copiedMixedField["matset_values/A/values"]);
    auto vfView = utils::make_array_view<float>(copiedMatset["volume_fractions/A/values"]);
    auto vfIndicesView = utils::make_array_view<int>(copiedMatset["volume_fractions/A/indices"]);

    axom::for_all<ExecSpace>(1, [=] AXOM_HOST_DEVICE(axom::IndexType) {
      resultsView[0] = xView[0];
      resultsView[1] = connView[0];
      resultsView[2] = sizesView[0];
      resultsView[3] = offsetsView[0];
      resultsView[4] = fieldView[0];
      resultsView[5] = mixedFieldView[0];
      resultsView[6] = vfView[0];
      resultsView[7] = vfIndicesView[0];
    });

    double resultsHost[8] = {};
    axom::copy(resultsHost, results.data(), sizeof(resultsHost));
    EXPECT_DOUBLE_EQ(resultsHost[0], 6.5);
    EXPECT_DOUBLE_EQ(resultsHost[1], 0.);
    EXPECT_DOUBLE_EQ(resultsHost[2], 1.);
    EXPECT_DOUBLE_EQ(resultsHost[3], 0.);
    EXPECT_DOUBLE_EQ(resultsHost[4], 7.5);
    EXPECT_DOUBLE_EQ(resultsHost[5], 8.5);
    EXPECT_DOUBLE_EQ(resultsHost[6], 1.);
    EXPECT_DOUBLE_EQ(resultsHost[7], 0.);
  }
};

TEST(bump_utilities, copy_single_value_blueprint_subtrees_seq)
{
  test_copy_single_value_blueprint_subtrees<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_utilities, copy_single_value_blueprint_subtrees_omp)
{
  test_copy_single_value_blueprint_subtrees<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_utilities, copy_single_value_blueprint_subtrees_cuda)
{
  test_copy_single_value_blueprint_subtrees<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_utilities, copy_single_value_blueprint_subtrees_hip)
{
  test_copy_single_value_blueprint_subtrees<hip_exec>::test();
}
#endif

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
