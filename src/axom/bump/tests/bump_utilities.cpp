// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit_blueprint_mesh_examples.hpp>

#include <iostream>
#include <algorithm>

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_conduit_allocate
{
  static void test()
  {
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    EXPECT_TRUE(c2a.getConduitAllocatorID() > 0);

    constexpr int nValues = 100;
    conduit::Node n;
    n.set_allocator(c2a.getConduitAllocatorID());
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
    n.set_allocator(c2a.getConduitAllocatorID());
    n.set(conduit::DataType::int32(0));
    EXPECT_EQ(n.dtype().number_of_elements(), 0);
  }
};

TEST(bump_blueprint_utilities, allocate_seq) { test_conduit_allocate<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_blueprint_utilities, allocate_omp) { test_conduit_allocate<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_blueprint_utilities, allocate_cuda) { test_conduit_allocate<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_blueprint_utilities, allocate_hip) { test_conduit_allocate<hip_exec>::test(); }
#endif

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

  static void create(conduit::Node &mesh)
  {
    const int d[3] = {10, 10, 10};
    conduit::blueprint::mesh::examples::braid("hexs", d[0], d[1], d[2], mesh);
  }
};

TEST(bump_blueprint_utilities, copy_seq) { test_copy_braid<seq_exec>::test(); }

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_blueprint_utilities, copy_omp) { test_copy_braid<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_blueprint_utilities, copy_cuda) { test_copy_braid<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_blueprint_utilities, copy_hip) { test_copy_braid<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
