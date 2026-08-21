// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_set_BivariateSet_device.cpp
 *
 * \brief Device-execution tests for the BivariateSet std::optional-returning queries
 *        and for std::optional itself inside device kernels.
 *
 * These accompany the host-side tests in slam_set_BivariateSet.cpp. 
 * The `AXOM_HOST_DEVICE findElementFlatIndexOptional(...)` wrappers construct a
 * `std::optional` in device code, so this file exercises two things on device:
 *
 *  1. that `std::optional<PositionType>` is constructible and queryable inside a RAJA kernel
 *  2. that `findElementFlatIndexOptional(...)` returns the correct value when evaluated in a kernel.
 *
 * The wrapper is only device-callable on the concrete (non-virtual) set instantiation (`ProductSet::ConcreteSet`).
 * The virtual `BivariateSet` interface would require a device-resident vtable, 
 * whereas on a concrete set captured by value the wrapper inlines
 * and its internal `findElementFlatIndex` call devirtualizes -- no vtable is touched.
 */

#include "gtest/gtest.h"

#include "axom/core/execution/execution_space.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <optional>

namespace
{
namespace slam = axom::slam;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

template <typename ExecSpace>
int getHostAccessibleAllocatorId()
{
#ifdef AXOM_USE_UMPIRE
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    return axom::detail::getAllocatorID<axom::MemorySpace::Unified>();
  }
#endif

  return axom::execution_space<ExecSpace>::allocatorID();
}

//------------------------------------------------------------------------------
template <typename ExecutionSpace>
class slam_set_bivariate_optional_device : public ::testing::Test
{
public:
  using ExecSpace = ExecutionSpace;

  // The device-usable (non-virtual) set instantiations.
  using ConcreteSetType = typename slam::RangeSet<SetPosition, SetElement>::ConcreteSet;
  using ProductSetType = typename slam::ProductSet<ConcreteSetType, ConcreteSetType>::ConcreteSet;

  // This translation unit is compiled by each enabled backend.
  // Keep these representative concept instantiations beside the device-captured types.
  static_assert(slam::SetLike<ConcreteSetType>);
  static_assert(slam::BivariateSetLike<ProductSetType>);
  static_assert(slam::DeviceCapturable<ConcreteSetType>);
  static_assert(slam::DeviceCapturable<ProductSetType>);
};

using MyTypes = ::testing::Types<
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  axom::OMP_EXEC,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  axom::CUDA_EXEC<256>,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  axom::HIP_EXEC<256>,
#endif
  axom::SEQ_EXEC>;

TYPED_TEST_SUITE(slam_set_bivariate_optional_device, MyTypes);

//------------------------------------------------------------------------------
// std::optional<PositionType> is constructible and queryable in a device kernel.
//------------------------------------------------------------------------------
AXOM_TYPED_TEST(slam_set_bivariate_optional_device, std_optional_in_kernel)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  const int allocatorId = getHostAccessibleAllocatorId<ExecSpace>();

  constexpr int N = 8;
  constexpr axom::IndexType SENTINEL = -1;

  axom::Array<axom::IndexType> results(N, N, allocatorId);
  axom::Array<int> consistent(N, N, allocatorId);
  auto results_v = results.view();
  auto consistent_v = consistent.view();

  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(int i) {
      std::optional<axom::IndexType> opt;
      if(i % 2 == 0)
      {
        opt = std::optional<axom::IndexType>(static_cast<axom::IndexType>(i * 10));
      }

      const bool has = opt.has_value();
      const bool asBool = static_cast<bool>(opt);
      const axom::IndexType viaValueOr = opt.value_or(SENTINEL);

      axom::IndexType decoded = SENTINEL;
      if(has)
      {
        decoded = *opt;
      }
      results_v[i] = decoded;

      const bool ok = (has == asBool) && (has ? (viaValueOr == decoded) : (viaValueOr == SENTINEL));
      consistent_v[i] = ok ? 1 : 0;
    });

  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(consistent[i], 1) << "std::optional query surface inconsistent on device at i=" << i;
    const axom::IndexType expected = (i % 2 == 0) ? static_cast<axom::IndexType>(i * 10) : SENTINEL;
    EXPECT_EQ(results[i], expected);
  }
}

//------------------------------------------------------------------------------
// findElementFlatIndexOptional on a concrete ProductSet, evaluated in a kernel.
//
// A ProductSet is dense: every (i, j) exists, so the returned optional must be
// engaged with FlatIndex == secondSetSize * i + j. The sets are placed in unified
// memory and the (pointer-holding) ProductSet is captured by value into the
// kernel, matching the BivariateMap device tests.
//------------------------------------------------------------------------------
AXOM_TYPED_TEST(slam_set_bivariate_optional_device, product_set_flat_index_optional_in_kernel)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using ConcreteSetType = typename TestFixture::ConcreteSetType;
  using ProductSetType = typename TestFixture::ProductSetType;

  const int allocatorId = getHostAccessibleAllocatorId<ExecSpace>();

  constexpr SetPosition SZ1 = 4;
  constexpr SetPosition SZ2 = 5;

  // Sets must be reachable from device code.
  // Place them in unified memory and point the ProductSet at them.
  axom::Array<ConcreteSetType> sets(2, 2, allocatorId);
  sets[0] = ConcreteSetType(SZ1);
  sets[1] = ConcreteSetType(SZ2);

  ProductSetType prodSet(&sets[0], &sets[1]);
  EXPECT_EQ(prodSet.size(), SZ1 * SZ2);
  EXPECT_TRUE(prodSet.isValid());

  const int totalSize = static_cast<int>(prodSet.size());
  axom::Array<int> ok(totalSize, totalSize, allocatorId);
  auto ok_v = ok.view();

  axom::for_all<ExecSpace>(
    prodSet.firstSetSize(),
    AXOM_LAMBDA(int i) {
      const auto sz2 = prodSet.secondSetSize();
      for(int j = 0; j < sz2; ++j)
      {
        const std::optional<SetPosition> flat = prodSet.findElementFlatIndexOptional(i, j);
        const SetPosition expected = sz2 * i + j;
        const bool good = flat.has_value() && (*flat == expected);
        ok_v[i * sz2 + j] = good ? 1 : 0;
      }
    });

  for(int idx = 0; idx < totalSize; ++idx)
  {
    EXPECT_EQ(ok[idx], 1) << "findElementFlatIndexOptional mismatch on device at flat index " << idx;
  }
}

}  // end anonymous namespace

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
