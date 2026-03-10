// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/execution/scans.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/Types.hpp"

#include "gtest/gtest.h"

#include <vector>
#include <type_traits>

using axom::IndexType;

// -----------------------------------------------------------------------------
// Reference helpers on host
// -----------------------------------------------------------------------------

template <typename InContainer, typename OutValue>
std::vector<OutValue> reference_exclusive_scan(const InContainer &input)
{
  std::vector<OutValue> result(input.size());
  OutValue total = 0;
  for(IndexType i = 0; i < static_cast<IndexType>(input.size()); ++i)
  {
    result[i] = total;
    total += static_cast<OutValue>(input[i]);
  }
  return result;
}

template <typename InContainer, typename OutValue>
std::vector<OutValue> reference_inclusive_scan(const InContainer &input)
{
  std::vector<OutValue> result(input.size());
  OutValue total = 0;
  for(IndexType i = 0; i < static_cast<IndexType>(input.size()); ++i)
  {
    total += static_cast<OutValue>(input[i]);
    result[i] = total;
  }
  return result;
}

// -----------------------------------------------------------------------------
// Data generation on host
// -----------------------------------------------------------------------------

// Fill a mask vector with 0/1 in the input type.
// alternating = true   -> 0,1,0,1, ...
// alternating = false  -> 1 every "stride" elements.
template <typename T>
std::vector<T> make_mask(IndexType n, bool alternating, int stride = 1)
{
  std::vector<T> mask(n);
  for(IndexType i = 0; i < n; ++i)
  {
    int val;
    if(alternating)
    {
      val = (i % 2 == 0) ? 0 : 1;
    }
    else
    {
      val = ((i % stride) == 0) ? 1 : 0;
    }
    mask[i] = static_cast<T>(val);
  }
  return mask;
}

// -----------------------------------------------------------------------------
// Common utilities for device/exec-space arrays
// -----------------------------------------------------------------------------

// Create an axom::Array<T> in the allocator for ExecSpace and fill it from host
template <typename ExecSpace, typename T>
axom::Array<T> create_exec_array_from_host(const std::vector<T> &hostData)
{
  const IndexType n = static_cast<IndexType>(hostData.size());
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  axom::Array<T> arr(n, n, allocatorID);

  if(n > 0)
  {
    axom::copy(arr.data(), hostData.data(), n * sizeof(T));
  }

  return arr;
}

// Copy an axom::Array<T> in exec-space back to a host std::vector<T>
template <typename T>
std::vector<T> copy_exec_array_to_host(const axom::Array<T> &arr)
{
  const IndexType n = arr.size();
  std::vector<T> hostData(n);

  if(n > 0)
  {
    axom::copy(hostData.data(), arr.data(), n * sizeof(T));
  }

  return hostData;
}

// -----------------------------------------------------------------------------
// Generic test runners: host <-> exec arrays, all type combinations
// -----------------------------------------------------------------------------

// exclusive_scan: InType host mask, OutType host result; using ExecSpace
template <typename ExecSpace, typename InType, typename OutType>
void run_exclusive_scan_test(IndexType n,
                             bool alternating_pattern,
                             int stride = 1)
{
  // 1. Host input
  auto hostMask = make_mask<InType>(n, alternating_pattern, stride);

  // 2. Exec-space Arrays (input + output)
  axom::Array<InType> devInput =
    create_exec_array_from_host<ExecSpace>(hostMask);

  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<OutType> devOutput(n, n, allocatorID);

  // 3. Run Axom scan in ExecSpace
  axom::exclusive_scan<ExecSpace>(devInput, devOutput);

  // 4. Copy output back to host and compute reference
  auto hostOutput = copy_exec_array_to_host(devOutput);
  auto hostRef = reference_exclusive_scan<std::vector<InType>, OutType>(hostMask);

  ASSERT_EQ(hostOutput.size(), hostRef.size());
  for(IndexType i = 0; i < n; ++i)
  {
    EXPECT_EQ(hostOutput[i], hostRef[i]) << "Mismatch at index " << i;
  }

  // 5. Check total sum is positive (since pattern is 0/1 and n > 0 for main tests)
  if(n > 0)
  {
    OutType total_axom;
    if(n == 1)
    {
      total_axom = static_cast<OutType>(hostMask[0]);
    }
    else
    {
      total_axom = hostOutput.back()
                   + static_cast<OutType>(hostMask.back());
    }

    OutType total_ref = 0;
    for(auto v : hostMask)
    {
      total_ref += static_cast<OutType>(v);
    }
    EXPECT_EQ(total_axom, total_ref);
    EXPECT_GT(total_axom, static_cast<OutType>(0));
  }
}

// inclusive_scan: InType host mask, OutType host result; using ExecSpace
template <typename ExecSpace, typename InType, typename OutType>
void run_inclusive_scan_test(IndexType n,
                             bool alternating_pattern,
                             int stride = 1)
{
  auto hostMask = make_mask<InType>(n, alternating_pattern, stride);

  axom::Array<InType> devInput =
    create_exec_array_from_host<ExecSpace>(hostMask);

  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<OutType> devOutput(n, n, allocatorID);

  axom::inclusive_scan<ExecSpace>(devInput, devOutput);

  auto hostOutput = copy_exec_array_to_host(devOutput);
  auto hostRef = reference_inclusive_scan<std::vector<InType>, OutType>(hostMask);

  ASSERT_EQ(hostOutput.size(), hostRef.size());
  for(IndexType i = 0; i < n; ++i)
  {
    EXPECT_EQ(hostOutput[i], hostRef[i]) << "Mismatch at index " << i;
  }

  if(n > 0)
  {
    OutType total_axom = hostOutput.back();
    OutType total_ref = 0;
    for(auto v : hostMask)
    {
      total_ref += static_cast<OutType>(v);
    }
    EXPECT_EQ(total_axom, total_ref);
    EXPECT_GT(total_axom, static_cast<OutType>(0));
  }
}

// exclusive_scan_inplace: InType host mask, OutType buffer in exec-space
// initialized from mask.
template <typename ExecSpace, typename InType, typename OutType>
void run_exclusive_scan_inplace_test(IndexType n,
                                     bool alternating_pattern,
                                     int stride = 1)
{
  auto hostMask = make_mask<InType>(n, alternating_pattern, stride);

  // Host buffer cast to OutType
  std::vector<OutType> hostBuf(n);
  for(IndexType i = 0; i < n; ++i)
  {
    hostBuf[i] = static_cast<OutType>(hostMask[i]);
  }

  axom::Array<OutType> devBuf =
    create_exec_array_from_host<ExecSpace>(hostBuf);

  axom::exclusive_scan_inplace<ExecSpace>(devBuf);

  auto hostOutput = copy_exec_array_to_host(devBuf);
  auto hostRef = reference_exclusive_scan<std::vector<InType>, OutType>(hostMask);

  ASSERT_EQ(hostOutput.size(), hostRef.size());
  for(IndexType i = 0; i < n; ++i)
  {
    EXPECT_EQ(hostOutput[i], hostRef[i]) << "Mismatch at index " << i;
  }

  if(n > 0)
  {
    OutType total_axom;
    if(n == 1)
    {
      total_axom = static_cast<OutType>(hostMask[0]);
    }
    else
    {
      total_axom = hostOutput.back()
                   + static_cast<OutType>(hostMask.back());
    }

    OutType total_ref = 0;
    for(auto v : hostMask)
    {
      total_ref += static_cast<OutType>(v);
    }
    EXPECT_EQ(total_axom, total_ref);
    EXPECT_GT(total_axom, static_cast<OutType>(0));
  }
}

// inclusive_scan_inplace: InType host mask, OutType buffer in exec-space
// initialized from mask.
template <typename ExecSpace, typename InType, typename OutType>
void run_inclusive_scan_inplace_test(IndexType n,
                                     bool alternating_pattern,
                                     int stride = 1)
{
  auto hostMask = make_mask<InType>(n, alternating_pattern, stride);

  std::vector<OutType> hostBuf(n);
  for(IndexType i = 0; i < n; ++i)
  {
    hostBuf[i] = static_cast<OutType>(hostMask[i]);
  }

  axom::Array<OutType> devBuf =
    create_exec_array_from_host<ExecSpace>(hostBuf);

  axom::inclusive_scan_inplace<ExecSpace>(devBuf);

  auto hostOutput = copy_exec_array_to_host(devBuf);
  auto hostRef = reference_inclusive_scan<std::vector<InType>, OutType>(hostMask);

  ASSERT_EQ(hostOutput.size(), hostRef.size());
  for(IndexType i = 0; i < n; ++i)
  {
    EXPECT_EQ(hostOutput[i], hostRef[i]) << "Mismatch at index " << i;
  }

  if(n > 0)
  {
    OutType total_axom = hostOutput.back();
    OutType total_ref = 0;
    for(auto v : hostMask)
    {
      total_ref += static_cast<OutType>(v);
    }
    EXPECT_EQ(total_axom, total_ref);
    EXPECT_GT(total_axom, static_cast<OutType>(0));
  }
}

// -----------------------------------------------------------------------------
// Type lists and execution space lists
// -----------------------------------------------------------------------------

// Input/output scalar types to test
using ScalarTypes = ::testing::Types<
  char,
  int,
  axom::IndexType>;

// Execution spaces to test (conditional on build)
using ExecSpaceTypes = ::testing::Types<
  axom::SEQ_EXEC
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  , axom::OMP_EXEC
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  , axom::CUDA_EXEC<256>
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  , axom::HIP_EXEC<256>
#endif
  >;

// A small helper struct to bundle ExecSpace and one scalar type
template <typename ExecSpace, typename Scalar>
struct ExecScalarPair
{
  using exec_space = ExecSpace;
  using scalar_type = Scalar;
};

// Build a combined list of (ExecSpace, Scalar) so we can iterate over
// scalar type combinations inside each test.
template <typename ExecSpace>
struct ExecScalarList
{
  using type =
    ::testing::Types<
      ExecScalarPair<ExecSpace, char>,
      ExecScalarPair<ExecSpace, int>,
      ExecScalarPair<ExecSpace, axom::IndexType>>;
};

// However, gtest does not directly support nested type lists nicely, so we
// instead parameterize only on ExecSpace and loop over scalar type combinations
// explicitly in each test.

// -----------------------------------------------------------------------------
// Test fixture parameterized by ExecSpace only
// -----------------------------------------------------------------------------

template <typename ExecSpace>
class ExecutionScanAllTypes : public ::testing::Test
{ };

TYPED_TEST_SUITE(ExecutionScanAllTypes, ExecSpaceTypes);

// -----------------------------------------------------------------------------
// Tests: for each ExecSpace, loop over all InType/OutType combinations
// among char, int, IndexType, and run all four scan variants.
// -----------------------------------------------------------------------------

TYPED_TEST(ExecutionScanAllTypes, ExclusiveScan_AllTypeCombos_128_Alternating)
{
  using ExecSpace = TypeParam;
  const IndexType N = 128;

  using InChar = char;
  using InInt = int;
  using InIdx = axom::IndexType;

  using OutChar = char;
  using OutInt = int;
  using OutIdx = axom::IndexType;

  // 3 x 3 input/output combinations
  run_exclusive_scan_test<ExecSpace, InChar, OutChar>(N, true);
  run_exclusive_scan_test<ExecSpace, InChar, OutInt>(N, true);
  run_exclusive_scan_test<ExecSpace, InChar, OutIdx>(N, true);

  run_exclusive_scan_test<ExecSpace, InInt, OutChar>(N, true);
  run_exclusive_scan_test<ExecSpace, InInt, OutInt>(N, true);
  run_exclusive_scan_test<ExecSpace, InInt, OutIdx>(N, true);

  run_exclusive_scan_test<ExecSpace, InIdx, OutChar>(N, true);
  run_exclusive_scan_test<ExecSpace, InIdx, OutInt>(N, true);
  run_exclusive_scan_test<ExecSpace, InIdx, OutIdx>(N, true);
}

TYPED_TEST(ExecutionScanAllTypes, InclusiveScan_AllTypeCombos_128_Alternating)
{
  using ExecSpace = TypeParam;
  const IndexType N = 128;

  using InChar = char;
  using InInt = int;
  using InIdx = axom::IndexType;

  using OutChar = char;
  using OutInt = int;
  using OutIdx = axom::IndexType;

  run_inclusive_scan_test<ExecSpace, InChar, OutChar>(N, true);
  run_inclusive_scan_test<ExecSpace, InChar, OutInt>(N, true);
  run_inclusive_scan_test<ExecSpace, InChar, OutIdx>(N, true);

  run_inclusive_scan_test<ExecSpace, InInt, OutChar>(N, true);
  run_inclusive_scan_test<ExecSpace, InInt, OutInt>(N, true);
  run_inclusive_scan_test<ExecSpace, InInt, OutIdx>(N, true);

  run_inclusive_scan_test<ExecSpace, InIdx, OutChar>(N, true);
  run_inclusive_scan_test<ExecSpace, InIdx, OutInt>(N, true);
  run_inclusive_scan_test<ExecSpace, InIdx, OutIdx>(N, true);
}

TYPED_TEST(ExecutionScanAllTypes, ExclusiveScanInplace_AllTypeCombos_128_Alternating)
{
  using ExecSpace = TypeParam;
  const IndexType N = 128;

  using InChar = char;
  using InInt = int;
  using InIdx = axom::IndexType;

  using OutChar = char;
  using OutInt = int;
  using OutIdx = axom::IndexType;

  run_exclusive_scan_inplace_test<ExecSpace, InChar, OutChar>(N, true);
  run_exclusive_scan_inplace_test<ExecSpace, InChar, OutInt>(N, true);
  run_exclusive_scan_inplace_test<ExecSpace, InChar, OutIdx>(N, true);

  run_exclusive_scan_inplace_test<ExecSpace, InInt, OutChar>(N, true);
  run_exclusive_scan_inplace_test<ExecSpace, InInt, OutInt>(N, true);
  run_exclusive_scan_inplace_test<ExecSpace, InInt, OutIdx>(N, true);

  run_exclusive_scan_inplace_test<ExecSpace, InIdx, OutChar>(N, true);
  run_exclusive_scan_inplace_test<ExecSpace, InIdx, OutInt>(N, true);
  run_exclusive_scan_inplace_test<ExecSpace, InIdx, OutIdx>(N, true);
}

TYPED_TEST(ExecutionScanAllTypes, InclusiveScanInplace_AllTypeCombos_128_Alternating)
{
  using ExecSpace = TypeParam;
  const IndexType N = 128;

  using InChar = char;
  using InInt = int;
  using InIdx = axom::IndexType;

  using OutChar = char;
  using OutInt = int;
  using OutIdx = axom::IndexType;

  run_inclusive_scan_inplace_test<ExecSpace, InChar, OutChar>(N, true);
  run_inclusive_scan_inplace_test<ExecSpace, InChar, OutInt>(N, true);
  run_inclusive_scan_inplace_test<ExecSpace, InChar, OutIdx>(N, true);

  run_inclusive_scan_inplace_test<ExecSpace, InInt, OutChar>(N, true);
  run_inclusive_scan_inplace_test<ExecSpace, InInt, OutInt>(N, true);
  run_inclusive_scan_inplace_test<ExecSpace, InInt, OutIdx>(N, true);

  run_inclusive_scan_inplace_test<ExecSpace, InIdx, OutChar>(N, true);
  run_inclusive_scan_inplace_test<ExecSpace, InIdx, OutInt>(N, true);
  run_inclusive_scan_inplace_test<ExecSpace, InIdx, OutIdx>(N, true);
}

// Sparse pattern, non power of two
TYPED_TEST(ExecutionScanAllTypes, ExclusiveScan_AllTypeCombos_257_Sparse)
{
  using ExecSpace = TypeParam;
  const IndexType N = 257;
  const int stride = 5;

  using InChar = char;
  using InInt = int;
  using InIdx = axom::IndexType;

  using OutChar = char;
  using OutInt = int;
  using OutIdx = axom::IndexType;

  run_exclusive_scan_test<ExecSpace, InChar, OutChar>(N, false, stride);
  run_exclusive_scan_test<ExecSpace, InChar, OutInt>(N, false, stride);
  run_exclusive_scan_test<ExecSpace, InChar, OutIdx>(N, false, stride);

  run_exclusive_scan_test<ExecSpace, InInt, OutChar>(N, false, stride);
  run_exclusive_scan_test<ExecSpace, InInt, OutInt>(N, false, stride);
  run_exclusive_scan_test<ExecSpace, InInt, OutIdx>(N, false, stride);

  run_exclusive_scan_test<ExecSpace, InIdx, OutChar>(N, false, stride);
  run_exclusive_scan_test<ExecSpace, InIdx, OutInt>(N, false, stride);
  run_exclusive_scan_test<ExecSpace, InIdx, OutIdx>(N, false, stride);
}

TYPED_TEST(ExecutionScanAllTypes, InclusiveScan_AllTypeCombos_257_Sparse)
{
  using ExecSpace = TypeParam;
  const IndexType N = 257;
  const int stride = 5;

  using InChar = char;
  using InInt = int;
  using InIdx = axom::IndexType;

  using OutChar = char;
  using OutInt = int;
  using OutIdx = axom::IndexType;

  run_inclusive_scan_test<ExecSpace, InChar, OutChar>(N, false, stride);
  run_inclusive_scan_test<ExecSpace, InChar, OutInt>(N, false, stride);
  run_inclusive_scan_test<ExecSpace, InChar, OutIdx>(N, false, stride);

  run_inclusive_scan_test<ExecSpace, InInt, OutChar>(N, false, stride);
  run_inclusive_scan_test<ExecSpace, InInt, OutInt>(N, false, stride);
  run_inclusive_scan_test<ExecSpace, InInt, OutIdx>(N, false, stride);

  run_inclusive_scan_test<ExecSpace, InIdx, OutChar>(N, false, stride);
  run_inclusive_scan_test<ExecSpace, InIdx, OutInt>(N, false, stride);
  run_inclusive_scan_test<ExecSpace, InIdx, OutIdx>(N, false, stride);
}

// Simple edge cases: empty and single element for one representative type combo
TYPED_TEST(ExecutionScanAllTypes, ExclusiveScan_Empty_IntToInt)
{
  using ExecSpace = TypeParam;

  std::vector<int> hostMask;
  axom::Array<int> devInput =
    create_exec_array_from_host<ExecSpace>(hostMask);

  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<int> devOutput(0, 0, allocatorID);

  axom::exclusive_scan<ExecSpace>(devInput, devOutput);

  auto hostOutput = copy_exec_array_to_host(devOutput);

  EXPECT_TRUE(hostMask.empty());
  EXPECT_TRUE(hostOutput.empty());
}

TYPED_TEST(ExecutionScanAllTypes, ExclusiveScan_SingleElement_CharToInt)
{
  using ExecSpace = TypeParam;

  std::vector<char> hostMask(1, 1);
  axom::Array<char> devInput =
    create_exec_array_from_host<ExecSpace>(hostMask);

  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<int> devOutput(1, 1, allocatorID);

  axom::exclusive_scan<ExecSpace>(devInput, devOutput);

  auto hostOutput = copy_exec_array_to_host(devOutput);

  ASSERT_EQ(hostOutput.size(), 1u);
  EXPECT_EQ(hostOutput[0], 0);
}

// If you do not have a separate test driver, you can add:
//
// int main(int argc, char** argv)
// {
//   ::testing::InitGoogleTest(&argc, argv);
//   return RUN_ALL_TESTS();
// }
