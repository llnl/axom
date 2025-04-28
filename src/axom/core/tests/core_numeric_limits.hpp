// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"  // for compile time definitions

#include "axom/core/NumericLimits.hpp"

// for gtest macros
#include "gtest/gtest.h"

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(core_NumericLimits, check_CPU)
{
  //
  // Tests to compare axom::numeric_limits to std::numeric_limits
  // to ensure that Axom type aliasing is correct.
  //
  EXPECT_TRUE(axom::numeric_limits<int>::lowest() == std::numeric_limits<int>::lowest());
  EXPECT_TRUE(axom::numeric_limits<int>::min() == std::numeric_limits<int>::min());
  EXPECT_TRUE(axom::numeric_limits<int>::max() == std::numeric_limits<int>::max());
  EXPECT_TRUE(axom::numeric_limits<int>::is_signed == std::numeric_limits<int>::is_signed);

  EXPECT_TRUE(axom::numeric_limits<float>::lowest() == std::numeric_limits<float>::lowest());
  EXPECT_TRUE(axom::numeric_limits<float>::min() == std::numeric_limits<float>::min());
  EXPECT_TRUE(axom::numeric_limits<float>::max() == std::numeric_limits<float>::max());

  EXPECT_TRUE(axom::numeric_limits<double>::lowest() == std::numeric_limits<double>::lowest());
  EXPECT_TRUE(axom::numeric_limits<double>::min() == std::numeric_limits<double>::min());
  EXPECT_TRUE(axom::numeric_limits<double>::max() == std::numeric_limits<double>::max());
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_CUDA)
//
// Tests to ensure axom::numeric_limits type alias does the correct thing
// in host and CUDA device code.
//

//
// Simple device kernel
//
__global__ void cuda_kernel(int* a, size_t* b, float* c, double* d)
{
  a[0] = axom::numeric_limits<int>::min();
  b[0] = axom::numeric_limits<size_t>::max();
  c[0] = axom::numeric_limits<float>::lowest();
  d[0] = axom::numeric_limits<double>::max();
}

TEST(core_NumericLimits, check_CUDA)
{
  //
  // Device memory allocation and initialiation for a few different types.
  //
  int* a;
  (void)cudaMalloc(&a, sizeof(int));
  (void)cudaMemset(a, 0, sizeof(int));

  size_t* b;
  (void)cudaMalloc(&b, sizeof(size_t));
  (void)cudaMemset(b, 0, sizeof(size_t));

  float* c;
  (void)cudaMalloc(&c, sizeof(float));
  (void)cudaMemset(c, 0, sizeof(float));

  double* d;
  (void)cudaMalloc(&d, sizeof(double));
  (void)cudaMemset(d, 0, sizeof(double));

  //
  // Set values in device code.
  //
  cuda_kernel<<<1, 1>>>(a, b, c, d);

  //
  // Copy device values back to host and compare with expectations....
  //
  int ha;
  size_t hb;
  float hc;
  double hd;
  (void)cudaMemcpy(&ha, a, sizeof(int), cudaMemcpyDeviceToHost);
  (void)cudaMemcpy(&hb, b, sizeof(size_t), cudaMemcpyDeviceToHost);
  (void)cudaMemcpy(&hc, c, sizeof(float), cudaMemcpyDeviceToHost);
  (void)cudaMemcpy(&hd, d, sizeof(double), cudaMemcpyDeviceToHost);

  EXPECT_TRUE(ha == axom::numeric_limits<int>::min());
  EXPECT_TRUE(hb == axom::numeric_limits<size_t>::max());
  EXPECT_TRUE(hc == axom::numeric_limits<float>::lowest());
  EXPECT_TRUE(hd == axom::numeric_limits<double>::max());
}
#endif

//------------------------------------------------------------------------------
#if defined(AXOM_USE_HIP)
//
// Tests to ensure axom::numeric_limits type alias does the correct thing
// in host and CUDA device code.
//

//
// Simple device kernel
//
__global__ void hip_kernel(int* a, size_t* b, float* c, double* d)
{
  a[0] = axom::numeric_limits<int>::min();
  b[0] = axom::numeric_limits<size_t>::max();
  c[0] = axom::numeric_limits<float>::lowest();
  d[0] = axom::numeric_limits<double>::max();
}

TEST(core_NumericLimits, check_HIP)
{
  //
  // Device memory allocation and initialiation for a few different types.
  //
  int* a;
  (void)hipMalloc(&a, sizeof(int));
  (void)hipMemset(a, 0, sizeof(int));

  size_t* b;
  (void)hipMalloc(&b, sizeof(size_t));
  (void)hipMemset(b, 0, sizeof(size_t));

  float* c;
  (void)hipMalloc(&c, sizeof(float));
  (void)hipMemset(c, 0, sizeof(float));

  double* d;
  (void)hipMalloc(&d, sizeof(double));
  (void)hipMemset(d, 0, sizeof(double));

  //
  // Set values in device code.
  //
  hip_kernel<<<1, 1>>>(a, b, c, d);

  //
  // Copy device values back to host and compare with expectations....
  //
  int ha;
  size_t hb;
  float hc;
  double hd;
  (void)hipMemcpy(&ha, a, sizeof(int), hipMemcpyDeviceToHost);
  (void)hipMemcpy(&hb, b, sizeof(size_t), hipMemcpyDeviceToHost);
  (void)hipMemcpy(&hc, c, sizeof(float), hipMemcpyDeviceToHost);
  (void)hipMemcpy(&hd, d, sizeof(double), hipMemcpyDeviceToHost);

  EXPECT_TRUE(ha == axom::numeric_limits<int>::min());
  EXPECT_TRUE(hb == axom::numeric_limits<size_t>::max());
  EXPECT_TRUE(hc == axom::numeric_limits<float>::lowest());
  EXPECT_TRUE(hd == axom::numeric_limits<double>::max());
}
#endif
