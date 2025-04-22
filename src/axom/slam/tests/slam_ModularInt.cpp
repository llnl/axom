// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_ModularInt.cpp
 *
 * Unit tests for the modular arithmetic class ModularInt
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/ModularInt.hpp"

namespace slam = axom::slam;

namespace
{
template <int N>
using IntSize = slam::policies::CompileTimeSize<int, N>;

using RTIntSize = slam::policies::RuntimeSize<int>;
}  // namespace

TEST(slam_modInt, runtime_modular_int_unitialized_and_full)
{
  using ModularIntType = slam::ModularInt<RTIntSize>;

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  SLIC_INFO("Checking that modular int with modulus zero fails.");
  SLIC_INFO("Note: Expecting a SLIC Failure: ");

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntType(0, 0), "")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntType(1, 0), "")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntType(), "") << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntType(1), "") << " SIZE of Modular int not allowed to be zero";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("Checking modular int with value set to "
            << " modulus equals (i.e. is equivalent to) 0 (runtime)");
  volatile int sz = 5;
  for(int i = 1; i < sz; ++i)
  {
    ModularIntType modIntFull(i, i);
    EXPECT_EQ(modIntFull, 0);
  }
}

TEST(slam_modInt, compile_modular_int_unitialized_and_full)
{
#ifdef AXOM_DEBUG
  using ModularIntZero = slam::ModularInt<IntSize<0>>;

  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  SLIC_INFO("Checking that modular int with modulus "
            << " zero fails.\nNote: Expecting a SLIC Failure: ");

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntZero(0, 0), "")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntZero(1, 0), "")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntZero(), "") << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(ModularIntZero(1), "") << " SIZE of Modular int not allowed to be zero";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("Checking modular int with value set to modulus"
            << " equals (i.e. is equivalent to) 0 for (compile time)");
  slam::ModularInt<IntSize<1>> m1(1);
  EXPECT_EQ(m1, 0);

  slam::ModularInt<IntSize<2>> m2(2);
  EXPECT_EQ(m2, 0);

  slam::ModularInt<IntSize<3>> m3(3);
  EXPECT_EQ(m3, 0);

  slam::ModularInt<IntSize<4>> m4(4);
  EXPECT_EQ(m4, 0);
}

TEST(slam_modInt, runtime_modular_int)
{
  SLIC_INFO("Checking modular int addition and subtraction"
            << " when supplying the max value at runtime");

  using ModularIntType = slam::ModularInt<RTIntSize>;

  volatile int sz = 937;

  for(int i = 0; i < sz; ++i)
  {
    ModularIntType modInt(i, sz);
    EXPECT_EQ(modInt, i);
    EXPECT_EQ(modInt, modInt + sz);
    EXPECT_EQ(modInt, modInt + 2 * sz);
    EXPECT_EQ(modInt, modInt - sz);
  }

  ModularIntType modIntUp(0, sz);
  ModularIntType modIntDn(0, sz);
  const int loopEnd = 3 * modIntUp.modulus();
  for(int i = 0; i < loopEnd; ++i)
  {
    EXPECT_EQ(modIntUp, modIntUp + sz);
    EXPECT_EQ(modIntUp, modIntUp + 2 * sz);
    EXPECT_EQ(modIntUp, modIntUp - sz);

    EXPECT_EQ(modIntDn, modIntDn + sz);
    EXPECT_EQ(modIntDn, modIntDn + 2 * sz);
    EXPECT_EQ(modIntDn, modIntDn - sz);

    ++modIntUp;
    --modIntDn;
  }
}

TEST(slam_modInt, equality)
{
  SLIC_INFO("Checking modular int equality");

  using RTModularInt = slam::ModularInt<RTIntSize>;

  // check equality when both have runtime size policies
  {
    volatile int sz = 7;

    RTModularInt a(4, sz);
    RTModularInt b(4, sz);
    EXPECT_EQ(a, b);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);

    // initially different values
    RTModularInt c(5, sz);
    EXPECT_NE(a, c);

    // ... but are equal after some manipulation
    EXPECT_EQ(a, c - 1);
    c--;
    EXPECT_EQ(a, c);

    // Not equal when modulus is different
    RTModularInt d(4, 2 * sz);
    EXPECT_NE(a, d);
  }

  // check equality when both have compile time size policies
  {
    const int sz = 7;

    slam::ModularInt<IntSize<sz>> a(4);
    slam::ModularInt<IntSize<sz>> b(4);
    EXPECT_EQ(a, b);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);

    // initially different values
    slam::ModularInt<IntSize<sz>> c(5);
    EXPECT_NE(a, c);

    // ... but are equal after some manipulation
    EXPECT_EQ(a, c - 1);
    c--;
    EXPECT_EQ(a, c);

    // Not equal when modulus is different
    slam::ModularInt<IntSize<2 * sz>> d(4);
    EXPECT_NE(a, d);
  }

  // check equality when they different size policies
  {
    const int sz = 7;

    RTModularInt a(4, sz);
    slam::ModularInt<IntSize<sz>> b(4);
    EXPECT_EQ(a, b);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);

    // check that we get the same results if we swap the operands
    EXPECT_TRUE(b == a);
    EXPECT_FALSE(b != a);

    // same modulus with initially different values
    slam::ModularInt<IntSize<sz>> c(5);
    EXPECT_NE(a, c);

    // ... but are equal after some manipulation
    EXPECT_EQ(a, c - 1);
    c--;
    EXPECT_EQ(a, c);

    // Not equal when modulus is different
    slam::ModularInt<IntSize<2 * sz>> d(4);
    EXPECT_NE(a, d);
  }
}

TEST(slam_modInt, copy_and_assign)
{
  SLIC_INFO("Checking modular int copy and assignment");

  using ModularIntType = axom::slam::ModularInt<RTIntSize>;

  volatile int sz = 7;

  // Note: copy assignment only affects its value, but not modulus

  // same modulus
  {
    ModularIntType a(5, sz);
    ModularIntType b(6, sz);
    EXPECT_NE(a, b);

    // copy constructor
    ModularIntType c(a);
    EXPECT_EQ(a, c);

    // copy assignment
    a = b;
    EXPECT_EQ(a, b);
  }

  // different modulus
  {
    ModularIntType a(5, sz);
    ModularIntType b(6, 2 * sz);
    EXPECT_NE(a, b);

    // copy constructor
    ModularIntType c(a);
    EXPECT_EQ(a, c);

    // copy assignment -- different modulus
    a = b;
    EXPECT_NE(a, b);
  }
}

TEST(slam_modInt, runtime_modular_int_mult)
{
  using ModularIntType = slam::ModularInt<RTIntSize>;

  volatile int sz = 10;

  SLIC_INFO("Checking modular int multiplication ");

  ModularIntType modInt5(5, sz);
  EXPECT_EQ(modInt5 * 2, 0);

  ModularIntType modInt2(2, sz);
  EXPECT_EQ(modInt2 * 5, 0);
  EXPECT_EQ(modInt2 * 4, 8);
  EXPECT_EQ(modInt2 * 6, 2);

  ModularIntType modInt3(3, sz);
  EXPECT_EQ(modInt3 * 0, 0);
  EXPECT_EQ(modInt3 * 1, 3);
  EXPECT_EQ(modInt3 * 2, 6);
  EXPECT_EQ(modInt3 * 3, 9);
  EXPECT_EQ(modInt3 * 4, 2);

  ModularIntType modInt13(13, sz);
  EXPECT_EQ(modInt13, 3);
  EXPECT_EQ(modInt13 * 2, 6);
}

TEST(slam_modInt, compiletime_modular_int)
{
  SLIC_INFO("Checking modular int addition and subtraction"
            << " when supplying the max value at compile time");

  const int SZ = 937;

  using ModularIntType = slam::ModularInt<IntSize<SZ>>;

  int sz = SZ;

  ModularIntType modIntZero(sz, sz);
  EXPECT_EQ(modIntZero, 0);

  for(int i = 0; i < sz; ++i)
  {
    ModularIntType modInt(i, sz);
    EXPECT_EQ(modInt, i);
    EXPECT_EQ(modInt, modInt + sz);
    EXPECT_EQ(modInt, modInt + 2 * sz);
    EXPECT_EQ(modInt, modInt - sz);
  }

  ModularIntType modIntUp(0, sz);
  ModularIntType modIntDn(0, sz);
  const int loopEnd = 3 * modIntUp.modulus();
  for(int i = 0; i < loopEnd; ++i)
  {
    EXPECT_EQ(modIntUp, modIntUp + sz);
    EXPECT_EQ(modIntUp, modIntUp + 2 * sz);
    EXPECT_EQ(modIntUp, modIntUp - sz);

    EXPECT_EQ(modIntDn, modIntDn + sz);
    EXPECT_EQ(modIntDn, modIntDn + 2 * sz);
    EXPECT_EQ(modIntDn, modIntDn - sz);

    ++modIntUp;
    --modIntDn;
  }
}

//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger

  result = RUN_ALL_TESTS();
  return result;
}
