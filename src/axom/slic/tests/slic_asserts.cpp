// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

/*!
 * \file slic_asserts.cpp
 *
 * The tests in this file check that SLIC macros properly output their message
 * and exit (when appropriate) when run from constructors and destructors.
 * They also exercise the SLIC SimpleLogger.
 */

namespace
{
/*!
 * A simple struct with an assert in the constructor
 */
struct AssertCtor
{
  AssertCtor() { SLIC_ASSERT_MSG(false, "Testing assert in .ctor"); }
};

/*!
 *  A simple struct with an assert in a method (foo)
 */
struct AssertMethod
{
  void foo() { SLIC_ASSERT_MSG(false, "Testing assert in class method"); }
};

/*!
 *  A simple struct with an assert in the destructor
 */
struct AssertDtor
{
  ~AssertDtor() { SLIC_ASSERT_MSG(false, "Testing assert in .dtor"); }
};

/*!
 *  A simple testing fixture with a SLIC_WARNING in the constructor.
 *  Note: gtest EXPECT_DEATH_IF_SUPPORTED has a return, so it cannot be
 *  used in a constructor.
 */
class SetFixtureC : public ::testing::Test
{
public:
  SetFixtureC()
  {
    SLIC_WARNING(
      "Testing warning in fixture .ctor -- this warning message should be "
      "logged");
  }
};

/*!
 *  A simple testing fixture with an assert in the SetUp function.
 */
class SetFixtureS : public ::testing::Test
{
public:
  void SetUp()
  {
#ifdef AXOM_DEBUG
    EXPECT_DEATH_IF_SUPPORTED(SLIC_ASSERT_MSG(false, "Testing assert in fixture setup"), "");
#else
    SLIC_WARNING("Testing warning in fixture setup");
#endif
  }
};

/*!
 *  A simple testing fixture with an assert in the TearDown function.
 */
class SetFixtureT : public ::testing::Test
{
public:
  void TearDown()
  {
#ifdef AXOM_DEBUG
    EXPECT_DEATH_IF_SUPPORTED(SLIC_ASSERT_MSG(false, "Testing assert in fixture teardown"), "");
#else
    SLIC_WARNING("Testing warning in fixture teardown");
#endif
  }
};

/*!
 *  A simple testing fixture with a SLIC_WARNING in the destructor.
 *  Note: gtest EXPECT_DEATH_IF_SUPPORTED has a return, so it cannot be used
 *  in a destructor.
 *
 */
class SetFixtureD : public ::testing::Test
{
public:
  ~SetFixtureD()
  {
    SLIC_WARNING(
      "Testing warning in fixture .dtor -- this warning message should be "
      "logged");
  }
};

}  // namespace

// ********  A series of tests exercising SLIC_ASSERT

TEST(slic_usage, in_test)
{
  SLIC_ASSERT_MSG(true, "Testing SLIC assert (true) in test body");
#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(SLIC_ASSERT_MSG(false, "Testing SLIC assert(false) in test body"), "")
    << "SLIC assert (false) from a test";
#else
  EXPECT_DEATH_IF_SUPPORTED(SLIC_ERROR_IF(true, "Testing SLIC error in test body for release mode"),
                            "")
    << "SLIC_ERROR_IF(false) from a test";

#endif
}

TEST(slic_usage, in_ctor)
{
#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(AssertCtor(), "") << " SLIC assert from class .ctor ";
#else
  AssertCtor();
#endif
}

TEST(slic_usage, in_method)
{
  AssertMethod am;
#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(am.foo(), "") << " SLIC assert from class method ";
#else
  am.foo();
#endif
}

TEST(slic_usage, in_dtor)
{
#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(AssertDtor(), "") << " SLIC assert from class .ctor ";
#else
  AssertDtor();
#endif
}

// A test using a test fixture with an assert in the setup phase
TEST_F(SetFixtureS, in_fixture_setup) { }

// A test using a test fixture with an assert in the teardown phase
TEST_F(SetFixtureT, in_fixture_teardown) { }

// Note (KW): the following two tests are warnings since ASSERT_DEATH does not
// work in a constructor or destructor. Specifically, the ASSERT_DEATH macro
// has a return statement which is not allowed in constructors or destructors

// A test using a test fixture with an assert in the constructor
TEST_F(SetFixtureC, in_fixture_ctor) { }

// A test using a test fixture with an assert in the destructor
TEST_F(SetFixtureD, in_fixture_dtor) { }

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  result = RUN_ALL_TESTS();

  return result;
}
