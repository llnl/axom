// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Timer.hpp"

#ifdef WIN32
  #include "windows.h"
void sleep(int numSeconds)
{
  int numMilliSecs = numSeconds * 1000;
  Sleep(numMilliSecs);
}
#else
  #include <unistd.h>  // for sleep()
#endif

TEST(utils_Timer, timer_check)
{
  axom::utilities::Timer t;

  std::cout << "Checking that new timer indicates 0 time elapsed" << std::endl;
  EXPECT_EQ(0., t.elapsed());

  t.start();

  sleep(1);

  t.stop();

  std::cout << "Simple test for elapsed time in different units." << std::endl;
  EXPECT_GT(t.elapsedTimeInMicroSec(), t.elapsedTimeInMilliSec());
  EXPECT_GT(t.elapsedTimeInMilliSec(), t.elapsedTimeInSec());
  EXPECT_EQ(t.elapsed(), t.elapsedTimeInSec());

  std::cout << "Testing that reset() indicates 0 elapsed time." << std::endl;
  t.reset();
  ASSERT_DOUBLE_EQ(0., t.elapsed());
}

TEST(utils_Timer, timer_check_duration)
{
  axom::utilities::Timer t;
  t.start();

  sleep(1);

  t.stop();
  double e = t.elapsed();
  std::cout << "Elapsed: " << e << std::endl;

  EXPECT_GE(e, 1.0);
  EXPECT_LT(e, 2.0);
}

TEST(utils_Timer, timer_check_sum)
{
  const int N = 3;
  axom::utilities::Timer t1(false);
  axom::utilities::Timer t2(false);
  for(int n = 0; n < N; ++n)
  {
    t2.start();
    t1.start();
    sleep(1);
    t1.stop();
    sleep(1);
    t2.stop();
  }

  std::cout << "t1 measured: " << t1.elapsed() << "s in " << t1.cycleCount() << " cycles"
            << std::endl;
  std::cout << "t2 measured: " << t2.elapsed() << "s in " << t2.cycleCount() << " cycles"
            << std::endl;

  EXPECT_EQ(t1.cycleCount(), N);
  EXPECT_EQ(t2.cycleCount(), N);

  // Inaccuracy per start-stop cycle, estimated from sample runs.
#if defined(_WIN32)
  const double tol = 0.03;
#elif defined(__OSX__) || defined(__APPLE__)
  const double tol = 0.3;
#else
  const double tol = 0.0004;
#endif

  EXPECT_GE(t1.elapsed() / N, 1 - tol);
  EXPECT_LE(t1.elapsed() / N, 1 + tol);
  EXPECT_GE(t2.elapsed() / N, 2 - tol);
  EXPECT_LE(t2.elapsed() / N, 2 + tol);
}
