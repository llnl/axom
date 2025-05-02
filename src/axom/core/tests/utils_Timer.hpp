// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Timer.hpp"

#ifdef WIN32
  #include "windows.h"  // for Sleep(), in milliseconds
#else
  #include <unistd.h>  // for usleep(), in microseconds
#endif

// cross-platform utility for calling sleep function, in milliseconds
void sleep_ms(std::uint32_t ms)
{
#ifdef WIN32
  Sleep(ms);  // already in milliseconds
#else
  usleep(ms * 1000);  // convert to microseconds
#endif
}

TEST(utils_Timer, timer_check)
{
  axom::utilities::Timer t;

  std::cout << "Checking that new timer indicates 0 time elapsed" << std::endl;
  EXPECT_EQ(0., t.elapsed());

  t.start();

  const std::uint32_t sleep_duration_ms = 100;
  sleep_ms(sleep_duration_ms);

  t.stop();

  std::cout << "Simple test for elapsed time in different units." << std::endl;
  EXPECT_DOUBLE_EQ(t.elapsedTimeInMicroSec(), 1000. * t.elapsedTimeInMilliSec());
  EXPECT_DOUBLE_EQ(t.elapsedTimeInMilliSec(), 1000. * t.elapsedTimeInSec());
  EXPECT_EQ(t.elapsed(), t.elapsedTimeInSec());

  std::cout << "Testing that reset() indicates 0 elapsed time." << std::endl;
  t.reset();
  ASSERT_DOUBLE_EQ(0., t.elapsed());
}

TEST(utils_Timer, timer_check_duration)
{
  const std::uint32_t sleep_duration_ms = 100;
  const double sleep_duration_s = sleep_duration_ms / 1000.;

  axom::utilities::Timer t;
  t.start();

  sleep_ms(sleep_duration_ms);

  t.stop();
  const double e = t.elapsed();
  std::cout << "Elapsed: " << e << std::endl;

  EXPECT_GE(e, sleep_duration_s);
  EXPECT_LT(e, sleep_duration_s + 1);
}

TEST(utils_Timer, timer_check_sum)
{
  const int N = 3;
  const std::uint32_t sleep_duration_ms = 100;
  const double sleep_duration_s = sleep_duration_ms / 1000.;

  axom::utilities::Timer t1(false);
  axom::utilities::Timer t2(false);
  for(int n = 0; n < N; ++n)
  {
    t2.start();
    t1.start();
    sleep_ms(sleep_duration_ms);
    t1.stop();
    sleep_ms(sleep_duration_ms);
    t2.stop();
  }

  std::cout << "t1 measured: " << t1.elapsed() << "s in " << t1.cycleCount() << " cycles\n";
  std::cout << "t2 measured: " << t2.elapsed() << "s in " << t2.cycleCount() << " cycles\n";

  EXPECT_EQ(t1.cycleCount(), N);
  EXPECT_EQ(t2.cycleCount(), N);

  EXPECT_GE(t1.elapsed(), N * sleep_duration_s);
  EXPECT_LT(t1.elapsed(), N * sleep_duration_s + 1);

  EXPECT_GE(t2.elapsed(), 2 * N * sleep_duration_s);
  EXPECT_LT(t2.elapsed(), 2 * N * sleep_duration_s + 1);
}
