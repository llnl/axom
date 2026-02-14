// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/Determinants.hpp"

#include "axom/fmt.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>

namespace
{

/**
 * \brief Returns the number of IEEE-754 double-precision representable values
 *        (ULPs) separating \p a and \p b.
 *
 * Computes the absolute difference in their ordered bit representations.
 * Assumes IEEE-754 64-bit floats
 * Returns max uint64_t if either argument is NaN; and 0 for equal params (e.g. +0 and -0).
 */
std::uint64_t ulps_apart(double a, double b)
{
  static_assert(sizeof(double) == sizeof(std::uint64_t), "unexpected double size");
  static_assert(std::numeric_limits<double>::is_iec559, "requires IEEE-754");

  if(std::isnan(a) || std::isnan(b))
  {
    return std::numeric_limits<std::uint64_t>::max();
  }

  // Treat +0 and -0 (and any equal values) as distance 0
  if(a == b)
  {
    return 0;
  }

  // Order-preserving bit transform for ULP comparison.
  auto to_ordered = [](double x) -> std::uint64_t {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof(bits));
    constexpr std::uint64_t sign = std::uint64_t {1} << 63;
    return (bits & sign) ? ~bits : (bits + sign);
  };

  const std::uint64_t oa = to_ordered(a);
  const std::uint64_t ob = to_ordered(b);
  return (oa > ob) ? (oa - ob) : (ob - oa);
}

// "raw" implementation of 2x2 determinant for comparisons to axom::utilities::determinant
// uses `volatile` in an attempt to allide compiler optimization that internally use fma calculations
double det2_raw(double a00, double a01, double a10, double a11)
{
  volatile double p0 = a00 * a11;
  volatile double p1 = a01 * a10;
  return p0 - p1;
}

// "raw" implementation of 2x2 determinant for comparisons to axom::utilities::determinant
// uses `volatile` in an attempt to allide compiler optimization that internally use fma calculations
double det3_raw(double a00,
                double a01,
                double a02,
                double a10,
                double a11,
                double a12,
                double a20,
                double a21,
                double a22)
{
  volatile double m01 = a11 * a22 - a12 * a21;
  volatile double m02 = a10 * a22 - a12 * a20;
  volatile double m12 = a10 * a21 - a11 * a20;

  volatile double t0 = a00 * m01;
  volatile double t1 = a01 * m02;
  volatile double t2 = a02 * m12;
  return (t0 - t1) + t2;
}
}  // namespace

TEST(numerics_determinants, determinant_of_In)
{
  const int N = 25;

  for(int i = 2; i < N; ++i)
  {
    axom::numerics::Matrix<double> In = axom::numerics::Matrix<double>::identity(i);
    double det = axom::numerics::determinant(In);
    EXPECT_DOUBLE_EQ(1.0, det);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_determinants, determinant5x5)
{
  const int N = 5;
  const double EPS = 1e-11;

  // clang-format off
  axom::numerics::Matrix<double> A(N, N);
  A(0, 0) = 1;  A(0, 1) = 2;   A(0, 2) = 4;   A(0, 3) = 3;  A(0, 4) = 0;
  A(1, 0) = 2;  A(1, 1) = 1;   A(1, 2) = -1;  A(1, 3) = 1;  A(1, 4) = 3;
  A(2, 0) = 4;  A(2, 1) = -1;  A(2, 2) = -2;  A(2, 3) = 5;  A(2, 4) = 1;  
  A(3, 0) = 7;  A(3, 1) = 3;   A(3, 2) = 6;   A(3, 3) = 2;  A(3, 4) = 1;
  A(4, 0) = 1;  A(4, 1) = 0;   A(4, 2) = -1;  A(4, 3) = 1;  A(4, 4) = 1;
  // clang-format on

  double computed_det = axom::numerics::determinant(A);
  double expected_det = -34.0;
  EXPECT_NEAR(expected_det, computed_det, EPS);
}

//------------------------------------------------------------------------------
TEST(numerics_determinants, check_ulps_apart_function)
{
  constexpr auto zero = std::uint64_t {0};
  constexpr auto one = std::uint64_t {1};
  constexpr auto two = std::uint64_t {2};
  constexpr auto three = std::uint64_t {3};

  // test identity/equality
  {
    EXPECT_EQ(ulps_apart(1., 1.), zero);
    EXPECT_EQ(ulps_apart(-0., +0.), zero);
  }

  // test values around 1
  {
    const double x = 1.;
    const double up = std::nextafter(x, x + 1.);    // next ulp in positive direction
    const double down = std::nextafter(x, x - 1.);  // nex ulp in negative direction
    EXPECT_EQ(ulps_apart(x, up), one);
    EXPECT_EQ(ulps_apart(x, down), one);
    EXPECT_EQ(ulps_apart(down, up), two);

    // check symmetric
    EXPECT_EQ(ulps_apart(x, up), ulps_apart(up, x));

    // check a few increasing values
    double y = x;
    for(int k = 1; k <= 10; ++k)
    {
      y = nextafter(y, y + 1.);
      EXPECT_EQ(ulps_apart(x, y), static_cast<uint64_t>(k));
    }
  }

  // test values around zero
  {
    const double ms = std::numeric_limits<double>::denorm_min();
    const double nms = -ms;
    EXPECT_EQ(ulps_apart(nms, -0.0), one);
    EXPECT_EQ(ulps_apart(+0.0, ms), one);
    EXPECT_EQ(ulps_apart(nms, ms), three);  // -0 to +0 still counts as one ulp

    const double up = std::nextafter(ms, 1.);
    EXPECT_EQ(ulps_apart(ms, up), one);

    const double down = std::nextafter(nms, -1.);
    EXPECT_EQ(ulps_apart(nms, down), one);
  }

  // test inf
  {
    const double mx = std::numeric_limits<double>::max();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_EQ(ulps_apart(mx, nextafter(mx, inf)), one);
    EXPECT_EQ(nextafter(mx, inf), inf);
    EXPECT_EQ(ulps_apart(mx, inf), one);

    const double ninf = -inf;
    const double nmx = -mx;
    EXPECT_EQ(nextafter(nmx, ninf), ninf);
    EXPECT_EQ(ulps_apart(nmx, ninf), one);
  }

  // test NaN
  {
    constexpr uint64_t UMAX = std::numeric_limits<std::uint64_t>::max();
    const double qnan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(ulps_apart(qnan, 1.0), UMAX);
    EXPECT_EQ(ulps_apart(1.0, qnan), UMAX);
    EXPECT_EQ(ulps_apart(qnan, qnan), UMAX);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_determinants, determinant_fma_reduces_ulp_for_cancellation)
{
  // Choose values that are exactly representable in double, but whose products are not.
  // This forces rounding in intermediate products, and the determinant computation is sensitive to cancellation.
  //
  // Let b=c=2^27, a=2^27+1, d=2^27-1.
  // Exact det = a*d - b*c = (2^54 - 1) - 2^54 = -1.
  const double base = std::ldexp(1.0, 27);
  const double a = base + 1.0;
  const double b = base;
  const double c = base;
  const double d = base - 1.0;
  constexpr double expected = -1.0;

  {
    const double det_fma = axom::numerics::determinant(a, b, c, d);
    const double det_raw = det2_raw(a, b, c, d);

    const auto ulp_fma = ulps_apart(det_fma, expected);
    const auto ulp_raw = ulps_apart(det_raw, expected);
    axom::fmt::print(
      "2x2 cancellation: expected={} det_fma={} det_raw={} (ulps from expected: fma={} raw={})\n",
      expected,
      det_fma,
      det_raw,
      ulp_fma,
      ulp_raw);

    EXPECT_DOUBLE_EQ(expected, det_fma);
    EXPECT_LE(ulp_fma, ulp_raw);
  }

  // Embed the same 2x2 cancellation into a 3x3 determinant.
  {
    const double det_fma = axom::numerics::determinant(a, b, 0.0, c, d, 0.0, 0.0, 0.0, 1.0);
    const double det_raw = det3_raw(a, b, 0.0, c, d, 0.0, 0.0, 0.0, 1.0);

    const auto ulp_fma = ulps_apart(det_fma, expected);
    const auto ulp_raw = ulps_apart(det_raw, expected);

    axom::fmt::print(
      "3x3 cancellation: expected={} det_fma={} det_raw={} (ulps from expected: fma={} raw={})\n",
      expected,
      det_fma,
      det_raw,
      ulp_fma,
      ulp_raw);

    EXPECT_DOUBLE_EQ(expected, det_fma);
    EXPECT_LE(ulp_fma, ulp_raw);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_determinants, determinant_fma_reduces_ulp_statistics_2x2)
{
#if !defined(__SIZEOF_INT128__)
  GTEST_SKIP() << "Test requires __int128 for an exact integer reference determinant.";
#else
  // Use many near-cancellation 2x2 determinants with integer inputs.
  // Compare the ULP error against the correctly-rounded double value of the exact integer determinant.

  constexpr std::int64_t base = (1LL << 27);  // exactly representable as double
  constexpr int SAMPLES = 5000;

  std::mt19937_64 rng(42ULL);
  std::uniform_int_distribution<std::int64_t> delta_dist(-2048, 2048);

  long double sum_raw_ulp = 0.0L;
  long double sum_fma_ulp = 0.0L;
  long double sum_delta = 0.0L;

  std::uint64_t max_raw_ulp = 0;
  std::uint64_t max_fma_ulp = 0;
  int better = 0;
  int worse = 0;
  int tied = 0;

  for(int i = 0; i < SAMPLES; ++i)
  {
    // use some explicit test cases, or random numbers otherwise
    std::int64_t da {}, dd {};
    switch(i)
    {
    case 0:
      da = 1.;
      dd = -1.;
      break;
    case 1:
      da = 2048;
      dd = -2048;
      break;
    case 2:
      da = 2048;
      dd = 2028;
      break;
    default:
      da = delta_dist(rng);
      dd = delta_dist(rng);
      break;
    }

    // actual determinant should be an integer: (a_i * d_i) - (b_i * c_i)
    const std::int64_t a_i = base + da;
    const std::int64_t b_i = base;
    const std::int64_t c_i = base;
    const std::int64_t d_i = base + dd;

    const __int128 exact =  //
      static_cast<__int128>(a_i) * static_cast<__int128>(d_i) -
      static_cast<__int128>(b_i) * static_cast<__int128>(c_i);
    const double ref = static_cast<double>(exact);

    const double a = static_cast<double>(a_i);
    const double b = static_cast<double>(b_i);
    const double c = static_cast<double>(c_i);
    const double d = static_cast<double>(d_i);

    const double det_fma = axom::numerics::determinant(a, b, c, d);
    const double det_raw = det2_raw(a, b, c, d);

    const std::uint64_t ulp_fma = ulps_apart(det_fma, ref);
    const std::uint64_t ulp_raw = ulps_apart(det_raw, ref);

    sum_fma_ulp += static_cast<long double>(ulp_fma);
    sum_raw_ulp += static_cast<long double>(ulp_raw);
    sum_delta += static_cast<long double>(ulp_raw) - static_cast<long double>(ulp_fma);

    max_fma_ulp = std::max(max_fma_ulp, ulp_fma);
    max_raw_ulp = std::max(max_raw_ulp, ulp_raw);

    if(ulp_fma < ulp_raw)
    {
      ++better;
    }
    else if(ulp_fma > ulp_raw)
    {
      ++worse;
    }
    else
    {
      ++tied;
    }
  }

  const long double mean_fma = sum_fma_ulp / static_cast<long double>(SAMPLES);
  const long double mean_raw = sum_raw_ulp / static_cast<long double>(SAMPLES);
  const long double mean_delta = sum_delta / static_cast<long double>(SAMPLES);

  axom::fmt::print(
    "2x2 ULP stats for {} samples\n"
    "\t mean_fma={} mean_raw={} mean_delta(raw-fma)={}\n"
    "\t max_fma={} max_raw={}\n"
    "\t better={} worse={} tied={}\n",
    SAMPLES,
    static_cast<double>(mean_fma),
    static_cast<double>(mean_raw),
    static_cast<double>(mean_delta),
    max_fma_ulp,
    max_raw_ulp,
    better,
    worse,
    tied);

  EXPECT_GE(mean_delta, 0.0L);
  EXPECT_GE(better, worse);
#endif
}
