// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic.hpp"

#ifdef AXOM_USE_CALIPER
  #include <caliper/cali.h>
#endif

#include "axom/core/Types.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"

#include "axom/primal/geometry/Point.hpp"

#include "axom/primal/operators/clip.hpp"

#include <limits>
#include <cmath>

#include "gtest/gtest.h"

/*
 * Clip various combinations of primitives and report time taken.
 */

namespace Primal3D
{
using PointType = axom::primal::Point<double, 3>;
using TetrahedronType = axom::primal::Tetrahedron<double, 3>;
}  // namespace Primal3D

template <typename ExecSpace>
double array_sum(axom::ArrayView<double> v)
{
  axom::ReduceSum<ExecSpace, double> sum(0);
  axom::for_all<ExecSpace>(v.size(), AXOM_LAMBDA(axom::IndexType i) { sum += v[i]; });
  double rval = sum.get();
  return rval;
}

constexpr double EPS = 1e-10;
constexpr bool tryFixOrientation = false;
constexpr axom::IndexType repCount = 10000;  // For reliable timings, try 1e6.

template <typename ExecSpace>
void time_repeat_clips(const Primal3D::TetrahedronType &a,
                       const Primal3D::TetrahedronType &b,
                       axom::IndexType count,
                       const std::string &caseName)
{
  using namespace Primal3D;
  const std::string timerName = caseName + axom::execution_space<ExecSpace>::name();
  auto poly = axom::primal::clip(a, b, EPS, tryFixOrientation);
  const double singleVol = poly.volume();

  SLIC_INFO(axom::fmt::format("{} with {} repetitions, volume {}", caseName, count, singleVol));

  const int allocId = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<double> vols(axom::ArrayOptions::Uninitialized(), count, count, allocId);

  axom::Array<TetrahedronType> as(axom::ArrayOptions::Uninitialized(), 0, count, allocId);
  axom::Array<TetrahedronType> bs(axom::ArrayOptions::Uninitialized(), 0, count, allocId);

  as.insert(0, count, a);
  bs.insert(0, count, b);

  auto asView = as.view();
  auto bsView = bs.view();
  auto volsView = vols.view();

  AXOM_ANNOTATE_BEGIN(timerName);
  axom::for_all<ExecSpace>(
    count,
    AXOM_LAMBDA(axom::IndexType i) {
      auto poly = axom::primal::clip(asView[i], bsView[i], EPS, tryFixOrientation);
      volsView[i] = poly.volume();
    });
  AXOM_ANNOTATE_END(timerName);

  // Verify correctness.
  double avgVol = array_sum<ExecSpace>(vols.view()) / count;
  EXPECT_NEAR(avgVol, singleVol, EPS);
}

void time_repeat_clips_all(const Primal3D::TetrahedronType &a,
                           const Primal3D::TetrahedronType &b,
                           axom::IndexType count,
                           const std::string &caseName)
{
  time_repeat_clips<axom::SEQ_EXEC>(a, b, count, caseName);

#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  {
    time_repeat_clips<axom::OMP_EXEC>(a, b, count, caseName);
  }
#endif

#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  {
    time_repeat_clips<axom::CUDA_EXEC<256>>(a, b, count, caseName);
  }
#endif

#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  {
    time_repeat_clips<axom::HIP_EXEC<256>>(a, b, count, caseName);
  }
#endif
}

//
// no overlap and should be as fast as possible.
//
TEST(primal_clip, fast_miss)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  /*
    backward and forward don't intersect and no plane from one intersects the other.
    so their clipping should be as fast as possible.
  */
  Primal3D::TetrahedronType backward({1, -1, -1}, {-1, -1, -1}, {0, 1, -1}, {0, 0, -2});
  Primal3D::TetrahedronType forward({1, -1, 1}, {0, 1, 1}, {-1, -1, 1}, {0, 0, 2});
  time_repeat_clips_all(backward, forward, repCount, name);
}

//
// the first is outside the second.
//
TEST(primal_clip, outside)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType in({1, -1, 0}, {0, 1, 0}, {-1, -1, 0}, {0, 0, 1});
  Primal3D::TetrahedronType out({2, -2, 0}, {0, 2, 0}, {-2, -2, 0}, {0, 0, 2});
  time_repeat_clips_all(in, out, repCount, name);
}

//
// the first is inside the second.
//
TEST(primal_clip, inside)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType in({1, -1, 0}, {0, 1, 0}, {-1, -1, 0}, {0, 0, 1});
  Primal3D::TetrahedronType out({2, -2, 0}, {0, 2, 0}, {-2, -2, 0}, {0, 0, 2});
  time_repeat_clips_all(out, in, repCount, name);
}

//
// one's tip pokes the other's side.
//
TEST(primal_clip, poked_side)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType sharp({1, -1, 0}, {0, 1, 0}, {-1, -1, 0}, {0, 0, 1});
  Primal3D::TetrahedronType poked({1, -1, 0.5}, {0, 1, 0.5}, {-1, -1, 0.5}, {0, 0, 1.5});
  time_repeat_clips_all(sharp, poked, repCount, name);
}

//
// one's tip impales the other's side.
//
TEST(primal_clip, impaled_side)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType sharp({1, -1, 0}, {0, 1, 0}, {-1, -1, 0}, {0, 0, 2});
  Primal3D::TetrahedronType impaled({1, -1, 0.5}, {0, 1, 0.5}, {-1, -1, 0.5}, {0, 0, 1.5});
  time_repeat_clips_all(sharp, impaled, repCount, name);
}

//
// each has an edge through the other.
//
TEST(primal_clip, edge_through)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType high({1, 0, -0.5}, {-1, 0, -0.5}, {0, -1, 1}, {0, 1, 1});
  Primal3D::TetrahedronType low({0, 1, 0.5}, {0, -1, 0.5}, {-1, 0, -1}, {1, 0, -1});
  time_repeat_clips_all(high, low, repCount, name);
}

//
// the tips of each crash head-on.
//
TEST(primal_clip, head_on)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType backward({1, -1, 0.5}, {-1, -1, 0.5}, {0, 1, 0.5}, {0, 0, -0.5});
  Primal3D::TetrahedronType forward({1, -1, -0.5}, {0, 1, -0.5}, {-1, -1, -0.5}, {0, 0, 0.5});
  time_repeat_clips_all(backward, forward, repCount, name);
}

//
// the tips of each crash head-on and one goes all
// the way through.
//
TEST(primal_clip, head_on2)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  /*
    backward and forward don't intersect and no plane from one intersects the other.
    so their clipping should be as fast as possible.
  */
  Primal3D::TetrahedronType backward({1, -1, 0.5}, {-1, -1, 0.5}, {0, 1, 0.5}, {0, 0, -0.5});
  Primal3D::TetrahedronType forward({1, -1, -0.5}, {0, 1, -0.5}, {-1, -1, -0.5}, {0, 0, 1.5});
  time_repeat_clips_all(backward, forward, repCount, name);
}

//
// the tips of each crash head-on and both go all
// the way through.
//
TEST(primal_clip, head_on3)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  /*
    backward and forward don't intersect and no plane from one intersects the other.
    so their clipping should be as fast as possible.
  */
  Primal3D::TetrahedronType backward({1, -1, 0.5}, {-1, -1, 0.5}, {0, 1, 0.5}, {0, 0, -1.5});
  Primal3D::TetrahedronType forward({1, -1, -0.5}, {0, 1, -0.5}, {-1, -1, -0.5}, {0, 0, 1.5});
  time_repeat_clips_all(backward, forward, repCount, name);
}

//
// the tets have the same center
// but the points of one stick out the sides of the other
// so the union looks like an 8-pointed ball.
//
TEST(primal_clip, eight_point)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType a({1, 0, 0}, {-1, 0, 0}, {0, -1, 1}, {0, 1, 1});
  Primal3D::TetrahedronType b({1, 0, 1}, {0, -1, 0}, {-1, 0, 1}, {0, 1, 0});
  time_repeat_clips_all(a, b, repCount, name);
}

//
// like eight_point, but shift the first 20% along +z
//
TEST(primal_clip, eight_point2)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType a({1, 0, 0.2}, {-1, 0, 0.2}, {0, -1, 1.2}, {0, 1, 1.2});
  Primal3D::TetrahedronType b({1, 0, 1}, {0, -1, 0}, {-1, 0, 1}, {0, 1, 0});
  time_repeat_clips_all(a, b, repCount, name);
}

//
// like eight_point, but shift the first half-way along +z
//
TEST(primal_clip, eight_point3)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType a({1, 0, 0.5}, {-1, 0, 0.5}, {0, -1, 1.5}, {0, 1, 1.5});
  Primal3D::TetrahedronType b({1, 0, 1}, {0, -1, 0}, {-1, 0, 1}, {0, 1, 0});
  time_repeat_clips_all(a, b, repCount, name);
}

//
// like eight_point, but shift the first 80% along +z
//
TEST(primal_clip, eight_point4)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType a({1, 0, 0.8}, {-1, 0, 0.8}, {0, -1, 1.8}, {0, 1, 1.8});
  Primal3D::TetrahedronType b({1, 0, 1}, {0, -1, 0}, {-1, 0, 1}, {0, 1, 0});
  time_repeat_clips_all(a, b, repCount, name);
}

//
// like eight_point, but shift the first 150% along +z, resulting in no overlap
//
TEST(primal_clip, eight_point5)
{
  std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  Primal3D::TetrahedronType a({1, 0, 1.5}, {-1, 0, 1.5}, {0, -1, 2.5}, {0, 1, 2.5});
  Primal3D::TetrahedronType b({1, 0, 1}, {0, -1, 0}, {-1, 0, 1}, {0, 1, 0});
  time_repeat_clips_all(a, b, repCount, name);
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

#ifdef AXOM_USE_CALIPER
  std::string caliAnnotationMode = "report";
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(caliAnnotationMode);
#else
  SLIC_INFO("No timer report.  Enable Caliper to get it.");
#endif

  int result = RUN_ALL_TESTS();

  return result;
}
