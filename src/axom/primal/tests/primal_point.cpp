// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Point.hpp"
#include "axom/core/NumericArray.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/slic.hpp"

using namespace axom;

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_point_policy()
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;

  double* coords = axom::allocate<double>(DIM, axom::execution_space<ExecSpace>::allocatorID());

  double coords_host[DIM];

  axom::for_all<ExecSpace>(
    1,
    AXOM_LAMBDA(int /*i*/) {
      PointType ones = PointType::ones();
      ones.to_array(coords);
    });

  axom::copy(&coords_host, coords, DIM * sizeof(double));

  EXPECT_EQ(PointType(coords_host), PointType::ones());
  axom::deallocate(coords);
}

//------------------------------------------------------------------------------
TEST(primal_point, point_default_constructor)
{
  static const int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;

  QPoint pt;

  EXPECT_EQ(pt[0], CoordType());
  EXPECT_EQ(pt[1], CoordType());
  EXPECT_EQ(pt[2], CoordType());
  EXPECT_EQ(pt.dimension(), DIM);
}

//------------------------------------------------------------------------------
TEST(primal_point, point_singleVal_constructor)
{
  static const int DIM = 5;
  using CoordType = int;
  using QPoint = primal::Point<CoordType, DIM>;
  const int singleVal = 10;

  //
  SLIC_INFO("\nprimal: testing constructor that sets all values to a singleVal."
            << "Default second parameter of constructor is DIM.");
  QPoint pt1(singleVal);
  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt1[dim], singleVal);
  }

  //
  SLIC_INFO("\nprimal: testing constructor that sets all values to a singleVal."
            << "Using explicit second parameter set to DIM.");
  QPoint pt2(singleVal, DIM);
  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt2[dim], singleVal);
  }

  //
  SLIC_INFO("\nprimal: testing constructor that sets all values to a singleVal."
            << "Using explicit second parameter set higher than DIM.");
  QPoint pt3(singleVal, DIM * 2);
  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt3[dim], singleVal);
  }

  //
  SLIC_INFO("\nprimal: testing constructor that sets *some* values a singleVal."
            << "Using explicit second parameter set less than DIM."
            << "Other values should be set to zero.");

  int numVals = DIM / 2;
  QPoint pt4(singleVal, numVals);
  for(int dim = 0; dim < numVals; ++dim)
  {
    EXPECT_EQ(pt4[dim], singleVal);
  }

  for(int dim = numVals + 1; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt4[dim], CoordType());
  }
}

//------------------------------------------------------------------------------
TEST(primal_point, point_array_constructor)
{
  static const int DIM = 5;
  using CoordType = int;
  using QPoint = primal::Point<CoordType, DIM>;

  // Set elt i of input array to i
  CoordType arr[DIM];
  for(int dim = 0; dim < DIM; ++dim)
  {
    arr[dim] = dim;
  }

  //
  SLIC_INFO("\nprimal: testing constructor that copies entire array. "
            << "Default second parameter of constructor is DIM.");
  QPoint pt1(arr);

  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt1[dim], arr[dim]);
  }

  //
  SLIC_INFO("\nprimal: testing constructor that copies entire array. "
            << "Using explicit second parameter set to DIM.");
  QPoint pt2(arr, DIM);

  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt2[dim], arr[dim]);
  }

  //
  SLIC_INFO("\nprimal: testing constructor that copies entire array. "
            << "Using explicit second parameter set higher than DIM.");
  QPoint pt3(arr, DIM * 2);

  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt3[dim], arr[dim]);
  }

  //
  SLIC_INFO("\nprimal: testing constructor that sets *some* values to a singleVal. "
            << "Using explicit second parameter set less than DIM. "
            << "Other values should be set to zero. ");

  int numVals = DIM / 2;
  QPoint pt4(arr, numVals);

  for(int dim = 0; dim < numVals; ++dim)
  {
    EXPECT_EQ(pt4[dim], arr[dim]);
  }

  for(int dim = numVals + 1; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt4[dim], CoordType());
  }
}

//------------------------------------------------------------------------------
TEST(primal_point, point_numericArray_constructor)
{
  static const int DIM = 5;
  using CoordType = int;
  using QArray = axom::NumericArray<CoordType, DIM>;
  using QPoint = primal::Point<CoordType, DIM>;

  // Set elt i of input array to i
  CoordType arr[DIM];
  for(int dim = 0; dim < DIM; ++dim)
  {
    arr[dim] = dim;
  }

  //
  SLIC_INFO("\nprimal: testing constructor that copies entire array. "
            << "Default second parameter of constructor is DIM.");
  QArray arr1(arr);
  QPoint pt(arr);

  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(pt[dim], arr[dim]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_point, point_initializerList_constructor)
{
  primal::Point<int, 3> fromInitializerList = {10, 20, 30};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerList[i]);
  }

  primal::Point<int, 4> listTooShort = {10, 20};
  for(int i = 0; i < 2; ++i)
  {
    EXPECT_EQ(10 * (i + 1), listTooShort[i]);
  }
  for(int i = 2; i < 4; ++i)
  {
    EXPECT_EQ(0, listTooShort[i]);
  }

  primal::Point<int, 3> listTooLong = {10, 20, 30, 40};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), listTooLong[i]);
  }

  primal::Point<int, 3> noEqualsSign {10, 20, 30};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), noEqualsSign[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_point, point_copy_and_assignment)
{
  static const int DIM = 5;
  using CoordType = int;
  using QPoint = primal::Point<CoordType, DIM>;

  // Set elt i of input array to i
  CoordType arr[DIM];
  for(int dim = 0; dim < DIM; ++dim)
  {
    arr[dim] = dim;
  }
  QPoint ptOrig(arr);

  //
  SLIC_INFO("\nprimal: testing copy constructor");
  QPoint ptCopy1(ptOrig);
  QPoint ptCopy2 = ptOrig;

  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(ptOrig[dim], ptCopy1[dim]);
    EXPECT_EQ(ptOrig[dim], ptCopy2[dim]);
  }

  //
  SLIC_INFO("\nprimal: testing assignment operator");
  QPoint ptAssign(5);
  ptAssign = ptOrig;
  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(ptOrig[dim], ptAssign[dim]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_point, point_equality)
{
  static const int DIM = 5;
  using CoordType = int;
  using QPoint = primal::Point<CoordType, DIM>;

  // Set elt i of input array to i
  CoordType arr[DIM];
  for(int dim = 0; dim < DIM; ++dim)
  {
    arr[dim] = dim;
  }
  QPoint ptOrig(arr);

  //
  SLIC_INFO("\nprimal: testing equality of same point");
  QPoint ptCopy1(ptOrig);

  EXPECT_TRUE(ptOrig == ptCopy1);
  EXPECT_FALSE(ptOrig != ptCopy1);

  //
  SLIC_INFO("\nprimal: testing inequality of different points");
  QPoint ptDifferent(7);
  EXPECT_FALSE(ptOrig == ptDifferent);
  EXPECT_TRUE(ptOrig != ptDifferent);

  //
  SLIC_INFO("\nprimal: Testing that zero() and ones()");
  EXPECT_EQ(QPoint::zero(), QPoint());
  EXPECT_EQ(QPoint::zero(), QPoint(0));
  EXPECT_EQ(QPoint::ones(), QPoint(1));
}

//------------------------------------------------------------------------------
TEST(primal_point, point_to_array)
{
  static const int DIM = 5;
  using CoordType = int;
  using QPoint = primal::Point<CoordType, DIM>;

  // Set elt i of input array to i
  CoordType arr[DIM];
  for(int dim = 0; dim < DIM; ++dim)
  {
    arr[dim] = dim;
  }
  QPoint pt(arr);

  CoordType outputArr[DIM];
  pt.to_array(outputArr);

  for(int dim = 0; dim < DIM; ++dim)
  {
    EXPECT_EQ(arr[dim], outputArr[dim]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_point, point_make_point)
{
  static const int DIM = 3;
  using CoordType = int;
  using QPoint = primal::Point<CoordType, DIM>;

  const int x = 10;
  const int y = 20;
  const int z = 30;

  SLIC_INFO("\nprimal: Testing make_point with two coordinates.");
  QPoint pt2 = QPoint::make_point(x, y);
  EXPECT_EQ(pt2[0], x);
  EXPECT_EQ(pt2[1], y);
  EXPECT_EQ(pt2[2], 0);

  //
  SLIC_INFO("\nprimal: Testing make_point with three coordinates.");
  QPoint pt3 = QPoint::make_point(x, y, z);
  EXPECT_EQ(pt3[0], x);
  EXPECT_EQ(pt3[1], y);
  EXPECT_EQ(pt3[2], z);
}

//------------------------------------------------------------------------------
TEST(primal_point, point_midpoint)
{
  static const int DIM = 3;
  using CoordType = int;
  using QPoint = primal::Point<CoordType, DIM>;

  QPoint p10(10);
  QPoint p30(30);
  QPoint p50(50);

  EXPECT_TRUE(p30 == QPoint::midpoint(p10, p50));
}

//------------------------------------------------------------------------------
TEST(primal_point, point_linear_interpolation)
{
  constexpr int DIM = 3;
  using QPoint = primal::Point<double, DIM>;

  QPoint p0;
  QPoint p1(100);

  EXPECT_TRUE(QPoint::lerp(p0, p1, 0) == p0);
  EXPECT_TRUE(QPoint::lerp(p0, p1, 1) == p1);
  EXPECT_TRUE(QPoint::lerp(p0, p1, 0.5) == QPoint::midpoint(p0, p1));
  EXPECT_TRUE(QPoint::lerp(p0, p1, .25) == QPoint(25));
  EXPECT_TRUE(QPoint::lerp(p0, p1, .75) == QPoint(75));
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(primal_point, point_check_policies)
{
  using seq_exec = axom::SEQ_EXEC;
  check_point_policy<seq_exec>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using omp_exec = axom::OMP_EXEC;
  check_point_policy<omp_exec>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_point_policy<cuda_exec>();
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_point_policy<hip_exec>();
#endif
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
