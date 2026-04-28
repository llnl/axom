// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/ArrayView.hpp"
#include "axom/core/numerics/matvecops.hpp"
#include "axom/core/numerics/transforms.hpp"

namespace
{
void expect_matrix_near(const axom::numerics::Matrix<double>& actual,
                        const axom::numerics::Matrix<double>& expected,
                        double tolerance = 1e-12)
{
  ASSERT_EQ(actual.getNumRows(), expected.getNumRows());
  ASSERT_EQ(actual.getNumColumns(), expected.getNumColumns());

  for(axom::IndexType i = 0; i < actual.getNumRows(); ++i)
  {
    for(axom::IndexType j = 0; j < actual.getNumColumns(); ++j)
    {
      EXPECT_NEAR(actual(i, j), expected(i, j), tolerance);
    }
  }
}

void expect_vector_near(const double* actual,
                        const double* expected,
                        int size,
                        double tolerance = 1e-12)
{
  for(int i = 0; i < size; ++i)
  {
    EXPECT_NEAR(actual[i], expected[i], tolerance);
  }
}
}  // namespace

TEST(numerics_transforms, scale_2d_about_center)
{
  double centerData[2] = {2.0, 3.0};
  axom::ArrayView<double> center(centerData, 2);

  auto actual = axom::numerics::transforms::scale(4.0, 5.0, center);

  axom::numerics::Matrix<double> expected = axom::numerics::Matrix<double>::identity(3);
  expected(0, 0) = 4.0;
  expected(1, 1) = 5.0;
  expected(0, 2) = -6.0;
  expected(1, 2) = -12.0;
  expect_matrix_near(actual, expected);

  double point[3] = {3.0, 4.0, 1.0};
  double result[3] = {0.0, 0.0, 0.0};
  double expectedPoint[3] = {6.0, 8.0, 1.0};
  axom::numerics::matrix_vector_multiply(actual, point, result);
  expect_vector_near(result, expectedPoint, 3);
}

TEST(numerics_transforms, scale_3d_about_center)
{
  double centerData[3] = {0.5, 0.5, 0.5};
  axom::ArrayView<double> center(centerData, 3);

  auto actual = axom::numerics::transforms::scale(2.0, 3.0, 4.0, center);

  axom::numerics::Matrix<double> expected = axom::numerics::Matrix<double>::identity(4);
  expected(0, 0) = 2.0;
  expected(1, 1) = 3.0;
  expected(2, 2) = 4.0;
  expected(0, 3) = -0.5;
  expected(1, 3) = -1.0;
  expected(2, 3) = -1.5;
  expect_matrix_near(actual, expected);

  double point[4] = {1.5, 1.5, 1.5, 1.0};
  double result[4] = {0.0, 0.0, 0.0, 0.0};
  double expectedPoint[4] = {2.5, 3.5, 4.5, 1.0};
  axom::numerics::matrix_vector_multiply(actual, point, result);
  expect_vector_near(result, expectedPoint, 4);
}

TEST(numerics_transforms, scale_about_zero_center_matches_origin_scale)
{
  double center2dData[2] = {0.0, 0.0};
  axom::ArrayView<double> center2d(center2dData, 2);
  expect_matrix_near(axom::numerics::transforms::scale(2.0, 3.0, center2d),
                     axom::numerics::transforms::scale(2.0, 3.0));

  double center3dData[3] = {0.0, 0.0, 0.0};
  axom::ArrayView<double> center3d(center3dData, 3);
  expect_matrix_near(axom::numerics::transforms::scale(2.0, 3.0, 4.0, center3d),
                     axom::numerics::transforms::scale(2.0, 3.0, 4.0, 4));
}
