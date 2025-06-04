// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/NumericArray.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/CoordinateTransformer.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Point.hpp"

#include <array>

namespace primal = axom::primal;

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_translation()
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  PointType a({1, 2, 3});
  VectorType d({15, 16, 17});
  primal::CoordinateTransformer<double> dt;
  dt.addTranslation(d);
  PointType b(dt.getTransformed(a.array()));
  VectorType diff(b.array() - a.array() - d.array());
  double diffNorm = diff.norm();
  EXPECT_TRUE(axom::utilities::isNearlyEqual(diffNorm, 1e-12));
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, translation)
{
  check_translation<axom::SEQ_EXEC>();
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_rotate_to_axis()
{
  // Check rotations from axis to axis.
  const int DIM = 3;
  using VectorType = primal::Vector<double, DIM>;

  VectorType x({1, 0, 0});
  VectorType y({0, 1, 0});
  VectorType z({0, 0, 1});
  VectorType oct1({1, 1, 1});

  const int n = 15; // Number of pairs in startsAndEnds
  VectorType startsAndEnds[2*n] = {x, x, y, y, z, z,
                                   x, y, y, z, z, x,
                                   x, -x, x, -y, x, -z,
                                   y, -x, y, -y, y, -z,
                                   z, -x, z, -y, z, -z};
  for(int i = 0; i < n; ++i)
  {
    VectorType startDir = startsAndEnds[2*i];
    VectorType endDir = startsAndEnds[2*i + 1];
    primal::CoordinateTransformer<double> rotation;
    rotation.addRotation(startDir, endDir);
    VectorType result(rotation.getTransformed(startDir.array()));
    VectorType diff(result.array() - endDir.array());
    std::cout << startDir << ' ' << endDir << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, rotate_to_axis)
{
  check_rotate_to_axis<axom::SEQ_EXEC>();
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_rotate_about_bisector()
{
  // Check rotations about octant bisectors.
  const int DIM = 3;
  using VectorType = primal::Vector<double, DIM>;

  VectorType x({1, 0, 0});
  VectorType y({0, 1, 0});
  VectorType z({0, 0, 1});
  VectorType oct1({1, 1, 1});
  VectorType oct2({-1, 1, 1});
  VectorType oct3({-1, -1, 1});
  VectorType oct4({1, -1, 1});
  VectorType oct5({1, 1, -1});
  VectorType oct6({-1, 1, -1});
  VectorType oct7({-1, -1, -1});
  VectorType oct8({1, -1, -1});

  double angle = 2*M_PI/3; // 1/3 full rotation goes from axis to axis.

  const int n = 8; // Number of octants and 7-tuples in vectors.
  // The 7-tuples in vectors are the octant bisector and 3 start-end pairs.
  VectorType vectors[7*n] = {oct1, x, y, y, z, z, x,
                             oct2, y, -x, -x, z, z, y,
                             oct3, z, -x, -x, -y, -y, z,
                             oct4, x, z, z, -y, -y, x,
                             oct5, x, -z, -z, y, y, x,
                             oct6, y, -z, -z, -x, -x, y,
                             oct7, -x, -z, -z, -y, -y, -x,
                             oct8, x, -y, -y, -z, -z, x};
  for(int i = 0; i < n; ++i)
  {
    int j = 7*i;
    VectorType rotAxis = vectors[j];
    for(int k=0; k<3; ++k)
    {
      auto startDir = vectors[j + 2*k + 1];
      auto endDir = vectors[j + 2*k + 2];
      primal::CoordinateTransformer<double> rotation;
      rotation.addRotation(rotAxis, angle);
      VectorType result(rotation.getTransformed(startDir.array()));
      VectorType diff(result.array() - endDir.array());
      std::cout << rotAxis << ' ' << startDir << ' ' << endDir << ' ' << result << ' ' << diff << std::endl;
      EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, rotate_about_bisector)
{
  check_rotate_about_bisector<axom::SEQ_EXEC>();
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_translate_rotate()
{
  const int DIM = 3;
  using VectorType = primal::Vector<double, DIM>;

  VectorType x({1, 0, 0});
  VectorType y({0, 1, 0});
  VectorType z({0, 0, 1});

  {
    VectorType pt({1, 1, 1});
    primal::CoordinateTransformer<double> transformer;
    transformer.addTranslation(VectorType{0, 1, 0});
    transformer.addRotation(x, y);
    VectorType correct({-2, 1, 1});
    VectorType result(transformer.getTransformed(pt.array()));
    VectorType diff(result.array() - correct.array());
    std::cout << pt << ' ' << result << ' ' << correct << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, translate_rotate)
{
  check_translate_rotate<axom::SEQ_EXEC>();
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_inverse()
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  {
    // Simple shift
    PointType start({1, 2, 3});
    VectorType shift({10, 20, 30});
    primal::CoordinateTransformer<double> transformer;
    transformer.addTranslation(shift);
    primal::CoordinateTransformer<double> inverse = transformer.getInverse();
    PointType changed = transformer.getTransformed(start);
    PointType undone = inverse.getTransformed(changed);
    VectorType diff(start.array() - undone.array());
    std::cout << start<< ' ' << changed << ' ' << undone << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }

  {
    // Simple rotate
    PointType start({1, 2, 30});
    VectorType axis({1, .5, .2});
    double angle = M_PI/180*10;
    primal::CoordinateTransformer<double> transformer;
    transformer.addRotation(axis, angle);
    primal::CoordinateTransformer<double> inverse = transformer.getInverse();
    PointType changed = transformer.getTransformed(start);
    PointType undone = inverse.getTransformed(changed);
    VectorType diff(start.array() - undone.array());
    std::cout << start<< ' ' << changed << ' ' << undone << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }

  {
    // A bunch of random translate + rotate transforms.
    const int n = 3;
    for(int i = 0; i < n; ++i)
    {
      VectorType startPt;
      startPt[0] = axom::utilities::random_real(-10., 10.);
      startPt[1] = axom::utilities::random_real(-10., 10.);
      startPt[2] = axom::utilities::random_real(-10., 10.);
      VectorType shift({axom::utilities::random_real(-100., 100.),
                        axom::utilities::random_real(-100., 100.),
                        axom::utilities::random_real(-100., 100.)});
      VectorType axis({axom::utilities::random_real(-100., 100.),
                       axom::utilities::random_real(-100., 100.),
                       axom::utilities::random_real(-100., 100.)});
      double angle = axom::utilities::random_real(-2*M_PI, 2*M_PI);

      primal::CoordinateTransformer<double> transformer;
      transformer.addTranslation(shift);
      transformer.addRotation(axis, angle);
      primal::CoordinateTransformer<double> inverse = transformer.getInverse();
      VectorType endPt = transformer.getTransformed(startPt);
      VectorType result = inverse.getTransformed(endPt);
      VectorType diff(result.array() - startPt.array());
      std::cout << startPt << ' ' << endPt << ' ' << result << ' ' << diff << std::endl;
      EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, inverse)
{
  check_inverse<axom::SEQ_EXEC>();
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
