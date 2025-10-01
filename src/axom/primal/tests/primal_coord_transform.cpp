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
#include <math.h>

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_coord_transform, matrix_consistency)
{
  axom::numerics::Matrix<double> m(4, 4);
  for(int r = 0; r < 3; ++r)
  {
    for(int c = 0; c < 4; ++c)
    {
      m(r, c) = 10 * r + c;
    }
  }
  m(3, 0) = m(3, 1) = m(3, 2) = 0.0;
  m(3, 3) = 1.0;

  primal::experimental::CoordinateTransformer<double> ct(m);
  EXPECT_EQ(ct.getMatrix(), m);

  primal::experimental::CoordinateTransformer<double> ct1(ct);
  EXPECT_EQ(ct1.getMatrix(), ct.getMatrix());
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, translation)
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  PointType a({1, 2, 3});
  VectorType d({15, 16, 17});
  primal::experimental::CoordinateTransformer<double> dt;
  dt.applyTranslation(d);
  PointType b(dt.getTransformed(a.array()));
  VectorType diff(b.array() - a.array() - d.array());
  double diffNorm = diff.norm();
  EXPECT_TRUE(axom::utilities::isNearlyEqual(diffNorm, 1e-12));
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, rotate_to_axis)
{
  // Check rotations from axis to axis.
  const int DIM = 3;
  using NumArrayType = axom::NumericArray<double, DIM>;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  // Points along axes
  NumArrayType x({1, 0, 0});
  NumArrayType y({0, 1, 0});
  NumArrayType z({0, 0, 1});

  const int n = 15;  // Number of pairs in startsAndEnds
  // clang-format off
  NumArrayType startsAndEnds[n][2] = {{x,  x}, {y,  y}, {z,  z},
                                      {x,  y}, {y,  z}, {z,  x},
                                      {x, -x}, {x, -y}, {x, -z},
                                      {y, -x}, {y, -y}, {y, -z},
                                      {z, -x}, {z, -y}, {z,  -z}};
  // clang-format on

  for(int i = 0; i < n; ++i)
  {
    VectorType startDir(startsAndEnds[i][0]);
    VectorType endDir(startsAndEnds[i][1]);
    primal::experimental::CoordinateTransformer<double> rotation;
    rotation.applyRotation(startDir, endDir);
    if(startDir == -endDir)
    {
      // Ill-defined rotation: check for invalid transformer.
      EXPECT_FALSE(rotation.isValid());
      continue;
    }
    PointType result(rotation.getTransformed(startDir.array()));
    VectorType diff(result.array() - endDir.array());
    std::cout << startDir << ' ' << endDir << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, rotate_about_bisector)
{
  // Check rotations about octant bisectors.
  const int DIM = 3;
  using NumArrayType = axom::NumericArray<double, DIM>;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  // Points along axes
  NumArrayType x({1, 0, 0});
  NumArrayType y({0, 1, 0});
  NumArrayType z({0, 0, 1});
  // Octant bisectors
  constexpr int nOct = 8;  // Number of octants
  VectorType oct1({1, 1, 1});
  VectorType oct2({-1, 1, 1});
  VectorType oct3({-1, -1, 1});
  VectorType oct4({1, -1, 1});
  VectorType oct5({1, 1, -1});
  VectorType oct6({-1, 1, -1});
  VectorType oct7({-1, -1, -1});
  VectorType oct8({1, -1, -1});

  double angle = 2 * M_PI / 3;  // 1/3 full rotation about bisectors goes from axis to axis.

  // Test in 8 octants.  Do 3 rotations per octant.  Each rotation has a {start, end} pair.
  VectorType rotAxes[nOct] = {oct1.unitVector(),
                              oct2.unitVector(),
                              oct3.unitVector(),
                              oct4.unitVector(),
                              oct5.unitVector(),
                              oct6.unitVector(),
                              oct7.unitVector(),
                              oct8.unitVector()};
  NumArrayType startEnds[nOct][3][2] = {{{x, y}, {y, z}, {z, x}},
                                        {{y, -x}, {-x, z}, {z, y}},
                                        {{z, -x}, {-x, -y}, {-y, z}},
                                        {{x, z}, {z, -y}, {-y, x}},
                                        {{x, -z}, {-z, y}, {y, x}},
                                        {{y, -z}, {-z, -x}, {-x, y}},
                                        {{-x, -z}, {-z, -y}, {-y, -x}},
                                        {{x, -y}, {-y, -z}, {-z, x}}};
  for(int i = 0; i < nOct; ++i)
  {
    const VectorType& rotAxis = rotAxes[i];
    for(int k = 0; k < 3; ++k)
    {
      PointType startPt(startEnds[i][k][0]);
      PointType endPt(startEnds[i][k][1]);
      primal::experimental::CoordinateTransformer<double> rotation;
      rotation.applyRotation(rotAxis, angle);
      PointType result(rotation.getTransformed(startPt.array()));
      VectorType diff(result.array() - endPt.array());
      std::cout << rotAxis << ' ' << startPt << ' ' << endPt << ' ' << result << ' ' << diff
                << std::endl;
      EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, translate_rotate)
{
  const int DIM = 3;
  using VectorType = primal::Vector<double, DIM>;

  VectorType x({1, 0, 0});
  VectorType y({0, 1, 0});
  VectorType z({0, 0, 1});

  {
    VectorType pt({1, 1, 1});
    primal::experimental::CoordinateTransformer<double> transformer;
    transformer.applyTranslation(VectorType {0, 1, 0});
    transformer.applyRotation(x, y);
    VectorType correct({-2, 1, 1});
    VectorType result(transformer.getTransformed(pt.array()));
    VectorType diff(result.array() - correct.array());
    std::cout << pt << ' ' << result << ' ' << correct << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, to_dest_pts)
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  PointType Ps[4] = {{1, 2, 3}, {2, 2, 3}, {1, 4, 3}, {1, 2, 6}};
  PointType Qs[4] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  primal::experimental::CoordinateTransformer<double> transformer;
  transformer.setByTerminusPts(Ps, Qs);
  for(int i = 0; i < 4; ++i)
  {
    PointType result = Ps[i];
    transformer.transform(result.array());
    const auto& correct = Qs[i];
    VectorType diff(result.array() - Qs[i].array());
    std::cout << Ps[i] << ' ' << result << ' ' << correct << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  };

  auto invTran = transformer.getInverse();
  for(int i = 0; i < 4; ++i)
  {
    PointType result = Qs[i];
    invTran.transform(result.array());
    const auto& correct = Ps[i];
    VectorType diff(result.array() - Ps[i].array());
    std::cout << Qs[i] << ' ' << result << ' ' << correct << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  };
}

//------------------------------------------------------------------------------
TEST(primal_coord_transform, inverse)
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  {
    // Simple shift
    PointType start({1, 2, 3});
    VectorType shift({10, 20, 30});
    primal::experimental::CoordinateTransformer<double> transformer;
    transformer.applyTranslation(shift);
    primal::experimental::CoordinateTransformer<double> inverse = transformer.getInverse();
    PointType changed = transformer.getTransformed(start);
    PointType undone = inverse.getTransformed(changed);
    VectorType diff(start.array() - undone.array());
    std::cout << start << ' ' << changed << ' ' << undone << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }

  {
    // Simple rotate
    PointType start({1, 2, 30});
    VectorType axis({1, .5, .2});
    double angle = M_PI / 180 * 10;
    primal::experimental::CoordinateTransformer<double> transformer;
    transformer.applyRotation(axis.unitVector(), angle);
    primal::experimental::CoordinateTransformer<double> inverse = transformer.getInverse();
    PointType changed = transformer.getTransformed(start);
    PointType undone = inverse.getTransformed(changed);
    VectorType diff(start.array() - undone.array());
    std::cout << start << ' ' << changed << ' ' << undone << ' ' << diff << std::endl;
    EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
  }

  {
    // A bunch of random translate + rotate transforms.
    const int n = 3;
    for(int i = 0; i < n; ++i)
    {
      PointType startPt;
      startPt[0] = axom::utilities::random_real(-10., 10.);
      startPt[1] = axom::utilities::random_real(-10., 10.);
      startPt[2] = axom::utilities::random_real(-10., 10.);
      VectorType shift({axom::utilities::random_real(-100., 100.),
                        axom::utilities::random_real(-100., 100.),
                        axom::utilities::random_real(-100., 100.)});
      VectorType axis({axom::utilities::random_real(-100., 100.),
                       axom::utilities::random_real(-100., 100.),
                       axom::utilities::random_real(-100., 100.)});
      double angle = axom::utilities::random_real(-2 * M_PI, 2 * M_PI);

      primal::experimental::CoordinateTransformer<double> transformer;
      transformer.applyTranslation(shift);
      transformer.applyRotation(axis.unitVector(), angle);
      primal::experimental::CoordinateTransformer<double> inverse = transformer.getInverse();
      PointType endPt = transformer.getTransformed(startPt);
      PointType result = inverse.getTransformed(endPt);
      VectorType diff(result.array() - startPt.array());
      std::cout << startPt << ' ' << endPt << ' ' << result << ' ' << diff << std::endl;
      EXPECT_TRUE(axom::utilities::isNearlyEqual(diff.norm(), 1e-12));
    }
  }
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
