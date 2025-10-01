// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal/geometry/Cone.hpp"
#include "axom/primal/geometry/Point.hpp"

#include <math.h>

#include "axom/slic.hpp"
#include "gtest/gtest.h"

namespace primal = axom::primal;

constexpr int NDIMS = 3;

using ConeType = primal::Cone<double, NDIMS>;
using PointType = primal::Point<double, NDIMS>;
using VectorType = primal::Vector<double, NDIMS>;

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
double cone_volume(double baseRad, double topRad, double len)
{
  return M_PI / 3 * len * (baseRad * baseRad + baseRad * topRad + topRad * topRad);
}

//------------------------------------------------------------------------------
TEST(primal_cone, default_constructor)
{
  const double defaultBaseRad = 1.0;
  const double defaultTopRad = 0.0;
  const double defaultLength = 1.0;
  const PointType defaultBaseCenter(0.0, NDIMS);
  const VectorType defaultDir({1.0, 0.0, 0.0});
  const double defaultVol = cone_volume(defaultBaseRad, defaultTopRad, defaultLength);
  ConeType defaultCone;
  EXPECT_EQ(defaultCone.getBaseRadius(), defaultBaseRad);
  EXPECT_EQ(defaultCone.getTopRadius(), defaultTopRad);
  EXPECT_EQ(defaultCone.getLength(), defaultLength);
  EXPECT_EQ(defaultCone.volume(), defaultVol);
  EXPECT_EQ(defaultCone.getBaseCenter(), defaultBaseCenter);
  EXPECT_EQ(defaultCone.getDirection(), defaultDir);

  const double a = .643245324;
  const double b = 1 - a;
  EXPECT_EQ(defaultCone.getRadiusAt(b * defaultLength), a * defaultBaseRad + b * defaultTopRad);
}

//------------------------------------------------------------------------------
TEST(primal_cone, full_constructor)
{
  const double aBaseRad = 3.5;
  const double aTopRad = 0.0;
  const double aLength = 1.0;
  const PointType aBaseCenter({1.0, 2.0, 3.0});
  const VectorType aDir({0.5, 0.6, 0.7});
  const double aVol = cone_volume(aBaseRad, aTopRad, aLength);
  ConeType coneA(aBaseRad, aTopRad, aLength, aDir, aBaseCenter);
  EXPECT_EQ(coneA.getBaseRadius(), aBaseRad);
  EXPECT_EQ(coneA.getTopRadius(), aTopRad);
  EXPECT_EQ(coneA.getLength(), aLength);
  EXPECT_EQ(coneA.volume(), aVol);
  EXPECT_EQ(coneA.getBaseCenter(), aBaseCenter);
  EXPECT_EQ(coneA.getDirection(), aDir.unitVector());

  const double a = .643245324;
  const double b = 1 - a;
  EXPECT_EQ(coneA.getRadiusAt(b * aLength), a * aBaseRad + b * aTopRad);
}

//------------------------------------------------------------------------------
TEST(primal_cone, copy_constructor)
{
  const double aBaseRad = 3.5;
  const double aTopRad = 0.0;
  const double aLength = 1.0;
  const PointType aBaseCenter({1.0, 2.0, 3.0});
  const VectorType aDir({0.5, 0.6, 0.7});
  ConeType coneA(aBaseRad, aTopRad, aLength, aDir, aBaseCenter);
  ConeType coneB(coneA);

  EXPECT_EQ(coneA.getBaseRadius(), coneB.getBaseRadius());
  EXPECT_EQ(coneA.getTopRadius(), coneB.getTopRadius());
  EXPECT_EQ(coneA.getLength(), coneB.getLength());
  EXPECT_EQ(coneA.getBaseCenter(), coneB.getBaseCenter());
  EXPECT_EQ(coneA.volume(), coneB.volume());
}

//------------------------------------------------------------------------------
TEST(primal_cone, assignment_operator)
{
  const double aBaseRad = 3.5;
  const double aTopRad = 0.0;
  const double aLength = 1.0;
  const PointType aBaseCenter({1.0, 2.0, 3.0});
  const VectorType aDir({0.5, 0.6, 0.7});
  ConeType coneA(aBaseRad, aTopRad, aLength, aDir, aBaseCenter);
  ConeType coneB;
  coneB = coneA;

  EXPECT_EQ(coneA.getBaseRadius(), coneB.getBaseRadius());
  EXPECT_EQ(coneA.getTopRadius(), coneB.getTopRadius());
  EXPECT_EQ(coneA.getLength(), coneB.getLength());
  EXPECT_EQ(coneA.getBaseCenter(), coneB.getBaseCenter());
  EXPECT_EQ(coneA.volume(), coneB.volume());
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
