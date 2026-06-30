// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"
#include "axom/klee/tests/KleeMatchers.hpp"
#include "axom/klee/tests/KleeTestUtils.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>

namespace axom
{
namespace klee
{
using test::AlmostEqPoint;
using test::MockOperator;
using ::testing::HasSubstr;
using ::testing::Return;

TEST(GeometryTest, dimensions_noOperators)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::mils};
  Geometry geometry {startProperties, "test format", "test path", nullptr};
  EXPECT_EQ(startProperties, geometry.getStartProperties());
  EXPECT_EQ(startProperties, geometry.getEndProperties());

  EXPECT_TRUE(geometry.hasGeometry());

  EXPECT_EQ(startProperties.dimensions, geometry.getInputDimensions());
  EXPECT_EQ(startProperties.dimensions, geometry.getOutputDimensions());
}

TEST(GeometryTest, dimensions_dimensionPreservingOperator)
{
  TransformableGeometryProperties startProperties {Dimensions::Two, LengthUnit::mils};
  TransformableGeometryProperties endProperties {Dimensions::Three, LengthUnit::cm};
  auto mockOperator = std::make_shared<MockOperator>(startProperties);
  Geometry geometry {startProperties, "test format", "test path", mockOperator};

  ON_CALL(*mockOperator, getEndProperties()).WillByDefault(Return(endProperties));
  EXPECT_CALL(*mockOperator, getEndProperties());

  EXPECT_EQ(startProperties.dimensions, geometry.getInputDimensions());
  EXPECT_EQ(endProperties.dimensions, geometry.getOutputDimensions());
}

TEST(GeometryTest, emptyPath)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::mils};
  Geometry geometry {startProperties, "none", "", nullptr};

  EXPECT_FALSE(geometry.hasGeometry());
}

TEST(GeometryTest, hasNonAffineOperators)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::cm};
  Geometry affineGeometry {
    startProperties,
    "test format",
    "test path",
    std::make_shared<Translation>(primal::Vector3D {1., 2., 3.}, startProperties)};
  EXPECT_FALSE(affineGeometry.hasNonAffineOperators());

  auto composite = std::make_shared<CompositeOperator>(startProperties);
  composite->addOperator(
    std::make_shared<Translation>(primal::Vector3D {1., 0., 0.}, startProperties));
  composite->addOperator(std::make_shared<PointTransform>(
    [](const Geometry::Point3D& p) { return Geometry::Point3D {p[0] * 2., p[1], p[2]}; },
    startProperties,
    "operators/2/transform",
    "Error evaluating callback for 'transform'"));
  Geometry nonAffineGeometry {startProperties, "test format", "test path", composite};
  EXPECT_TRUE(nonAffineGeometry.hasNonAffineOperators());
}

TEST(GeometryTest, applyTransformPreservesOperatorOrder)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::cm};
  auto composite = std::make_shared<CompositeOperator>(startProperties);
  composite->addOperator(
    std::make_shared<Translation>(primal::Vector3D {1., 0., 0.}, startProperties));
  composite->addOperator(std::make_shared<PointTransform>(
    [](const Geometry::Point3D& p) { return Geometry::Point3D {p[0] * 2., p[1] + 3., p[2]}; },
    startProperties,
    "operators/2/transform",
    "Error evaluating callback for 'transform'"));
  Geometry geometry {startProperties, "test format", "test path", composite};

  EXPECT_THAT(geometry.applyTransform({1., 2., 3.}), AlmostEqPoint(Geometry::Point3D {4., 5., 3.}));
}

TEST(GeometryTest, getTransformRejectsNonAffineOperators)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::cm};
  auto transform = std::make_shared<PointTransform>([](const Geometry::Point3D& p) { return p; },
                                                    startProperties,
                                                    "operators/1/transform",
                                                    "Error evaluating callback for 'transform'");
  Geometry geometry {startProperties, "test format", "test path", transform};

  try
  {
    geometry.getTransform();
    FAIL() << "Expected non-affine getTransform() to throw";
  }
  catch(const KleeError& err)
  {
    EXPECT_THAT(err.what(), HasSubstr("Cannot convert geometry to matrix"));
    EXPECT_THAT(err.what(), HasSubstr("transform"));
  }
}

}  // namespace klee
}  // namespace axom
