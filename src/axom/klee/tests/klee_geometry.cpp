// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/tests/KleeMatchers.hpp"
#include "axom/klee/tests/KleeTestUtils.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>

namespace axom
{
namespace klee
{
using test::AlmostEqMatrix;
using test::MockOperator;
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

TEST(GeometryTest, getTransform_sliceOperator)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::cm};
  auto slice = std::make_shared<SliceOperator>(primal::Point3D {10, 20, 30},
                                               primal::Vector3D {1, 4, 8},
                                               primal::Vector3D {-4, 1, 0},
                                               startProperties);
  Geometry geometry {startProperties, "test format", "test path", slice};

  EXPECT_THAT(geometry.getTransform(), AlmostEqMatrix(slice->toMatrix()));
}

}  // namespace klee
}  // namespace axom
