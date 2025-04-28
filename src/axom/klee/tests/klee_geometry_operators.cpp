// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/klee/GeometryOperators.hpp"

#include "axom/klee/tests/KleeMatchers.hpp"
#include "axom/klee/tests/KleeTestUtils.hpp"

#include <math.h>
#include <stdexcept>
#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace axom
{
namespace klee
{
namespace
{
using test::affine;
using test::AlmostEqMatrix;
using test::AlmostEqPoint;
using test::AlmostEqVector;
using test::MockOperator;

using ::testing::AtLeast;
using ::testing::ElementsAre;
using ::testing::Invoke;
using ::testing::IsEmpty;
using ::testing::Matcher;
using ::testing::Ref;
using ::testing::Return;

using primal::Point3D;
using primal::Vector3D;

template <typename ColumnVector>
ColumnVector operator*(const numerics::Matrix<double> &matrix, const ColumnVector &rhs)
{
  if(matrix.getNumRows() != matrix.getNumRows() || matrix.getNumRows() != rhs.dimension())
  {
    throw std::logic_error("Can't multiply entities of this size");
  }
  ColumnVector result;
  matrix_vector_multiply(matrix, rhs.data(), result.data());
  return result;
}

primal::Vector<double, 4> affineVec(const Vector3D &vec3d)
{
  primal::Vector<double, 4> vector {vec3d.data(), 3};
  vector[3] = 0;
  return vector;
}

primal::Point<double, 4> affinePoint(const Point3D &point3d)
{
  primal::Point<double, 4> point {point3d.data(), 3};
  point[3] = 1;
  return point;
}

Dimensions ALL_DIMS[] = {Dimensions::Two, Dimensions::Three};

class MockVisitor : public GeometryOperatorVisitor
{
public:
  MOCK_METHOD(void, visit, (const Translation &translation), (override));
  MOCK_METHOD(void, visit, (const Rotation &rotation), (override));
  MOCK_METHOD(void, visit, (const Scale &scale), (override));
  MOCK_METHOD(void, visit, (const UnitConverter &converter), (override));
  MOCK_METHOD(void, visit, (const CompositeOperator &op), (override));
  MOCK_METHOD(void, visit, (const SliceOperator &op), (override));
};

TEST(GeometryOperator, getProperties)
{
  TransformableGeometryProperties constructorProps {Dimensions::Three, LengthUnit::mm};
  MockOperator op {constructorProps};

  auto startProps = op.getStartProperties();
  EXPECT_EQ(constructorProps.dimensions, startProps.dimensions);
  EXPECT_EQ(constructorProps.units, startProps.units);

  ON_CALL(op, getEndProperties()).WillByDefault(Invoke(&op, &MockOperator::getBaseEndProperties));
  EXPECT_CALL(op, getEndProperties());
  auto endProperties = op.getEndProperties();
  EXPECT_EQ(constructorProps.dimensions, endProperties.dimensions);
  EXPECT_EQ(constructorProps.units, endProperties.units);
}

TEST(Translation, basics)
{
  Vector3D offset {10, 20, 30};
  Translation translation {offset, {Dimensions::Two, LengthUnit::cm}};
  EXPECT_THAT(translation.getOffset(), AlmostEqVector(offset));
}

TEST(Translation, toMatrix)
{
  for(Dimensions dims : ALL_DIMS)
  {
    Vector3D offset {10, 20, 30};
    Translation translation {offset, {dims, LengthUnit::cm}};
    EXPECT_THAT(translation.toMatrix(),
                AlmostEqMatrix(affine({{{1, 0, 0, 10}, {0, 1, 0, 20}, {0, 0, 1, 30}}})));
  }
}

TEST(Translation, accept)
{
  Translation translation {{10, 20, 30}, {Dimensions::Two, LengthUnit::cm}};
  MockVisitor visitor;
  EXPECT_CALL(visitor, visit(Matcher<const Translation &>(Ref(translation))));
  translation.accept(visitor);
}

TEST(Rotation, basics)
{
  Vector3D axis {10, 20, 30};
  Point3D center {40, 50, 60};
  double angle = 45;
  Rotation rotation {angle, center, axis, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_DOUBLE_EQ(angle, rotation.getAngle());
  EXPECT_THAT(rotation.getCenter(), AlmostEqPoint(center));
  EXPECT_THAT(rotation.getAxis(), AlmostEqVector(axis));
}

TEST(Rotation, rotate2d_with_center)
{
  Vector3D axis {0, 0, 1};
  Point3D center {10, 20, 0};
  double angle = 30;
  Rotation rotation {angle, center, axis, {Dimensions::Two, LengthUnit::cm}};

  double sin30 = 0.5;
  double cos30 = std::sqrt(3) / 2;

  double xOffset = 10 - 10 * cos30 + 20 * sin30;
  double yOffset = 20 - 10 * sin30 - 20 * cos30;

  EXPECT_THAT(rotation.toMatrix(),
              AlmostEqMatrix(
                affine({{{cos30, -sin30, 0, xOffset}, {sin30, cos30, 0, yOffset}, {0, 0, 1, 0}}})));

  // Sanity check to make sure things work with nice rotation since
  // the test and implementation use similar approaches (though the equations
  // were verified by hand).
  Rotation rotate90 {90, {10, 20, 0}, {0, 0, 1}, {Dimensions::Two, LengthUnit::cm}};

  EXPECT_THAT(rotate90.toMatrix() * affinePoint({15, 20, 0}),
              AlmostEqPoint(affinePoint({10, 25, 0})));
}

TEST(Rotation, rotate3d_axis_aligned)
{
  double sin30 = 0.5;
  double cos30 = std::sqrt(3) / 2;

  Point3D origin {0, 0, 0};

  Rotation rotatedAboutXAxis {30, origin, {1, 0, 0}, {Dimensions::Three, LengthUnit::cm}};

  EXPECT_THAT(rotatedAboutXAxis.toMatrix(),
              AlmostEqMatrix(affine({{{1, 0, 0, 0}, {0, cos30, -sin30, 0}, {0, sin30, cos30, 0}}})));

  Rotation rotatedAboutYAxis {30, origin, {0, 1, 0}, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(rotatedAboutYAxis.toMatrix(),
              AlmostEqMatrix(affine({{{cos30, 0, sin30, 0}, {0, 1, 0, 0}, {-sin30, 0, cos30, 0}}})));

  Rotation rotatedAboutZAxis {30, origin, {0, 0, 1}, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(rotatedAboutZAxis.toMatrix(),
              AlmostEqMatrix(affine({{{cos30, -sin30, 0, 0}, {sin30, cos30, 0, 0}, {0, 0, 1, 0}}})));
}

TEST(Rotation, rotate3d_with_center)
{
  Rotation rotation {90, {10, 20, 30}, {1, 1, 0}, {Dimensions::Three, LengthUnit::cm}};

  double halfRoot2 = std::sqrt(2) / 2;
  numerics::Matrix<double> expected =
    affine({{{0.5, 0.5, halfRoot2, 0}, {0.5, 0.5, -halfRoot2, 0}, {-halfRoot2, halfRoot2, 0, 0}}});
  expected(0, 3) = 10 - 10 * expected(0, 0) - 20 * expected(0, 1) - 30 * expected(0, 2);
  expected(1, 3) = 20 - 10 * expected(1, 0) - 20 * expected(1, 1) - 30 * expected(1, 2);
  expected(2, 3) = 30 - 10 * expected(2, 0) - 20 * expected(2, 1) - 30 * expected(2, 2);

  EXPECT_THAT(rotation.toMatrix(), AlmostEqMatrix(expected));

  EXPECT_THAT(rotation.toMatrix() * affinePoint({11, 20, 30}),
              AlmostEqPoint(affinePoint({10.5, 20.5, 30 - 1 / std::sqrt(2)})));

  // Use the pythagorean quadruple (1, 4, 8, 9) for easier manual checking
  rotation = Rotation {90, {10, 20, 30}, {1, -4, 8}, {Dimensions::Three, LengthUnit::cm}};

  expected = affine({{{1, -76, -28, 0}, {68, 16, -41, 0}, {44, -23, 64, 0}}});
  matrix_scalar_multiply(expected, 1 / 81.0);
  expected(3, 3) = 1;
  expected(0, 3) = 10 - 10 * expected(0, 0) - 20 * expected(0, 1) - 30 * expected(0, 2);
  expected(1, 3) = 20 - 10 * expected(1, 0) - 20 * expected(1, 1) - 30 * expected(1, 2);
  expected(2, 3) = 30 - 10 * expected(2, 0) - 20 * expected(2, 1) - 30 * expected(2, 2);
  EXPECT_THAT(rotation.toMatrix(), AlmostEqMatrix(expected));
}

TEST(Rotation, accept)
{
  Rotation rotation {90, {0, 0, 0}, {1, 2, 3}, {Dimensions::Three, LengthUnit::cm}};
  MockVisitor visitor;
  EXPECT_CALL(visitor, visit(Matcher<const Rotation &>(Ref(rotation))));
  rotation.accept(visitor);
}

TEST(Scale, basics)
{
  Scale scale {2, 3, 4, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_DOUBLE_EQ(2, scale.getXFactor());
  EXPECT_DOUBLE_EQ(3, scale.getYFactor());
  EXPECT_DOUBLE_EQ(4, scale.getZFactor());
}

TEST(Scale, toMatrix)
{
  Scale scale {2, 3, 4, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(scale.toMatrix(), AlmostEqMatrix(affine({{{2, 0, 0, 0}, {0, 3, 0, 0}, {0, 0, 4, 0}}})));
}

TEST(Scale, accept)
{
  Scale scale {1, 2, 3, {Dimensions::Three, LengthUnit::cm}};
  MockVisitor visitor;
  EXPECT_CALL(visitor, visit(Matcher<const Scale &>(Ref(scale))));
  scale.accept(visitor);
}

TEST(UnitConverter, basics)
{
  UnitConverter converter {LengthUnit::m, {Dimensions::Three, LengthUnit::cm}};
  TransformableGeometryProperties expectedEndProperties {Dimensions::Three, LengthUnit::m};
  EXPECT_EQ(expectedEndProperties, converter.getEndProperties());
  EXPECT_DOUBLE_EQ(0.01, converter.getConversionFactor());
}

TEST(UnitConverter, toMatrix)
{
  UnitConverter converter {LengthUnit::m, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(converter.toMatrix(),
              AlmostEqMatrix(affine({{{0.01, 0, 0, 0}, {0, 0.01, 0, 0}, {0, 0, 0.01, 0}}})));
}

TEST(UnitConverter, accept)
{
  UnitConverter converter {LengthUnit::m, {Dimensions::Three, LengthUnit::cm}};
  MockVisitor visitor;
  EXPECT_CALL(visitor, visit(Matcher<const UnitConverter &>(Ref(converter))));
  converter.accept(visitor);
}

TEST(CompositeOperator, empty)
{
  CompositeOperator op {{Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(op.getOperators(), IsEmpty());
}

TEST(CompositeOperator, addOperator)
{
  TransformableGeometryProperties startProps {Dimensions::Three, LengthUnit::cm};
  CompositeOperator compositeOp {startProps};
  EXPECT_EQ(startProps, compositeOp.getStartProperties());
  EXPECT_EQ(startProps, compositeOp.getEndProperties());

  auto startMock = std::make_shared<MockOperator>(startProps);
  TransformableGeometryProperties startMockEndProps {Dimensions::Three, LengthUnit::cm};
  ON_CALL(*startMock, getEndProperties()).WillByDefault(Return(startMockEndProps));
  compositeOp.addOperator(startMock);
  EXPECT_CALL(*startMock, getEndProperties()).Times(AtLeast(1));
  EXPECT_EQ(startProps, compositeOp.getStartProperties());
  EXPECT_EQ(startMockEndProps, compositeOp.getEndProperties());
  EXPECT_THAT(compositeOp.getOperators(), ElementsAre(startMock));

  TransformableGeometryProperties midMockEndProps {Dimensions::Three, LengthUnit::inches};
  auto midMock = std::make_shared<MockOperator>(startMockEndProps);
  ON_CALL(*midMock, getEndProperties()).WillByDefault(Return(midMockEndProps));
  compositeOp.addOperator(midMock);
  EXPECT_CALL(*midMock, getEndProperties()).Times(AtLeast(1));
  EXPECT_EQ(startProps, compositeOp.getStartProperties());
  EXPECT_EQ(midMockEndProps, compositeOp.getEndProperties());
  EXPECT_THAT(compositeOp.getOperators(), ElementsAre(startMock, midMock));

  TransformableGeometryProperties endMockEndProps {Dimensions::Two, LengthUnit::angstrom};
  auto endMock = std::make_shared<MockOperator>(midMockEndProps);
  ON_CALL(*endMock, getEndProperties()).WillByDefault(Return(endMockEndProps));
  compositeOp.addOperator(endMock);
  EXPECT_CALL(*endMock, getEndProperties()).Times(AtLeast(1));
  EXPECT_EQ(startProps, compositeOp.getStartProperties());
  EXPECT_EQ(endMockEndProps, compositeOp.getEndProperties());
  EXPECT_THAT(compositeOp.getOperators(), ElementsAre(startMock, midMock, endMock));
}

TEST(CompositeOperator, addOperator_doesNotMatchInitial)
{
  TransformableGeometryProperties startProps {Dimensions::Three, LengthUnit::cm};
  CompositeOperator compositeOp {startProps};

  auto startMock =
    std::make_shared<MockOperator>(TransformableGeometryProperties {Dimensions::Two, LengthUnit::cm});
  EXPECT_THROW(compositeOp.addOperator(startMock), std::invalid_argument);
}

TEST(CompositeOperator, addOperator_doesNotMatchLast)
{
  TransformableGeometryProperties startProps {Dimensions::Three, LengthUnit::cm};
  CompositeOperator compositeOp {startProps};

  auto startMock = std::make_shared<MockOperator>(
    TransformableGeometryProperties {Dimensions::Three, LengthUnit::cm});
  compositeOp.addOperator(startMock);

  auto midMock =
    std::make_shared<MockOperator>(TransformableGeometryProperties {Dimensions::Two, LengthUnit::mm});
  EXPECT_CALL(*startMock, getEndProperties());
  EXPECT_THROW(compositeOp.addOperator(midMock), std::invalid_argument);
}

TEST(CompositeOperator, accept)
{
  CompositeOperator composite {{Dimensions::Three, LengthUnit::cm}};
  MockVisitor visitor;
  EXPECT_CALL(visitor, visit(Matcher<const CompositeOperator &>(Ref(composite))));
  composite.accept(visitor);
}

TEST(Slice, basics)
{
  Point3D origin {1, 2, 3};
  Vector3D normal {4, 5, 6};
  Vector3D up {-5, 4, 0};
  SliceOperator slice {origin, normal, up, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(slice.getOrigin(), AlmostEqPoint(origin));
  EXPECT_THAT(slice.getNormal(), AlmostEqVector(normal));
  EXPECT_THAT(slice.getUp(), AlmostEqVector(up));

  TransformableGeometryProperties expectedEndProperties {Dimensions::Two, LengthUnit::cm};
  EXPECT_EQ(expectedEndProperties, slice.getEndProperties());
}

TEST(Slice, toMatrix_translateOnly)
{
  Point3D origin {10, 20, 30};
  Vector3D normal {0, 0, 1};
  Vector3D up {0, 1, 0};
  SliceOperator slice {origin, normal, up, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(slice.toMatrix(),
              AlmostEqMatrix(affine({{{1, 0, 0, -10}, {0, 1, 0, -20}, {0, 0, 1, -30}}})));
}

TEST(Slice, toMatrix_vectorScaleIgnored)
{
  Point3D origin {10, 20, 30};
  Vector3D normal {0, 0, 100};
  Vector3D up {0, 200, 0};
  SliceOperator slice {origin, normal, up, {Dimensions::Three, LengthUnit::cm}};
  EXPECT_THAT(slice.toMatrix(),
              AlmostEqMatrix(affine({{{1, 0, 0, -10}, {0, 1, 0, -20}, {0, 0, 1, -30}}})));
}

TEST(Slice, toMatrix_general)
{
  Point3D origin {10, 20, 30};
  Vector3D normal {1, 4, 8};
  Vector3D up {-4, 1, 0};
  SliceOperator slice {origin, normal, up, {Dimensions::Three, LengthUnit::cm}};

  auto sliceMatrix = slice.toMatrix();

  // Expected matrix calculated with Mathematica
  double root17 = std::sqrt(17);
  EXPECT_THAT(
    sliceMatrix,
    AlmostEqMatrix(affine({{{8 / (9 * root17), 32 / (9 * root17), -root17 / 9, -70 / (3 * root17)},
                            {-4 / root17, 1 / root17, 0, 20 / root17},
                            {1 / 9.0, 4 / 9.0, 8 / 9.0, -110 / 3.0}}})));

  // Extra checks to make sure we input things right both here and in
  // Mathematica. Verify expected operations on vectors and points.
  EXPECT_THAT(sliceMatrix * affinePoint(origin), AlmostEqPoint(affinePoint(primal::Point3D {0.0})));
  EXPECT_THAT(sliceMatrix * affineVec(normal.unitVector()), AlmostEqVector(affineVec({0, 0, 1})));
  EXPECT_THAT(sliceMatrix * affineVec(up.unitVector()), AlmostEqVector(affineVec({0, 1, 0})));
}

TEST(Slice, accept)
{
  SliceOperator slice {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {Dimensions::Three, LengthUnit::cm}};
  MockVisitor visitor;
  EXPECT_CALL(visitor, visit(Matcher<const SliceOperator &>(Ref(slice))));
  slice.accept(visitor);
}
}  // namespace
}  // namespace klee
}  // namespace axom
