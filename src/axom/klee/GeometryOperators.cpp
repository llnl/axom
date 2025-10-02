// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/numerics/matvecops.hpp"

#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/Units.hpp"

#include <cmath>
#include <stdexcept>

namespace axom
{
namespace klee
{
GeometryOperator::GeometryOperator(const TransformableGeometryProperties &startProperties)
  : m_startProperties(startProperties)
{ }

void CompositeOperator::accept(GeometryOperatorVisitor &visitor) const { visitor.visit(*this); }

void CompositeOperator::addOperator(const OpPtr &op)
{
  if(getEndProperties() != op->getStartProperties())
  {
    throw std::invalid_argument("Start and end properties don't match");
  }
  m_operators.emplace_back(op);
}

TransformableGeometryProperties CompositeOperator::getEndProperties() const
{
  if(m_operators.empty())
  {
    return GeometryOperator::getEndProperties();
  }
  return (*m_operators.rbegin())->getEndProperties();
}

Translation::Translation(const primal::Vector3D &offset,
                         const TransformableGeometryProperties &startProperties)
  : MatrixOperator {startProperties}
  , m_offset {offset}
{ }

numerics::Matrix<double> Translation::toMatrix() const
{
  auto transformation = numerics::Matrix<double>::identity(4);
  for(int i = 0; i < 3; ++i)
  {
    transformation(i, 3) = m_offset[i];
  }
  return transformation;
}

void Translation::accept(GeometryOperatorVisitor &visitor) const { visitor.visit(*this); }

Rotation::Rotation(double angle,
                   const primal::Point3D &center,
                   const primal::Vector3D &axis,
                   const TransformableGeometryProperties &startProperties)
  : MatrixOperator {startProperties}
  , m_angle {angle}
  , m_center {center}
  , m_axis {axis}
{ }

numerics::Matrix<double> Rotation::toMatrix() const
{
  // Use the matrix listed at
  // https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
  // in the section "Rotation matrix from angle and axis"
  double angleInRadians = m_angle * M_PI / 180;
  double sinT = std::sin(angleInRadians);
  double cosT = std::cos(angleInRadians);

  auto unitAxis = m_axis.unitVector();

  // Match names in wikipedia article for ease of implementation and
  // to keep expressions shorter.
  double ux = unitAxis[0];
  double uy = unitAxis[1];
  double uz = unitAxis[2];

  numerics::Matrix<double> transformation = numerics::Matrix<double>::zeros(4, 4);
  transformation(0, 0) = cosT + ux * ux * (1 - cosT);
  transformation(0, 1) = ux * uy * (1 - cosT) - uz * sinT;
  transformation(0, 2) = ux * uz * (1 - cosT) + uy * sinT;
  transformation(1, 0) = uy * ux * (1 - cosT) + uz * sinT;
  transformation(1, 1) = cosT + uy * uy * (1 - cosT);
  transformation(1, 2) = uy * uz * (1 - cosT) - ux * sinT;
  transformation(2, 0) = uz * ux * (1 - cosT) - uy * sinT;
  transformation(2, 1) = uz * uy * (1 - cosT) + ux * sinT;
  transformation(2, 2) = cosT + uz * uz * (1 - cosT);
  transformation(3, 3) = 1;

  for(int i = 0; i < 3; ++i)
  {
    transformation(i, 3) = m_center[i];
    for(int j = 0; j < 3; ++j)
    {
      transformation(i, 3) -= m_center[j] * transformation(i, j);
    }
  }

  return transformation;
}

void Rotation::accept(GeometryOperatorVisitor &visitor) const { visitor.visit(*this); }

Scale::Scale(double xFactor,
             double yFactor,
             double zFactor,
             const TransformableGeometryProperties &startProperties)
  : MatrixOperator {startProperties}
  , m_xFactor {xFactor}
  , m_yFactor {yFactor}
  , m_zFactor {zFactor}
{ }

numerics::Matrix<double> Scale::toMatrix() const
{
  auto transformation = numerics::Matrix<double>::zeros(4, 4);
  transformation(0, 0) = m_xFactor;
  transformation(1, 1) = m_yFactor;
  transformation(2, 2) = m_zFactor;
  transformation(3, 3) = 1;
  return transformation;
}

void Scale::accept(GeometryOperatorVisitor &visitor) const { visitor.visit(*this); }

UnitConverter::UnitConverter(LengthUnit endUnits,
                             const TransformableGeometryProperties &startProperties)
  : MatrixOperator {startProperties}
  , m_endUnits {endUnits}
{ }

TransformableGeometryProperties UnitConverter::getEndProperties() const
{
  return TransformableGeometryProperties {getStartProperties().dimensions, m_endUnits};
};

numerics::Matrix<double> UnitConverter::toMatrix() const
{
  double factor = getConversionFactor();
  Scale scale(factor, factor, factor, getStartProperties());
  return scale.toMatrix();
}

void UnitConverter::accept(GeometryOperatorVisitor &visitor) const { visitor.visit(*this); };

double UnitConverter::getConversionFactor() const
{
  return klee::getConversionFactor(getStartProperties().units, m_endUnits);
};

SliceOperator::SliceOperator(const primal::Point3D &origin,
                             const primal::Vector3D &normal,
                             const primal::Vector3D &up,
                             const TransformableGeometryProperties &startProperties)
  : MatrixOperator {startProperties}
  , m_origin {origin}
  , m_normal {normal}
  , m_up {up}
{ }

numerics::Matrix<double> SliceOperator::toMatrix() const
{
  auto rotation = createRotation();
  auto translation = createTranslationToOrigin();

  numerics::Matrix<double> transformation = numerics::Matrix<double>::identity(4);
  numerics::matrix_multiply(rotation, translation, transformation);

  return transformation;
}

numerics::Matrix<double> SliceOperator::createRotation() const
{
  // The basic idea of our implementation is express the coordinate system of
  // the slice as a normal basis, and then do a change of basis to the
  // standard basis.
  auto unitNormal = m_normal.unitVector();
  auto unitUp = m_up.unitVector();
  auto unitRight = calculateRightVector();

  numerics::Matrix<double> rotation = numerics::Matrix<double>::identity(4);

  // Since our matrix is orthogonal, we can just transpose it to get the
  // inverse. We just do it in-place since we don't already have it set
  // up as a matrix.
  for(int i = 0; i < 3; ++i)
  {
    rotation(0, i) = unitRight[i];
    rotation(1, i) = unitUp[i];
    rotation(2, i) = unitNormal[i];
  }
  return rotation;
}

numerics::Matrix<double> SliceOperator::createTranslationToOrigin() const
{
  numerics::Matrix<double> translation = numerics::Matrix<double>::identity(4);
  for(int i = 0; i < 3; ++i)
  {
    translation(i, 3) = -1 * m_origin[i];
  }

  return translation;
}

primal::Vector3D SliceOperator::calculateRightVector() const
{
  Rotation rotation {270, primal::Point3D {0.0}, m_normal.unitVector(), getStartProperties()};
  primal::Vector<double, 4> unitRightAffine;
  auto unitUp = m_up.unitVector();
  primal::Vector<double, 4> upAffine {unitUp.data(), 3};
  numerics::matrix_vector_multiply(rotation.toMatrix(), upAffine.data(), unitRightAffine.data());
  return primal::Vector3D {unitRightAffine.data()};
}

void SliceOperator::accept(GeometryOperatorVisitor &visitor) const { visitor.visit(*this); }

TransformableGeometryProperties SliceOperator::getEndProperties() const
{
  TransformableGeometryProperties result = getStartProperties();
  result.dimensions = Dimensions::Two;
  return result;
}

}  // namespace klee
}  // namespace axom
