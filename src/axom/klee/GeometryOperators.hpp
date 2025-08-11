// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_GEOMETRYOPERATOR_HPP_
#define AXOM_KLEE_GEOMETRYOPERATOR_HPP_

#include "axom/config.hpp"
#include "axom/core/numerics/Matrix.hpp"

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Geometry.hpp"
#include "axom/klee/Units.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <memory>
#include <vector>

namespace axom
{
namespace klee
{
class GeometryOperatorVisitor;

/**
 * \brief A GeometryOperator describes an operation to perform on the Geometry of a Shape
 *
 * There is a subclass of GeometryOperator for each operator defined in
 * the format specification. You can figure out which one you're working
 * with via the accept() method. This follows the standard visitor pattern.
 */
class GeometryOperator
{
public:
  /**
   * Create an operator with the given start properties
   * \param startProperties the properties before the operator is applied
   */
  explicit GeometryOperator(const TransformableGeometryProperties &startProperties);

  virtual ~GeometryOperator() = default;

  /**
   * Get the properties that the operator expects to start in
   *
   * \return the properties which must be true before this operator is applied
   */
  const TransformableGeometryProperties &getStartProperties() const { return m_startProperties; }

  /**
   * Get the properties after this operator is applied
   *
   * \return the properties which are true after this operator is applied
   */
  virtual TransformableGeometryProperties getEndProperties() const
  {
    // By default, end properties are the same as start properties, so this is correct.
    return m_startProperties;
  }

  /**
   * Accept the given visitor. The appropriate visit() method will
   * be called on the visitor based on the run-time type of this operator.
   *
   * \param visitor the visitor to accept.
   */
  virtual void accept(GeometryOperatorVisitor &visitor) const = 0;

private:
  TransformableGeometryProperties m_startProperties;
};

/**
 * A MatrixOperator is a type of GeometryOperator whose operation can be
 * expressed as a 4x4 affine transformation matrix.
 */
class MatrixOperator : public GeometryOperator
{
public:
  using GeometryOperator::GeometryOperator;

  /**
   * Convert this operator to its matrix representation.
   *
   * \return a 4x4 affine transformation matrix
   */
  virtual numerics::Matrix<double> toMatrix() const = 0;
};

/// A CompositeOperator contains a list of GeometryOperators that are applied in order
class CompositeOperator : public GeometryOperator
{
public:
  using GeometryOperator::GeometryOperator;

  using OpPtr = std::shared_ptr<const GeometryOperator>;

  void accept(GeometryOperatorVisitor &visitor) const override;

  /**
   * Add the given operator to the end of the list of operators in this composite.
   *
   * \param op the operator to add
   */
  void addOperator(const OpPtr &op);

  /**
   * Get a list of all the operators. They should be applied in order.
   *
   * \return the list of operators
   */
  const std::vector<OpPtr> &getOperators() const { return m_operators; }

  TransformableGeometryProperties getEndProperties() const override;

private:
  std::vector<OpPtr> m_operators;
};

/// A Translation is a GeometryOperator which translates points.
class Translation : public MatrixOperator
{
public:
  /**
   * Create a Translation.
   *
   * \param offset the amount by which to offset points
   * \param startProperties the initial properties, as in the parent class. 
   * If the number of dimensions is 2, the 3rd entry in the offset should be zero, but this is not checked.
   */
  Translation(const primal::Vector3D &offset, const TransformableGeometryProperties &startProperties);

  /**
   * Get the amount by which to offset points.
   *
   * \return a vector by which points should be offset
   */
  const primal::Vector3D &getOffset() const { return m_offset; }

  numerics::Matrix<double> toMatrix() const override;

  void accept(GeometryOperatorVisitor &visitor) const override;

private:
  primal::Vector3D m_offset;
};

/// A Rotation is a GeometryOperator which rotates points about a given axis.
class Rotation : public MatrixOperator
{
public:
  /**
   * Create a Rotation.
   *
   * \param angle the angle, in degrees, by which to rotate. Rotations are counter-clockwise.
   * \param center the center of rotation
   * \param axis the axis about which to rotate points
   * \param startProperties the initial properties, as in the parent class. 
   * If the number of dimensions is 2, the axis should be [0, 0, 1], but this is not checked.
   */
  Rotation(double angle,
           const primal::Point3D &center,
           const primal::Vector3D &axis,
           const TransformableGeometryProperties &startProperties);

  /**
   * Get the angle of rotation.
   *
   * \return the amount by which to rotate in degrees.
   */
  double getAngle() const { return m_angle; }

  /**
   * Get the center of rotation.
   *
   * \return the point about which to rotate in 2D, and in 3D, the point
   * which defines the axis of rotation along with getAxis().
   */
  const primal::Point3D &getCenter() const { return m_center; }

  /**
   * The direction of the axis of rotation.
   *
   * \return the vector, which when combined with the center, defines the axis of rotation.
   */
  const primal::Vector3D &getAxis() const { return m_axis; }

  numerics::Matrix<double> toMatrix() const override;

  void accept(GeometryOperatorVisitor &visitor) const override;

private:
  double m_angle;
  primal::Point3D m_center;
  primal::Vector3D m_axis;
};

/// A GeometryOperator for scaling shapes
class Scale : public MatrixOperator
{
public:
  /**
   * Create a new Scale operator.
   *
   * \param xFactor the amount by which to scale in the x direction
   * \param yFactor the amount by which to scale in the y direction
   * \param zFactor the amount by which to scale in the z direction
   * \param startProperties the initial properties, as in the parent class. 
   * If the number of dimensions is 2, zFactor should be 1.0, but this is not enforced.
   */
  Scale(double xFactor,
        double yFactor,
        double zFactor,
        const TransformableGeometryProperties &startProperties);

  /**
   * Get the scale factor in the x direction.
   *
   * \return the x scale factor
   */
  double getXFactor() const { return m_xFactor; }

  /**
   * Get the scale factor in the y direction.
   *
   * \return the y scale factor
   */
  double getYFactor() const { return m_yFactor; }

  /**
   * Get the scale factor in the z direction.
   *
   * \return the z scale factor
   */
  double getZFactor() const { return m_zFactor; }

  numerics::Matrix<double> toMatrix() const override;

  void accept(GeometryOperatorVisitor &visitor) const override;

private:
  double m_xFactor;
  double m_yFactor;
  double m_zFactor;
};

/// An operator for converting units
class UnitConverter : public MatrixOperator
{
public:
  /**
   * Convert from the units specified in the start properties to the given end units
   * \param endUnits the units at the end of the operation
   * \param startProperties the properties before the operation
   */
  UnitConverter(LengthUnit endUnits, const TransformableGeometryProperties &startProperties);

  TransformableGeometryProperties getEndProperties() const override;

  numerics::Matrix<double> toMatrix() const override;

  void accept(GeometryOperatorVisitor &visitor) const override;

  /**
   * Get the conversion factor used to convert from the start units to the end units
   * \return the unit conversion factor
   */
  double getConversionFactor() const;

private:
  LengthUnit m_endUnits;
};

/**
 * A SliceOperator takes a 3D shape and makes it 2D by defining a cut plane.
 * This plane is augmented with an origin and an "up" direction to orient the slice to the x-y plane.
 */
class SliceOperator : public MatrixOperator
{
public:
  /**
   * Create a new Slice.
   *
   * \param origin the origin of the coordinate system
   * \param normal a vector normal to the slice plane. Cannot be a zero vector.
   * \param up the direction of the positive Y axis. Must be normal to the "normal" vector.
   * \param startProperties the initial properties, as in the parent class. 
   * The number of dimensions should be 3, though this is not checked.
   */
  SliceOperator(const primal::Point3D &origin,
                const primal::Vector3D &normal,
                const primal::Vector3D &up,
                const TransformableGeometryProperties &startProperties);

  /**
   * Get the origin of the coordinate system.
   *
   * \return the system's origin
   */
  const primal::Point3D &getOrigin() const { return m_origin; }

  /**
   * Get a vector normal to the slice plane.
   *
   * \return a vector normal to the slice plane
   */
  const primal::Vector3D &getNormal() const { return m_normal; }

  /**
   * Get a vector in the direction of the positive Y axis.
   *
   * \return the direction of the positive Y axis
   */
  const primal::Vector3D &getUp() const { return m_up; }

  numerics::Matrix<double> toMatrix() const override;

  void accept(GeometryOperatorVisitor &visitor) const override;

  TransformableGeometryProperties getEndProperties() const override;

private:
  numerics::Matrix<double> createRotation() const;
  numerics::Matrix<double> createTranslationToOrigin() const;
  primal::Vector3D calculateRightVector() const;

  primal::Point3D m_origin;
  primal::Vector3D m_normal;
  primal::Vector3D m_up;
};

/**
 * A GeometryOperatorVisitor defines visitor for interacting with
 * instances of GeometryOperator. It defines a "visit()" method for each
 * type of operator that a user can specify in the input file, as well
 * as one for CompositeOperator to work with a whole list of operators.
 */
class GeometryOperatorVisitor
{
public:
  virtual ~GeometryOperatorVisitor() = default;

  virtual void visit(const Translation &translation) = 0;

  virtual void visit(const Rotation &rotation) = 0;

  virtual void visit(const Scale &scale) = 0;

  virtual void visit(const UnitConverter &converter) = 0;

  virtual void visit(const CompositeOperator &composite) = 0;

  virtual void visit(const SliceOperator &slice) = 0;
};

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_GEOMETRYOPERATOR_HPP_
