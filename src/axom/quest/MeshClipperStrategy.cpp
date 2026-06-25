// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/klee/GeometryOperators.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{
namespace internal
{
/*!
 * \brief Implementation of a GeometryOperatorVisitor for processing klee shape operators
 *
 * This class extracts the matrix form of supported operators.
 * For other operators, it sets the valid flag to false.
 * To use, check the \a isValid() function after visiting and then call the \a getMatrix() function.
 */
class AffineMatrixVisitor : public klee::GeometryOperatorVisitor
{
public:
  AffineMatrixVisitor() : m_matrix(4, 4) { }

  void visit(const klee::Translation& translation) override
  {
    m_matrix = translation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Rotation& rotation) override
  {
    m_matrix = rotation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Scale& scale) override
  {
    m_matrix = scale.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::UnitConverter& converter) override
  {
    m_matrix = converter.toMatrix();
    m_isValid = true;
  }

  void visit(const klee::CompositeOperator&) override
  {
    SLIC_WARNING("CompositeOperator not supported for Shaper query");
    m_isValid = false;
  }
  void visit(const klee::SliceOperator&) override
  {
    SLIC_WARNING("SliceOperator not yet supported for Shaper query");
    m_isValid = false;
  }

  const numerics::Matrix<double>& getMatrix() const { return m_matrix; }
  bool isValid() const { return m_isValid; }

private:
  bool m_isValid {false};
  numerics::Matrix<double> m_matrix;
};

}  // end namespace internal

MeshClipperStrategy::MeshClipperStrategy(const klee::Geometry& kGeom)
  : m_info(kGeom.asHierarchy())
  , m_extTrans(kGeom.getTransform())
{ }

const std::string& MeshClipperStrategy::name() const
{
  static const std::string n = "UNNAMED";
  return n;
}

const axom::primal::BoundingBox<double, 2>& MeshClipperStrategy::getBoundingBox2D() const
{
  static const axom::primal::BoundingBox<double, 2> invalidBb2d;
  return invalidBb2d;
}

const axom::primal::BoundingBox<double, 3>& MeshClipperStrategy::getBoundingBox3D() const
{
  static const axom::primal::BoundingBox<double, 3> invalidBb3d;
  return invalidBb3d;
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
