// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/klee/AffineMatrixVisitor.hpp"

namespace axom::klee
{

AffineMatrixVisitor::AffineMatrixVisitor()
  : GeometryOperatorVisitor()
  , m_isValid(false)
  , m_matrix(4, 4)
{ }

void AffineMatrixVisitor::visit(const klee::Translation& translation)
{
  m_matrix = translation.toMatrix();
  m_isValid = true;
}
void AffineMatrixVisitor::visit(const klee::Rotation& rotation)
{
  m_matrix = rotation.toMatrix();
  m_isValid = true;
}
void AffineMatrixVisitor::visit(const klee::Scale& scale)
{
  m_matrix = scale.toMatrix();
  m_isValid = true;
}
void AffineMatrixVisitor::visit(const klee::UnitConverter& converter)
{
  m_matrix = converter.toMatrix();
  m_isValid = true;
}

void AffineMatrixVisitor::visit(const klee::CompositeOperator&)
{
  SLIC_WARNING_ROOT("CompositeOperator not supported for Shaper query");
  m_isValid = false;
}
void AffineMatrixVisitor::visit(const klee::SliceOperator&)
{
  SLIC_WARNING_ROOT("SliceOperator not yet supported for Shaper query");
  m_isValid = false;
}

}  // end namespace axom::klee
