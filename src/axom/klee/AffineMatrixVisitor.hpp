// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_AFFINE_MATRIX_VISITOR_HPP_
#define AXOM_KLEE_AFFINE_MATRIX_VISITOR_HPP_
#include "axom/klee/GeometryOperators.hpp"

namespace axom::klee
{

/*!
 * \brief Implementation of a GeometryOperatorVisitor for processing klee shape operators
 *
 * This class extracts the matrix form of supported operators and marks the operator as unvalid otherwise
 * To use, check the \a isValid() function after visiting and then call the \a getMatrix() function.
 */
class AffineMatrixVisitor : public GeometryOperatorVisitor
{
public:
  AffineMatrixVisitor();

  void visit(const klee::Translation& translation) override;
  void visit(const klee::Rotation& rotation) override;
  void visit(const klee::Scale& scale) override;
  void visit(const klee::UnitConverter& converter) override;

  void visit(const klee::CompositeOperator&) override;
  void visit(const klee::SliceOperator&) override;

  const numerics::Matrix<double>& getMatrix() const { return m_matrix; }
  bool isValid() const { return m_isValid; }

private:
  bool m_isValid;
  numerics::Matrix<double> m_matrix;
};

}  // end namespace axom::klee

#endif
