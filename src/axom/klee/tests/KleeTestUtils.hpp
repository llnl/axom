// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_TESTUTILS_HPP_
#define AXOM_KLEE_TESTUTILS_HPP_

#include "axom/klee/GeometryOperators.hpp"

#include "axom/core/numerics/Matrix.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "gmock/gmock.h"

#include <array>

namespace axom
{
namespace klee
{
namespace test
{
/**
 * Create a 3D affine transformation matrix
 *
 * \param values the values of the matrix, in row-major order
 * \return the affine transformation matrix represented by the rows
 */
numerics::Matrix<double> affine(const std::array<std::array<double, 4>, 3> &values);

class MockOperator : public GeometryOperator
{
public:
  using GeometryOperator::GeometryOperator;
  MOCK_METHOD(TransformableGeometryProperties, getEndProperties, (), (const));
  MOCK_METHOD(void, accept, (GeometryOperatorVisitor &), (const));
  TransformableGeometryProperties getBaseEndProperties() const
  {
    return GeometryOperator::getEndProperties();
  }
};

}  // namespace test
}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_TESTUTILS_HPP_
