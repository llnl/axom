// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file GregoryPatch.hpp
 *
 * \brief A Gregory patch primitive
 */

#ifndef AXOM_PRIMAL_GREGORY_PATCH_HPP_
#define AXOM_PRIMAL_GREGORY_PATCH_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include <ostream>
#include <math.h>

#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class GregoryPatch;

/*! \brief Overloaded output operator for Gregory Patches*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const GregoryPatch<T, NDIMS>& nPatch);

/*!
 * \class GregoryPatch
 *
 * \brief Represents a 3D Gregory patch defined by ...
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class GregoryPatch
{
public:
  // The number of control points for a bicubic Gregory patch is fixed:
  //  - 12 exterior control points
  //  - 2 interior control points for each of 4 boundary curves
  static constexpr int NPTS = 20;

  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  using CoordsVec = axom::StackArray<PointType, NPTS>;

  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;

  AXOM_STATIC_ASSERT_MSG((NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
                         "A Gregory Patch object may be defined in 1-, 2-, or 3-D");

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A Gregory Patch must be defined using an arithmetic type");

public:
  GregoryPatch() = default;

  explicit GregoryPatch(ArrayView<const PointType> controlPoints)
  {
    SLIC_ASSERT(controlPoints.size() == NPTS);
    for(int i = 0; i < NPTS; ++i)
    {
      m_controlPoints[i] = controlPoints[i];
    }
  }

  GregoryPatch(ArrayView<const PointType> nodePositions, ArrayView<const VectorType> nodeVectors)
  {
    // Store the position and orthogonal unit vector at each corner
    SLIC_ASSERT(nodePositions.size() == 4);
    SLIC_ASSERT(nodeVectors.size() == 4);

    axom::Array<VectorType> v(4);
    for(int i = 0; i < 4; ++i)
    {
      getCorner(i) = nodePositions[i];
      v[i] = nodeVectors[i].unitVector();
    }

    axom::Array<VectorType> c0(4), c2(4);
    axom::Array<VectorType> a0(4), a3(4);
    for(int i = 0; i < 4; ++i)
    {
      const int ip1 = (i + 1) % 4;
      const int im1 = (i + 3) % 4;

      const VectorType dx(getCorner(i), getCorner(ip1));

      c0[i] = (dx - dx.dot(v[i]) * v[i]) / 3;
      a0[i] = VectorType::cross_product(v[i], dx).unitVector();
      c2[i] = (dx - dx.dot(v[ip1]) * v[ip1]) / 3;
      a3[i] = VectorType::cross_product(v[ip1], dx).unitVector();

      // Other stuff for adjusted normals

      // Use Chiyokura algorithm to define the interior control points
      const PointType x1(getCorner(i).array() + c0[i].array());
      const PointType x2(getCorner(ip1).array() - c2[i].array());

      getBoundaryPoint(i, 1) = x1;
      getBoundaryPoint(i, 2) = x2;

      const VectorType c1(x1, x2);
      const VectorType b0 = -c2[im1];
      const VectorType b3 = c0[ip1];

      const double k0 = a0[i].dot(b0);
      const double k1 = a3[i].dot(b3);
      const double h0 = c0[i].dot(b0) / c0[i].dot(c0[i]);
      const double h1 = c2[i].dot(b3) / c2[i].dot(c2[i]);

      const VectorType b1 = ((k0 + k1) * a0[i] + k0 * a3[i] + 2.0 * h0 * c1 + h1 * c0[i]) / 3.0;
      const VectorType b2 = ((k0 + k1) * a3[i] + k1 * a0[i] + 2.0 * h1 * c1 + h0 * c2[i]) / 3.0;

      getTangent(i, 0) = PointType(x1.array() + b1.array());
      getTangent(i, 1) = PointType(x2.array() + b2.array());
    }
  }

  PointType& getCorner(int i) { return m_controlPoints[i]; }
  const PointType& getCorner(int i) const { return m_controlPoints[i]; }

  PointType& getTangent(int e, int t) { return m_controlPoints[12 + 2 * e + t]; }
  const PointType& getTangent(int e, int t) const { return m_controlPoints[12 + 2 * e + t]; }

  PointType& getBoundaryPoint(int e, int k) { return m_controlPoints[m_index_map[e][k]]; }
  const PointType& getBoundaryPoint(int e, int k) const
  {
    return m_controlPoints[m_index_map[e][k]];
  }

  // Evaluate the patch by constructing the equivalent Bezier patch with interior control nodes
  //  defined in terms of the tangent vectors and the evaluation parameters
  PointType evaluate(T u, T v) const
  {
    auto gregoryBlend = [&](const PointType& a, const PointType& b, T wa, T wb) -> PointType {
      const T denom = wa + wb;
      if(axom::utilities::isNearlyEqual(denom, T(0)))
      {
        return a;
      }
      return PointType((wa * a.array() + wb * b.array()) / denom);
    };

    const T um = T(1) - u;
    const T vm = T(1) - v;

    const PointType Q11 = gregoryBlend(getTangent(0, 0), getTangent(3, 1), u, v);
    const PointType Q21 = gregoryBlend(getTangent(0, 1), getTangent(1, 0), um, v);
    const PointType Q12 = gregoryBlend(getTangent(2, 1), getTangent(3, 0), u, vm);
    const PointType Q22 = gregoryBlend(getTangent(2, 0), getTangent(1, 1), um, vm);

    auto bpatch = get_bezier_boundary();

    bpatch(1, 1) = Q11;
    bpatch(2, 1) = Q21;
    bpatch(1, 2) = Q12;
    bpatch(2, 2) = Q22;

    return bpatch.evaluate(u, v);
  }

  void print(std::ostream& os) const
  {
    os << "GregoryPatch<" << NDIMS << "D>(";
    for(int i = 0; i < NPTS; ++i)
    {
      os << m_controlPoints[i];
      it if(i + 1 < NPTS) { os << ", "; }
    }
    os << ")";
  }

private:
  // Copies over the boundary points to a BezierPatch object,
  //  leaving the 4 interior control points uninitialized
  primal::BezierPatch<double, 3> get_bezier_boundary() const
  {
    BezierPatch<T, NDIMS> bpatch(3, 3);

    bpatch(0, 0) = getCorner(0);
    bpatch(1, 0) = getBoundaryPoint(0, 1);
    bpatch(2, 0) = getBoundaryPoint(0, 2);
    bpatch(3, 0) = getCorner(1);

    bpatch(3, 1) = getBoundaryPoint(1, 1);
    bpatch(3, 2) = getBoundaryPoint(1, 2);
    bpatch(3, 3) = getCorner(2);

    bpatch(2, 3) = getBoundaryPoint(2, 1);
    bpatch(1, 3) = getBoundaryPoint(2, 2);
    bpatch(0, 3) = getCorner(3);

    bpatch(0, 2) = getBoundaryPoint(3, 1);
    bpatch(0, 1) = getBoundaryPoint(3, 2);

    return bpatch;
  }

  CoordsVec m_controlPoints;

  // Map of boundary curve control points into internal storage
  static constexpr int m_index_map[4][4] = {{/*V0*/ 0, /*E01*/ 4, /*E02*/ 5, /*V1*/ 1},
                                            {/*V1*/ 1, /*E11*/ 6, /*E12*/ 7, /*V2*/ 2},
                                            {/*V2*/ 2, /*E21*/ 8, /*E22*/ 9, /*V3*/ 3},
                                            {/*V3*/ 3, /*E31*/ 10, /*E32*/ 11, /*V0*/ 0}};
};

//------------------------------------------------------------------------------
/// Free functions related to GregoryPatch
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const GregoryPatch<T, NDIMS>& nPatch)
{
  nPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::GregoryPatch using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::GregoryPatch<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_GREGORY_PATCH_HPP_
