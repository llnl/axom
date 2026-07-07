// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file GregoryTriangle.hpp
 *
 * \brief A bicubic Gregory triangle primitive
 */

#ifndef AXOM_PRIMAL_GREGORY_TRIANGLE_HPP_
#define AXOM_PRIMAL_GREGORY_TRIANGLE_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BezierTriangle.hpp"
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
template <typename T>
class GregoryTriangle;

/*! \brief Overloaded output operator for Gregory Triangles*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const GregoryTriangle<T>& nTri);

/*!
 * \class GregoryTriangle
 *
 * \brief Represents a 3D Gregory triangle defined by the control points of 3 degree-elevated
 *          cubic Bezier curves (i.e. quartic curves with identical geometry to cubics), and
 *          an additional two "Gregory points" for each edge which determine the surface.
 * 
 * Degree elevation of the boundary curves is necessary to provide sufficient degrees of freedom
 *   for the blending of the internal nodes.
 * 
 * \tparam T the coordinate type, e.g., double, float, etc.
 */
template <typename T>
class GregoryTriangle
{
public:
  // The number of control points for a hybrid quartic-cubic Gregory triangle is fixed:
  //  - 12 exterior control points for each of three degree-elevated cubic curves
  //  - 6 interior control points for each of 3 boundary curves
  static constexpr int NPTS = 18;

  using PointType = Point<T, 3>;
  using VectorType = Vector<T, 3>;

  using CoordsVec = axom::StackArray<PointType, NPTS>;
  using BezierType = BezierCurve<T, 3>;

  using BoundingBoxType = BoundingBox<T, 3>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, 3>;

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A Gregory Triangle must be defined using an arithmetic type");

public:
  GregoryTriangle() = default;

  explicit GregoryTriangle(ArrayView<const PointType> controlPoints)
  {
    SLIC_ASSERT(controlPoints.size() == NPTS);
    for(int i = 0; i < NPTS; ++i)
    {
      m_controlPoints[i] = controlPoints[i];
    }
  }

  GregoryTriangle(ArrayView<const PointType> nodePositions, ArrayView<const VectorType> nodeVectors)
  {
    // Store the position and orthogonal unit vector at each corner
    SLIC_ASSERT(nodePositions.size() == 3);
    SLIC_ASSERT(nodeVectors.size() == 3);

    axom::Array<VectorType> v(4);
    for(int i = 0; i < 3; ++i)
    {
      getCorner(i) = nodePositions[i];
      v[i] = nodeVectors[i].unitVector();
    }

    // Initialize the boundary points for the quartic Bezier triangle
    axom::StackArray<axom::StackArray<VectorType, 3>, 3> cubic_deriv_cp;
    for(int k = 0; k < 3; ++k)  // Loop over edges
    {
      const int kp1 = (k + 1) % 3;
      const VectorType dx(nodePositions[k], nodePositions[kp1]);

      const VectorType c0 = (dx - dx.dot(v[k]) * v[k]) / 3.0;
      const VectorType c2 = (dx - dx.dot(v[kp1]) * v[kp1]) / 3.0;

      // Define the cubic Bezier which represents the boundary of the curve
      BezierType cubic(axom::Array {nodePositions[k],
                                    PointType {nodePositions[k].array() + c0.array()},
                                    PointType {nodePositions[kp1].array() - c2.array()},
                                    nodePositions[kp1]},
                       3);

      // Store the control points of the derivative of this cubic for later
      cubic_deriv_cp[k][0] = VectorType(cubic[0], cubic[1]);
      cubic_deriv_cp[k][1] = VectorType(cubic[1], cubic[2]);
      cubic_deriv_cp[k][2] = VectorType(cubic[2], cubic[3]);

      // Do degree elevation on the cubic curve, which defines the quartic boundary control points
      cubic.degreeElevate();
      getBoundaryPoint(k, 0) = cubic[0];
      getBoundaryPoint(k, 1) = cubic[1];
      getBoundaryPoint(k, 2) = cubic[2];
      getBoundaryPoint(k, 3) = cubic[3];
      // getBoundaryPoint(i, 4) will be set for the next edge
    }

    // Define the interior control nodse at each vertex
    for(int k = 0; k < 3; ++k)
    {
      const int km1 = (k + 2) % 3;
      const int kp1 = (k + 1) % 3;

      // Get control points at and around the edge's starting vertex
      const auto& p0 = getBoundaryPoint(km1, 3);  // Vertex - 1
      const auto& q0 = getBoundaryPoint(k, 0);    // Vertex
      const auto& q1 = getBoundaryPoint(k, 1);    // Vertex + 1
      const auto deriv0 = 0.5 * VectorType(p0, q0) + 0.5 * VectorType(p0, q1);

      // Get control points at and around the edge's ending vertex
      const auto& q3 = getBoundaryPoint(k, 3);    // Vertex - 1
      const auto& q4 = getBoundaryPoint(k, 4);    // Vertex
      const auto& p3 = getBoundaryPoint(kp1, 1);  // Vertex + 1
      const auto deriv1 = 0.5 * VectorType(p3, q3) + 0.5 * VectorType(p3, q4);

      // Compute a boundary cross derivative that varies across the edge
      const VectorType dx(getCorner(k), getCorner(kp1));
      const VectorType a0 = VectorType::cross_product(nodeVectors[k], dx).unitVector();
      const VectorType a3 = VectorType::cross_product(nodeVectors[kp1], dx).unitVector();

      // Elevate the linear cross derivative a(t) = (1-t)a0 + t*a3 to quadratic
      const axom::StackArray<VectorType, 3> aHat = {a0, 0.5 * (a0 + a3), a3};

      // Compute the CPs for blending functions k(t) = (1-t)*k0 + t*k1 and
      //                                        h(t) = (1-t)*h- + t*h1
      auto& c = cubic_deriv_cp[k];
      const double k0 = aHat[0].dot(deriv0);
      const double k1 = aHat[2].dot(deriv1);

      const double h0 = c[0].dot(deriv0) / c[0].dot(c[0]);
      const double h1 = c[2].dot(deriv1) / c[2].dot(c[2]);

      // Compute cross derivatives at edge interior points and use them for interior CP
      axom::StackArray<VectorType, 2> deriv;
      for(int j = 1; j < 3; ++j)
      {
        const double fac = j / 3.0;
        deriv[j - 1] =
          (1. - fac) * (k0 * aHat[j] + h0 * c[j]) + fac * (k1 * aHat[j - 1] + h1 * c[j - 1]);
      }

      const auto& q2 = getBoundaryPoint(k, 2);
      getTangent(k, 0) = PointType(PointType::lerp(q1, q2, 0.5).array() - deriv[0].array());
      getTangent(k, 1) = PointType(PointType::lerp(q2, q3, 0.5).array() - deriv[1].array());
    }
  }

  PointType& getCorner(int i) { return m_controlPoints[i]; }
  const PointType& getCorner(int i) const { return m_controlPoints[i]; }

  PointType& getTangent(int e, int t) { return m_controlPoints[12 + 2 * e + t]; }
  const PointType& getTangent(int e, int t) const { return m_controlPoints[12 + 2 * e + t]; }

  void getTangentsByCorner(int i, PointType& v0, PointType& v1) const
  {
    v0 = getTangent((i + 2) % 3, 1);
    v1 = getTangent(i, 0);
  }

  PointType& getBoundaryPoint(int e, int k) { return m_controlPoints[s_edge_index_map[e][k]]; }
  const PointType& getBoundaryPoint(int e, int k) const
  {
    return m_controlPoints[s_edge_index_map[e][k]];
  }

  // Evaluate the triangle by constructing the equivalent Bezier triangle with interior control nodes
  //  defined in terms of the tangent vectors and the evaluation parameters
  PointType evaluate(T u0, T v0) const
  {
    const auto intermediate = setup_intermediate_bezier(u0, v0, 0);
    return intermediate.btri.evaluate(u0, v0);
  }

  /*!
   * \brief Evaluates all first derivatives of the Gregory patch at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the Gregory patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
  */
  void evaluateFirstDerivatives(T u0, T v0, PointType& eval, VectorType& Du, VectorType& Dv) const
  {
    const auto intermediate = setup_intermediate_bezier(u0, v0, 1);
    intermediate.btri.evaluateFirstDerivatives(u0, v0, eval, Du, Dv);

    // Chain rule correction due to (u0,v0)-dependent interior control points.
    // For a quartic Bezier triangle, the basis weight of control point (i,j)
    // is:  4!/(i! j! k!) * u^i v^j w^k  with k = 4-i-j and barycentric weights
    // (lambda0, lambda1, lambda2) = (1-u0-v0, v0, u0) corresponding to vertices
    // (0,0), (0,order), (order,0), respectively.
    const T u = u0;
    const T v = v0;
    const T w = T(1) - u0 - v0;

    const T w11 = T(12) * u * v * w * w;  // (i,j)=(1,1), k=2
    const T w21 = T(12) * u * u * v * w;  // (i,j)=(2,1), k=1
    const T w12 = T(12) * u * v * v * w;  // (i,j)=(1,2), k=1

    Du += w11 * intermediate.Q_u[0] + w21 * intermediate.Q_u[1] + w12 * intermediate.Q_u[2];
    Dv += w11 * intermediate.Q_v[0] + w21 * intermediate.Q_v[1] + w12 * intermediate.Q_v[2];
  }

  // Not yet completed
  // void evaluateSecondDerivatives(T u,
  //                                T v,
  //                                PointType& eval,
  //                                VectorType& Du,
  //                                VectorType& Dv,
  //                                VectorType& DuDu,
  //                                VectorType& DvDv,
  //                                VectorType& DuDv) const
  // {
  // }

  VectorType du(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return Du;
  }

  VectorType dv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return Dv;
  }

  // Not yet finished
  // VectorType dudu(T u, T v) const
  // {
  //   PointType eval;
  //   VectorType Du, Dv, DuDu, DvDv, DuDv;
  //   evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
  //   return DuDu;
  // }

  // VectorType dvdv(T u, T v) const
  // {
  //   PointType eval;
  //   VectorType Du, Dv, DuDu, DvDv, DuDv;
  //   evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
  //   return DvDv;
  // }

  // VectorType dudv(T u, T v) const
  // {
  //   PointType eval;
  //   VectorType Du, Dv, DuDu, DvDv, DuDv;
  //   evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
  //   return DuDv;
  // }

  /// \brief Returns an axis-aligned bounding box containing the patch
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  /// \brief Returns an oriented bounding box containing the patch
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  void print(std::ostream& os) const
  {
    os << "GregoryTriangle(";
    for(int i = 0; i < NPTS; ++i)
    {
      os << m_controlPoints[i];
      if(i + 1 < NPTS)
      {
        os << ", ";
      }
    }
    os << ")";
  }

private:
  struct IntermediateBlendingDerivatives
  {
    BezierTriangle<T, 3> btri;
    PointType Q[3];
    VectorType Q_u[3];
    VectorType Q_v[3];
  };

  IntermediateBlendingDerivatives setup_intermediate_bezier(T u0, T v0, int derivative_order) const
  {
    IntermediateBlendingDerivatives out;
    out.btri = get_bezier_boundary();

    // Parameter convention matches BezierTriangle::evaluate():
    // barycentric weights (lambda0, lambda1, lambda2) are (1-u0-v0, v0, u0)
    const T u = u0;
    const T v = v0;
    const T w = T(1) - u0 - v0;

    auto blend = [&](const PointType& A,
                     const PointType& B,
                     T wa,
                     T wb,
                     T wa_u0,
                     T wb_u0,
                     T wa_v0,
                     T wb_v0,
                     PointType& Q,
                     VectorType& Q_u0,
                     VectorType& Q_v0) {
      const T denom = wa + wb;
      if(axom::utilities::isNearlyEqual(denom, T(0)))
      {
        Q = A;
        Q_u0 = VectorType(T(0));
        Q_v0 = VectorType(T(0));
        return;
      }

      Q = PointType((wa * A.array() + wb * B.array()) / denom);

      if(derivative_order >= 1)
      {
        const auto dQ_u0 =
          (wa_u0 * (A.array() - Q.array()) + wb_u0 * (B.array() - Q.array())) / denom;
        const auto dQ_v0 =
          (wa_v0 * (A.array() - Q.array()) + wb_v0 * (B.array() - Q.array())) / denom;
        Q_u0 = VectorType(dQ_u0);
        Q_v0 = VectorType(dQ_v0);
      }
      else
      {
        Q_u0 = VectorType(T(0));
        Q_v0 = VectorType(T(0));
      }
    };

    // These are the three (u0,v0)-dependent interior control points for the
    // equivalent quartic Bezier triangle. Each is a blend of the two Gregory
    // points adjacent to the corresponding vertex.
    //
    // By convention, edge k connects vertex k -> vertex (k+1)%3.
    //  - getTangent(k,0) is the Gregory point near the starting vertex k
    //  - getTangent(k,1) is the Gregory point near the ending   vertex (k+1)%3
    //
    // Control net slots:
    //  Q[0] -> btri(1,1)
    //  Q[1] -> btri(2,1)
    //  Q[2] -> btri(1,2)
    blend(getTangent(2, 1),
          getTangent(0, 0),
          u,
          v,
          T(1),
          T(0),
          T(0),
          T(1),
          out.Q[0],
          out.Q_u[0],
          out.Q_v[0]);

    blend(getTangent(1, 1),
          getTangent(2, 0),
          v,
          w,
          T(0),
          T(-1),
          T(1),
          T(-1),
          out.Q[1],
          out.Q_u[1],
          out.Q_v[1]);

    blend(getTangent(0, 1),
          getTangent(1, 0),
          w,
          u,
          T(-1),
          T(1),
          T(-1),
          T(0),
          out.Q[2],
          out.Q_u[2],
          out.Q_v[2]);

    setBezierInterior(out.btri, out.Q);
    return out;
  }

  static void setBezierInterior(BezierTriangle<T, 3>& btri, const PointType Q[3])
  {
    btri(1, 1) = Q[0];
    btri(2, 1) = Q[1];
    btri(1, 2) = Q[2];
  }

  // Copies over the boundary points to a BezierPatch object,
  //  leaving the 4 interior control points uninitialized
  BezierTriangle<T, 3> get_bezier_boundary() const
  {
    BezierTriangle<T, 3> btri(4);

    // Edge 0
    btri(0, 0) = getBoundaryPoint(0, 0);
    btri(0, 1) = getBoundaryPoint(0, 1);
    btri(0, 2) = getBoundaryPoint(0, 2);
    btri(0, 3) = getBoundaryPoint(0, 3);
    btri(0, 4) = getBoundaryPoint(0, 4);

    // Edge 1
    // btri(0, 4) = getBoundaryPOint(1, 0);
    btri(1, 3) = getBoundaryPoint(1, 1);
    btri(2, 2) = getBoundaryPoint(1, 2);
    btri(3, 1) = getBoundaryPoint(1, 3);
    btri(4, 0) = getBoundaryPoint(1, 4);

    // Edge 2
    // btri(4, 0) = getBoundaryPoint(2, 0);
    btri(3, 0) = getBoundaryPoint(2, 1);
    btri(2, 0) = getBoundaryPoint(2, 2);
    btri(1, 0) = getBoundaryPoint(2, 3);
    // btri(0, 0) = getBoundaryPoint(2, 4);

    return btri;
  }

  CoordsVec m_controlPoints;

  // Map of boundary curve control points into internal storage
  static constexpr int s_edge_index_map[3][5] = {
    {/*V0*/ 0, /*E01*/ 3, /*E02*/ 4, /*E03*/ 5, /*V1*/ 1},
    {/*V0*/ 1, /*E01*/ 6, /*E02*/ 7, /*E03*/ 8, /*V1*/ 2},
    {/*V0*/ 2, /*E01*/ 9, /*E02*/ 10, /*E03*/ 11, /*V1*/ 0}};
};

//------------------------------------------------------------------------------
/// Free functions related to GregoryTriangle
//------------------------------------------------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const GregoryTriangle<T>& nPatch)
{
  nPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::GregoryTriangle using fmt
template <typename T>
struct axom::fmt::formatter<axom::primal::GregoryTriangle<T>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_GREGORY_TRIANGLE_HPP_
