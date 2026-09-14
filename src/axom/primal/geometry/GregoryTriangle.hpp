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
    SLIC_ASSERT(controlPoints.data() != nullptr);
    for(int i = 0; i < NPTS; ++i)
    {
      m_controlPoints[i] = controlPoints[i];
    }
  }

  explicit GregoryTriangle(ArrayView<PointType> controlPoints)
    : GregoryTriangle(ArrayView<const PointType>(controlPoints.data(), controlPoints.size()))
  { }

  explicit GregoryTriangle(const PointType* pts)
    : GregoryTriangle(ArrayView<const PointType>(pts, NPTS))
  { }

  explicit GregoryTriangle(PointType* pts) : GregoryTriangle(ArrayView<const PointType>(pts, NPTS))
  { }

  explicit GregoryTriangle(const CoordsVec& pts)
    : GregoryTriangle(ArrayView<const PointType>(pts.data(), pts.size()))
  { }

  explicit GregoryTriangle(const BezierTriangle<T, 3>& bTri)
  {
    SLIC_ASSERT(bTri.getOrder() == 4);
    SLIC_ASSERT(!bTri.isRational());

    getCorner(0) = bTri(0, 0);
    getCorner(1) = bTri(0, 4);
    getCorner(2) = bTri(4, 0);

    getBoundaryPoint(0, 1) = bTri(1, 3);
    getBoundaryPoint(0, 2) = bTri(2, 2);
    getBoundaryPoint(0, 3) = bTri(3, 1);
    getBoundaryPoint(1, 1) = bTri(3, 0);
    getBoundaryPoint(1, 2) = bTri(2, 0);
    getBoundaryPoint(1, 3) = bTri(1, 0);
    getBoundaryPoint(2, 1) = bTri(0, 1);
    getBoundaryPoint(2, 2) = bTri(0, 2);
    getBoundaryPoint(2, 3) = bTri(0, 3);

    getTangent(0, 0) = bTri(1, 2);
    getTangent(0, 1) = bTri(2, 1);
    getTangent(1, 0) = bTri(2, 1);
    getTangent(1, 1) = bTri(1, 1);
    getTangent(2, 0) = bTri(1, 1);
    getTangent(2, 1) = bTri(1, 2);
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
      const int start = (k + 1) % 3;
      const int end = (k + 2) % 3;
      const VectorType dx(nodePositions[start], nodePositions[end]);

      const VectorType c0 = (dx - dx.dot(v[start]) * v[start]) / 3.0;
      const VectorType c2 = (dx - dx.dot(v[end]) * v[end]) / 3.0;

      // Define the cubic Bezier which represents the boundary of the curve
      BezierType cubic(axom::Array {nodePositions[start],
                                    PointType {nodePositions[start].array() + c0.array()},
                                    PointType {nodePositions[end].array() - c2.array()},
                                    nodePositions[end]},
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

    // Define the interior control nodes at each vertex
    for(int k = 0; k < 3; ++k)
    {
      const int start = (k + 1) % 3;
      const int end = (k + 2) % 3;
      const int prev_edge = (k + 2) % 3;
      const int next_edge = (k + 1) % 3;

      // Get control points at and around the edge's starting vertex
      const auto& p0 = getBoundaryPoint(prev_edge, 3);  // Previous edge
      const auto& q0 = getBoundaryPoint(k, 0);          // Edge start
      const auto& q1 = getBoundaryPoint(k, 1);          // Current edge
      const auto deriv0 = 0.5 * VectorType(p0, q0) + 0.5 * VectorType(p0, q1);

      // Get control points at and around the edge's ending vertex
      const auto& q3 = getBoundaryPoint(k, 3);          // Current edge
      const auto& q4 = getBoundaryPoint(k, 4);          // Edge end
      const auto& p3 = getBoundaryPoint(next_edge, 1);  // Next edge
      const auto deriv1 = 0.5 * VectorType(p3, q3) + 0.5 * VectorType(p3, q4);

      // Compute a boundary cross derivative that varies across the edge
      const VectorType dx(getCorner(start), getCorner(end));
      const VectorType a0 = VectorType::cross_product(nodeVectors[start], dx).unitVector();
      const VectorType a3 = VectorType::cross_product(nodeVectors[end], dx).unitVector();

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
    v0 = getTangent((i + 1) % 3, 1);
    v1 = getTangent((i + 2) % 3, 0);
  }

  PointType& getBoundaryPoint(int e, int k) { return m_controlPoints[s_edge_index_map[e][k]]; }
  const PointType& getBoundaryPoint(int e, int k) const
  {
    return m_controlPoints[s_edge_index_map[e][k]];
  }

  CoordsVec& getControlPoints() { return m_controlPoints; }
  const CoordsVec& getControlPoints() const { return m_controlPoints; }

  // Evaluate the triangle by constructing the equivalent Bezier triangle with interior control nodes
  //  defined in terms of the tangent vectors and the evaluation parameters
  PointType evaluate(T u0, T v0) const
  {
    const auto intermediate = setup_intermediate_bezier(u0, v0, 0);
    return intermediate.btri.evaluate(u0, v0);
  }

  /*!
   * \brief Evaluates all first derivatives of the Gregory triangle at (\a u0, \a v0)
   *
   * \param [in] u0 Parameter value at which to evaluate on the first axis
   * \param [in] v0 Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the Gregory triangle at (u0, v0)
   * \param [out] Du The vector value of S_u(u0, v0)
   * \param [out] Dv The vector value of S_v(u0, v0)
  */
  void evaluateFirstDerivatives(T u0, T v0, PointType& eval, VectorType& Du, VectorType& Dv) const
  {
    const auto intermediate = setup_intermediate_bezier(u0, v0, 1);
    intermediate.btri.evaluateFirstDerivatives(u0, v0, eval, Du, Dv);

    // Chain rule correction due to (u0,v0)-dependent interior control points.
    // This follows BezierTriangle's barycentric convention:
    //   {u, v, w} = {1 - u0 - v0, v0, u0}
    const T u = T(1) - u0 - v0;
    const T v = v0;
    const T w = u0;

    axom::StaticArray<T, 3> B, B_u0, B_v0;
    evaluate_quartic_interior_basis(u, v, w, B, B_u0, B_v0);

    Du += B[0] * intermediate.Q_u[0] + B[1] * intermediate.Q_u[1] + B[2] * intermediate.Q_u[2];
    Dv += B[0] * intermediate.Q_v[0] + B[1] * intermediate.Q_v[1] + B[2] * intermediate.Q_v[2];
  }

  /*!
   * \brief Evaluates all second derivatives of the Gregory triangle at (\a u0, \a v0)
   *
   * \param [in] u0 Parameter value at which to evaluate on the first axis
   * \param [in] v0 Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the Gregory triangle at (u0, v0)
   * \param [out] Du The vector value of S_u(u0, v0)
   * \param [out] Dv The vector value of S_v(u0, v0)
   * \param [out] DuDu The vector value of S_uu(u0, v0)
   * \param [out] DvDv The vector value of S_vv(u0, v0)
   * \param [out] DuDv The vector value of S_uv(u0, v0) == S_vu(u0, v0)
   */
  void evaluateSecondDerivatives(T u0,
                                 T v0,
                                 PointType& eval,
                                 VectorType& Du,
                                 VectorType& Dv,
                                 VectorType& DuDu,
                                 VectorType& DvDv,
                                 VectorType& DuDv) const
  {
    const auto intermediate = setup_intermediate_bezier(u0, v0, 2);
    intermediate.btri.evaluateSecondDerivatives(u0, v0, eval, Du, Dv, DuDu, DvDv, DuDv);

    const T u = T(1) - u0 - v0;
    const T v = v0;
    const T w = u0;

    axom::StaticArray<T, 3> B, B_u0, B_v0;
    evaluate_quartic_interior_basis(u, v, w, B, B_u0, B_v0);

    // First derivative corrections
    Du += B[0] * intermediate.Q_u[0] + B[1] * intermediate.Q_u[1] + B[2] * intermediate.Q_u[2];
    Dv += B[0] * intermediate.Q_v[0] + B[1] * intermediate.Q_v[1] + B[2] * intermediate.Q_v[2];

    // Second derivative corrections
    DuDu += T(2) *
        (B_u0[0] * intermediate.Q_u[0] + B_u0[1] * intermediate.Q_u[1] +
         B_u0[2] * intermediate.Q_u[2]) +
      (B[0] * intermediate.Q_uu[0] + B[1] * intermediate.Q_uu[1] + B[2] * intermediate.Q_uu[2]);

    DvDv += T(2) *
        (B_v0[0] * intermediate.Q_v[0] + B_v0[1] * intermediate.Q_v[1] +
         B_v0[2] * intermediate.Q_v[2]) +
      (B[0] * intermediate.Q_vv[0] + B[1] * intermediate.Q_vv[1] + B[2] * intermediate.Q_vv[2]);

    DuDv += (B_u0[0] * intermediate.Q_v[0] + B_u0[1] * intermediate.Q_v[1] +
             B_u0[2] * intermediate.Q_v[2]) +
      (B_v0[0] * intermediate.Q_u[0] + B_v0[1] * intermediate.Q_u[1] + B_v0[2] * intermediate.Q_u[2]) +
      (B[0] * intermediate.Q_uv[0] + B[1] * intermediate.Q_uv[1] + B[2] * intermediate.Q_uv[2]);
  }

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

  VectorType dudu(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDu;
  }

  VectorType dvdv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DvDv;
  }

  VectorType dudv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDv;
  }

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
    VectorType Q_uu[3];
    VectorType Q_vv[3];
    VectorType Q_uv[3];
  };

  IntermediateBlendingDerivatives setup_intermediate_bezier(T u0, T v0, int derivative_order) const
  {
    IntermediateBlendingDerivatives out;
    out.btri = get_bezier_boundary();

    // Parameter convention matches BezierTriangle::evaluate():
    // barycentric coordinates {u, v, w} are {1-u0-v0, v0, u0}
    const T u = T(1) - u0 - v0;
    const T v = v0;
    const T w = u0;

    // clang-format off
    auto blend = [&](const PointType& A, const PointType& B, // Internal Gregory points
                     T wa, T wb,                             // Barycentric coordinates for eval
                     T wa_u0, T wb_u0, T wa_v0, T wb_v0,     // Derivatives of Barycentric coords
                     PointType& Q,
                     VectorType& Q_u0, VectorType& Q_v0,
                     VectorType& Q_u0u0, VectorType& Q_v0v0, VectorType& Q_u0v0) {
      // clang-format on
      const T denom = wa + wb;
      if(axom::utilities::isNearlyEqual(denom, T(0)))
      {
        Q = A;
        Q_u0 = VectorType(T(0));
        Q_v0 = VectorType(T(0));
        Q_u0u0 = VectorType(T(0));
        Q_v0v0 = VectorType(T(0));
        Q_u0v0 = VectorType(T(0));
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

      if(derivative_order >= 2)
      {
        const T denom_u0 = wa_u0 + wb_u0;
        const T denom_v0 = wa_v0 + wb_v0;
        Q_u0u0 = (-T(2) * denom_u0 / denom) * Q_u0;
        Q_v0v0 = (-T(2) * denom_v0 / denom) * Q_v0;
        Q_u0v0 = (-(denom_u0 * Q_v0 + denom_v0 * Q_u0)) / denom;
      }
      else
      {
        Q_u0u0 = VectorType(T(0));
        Q_v0v0 = VectorType(T(0));
        Q_u0v0 = VectorType(T(0));
      }
    };

    // Get the three (u0, v0)-dependent interior control points for the equivalent
    //  quartic Bezier triangle. Each is blended from the two Gregory points adjacent
    //  to the corresponding vertex,

    // clang-format off
    blend(getTangent(1, 1), getTangent(2, 0),
          w, v,
          T(1), T(0), T(0), T(1),
          out.Q[0],
          out.Q_u[0], out.Q_v[0],
          out.Q_uu[0], out.Q_vv[0], out.Q_uv[0]);

    blend(getTangent(0, 1), getTangent(1, 0),
          v, u,
          T(0), T(-1), T(1), T(-1),
          out.Q[1],
          out.Q_u[1], out.Q_v[1],
          out.Q_uu[1], out.Q_vv[1], out.Q_uv[1]);

    blend(getTangent(2, 1), getTangent(0, 0),
          u, w,
          T(-1), T(1), T(-1), T(0),
          out.Q[2],
          out.Q_u[2], out.Q_v[2],
          out.Q_uu[2], out.Q_vv[2], out.Q_uv[2]);
    // clang-format on

    set_bezier_interior(out.btri, out.Q);
    return out;
  }

  /*!
   * \brief Assigns the three interior control points of a biquartic Bezier triangle
   *
   * \param [in,out] btri The biquartic Bezier triangle to update
   * \param [in] Q The 3 interior control points
   */
  static void set_bezier_interior(BezierTriangle<T, 3>& btri, const PointType Q[3])
  {
    btri(1, 1) = Q[0];
    btri(2, 1) = Q[1];
    btri(1, 2) = Q[2];
  }

  /*!
   * \brief Evaluates triangular Bernstein basis functions and their first derivatives
   *
   * \param [in] u First standard barycentric coordinate, equal to `1-u0-v0`
   * \param [in] v Second standard barycentric coordinate, equal to `v0`
   * \param [in] w Third standard barycentric coordinate, equal to `u0`
   * \param [out] B Basis values for the three interior control points
   * \param [out] B_u0 First derivatives of \a B with respect to `u0`
   * \param [out] B_v0 First derivatives of \a B with respect to `v0`
   */
  static void evaluate_quartic_interior_basis(T u,
                                              T v,
                                              T w,
                                              axom::StaticArray<T, 3>& B,
                                              axom::StaticArray<T, 3>& B_u0,
                                              axom::StaticArray<T, 3>& B_v0)
  {
    B[0] = T(12) * w * v * u * u;  // (i,j)=(1,1), k=2
    B[1] = T(12) * w * w * v * u;  // (i,j)=(2,1), k=1
    B[2] = T(12) * w * v * v * u;  // (i,j)=(1,2), k=1

    B_u0[0] = T(12) * v * u * (u - T(2) * w);
    B_u0[1] = T(12) * w * v * (T(2) * u - w);
    B_u0[2] = T(12) * v * v * (u - w);

    B_v0[0] = T(12) * w * u * (u - T(2) * v);
    B_v0[1] = T(12) * w * w * (u - v);
    B_v0[2] = T(12) * w * v * (T(2) * u - v);
  }

  // Copies over the boundary points to a BezierTriangle object,
  //  leaving the 3 interior control points uninitialized
  BezierTriangle<T, 3> get_bezier_boundary() const
  {
    BezierTriangle<T, 3> btri(4);

    // Edge 0
    btri(0, 4) = getBoundaryPoint(0, 0);
    btri(1, 3) = getBoundaryPoint(0, 1);
    btri(2, 2) = getBoundaryPoint(0, 2);
    btri(3, 1) = getBoundaryPoint(0, 3);
    btri(4, 0) = getBoundaryPoint(0, 4);

    // Edge 1
    // btri(4, 0) = getBoundaryPoint(1, 0);
    btri(3, 0) = getBoundaryPoint(1, 1);
    btri(2, 0) = getBoundaryPoint(1, 2);
    btri(1, 0) = getBoundaryPoint(1, 3);
    btri(0, 0) = getBoundaryPoint(1, 4);

    // Edge 2
    // btri(0, 0) = getBoundaryPoint(2, 0);
    btri(0, 1) = getBoundaryPoint(2, 1);
    btri(0, 2) = getBoundaryPoint(2, 2);
    btri(0, 3) = getBoundaryPoint(2, 3);
    // btri(0, 4) = getBoundaryPoint(2, 4);

    return btri;
  }

  CoordsVec m_controlPoints;

  // Map of BezierTriangle-style boundary curve control points into internal storage
  static constexpr int s_edge_index_map[3][5] = {
    {/*V1*/ 1, /*E01*/ 6, /*E02*/ 7, /*E03*/ 8, /*V2*/ 2},
    {/*V2*/ 2, /*E11*/ 9, /*E12*/ 10, /*E13*/ 11, /*V0*/ 0},
    {/*V0*/ 0, /*E21*/ 3, /*E22*/ 4, /*E23*/ 5, /*V1*/ 1}};
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
