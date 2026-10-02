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
  //  - 12 exterior control points across three degree-elevated cubic curves
  //  - 6 interior control points, two for each boundary curve
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
  ///@{
  /**
   * @name Constructors for GregoryTriangle
   *
   * The constructors allow initialization from:
   * - the 18 Gregory triangle control points,
   * - a polynomial quartic Bezier triangle,
   * - C-style arrays, Axom StackArrays, or Axom ArrayViews,
   * - three corner positions with associated corner normal vectors.
   *
   * The 18-point control net is stored as:
   * - indices 0-2: corners,
   * - indices 3-11: three boundary control points for each edge,
   * - indices 12-17: two Gregory tangent points for each edge.
   *
   * Boundary edge \a e is directed from corner `(e+1)%3` to corner `(e+2)%3`.
   */

  /*!
   * \brief Default constructor for a Gregory triangle
   *
   * The fixed-size control net is default-initialized.
   */
  GregoryTriangle() = default;

  /*!
   * \brief Constructor from an ArrayView over the control points
   *
   * \param [in] controlPoints ArrayView of the 18 Gregory triangle control points
   * \pre \a controlPoints must contain exactly `NPTS` points
   */
  explicit GregoryTriangle(ArrayView<const PointType> controlPoints)
  {
    SLIC_ASSERT(controlPoints.size() == NPTS);
    SLIC_ASSERT(controlPoints.data() != nullptr);
    for(int i = 0; i < NPTS; ++i)
    {
      m_controlPoints[i] = controlPoints[i];
    }
  }

  /*!
   * \brief Constructor from a non-const ArrayView over the control points
   *
   * \param [in] controlPoints ArrayView of the 18 Gregory triangle control points
   * \pre \a controlPoints must contain exactly `NPTS` points
   */
  explicit GregoryTriangle(ArrayView<PointType> controlPoints)
    : GregoryTriangle(ArrayView<const PointType>(controlPoints.data(), controlPoints.size()))
  { }

  /*!
   * \brief Constructor from a C-style array of control points
   *
   * \param [in] pts A C-style array of 18 Gregory triangle control points
   * \pre \a pts must be non-null and contain at least `NPTS` points
   */
  explicit GregoryTriangle(const PointType* pts)
    : GregoryTriangle(ArrayView<const PointType>(pts, NPTS))
  { }

  /*!
   * \brief Constructor from a C-style array of control points
   *
   * \param [in] pts A C-style array of 18 Gregory triangle control points
   * \pre \a pts must be non-null and contain at least `NPTS` points
   */
  explicit GregoryTriangle(PointType* pts) : GregoryTriangle(ArrayView<const PointType>(pts, NPTS))
  { }

  /*!
   * \brief Constructor from an Axom StackArray of control points
   *
   * \param [in] pts StackArray containing the 18 Gregory triangle control points
   */
  explicit GregoryTriangle(const CoordsVec& pts)
    : GregoryTriangle(ArrayView<const PointType>(pts.data(), pts.size()))
  { }

  /*!
   * \brief Constructor from a polynomial quartic Bezier triangle
   *
   * \param [in] bTri A polynomial Bezier triangle of order 4
   *
   * This creates a Gregory triangle that exactly reproduces the input quartic Bezier triangle.
   * The Gregory tangent pairs are duplicated from the three Bezier interior control points,
   * causing the parameter-dependent Gregory blends to collapse to fixed Bezier points.
   *
   * \pre \a bTri must have order 4
   * \pre \a bTri must be polynomial, not rational
   */
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

  /*!
   * \brief Constructor from vertex points and corner normal vectors
   *
   * \param [in] nodePositions ArrayView of the three corner positions
   * \param [in] nodeVectors ArrayView of the three corner normal vectors
   *
   * Deterministically compute hybrid cubic-quartic boundary control points and 
   * Gregory tangent points using local vertex information.
   *
   * Algorithm derived from the rectangular analog in
   *  A 3D contact smoothing method using Gregory patches, 
   *  Michael Anthony Puso, Tod A. Laursen, International Journal for Numerical Methods in Engineering
   *  Volume 54, Issue 8 (June 2002)
   
   * \pre \a nodePositions and \a nodeVectors must each contain exactly 3 entries
   */
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
      cubic.degreeElevate(4);
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

  ///@}

  /*!
   * \brief Returns the \a i-th corner point, oriented ccw
   *
   * \param [in] i Corner index in `[0, 2]`
   */
  PointType& getCorner(int i)
  {
    SLIC_ASSERT(i >= 0 && i < 3);
    return m_controlPoints[i];
  }

  /*!
   * \brief Returns the \a i-th corner point, oriented ccw
   *
   * \param [in] i Corner index in `[0, 2]`
   */
  const PointType& getCorner(int i) const
  {
    SLIC_ASSERT(i >= 0 && i < 3);
    return m_controlPoints[i];
  }

  /*!
   * \brief Returns a Gregory tangent point for an edge
   *
   * \param [in] e Edge index in `[0, 2]`, oriented ccw
   * \param [in] t Tangent point index in `[0, 1]`
   */
  PointType& getTangent(int e, int t)
  {
    SLIC_ASSERT(e >= 0 && e < 3);
    SLIC_ASSERT(t >= 0 && t < 2);
    return m_controlPoints[12 + 2 * e + t];
  }

  /*!
   * \brief Returns a Gregory tangent point for an edge
   *
   * \param [in] e Edge index in `[0, 2]`, oriented ccw
   * \param [in] t Tangent point index in `[0, 1]`
   */
  const PointType& getTangent(int e, int t) const
  {
    SLIC_ASSERT(e >= 0 && e < 3);
    SLIC_ASSERT(t >= 0 && t < 2);
    return m_controlPoints[12 + 2 * e + t];
  }

  /*!
   * \brief Returns the two Gregory tangent points adjacent to a corner
   *
   * \param [in] i Corner index in `[0, 2]`, oriented ccw
   * \param [out] v0 Tangent point from the preceding edge
   * \param [out] v1 Tangent point from the following edge
   */
  void getTangentsByCorner(int i, PointType& v0, PointType& v1) const
  {
    SLIC_ASSERT(i >= 0 && i < 3);
    v0 = getTangent((i + 1) % 3, 1);
    v1 = getTangent((i + 2) % 3, 0);
  }

  /*!
   * \brief Returns a control point on a boundary edge
   *
   * \param [in] e Edge index in `[0, 2]`
   * \param [in] k Boundary point index in `[0, 4]`
   */
  PointType& getBoundaryPoint(int e, int k)
  {
    SLIC_ASSERT(e >= 0 && e < 3);
    SLIC_ASSERT(k >= 0 && k < 5);
    return m_controlPoints[s_edge_index_map[e][k]];
  }

  /*!
   * \brief Returns a control point on a boundary edge
   *
   * \param [in] e Edge index in `[0, 2]`
   * \param [in] k Boundary point index in `[0, 4]`
   */
  const PointType& getBoundaryPoint(int e, int k) const
  {
    SLIC_ASSERT(e >= 0 && e < 3);
    SLIC_ASSERT(k >= 0 && k < 5);
    return m_controlPoints[s_edge_index_map[e][k]];
  }

  /*!
   * \brief Returns a reference to the triangle's control points
   */
  CoordsVec& getControlPoints() { return m_controlPoints; }

  /// \brief Returns a reference to the triangle's control points
  const CoordsVec& getControlPoints() const { return m_controlPoints; }

  /*!
   * \brief Evaluates the Gregory triangle at the given parameter values
   *
   * \param [in] u0 Parameter value on the first axis
   * \param [in] v0 Parameter value on the second axis
   *
   * A Gregory triangle is evaluated by constructing the equivalent quartic Bezier triangle whose
   * interior control points are blended from the Gregory tangent points at (\a u0, \a v0).
   */
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

  /*!
   * \brief Evaluates the first derivative in the u direction
   *
   * \param [in] u Parameter value on the first axis
   * \param [in] v Parameter value on the second axis
   */
  VectorType du(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return Du;
  }

  /*!
   * \brief Evaluates the first derivative in the v direction
   *
   * \param [in] u Parameter value on the first axis
   * \param [in] v Parameter value on the second axis
   */
  VectorType dv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return Dv;
  }

  /*!
   * \brief Evaluates the second derivative in the u direction
   *
   * \param [in] u Parameter value on the first axis
   * \param [in] v Parameter value on the second axis
   */
  VectorType dudu(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDu;
  }

  /*!
   * \brief Evaluates the second derivative in the v direction
   *
   * \param [in] u Parameter value on the first axis
   * \param [in] v Parameter value on the second axis
   */
  VectorType dvdv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DvDv;
  }

  /*!
   * \brief Evaluates the mixed second derivative
   *
   * \param [in] u Parameter value on the first axis
   * \param [in] v Parameter value on the second axis
   */
  VectorType dudv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDv;
  }

  /// \brief Returns an axis-aligned bounding box containing the triangle
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  /// \brief Returns an oriented bounding box containing the triangle
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Simple formatted print of a Gregory Triangle instance
   *
   * \param os The output stream to write to
   */
  void print(std::ostream& os) const
  {
    os << "GregoryTriangle(vertices [";
    for(int i = 0; i < 3; ++i)
    {
      os << getCorner(i) << (i < 2 ? ", " : "]");
    }

    os << ", edge points [";
    for(int e = 0; e < 3; ++e)
    {
      for(int k = 1; k < 4; ++k)
      {
        os << getBoundaryPoint(e, k) << (e < 2 || k < 3 ? ", " : "]");
      }
    }

    os << ", tangent points [";
    for(int e = 0; e < 3; ++e)
    {
      for(int t = 0; t < 2; ++t)
      {
        os << getTangent(e, t) << (e < 2 || t < 1 ? ", " : "]");
      }
    }

    os << ")";
  }

private:
  /*!
   * \brief Stores the temporary Bezier triangle and blended interior point derivatives
   *
   * The Gregory triangle evaluation converts the control net to a quartic Bezier triangle at a
   * specific parameter value. The three interior Bezier points, `Q`, depend on the evaluation
   * parameters, so derivative evaluation also requires their first and second derivatives.
   */
  struct IntermediateBlendingDerivatives
  {
    /// \brief Equivalent quartic Bezier triangle for the requested parameter value
    BezierTriangle<T, 3> btri;

    /// \brief Blended interior Bezier control points and derivatives
    PointType Q[3];
    VectorType Q_u[3];
    VectorType Q_v[3];
    VectorType Q_uu[3];
    VectorType Q_vv[3];
    VectorType Q_uv[3];
  };

  /*!
   * \brief Constructs the equivalent Bezier triangle and parameter-dependent interior data
   *
   * \param [in] u0 Parameter value on the first axis
   * \param [in] v0 Parameter value on the second axis
   * \param [in] derivative_order Highest derivative order to compute, in `[0, 2]`
   *
   * The returned quartic Bezier triangle has the Gregory boundary control points copied directly
   * and the three interior control points blended from the Gregory tangent points.
   */
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
   * \brief Copies the Gregory boundary into a quartic BezierTriangle object
   *
   * The returned triangle has its 12 exterior control points initialized from the Gregory
   * triangle boundary. The three interior control points are intentionally left uninitialized.
   * 
   * \sa set_bezier_interior(bezierTriangle<T, 3>&, const PointType[3])
   */
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

  /*!
   * \brief Maps boundary curve control point indices to control net storage indices
   *
   * The first index selects a directed edge. The second index selects one of the five
   * degree-elevated cubic boundary control points on that edge.
   */
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
