// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file BezierTriangle.hpp
 *
 * \brief A BezierTriangle primitive
 */

#ifndef AXOM_PRIMAL_BEZIERTRIANGLE_HPP_
#define AXOM_PRIMAL_BEZIERTRIANGLE_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/squared_distance.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class BezierTriangle;

/*! \brief Overloaded output operator for Bezier Triangles*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierTriangle<T, NDIMS>& bTri);

/*!
 * \class BezierTriangle
 *
 * \brief Represents a Bezier triangle defined by a triangular array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * A Bezier triangle of order \a N has \f$ (N+1)(N+2)/2 \f$ control points.
 * It is parametrized over the domain \f$ u \ge 0, v \ge 0, u+v \le 1 \f$.
 *
 * Control points are indexed using integer coordinates \f$ (i,j) \f$ with
 * \f$ 0 \le i \le N \f$ and \f$ 0 \le j \le N-i \f$ and accessed via `operator()(i,j)`.
 * Internally, the triangular control net is stored in a 1D array and `triIndex(N,i,j)`
 * maps \f$ (i,j) \f$ to that linear storage index.
 *
 * Rational triangles are represented by an additional set of positive weights.
 * Polynomial (nonrational) Bezier triangles are identified by an empty weights array.
 *
 * \note This triangle uses permuted barycentric coordinates for evaluation such that, when
 * `getOrder()==1`, the parameter values correspond to the triangle vertices:
 * - `evaluate(0,0) == (*this)(0,0)`
 * - `evaluate(0,1) == (*this)(0,1)`
 * - `evaluate(1,0) == (*this)(1,0)`
 */
template <typename T, int NDIMS>
class BezierTriangle
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  using CoordsVec = axom::Array<PointType, 1>;
  using WeightsVec = axom::Array<T, 1>;

  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;
  using BezierCurveType = primal::BezierCurve<T, NDIMS>;

  AXOM_STATIC_ASSERT_MSG((NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
                         "A Bezier Triangle object may be defined in 1-, 2-, or 3-D");

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A Bezier Triangle must be defined using an arithmetic type");

public:
  ///@{
  /**
   * @name Constructors for BezierTriangle
   *
   * The constructors allow for flexible initialization of BezierTriangle objects from:
   * - 1D Axom arrays/views of control points and weights,
   * - C-style arrays of control points and weights,
   * - a specified polynomial order,
   * - rational or polynomial (nonrational) triangles, depending on the presence of weights.
   *
   * The triangle is parametrized over the domain \f$ u \ge 0, v \ge 0, u+v \le 1 \f$.
   *
   * Rational triangles are identified by a non-empty weights array,
   * and nonrational triangles by an empty weights array.
   * All weights must be greater than 0 in a rational triangle.
   *
   * For 1D control point/weight arrays, the expected layout corresponds to the indexing
   * used by `operator()(i,j)`:
     \verbatim
      pts[ triIndex(N,i,j) ]  <->  (*this)(i,j)   for 0<=i<=N and 0<=j<=N-i
     \endverbatim
   */

  /**
   * \brief Constructor from ArrayViews of control points and weights
   *
   * \param [in] controlPoints ArrayView of control points (size: (ord+1)(ord+2)/2 or 0)
   * \param [in] weights ArrayView of weights (size: (ord+1)(ord+2)/2 or 0)
   * \param [in] ord The triangle's polynomial order
   *
   * If \a controlPoints is empty, we still allocate space for the control points.
   * \pre ord must be greater than or equal to -1
   */
  BezierTriangle(axom::ArrayView<const PointType> controlPoints,
                 axom::ArrayView<const T> weights,
                 int ord)
    : m_ord(ord)
  {
    SLIC_ASSERT(ord >= -1);

    const int SZ = (m_ord >= 0) ? triSize(m_ord) : 0;
    SLIC_ASSERT(controlPoints.size() >= weights.size());

    // note: always allocate space for control points
    if(controlPoints.empty())
    {
      m_controlPoints.resize(SZ);
    }
    else
    {
      SLIC_ASSERT(controlPoints.data() != nullptr);
      SLIC_ASSERT(controlPoints.size() == SZ);
      m_controlPoints = controlPoints;
    }

    // note: only allocate space for weights when they are supplied
    if(!weights.empty())
    {
      SLIC_ASSERT(weights.data() != nullptr);
      SLIC_ASSERT(weights.size() == SZ);
      m_weights = weights;
      SLIC_ASSERT(is_valid_rational());
    }
  }

  /// Constructor from ArrayViews of (non-const) control points and weights
  BezierTriangle(axom::ArrayView<PointType> controlPoints, axom::ArrayView<T> weights, int ord)
    : BezierTriangle(axom::ArrayView<const PointType>(controlPoints.data(), controlPoints.size()),
                     axom::ArrayView<const T>(weights.data(), weights.size()),
                     ord)
  { }

  /*!
   * \brief Constructor for a polynomial (nonrational) Bezier Triangle that reserves space
   *
   * \param [in] ord The triangle's polynomial order
   * \pre ord must be greater than or equal to -1
   */
  explicit BezierTriangle(int ord = -1)
    : BezierTriangle(axom::ArrayView<const PointType>(nullptr, 0),
                     axom::ArrayView<const T>(nullptr, 0),
                     ord)
  { }

  /*!
   * \brief Constructor for a polynomial Bezier Triangle from an array of coordinates
   *
   * \param [in] pts A 1D C-style array of (ord+1)(ord+2)/2 control points
   * \param [in] ord The triangle's polynomial order
   * \pre ord is greater than or equal to zero
   */
  BezierTriangle(const PointType* pts, int ord)
    : BezierTriangle(axom::ArrayView<const PointType>(pts, triSize(ord)),
                     axom::ArrayView<const T>(nullptr, 0),
                     ord)
  { }

  /*!
   * \brief Constructor for a rational Bezier Triangle from arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of (ord+1)(ord+2)/2 control points
   * \param [in] weights A 1D C-style array of (ord+1)(ord+2)/2 positive weights
   * \param [in] ord The triangle's polynomial order
   * \pre ord is greater than or equal to zero
   *
   * If \a weights is the null pointer, creates a nonrational triangle.
   */
  BezierTriangle(const PointType* pts, const T* weights, int ord)
    : BezierTriangle(axom::ArrayView<const PointType>(pts, triSize(ord)),
                     axom::ArrayView<const T>(weights, weights ? triSize(ord) : 0),
                     ord)
  { }

  /*!
   * \brief Constructor from an Axom array of control points
   *
   * \param [in] pts A 1D Axom array of (ord+1)(ord+2)/2 control points
   * \param [in] ord The triangle's polynomial order (>= 0)
   */
  BezierTriangle(const CoordsVec& pts, int ord)
    : BezierTriangle(pts.view(), axom::ArrayView<const T>(nullptr, 0), ord)
  { }

  /*!
   * \brief Constructor from Axom arrays of control points and weights
   *
   * \param [in] pts A 1D Axom array of (ord+1)(ord+2)/2 control points
   * \param [in] weights A 1D Axom array of (ord+1)(ord+2)/2 positive weights
   * \param [in] ord The triangle's polynomial order (>= 0)
   */
  BezierTriangle(const CoordsVec& pts, const WeightsVec& weights, int ord)
    : BezierTriangle(pts.view(), weights.view(), ord)
  { }

  ///@}

  /*!
   * \brief Returns true when this triangle is rational
   *
   * A rational triangle has a weight per control point; polynomial (nonrational) triangles
   * are identified by an empty weight array.
   */
  bool isRational() const { return !m_weights.empty(); }

  /*!
   * \brief Returns a reference to the triangle's control points
   *
   * The control net contains `triSize(getOrder())` points (or 0 when `getOrder()<0`).
   */
  CoordsVec& getControlPoints() { return m_controlPoints; }

  /// \overload
  const CoordsVec& getControlPoints() const { return m_controlPoints; }

  /*!
   * \brief Returns a reference to the triangle's weights
   *
   * The weight array is empty for polynomial triangles. For rational triangles it contains
   * `triSize(getOrder())` positive weights.
   */
  WeightsVec& getWeights() { return m_weights; }

  /// \overload
  const WeightsVec& getWeights() const { return m_weights; }

  /*!
   * \brief Returns an axis-aligned bounding box containing the Bezier triangle
   *
   * \note The returned box is computed from the control points.
   */
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Returns an oriented bounding box containing the Bezier triangle
   *
   * \note The returned box is computed from the control points.
   */
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Sets the order of the Bezier triangle and resizes internal storage
   *
   * \param [in] ord The polynomial order 
   * 
   * \pre ord must be greater than or equal to -1
   *
   * \note This function only resizes the control point and weight arrays and does not
   * initialize their values. If the triangle is rational (`isRational()==true`), the
   * weight array is also resized to match the number of control points.
   */
  void setOrder(int ord)
  {
    SLIC_ASSERT(ord >= -1);

    m_ord = ord;

    const int SZ = (m_ord >= 0) ? triSize(m_ord) : 0;
    m_controlPoints.resize(SZ);
    if(isRational())
    {
      m_weights.resize(SZ);
    }
  }

  /*!
   * \brief Returns the polynomial order of the triangle
   *
   * \note A default constructed triangle has order -1 and no control points.
   */
  int getOrder() const { return m_ord; }

  /*!
   * \brief Access a control point in the triangular control net
   *
   * \param [in] i The first index (0 <= i <= getOrder())
   * \param [in] j The second index (0 <= j <= getOrder()-i)
   *
   * \pre \a i and \a j are in range and \a i+\a j <= getOrder()
   */
  PointType& operator()(int i, int j)
  {
    SLIC_ASSERT(i >= 0);
    SLIC_ASSERT(j >= 0);
    SLIC_ASSERT(i + j <= m_ord);
    return m_controlPoints[triIndex(m_ord, i, j)];
  }

  const PointType& operator()(int i, int j) const
  {
    SLIC_ASSERT(i >= 0);
    SLIC_ASSERT(j >= 0);
    SLIC_ASSERT(i + j <= m_ord);
    return m_controlPoints[triIndex(m_ord, i, j)];
  }

  /*!
   * \brief Evaluates the Bezier triangle at \a (u,v)
   *
   * \param [in] u Parameter value along the \a u axis
   * \param [in] v Parameter value along the \a v axis
   *
   * \pre getOrder() >= 0
   * \pre u >= 0, v >= 0, and u+v <= 1
   *
   * \return Point value S(u,v)
   *
   * \note In the rational case, evaluation is performed in projective space and divided
   * by the evaluated weight.
   */
  PointType evaluate(T u, T v) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(u >= T(0));
    SLIC_ASSERT(v >= T(0));
    SLIC_ASSERT(u + v <= T(1));

    if(!isRational())
    {
      PointType ptval;

      const int npts = m_controlPoints.size();
      axom::Array<T> dCarray(npts);

      // Run de Casteljau algorithm on each dimension
      for(int N = 0; N < NDIMS; ++N)
      {
        for(int n = 0; n < npts; ++n)
        {
          dCarray[n] = m_controlPoints[n][N];
        }

        for(int p = 1; p <= m_ord; ++p)
        {
          const int end = m_ord - p + 1;
          for(int i = 0; i < end; ++i)
          {
            for(int j = 0; j < end - i; ++j)
            {
              const auto& A = dCarray[triIndex(end, i, j)];
              const auto& B = dCarray[triIndex(end, i, j + 1)];
              const auto& C = dCarray[triIndex(end, i + 1, j)];
              dCarray[triIndex(end - 1, i, j)] = A + u * (C - A) + v * (B - A);
            }
          }
        }

        ptval[N] = dCarray[0];
      }

      return ptval;
    }
    else
    {
      // Rational case: evaluate in projective space, then divide by weight
      BezierTriangle<T, NDIMS> projective(m_ord);
      BezierTriangle<T, 1> weights(m_ord);
      fill_projective_triangles(projective, weights);

      const Point<T, NDIMS> P = projective.evaluate(u, v);
      const Point<T, 1> W = weights.evaluate(u, v);

      PointType eval;
      for(int N = 0; N < NDIMS; ++N)
      {
        eval[N] = P[N] / W[0];
      }
      return eval;
    }
  }

  /*!
   * \brief Evaluates first derivatives of the Bezier triangle at \a (u,v)
   *
   * \param [in] u Parameter value along the \a u axis
   * \param [in] v Parameter value along the \a v axis
   * \param [out] eval Point value S(u,v)
   * \param [out] Du First derivative S_u(u,v)
   * \param [out] Dv First derivative S_v(u,v)
   *
   * \pre getOrder() >= 0
   * \pre u >= 0, v >= 0, and u+v <= 1
   */
  void evaluateFirstDerivatives(T u,
                                T v,
                                Point<T, NDIMS>& eval,
                                Vector<T, NDIMS>& Du,
                                Vector<T, NDIMS>& Dv) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(u >= T(0));
    SLIC_ASSERT(v >= T(0));
    SLIC_ASSERT(u + v <= T(1));

    if(m_ord == 0)
    {
      eval = m_controlPoints[0];
      for(int N = 0; N < NDIMS; ++N)
      {
        Du[N] = T(0);
        Dv[N] = T(0);
      }
      return;
    }

    if(!isRational())
    {
      const int npts = m_controlPoints.size();
      axom::Array<T> dCarray(npts);

      // Run de Casteljau algorithm on each dimension
      for(int N = 0; N < NDIMS; ++N)
      {
        for(int n = 0; n < npts; ++n)
        {
          dCarray[n] = m_controlPoints[n][N];
        }

        for(int p = 1; p <= m_ord - 1; ++p)
        {
          const int end = m_ord - p + 1;
          for(int i = 0; i < end; ++i)
          {
            for(int j = 0; j < end - i; ++j)
            {
              const auto& A = dCarray[triIndex(end, i, j)];
              const auto& B = dCarray[triIndex(end, i, j + 1)];
              const auto& C = dCarray[triIndex(end, i + 1, j)];
              dCarray[triIndex(end - 1, i, j)] = A + u * (C - A) + v * (B - A);
            }
          }
        }

        // The last reduction yields a linear triangle:
        //   S(u,v) = A + u(C-A) + v(B-A)
        Du[N] = (dCarray[2] - dCarray[0]);
        Dv[N] = (dCarray[1] - dCarray[0]);
        eval[N] = dCarray[0] + u * Du[N] + v * Dv[N];

        Du[N] *= m_ord;
        Dv[N] *= m_ord;
      }
    }
    else
    {
      BezierTriangle<T, NDIMS> projective(m_ord);
      BezierTriangle<T, 1> weights(m_ord);
      fill_projective_triangles(projective, weights);

      Point<T, NDIMS> P;
      Vector<T, NDIMS> P_u, P_v;

      Point<T, 1> W;
      Vector<T, 1> W_u, W_v;

      projective.evaluateFirstDerivatives(u, v, P, P_u, P_v);
      weights.evaluateFirstDerivatives(u, v, W, W_u, W_v);

      for(int N = 0; N < NDIMS; ++N)
      {
        eval[N] = P[N] / W[0];
        Du[N] = (P_u[N] - eval[N] * W_u[0]) / W[0];
        Dv[N] = (P_v[N] - eval[N] * W_v[0]) / W[0];
      }
    }
  }

  /*!
   * \brief Evaluates all linear derivatives of a Bezier triangle at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate along the u axis
   * \param [in] v Parameter value at which to evaluate along the v axis
   * \param [out] eval The point value of the Bezier triangle at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   */
  void evaluateLinearDerivatives(T u,
                                 T v,
                                 Point<T, NDIMS>& eval,
                                 Vector<T, NDIMS>& Du,
                                 Vector<T, NDIMS>& Dv,
                                 Vector<T, NDIMS>& DuDv) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(u >= T(0));
    SLIC_ASSERT(v >= T(0));
    SLIC_ASSERT(u + v <= T(1));

    if(!isRational())
    {
      if(m_ord < 2)
      {
        evaluateFirstDerivatives(u, v, eval, Du, Dv);
        for(int N = 0; N < NDIMS; ++N)
        {
          DuDv[N] = T(0);
        }
        return;
      }

      const int npts = m_controlPoints.size();
      axom::Array<T> dCarray(npts);

      const T n_ord = static_cast<T>(m_ord);
      const T n_ord_nm1 = static_cast<T>(m_ord) * static_cast<T>(m_ord - 1);

      for(int N = 0; N < NDIMS; ++N)
      {
        for(int n = 0; n < npts; ++n)
        {
          dCarray[n] = m_controlPoints[n][N];
        }

        // Reduce to a quadratic triangle (order 2)
        for(int p = 1; p <= m_ord - 2; ++p)
        {
          const int end = m_ord - p + 1;
          for(int i = 0; i < end; ++i)
          {
            for(int j = 0; j < end - i; ++j)
            {
              const auto& A = dCarray[triIndex(end, i, j)];
              const auto& B = dCarray[triIndex(end, i, j + 1)];
              const auto& C = dCarray[triIndex(end, i + 1, j)];
              dCarray[triIndex(end - 1, i, j)] = A + u * (C - A) + v * (B - A);
            }
          }
        }

        // Extract the reduced quadratic triangle values
        const T Q00 = dCarray[triIndex(2, 0, 0)];
        const T Q01 = dCarray[triIndex(2, 0, 1)];
        const T Q02 = dCarray[triIndex(2, 0, 2)];
        const T Q10 = dCarray[triIndex(2, 1, 0)];
        const T Q11 = dCarray[triIndex(2, 1, 1)];
        const T Q20 = dCarray[triIndex(2, 2, 0)];

        // One reduction yields a linear triangle (order 1)
        const T L00 = Q00 + u * (Q10 - Q00) + v * (Q01 - Q00);
        const T L01 = Q01 + u * (Q11 - Q01) + v * (Q02 - Q01);
        const T L10 = Q10 + u * (Q20 - Q10) + v * (Q11 - Q10);

        eval[N] = L00 + u * (L10 - L00) + v * (L01 - L00);
        Du[N] = n_ord * (L10 - L00);
        Dv[N] = n_ord * (L01 - L00);
        DuDv[N] = n_ord_nm1 * (Q11 - Q10 - Q01 + Q00);
      }
    }
    else
    {
      BezierTriangle<T, NDIMS> projective(m_ord);
      BezierTriangle<T, 1> weights(m_ord);
      fill_projective_triangles(projective, weights);

      Point<T, NDIMS> P;
      Vector<T, NDIMS> P_u, P_v, P_uv;

      Point<T, 1> W;
      Vector<T, 1> W_u, W_v, W_uv;

      projective.evaluateLinearDerivatives(u, v, P, P_u, P_v, P_uv);
      weights.evaluateLinearDerivatives(u, v, W, W_u, W_v, W_uv);

      for(int N = 0; N < NDIMS; ++N)
      {
        eval[N] = P[N] / W[0];
        Du[N] = (P_u[N] - eval[N] * W_u[0]) / W[0];
        Dv[N] = (P_v[N] - eval[N] * W_v[0]) / W[0];
        DuDv[N] = (P_uv[N] - Du[N] * W_v[0] - Dv[N] * W_u[0] - eval[N] * W_uv[0]) / W[0];
      }
    }
  }

  /*!
   * \brief Evaluates all second derivatives of a Bezier triangle at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate along the u axis
   * \param [in] v Parameter value at which to evaluate along the v axis
   * \param [out] eval The point value of the Bezier triangle at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDu The vector value of S_uu(u, v)
   * \param [out] DvDv The vector value of S_vv(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   */
  void evaluateSecondDerivatives(T u,
                                 T v,
                                 Point<T, NDIMS>& eval,
                                 Vector<T, NDIMS>& Du,
                                 Vector<T, NDIMS>& Dv,
                                 Vector<T, NDIMS>& DuDu,
                                 Vector<T, NDIMS>& DvDv,
                                 Vector<T, NDIMS>& DuDv) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(u >= T(0));
    SLIC_ASSERT(v >= T(0));
    SLIC_ASSERT(u + v <= T(1));

    if(m_ord == 0)
    {
      eval = m_controlPoints[0];
      for(int N = 0; N < NDIMS; ++N)
      {
        Du[N] = T(0);
        Dv[N] = T(0);
        DuDu[N] = T(0);
        DvDv[N] = T(0);
        DuDv[N] = T(0);
      }
      return;
    }

    if(m_ord == 1)
    {
      evaluateFirstDerivatives(u, v, eval, Du, Dv);
      for(int N = 0; N < NDIMS; ++N)
      {
        DuDu[N] = T(0);
        DvDv[N] = T(0);
        DuDv[N] = T(0);
      }
      return;
    }

    if(!isRational())
    {
      const int npts = m_controlPoints.size();
      axom::Array<T> dCarray(npts);

      const T n_ord = static_cast<T>(m_ord);
      const T n_ord_nm1 = static_cast<T>(m_ord) * static_cast<T>(m_ord - 1);

      for(int N = 0; N < NDIMS; ++N)
      {
        for(int n = 0; n < npts; ++n)
        {
          dCarray[n] = m_controlPoints[n][N];
        }

        // Reduce to a quadratic triangle (order 2)
        for(int p = 1; p <= m_ord - 2; ++p)
        {
          const int end = m_ord - p + 1;
          for(int i = 0; i < end; ++i)
          {
            for(int j = 0; j < end - i; ++j)
            {
              const auto& A = dCarray[triIndex(end, i, j)];
              const auto& B = dCarray[triIndex(end, i, j + 1)];
              const auto& C = dCarray[triIndex(end, i + 1, j)];
              dCarray[triIndex(end - 1, i, j)] = A + u * (C - A) + v * (B - A);
            }
          }
        }

        // Extract the reduced quadratic triangle values
        const T Q00 = dCarray[triIndex(2, 0, 0)];
        const T Q01 = dCarray[triIndex(2, 0, 1)];
        const T Q02 = dCarray[triIndex(2, 0, 2)];
        const T Q10 = dCarray[triIndex(2, 1, 0)];
        const T Q11 = dCarray[triIndex(2, 1, 1)];
        const T Q20 = dCarray[triIndex(2, 2, 0)];

        // One reduction yields a linear triangle (order 1)
        const T L00 = Q00 + u * (Q10 - Q00) + v * (Q01 - Q00);
        const T L01 = Q01 + u * (Q11 - Q01) + v * (Q02 - Q01);
        const T L10 = Q10 + u * (Q20 - Q10) + v * (Q11 - Q10);

        eval[N] = L00 + u * (L10 - L00) + v * (L01 - L00);
        Du[N] = n_ord * (L10 - L00);
        Dv[N] = n_ord * (L01 - L00);

        // Second derivatives from second differences of the reduced quadratic triangle
        DuDu[N] = n_ord_nm1 * (Q20 - T(2) * Q10 + Q00);
        DvDv[N] = n_ord_nm1 * (Q02 - T(2) * Q01 + Q00);
        DuDv[N] = n_ord_nm1 * (Q11 - Q10 - Q01 + Q00);
      }
    }
    else
    {
      BezierTriangle<T, NDIMS> projective(m_ord);
      BezierTriangle<T, 1> weights(m_ord);
      fill_projective_triangles(projective, weights);

      Point<T, NDIMS> P;
      Vector<T, NDIMS> P_u, P_v, P_uu, P_vv, P_uv;

      Point<T, 1> W;
      Vector<T, 1> W_u, W_v, W_uu, W_vv, W_uv;

      projective.evaluateSecondDerivatives(u, v, P, P_u, P_v, P_uu, P_vv, P_uv);
      weights.evaluateSecondDerivatives(u, v, W, W_u, W_v, W_uu, W_vv, W_uv);

      for(int N = 0; N < NDIMS; ++N)
      {
        eval[N] = P[N] / W[0];
        Du[N] = (P_u[N] - eval[N] * W_u[0]) / W[0];
        Dv[N] = (P_v[N] - eval[N] * W_v[0]) / W[0];
        DuDu[N] = (P_uu[N] - T(2) * W_u[0] * Du[N] - eval[N] * W_uu[0]) / W[0];
        DvDv[N] = (P_vv[N] - T(2) * W_v[0] * Dv[N] - eval[N] * W_vv[0]) / W[0];
        DuDv[N] = (P_uv[N] - Du[N] * W_v[0] - Dv[N] * W_u[0] - eval[N] * W_uv[0]) / W[0];
      }
    }
  }

  /*!
   * \brief Computes a tangent of a Bezier triangle at (\a u, \a v) along the u axis
   */
  VectorType du(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return Du;
  }

  /*!
   * \brief Computes a tangent of a Bezier triangle at (\a u, \a v) along the v axis
   */
  VectorType dv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return Dv;
  }

  /*!
   * \brief Computes the second derivative of a Bezier triangle at (\a u, \a v) along the u axis
   */
  VectorType dudu(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDu;
  }

  /*!
   * \brief Computes the second derivative of a Bezier triangle at (\a u, \a v) along the v axis
   */
  VectorType dvdv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DvDv;
  }

  /*!
   * \brief Computes the mixed second derivative of a Bezier triangle at (\a u, \a v)
   */
  VectorType dudv(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDv;
  }

  /// Convenience alias for S_vu(u,v), which equals S_uv(u,v) for polynomial triangles
  VectorType dvdu(T u, T v) const { return dudv(u, v); }

  /*!
   * \brief Computes the normal vector of a Bezier triangle at (\a u, \a v)
   *
   * \note Only meaningful for NDIMS==3.
   */
  VectorType normal(T u, T v) const
  {
    Point<T, NDIMS> eval;
    Vector<T, NDIMS> Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return VectorType::cross_product(Du, Dv);
  }

  /*!
   * \brief Simple formatted print of a Bezier Triangle instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{ order " << m_ord << " Bezier Triangle ";

    for(int i = 0; i <= m_ord; ++i)
    {
      for(int j = 0; j <= m_ord - i; ++j)
      {
        os << (*this)(i, j) << ((i < m_ord || j < m_ord - i) ? "," : "");
      }
    }

    if(isRational())
    {
      os << ", weights [";
      for(int i = 0; i <= m_ord; ++i)
      {
        for(int j = 0; j <= m_ord - i; ++j)
        {
          os << m_weights[triIndex(m_ord, i, j)] << ((i < m_ord || j < m_ord - i) ? "," : "");
        }
      }
      os << "]";
    }

    os << "}";
    return os;
  }

  /*!
   * \brief Returns the number of control points for a triangle of order \a ord
   *
   * \param [in] ord Triangle order
   */
  static constexpr size_t triSize(int ord)
  {
    return (ord >= 0) ? static_cast<size_t>((ord + 1) * (ord + 2) / 2) : size_t {0};
  }

  /*!
   * \brief Maps triangular indices \a (i,j) to the linear storage index
   *
   * \param [in] ord Triangle order
   * \param [in] i First control net index
   * \param [in] j Second control net index
   *
   * \pre ord >= 0, i >= 0, j >= 0, and i+j <= ord
   */
  static constexpr size_t triIndex(int ord, int i, int j)
  {
    return static_cast<size_t>(i * (2 * ord + 3 - i) / 2 + j);
  }

private:
  /// Check that the weights used are positive, and
  ///  that there is one for each control node
  bool is_valid_rational() const
  {
    if(!isRational())
    {
      return true;
    }

    if(m_weights.size() != m_controlPoints.size())
    {
      return false;
    }

    for(int p = 0; p < m_weights.size(); ++p)
    {
      if(m_weights[p] <= 0)
      {
        return false;
      }
    }

    return true;
  }

  void fill_projective_triangles(BezierTriangle<T, NDIMS>& projective,
                                 BezierTriangle<T, 1>& weights) const
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(is_valid_rational());

    for(int i = 0; i <= m_ord; ++i)
    {
      for(int j = 0; j <= m_ord - i; ++j)
      {
        const T w = m_weights[triIndex(m_ord, i, j)];
        weights(i, j)[0] = w;

        for(int N = 0; N < NDIMS; ++N)
        {
          projective(i, j)[N] = (*this)(i, j)[N] * w;
        }
      }
    }
  }

private:
  int m_ord;

  CoordsVec m_controlPoints;
  WeightsVec m_weights;
};

//------------------------------------------------------------------------------
/// Free functions related to BezierTriangle
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierTriangle<T, NDIMS>& bTri)
{
  bTri.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_BEZIERTRIANGLE_HPP_
