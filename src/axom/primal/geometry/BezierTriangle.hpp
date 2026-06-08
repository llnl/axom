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
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
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
 * \tparam NDIMS The dimension of each control point, e.g. 1 for rational weights
 *
 * A Bezier triangle of order \a N has \f$ (N+1)(N+2)/2 \f$ control points.
 * It is parametrized over the domain \f$ u0 \ge 0, v0 \ge 0, u0+v0 \le 1 \f$.
 *
 * Control points are indexed using integer coordinates \f$ (i,j) \f$ with
 * \f$ 0 \le i \le N \f$ and \f$ 0 \le j \le N-i \f$ and accessed via `operator()(i,j)`.
 * Internally, the triangular control net is stored in a 1D array and `triIndex(N,i,j)`
 * maps \f$ (i,j) \f$ to that linear storage index.
 *
 * Rational triangles are represented by an additional set of positive weights.
 * Polynomial (nonrational) Bezier triangles are identified by an empty weights array.
 *
 * A default-constructed triangle will have order -1, and is "invalid".
 * Arrays of nodes and weights will be empty, and most methods are invalid 
 * 
 * \note This triangle uses permuted barycentric coordinates (u0, v0) for evaluation such that, when
 * `getOrder()==1`, the parameter values correspond to the triangle vertices:
 * - `evaluate(0,0) == (*this)(0,0)`
 * - `evaluate(0,1) == (*this)(0,1)`
 * - `evaluate(1,0) == (*this)(1,0)`
 * 
 * These are mapped to standard Barycentric coordinates {u,v,w} through (u0, v0) = {1 - u0 - v0, v0, u0}:
 *
 *   Parametric (u0, v0):       Barycentric {u,v,w}:
 *        (1, 0)                      {0,0,1}
 *          /\                          /\
 *         /  \                        /  \
 *   ^    /    \         <--->        /    \
 *   |   /      \                    /      \
 *   |  /        \                  /        \
 *  v0 /__________\                /__________\
 * (0, 0) u0 ---> (0, 1)       {1,0,0}        {0,1,0}
 * 
 */
template <typename T, int NDIMS>
class BezierTriangle
{
public:
  using PointType = Point<T, NDIMS>;
  using Barycentric = Point<T, 3>;
  using VectorType = Vector<T, NDIMS>;

  using CoordsVec = axom::Array<PointType, 1>;
  using WeightsVec = axom::Array<T, 1>;

  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;
  using BezierCurveType = primal::BezierCurve<T, NDIMS>;
  using EdgesVec = axom::Array<BezierCurveType>;

  AXOM_STATIC_ASSERT_MSG((NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
                         "A Bezier Triangle object may be defined in 1-, 2-, or 3-D");

  // Allows template types to access private member data of other allowable templates
  friend class BezierTriangle<T, 1>;
  friend class BezierTriangle<T, 2>;
  friend class BezierTriangle<T, 3>;

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
   *    pts[ triIndex(N,i,j) ]  <->  (*this)(i,j)   for 0<=i<=N and 0<=j<=N-i
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

  /// \brief Constructor from ArrayViews of (non-const) control points and weights
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

  /// Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      m_weights.resize(triSize(m_ord));
      m_weights.fill(1.0);
    }
  }

  /// Make nonrational by shrinking array of weights
  void makeNonrational() { m_weights.clear(); }

  /*!
   * \brief Returns a reference to the triangle's control points
   *
   * The control net contains `triSize(getOrder())` points (or 0 when `getOrder()<0`).
   */
  CoordsVec& getControlPoints() { return m_controlPoints; }

  /// \brief Returns a reference to the triangle's control points
  const CoordsVec& getControlPoints() const { return m_controlPoints; }

  /*!
   * \brief Returns a reference to the triangle's weights
   *
   * The weight array is empty for polynomial triangles. For rational triangles it contains
   * `triSize(getOrder())` positive weights.
   */
  WeightsVec& getWeights() { return m_weights; }

  /// \brief Returns a reference to the triangle's weights
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
   * \brief Return one vertex from the Bezier triangle
   *
   * \param [in] vertIdx Index of the requested vertex
   *
   * The vertices are returned in counter-clockwise order with respect to
   * the first control point (*this)(0, 0) == evaluate(0, 0)
   *
   * \return The PointType object at the vertex
   */
  PointType getVertex(int vertIdx) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(vertIdx >= 0 && vertIdx < 3);

    switch(vertIdx)
    {
    case 0:
      return (*this)(0, 0);
    case 1:
      return (*this)(0, m_ord);
    default:
      return (*this)(m_ord, 0);
    }
  }

  /*!
   * \brief Returns a copy of one of the boundary edges of the Bezier triangle
   *
   * \param [in] edgeIdx Index of the requested edge in \a [0,2]
   *
   * The edges are returned in counter-clockwise order with respect to the
   * parameter domain (u,v), with edge 0 across from corner 0 (i.e. evaluate(0,0)):
   * - \a edgeIdx = 0: u+v = 1 from `evaluate(1,0)` to `evaluate(0,1)`
   * - \a edgeIdx = 1: u = 0 from `evaluate(0,1)` to `evaluate(0,0)`
   * - \a edgeIdx = 2: v = 0 from `evaluate(0,0)` to `evaluate(1,0)`
   *
   * \return A copy of the Bezier curve representing the requested edge.
   */
  BezierCurveType getEdge(int edgeIdx) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(edgeIdx >= 0 && edgeIdx < 3);

    axom::Array<PointType> pts(m_ord + 1);
    axom::Array<T> wts;
    if(isRational())
    {
      wts.resize(m_ord + 1);
    }

    switch(edgeIdx)
    {
    case 0:
      for(int k = 0; k <= m_ord; ++k)
      {
        const int idx = triIndex(m_ord, k, m_ord - k);
        pts[k] = m_controlPoints[idx];
        if(isRational())
        {
          wts[k] = m_weights[idx];
        }
      }
      break;
    case 1:
      for(int k = 0; k <= m_ord; ++k)
      {
        const int idx = triIndex(m_ord, m_ord - k, 0);
        pts[k] = m_controlPoints[idx];
        if(isRational())
        {
          wts[k] = m_weights[idx];
        }
      }
      break;
    case 2:
      for(int k = 0; k <= m_ord; ++k)
      {
        const int idx = triIndex(m_ord, 0, k);
        pts[k] = m_controlPoints[idx];
        if(isRational())
        {
          wts[k] = m_weights[idx];
        }
      }
      break;
    default:
      break;
    }

    return isRational() ? BezierCurveType(pts, wts, m_ord) : BezierCurveType(pts, m_ord);
  }

  /*!
   * \brief Returns all three boundary edges of the Bezier triangle
   *
   * \return An array of three Bezier curves, ordered the same as `getEdge(int)`.
   */
  EdgesVec getEdges() const
  {
    SLIC_ASSERT(m_ord >= 0);
    EdgesVec edges;
    edges.reserve(3);
    edges.push_back(getEdge(0));
    edges.push_back(getEdge(1));
    edges.push_back(getEdge(2));
    return edges;
  }

  /*!
   * \brief Restricts this Bezier triangle to a subtriangle of the parameter domain
   *
   * See overload returning a `BezierTriangle` for the vertex mapping convention.
   *
   * \param [in] Qa Barycentric coordinates of the first subtriangle vertex `(u,v,w)`
   * \param [in] Qb Barycentric coordinates of the second subtriangle vertex `(u,v,w)`
   * \param [in] Qc Barycentric coordinates of the third subtriangle vertex `(u,v,w)`
   * \param [out] out Output restricted Bezier triangle
   *
   * \pre getOrder() >= 0
   *
   * \note The barycentric inputs \a Qa, \a Qb, \a Qc are standard Barycentric coordiantes
   *       related to parameter convention through (u0 = Qc, v0 = Qb)
   */
  void restrictToSubtriangle(const Barycentric& Qa,
                             const Barycentric& Qb,
                             const Barycentric& Qc,
                             BezierTriangle& out) const
  {
    using TriangularArray = axom::Array<PointType>;

    SLIC_ASSERT(m_ord >= 0);
    if(isRational())
    {
      // Rational case: restrict in projective space, then convert back.
      BezierTriangle<T, NDIMS> projective(m_ord);
      BezierTriangle<T, 1> weights(m_ord);
      fill_projective_triangles(projective, weights);

      BezierTriangle<T, NDIMS> proj_out(m_ord);
      BezierTriangle<T, 1> w_out(m_ord);

      projective.restrictToSubtriangle(Qa, Qb, Qc, proj_out);
      weights.restrictToSubtriangle(Qa, Qb, Qc, w_out);

      set_from_projective_triangles(proj_out, w_out, out);
      return;
    }

    const int n = m_ord;

    // Given a degree `deg` control net in `prev`, perform one reduction step so that
    // `next` becomes the degree `deg-1` control net at barycentric point `Q`.
    auto reduce_once =
      [&](const TriangularArray& prev, TriangularArray& next, int deg, const Barycentric& Q) {
        const int newDeg = deg - 1;

        for(int ii = 0; ii <= newDeg; ++ii)
        {
          for(int jj = 0; jj <= newDeg - ii; ++jj)
          {
            const auto& A0 = prev[triIndex(deg, ii, jj)];
            const auto& B0 = prev[triIndex(deg, ii, jj + 1)];
            const auto& C0 = prev[triIndex(deg, ii + 1, jj)];

            next[triIndex(newDeg, ii, jj)] = triInterpolate(A0, B0, C0, Q);
          }
        }
      };

    // The restricted control net over (Qa,Qb,Qc) is given by blossom values:
    //   out(i,j) = b(Qc^i, Qb^j, Qa^(n-i-j))  for 0<=i<=n and 0<=j<=n-i
    //
    // At layer p, we maintain all intermediate nets with exactly p fixed arguments,
    //  i.e. b(Qc^i0, Qb^j0, Qc^k0) with i0 + j0 + k0 = p, and use them to compute nets
    //  for blossom values with one more fixed argument until p == n
    axom::Array<TriangularArray> prevNets(1);
    prevNets[0] = m_controlPoints;

    for(int p = 1; p <= n; ++p)
    {
      const int degPrev = n - (p - 1);
      const int degCurr = degPrev - 1;

      // Allocate triangular arrays for each layer-p net.
      axom::Array<TriangularArray> currNets(triSize(p));
      for(auto& net : currNets)
      {
        net.resize(triSize(degCurr));
      }

      // Iterate over count triples with i0 + j0 + k0 == p.
      for(int i0 = 0; i0 <= p; ++i0)
      {
        for(int j0 = 0; j0 <= p - i0; ++j0)
        {
          const int k0 = p - i0 - j0;
          const int idx = triIndex(p, i0, j0);

          // Increase Va, then Vb, then Vc.
          // The order is arbitrary, since the blossom is symmetric.
          if(k0 > 0)
          {
            const int predIdx = triIndex(p - 1, i0, j0);
            reduce_once(prevNets[predIdx], currNets[idx], degPrev, Qa);
          }
          else if(j0 > 0)
          {
            const int predIdx = triIndex(p - 1, i0, j0 - 1);
            reduce_once(prevNets[predIdx], currNets[idx], degPrev, Qb);
          }
          else
          {
            SLIC_ASSERT(i0 > 0);
            const int predIdx = triIndex(p - 1, i0 - 1, j0);
            reduce_once(prevNets[predIdx], currNets[idx], degPrev, Qc);
          }
        }
      }

      // Discard the previous nets, since layer p determines layer p+1
      prevNets.swap(currNets);
    }

    // Each net in the last layer (p==n) has degree 0 and contains a single control point.
    out.setOrder(n);
    out.makeNonrational();
    for(int i = 0; i <= n; ++i)
    {
      for(int j = 0; j <= n - i; ++j)
      {
        const int idx = triIndex(n, i, j);
        out.m_controlPoints[idx] = prevNets[idx][0];
      }
    }
  }

  /*!
   * \brief Splits a Bezier triangle into three subtriangles by connecting an
   * interior parameter point \a (u0,v0) to the triangle's three vertices
   *
   * \param [in] u0 Parameter value along the \a u axis for the split point
   * \param [in] v0 Parameter value along the \a v axis for the split point
   * \param [out] t0 Subtriangle over the parameter triangle with vertices
   *  `(0,1)`, `(1,0)`, and `(u0,v0)` (preserves edge 0)
   * \param [out] t1 Subtriangle over the parameter triangle with vertices
   *  `(1,0)`, `(0,0)`, and `(u0,v0)` (preserves edge 1)
   * \param [out] t2 Subtriangle over the parameter triangle with vertices
   *  `(0,0)`, `(0,1)`, and `(u0,v0)` (preserves edge 2)
   *
   * \pre \a u0 > 0, \a v0 > 0, and \a u0 + \a v0 < 1
   *
   *            A
   *           /|\
   *          / | \
   *         /t2|t1\
   *        /  /Q\  \
   *       / /  t0 \ \
   *      //_________\\
   *     B             C
   *
   * \return t0 = Tri(B,C,Q), t1 = Tri(C,A,Q), t2 = Tri(A,B,Q)
   */
  void split(T u0, T v0, BezierTriangle& t0, BezierTriangle& t1, BezierTriangle& t2) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(is_valid_interior_parameter(u0, v0));

    if(isRational())
    {
      // Rational case: split in projective space and for weights, then convert back.
      BezierTriangle<T, NDIMS> projective(m_ord);
      BezierTriangle<T, 1> weights(m_ord);
      fill_projective_triangles(projective, weights);

      BezierTriangle<T, NDIMS> p0, p1, p2;
      BezierTriangle<T, 1> w0, w1, w2;
      projective.split(u0, v0, p0, p1, p2);
      weights.split(u0, v0, w0, w1, w2);

      set_from_projective_triangles(p0, w0, t0);
      set_from_projective_triangles(p1, w1, t1);
      set_from_projective_triangles(p2, w2, t2);
      return;
    }

    // Q is the split point in barycentric coordinates over the reference parameter triangle:
    const Barycentric Q {T(1) - u0 - v0, v0, u0};

    const int n = m_ord;
    axom::Array<axom::Array<PointType>> net(n + 1);
    net[0] = m_controlPoints;

    for(int p = 1; p <= n; ++p)
    {
      const int deg = n - p + 1;
      const int newDeg = deg - 1;

      net[p].resize(triSize(newDeg));
      const auto& prev = net[p - 1];

      for(int ii = 0; ii <= newDeg; ++ii)
      {
        for(int jj = 0; jj <= newDeg - ii; ++jj)
        {
          const auto& A0 = prev[triIndex(deg, ii, jj)];
          const auto& B0 = prev[triIndex(deg, ii, jj + 1)];
          const auto& C0 = prev[triIndex(deg, ii + 1, jj)];

          net[p][triIndex(newDeg, ii, jj)] = triInterpolate(A0, B0, C0, Q);
        }
      }
    }

    // Subtriangle (B,C,Q): (0,1), (1,0), (u0,v0)
    t0.setOrder(n);
    t0.makeNonrational();
    for(int i = 0; i <= n; ++i)
    {
      const int deg = n - i;
      for(int j = 0; j <= n - i; ++j)
      {
        t0(i, j) = net[i][triIndex(deg, j, deg - j)];
      }
    }

    // Subtriangle (C,A,Q): (1,0), (0,0), (u0,v0)
    t1.setOrder(n);
    t1.makeNonrational();
    for(int i = 0; i <= n; ++i)
    {
      const int deg = n - i;
      for(int j = 0; j <= n - i; ++j)
      {
        t1(i, j) = net[i][triIndex(deg, deg - j, 0)];
      }
    }

    // Subtriangle (A,B,Q): (0,0), (0,1), (u0,v0)
    t2.setOrder(n);
    t2.makeNonrational();
    for(int i = 0; i <= n; ++i)
    {
      const int deg = n - i;
      for(int j = 0; j <= n - i; ++j)
      {
        t2(i, j) = net[i][triIndex(deg, 0, j)];
      }
    }
  }

  /*!
   * \brief Splits a Bezier triangle into two subtriangles by connecting a
   * point on a boundary edge to the opposite vertex
   *
   * \param [in] edgeIdx Index of the boundary edge to split (same convention as `getEdge(int)`)
   * \param [in] s Parameter in \a (0,1) locating the split point along the chosen edge
   * \param [out] t0 First output subtriangle
   * \param [out] t1 Second output subtriangle
   *
   * \pre edgeIdx is 0, 1, or 2
   * \pre \a s is in (0,1)
   *
   * Taking P0 as the vertex opposite edge `edgeIdx`:
   *
   *            P0
   *           /|\
   *          / | \
   *         /  |  \
   *        /   |   \
   *       / t0 | t1 \
   *      /_____|_____\
   *   P1       Q      P2
   *   s=0      s      s=1
   *
   * \return t0 = Tri( P0, P1, Q ), t1 = Tri( P2, P0, Q )
   */
  void split(int edgeIdx, T s, BezierTriangle& t0, BezierTriangle& t1) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(edgeIdx >= 0 && edgeIdx < 3);
    SLIC_ASSERT(s > T(0) && s < T(1));

    if(isRational())
    {
      // Rational case: split in projective space and for weights, then convert back.
      BezierTriangle<T, NDIMS> projective(m_ord);
      BezierTriangle<T, 1> weights(m_ord);
      fill_projective_triangles(projective, weights);

      BezierTriangle<T, NDIMS> p0, p1;
      BezierTriangle<T, 1> w0, w1;
      projective.split(edgeIdx, s, p0, p1);
      weights.split(edgeIdx, s, w0, w1);

      set_from_projective_triangles(p0, w0, t0);
      set_from_projective_triangles(p1, w1, t1);
      return;
    }

    using TriangularArray = axom::Array<PointType>;

    Barycentric Q;
    switch(edgeIdx)
    {
    case 0:  // (u=0, v=1-s, w=s)
      Q = Barycentric {T(0), T(1) - s, s};
      break;
    case 1:  // (u=s, v=0, w=1-s)
      Q = Barycentric {s, T(0), T(1) - s};
      break;
    case 2:  // (u=1-s, v=s, w=0)
      Q = Barycentric {T(1) - s, s, T(0)};
      break;
    default:
      break;
    }

    const int n = m_ord;
    axom::Array<TriangularArray> net(n + 1);
    net[0] = m_controlPoints;

    for(int p = 1; p <= n; ++p)
    {
      const int deg = n - p + 1;
      const int newDeg = deg - 1;

      net[p].resize(triSize(newDeg));
      const auto& prev = net[p - 1];

      for(int ii = 0; ii <= newDeg; ++ii)
      {
        for(int jj = 0; jj <= newDeg - ii; ++jj)
        {
          const auto& A0 = prev[triIndex(deg, ii, jj)];
          const auto& B0 = prev[triIndex(deg, ii, jj + 1)];
          const auto& C0 = prev[triIndex(deg, ii + 1, jj)];

          net[p][triIndex(newDeg, ii, jj)] = triInterpolate(A0, B0, C0, Q);
        }
      }
    }

    auto fill_from_net = [&](int keptEdge, BezierTriangle& out) {
      out.setOrder(n);
      out.makeNonrational();
      for(int i = 0; i <= n; ++i)
      {
        const int deg = n - i;
        for(int j = 0; j <= n - i; ++j)
        {
          switch(keptEdge)
          {
          case 0:
            // Keeps original BC line: (j,deg-j)
            out(i, j) = net[i][triIndex(deg, j, deg - j)];
            break;
          case 1:
            // Keeps original CA line: (deg-j,0)
            out(i, j) = net[i][triIndex(deg, deg - j, 0)];
            break;
          case 2:
            // Keeps original AB line: (0,j)
            out(i, j) = net[i][triIndex(deg, 0, j)];
            break;
          default:
            break;
          }
        }
      }
    };

    switch(edgeIdx)
    {
    case 0:
      // Q on BC -> keep AB and CA
      fill_from_net(2, t0);  // (A,B,Q)
      fill_from_net(1, t1);  // (C,A,Q)
      break;
    case 1:
      // Q on CA -> keep BC and AB
      fill_from_net(0, t0);  // (B,C,Q)
      fill_from_net(2, t1);  // (A,B,Q)
      break;
    case 2:
      // Q on AB -> keep CA and BC
      fill_from_net(1, t0);  // (C,A,Q)
      fill_from_net(0, t1);  // (B,C,Q)
      break;
    default:
      break;
    }
  }

  /*!
   * \brief Uniform 4-way split at edge midpoints
   *
   * \param [out] t0 Subtriangle near vertex `(0,0)`
   * \param [out] t1 Subtriangle near vertex `(0,1)`
   * \param [out] t2 Subtriangle near vertex `(1,0)`
   * \param [out] t3 Central subtriangle
   *
   * This is equivalent to `split(0.5, 0.5, 0.5, ...)`, but with optimizations to reduce
   *  redundant computations by sharing intermediate control nets.
   * We also separate the implementation based on triangle rationality to improve performance
   *
   *           C
   *           /\
   *          /t2\
   *      P1 /____\ P0
   *        /\ t3 /\
   *       /t0\  /t1\
   *      /____\/____\
   *     A     P2      B
   *
   * \return t0 = Tri( A, P2, P1 ), t1 = Tri( B, P0, P2 ), t2 = Tri( C, P1, P0 ), t3 = Tri( P1, P0, P2 )
   */
  void uniformSplit(BezierTriangle& t0, BezierTriangle& t1, BezierTriangle& t2, BezierTriangle& t3) const
  {
    SLIC_ASSERT(m_ord >= 0);
    const int n = m_ord;

    t0.setOrder(n);
    t1.setOrder(n);
    t2.setOrder(n);
    t3.setOrder(n);

    if(!isRational())
    {
      t0.makeNonrational();
      t1.makeNonrational();
      t2.makeNonrational();
      t3.makeNonrational();

      // For polynomial triangles, these accessors just use the regular control points
      auto get_point = [&](int i, int j) -> PointType { return (*this)(i, j); };
      auto set_point = [&](BezierTriangle& out, int i, int j, const PointType& pt) {
        out(i, j) = pt;
      };
      uniform_split_impl<PointType>(get_point, set_point, t0, t1, t2, t3);
    }
    else
    {
      using HomogeneousPoint = Point<T, NDIMS + 1>;
      t0.makeRational();
      t1.makeRational();
      t2.makeRational();
      t3.makeRational();

      // For rational triangles, these accessors generate homogeneous control points
      auto get_hom_point = [&](int i, int j) -> HomogeneousPoint {
        const int idx = triIndex(n, i, j);

        HomogeneousPoint hp;
        hp[NDIMS] = m_weights[idx];
        for(int d = 0; d < NDIMS; ++d)
        {
          hp[d] = m_controlPoints[idx][d] * m_weights[idx];
        }
        return hp;
      };
      auto set_hom_point = [&](BezierTriangle& out, int i, int j, const HomogeneousPoint& hp) {
        const int idx = triIndex(n, i, j);

        const T w = hp[NDIMS];
        SLIC_ASSERT(w > T(0));

        out.m_weights[idx] = w;
        for(int d = 0; d < NDIMS; ++d)
        {
          out.m_controlPoints[idx][d] = hp[d] / w;
        }
      };
      uniform_split_impl<HomogeneousPoint>(get_hom_point, set_hom_point, t0, t1, t2, t3);
    }
  }

  /*!
   * \brief Splits a Bezier triangle into four subtriangles by inserting one
   * split point on each boundary edge and connecting the split points pairwise
   *
   * \param [in] s1 Parameter in \a (0,1) locating the split point on edge 0 (same convention as `getEdge(0)`)
   * \param [in] s2 Parameter in \a (0,1) locating the split point on edge 1 (same convention as `getEdge(1)`)
   * \param [in] s3 Parameter in \a (0,1) locating the split point on edge 2 (same convention as `getEdge(2)`)
   * \param [out] t1 Subtriangle near vertex `(0,0)` with vertices `(0,0)`, point on edge 2, point on edge 1
   * \param [out] t2 Subtriangle near vertex `(0,1)` with vertices `(0,1)`, point on edge 0, point on edge 2
   * \param [out] t3 Subtriangle near vertex `(1,0)` with vertices `(1,0)`, point on edge 1, point on edge 0
   * \param [out] t4 Central subtriangle with vertices point on edge 1, point on edge 0, point on edge 2
   *
   * \pre \a s1, \a s2, \a s3 are all in (0,1)
   *           C
   *           /\
   *          /t3\
   *      P1 /____\ P0
   *        /\ t4 /\
   *       /t1\  /t2\
   *      /____\/____\
   *     A     P2      B
   *
   * \return t1 = Tri( A, P2, P1 ), t2 = Tri( B, P0, P2 ), t3 = Tri( C, P1, P0 ), t4 = Tri( P1, P0, P2 )
   */
  void split(T s1,
             T s2,
             T s3,
             BezierTriangle& t1,
             BezierTriangle& t2,
             BezierTriangle& t3,
             BezierTriangle& t4) const
  {
    SLIC_ASSERT(m_ord >= 0);
    SLIC_ASSERT(s1 > T(0) && s1 < T(1));
    SLIC_ASSERT(s2 > T(0) && s2 < T(1));
    SLIC_ASSERT(s3 > T(0) && s3 < T(1));

    // Standard Barycentric coordinates in {u,v,w} where w = 1-u-v
    // and the triangle vertices are: A={1,0,0}, B={0,1,0}, C={0,0,1}
    const Barycentric A {T(1), T(0), T(0)};
    const Barycentric B {T(0), T(1), T(0)};
    const Barycentric C {T(0), T(0), T(1)};

    // Edge points {u,v,w} following the same orientation as getEdge(0..2)
    // edge0: B->C  (w == 0)
    // edge1: C->A  (v == 0)
    // edge2: A->B  (u == 0)
    const Barycentric P0 {T(0), T(1) - s1, s1};
    const Barycentric P1 {s2, T(0), T(1) - s2};
    const Barycentric P2 {T(1) - s3, s3, T(0)};

    // Corner triangles (ordered to match the diagram and preserve interior-edge orientation)
    restrictToSubtriangle(A, P2, P1, t1);
    restrictToSubtriangle(B, P0, P2, t2);
    restrictToSubtriangle(C, P1, P0, t3);
    restrictToSubtriangle(P1, P0, P2, t4);
  }

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
    SLIC_ASSERT(isValidIndex(m_ord, i, j));
    return m_controlPoints[triIndex(m_ord, i, j)];
  }

  /// \brief Access a control point in the triangular control net
  const PointType& operator()(int i, int j) const
  {
    SLIC_ASSERT(isValidIndex(m_ord, i, j));
    return m_controlPoints[triIndex(m_ord, i, j)];
  }

  /*!
   * \brief Get a specific weight from a rational Bezier triangle
   *
   * \param [in] i First control net index
   * \param [in] j Second control net index
   * \pre Requires that the triangle be rational
   */
  const T& getWeight(int i, int j) const
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(isValidIndex(m_ord, i, j));
    return m_weights[triIndex(m_ord, i, j)];
  }

  /*!
   * \brief Set the weight at a specific index for a rational Bezier triangle
   *
   * \param [in] i First control net index
   * \param [in] j Second control net index
   * \param [in] weight The updated value of the weight
   * \pre Requires that the triangle be rational
   * \pre Requires that the weight be positive
   */
  void setWeight(int i, int j, T weight)
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(weight > T {0});
    SLIC_ASSERT(isValidIndex(m_ord, i, j));
    m_weights[triIndex(m_ord, i, j)] = weight;
  }

  /*!
   * \brief Evaluates the Bezier triangle at \a (u0,v0)
   *
   * \param [in] u0 Parameter value along the \a u axis
   * \param [in] v0 Parameter value along the \a v axis
   *
   * \note Evaluation uses permuted barycentric coordinates such that
   *   parameter values (u0, v0) correspond to the triangle vertices:
   *     - `evaluate(0,0) == (*this)(0,0)`
   *     - `evaluate(0,1) == (*this)(0,order)`
   *     - `evaluate(1,0) == (*this)(order,0)`
   *
   * \warning Will automatically extrapolate if (u0, v0) are outside the triangular
   *            domain (u0 >= 0, v0 >= 0, and u0+v0 <= 1)
   * 
   * \return Point value S(u0,v0)
   */
  PointType evaluate(T u0, T v0) const
  {
    SLIC_ASSERT(m_ord >= 0);

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
              dCarray[triIndex(end - 1, i, j)] =
                triInterpolate(A, B, C, Barycentric {1.0 - u0 - v0, v0, u0});
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

      const Point<T, NDIMS> P = projective.evaluate(u0, v0);
      const Point<T, 1> W = weights.evaluate(u0, v0);

      PointType eval;
      for(int N = 0; N < NDIMS; ++N)
      {
        eval[N] = P[N] / W[0];
      }
      return eval;
    }
  }

  /*!
   * \brief Evaluates first derivatives of the Bezier triangle at \a (u0,v0)
   *
   * \param [in] u0 Parameter value along the \a u axis
   * \param [in] v0 Parameter value along the \a v axis
   * \param [out] eval Point value S(u,v)
   * \param [out] Du First derivative S_u(u,v)
   * \param [out] Dv First derivative S_v(u,v)
   *
   * \warning Will automatically extrapolate if (u0, v0) are outside the triangular
   *            domain (u0 >= 0, v0 >= 0, and u0+v0 <= 1)
   */
  void evaluateFirstDerivatives(T u0,
                                T v0,
                                Point<T, NDIMS>& eval,
                                Vector<T, NDIMS>& Du,
                                Vector<T, NDIMS>& Dv) const
  {
    SLIC_ASSERT(m_ord >= 0);

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
              dCarray[triIndex(end - 1, i, j)] =
                triInterpolate(A, B, C, Barycentric {1.0 - u0 - v0, v0, u0});
            }
          }
        }

        // The last reduction yields a linear triangle:
        //   S(u,v) = A + u0(C-A) + v0(B-A)
        Du[N] = (dCarray[2] - dCarray[0]);
        Dv[N] = (dCarray[1] - dCarray[0]);
        eval[N] = dCarray[0] + u0 * Du[N] + v0 * Dv[N];

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

      projective.evaluateFirstDerivatives(u0, v0, P, P_u, P_v);
      weights.evaluateFirstDerivatives(u0, v0, W, W_u, W_v);

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
   * \param [in] u0 Parameter value at which to evaluate along the u axis
   * \param [in] v0 Parameter value at which to evaluate along the v axis
   * \param [out] eval The point value of the Bezier triangle at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   *
   * \warning Will automatically extrapolate if (u0, v0) are outside the triangular
   *            domain (u0 >= 0, v0 >= 0, and u0+v0 <= 1)
   */
  void evaluateLinearDerivatives(T u0,
                                 T v0,
                                 Point<T, NDIMS>& eval,
                                 Vector<T, NDIMS>& Du,
                                 Vector<T, NDIMS>& Dv,
                                 Vector<T, NDIMS>& DuDv) const
  {
    SLIC_ASSERT(m_ord >= 0);

    if(!isRational())
    {
      if(m_ord < 2)
      {
        evaluateFirstDerivatives(u0, v0, eval, Du, Dv);
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
              dCarray[triIndex(end - 1, i, j)] =
                triInterpolate(A, B, C, Barycentric {1.0 - u0 - v0, v0, u0});
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
        const T L00 = Q00 + u0 * (Q10 - Q00) + v0 * (Q01 - Q00);
        const T L01 = Q01 + u0 * (Q11 - Q01) + v0 * (Q02 - Q01);
        const T L10 = Q10 + u0 * (Q20 - Q10) + v0 * (Q11 - Q10);

        eval[N] = L00 + u0 * (L10 - L00) + v0 * (L01 - L00);
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

      projective.evaluateLinearDerivatives(u0, v0, P, P_u, P_v, P_uv);
      weights.evaluateLinearDerivatives(u0, v0, W, W_u, W_v, W_uv);

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
   * \brief Evaluates all second derivatives of a Bezier triangle at (\a u0, \a v0)
   *
   * \param [in] u0 Parameter value at which to evaluate along the u axis
   * \param [in] v0 Parameter value at which to evaluate along the v axis
   * \param [out] eval The point value of the Bezier triangle at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDu The vector value of S_uu(u, v)
   * \param [out] DvDv The vector value of S_vv(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   * 
   * \warning Will automatically extrapolate if (u0, v0) are outside the triangular
   *            domain (u0 >= 0, v0 >= 0, and u0+v0 <= 1)
   */
  void evaluateSecondDerivatives(T u0,
                                 T v0,
                                 Point<T, NDIMS>& eval,
                                 Vector<T, NDIMS>& Du,
                                 Vector<T, NDIMS>& Dv,
                                 Vector<T, NDIMS>& DuDu,
                                 Vector<T, NDIMS>& DvDv,
                                 Vector<T, NDIMS>& DuDv) const
  {
    SLIC_ASSERT(m_ord >= 0);

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
      evaluateFirstDerivatives(u0, v0, eval, Du, Dv);
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
              dCarray[triIndex(end - 1, i, j)] =
                triInterpolate(A, B, C, Barycentric {1.0 - u0 - v0, v0, u0});
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
        const T L00 = Q00 + u0 * (Q10 - Q00) + v0 * (Q01 - Q00);
        const T L01 = Q01 + u0 * (Q11 - Q01) + v0 * (Q02 - Q01);
        const T L10 = Q10 + u0 * (Q20 - Q10) + v0 * (Q11 - Q10);

        eval[N] = L00 + u0 * (L10 - L00) + v0 * (L01 - L00);
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

      projective.evaluateSecondDerivatives(u0, v0, P, P_u, P_v, P_uu, P_vv, P_uv);
      weights.evaluateSecondDerivatives(u0, v0, W, W_u, W_v, W_uu, W_vv, W_uv);

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

  /// \brief Computes a tangent of a Bezier triangle at (\a u0, \a v0) along the u axis
  VectorType du(T u0, T v0) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u0, v0, eval, Du, Dv);
    return Du;
  }

  /// \brief Computes a tangent of a Bezier triangle at (\a u0, \a v0) along the v axis
  VectorType dv(T u0, T v0) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u0, v0, eval, Du, Dv);
    return Dv;
  }

  /// \brief Computes the second derivative of a Bezier triangle at (\a u0, \a v0) along the u axis
  VectorType dudu(T u0, T v0) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u0, v0, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDu;
  }

  /// \brief Computes the second derivative of a Bezier triangle at (\a u0, \a v0) along the v axis
  VectorType dvdv(T u0, T v0) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u0, v0, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DvDv;
  }

  /// \brief Computes the mixed second derivative of a Bezier triangle at (\a u0, \a v0)
  VectorType dudv(T u0, T v0) const
  {
    PointType eval;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    evaluateSecondDerivatives(u0, v0, eval, Du, Dv, DuDu, DvDv, DuDv);
    return DuDv;
  }

  /// \brief Convenience alias for S_vu(u0,v0), which equals S_uv(u0,v0) for polynomial triangles
  VectorType dvdu(T u0, T v0) const { return dudv(u0, v0); }

  /*!
   * \brief Computes the normal vector of a Bezier triangle at (\a u0, \a v0)
   *
   * \note Only meaningful for NDIMS==3.
   */
  VectorType normal(T u0, T v0) const
  {
    Point<T, NDIMS> eval;
    Vector<T, NDIMS> Du, Dv;
    evaluateFirstDerivatives(u0, v0, eval, Du, Dv);
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
  static constexpr int triSize(int ord) { return (ord >= 0) ? ((ord + 1) * (ord + 2) / 2) : 0; }

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

  /// \brief Check if a given index is valid in the triangular array
  static constexpr bool isValidIndex(int ord, int i, int j)
  {
    return (i >= 0) && (j >= 0) && (i + j <= ord);
  }

  /// \brief Do a triangular interpolation from three points and a barycentric coordinate
  static constexpr PointType triInterpolate(const PointType& A,
                                            const PointType& B,
                                            const PointType& C,
                                            const Barycentric& Q)
  {
    return PointType {Q[0] * A.array() + Q[1] * B.array() + Q[2] * C.array()};
  }

  /// \brief Do a triangular interpolation from three coordinates and a barycentric coordinate
  static constexpr T triInterpolate(const T& A, const T& B, const T& C, const Barycentric& Q)
  {
    return Q[0] * A + Q[1] * B + Q[2] * C;
  }

private:
  /// \brief Check that each weight is positive and the size matches control nodes
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

  /// \brief Check if the given coordinates are interior to the triangle
  bool is_valid_interior_parameter(T u0, T v0) const
  {
    return (u0 > T {0}) && (v0 > T {0}) && (u0 + v0 < T {1});
  }

  /*!
   * \brief Private function to evaluate the uniform split algorithm
   *
   * \param [in] get_eval_point Lambda that returns an `EvalPointType` for control point `(i,j)`
   * \param [in] set_eval_point Lambda that assigns an `EvalPointType` to output triangle control point `(i,j)`
   *              If the triangle is polynomial, these should access the control point as-is.
   *              If the triangle is rational, these should access the homogeneous control point.
   *
   * Implements the algorithm from Kenneth I. Joy, "A Uniform Subdivision Method for Triangular Bezier Patches"
   *
   * We construct a degenerate "Bezier tetrahedron" of intermediate points via repeated
   *  midpoint averaging in the control net, then extract the four subtriangle control nets
   *  from its four faces (three corner triangles + one central triangle)
   *
   * We separate this templated implementation so that we can more efficiently process rational triangles.
   */
  template <typename EvalPointType, typename GetEvalPointFn, typename SetEvalPointFn>
  void uniform_split_impl(GetEvalPointFn get_eval_point,
                          SetEvalPointFn set_eval_point,
                          BezierTriangle& t0,
                          BezierTriangle& t1,
                          BezierTriangle& t2,
                          BezierTriangle& t3) const
  {
    const int n = m_ord;
    const int triN = triSize(n);
    constexpr int evalDims = EvalPointType::dimension();

    axom::Array<int> offsets(n + 2);
    offsets[0] = 0;
    for(int p = 0; p <= n; ++p)
    {
      offsets[p + 1] = offsets[p] + triSize(p);
    }
    const int tetN = offsets[n + 1];

    // tetSize: Number of (n1,n2,n3) triples with n1+n2+n3 <= m_ord.
    SLIC_ASSERT(tetN == (n + 1) * (n + 2) * (n + 3) / 6);

    // Storage for the degenerate tetrahedron:
    // - tetIdx indexes a (n1,n2,n3) triple with fixed sum p via:
    //     tetIdx(p,n1,n2) = offsets[p] + triIndex(p,n1,n2), with n3 = p-n1-n2
    // - within each state, we store a full triangle mesh at degree n indexed by triIndex(n,i,j)
    //   (even though only a subset of indices are valid for a given (n1,n2,n3)).
    axom::Array<axom::Array<EvalPointType>> tet(tetN);
    for(int s = 0; s < tetN; ++s)
    {
      tet[s].resize(triN);
    }

    auto tetIdx = [&](int p, int n1, int n2) -> int {
      SLIC_ASSERT(p >= 0 && p <= n);
      SLIC_ASSERT(n1 >= 0 && n1 <= p);
      SLIC_ASSERT(n2 >= 0 && n2 <= p - n1);
      return offsets[p] + triIndex(p, n1, n2);
    };
    auto getTetPt = [&](int state, int i, int j) -> const EvalPointType& {
      return tet[static_cast<std::size_t>(state)][triIndex(n, i, j)];
    };
    auto setTetPt = [&](int state, int i, int j, const EvalPointType& value) {
      tet[static_cast<std::size_t>(state)][triIndex(n, i, j)] = value;
    };

    const int base = tetIdx(0, 0, 0);
    for(int i = 0; i <= n; ++i)
    {
      for(int j = 0; j <= n - i; ++j)
      {
        setTetPt(base, i, j, get_eval_point(i, j));
      }
    }

    // Build the tetrahedron states in increasing p = n1+n2+n3.
    for(int p = 1; p <= n; ++p)
    {
      for(int n1 = 0; n1 <= p; ++n1)
      {
        for(int n2 = 0; n2 <= p - n1; ++n2)
        {
          const int n3 = p - n1 - n2;
          const int curr = tetIdx(p, n1, n2);

          // Pick the predecessor state and the second point used in the midpoint average.
          int prev = -1, di = 0, dj = 0;
          if(n1 > 0)
          {
            prev = tetIdx(p - 1, n1 - 1, n2);
            di = -1;
            dj = +1;
          }
          else if(n2 > 0)
          {
            prev = tetIdx(p - 1, n1, n2 - 1);
            di = 0;
            dj = -1;
          }
          else
          {
            // n1==0 && n2==0 so n3>0 here.
            prev = tetIdx(p - 1, 0, 0);
            di = +1;
            dj = 0;
          }

          // Valid domain: i>=n1, j>=n2, k>=n3.
          for(int i = n1; i <= n; ++i)
          {
            const int jMax = n - n3 - i;
            if(jMax < n2)
            {
              continue;
            }
            for(int j = n2; j <= jMax; ++j)
            {
              const auto& a = getTetPt(prev, i, j);
              const auto& b = getTetPt(prev, i + di, j + dj);

              EvalPointType out;
              for(int d = 0; d < evalDims; ++d)
              {
                out[d] = T(0.5) * (a[d] + b[d]);
              }
              setTetPt(curr, i, j, out);
            }
          }
        }
      }
    }

    // Extract the four faces Q1..Q4, given by
    //   Q1 = { P^{[m,0,k]}_{i,0,k} : i+k=n, m=0..i }
    //   Q2 = { P^{[i,m,0]}_{i,j,0} : i+j=n, m=0..j }
    //   Q3 = { P^{[0,j,m]}_{0,j,k} : j+k=n, m=0..k }
    //   Q4 = { P^{[i,j,k]}_{i,j,k} : i+j+k=n }

    // Q3 -> t0 : P^{[0,j,m]}_{0,j,k} with j+k=n and m<=k maps to local (m,j)
    // so that (0,1) is the midpoint on AB and (1,0) is the midpoint on AC.
    for(int j = 0; j <= n; ++j)
    {
      const int i = 0;
      for(int m = 0; m <= n - j; ++m)
      {
        const int p = j + m;
        const int state = tetIdx(p, 0, j);
        set_eval_point(t0, m, j, getTetPt(state, i, j));
      }
    }

    // Q2 -> t1 : P^{[i,m,0]}_{i,j,0} with i+j=n and m<=j maps to local (m,i)
    // so that (0,1) is the midpoint on BC and (1,0) is the midpoint on AB.
    for(int i = 0; i <= n; ++i)
    {
      const int j = n - i;
      for(int m = 0; m <= j; ++m)
      {
        const int p = i + m;
        const int state = tetIdx(p, i, m);
        set_eval_point(t1, m, i, getTetPt(state, i, j));
      }
    }

    // Q1 -> t2 : P^{[m,0,k]}_{i,0,k} with i+k=n and m<=i maps to local (m,k)
    // so that (0,1) is the midpoint on AC and (1,0) is the midpoint on BC.
    for(int k = 0; k <= n; ++k)
    {
      const int i = n - k;
      const int j = 0;
      for(int m = 0; m <= i; ++m)
      {
        const int p = m + k;
        const int state = tetIdx(p, m, 0);
        set_eval_point(t2, m, k, getTetPt(state, i, j));
      }
    }

    // Q4 -> t3 : P^{[i,j,k]}_{i,j,k} with i+j+k=n maps to local (j,i)
    // so that (0,1) is the midpoint on BC and (1,0) is the midpoint on AB.
    for(int i = 0; i <= n; ++i)
    {
      for(int j = 0; j <= n - i; ++j)
      {
        const int state = offsets[n] + triIndex(n, j, i);
        set_eval_point(t3, i, j, getTetPt(state, j, i));
      }
    }
  }

  /// \brief Fill a (possibly non-rational) triangle from projective control points and weights
  static void set_from_projective_triangles(const BezierTriangle<T, NDIMS>& projective,
                                            const BezierTriangle<T, 1>& weights,
                                            BezierTriangle& out)
  {
    const int ord = projective.getOrder();
    SLIC_ASSERT(ord == weights.getOrder());

    out.setOrder(ord);
    out.makeRational();

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord - i; ++j)
      {
        const int idx = triIndex(ord, i, j);
        const T w = weights.m_controlPoints[idx][0];

        SLIC_ASSERT(w > T(0));
        out.m_weights[idx] = w;

        for(int N = 0; N < NDIMS; ++N)
        {
          out.m_controlPoints[idx][N] = projective.m_controlPoints[idx][N] / w;
        }
      }
    }
  }

  /// \brief For a rational triangle, fill triangle objects with weighted control points and weights
  void fill_projective_triangles(BezierTriangle<T, NDIMS>& projective,
                                 BezierTriangle<T, 1>& weights) const
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(is_valid_rational());

    for(int i = 0; i <= m_ord; ++i)
    {
      for(int j = 0; j <= m_ord - i; ++j)
      {
        const int idx = triIndex(m_ord, i, j);
        const T w = m_weights[idx];
        weights.m_controlPoints[idx][0] = w;

        for(int N = 0; N < NDIMS; ++N)
        {
          projective.m_controlPoints[idx][N] = m_controlPoints[idx][N] * w;
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
