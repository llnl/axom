// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file GregoryPatch.hpp
 *
 * \brief A bicubic Gregory patch primitive
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
template <typename T>
class GregoryPatch;

/*! \brief Overloaded output operator for Gregory Patches*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const GregoryPatch<T>& nPatch);

/*!
 * \class GregoryPatch
 *
 * \brief Represents a 3D Gregory patch defined by the control points of 4 cubic Bezier curves
 *          for each boundary, and an additional two "Gregory points" for each edge which 
 *          determine the internal geometry of the surface.
 * \tparam T the coordinate type, e.g., double, float, etc.
 */
template <typename T>
class GregoryPatch
{
public:
  // The number of control points for a bicubic Gregory patch is fixed:
  //  - 12 exterior control points
  //  - 2 interior control points for each of 4 boundary curves
  static constexpr int NPTS = 20;

  using PointType = Point<T, 3>;
  using VectorType = Vector<T, 3>;

  using CoordsVec = axom::StackArray<PointType, NPTS>;

  using BoundingBoxType = BoundingBox<T, 3>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, 3>;

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A Gregory Patch must be defined using an arithmetic type");

public:
  ///@{
  /**
   * @name Constructors for GregoryPatch
   *
   * The constructors allow initialization from:
   * - the 20 Gregory patch control points,
   * - a polynomial bicubic Bezier patch,
   * - C-style arrays, Axom StackArrays, or Axom ArrayViews,
   * - four corner positions with associated corner normal vectors.
   *
   * The 20-point control net is stored as:
   * - indices 0-3: corners,
   * - indices 4-11: two boundary control points for each edge,
   * - indices 12-19: two Gregory tangent points for each edge.
   *
   * Boundary edge \a e is directed from corner \a e to corner `(e+1)%4`.
   */

  /*!
   * \brief Default constructor for a Gregory patch
   *
   * The fixed-size control net is default-initialized.
   */
  GregoryPatch() = default;

  /*!
   * \brief Constructor from an ArrayView over the control points
   *
   * \param [in] controlPoints ArrayView of the 20 Gregory patch control points
   * \pre \a controlPoints must contain exactly `NPTS` points
   */
  explicit GregoryPatch(ArrayView<const PointType> controlPoints)
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
   * \param [in] controlPoints ArrayView of the 20 Gregory patch control points
   * \pre \a controlPoints must contain exactly `NPTS` points
   */
  explicit GregoryPatch(ArrayView<PointType> controlPoints)
    : GregoryPatch(ArrayView<const PointType>(controlPoints.data(), controlPoints.size()))
  { }

  /*!
   * \brief Constructor from a C-style array of control points
   *
   * \param [in] pts A C-style array of 20 Gregory patch control points
   * \pre \a pts must be non-null and contain at least `NPTS` points
   */
  explicit GregoryPatch(const PointType* pts) : GregoryPatch(ArrayView<const PointType>(pts, NPTS))
  { }

  /*!
   * \brief Constructor from a C-style array of control points
   *
   * \param [in] pts A C-style array of 20 Gregory patch control points
   * \pre \a pts must be non-null and contain at least `NPTS` points
   */
  explicit GregoryPatch(PointType* pts) : GregoryPatch(ArrayView<const PointType>(pts, NPTS)) { }

  /*!
   * \brief Constructor from an Axom StackArray of control points
   *
   * \param [in] pts StackArray containing the 20 Gregory patch control points
   */
  explicit GregoryPatch(const CoordsVec& pts)
    : GregoryPatch(ArrayView<const PointType>(pts.data(), pts.size()))
  { }

  /*!
   * \brief Constructor from a polynomial bicubic Bezier patch
   *
   * \param [in] bPatch A polynomial Bezier patch of order (3, 3)
   *
   * This creates a Gregory patch that exactly reproduces the input bicubic Bezier patch.
   * The Gregory tangent pairs are duplicated from the four Bezier interior control points,
   * causing the parameter-dependent Gregory blends to collapse to fixed Bezier points.
   *
   * \pre \a bPatch must have order (3, 3)
   * \pre \a bPatch must be polynomial, not rational
   */
  explicit GregoryPatch(const BezierPatch<T, 3>& bPatch)
  {
    SLIC_ASSERT(bPatch.getOrder_u() == 3);
    SLIC_ASSERT(bPatch.getOrder_v() == 3);
    SLIC_ASSERT(!bPatch.isRational());

    getCorner(0) = bPatch(0, 0);
    getCorner(1) = bPatch(3, 0);
    getCorner(2) = bPatch(3, 3);
    getCorner(3) = bPatch(0, 3);

    getBoundaryPoint(0, 1) = bPatch(1, 0);
    getBoundaryPoint(0, 2) = bPatch(2, 0);
    getBoundaryPoint(1, 1) = bPatch(3, 1);
    getBoundaryPoint(1, 2) = bPatch(3, 2);
    getBoundaryPoint(2, 1) = bPatch(2, 3);
    getBoundaryPoint(2, 2) = bPatch(1, 3);
    getBoundaryPoint(3, 1) = bPatch(0, 2);
    getBoundaryPoint(3, 2) = bPatch(0, 1);

    getTangent(0, 0) = bPatch(1, 1);
    getTangent(0, 1) = bPatch(2, 1);
    getTangent(1, 0) = bPatch(2, 1);
    getTangent(1, 1) = bPatch(2, 2);
    getTangent(2, 0) = bPatch(2, 2);
    getTangent(2, 1) = bPatch(1, 2);
    getTangent(3, 0) = bPatch(1, 2);
    getTangent(3, 1) = bPatch(1, 1);
  }

  /*!
   * \brief Constructor from corner nodes and corner normal vectors
   *
   * \param [in] nodePositions ArrayView of the four corner positions
   * \param [in] nodeVectors ArrayView of the four corner normal vectors
   *
   * Deterministically compute cubic boundary control points and Gregory tangent points
   * using local corner information.
   *
   * \pre \a nodePositions and \a nodeVectors must each contain exactly 4 entries
   */
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

  ///@}

  /// \brief Returns the \a i-th corner point, oriented ccw
  PointType& getCorner(int i) { return m_controlPoints[i]; }

  /// \brief Returns the \a i-th corner point, oriented ccw
  const PointType& getCorner(int i) const { return m_controlPoints[i]; }

  /*!
   * \brief Returns a Gregory tangent point for an edge
   *
   * \param [in] e Edge index, oriented ccw
   * \param [in] t Tangent point index (either 0 or 1)
   */
  PointType& getTangent(int e, int t) { return m_controlPoints[12 + 2 * e + t]; }

  /*!
   * \brief Returns a Gregory tangent point for an edge
   *
   * \param [in] e Edge index, oriented ccw
   * \param [in] t Tangent point index (either 0 or 1)
   */
  const PointType& getTangent(int e, int t) const { return m_controlPoints[12 + 2 * e + t]; }

  /*!
   * \brief Returns the two Gregory tangent points adjacent to a corner
   *
   * \param [in] i Corner index, oriented ccw
   * \param [out] v0 Tangent point from the preceding edge
   * \param [out] v1 Tangent point from the following edge
   */
  void getTangentsByCorner(int i, PointType& v0, PointType& v1) const
  {
    v0 = getTangent((i + 3) % 4, 1);
    v1 = getTangent(i, 0);
  }

  /*!
   * \brief Returns a control point on a boundary edge
   *
   * \param [in] e Edge index in `[0, 3]`
   * \param [in] k Boundary point index in `[0, 3]`
   *
   * Values `k=0` and `k=3` are the edge's corner points; values `k=1` and `k=2`
   * are the cubic boundary control points.
   */
  PointType& getBoundaryPoint(int e, int k) { return m_controlPoints[s_edge_index_map[e][k]]; }

  /*!
   * \brief Returns a control point on a boundary edge
   *
   * \param [in] e Edge index in `[0, 3]`
   * \param [in] k Boundary point index in `[0, 3]`
   *
   * Values `k=0` and `k=3` are the edge's corner points; values `k=1` and `k=2`
   * are the cubic boundary control points.
   */
  const PointType& getBoundaryPoint(int e, int k) const
  {
    return m_controlPoints[s_edge_index_map[e][k]];
  }

  /*!
   * \brief Returns a reference to the patch's control points
   */
  CoordsVec& getControlPoints() { return m_controlPoints; }

  /// \brief Returns a reference to the patch's control points
  const CoordsVec& getControlPoints() const { return m_controlPoints; }

  /*!
   * \brief Evaluates the Gregory patch at the given parameter values
   *
   * \param [in] u Parameter value on the first axis
   * \param [in] v Parameter value on the second axis
   *
   * A cubic Gregory patch is evaluated by constructing the equivalent bicubic Bezier patch
   * whose interior control points are blended from the Gregory tangent points at (\a u, \a v).
   */
  PointType evaluate(T u, T v) const
  {
    const auto intermediate = setup_intermediate_bezier(u, v, 0);
    return intermediate.bpatch.evaluate(u, v);
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
  void evaluateFirstDerivatives(T u, T v, PointType& eval, VectorType& Du, VectorType& Dv) const
  {
    const auto intermediate = setup_intermediate_bezier(u, v, 1);
    intermediate.bpatch.evaluateFirstDerivatives(u, v, eval, Du, Dv);

    // Chain rule correction due to (u,v)-dependent interior control points
    axom::StaticArray<T, 4> Bu, dBu, Bv, dBv;
    evaluate_cubic_bernstein(u, Bu, dBu);
    evaluate_cubic_bernstein(v, Bv, dBv);

    const T w11 = Bu[1] * Bv[1];
    const T w21 = Bu[2] * Bv[1];
    const T w12 = Bu[1] * Bv[2];
    const T w22 = Bu[2] * Bv[2];

    Du += w11 * intermediate.Q_u[0][0] + w21 * intermediate.Q_u[1][0] +
      w12 * intermediate.Q_u[0][1] + w22 * intermediate.Q_u[1][1];
    Dv += w11 * intermediate.Q_v[0][0] + w21 * intermediate.Q_v[1][0] +
      w12 * intermediate.Q_v[0][1] + w22 * intermediate.Q_v[1][1];
  }

  /*!
   * \brief Evaluates all second derivatives of the Gregory patch at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the Gregory patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDu The vector value of S_uu(u, v)
   * \param [out] DvDv The vector value of S_vv(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   */
  void evaluateSecondDerivatives(T u,
                                 T v,
                                 PointType& eval,
                                 VectorType& Du,
                                 VectorType& Dv,
                                 VectorType& DuDu,
                                 VectorType& DvDv,
                                 VectorType& DuDv) const
  {
    const auto intermediate = setup_intermediate_bezier(u, v, 2);
    intermediate.bpatch.evaluateSecondDerivatives(u, v, eval, Du, Dv, DuDu, DvDv, DuDv);

    // Chain rule correction due to (u,v)-dependent interior control points
    axom::StaticArray<T, 4> Bu, dBu, Bv, dBv;
    evaluate_cubic_bernstein(u, Bu, dBu);
    evaluate_cubic_bernstein(v, Bv, dBv);

    const T w11 = Bu[1] * Bv[1];
    const T w21 = Bu[2] * Bv[1];
    const T w12 = Bu[1] * Bv[2];
    const T w22 = Bu[2] * Bv[2];

    const T wu11 = dBu[1] * Bv[1];
    const T wu21 = dBu[2] * Bv[1];
    const T wu12 = dBu[1] * Bv[2];
    const T wu22 = dBu[2] * Bv[2];

    const T wv11 = Bu[1] * dBv[1];
    const T wv21 = Bu[2] * dBv[1];
    const T wv12 = Bu[1] * dBv[2];
    const T wv22 = Bu[2] * dBv[2];

    // First derivative corrections
    Du += w11 * intermediate.Q_u[0][0] + w21 * intermediate.Q_u[1][0] +
      w12 * intermediate.Q_u[0][1] + w22 * intermediate.Q_u[1][1];
    Dv += w11 * intermediate.Q_v[0][0] + w21 * intermediate.Q_v[1][0] +
      w12 * intermediate.Q_v[0][1] + w22 * intermediate.Q_v[1][1];

    // Second derivative corrections
    DuDu += T(2) *
        (wu11 * intermediate.Q_u[0][0] + wu21 * intermediate.Q_u[1][0] +
         wu12 * intermediate.Q_u[0][1] + wu22 * intermediate.Q_u[1][1]) +
      (w11 * intermediate.Q_uu[0][0] + w21 * intermediate.Q_uu[1][0] +
       w12 * intermediate.Q_uu[0][1] + w22 * intermediate.Q_uu[1][1]);

    DvDv += T(2) *
        (wv11 * intermediate.Q_v[0][0] + wv21 * intermediate.Q_v[1][0] +
         wv12 * intermediate.Q_v[0][1] + wv22 * intermediate.Q_v[1][1]) +
      (w11 * intermediate.Q_vv[0][0] + w21 * intermediate.Q_vv[1][0] +
       w12 * intermediate.Q_vv[0][1] + w22 * intermediate.Q_vv[1][1]);

    DuDv += (wu11 * intermediate.Q_v[0][0] + wu21 * intermediate.Q_v[1][0] +
             wu12 * intermediate.Q_v[0][1] + wu22 * intermediate.Q_v[1][1]) +
      (wv11 * intermediate.Q_u[0][0] + wv21 * intermediate.Q_u[1][0] +
       wv12 * intermediate.Q_u[0][1] + wv22 * intermediate.Q_u[1][1]) +
      (w11 * intermediate.Q_uv[0][0] + w21 * intermediate.Q_uv[1][0] +
       w12 * intermediate.Q_uv[0][1] + w22 * intermediate.Q_uv[1][1]);
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

  /*!
   * \brief Simple formatted print of a Gregory Patch instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  void print(std::ostream& os) const
  {
    os << "GregoryPatch(";
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
  /*!
   * \brief Stores the temporary Bezier patch and blended interior point derivatives
   *
   * The Gregory patch evaluation converts the control net to a bicubic Bezier patch at
   * a specific parameter value. The four interior Bezier points, `Q`, depend on the
   * evaluation parameters, so derivative evaluation also requires their first and second
   * derivatives.
   */
  struct IntermediateBlendingDerivatives
  {
    /// \brief Equivalent bicubic Bezier patch for the requested parameter value
    BezierPatch<T, 3> bpatch;

    /// \brief Blended interior Bezier control points and derivatives
    PointType Q[2][2];
    VectorType Q_u[2][2];
    VectorType Q_v[2][2];
    VectorType Q_uu[2][2];
    VectorType Q_vv[2][2];
    VectorType Q_uv[2][2];
  };

  /*!
   * \brief Constructs the equivalent Bezier patch and parameter-dependent interior data
   *
   * \param [in] u Parameter value on the first axis
   * \param [in] v Parameter value on the second axis
   * \param [in] derivative_order Highest derivative order to compute, in `[0, 2]`
   *
   * The returned bicubic Bezier patch has the Gregory boundary control points copied
   * directly and the four interior control points blended from the Gregory tangent points.
   */
  IntermediateBlendingDerivatives setup_intermediate_bezier(T u, T v, int derivative_order) const
  {
    IntermediateBlendingDerivatives out;
    out.bpatch = get_bezier_boundary();

    const T um = T(1) - u;
    const T vm = T(1) - v;

    // For each interior point (i,j), select which corner's tangents to use.
    // This mapping matches the explicit construction used in evaluate():
    //   (0,0)->corner0, (1,0)->corner1, (0,1)->corner3, (1,1)->corner2
    static constexpr int corner_of[2][2] = {{0, 3}, {1, 2}};

    for(int j = 0; j < 2; ++j)
    {
      const T wb = (j == 0) ? v : vm;
      const T wb_v = (j == 0) ? T(1) : T(-1);

      for(int i = 0; i < 2; ++i)
      {
        const T wa = (i == 0) ? u : um;
        const T wa_u = (i == 0) ? T(1) : T(-1);

        const int corner = corner_of[i][j];

        PointType tPrev, tNext;
        getTangentsByCorner(corner, tPrev, tNext);

        // Corner parity determines which tangent is blended with which weight
        // (see explicit mapping in evaluate()).
        const bool swap = (corner % 2 == 0);  // corners 0 and 2
        const PointType& A = swap ? tNext : tPrev;
        const PointType& B = swap ? tPrev : tNext;

        const T denom = wa + wb;
        if(axom::utilities::isNearlyEqual(denom, T(0)))
        {
          out.Q[i][j] = A;
          out.Q_u[i][j] = VectorType(T(0));
          out.Q_v[i][j] = VectorType(T(0));
          out.Q_uu[i][j] = VectorType(T(0));
          out.Q_vv[i][j] = VectorType(T(0));
          out.Q_uv[i][j] = VectorType(T(0));
          continue;
        }

        out.Q[i][j] = PointType((wa * A.array() + wb * B.array()) / denom);

        if(derivative_order >= 1)
        {
          // Since wa depends only on u and wb depends only on v:
          //   dQ/du uses only wa_u; dQ/dv uses only wb_v.
          out.Q_u[i][j] = VectorType((wa_u * (A.array() - out.Q[i][j].array())) / denom);
          out.Q_v[i][j] = VectorType((wb_v * (B.array() - out.Q[i][j].array())) / denom);
        }
        else
        {
          out.Q_u[i][j] = VectorType(T(0));
          out.Q_v[i][j] = VectorType(T(0));
        }

        if(derivative_order >= 2)
        {
          // Second derivatives for linear weights (wa, wb):
          //   Q_uu = -(2*wa_u/denom) * Q_u
          //   Q_vv = -(2*wb_v/denom) * Q_v
          //   Q_uv = -(wa_u*Q_v + wb_v*Q_u)/denom
          out.Q_uu[i][j] = (-T(2) * wa_u / denom) * out.Q_u[i][j];
          out.Q_vv[i][j] = (-T(2) * wb_v / denom) * out.Q_v[i][j];
          out.Q_uv[i][j] = (-(wa_u * out.Q_v[i][j] + wb_v * out.Q_u[i][j])) / denom;
        }
        else
        {
          out.Q_uu[i][j] = VectorType(T(0));
          out.Q_vv[i][j] = VectorType(T(0));
          out.Q_uv[i][j] = VectorType(T(0));
        }
      }
    }

    set_bezier_interior(out.bpatch, out.Q);
    return out;
  }

  /*!
   * \brief Assigns the four interior control points of a bicubic Bezier patch
   *
   * \param [in,out] bpatch The bicubic Bezier patch to update
   * \param [in] Q The four interior control points, indexed by interior u and v position
   */
  static void set_bezier_interior(BezierPatch<T, 3>& bpatch, const PointType Q[2][2])
  {
    bpatch(1, 1) = Q[0][0];
    bpatch(2, 1) = Q[1][0];
    bpatch(1, 2) = Q[0][1];
    bpatch(2, 2) = Q[1][1];
  }

  /*!
   * \brief Copies the Gregory boundary into a bicubic Bezier patch
   *
   * The returned patch has its 12 exterior control points initialized from the Gregory
   * patch boundary. The four interior control points are intentionally left uninitialized.
   */
  BezierPatch<T, 3> get_bezier_boundary() const
  {
    BezierPatch<T, 3> bpatch(3, 3);

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

  /*!
   * \brief Evaluates cubic Bernstein basis functions and their first derivatives
   *
   * \param [in] t Parameter value
   * \param [out] B Cubic Bernstein basis values at \a t
   * \param [out] dB First derivative values of the cubic Bernstein basis at \a t
   */
  static void evaluate_cubic_bernstein(T t, axom::StaticArray<T, 4>& B, axom::StaticArray<T, 4>& dB)
  {
    const T tm = T(1) - t;
    const T tm2 = tm * tm;
    const T t2 = t * t;

    B[0] = tm2 * tm;
    B[1] = T(3) * t * tm2;
    B[2] = T(3) * t2 * tm;
    B[3] = t2 * t;

    dB[0] = -T(3) * tm2;
    dB[1] = T(3) * tm2 - T(6) * t * tm;
    dB[2] = T(6) * t * tm - T(3) * t2;
    dB[3] = T(3) * t2;
  }

  CoordsVec m_controlPoints;

  /*!
   * \brief Maps boundary curve control point indices to control net storage indices
   *
   * The first index selects a directed edge. The second index selects one of the four
   * cubic boundary control points on that edge.
   */
  static constexpr int s_edge_index_map[4][4] = {{/*V0*/ 0, /*E01*/ 4, /*E02*/ 5, /*V1*/ 1},
                                                 {/*V1*/ 1, /*E11*/ 6, /*E12*/ 7, /*V2*/ 2},
                                                 {/*V2*/ 2, /*E21*/ 8, /*E22*/ 9, /*V3*/ 3},
                                                 {/*V3*/ 3, /*E31*/ 10, /*E32*/ 11, /*V0*/ 0}};
};

//------------------------------------------------------------------------------
/// Free functions related to GregoryPatch
//------------------------------------------------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const GregoryPatch<T>& nPatch)
{
  nPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::GregoryPatch using fmt
template <typename T>
struct axom::fmt::formatter<axom::primal::GregoryPatch<T>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_GREGORY_PATCH_HPP_
