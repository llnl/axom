// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file NURBSPatch.hpp
 *
 * \brief A (trimmed) NURBSPatch primitive
 */

#ifndef AXOM_PRIMAL_NURBSPATCH_HPP_
#define AXOM_PRIMAL_NURBSPATCH_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/detail/winding_number_2d_impl.hpp"
#include "axom/primal/operators/detail/intersect_bezier_impl.hpp"

#include <ostream>
#include <math.h>

#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class NURBSPatch;

template <typename T>
class TrimmingCurveQuadratureData;

/*! \brief Overloaded output operator for NURBS Patches*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSPatch<T, NDIMS>& nPatch);

/*!
 * \class NURBSPatch
 *
 * \brief Represents a NURBS patch defined by a 2D array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * A NURBS patch has degrees `p` and `q` with knot vectors of length 
 * `r+1` and `s+1` respectively. There is a control net of (n + 1) * (m + 1) points
 * with r+1 = n+p+2 and s+1 = m+q+2. 
 * Optionally has an equal number of weights for rational patches.
 * 
 * The patch must be open (clamped on all boundaries) 
 *   and continuous (unless p = 0 or q = 0)
 * 
 * Nonrational NURBS patches are identified by an empty weights array.
 * 
 * A NURBS surface is identified as trimmed with an internal flag,
 *  as an empty trimming curve vector indicates a patch with no visible surface.
 */
template <typename T, int NDIMS>
class NURBSPatch
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using PlaneType = Plane<T, NDIMS>;

  using CoordsVec = axom::Array<PointType, 1>;
  using CoordsMat = axom::Array<PointType, 2>;
  using WeightsVec = axom::Array<T, 1>;
  using WeightsMat = axom::Array<T, 2>;
  using KnotVectorType = KnotVector<T>;

  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;
  using NURBSCurveType = primal::NURBSCurve<T, NDIMS>;

  using TrimmingCurveType = primal::NURBSCurve<T, 2>;
  using TrimmingCurveVec = axom::Array<TrimmingCurveType>;
  using ParameterPointType = Point<T, 2>;
  using ParameterBoundingBoxType = BoundingBox<T, 2>;

  AXOM_STATIC_ASSERT_MSG((NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
                         "A NURBS Patch object may be defined in 1-, 2-, or 3-D");

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A NURBS Patch must be defined using an arithmetic type");

public:
  ///@{
  /**
   * @name Constructors for NURBSPatch
   *
   * The NURBSPatch class provides a variety of constructors to support flexible initialization.
   * The default constructor creates an empty, invalid patch. Other constructors allow you to specify
   * degrees, control point counts, knot vectors, and weights, either using C-style arrays or axom::Array
   * containers, in both 1D and 2D forms. You can also construct a NURBSPatch from a BezierPatch.
   *
   * Depending on the constructor used, the resulting instance will have:
   * - Control points and weights arrays sized according to the provided parameters.
   * - Knot vectors initialized either as uniform (for degree/size-based constructors) or from provided arrays/objects.
   * - Rationality determined by the presence of weights (patches are rational if weights are provided, nonrational otherwise).
   * - Trimming curves are always empty and the patch is untrimmed by default.
   *
   * For 1D arrays, the mapping of control points and weights to the patch is lexicographical, i.e.
     \verbatim
      pts[0]                 -> nodes[0, 0],      ..., pts[npts_v]        -> nodes[0, npts_v]
      pts[npts_v+1]          -> nodes[1, 0],      ..., pts[2*npts_v]      -> nodes[1, npts_v]
                                                  ...
      pts[npts_u*(npts_v-1)] -> nodes[npts_u, 0], ..., pts[npts_u*npts_v] -> nodes[npts_u, npts_v]
     \endverbatim
   * 
   * All constructors ensure that the patch is internally consistent and valid, provided the input parameters are valid.
   */

  /*! 
   * \brief Default constructor for an empty (invalid) NURBS patch
   *
   * \note An empty NURBS patch is not valid
   */
  NURBSPatch() : NURBSPatch(0, 0, -1, -1) { }

  /*!
   * \brief Constructor for a simple NURBS surface that reserves space for
   *  the minimum (sensible) number of points for the given degrees
   * 
   * \param [in] deg_u, deg_v The patch's degrees on the first and second axis
   * \pre deg_u, deg_v both greater than or equal to 0, or both -1
   */
  NURBSPatch(int deg_u, int deg_v) : NURBSPatch(deg_u + 1, deg_v + 1, deg_u, deg_v) { }

  /*!
   * \brief Constructor for a simple NURBS surface that reserves space for
   *   \a npts_u * npts_v control points
   *
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] deg_u, deg_v The patch's degrees on the first and second axis
   * \pre Requires npts_d > deg_d and deg_d >= 0 for d = u, v 
   */
  NURBSPatch(int npts_u, int npts_v, int deg_u, int deg_v)
    : m_knotvec_u(npts_u, deg_u)
    , m_knotvec_v(npts_v, deg_v)
  {
    if(const bool is_empty = (deg_u == -1 && deg_v == -1); is_empty)
    {
      SLIC_ASSERT(npts_u == 0 && npts_v == 0);
    }
    else
    {
      SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);
      SLIC_ASSERT(npts_u > deg_u && npts_v > deg_v);
      m_controlPoints.resize(npts_u, npts_v);
    }
  }

  /*!
   * \brief Constructor for a NURBS surface from a Bezier surface
   *
   * \param [in] bezierPatch the Bezier patch to convert to a NURBS patch 
   */
  explicit NURBSPatch(const BezierPatch<T, NDIMS>& bezierPatch)
    : NURBSPatch(bezierPatch.getControlPoints(),
                 bezierPatch.getWeights(),
                 KnotVectorType(bezierPatch.getOrder_u() + 1, bezierPatch.getOrder_u()),
                 KnotVectorType(bezierPatch.getOrder_v() + 1, bezierPatch.getOrder_v()))
  { }

  /*!
   * \brief Constructor for a NURBSPatch from 2D ArrayViews of control points and weights 
   *        and KnotVectors for the u- and v- directions
   *
   * \param [in] controlPoints 2D ArrayView of control points
   * \param [in] weights 2D ArrayView of weights
   * \param [in] knotVector_u, knotVector_v The knot vectors in the u- and v- directions
   *
   * \pre If controlPoints is not empty, its sizes must match the number of control points 
   * implied by the knot vectors: knotVector_u.getNumControlPoints() * knotVector_v.getNumberControlPoints()
   * \pre If weights is not empty, its dimensions must match that of the controlPoints
   * \pre The KnotVector degrees must be valid: \a knotVector_u.getDegree() >= -1 
   * and \a knotVector_v.getDegree() >= -1. 
   * \pre If the degrees are not both -1 (i.e., not empty), then:
   * They must both be non-negative and the number of control points must be
   * greature than the degree for both axes.
   */
  NURBSPatch(ArrayView<const PointType, 2> controlPoints,
             ArrayView<const T, 2> weights,
             const KnotVectorType& knotVector_u,
             const KnotVectorType& knotVector_v)
    : m_knotvec_u(knotVector_u)
    , m_knotvec_v(knotVector_v)
  {
    const int knot_deg_u = m_knotvec_u.getDegree();
    const int knot_deg_v = m_knotvec_v.getDegree();
    SLIC_ASSERT(knot_deg_u >= -1 && knot_deg_v >= -1);

    if(const bool is_empty = (knot_deg_u == -1 && knot_deg_v == -1); is_empty)
    {
      SLIC_ASSERT(controlPoints.empty());
      SLIC_ASSERT(weights.empty());
    }
    else
    {
      SLIC_ASSERT(knot_deg_u >= 0 && knot_deg_v >= 0);
      AXOM_MAYBE_UNUSED const int deg_u = utilities::max(0, knot_deg_u);
      AXOM_MAYBE_UNUSED const int deg_v = utilities::max(0, knot_deg_v);
      const int npts_u = knotVector_u.getNumControlPoints();
      const int npts_v = knotVector_v.getNumControlPoints();
      SLIC_ASSERT(npts_u > deg_u && npts_v > deg_v);

      if(controlPoints.empty())
      {
        m_controlPoints.resize(npts_u, npts_v);
      }
      else
      {
        SLIC_ASSERT(controlPoints.data() != nullptr);
        SLIC_ASSERT(controlPoints.shape()[0] == npts_u);
        SLIC_ASSERT(controlPoints.shape()[1] == npts_v);
        m_controlPoints = controlPoints;
      }

      if(!weights.empty())
      {
        SLIC_ASSERT(weights.data() != nullptr);
        SLIC_ASSERT(weights.shape()[0] == npts_u);
        SLIC_ASSERT(weights.shape()[1] == npts_v);
        m_weights = weights;
      }

      SLIC_ASSERT(isValidNURBS());
    }
  }

  /*!
   * \brief Constructor for a NURBSPatch from 2D ArrayViews of control points
   *        and KnotVectors for the u- and v- directions
   *
   * \param [in] controlPoints 2D ArrayView of control points
   * \param [in] weights 2D ArrayView of weights
   * \param [in] knotVector_u, knotVector_v The knot vectors in the u- and v- directions
   */
  NURBSPatch(ArrayView<const PointType, 2> controlPoints,
             const KnotVectorType& knotVector_u,
             const KnotVectorType& knotVector_v)
    : NURBSPatch(controlPoints, axom::ArrayView<const T, 2>(nullptr, {{0, 0}}), knotVector_u, knotVector_v)
  { }

  /*!
   * \brief Constructor for a NURBSPatch from 2D ArrayViews of control points
   *        and KnotVectors for the u- and v- directions
   * \overload Overload for non-const PointType
   */
  NURBSPatch(ArrayView<PointType, 2> controlPoints,
             const KnotVectorType& knotVector_u,
             const KnotVectorType& knotVector_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(controlPoints.data(), controlPoints.shape()),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 knotVector_u,
                 knotVector_v)
  { }

  /*!
   * \brief Constructor for a NURBSPatch from 2D ArrayViews of control points and weights
   *        and KnotVectors for the u- and v- directions
   * \overload Overload for non-const PointType and weights
   */
  NURBSPatch(ArrayView<PointType, 2> controlPoints,
             ArrayView<T, 2> weights,
             const KnotVectorType& knotVector_u,
             const KnotVectorType& knotVector_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(controlPoints.data(), controlPoints.shape()),
                 axom::ArrayView<const T, 2>(weights.data(), weights.shape()),
                 knotVector_u,
                 knotVector_v)
  { }

  /*!
   * \brief Constructor for a NURBS Patch from an array of coordinates and degrees
   *
   * \param [in] pts A 1D C-style array of npts_u*npts_v control points
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] deg_u, deg_v The patch's degree on the first and second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   */
  NURBSPatch(const PointType* pts, int npts_u, int npts_v, int deg_u, int deg_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts, {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 KnotVectorType(npts_u, deg_u),
                 KnotVectorType(npts_v, deg_v))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D C-style array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] deg_u, deg_v The patch's degree on the first and second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   */
  NURBSPatch(const PointType* pts, const T* weights, int npts_u, int npts_v, int deg_u, int deg_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts, {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(weights, {{weights ? npts_u : 0, weights ? npts_v : 0}}),
                 KnotVectorType(npts_u, deg_u),
                 KnotVectorType(npts_v, deg_v))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 1D arrays of coordinates and degrees
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] deg_u, deg_v The patch's degree on the first and second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   */
  NURBSPatch(const CoordsVec& pts, int npts_u, int npts_v, int deg_u, int deg_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts.data(), {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 KnotVectorType(npts_u, deg_u),
                 KnotVectorType(npts_v, deg_v))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 1D arrays of coordinates and weights
   *
   * \param [in] pts A 1D axom::Array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D axom::Array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] deg_u, deg_v The patch's degree on the first and second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   */
  NURBSPatch(const CoordsVec& pts, const WeightsVec& weights, int npts_u, int npts_v, int deg_u, int deg_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts.data(), {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(weights.data(), {{npts_u, npts_v}}),
                 KnotVectorType(npts_u, deg_u),
                 KnotVectorType(npts_v, deg_v))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 2D arrays of coordinates and degrees
   *
   * \param [in] pts A 2D axom::Array of (npts_u, npts_v) control points
   * \param [in] deg_u, deg_v The patch's degree on the first and second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   */
  NURBSPatch(const CoordsMat& pts, int deg_u, int deg_v)
    : NURBSPatch(pts.view(),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 KnotVectorType(pts.shape()[0], deg_u),
                 KnotVectorType(pts.shape()[1], deg_v))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 2D arrays of coordinates and weights
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D axom::Array of (ord_u+1, ord_v+1) positive weights
   * \param [in] deg_u, deg_v The patch's degree on the first and second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   */
  NURBSPatch(const CoordsMat& pts, const WeightsMat& weights, int deg_u, int deg_v)
    : NURBSPatch(pts.view(),
                 weights.view(),
                 KnotVectorType(pts.shape()[0], deg_u),
                 KnotVectorType(pts.shape()[1], deg_v))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from C-style arrays of coordinates and knot vectors
   *
   * \param [in] pts A 1D C-style array of npts_u*npts_v control points
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] knots_u A 1D C-style array of npts_u + deg_u + 1 knots
   * \param [in] nkts_u The number of knots in the u direction
   * \param [in] knots_v A 1D C-style array of npts_v + deg_v + 1 knots
   * \param [in] nkts_v The number of knots in the v direction
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires valid pointers and knot vectors
   */
  NURBSPatch(const PointType* pts,
             int npts_u,
             int npts_v,
             const T* knots_u,
             int nkts_u,
             const T* knots_v,
             int nkts_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts, {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 KnotVectorType(axom::ArrayView<const T>(knots_u, nkts_u), nkts_u - npts_u - 1),
                 KnotVectorType(axom::ArrayView<const T>(knots_v, nkts_v), nkts_v - npts_v - 1))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from C-style arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of npts_u*npts_v control points
   * \param [in] weights A 1D C-style array of npts_u*npts_v positive weights
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] knots_u A 1D C-style array of npts_u + deg_u + 1 knots
   * \param [in] nkts_u The number of knots in the u direction
   * \param [in] knots_v A 1D C-style array of npts_v + deg_v + 1 knots
   * \param [in] nkts_v The number of knots in the v direction
   * 
   * For clamped and continuous patch axes, npts and the knot vector  uniquely determine the degree
   * \pre Requires valid pointers and knot vectors
   */
  NURBSPatch(const PointType* pts,
             const T* weights,
             int npts_u,
             int npts_v,
             const T* knots_u,
             int nkts_u,
             const T* knots_v,
             int nkts_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts, {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(weights, {{weights ? npts_u : 0, weights ? npts_v : 0}}),
                 KnotVectorType(axom::ArrayView<const T>(knots_u, nkts_u), nkts_u - npts_u - 1),
                 KnotVectorType(axom::ArrayView<const T>(knots_v, nkts_v), nkts_v - npts_v - 1))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates and knots
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             int npts_u,
             int npts_v,
             const axom::Array<T>& knots_u,
             const axom::Array<T>& knots_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts.data(), {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 KnotVectorType(knots_u.view(), knots_u.size() - npts_u - 1),
                 KnotVectorType(knots_v.view(), knots_v.size() - npts_v - 1))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates, weights, and knots
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] weights A 1D axom::Array of npts_u*npts_v positive weights
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             const WeightsVec& weights,
             int npts_u,
             int npts_v,
             const axom::Array<T>& knots_u,
             const axom::Array<T>& knots_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts.data(), {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(weights.data(), {{npts_u, npts_v}}),
                 KnotVectorType(knots_u.view(), knots_u.size() - npts_u - 1),
                 KnotVectorType(knots_v.view(), knots_v.size() - npts_v - 1))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates and KnotVectors
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] knotvec_u, knotvec_v  KnotVector objects for the first and second axis
   * 
   * For clamped and continuous patch axes, npts and the knot vectoruniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             int npts_u,
             int npts_v,
             const KnotVectorType& knotvec_u,
             const KnotVectorType& knotvec_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts.data(), {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 knotvec_u,
                 knotvec_v)
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates, weights, and KnotVectors
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] weights A 1D axom::Array of npts_u*npts_v positive weights
   * \param [in] npts_u, npts_v The number of control points on the first and second axis
   * \param [in] knotvec_u, knotvec_v KnotVector objects for the first and second axis
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             const WeightsVec& weights,
             int npts_u,
             int npts_v,
             const KnotVectorType& knotvec_u,
             const KnotVectorType& knotvec_v)
    : NURBSPatch(axom::ArrayView<const PointType, 2>(pts.data(), {{npts_u, npts_v}}),
                 axom::ArrayView<const T, 2>(weights.data(), {{npts_u, npts_v}}),
                 knotvec_u,
                 knotvec_v)
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 2D axom::Array array of coordinates and array of knots
   *
   * \param [in] pts A 2D axom::Array of (npts_u, npts_v) control points
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts, const axom::Array<T>& knots_u, const axom::Array<T>& knots_v)
    : NURBSPatch(pts.view(),
                 axom::ArrayView<const T, 2>(nullptr, {{0, 0}}),
                 KnotVectorType(knots_u.view(), knots_u.size() - pts.shape()[0] - 1),
                 KnotVectorType(knots_v.view(), knots_v.size() - pts.shape()[1] - 1))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 2D axom::Array array of coordinates, weights, and array of knots
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D axom::Array of (ord_u+1, ord_v+1) positive weights
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts,
             const WeightsMat& weights,
             const axom::Array<T>& knots_u,
             const axom::Array<T>& knots_v)
    : NURBSPatch(pts.view(),
                 weights.view(),
                 KnotVectorType(knots_u.view(), knots_u.size() - pts.shape()[0] - 1),
                 KnotVectorType(knots_v.view(), knots_v.size() - pts.shape()[1] - 1))
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array array of coordinates and KnotVector objects
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] knotvec_u, knotvec_v KnotVector objects for the first and second axis
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts, const KnotVectorType& knotvec_u, const KnotVectorType& knotvec_v)
    : NURBSPatch(pts.view(), axom::ArrayView<const T, 2>(nullptr, {{0, 0}}), knotvec_u, knotvec_v)
  { }

  /*!
   * \brief Constructor for a NURBS Patch from 2D axom::Array array of coordinates, weights, and KnotVector objects
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D axom::Array of (ord_u+1, ord_v+1) positive weights
   * \param [in] knotvec_u, knotvec_v KnotVector objects for the first and second axis
   * 
   * For clamped and continuous patch axes, npts and the knot vector uniquely determine the degree
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts,
             const WeightsMat& weights,
             const KnotVectorType& knotvec_u,
             const KnotVectorType& knotvec_v)
    : NURBSPatch(pts.view(), weights.view(), knotvec_u, knotvec_v)
  { }

  ///@}

  ///@{
  /// \name Query/modify patch properties (degree, rationality, ...)

  /*!
   * \brief Reset the degree and resize arrays of points (and weights)
   * 
   * \param [in] npts_u, npts_v The target number of control points on the first and second axis
   * \param [in] deg_u, deg_v The target degrees on the first and second axis
   * 
   * \warning This method will replace existing knot vectors with a uniform one.
   */
  void setParameters(int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(npts_u > deg_u && npts_v > deg_v);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);

    if(isRational())
    {
      m_weights.resize(npts_u, npts_v);
    }

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);

    makeNonrational();
  }

  /*!
   * \brief Reset the knot vector in u
   *
   * \param [in] deg The target degree
   * 
   * \warning This method does NOT change the existing control points, 
   *  i.e. does not perform degree elevation/reduction. 
   *  Will replace existing knot vector with a uniform one.
   *  
   * \pre Requires deg_u < npts_u and deg >= 0
   */
  void setDegree_u(int deg)
  {
    SLIC_ASSERT(0 <= deg && deg < getNumControlPoints_u());

    m_knotvec_u.makeUniform(getNumControlPoints_u(), deg);
  }

  /*!
   * \brief Reset the knot vector in v
   *
   * \param [in] deg The target degree
   * 
   * \warning This method does NOT change the existing control points, 
   *  i.e. does not perform degree elevation/reduction. 
   *  Will replace existing knot vector with a uniform one.
   *  
   * \pre Requires deg_v < npts_v and deg >= 0
   */
  void setDegree_v(int deg)
  {
    SLIC_ASSERT(0 <= deg && deg < getNumControlPoints_v());

    m_knotvec_v.makeUniform(getNumControlPoints_v(), deg);
  }

  /*!
   * \brief Reset the knot vector and increase the number of control points
   *
   * \param [in] deg_u The target degree in u
   * \param [in] deg_v The target degree in v
   *
   * \warning This method does NOT change the existing control points,
   *  i.e. is not performing degree elevation or reduction.
   * \pre Requires deg_u < npts_u and deg_v < npts_v
   */
  void setDegree(int deg_u, int deg_v)
  {
    setDegree_u(deg_u);
    setDegree_v(deg_v);
  }

  /// \brief Returns the degree of the NURBS Patch on the first axis
  int getDegree_u() const { return m_knotvec_u.getDegree(); }

  /// \brief Returns the degree of the NURBS Patch on the second axis
  int getDegree_v() const { return m_knotvec_v.getDegree(); }

  /// \brief Returns the order (degree + 1) of the NURBS Patch on the first axis
  int getOrder_u() const { return m_knotvec_u.getDegree() + 1; }

  /// \brief Returns the order of the NURBS Patch on the second axis
  int getOrder_v() const { return m_knotvec_v.getDegree() + 1; }

  /*!
   * \brief Set the number control points in u
   *
   * \param [in] npts The target number of control points
   * 
   * \warning This method does NOT maintain the patch shape,
   *  i.e. is not performing knot insertion/removal.
   *  Will replace existing knot vectots with uniform ones.
   */
  void setNumControlPoints(int npts_u, int npts_v)
  {
    SLIC_ASSERT(npts_u > getDegree_u());
    SLIC_ASSERT(npts_v > getDegree_v());

    m_controlPoints.resize(npts_u, npts_v);

    if(isRational())
    {
      m_weights.resize(npts_u, npts_v);
    }

    m_knotvec_u.makeUniform(npts_u, getDegree_u());
    m_knotvec_v.makeUniform(npts_v, getDegree_v());
  }

  /// \brief Returns the number of control points in the NURBS Patch on the first axis
  int getNumControlPoints_u() const { return static_cast<int>(m_controlPoints.shape()[0]); }

  /// \brief Returns the number of control points in the NURBS Patch on the second axis
  int getNumControlPoints_v() const { return static_cast<int>(m_controlPoints.shape()[1]); }

  /*!
   * \brief Set the number control points in u
   *
   * \param [in] npts The target number of control points
   * 
   * \warning This method does NOT maintain the patch shape,
   *  i.e. is not performing knot insertion/removal.
   */
  void setNumControlPoints_u(int npts)
  {
    SLIC_ASSERT(npts > getDegree_u());

    m_controlPoints.resize(npts, getNumControlPoints_v());

    if(isRational())
    {
      m_weights.resize(npts, getNumControlPoints_v());
    }

    m_knotvec_u.makeUniform(npts, getDegree_u());
  }

  /*!
   * \brief Set the number control points in v
   *
   * \param [in] npts The target number of control points
   * 
   * \warning This method does NOT maintain the patch shape,
   *  i.e. is not performing knot insertion/removal.
   */
  void setNumControlPoints_v(int npts)
  {
    SLIC_ASSERT(npts > getDegree_v());

    m_controlPoints.resize(getNumControlPoints_u(), npts);

    if(isRational())
    {
      m_weights.resize(getNumControlPoints_u(), npts);
    }

    m_knotvec_v.makeUniform(npts, getDegree_v());
  }

  /// Clears the list of control points, make nonrational
  void clear()
  {
    m_controlPoints.clear();
    m_knotvec_u.clear();
    m_knotvec_v.clear();
    m_trimmingCurves.clear();
    makeNonrational();
    makeUntrimmed();
  }

  /// \brief Use array size as flag for rationality
  bool isRational() const { return !m_weights.empty(); }

  /// \brief Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      auto patch_shape = m_controlPoints.shape();
      m_weights.resize(patch_shape[0], patch_shape[1]);
      m_weights.fill(1.0);
    }
  }

  /// \brief Make nonrational by shrinking array of weights
  void makeNonrational() { m_weights.clear(); }

  /// \brief Function to check if the NURBS surface is valid
  bool isValidNURBS() const
  {
    // Check monotonicity, open-ness, continuity of each knot vector
    if(!m_knotvec_u.isValid() || !m_knotvec_v.isValid())
    {
      return false;
    }

    // Number of knots must match the number of control points
    int deg_u = m_knotvec_u.getDegree();
    int deg_v = m_knotvec_v.getDegree();

    // Number of knots must match the number of control points
    auto patch_shape = m_controlPoints.shape();
    if(m_knotvec_u.getNumKnots() != patch_shape[0] + deg_u + 1 ||
       m_knotvec_v.getNumKnots() != patch_shape[1] + deg_v + 1)
    {
      return false;
    }

    if(isRational())
    {
      // Number of control points must match number of weights
      auto weights_shape = m_weights.shape();
      if(weights_shape[0] != patch_shape[0] || weights_shape[1] != patch_shape[1])
      {
        return false;
      }

      // Weights must be positive
      for(int i = 0; i < weights_shape[0]; ++i)
      {
        for(int j = 0; j < weights_shape[1]; ++j)
        {
          if(m_weights(i, j) <= 0.0)
          {
            return false;
          }
        }
      }
    }

    return true;
  }

  ///@}

  ///@{
  /// \name Query/modify patch's geometry (control points, weights, bounding box, ...)

  /// Retrieves the control point at index \a (idx_p, idx_q)
  PointType& operator()(int ui, int vi) { return m_controlPoints(ui, vi); }

  /// Retrieves the vector of control points at index \a idx
  const PointType& operator()(int ui, int vi) const { return m_controlPoints(ui, vi); }

  /// Returns a reference to the NURBS patch's control points
  CoordsMat& getControlPoints() { return m_controlPoints; }

  /// Returns a reference to the NURBS patch's control points
  const CoordsMat& getControlPoints() const { return m_controlPoints; }

  /*!
   * \brief Get a specific weight
   *
   * \param [in] ui The index of the weight on the first axis
   * \param [in] vi The index of the weight on the second axis
   * \pre Requires that the surface be rational
   */
  const T& getWeight(int ui, int vi) const
  {
    SLIC_ASSERT(isRational());
    return m_weights(ui, vi);
  }

  /*!
   * \brief Set the weight at a specific index
   *
   * \param [in] ui The index of the weight in on the first axis
   * \param [in] vi The index of the weight in on the second axis
   * \param [in] weight The updated value of the weight
   * \pre Requires that the surface be rational
   * \pre Requires that the weight be positive
   */
  void setWeight(int ui, int vi, T weight)
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(weight > 0);

    m_weights(ui, vi) = weight;
  }

  /// Returns a reference to the NURBS patch's weights
  WeightsMat& getWeights() { return m_weights; }

  /// Returns a const reference to the NURBS patch's weights
  const WeightsMat& getWeights() const { return m_weights; }

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

  ///@}

  ///@{
  //!  @name Methods for (untrimmed) Patch Parameterization.
  //!
  //! These methods operate on the parameterization of the patch leaving the geometry unchanged.

  /*!
   * \brief Set the knot value in the u vector at a specific index
   *
   * \param [in] idx The index of the knot
   * \param [in] knot The updated value of the knot
   */
  void setKnot_u(int idx, T knot) { m_knotvec_u[idx] = knot; }

  /*!
   * \brief Set the knot value in the v vector at a specific index
   *
   * \param [in] idx The index of the knot
   * \param [in] knot The updated value of the knot
   */
  void setKnot_v(int idx, T knot) { m_knotvec_v[idx] = knot; }

  /*! 
   * \brief Set the u knot vector by an axom::Array
   *
   * \param [in] knots The new knot vector
   */
  void setKnots_u(const axom::Array<T>& knots, int degree)
  {
    m_knotvec_u = KnotVectorType(knots, degree);
  }

  /*! 
   * \brief Set the v knot vector by an axom::Array
   *
   * \param [in] knots The new knot vector
   */
  void setKnots_v(const axom::Array<T>& knots, int degree)
  {
    m_knotvec_v = KnotVectorType(knots, degree);
  }

  /*! 
   * \brief Set the u knot vector by a KnotVector object
   *
   * \param [in] knotVector The new knot vector
   */
  void setKnots_u(const KnotVectorType& knotVector) { m_knotvec_u = knotVector; }

  /*! 
   * \brief Set the v knot vector by a KnotVector object
   *
   * \param [in] knotVector The new knot vector
   */
  void setKnots_v(const KnotVectorType& knotVector) { m_knotvec_v = knotVector; }

  /// \brief Return a reference to the KnotVector instance on the first axis
  KnotVectorType& getKnots_u() { return m_knotvec_u; }

  /// \brief Return a const reference to the KnotVector instance on the first axis
  const KnotVectorType& getKnots_u() const { return m_knotvec_u; }

  /// \brief Get the minimum knot value in the u-axis
  T getMinKnot_u() const { return m_knotvec_u.getMinKnot(); }

  /// \brief Get the maximum knot value in the u-axis
  T getMaxKnot_u() const { return m_knotvec_u.getMaxKnot(); }

  /// \brief Get the length of the parameter space bounding box
  T getParameterSpaceDiagonal() const
  {
    T u_length = getMaxKnot_u() - getMinKnot_u();
    T v_length = getMaxKnot_v() - getMinKnot_v();

    return std::sqrt(u_length * u_length + v_length * v_length);
  }

  /// \brief Return a reference to the KnotVector instance on the second axis
  KnotVectorType& getKnots_v() { return m_knotvec_v; }

  /// \brief Return a const reference to the KnotVector instance on the second axis
  const KnotVectorType& getKnots_v() const { return m_knotvec_v; }

  /// \brief Get the minimum knot value in the v-axis
  T getMinKnot_v() const { return m_knotvec_v.getMinKnot(); }

  /// \brief Get the maximum knot value in the v-axis
  T getMaxKnot_v() const { return m_knotvec_v.getMaxKnot(); }

  /// \brief Return the length of the knot vector on the first axis
  int getNumKnots_u() const { return m_knotvec_u.getNumKnots(); }

  /// \brief Return the length of the knot vector on the second axis
  int getNumKnots_v() const { return m_knotvec_v.getNumKnots(); }

  /*! 
   * \brief Insert a knot to the u knot vector to have the given multiplicity
   *
   * \param [in] u The parameter value of the knot to insert
   * \param [in] target_multiplicity The multiplicity of the knot to insert
   * \return The index of the new knot
   * 
   * Algorithm A5.3 on p. 155 of "The NURBS Book"
   * 
   * \note If the knot is already present, it will be inserted
   *  up to the given multiplicity, or the maximum permitted by the degree
   * 
   * \pre Requires \a u in the span of the knots (up to a small tolerance)
   * 
   * \note If u is outside the knot span up this tolerance, it is clamped to the span
   * 
   * \return The (maximum) index of the new knot
   */
  axom::IndexType insertKnot_u(T u, int target_multiplicity = 1)
  {
    SLIC_ASSERT_MSG(isValidParameter_u(u, 1e-5),
                    axom::fmt::format("Requested u-parameter {} for knot insertion is outside "
                                      "valid range [{},{}] with tolerance {}",
                                      u,
                                      getMinKnot_u(),
                                      getMaxKnot_u(),
                                      1e-5));

    u = axom::utilities::clampVal(u, getMinKnot_u(), getMaxKnot_u());

    SLIC_ASSERT(target_multiplicity > 0);

    const bool isRationalPatch = isRational();

    const int np = getNumControlPoints_u() - 1;
    const int p = getDegree_u();

    const int nq = getNumControlPoints_v() - 1;

    // Find the span and initial multiplicity of the knot
    int s = 0;
    const auto k = m_knotvec_u.findSpan(u, s);

    // Find how many knots we need to insert
    int r = axom::utilities::min(target_multiplicity - s, p - s);
    if(r <= 0)
    {
      return k;
    }

    // Temp variable
    axom::IndexType L;

    // Compute the alphas, which depend only on the knot vector
    axom::Array<T, 2> alpha(p - s, r + 1);
    for(int j = 1; j <= r; ++j)
    {
      L = k - p + j;
      for(int i = 0; i <= p - j - s; ++i)
      {
        alpha[i][j] = (u - m_knotvec_u[L + i]) / (m_knotvec_u[i + k + 1] - m_knotvec_u[L + i]);
      }
    }

    // Store the new control points and weights
    CoordsMat newControlPoints(np + 1 + r, nq + 1);
    WeightsMat newWeights(0, 0);
    if(isRationalPatch)
    {
      newWeights.resize(np + 1 + r, nq + 1);
    }

    // Store a temporary array of points and weights
    CoordsVec tempControlPoints(p + 1);
    WeightsVec tempWeights(isRationalPatch ? p + 1 : 0);

    // Insert the knot for each row
    for(int row = 0; row <= nq; ++row)
    {
      // Save unaltered control points
      for(int i = 0; i <= k - p; ++i)
      {
        newControlPoints(i, row) = m_controlPoints(i, row);
        if(isRationalPatch)
        {
          newWeights(i, row) = m_weights(i, row);
        }
      }

      for(int i = k - s; i <= np; ++i)
      {
        newControlPoints(i + r, row) = m_controlPoints(i, row);
        if(isRationalPatch)
        {
          newWeights(i + r, row) = m_weights(i, row);
        }
      }

      // Load auxiliary control points
      for(int i = 0; i <= p - s; ++i)
      {
        for(int n = 0; n < NDIMS; ++n)
        {
          tempControlPoints[i][n] =
            m_controlPoints(k - p + i, row)[n] * (isRationalPatch ? m_weights(k - p + i, row) : 1.0);
        }

        if(isRationalPatch)
        {
          tempWeights[i] = m_weights(k - p + i, row);
        }
      }

      // Insert the knot r times
      for(int j = 1; j <= r; ++j)
      {
        L = k - p + j;
        for(int i = 0; i <= p - j - s; ++i)
        {
          tempControlPoints[i].array() = alpha(i, j) * tempControlPoints[i + 1].array() +
            (1.0 - alpha(i, j)) * tempControlPoints[i].array();

          if(isRationalPatch)
          {
            tempWeights[i] = alpha(i, j) * tempWeights[i + 1] + (1.0 - alpha(i, j)) * tempWeights[i];
          }
        }

        for(int n = 0; n < NDIMS; ++n)
        {
          newControlPoints(L, row)[n] =
            tempControlPoints[0][n] / (isRationalPatch ? tempWeights[0] : 1.0);
          newControlPoints(k + r - j - s, row)[n] =
            tempControlPoints[p - j - s][n] / (isRationalPatch ? tempWeights[p - j - s] : 1.0);
        }

        if(isRationalPatch)
        {
          newWeights(L, row) = tempWeights[0];
          newWeights(k + r - j - s, row) = tempWeights[p - j - s];
        }
      }

      // Load the remaining control points
      for(int i = L + 1; i < k - s; ++i)
      {
        for(int n = 0; n < NDIMS; ++n)
        {
          newControlPoints(i, row)[n] =
            tempControlPoints[i - L][n] / (isRationalPatch ? tempWeights[i - L] : 1.0);
        }

        if(isRationalPatch)
        {
          newWeights(i, row) = tempWeights[i - L];
        }
      }
    }

    // Update the knot vector and control points
    m_knotvec_u.insertKnotBySpan(k, u, r);
    m_controlPoints = newControlPoints;
    m_weights = newWeights;

    return k + r;
  }

  /*! 
   * \brief Insert a knot to the v knot vector to have the given multiplicity
   *
   * \param [in] v The parameter value of the knot to insert
   * \param [in] target_multiplicity The multiplicity of the knot to insert
   * \return The index of the new knot
   * 
   * Algorithm A5.3 on p. 155 of "The NURBS Book"
   * 
   * \note If the knot is already present, it will be inserted
   *  up to the given multiplicity, or the maximum permitted by the degree
   * 
   * \pre Requires \a v in the span of the knots (up to a small tolerance)
   * 
   * \note If v is outside the knot span up this tolerance, it is clamped to the span
   * 
   * \return The (maximum) index of the new knot
   */
  axom::IndexType insertKnot_v(T v, int target_multiplicity = 1)
  {
    SLIC_ASSERT_MSG(isValidParameter_v(v, 1e-5),
                    axom::fmt::format("Requested v-parameter {} for knot insertion is outside "
                                      "valid range [{},{}] with tolerance {}",
                                      v,
                                      getMinKnot_v(),
                                      getMaxKnot_v(),
                                      1e-5));

    v = axom::utilities::clampVal(v, getMinKnot_v(), getMaxKnot_v());

    SLIC_ASSERT(target_multiplicity > 0);

    const bool isRationalPatch = isRational();

    const int np = getNumControlPoints_u() - 1;

    const int nq = getNumControlPoints_v() - 1;
    const int q = getDegree_v();

    // Find the span and initial multiplicity of the knot
    int s = 0;
    const auto k = m_knotvec_v.findSpan(v, s);

    // Find how many knots we need to insert
    int r = axom::utilities::min(target_multiplicity - s, q - s);
    if(r <= 0)
    {
      return k;
    }

    // Temp variable
    axom::IndexType L;

    // Compute the alphas, which depend only on the knot vector
    axom::Array<T, 2> alpha(q - s, r + 1);
    for(int j = 1; j <= r; ++j)
    {
      L = k - q + j;
      for(int i = 0; i <= q - j - s; ++i)
      {
        alpha[i][j] = (v - m_knotvec_v[L + i]) / (m_knotvec_v[i + k + 1] - m_knotvec_v[L + i]);
      }
    }

    // Store the new control points and weights
    CoordsMat newControlPoints(np + 1, nq + 1 + r);
    WeightsMat newWeights(0, 0);
    if(isRationalPatch)
    {
      newWeights.resize(np + 1, nq + 1 + r);
    }

    // Store a temporary array of points and weights
    CoordsVec tempControlPoints(q + 1);
    WeightsVec tempWeights(isRationalPatch ? q + 1 : 0);

    // Insert the knot for each row
    for(int col = 0; col <= np; ++col)
    {
      // Save unaltered control points
      for(int i = 0; i <= k - q; ++i)
      {
        newControlPoints(col, i) = m_controlPoints(col, i);
        if(isRationalPatch)
        {
          newWeights(col, i) = m_weights(col, i);
        }
      }

      for(int i = k - s; i <= nq; ++i)
      {
        newControlPoints(col, i + r) = m_controlPoints(col, i);
        if(isRationalPatch)
        {
          newWeights(col, i + r) = m_weights(col, i);
        }
      }

      // Load auxiliary control points
      for(int i = 0; i <= q - s; ++i)
      {
        for(int n = 0; n < NDIMS; ++n)
        {
          tempControlPoints[i][n] =
            m_controlPoints(col, k - q + i)[n] * (isRationalPatch ? m_weights(col, k - q + i) : 1.0);
        }

        if(isRationalPatch)
        {
          tempWeights[i] = m_weights(col, k - q + i);
        }
      }

      // Insert the knot r times
      for(int j = 1; j <= r; ++j)
      {
        L = k - q + j;
        for(int i = 0; i <= q - j - s; ++i)
        {
          tempControlPoints[i].array() = alpha(i, j) * tempControlPoints[i + 1].array() +
            (1.0 - alpha(i, j)) * tempControlPoints[i].array();

          if(isRationalPatch)
          {
            tempWeights[i] = alpha(i, j) * tempWeights[i + 1] + (1.0 - alpha(i, j)) * tempWeights[i];
          }
        }

        for(int n = 0; n < NDIMS; ++n)
        {
          newControlPoints(col, L)[n] =
            tempControlPoints[0][n] / (isRationalPatch ? tempWeights[0] : 1.0);
          newControlPoints(col, k + r - j - s)[n] =
            tempControlPoints[q - j - s][n] / (isRationalPatch ? tempWeights[q - j - s] : 1.0);
        }

        if(isRationalPatch)
        {
          newWeights(col, L) = tempWeights[0];
          newWeights(col, k + r - j - s) = tempWeights[q - j - s];
        }
      }

      // Load the remaining control points
      for(int i = L + 1; i < k - s; ++i)
      {
        for(int n = 0; n < NDIMS; ++n)
        {
          newControlPoints(col, i)[n] =
            tempControlPoints[i - L][n] / (isRationalPatch ? tempWeights[i - L] : 1.0);
        }

        if(isRationalPatch)
        {
          newWeights(col, i) = tempWeights[i - L];
        }
      }
    }

    // Update the knot vector and control points
    m_knotvec_v.insertKnotBySpan(k, v, r);
    m_controlPoints = newControlPoints;
    m_weights = newWeights;

    return k + r;
  }

  /// \brief Reverses all trimming curves in the patch
  /// \warning Trimming curves should be oriented CCW in parameter space,
  ///   and this method may make them CW
  void reverseTrimmingCurves()
  {
    for(auto& curve : m_trimmingCurves)
    {
      curve.reverseOrientation();
    }
  }

  /*!
   * \brief Reverses the order of one direction of the NURBS patch's control points and weights
   *
   * This method does not affect the position of the patch in space, or its 
   *  trimming curves, but it does reverse the patch's normal vectors. 
   * 
   * \param [in] axis orientation of patch. 0 to reverse in u, 1 for reverse in v
   */
  void reverseOrientation(int axis)
  {
    if(axis == 0)
    {
      reverseOrientation_u();
    }
    else
    {
      reverseOrientation_v();
    }
  }

  /// \brief Reverses the order of the control points, weights, and knots on the first axis
  void reverseOrientation_u()
  {
    auto patch_shape = m_controlPoints.shape();
    const int npts_u_mid = (patch_shape[0] + 1) / 2;

    for(int q = 0; q < patch_shape[1]; ++q)
    {
      for(int i = 0; i < npts_u_mid; ++i)
      {
        axom::utilities::swap(m_controlPoints(i, q), m_controlPoints(patch_shape[0] - i - 1, q));
      }

      if(isRational())
      {
        for(int i = 0; i < npts_u_mid; ++i)
        {
          axom::utilities::swap(m_weights(i, q), m_weights(patch_shape[0] - i - 1, q));
        }
      }
    }

    m_knotvec_u.reverse();

    // Mirror the trimming curves on the u-axis
    auto min_u = m_knotvec_u[0];
    auto max_u = m_knotvec_u[m_knotvec_u.getNumKnots() - 1];

    for(auto& curve : m_trimmingCurves)
    {
      for(int i = 0; i < curve.getNumControlPoints(); ++i)
      {
        curve[i][0] = min_u + max_u - curve[i][0];
      }

      curve.reverseOrientation();
    }
  }

  /// \brief Reverses the order of the control points, weights, and knots on the second axis
  void reverseOrientation_v()
  {
    auto patch_shape = m_controlPoints.shape();
    const int npts_v_mid = (patch_shape[1] + 1) / 2;

    for(int p = 0; p < patch_shape[0]; ++p)
    {
      for(int i = 0; i < npts_v_mid; ++i)
      {
        axom::utilities::swap(m_controlPoints(p, i), m_controlPoints(p, patch_shape[1] - i - 1));
      }

      if(isRational())
      {
        for(int i = 0; i < npts_v_mid; ++i)
        {
          axom::utilities::swap(m_weights(p, i), m_weights(p, patch_shape[1] - i - 1));
        }
      }
    }

    m_knotvec_v.reverse();

    // Mirror the trimming curves on the v-axis
    auto min_v = m_knotvec_v[0];
    auto max_v = m_knotvec_v[m_knotvec_v.getNumKnots() - 1];

    for(auto& curve : m_trimmingCurves)
    {
      for(int i = 0; i < curve.getNumControlPoints(); ++i)
      {
        curve[i][1] = min_v + max_v - curve[i][1];
      }

      curve.reverseOrientation();
    }
  }

  /*!
   * \brief Swap the axes such that s(u, v) becomes s(v, u)
   *
   * This method does not affect the position of the patch in space,
   *  or its trimming curves.
   */
  void swapAxes()
  {
    auto patch_shape = m_controlPoints.shape();

    CoordsMat new_controlPoints(patch_shape[1], patch_shape[0]);

    for(int p = 0; p < patch_shape[0]; ++p)
    {
      for(int q = 0; q < patch_shape[1]; ++q)
      {
        new_controlPoints(q, p) = m_controlPoints(p, q);
      }
    }

    m_controlPoints = new_controlPoints;

    if(isRational())
    {
      WeightsMat new_weights(patch_shape[1], patch_shape[0]);

      for(int p = 0; p < patch_shape[0]; ++p)
      {
        for(int q = 0; q < patch_shape[1]; ++q)
        {
          new_weights(q, p) = m_weights(p, q);
        }
      }

      m_weights = new_weights;
    }

    axom::utilities::swap(m_knotvec_u, m_knotvec_v);

    for(auto& curve : m_trimmingCurves)
    {
      for(int j = 0; j < curve.getNumControlPoints(); ++j)
      {
        axom::utilities::swap(curve[j][0], curve[j][1]);
      }

      curve.reverseOrientation();
    }
  }

  /// \brief Normalize the knot vectors to the span [0, 1]
  void normalize()
  {
    rescaleTrimmingCurves_u(getMinKnot_u(), getMaxKnot_u(), 0.0, 1.0);
    rescaleTrimmingCurves_v(getMinKnot_v(), getMaxKnot_v(), 0.0, 1.0);

    m_knotvec_u.normalize();
    m_knotvec_v.normalize();
  }

  /// \brief Normalize the knot vector in u to the span [0, 1]
  void normalize_u()
  {
    rescaleTrimmingCurves_u(getMinKnot_u(), getMaxKnot_u(), 0.0, 1.0);

    m_knotvec_u.normalize();
  }

  /// \brief Normalize the knot vector in v to the span [0, 1]
  void normalize_v()
  {
    rescaleTrimmingCurves_v(getMinKnot_v(), getMaxKnot_v(), 0.0, 1.0);

    m_knotvec_v.normalize();
  }

  /// \brief Normalize to the span [0, N] x [0, M] where N and M are the number of spans in u and v
  void normalizeBySpan()
  {
    auto n = m_knotvec_u.getNumKnotSpans();
    auto m = m_knotvec_v.getNumKnotSpans();

    rescaleTrimmingCurves_u(getMinKnot_u(), getMaxKnot_u(), 0.0, n);
    rescaleTrimmingCurves_v(getMinKnot_v(), getMaxKnot_v(), 0.0, m);

    m_knotvec_u.rescale(0, n);
    m_knotvec_v.rescale(0, m);
  }

  /*!
   * \brief Rescale both knot vectors to the span of [a, b]
   * 
   * \param [in] a The lower bound of the new knot vector
   * \param [in] b The upper bound of the new knot vector
   * 
   * \pre Requires a < b
   */
  void rescale(T a, T b)
  {
    SLIC_ASSERT(a < b);

    rescaleTrimmingCurves_u(getMinKnot_u(), getMaxKnot_u(), a, b);
    rescaleTrimmingCurves_v(getMinKnot_v(), getMaxKnot_v(), a, b);

    m_knotvec_u.rescale(a, b);
    m_knotvec_v.rescale(a, b);
  }

  /*!
   * \brief Rescale the knot vector in u to the span of [a, b]
   * 
   * \param [in] a The lower bound of the new knot vector
   * \param [in] b The upper bound of the new knot vector
   * 
   * \pre Requires a < b
   */
  void rescale_u(T a, T b)
  {
    SLIC_ASSERT(a < b);

    rescaleTrimmingCurves_u(getMinKnot_u(), getMaxKnot_u(), a, b);

    m_knotvec_u.rescale(a, b);
  }

  /*!
   * \brief Rescale the knot vector in v to the span of [a, b]
   * 
   * \param [in] a The lower bound of the new knot vector
   * \param [in] b The upper bound of the new knot vector
   * 
   * \pre Requires a < b
   */
  void rescale_v(T a, T b)
  {
    SLIC_ASSERT(a < b);

    rescaleTrimmingCurves_v(getMinKnot_v(), getMaxKnot_v(), a, b);

    m_knotvec_v.rescale(a, b);
  }

  /// \brief Function to check if the u parameter is within the knot span
  bool isValidParameter_u(T u, T EPS = 1e-8) const
  {
    return u >= m_knotvec_u[0] - EPS && u <= m_knotvec_u[m_knotvec_u.getNumKnots() - 1] + EPS;
  }

  /// \brief Function to check if the v parameter is within the knot span
  bool isValidParameter_v(T v, T EPS = 1e-8) const
  {
    return v >= m_knotvec_v[0] - EPS && v <= m_knotvec_v[m_knotvec_v.getNumKnots() - 1] + EPS;
  }

  /// \brief Checks if given u parameter is *interior* to the knot span
  bool isValidInteriorParameter_u(T t) const { return m_knotvec_u.isValidInteriorParameter(t); }

  /// \brief Checks if given v parameter is *interior* to the knot span
  bool isValidInteriorParameter_v(T t) const { return m_knotvec_v.isValidInteriorParameter(t); }

  /*!
   * \brief Scale the parameter space of the NURBS patch geometry 
   *         linearly (by tangents) in all directions
   *
   * \param [in] scaleFactor The multiplicative factor to expand each knot vector by
   * \param [in] removeTrimmingCurves If true, the resulting patch has no trimming curves
   *
   * Algorithm from Wolters, Hans J., "Extensions: Extrapolation Methods for CAD", 1999
   * 
   * \note This function only affects the geometry of the untrimmed NURBS patch, 
   *       and does not affect any existing trimming curves (unless explicitly removed)
   * 
   * \post If removeTrimmingCurves is false, the resulting patch will be trimmed.
   * 
   * \warning Method becomes numerically unstable for large values of scaleFactor,
   *           or for rational patches with a large range of weights.
   */
  void scaleParameterSpace(double scaleFactor, bool removeTrimmingCurves = false)
  {
    SLIC_ASSERT(scaleFactor >= 1.0);
    SLIC_WARNING_IF(scaleFactor > 1.15,
                    "Expanding patch parameter space is numerically unstable "
                    "for large values of scaleFactor.");

    double expansionAmount_u = (getMaxKnot_u() - getMinKnot_u()) * (scaleFactor - 1.0);
    double expansionAmount_v = (getMaxKnot_v() - getMinKnot_v()) * (scaleFactor - 1.0);

    expandParameterSpace(expansionAmount_u, expansionAmount_v, removeTrimmingCurves);
  }

  /*!
   * \brief Expand the parameter space of the NURBS patch geometry 
   *         linearly (by tangents) in all directions by a fixed amount
   *
   * \param [in] expansionAmount_u The absolute additive amount by which the u is expanded
   * \param [in] expansionAmount_v The absolute additive amount by which the v is expanded
   * \param [in] removeTrimmingCurves If true, the resulting patch has no trimming curves
   *
   * Algorithm from Wolters, Hans J., "Extensions: Extrapolation Methods for CAD", 1999
   * 
   * \note This function only affects the geometry of the untrimmed NURBS patch, 
   *       and does not affect any existing trimming curves (unless explicitly removed)
   * 
   * \post If removeTrimmingCurves is false, the resulting patch will be trimmed.
   * 
   * \warning Method becomes numerically unstable for values of expansionAmount that are
   *           a large fraction of the existing parameter size, or for rational patches 
   *           with a large range of weights.
   */
  void expandParameterSpace(double expansionAmount_u,
                            double expansionAmount_v,
                            bool removeTrimmingCurves = false)
  {
    SLIC_ASSERT(expansionAmount_u > 0.0 && expansionAmount_v > 0.0);
    SLIC_WARNING_IF(expansionAmount_u > 1.15 * (getMaxKnot_u() - getMinKnot_u()) ||
                      expansionAmount_v > 1.15 * (getMaxKnot_v() - getMinKnot_v()),
                    "Expanding patch parameter space is numerically unstable "
                    "for values of expansionAmount that are a large fraction of "
                    "the patch's length in parameter space.");

    if(removeTrimmingCurves)
    {
      m_trimmingCurves.clear();
    }
    else if(!isTrimmed())
    {
      // If the patch is untrimmed, we need to create new trimming curves
      //  to match the original parameter space
      makeTriviallyTrimmed();
    }

    auto n = getNumControlPoints_u();
    auto m = getNumControlPoints_v();

    if(n <= 1 || m <= 1)
    {
      return;
    }

    // When the patch is expanded in homogeneous space,
    //  weights may become negative.
    // We gurantee no negative weights by restricting this expansion
    //  to w > min_weight in homogeous space
    double min_weight = 0.0;
    if(isRational())
    {
      min_weight = m_weights(0, 0);

      for(int i = 0; i < n; ++i)
      {
        for(int j = 0; j < m; ++j)
        {
          min_weight = std::min(min_weight, m_weights(i, j));
        }
      }

      min_weight *= 0.5;
    }

    auto deg_u = getDegree_u();
    auto deg_v = getDegree_v();

    CoordsMat newControlPoints(n + 2 * deg_u, m + 2 * deg_v);
    WeightsMat newWeights(0, 0);
    if(isRational())
    {
      newWeights.resize(n + 2 * deg_u, m + 2 * deg_v);
      newWeights.fill(1.0);
    }

    axom::Array<T> newKnotVec_u, newKnotVec_v;

    // Copy the original control points
    for(int i = 0; i < n; ++i)
    {
      for(int j = 0; j < m; ++j)
      {
        newControlPoints(i + deg_u, j + deg_v) = m_controlPoints(i, j);
        if(isRational())
        {
          newWeights(i + deg_u, j + deg_v) = m_weights(i, j);
        }
      }
    }

    int nkts_v = m_knotvec_v.getNumKnots();
    int nkts_u = m_knotvec_u.getNumKnots();

    // Add the control points on the v direction
    for(int i = 0; i < n; ++i)
    {
      if(!isRational())
      {
        Vector<T, 3> v(m_controlPoints(i, 1), m_controlPoints(i, 0));
        double alpha = deg_v * expansionAmount_v / (m_knotvec_v[0] - m_knotvec_v[deg_v + 1]);

        for(int j = 0; j < deg_v; ++j)
        {
          newControlPoints(i + deg_u, j).array() =
            m_controlPoints(i, 0).array() + static_cast<T>(j - deg_v) / (deg_v)*alpha * v.array();
        }

        v = Vector<T, 3>(m_controlPoints(i, m - 2), m_controlPoints(i, m - 1));
        alpha =
          deg_v * expansionAmount_v / (m_knotvec_v[nkts_v - 1] - m_knotvec_v[nkts_v - deg_v - 2]);

        for(int j = 0; j < deg_v; ++j)
        {
          newControlPoints(i + deg_u, m + deg_v + j).array() =
            m_controlPoints(i, m - 1).array() + static_cast<T>(j + 1) / (deg_v)*alpha * v.array();
        }
      }
      else
      {
        Vector<T, 3> v(Point<T, 3>(m_controlPoints(i, 1).array() * m_weights(i, 1)),
                       Point<T, 3>(m_controlPoints(i, 0).array() * m_weights(i, 0)));
        double d_weight = m_weights(i, 0) - m_weights(i, 1);
        double alpha = deg_v * expansionAmount_v / (m_knotvec_v[0] - m_knotvec_v[deg_v + 1]);

        // New weights can't be less than min_weight
        if(d_weight != 0 && (m_weights(i, 0) - alpha * d_weight < min_weight))
        {
          alpha = (m_weights(i, 0) - min_weight) / d_weight;
        }

        for(int j = 0; j < deg_v; ++j)
        {
          newWeights(i + deg_u, j) =
            m_weights(i, 0) + static_cast<T>(j - deg_v) / (deg_v)*alpha * d_weight;

          newControlPoints(i + deg_u, j).array() =
            (m_controlPoints(i, 0).array() * m_weights(i, 0) +
             static_cast<T>(j - deg_v) / (deg_v)*alpha * v.array()) /
            newWeights(i + deg_u, j);
        }

        v = Vector<T, 3>(Point<T, 3>(m_controlPoints(i, m - 2).array() * m_weights(i, m - 2)),
                         Point<T, 3>(m_controlPoints(i, m - 1).array() * m_weights(i, m - 1)));
        d_weight = m_weights(i, m - 1) - m_weights(i, m - 2);
        alpha =
          deg_v * expansionAmount_v / (m_knotvec_v[nkts_v - 1] - m_knotvec_v[nkts_v - deg_v - 2]);

        // New weights can't be less than min_weight
        if(d_weight != 0 && (m_weights(i, m - 1) + alpha * d_weight < min_weight))
        {
          alpha = (min_weight - m_weights(i, m - 1)) / d_weight;
        }

        for(int j = 0; j < deg_v; ++j)
        {
          newWeights(i + deg_u, m + deg_v + j) =
            m_weights(i, m - 1) + static_cast<T>(j + 1) / (deg_v)*alpha * d_weight;

          newControlPoints(i + deg_u, m + deg_v + j).array() =
            (m_controlPoints(i, m - 1).array() * m_weights(i, m - 1) +
             static_cast<T>(j + 1) / (deg_v)*alpha * v.array()) /
            newWeights(i + deg_u, m + deg_v + j);
        }
      }
    }

    // Add the control points on the u direction
    //  Note that this method uses values added in the v direction,
    //  making it slightly anisotropic at the corners
    for(int j = 0; j < m + 2 * deg_v; ++j)
    {
      if(!isRational())
      {
        Vector<T, 3> v(newControlPoints(deg_u + 1, j), newControlPoints(deg_u, j));
        double alpha = deg_u * expansionAmount_v / (m_knotvec_u[0] - m_knotvec_u[deg_u + 1]);
        for(int i = 0; i < deg_u; ++i)
        {
          newControlPoints(i, j).array() = newControlPoints(deg_u, j).array() +
            static_cast<T>(i - deg_u) / (deg_u)*alpha * v.array();
        }

        v = Vector<T, 3>(newControlPoints(n - 2 + deg_u, j), newControlPoints(n + deg_u - 1, j));
        alpha =
          deg_u * expansionAmount_v / (m_knotvec_u[nkts_u - 1] - m_knotvec_u[nkts_u - deg_u - 2]);
        for(int i = 0; i < deg_u; ++i)
        {
          newControlPoints(n + deg_u + i, j).array() = newControlPoints(n + deg_u - 1, j).array() +
            static_cast<T>(i + 1) / (deg_u)*alpha * v.array();
        }
      }
      else
      {
        Vector<T, 3> v(Point<T, 3>(newControlPoints(deg_u + 1, j).array() * newWeights(deg_u + 1, j)),
                       Point<T, 3>(newControlPoints(deg_u, j).array() * newWeights(deg_u, j)));
        double d_weight = newWeights(deg_u, j) - newWeights(deg_u + 1, j);
        double alpha = deg_u * expansionAmount_v / (m_knotvec_u[0] - m_knotvec_u[deg_u + 1]);

        // New weights can't be less than min_weight
        if(d_weight != 0 && (newWeights(deg_u, j) - alpha * d_weight < min_weight))
        {
          alpha = (newWeights(deg_u, j) - min_weight) / d_weight;
        }

        for(int i = 0; i < deg_u; ++i)
        {
          newWeights(i, j) =
            newWeights(deg_u, j) + static_cast<T>(i - deg_u) / (deg_u)*alpha * d_weight;

          newControlPoints(i, j).array() =
            (newControlPoints(deg_u, j).array() * newWeights(deg_u, j) +
             static_cast<T>(i - deg_u) / (deg_u)*alpha * v.array()) /
            newWeights(i, j);
        }

        v = Vector<T, 3>(
          Point<T, 3>(newControlPoints(n - 2 + deg_u, j).array() * newWeights(n - 2 + deg_u, j)),
          Point<T, 3>(newControlPoints(n + deg_u - 1, j).array() * newWeights(n + deg_u - 1, j)));
        d_weight = newWeights(n + deg_u - 1, j) - newWeights(n - 2 + deg_u, j);
        alpha =
          deg_u * expansionAmount_v / (m_knotvec_u[nkts_u - 1] - m_knotvec_u[nkts_u - deg_u - 2]);

        // New weights can't be less than min_weight
        if(d_weight != 0 && (newWeights(n + deg_u - 1, j) + alpha * d_weight < min_weight))
        {
          alpha = (min_weight - newWeights(n + deg_u - 1, j)) / d_weight;
        }

        for(int i = 0; i < deg_u; ++i)
        {
          newWeights(n + deg_u + i, j) =
            newWeights(n + deg_u - 1, j) + static_cast<T>(i + 1) / (deg_u)*d_weight * alpha;

          newControlPoints(n + deg_u + i, j).array() =
            (newControlPoints(n + deg_u - 1, j).array() * newWeights(n + deg_u - 1, j) +
             static_cast<T>(i + 1) / (deg_u)*alpha * v.array()) /
            newWeights(n + deg_u + i, j);
        }
      }
    }

    // Fix the u knot vector
    newKnotVec_u.resize(m_knotvec_u.getNumKnots() + 2 * m_knotvec_u.getDegree());
    for(int i = 0; i <= deg_u; ++i)
    {
      newKnotVec_u[i] = m_knotvec_u[0] - expansionAmount_u;
    }
    for(int i = 0; i < m_knotvec_u.getNumKnots() - 2; ++i)
    {
      newKnotVec_u[i + deg_u + 1] = m_knotvec_u[i + 1];
    }
    for(int i = 0; i <= deg_u; ++i)
    {
      newKnotVec_u[i + deg_u + m_knotvec_u.getNumKnots() - 1] =
        m_knotvec_u[m_knotvec_u.getNumKnots() - 1] + expansionAmount_u;
    }

    // Fix the v knot vector
    newKnotVec_v.resize(m_knotvec_v.getNumKnots() + 2 * m_knotvec_v.getDegree());
    for(int i = 0; i <= deg_v; ++i)
    {
      newKnotVec_v[i] = m_knotvec_v[0] - expansionAmount_v;
    }
    for(int i = 0; i < m_knotvec_v.getNumKnots() - 2; ++i)
    {
      newKnotVec_v[i + deg_v + 1] = m_knotvec_v[i + 1];
    }
    for(int i = 0; i <= deg_v; ++i)
    {
      newKnotVec_v[i + deg_v + m_knotvec_v.getNumKnots() - 1] =
        m_knotvec_v[m_knotvec_v.getNumKnots() - 1] + expansionAmount_v;
    }

    m_controlPoints = newControlPoints;
    m_weights = newWeights;

    m_knotvec_u = KnotVectorType(newKnotVec_u, deg_u);
    m_knotvec_v = KnotVectorType(newKnotVec_v, deg_v);
  }
  ///@}

  ///@{
  //! \name Functions to evaluate (untrimmed) patch and its derivatives and normals along points or lines
  //!
  //! These methods operate only on the geometry of the patch
  //!  as defined by the control points and weights.
  //! They do not modify nor interact with the trimming curves.

  /*!
   * \brief Evaluate the NURBS patch geometry at a particular parameter value \a t
   *
   * \param [in] u The parameter value on the first axis
   * \param [in] v The parameter value on the second axis
   * 
   * Adapted from Algorithm A3.5 on page 103 of "The NURBS Book"
   * 
   * \pre Requires \a u, v in the span of each knot vector (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  PointType evaluate(T u, T v) const
  {
    u = axom::utilities::clampVal(u, getMinKnot_u(), getMaxKnot_u());
    v = axom::utilities::clampVal(v, getMinKnot_v(), getMaxKnot_v());

    const auto span_u = m_knotvec_u.findSpan(u);
    const auto span_v = m_knotvec_v.findSpan(v);

    const auto basis_funs_u = m_knotvec_u.calculateBasisFunctionsBySpan(span_u, u);
    const auto basis_funs_v = m_knotvec_v.calculateBasisFunctionsBySpan(span_v, v);

    const auto deg_u = getDegree_u();
    const auto deg_v = getDegree_v();

    int ind_u = span_u - deg_u;

    PointType S = PointType::zero();

    if(isRational())
    {
      // Evaluate the homogeneous point
      Point<T, NDIMS + 1> Sw = Point<T, NDIMS + 1>::zero();
      for(int l = 0; l <= deg_v; ++l)
      {
        Point<T, NDIMS + 1> temp = Point<T, NDIMS + 1>::zero();
        int ind_v = span_v - deg_v + l;
        for(int k = 0; k <= deg_u; ++k)
        {
          auto& the_weight = m_weights(ind_u + k, ind_v);
          auto& the_pt = m_controlPoints(ind_u + k, ind_v);

          for(int i = 0; i < NDIMS; ++i)
          {
            temp[i] += basis_funs_u[k] * the_weight * the_pt[i];
          }
          temp[NDIMS] += basis_funs_u[k] * the_weight;
        }

        for(int i = 0; i < NDIMS; ++i)
        {
          Sw[i] += basis_funs_v[l] * temp[i];
        }
        Sw[NDIMS] += basis_funs_v[l] * temp[NDIMS];
      }

      // Project the point back to coordinate space
      for(int i = 0; i < NDIMS; ++i)
      {
        S[i] = Sw[i] / Sw[NDIMS];
      }
    }
    else
    {
      for(int l = 0; l <= deg_v; ++l)
      {
        PointType temp = PointType::zero();
        int ind_v = span_v - deg_v + l;
        for(int k = 0; k <= deg_u; ++k)
        {
          for(int i = 0; i < NDIMS; ++i)
          {
            temp[i] += basis_funs_u[k] * m_controlPoints(ind_u + k, ind_v)[i];
          }
        }
        for(int i = 0; i < NDIMS; ++i)
        {
          S[i] += basis_funs_v[l] * temp[i];
        }
      }
    }

    return S;
  }

  /*!
   * \brief Returns a NURBS patch isocurve for a fixed parameter value of \a u or \a v
   *
   * \param [in] uv parameter value at which to construct the isocurve
   * \param [in] axis orientation of curve. 0 for fixed u, 1 for fixed v
   * \return c The isocurve C(v) = S(u, v) for fixed u or C(u) = S(u, v) for fixed v
   *
   * \pre Requires \a uv be in the span of the relevant knot vector
   */
  NURBSCurveType isocurve(T uv, int axis) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));

    if(axis == 0)
    {
      return isocurve_u(uv);
    }
    else
    {
      return isocurve_v(uv);
    }
  }

  /*!
   * \brief Returns a NURBS patch isocurve with a fixed value of u
   *
   * \param [in] u Parameter value fixed in the isocurve
   * \return c The isocurve C(v) = S(u, v) for fixed u
   * 
   * \pre Requires \a u be in the span of the knot vector (up to a small tolerance)
   * 
   * \note If u is outside the knot span up this tolerance, it is clamped to the span
   */
  NURBSCurveType isocurve_u(T u) const
  {
    SLIC_ASSERT_MSG(isValidParameter_u(u, 1e-5),
                    axom::fmt::format("Requested u-parameter {} for isocurve evaluation is "
                                      "outside valid range [{},{}] with tolerance {}",
                                      u,
                                      getMinKnot_u(),
                                      getMaxKnot_u(),
                                      1e-5));

    u = axom::utilities::clampVal(u, m_knotvec_u[0], m_knotvec_u[m_knotvec_u.getNumKnots() - 1]);

    using axom::utilities::lerp;

    bool isRationalPatch = isRational();

    auto patch_shape = m_controlPoints.shape();
    const int deg_u = m_knotvec_u.getDegree();
    const int deg_v = m_knotvec_v.getDegree();

    NURBSCurveType c(patch_shape[1], deg_v);
    if(isRationalPatch)
    {
      c.makeRational();
    }

    // Find the control points by evaluating each column of the patch
    const auto span_u = m_knotvec_u.findSpan(u);
    const auto N_evals_u = m_knotvec_u.calculateBasisFunctionsBySpan(span_u, u);
    for(int q = 0; q < patch_shape[1]; ++q)
    {
      Point<T, NDIMS + 1> H;
      for(int j = 0; j <= deg_u; ++j)
      {
        const auto offset = span_u - deg_u + j;
        const T weight = isRationalPatch ? m_weights(offset, q) : 1.0;
        const auto& controlPoint = m_controlPoints(offset, q);

        for(int i = 0; i < NDIMS; ++i)
        {
          H[i] += N_evals_u[j] * weight * controlPoint[i];
        }
        H[NDIMS] += N_evals_u[j] * weight;
      }

      for(int i = 0; i < NDIMS; ++i)
      {
        c[q][i] = H[i] / H[NDIMS];
      }

      if(isRationalPatch)
      {
        c.setWeight(q, H[NDIMS]);
      }
    }

    c.setKnots(m_knotvec_v);

    return c;
  }

  /*!
   * \brief Returns a NURBS patch isocurve with a fixed value of v
   *
   * \param [in] v Parameter value fixed in the isocurve
   * \return c The isocurve C(u) = S(u, v) for fixed v
   * 
   * \pre Requires \a v be in the span of the knot vector (up to a small tolerance)
   * 
   * \note If v is outside the knot span up this tolerance, it is clamped to the span
   */
  NURBSCurveType isocurve_v(T v) const
  {
    SLIC_ASSERT_MSG(isValidParameter_v(v, 1e-5),
                    axom::fmt::format("Requested v-parameter {} for isocurve evaluation is "
                                      "outside valid range [{},{}] with tolerance {}",
                                      v,
                                      getMinKnot_v(),
                                      getMaxKnot_v(),
                                      1e-5));

    v = axom::utilities::clampVal(v, m_knotvec_v[0], m_knotvec_v[m_knotvec_v.getNumKnots() - 1]);

    using axom::utilities::lerp;

    bool isRationalPatch = isRational();

    auto patch_shape = m_controlPoints.shape();
    const int deg_u = m_knotvec_u.getDegree();
    const int deg_v = m_knotvec_v.getDegree();

    NURBSCurveType c(patch_shape[0], deg_u);
    if(isRationalPatch)
    {
      c.makeRational();
    }

    // Find the control points by evaluating each row of the patch
    const auto span_v = m_knotvec_v.findSpan(v);
    const auto N_evals_v = m_knotvec_v.calculateBasisFunctionsBySpan(span_v, v);
    for(int p = 0; p < patch_shape[0]; ++p)
    {
      Point<T, NDIMS + 1> H;
      for(int i = 0; i <= deg_v; ++i)
      {
        const auto offset = span_v - deg_v + i;
        const T weight = isRationalPatch ? m_weights(p, offset) : 1.0;
        const auto& controlPoint = m_controlPoints(p, offset);

        for(int j = 0; j < NDIMS; ++j)
        {
          H[j] += N_evals_v[i] * weight * controlPoint[j];
        }
        H[NDIMS] += N_evals_v[i] * weight;
      }

      for(int j = 0; j < NDIMS; ++j)
      {
        c[p][j] = H[j] / H[NDIMS];
      }

      if(isRationalPatch)
      {
        c.setWeight(p, H[NDIMS]);
      }
    }

    c.setKnots(m_knotvec_u);

    return c;
  }

  /*!
   * \brief Evaluate the NURBS patch geometry and the first \a d derivatives at parameter \a u, \a v
   *
   * \param [in] u The parameter value on the first axis
   * \param [in] v The parameter value on the second axis
   * \param [in] d The number of derivatives to evaluate
   * \param [out] ders A matrix of size d+1 x d+1 containing the derivatives
   * 
   * ders[i][j] is the derivative of S with respect to u i times and v j times.
   *  For consistency, ders[0][0] contains the evaluation point stored as a vector
   *  
   * Implementation adapted from Algorithm A3.6 on p. 111 of "The NURBS Book".
   * Rational derivatives from Algorithm A4.4 on p. 137 of "The NURBS Book".
   * 
   * \pre Requires \a u, v be in the span of the knots (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  void evaluateDerivatives(T u, T v, int d, axom::Array<VectorType, 2>& ders) const
  {
    u = axom::utilities::clampVal(u, getMinKnot_u(), getMaxKnot_u());
    v = axom::utilities::clampVal(v, getMinKnot_v(), getMaxKnot_v());

    const int deg_u = getDegree_u();
    const int du = axom::utilities::min(d, deg_u);

    const int deg_v = getDegree_v();
    const int dv = axom::utilities::min(d, deg_v);

    // Matrix for derivatives
    ders.resize(d + 1, d + 1);
    ders.fill(VectorType(0.0));

    // Matrix for derivatives of homogeneous surface
    //  Store w_{ui, uj} in Awders[i][j][NDIMS]
    axom::Array<Point<T, NDIMS + 1>, 2> Awders(d + 1, d + 1);
    Awders.fill(Point<T, NDIMS + 1>::zero());

    const bool isCurveRational = isRational();

    // Find the span of the knot vectors and basis function derivatives
    const auto span_u = m_knotvec_u.findSpan(u);
    const auto N_evals_u = m_knotvec_u.derivativeBasisFunctionsBySpan(span_u, u, du);

    const auto span_v = m_knotvec_v.findSpan(v);
    const auto N_evals_v = m_knotvec_v.derivativeBasisFunctionsBySpan(span_v, v, dv);

    for(int k = 0; k <= du; ++k)
    {
      axom::Array<Point<T, NDIMS + 1>> temp(deg_v + 1);

      for(int s = 0; s <= deg_v; ++s)
      {
        temp[s] = Point<T, NDIMS + 1>::zero();
        for(int r = 0; r <= deg_u; ++r)
        {
          auto the_weight = isCurveRational ? m_weights(span_u - deg_u + r, span_v - deg_v + s) : 1.0;
          auto& the_pt = m_controlPoints(span_u - deg_u + r, span_v - deg_v + s);

          for(int n = 0; n < NDIMS; ++n)
          {
            temp[s][n] += N_evals_u[k][r] * the_weight * the_pt[n];
          }
          temp[s][NDIMS] += N_evals_u[k][r] * the_weight;
        }
      }

      int dd = axom::utilities::min(d - k, dv);
      for(int l = 0; l <= dd; ++l)
      {
        for(int s = 0; s <= deg_v; ++s)
        {
          for(int n = 0; n < NDIMS + 1; ++n)
          {
            Awders[k][l][n] += N_evals_v[l][s] * temp[s][n];
          }
        }
      }
    }

    // Compute the derivatives of the homogeneous surface
    for(int k = 0; k <= d; ++k)
    {
      for(int l = 0; l <= d - k; ++l)
      {
        auto v = Awders[k][l];

        for(int j = 0; j <= l; ++j)
        {
          auto bin = axom::utilities::binomialCoefficient(l, j);
          for(int n = 0; n < NDIMS; ++n)
          {
            v[n] -= bin * Awders[0][j][NDIMS] * ders[k][l - j][n];
          }
        }

        for(int i = 1; i <= k; ++i)
        {
          auto bin = axom::utilities::binomialCoefficient(k, i);
          for(int n = 0; n < NDIMS; ++n)
          {
            v[n] -= bin * Awders[i][0][NDIMS] * ders[k - i][l][n];
          }

          auto v2 = Point<T, NDIMS + 1>::zero();
          for(int j = 1; j <= l; ++j)
          {
            auto bin = axom::utilities::binomialCoefficient(l, j);
            for(int n = 0; n < NDIMS; ++n)
            {
              v2[n] += bin * Awders[i][j][NDIMS] * ders[k - i][l - j][n];
            }
          }

          for(int n = 0; n < NDIMS; ++n)
          {
            v[n] -= bin * v2[n];
          }
        }

        for(int n = 0; n < NDIMS; ++n)
        {
          ders[k][l][n] = v[n] / Awders[0][0][NDIMS];
        }
      }
    }
  }

  /*!
   * \brief Evaluates all first derivatives of the NURBS patch geometry at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the NURBS patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   *
   * \pre We require evaluation of the patch at \a u and \a v between 0 and 1
   */
  void evaluateFirstDerivatives(T u, T v, PointType& eval, VectorType& Du, VectorType& Dv) const
  {
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 1, ders);

    eval = PointType(ders[0][0].array());
    Du = ders[1][0];
    Dv = ders[0][1];
  }

  /*!
   * \brief Evaluates all linear derivatives of the NURBS patch geometry at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the NURBS patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   *
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  void evaluateLinearDerivatives(T u,
                                 T v,
                                 PointType& eval,
                                 VectorType& Du,
                                 VectorType& Dv,
                                 VectorType& DuDv) const
  {
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 1, ders);

    eval = PointType(ders[0][0]);
    Du = ders[1][0];
    Dv = ders[0][1];
    DuDv = ders[1][1];
  }

  /*!
   * \brief Evaluates all second derivatives of the NURBS patch geometry at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the NURBS patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDu The vector value of S_uu(u, v)
   * \param [out] DvDv The vector value of S_vv(u, v)
   * \param [out] DuDv The vector value of S_uu(u, v)
   *
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
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
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 2, ders);

    eval = PointType(ders[0][0]);
    Du = ders[1][0];
    Dv = ders[0][1];
    DuDu = ders[2][0];
    DvDv = ders[0][2];
    DuDv = ders[1][1];
  }

  /*!
   * \brief Computes a tangent in u of the NURBS patch geometry at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   *
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  VectorType du(T u, T v) const
  {
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 1, ders);

    return ders[1][0];
  }

  /*!
   * \brief Computes a tangent in v of the NURBS patch geometry at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   *
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  VectorType dv(T u, T v) const
  {
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 1, ders);

    return ders[0][1];
  }

  /*!
   * \brief Computes the second derivative in u of the NURBS patch geometry patch at (\a u, \a v)
   * 
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * 
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  VectorType dudu(T u, T v) const
  {
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 2, ders);

    return ders[2][0];
  }

  /*!
   * \brief Computes the second derivative in v of the NURBS patch geometry patch at (\a u, \a v)
   * 
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * 
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  VectorType dvdv(T u, T v) const
  {
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 2, ders);

    return ders[0][2];
  }

  /*!
   * \brief Computes the mixed second derivative in u and v of the NURBS patch geometry at (\a u, \a v)
   * 
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * 
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  VectorType dudv(T u, T v) const
  {
    axom::Array<VectorType, 2> ders;
    evaluateDerivatives(u, v, 2, ders);

    return ders[1][1];
  }

  /*!
   * \brief Computes the mixed second derivative in u and v of the NURBS patch geometry at (\a u, \a v)
   * 
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * 
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  VectorType dvdu(T u, T v) const { return dudv(u, v); }

  /*!
   * \brief Computes the normal vector to the NURBS patch geometry at (\a u, \a v)
   * 
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * 
   * \pre We require evaluation of the patch at \a u and \a v (up to a small tolerance)
   * 
   * \note If u/v is outside the knot span up this tolerance, it is clamped to the span
   */
  VectorType normal(T u, T v) const
  {
    PointType eval;
    VectorType Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);

    return VectorType::cross_product(Du, Dv);
  }

  /*!
   * \brief Calculate the average normal for the (untrimmed) patch
   * 
   * \param [in] npts The number of quadrature nodes used in each component integral
   *
   * Algorithm from "Mean normal vector to a surface bounded by Bézier curves"
   *  by Kenji Ueda, 1996
   * 
   * Projects the 4 boundary curves of the patch along each coordinate axis, 
   *  then computes the 2D area of that projection to get the corresponding
   *  component of the average surface normal.
   *  
   * Evaluates the integral with Gauss-Legendre quadrature on each boundary curve
   * 
   * \return The calculated mean surface normal
   */
  VectorType calculateUntrimmedPatchNormal(int npts = 20) const
  {
    SLIC_ASSERT(NDIMS == 3);

    VectorType ret_vec;
    auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };

    // Set up the correct sizes and weights of the bounding curves
    axom::Array<NURBSCurve<T, 2>> boundingPoly(4);

    const int npts_u = getNumControlPoints_u();
    const int npts_v = getNumControlPoints_v();

    boundingPoly[0].setParameters(npts_v, getDegree_v());
    boundingPoly[0].setKnots(getKnots_v());

    boundingPoly[1].setParameters(npts_u, getDegree_u());
    boundingPoly[1].setKnots(getKnots_u());

    boundingPoly[2].setParameters(npts_v, getDegree_v());
    boundingPoly[2].setKnots(getKnots_v());

    boundingPoly[3].setParameters(npts_u, getDegree_u());
    boundingPoly[3].setKnots(getKnots_u());

    if(isRational())
    {
      for(int i = 0; i < 4; ++i)
      {
        boundingPoly[i].makeRational();
      }

      for(int m = 0; m < npts_v; ++m)
      {
        boundingPoly[0].setWeight(m, m_weights(0, npts_v - 1 - m));
        boundingPoly[2].setWeight(m, m_weights(npts_u - 1, m));
      }

      for(int n = 0; n < npts_u; ++n)
      {
        boundingPoly[1].setWeight(n, m_weights(npts_u - 1 - n, npts_v - 1));
        boundingPoly[3].setWeight(n, m_weights(n, 0));
      }
    }

    // Get each component by projecting boundaries onto the other coordinate axes
    for(int N = 0; N < 3; ++N)
    {
      int ind_3d = 0;
      for(int i = 0; i < 2; ++i, ++ind_3d)
      {
        // Skip the corresponding coordinate used to access the 3D point
        if(ind_3d == N)
        {
          --i;
          continue;
        }

        for(int m = 0; m < npts_v; ++m)
        {
          boundingPoly[0][m][i] = m_controlPoints(0, npts_v - 1 - m)[ind_3d];
          boundingPoly[2][m][i] = m_controlPoints(npts_u - 1, m)[ind_3d];
        }

        for(int n = 0; n < npts_u; ++n)
        {
          boundingPoly[1][n][i] = m_controlPoints(npts_u - 1 - n, npts_v - 1)[ind_3d];
          boundingPoly[3][n][i] = m_controlPoints(n, 0)[ind_3d];
        }
      }

      // Find the area of the resulting projection
      ret_vec[N] = evaluate_area_integral(boundingPoly, const_integrand, npts);
    }

    // Need to flip the y-component to account for the flipped projection
    ret_vec[1] = -ret_vec[1];
    return ret_vec;
  }
  ///@}

  ///@{
  //!  @name Methods for (trimmed) Surface Geometry.
  //!
  //! These methods operate on the visible geometry of the NURBS surface
  //!  as defined by the patch geometry and the trimming curves
  //!  which define visibility.

  /*!
   * \brief Get array of trimming curves
   */
  const TrimmingCurveVec& getTrimmingCurves() const { return m_trimmingCurves; }

  /// \brief Get mutable array of trimming curves
  TrimmingCurveVec& getTrimmingCurves() { return m_trimmingCurves; }

  /// \brief Get a trimming curve by index
  const TrimmingCurveType& getTrimmingCurve(int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < m_trimmingCurves.size());
    return m_trimmingCurves[idx];
  }

  /// \brief Add a trimming curve
  void addTrimmingCurve(const TrimmingCurveType& curve)
  {
    m_isTrimmed = true;
    m_trimmingCurves.push_back(curve);
  }

  /// \brief Add array of trimming curves
  void addTrimmingCurves(const TrimmingCurveVec& curves)
  {
    m_isTrimmed = true;
    for(int i = 0; i < curves.size(); ++i)
    {
      m_trimmingCurves.push_back(curves[i]);
    }
  }

  /// \brief Set the array of trimming curves
  void setTrimmingCurves(const TrimmingCurveVec& curves)
  {
    m_isTrimmed = true;
    m_trimmingCurves = curves;
  }

  /// \brief Clear trimming curves, but DON'T mark as untrimmed
  void clearTrimmingCurves() { m_trimmingCurves.clear(); }

  /// \brief Get number of trimming curves
  int getNumTrimmingCurves() const { return m_trimmingCurves.size(); }

  /// \brief use boolean flag for trimmed-ness
  bool isTrimmed() const { return m_isTrimmed; }

  /// \brief Mark as trimmed
  void markAsTrimmed() { m_isTrimmed = true; }

  /// \brief Delete all trimming curves
  void makeUntrimmed()
  {
    m_isTrimmed = false;
    m_trimmingCurves.clear();
  }

  /// \brief Make trivially trimmed by adding trimming curves at each boundary
  void makeTriviallyTrimmed()
  {
    if(isTrimmed())
    {
      m_trimmingCurves.clear();
    }

    const double min_u = m_knotvec_u[0];
    const double max_u = m_knotvec_u[m_knotvec_u.getNumKnots() - 1];

    const double min_v = m_knotvec_v[0];
    const double max_v = m_knotvec_v[m_knotvec_v.getNumKnots() - 1];

    // For each min/max u/v, add a straight trimming curve along the boundary

    // Bottom
    addTrimmingCurve(TrimmingCurveType::make_linear_segment_nurbs({min_u, min_v}, {max_u, min_v}));

    // Top
    addTrimmingCurve(TrimmingCurveType::make_linear_segment_nurbs({max_u, max_v}, {min_u, max_v}));

    // Left
    addTrimmingCurve(TrimmingCurveType::make_linear_segment_nurbs({min_u, max_v}, {min_u, min_v}));

    // Right
    addTrimmingCurve(TrimmingCurveType::make_linear_segment_nurbs({max_u, min_v}, {max_u, max_v}));

    markAsTrimmed();
  }

  /*!
   * \brief Check if a parameter point is visible on the NURBS patch via a trim test
   *
   * \param [in] u The parameter value on the first axis
   * \param [in] v The parameter value on the second axis
   * 
   * Checks for containment of the parameter point in 
   *  the collection of trimming curves via an even-odd rule.
   * 
   * If the collection of trimming curves does not form closed loops, 
   *  then the (now fractional) generalized winding number is rounded to the
   *  nearest integer before the even-odd rule is applied.
   */
  bool isVisible(T u, T v) const
  {
    if(!isTrimmed())
    {
      return (m_knotvec_u.isValidParameter(u) && m_knotvec_v.isValidParameter(v));
    }

    ParameterPointType uv = {u, v};

    double gwn = 0.0;
    for(const auto& curve : m_trimmingCurves)
    {
      bool isOnThisCurve = false;
      gwn += detail::nurbs_winding_number(uv, curve, isOnThisCurve);

      if(isOnThisCurve)
      {
        return true;
      }
    }

    return std::lround(gwn) % 2 != 0;
  }

  /*!
   * \brief Restrict the edges of a NURBS surface to the given parameter values, if necessary
   *
   * \param [in] min_u The minimum value of the u parameter
   * \param [in] max_u The maximum value of the u parameter
   * \param [in] min_v The minimum value of the v parameter
   * \param [in] max_v The maximum value of the v parameter
   * \param [in] normalize If true, normalize the patch to the range [0, 1]^2
   * 
   * If a given parameter is less/greater than the minimum/maximum knot value,
   *  then the patch is not changed along that direction and axis
   * 
   * \pre Requires that min_u < max_u and min_v < max_v  
   * \post A patch whose parameter space is a subset of [min_u, max_u] x [min_v, max_v]
   * 
   * \sa NURBSPatch::split()
   */
  void clip(T min_u, T max_u, T min_v, T max_v, bool normalizeParameters = false)
  {
    SLIC_ASSERT(min_u < max_u);
    SLIC_ASSERT(min_v < max_v);
    NURBSPatch dummy_patch;

    if(min_u > getMinKnot_u())
    {
      this->split_u(min_u, dummy_patch, *this);
    }
    if(min_v > getMinKnot_v())
    {
      this->split_v(min_v, dummy_patch, *this);
    }
    if(max_u < getMaxKnot_u())
    {
      this->split_u(max_u, *this, dummy_patch);
    }
    if(max_v < getMaxKnot_v())
    {
      this->split_v(max_v, *this, dummy_patch);
    }

    if(normalizeParameters)
    {
      normalize();
    }
  }

  ///@}

  ///@{
  /// \name Functions dealing with (untrimmed) patch subdivision

  /*!
   * \brief Splits the NURBS patch geometry (at each internal knot) into several Bezier patches
   * 
   * If either degree_u or degree_v is zero, the resulting Bezier patches along 
   *  that axis will be disconnected and order 0
   * 
   * This method ignores any trimming curves in the patch, 
   *  and returns all extracted patches of the untrimmed patch.
   * 
   * Algorithm A5.7 on p. 177 of "The NURBS Book"
   * 
   * \return An array of Bezier patches ordered lexicographically (in v, then u)
   */
  axom::Array<BezierPatch<T, NDIMS>> extractBezier() const
  {
    const bool isRationalPatch = isRational();

    const auto n = getNumControlPoints_u() - 1;
    const auto p = getDegree_u();
    const auto kp = m_knotvec_u.getNumKnotSpans();

    const auto m = getNumControlPoints_v() - 1;
    const auto q = getDegree_v();
    const auto kq = m_knotvec_v.getNumKnotSpans();

    axom::Array<NURBSPatch<T, NDIMS>> strips(kp);
    for(int i = 0; i < strips.size(); ++i)
    {
      strips[i].setParameters(p + 1, m + 1, p, q);
      if(isRationalPatch)
      {
        strips[i].makeRational();
      }
    }

    axom::Array<T> alphas(axom::utilities::max(0, axom::utilities::max(p - 1, q - 1)));

    // Do Bezier extraction on the u-axis, which returns a collection of Bezier strips
    if(p == 0)
    {
      for(int i = 0; i < n + 1; ++i)
      {
        for(int row = 0; row < m + 1; ++row)
        {
          strips[i](0, row) = m_controlPoints(i, row);
          if(isRationalPatch)
          {
            strips[i].setWeight(0, row, m_weights(i, row));
          }
        }
      }
    }
    else
    {
      int a = p;
      int b = p + 1;
      int ns = 0;

      for(int i = 0; i <= p; ++i)
      {
        for(int row = 0; row <= m; ++row)
        {
          strips[ns](i, row) = m_controlPoints(i, row);
          if(isRationalPatch)
          {
            strips[ns].setWeight(i, row, m_weights(i, row));
          }
        }
      }

      while(b < n + p + 1)
      {
        // Get multiplicity of the knot
        int i = b;
        while(b < n + p + 1 && m_knotvec_u[b] == m_knotvec_u[b + 1])
        {
          ++b;
        }
        int mult = b - i + 1;

        if(mult < p)
        {
          // Get the numerator and the alphas
          T numer = m_knotvec_u[b] - m_knotvec_u[a];

          for(int j = p; j > mult; --j)
          {
            alphas[j - mult - 1] = numer / (m_knotvec_u[a + j] - m_knotvec_u[a]);
          }

          // Do the knot insertion in place
          for(int j = 1; j <= p - mult; ++j)
          {
            int save = p - mult - j;
            int s = mult + j;

            for(int k = p; k >= s; --k)
            {
              T alpha = alphas[k - s];
              for(int row = 0; row <= m; ++row)
              {
                T weight_k = isRationalPatch ? strips[ns].getWeight(k, row) : 1.0;
                T weight_km1 = isRationalPatch ? strips[ns].getWeight(k - 1, row) : 1.0;

                if(isRationalPatch)
                {
                  strips[ns].setWeight(k, row, alpha * weight_k + (1.0 - alpha) * weight_km1);
                }

                for(int N = 0; N < NDIMS; ++N)
                {
                  strips[ns](k, row)[N] = (alpha * strips[ns](k, row)[N] * weight_k +
                                           (1.0 - alpha) * strips[ns](k - 1, row)[N] * weight_km1) /
                    (isRationalPatch ? strips[ns].getWeight(k, row) : 1.0);
                }
              }
            }

            if(b < n + p + 1)
            {
              for(int row = 0; row <= m; ++row)
              {
                strips[ns + 1](save, row) = strips[ns](p, row);
                if(isRationalPatch)
                {
                  strips[ns + 1].setWeight(save, row, strips[ns].getWeight(p, row));
                }
              }
            }
          }
        }

        ++ns;

        if(b < n + p + 1)
        {
          for(int j = p - mult; j <= p; ++j)
          {
            for(int row = 0; row <= m; ++row)
            {
              strips[ns](j, row) = m_controlPoints(b - p + j, row);
              if(isRationalPatch)
              {
                strips[ns].setWeight(j, row, m_weights(b - p + j, row));
              }
            }
          }
          a = b;
          b++;
        }
      }
    }

    // For each strip, do Bezier extraction on the v-axis
    axom::Array<BezierPatch<T, NDIMS>> beziers(kp * kq);
    for(int i = 0; i < beziers.size(); ++i)
    {
      beziers[i].setOrder(p, q);
      if(isRationalPatch)
      {
        beziers[i].makeRational();
      }
    }

    for(int s_i = 0; s_i < strips.size(); ++s_i)
    {
      auto& strip = strips[s_i];
      int n_i = strip.getNumControlPoints_u() - 1;
      int nb = s_i * m_knotvec_v.getNumKnotSpans();

      // Handle this case separately
      if(q == 0)
      {
        for(int i = 0; i < m + 1; ++i)
        {
          for(int col = 0; col < n_i + 1; ++col)
          {
            beziers[nb](col, 0) = strip(col, i);
            if(isRationalPatch)
            {
              beziers[nb].setWeight(col, 0, strip.getWeight(col, i));
            }
          }

          ++nb;
        }

        continue;
      }

      int a = q;
      int b = q + 1;

      for(int i = 0; i <= q; ++i)
      {
        for(int col = 0; col <= n_i; ++col)
        {
          beziers[nb](col, i) = strip(col, i);
          if(isRationalPatch)
          {
            beziers[nb].setWeight(col, i, strip.getWeight(col, i));
          }
        }
      }

      while(b < m + q + 1)
      {
        // Get multiplicity of the knot
        int i = b;
        while(b < m + q + 1 && m_knotvec_v[b] == m_knotvec_v[b + 1])
        {
          ++b;
        }
        int mult = b - i + 1;

        if(mult < q)
        {
          // Get the numerator and the alphas
          T numer = m_knotvec_v[b] - m_knotvec_v[a];

          for(int j = q; j > mult; --j)
          {
            alphas[j - mult - 1] = numer / (m_knotvec_v[a + j] - m_knotvec_v[a]);
          }

          // Do the knot insertion in place
          for(int j = 1; j <= q - mult; ++j)
          {
            int save = q - mult - j;
            int s = mult + j;

            for(int k = q; k >= s; --k)
            {
              T alpha = alphas[k - s];
              for(int col = 0; col <= n_i; ++col)
              {
                T weight_k = isRationalPatch ? beziers[nb].getWeight(col, k) : 1.0;
                T weight_km1 = isRationalPatch ? beziers[nb].getWeight(col, k - 1) : 1.0;

                if(isRationalPatch)
                {
                  beziers[nb].setWeight(col, k, alpha * weight_k + (1.0 - alpha) * weight_km1);
                }

                for(int N = 0; N < NDIMS; ++N)
                {
                  beziers[nb](col, k)[N] = (alpha * beziers[nb](col, k)[N] * weight_k +
                                            (1.0 - alpha) * beziers[nb](col, k - 1)[N] * weight_km1) /
                    (isRationalPatch ? beziers[nb].getWeight(col, k) : 1.0);
                }
              }
            }

            if(b < m + q + 1)
            {
              for(int col = 0; col <= n_i; ++col)
              {
                beziers[nb + 1](col, save) = beziers[nb](col, q);
                if(isRationalPatch)
                {
                  beziers[nb + 1].setWeight(col, save, beziers[nb].getWeight(col, q));
                }
              }
            }
          }
        }

        ++nb;

        if(b < m + q + 1)
        {
          for(int j = q - mult; j <= q; ++j)
          {
            for(int col = 0; col <= n_i; ++col)
            {
              beziers[nb](col, j) = strip(col, b - q + j);
              if(isRationalPatch)
              {
                beziers[nb].setWeight(col, j, strip.getWeight(col, b - q + j));
              }
            }
          }
          a = b;
          b++;
        }
      }
    }

    return beziers;
  }

  /*!
   * \brief Splits the NURBS patch geometry into several one-span trimmed NURBS surfaces
   *  
   * \return An array of Bezier patches ordered lexicographically (in v, then u)
   * 
   * \sa extractBezier
   */
  axom::Array<NURBSPatch<T, NDIMS>> extractTrimmedBezier() const
  {
    // Loop over the original set of trimming curves and split them over the knot vectors
    axom::Array<T> knot_vals_u = getKnots_u().getUniqueKnots();
    axom::Array<T> knot_vals_v = getKnots_v().getUniqueKnots();

    const auto num_knot_span_u = knot_vals_u.size() - 1;
    const auto num_knot_span_v = knot_vals_v.size() - 1;

    axom::Array<NURBSPatch<T, 3>> split_patches(num_knot_span_u * num_knot_span_v);

    // This method is nominally faster if the Bezier extraction routine is called separately,
    //  as this avoids a handful of repeated calculations
    split_patches[0] = *this;
    for(int i = 0; i < num_knot_span_u - 1; ++i)
    {
      split_patches[i * num_knot_span_v].split_u(knot_vals_u[i + 1],
                                                 split_patches[i * num_knot_span_v],
                                                 split_patches[(i + 1) * num_knot_span_v]);
    }

    for(int i = 0; i < num_knot_span_u; ++i)
    {
      for(int j = 0; j < num_knot_span_v - 1; ++j)
      {
        split_patches[i * num_knot_span_v + j].split_v(knot_vals_v[j + 1],
                                                       split_patches[i * num_knot_span_v + j],
                                                       split_patches[i * num_knot_span_v + j + 1]);
      }
    }

    return split_patches;
  }

  /*!
   * \brief Calculate the average normal for the trimmed patch
   * 
   * \param [in] npts The number of quadrature nodes used in each component integral
   *
   * Decomposes the NURBS surface into trimmed Bezier components (to ensure differentiability of the integrand) 
   *  and evaluates the integral numerically on each component using trimming curves
   * 
   * Evaluates the integral with Gauss-Legendre quadrature on each boundary curve
   * 
   * \return The calculated mean surface normal
   */
  VectorType calculateTrimmedPatchNormal(int npts = 20) const
  {
    SLIC_ASSERT(NDIMS == 3);

    VectorType ret_vec;

    // Split the patch along the unique knot values to improve convergence
    for(const auto& nPatch : extractTrimmedBezier())
    {
      // Integrate the surface normal over the patches
      ret_vec += evaluate_area_integral(
        nPatch.getTrimmingCurves(),
        [&nPatch](Point2D x) -> Vector<T, 3> { return nPatch.normal(x[0], x[1]); },
        npts);
    }

    return ret_vec;
  }
  //@}

  ///@{
  /// \name Functions dealing with trimmed patch subdivision

  /*!
    * \brief Splits a NURBS surface into four NURBS patches
    *
    * \param [in] u parameter value at which to bisect on the first axis
    * \param [in] v parameter value at which to bisect on the second axis
    * \param [out] p1 First output NURBS patch
    * \param [out] p2 Second output NURBS patch
    * \param [out] p3 Third output NURBS patch
    * \param [out] p4 Fourth output NURBS patch
    * \param [in] normalizeParameters If true, normalize the knot span of the
    *                                   output patches to the range [0, 1]^2
    *
    *          v = v_max
    *          ----------------------
    *          |         |          |
    *          |   p3    |    p4    |
    *          |         |          |
    *          --------(u,v)---------
    *          |         |          |
    *          |   p1    |    p2    |
    *          |         |          |
    *          ---------------------- u = u_max
    *  u/v_min
    * 
    * \note If u/v is not strictly interior to the knot span, will return an invalid NURBS
    *  for the invalid portion and the original surface for the rest
    *
    * \return True if and only if the patch was split (i.e., u, v is in the knot span)
    */
  bool split(T u,
             T v,
             NURBSPatch& p1,
             NURBSPatch& p2,
             NURBSPatch& p3,
             NURBSPatch& p4,
             bool normalizeParameters = false) const
  {
    bool wasSplit = true;

    // Bisect the patch along the u direction
    wasSplit = split_u(u, p1, p2, false) && wasSplit;

    // Temporarily store the result in each half and split again
    NURBSPatch p0(p1);
    wasSplit = p0.split_v(v, p1, p3, false) && wasSplit;

    p0 = p2;
    wasSplit = p0.split_v(v, p2, p4, false) && wasSplit;

    if(normalizeParameters)
    {
      p1.normalize();
      p2.normalize();
      p3.normalize();
      p4.normalize();
    }

    return wasSplit;
  }

  /*!
   * \brief Split the NURBS surface in two along the u direction
   *
   * \note If u is not strictly interior to the knot span, will return an invalid NURBS
   *  for the invalid portion and the original surface for the rest
   * 
   * \return True if and only if the patch was split (i.e., u is in the knot span)
   */
  bool split_u(T u, NURBSPatch& p1, NURBSPatch& p2, bool normalizeParameters = false) const
  {
    // If the patch is not valid, return two invalid patches
    if(m_controlPoints.size() == 0)
    {
      p1 = NURBSPatch();
      p2 = NURBSPatch();
      return false;
    }

    // If u is outside hte knot span, return the original patch
    //  and an invalid NURBS patch
    if(!isValidInteriorParameter_u(u))
    {
      if(u <= getMinKnot_u())
      {
        p1 = NURBSPatch();
        p2 = *this;
      }
      else if(u >= getMaxKnot_u())
      {
        p1 = *this;
        p2 = NURBSPatch();
      }

      if(normalizeParameters)
      {
        p1.normalize_u();
        p2.normalize_u();
      }

      return false;
    }

    // Split the untrimmed geometry
    uncheckedSplit_u(u, p1, p2);

    // Split the trimming curves if necessary
    if(isTrimmed())
    {
      constexpr bool splitInU = true;
      splitTrimmingCurves(u, splitInU, p1.getTrimmingCurves(), p2.getTrimmingCurves());
    }

    if(normalizeParameters)
    {
      p1.normalize_u();
      p2.normalize_u();
    }

    return true;
  }

  /*!
   * \brief Split the NURBS surface in two along the v direction
   *
   * \note If v is not strictly interior to the knot span, will return an invalid NURBS
   *  for the invalid portion and the original surface for the rest
   * 
   * \return True if and only if the patch was split (i.e., v is in the knot span)
   */
  bool split_v(T v, NURBSPatch& p1, NURBSPatch& p2, bool normalizeParameters = false) const
  {
    // If the patch is not valid, return two invalid patches
    if(m_controlPoints.size() == 0)
    {
      p1 = NURBSPatch();
      p2 = NURBSPatch();
      return false;
    }

    // If v is outside the knot span, return the original patch
    //  and an invalid NURBS patch
    if(!isValidInteriorParameter_v(v))
    {
      if(v <= getMinKnot_v())
      {
        p1 = NURBSPatch();
        p2 = *this;
      }
      else if(v >= getMaxKnot_v())
      {
        p1 = *this;
        p2 = NURBSPatch();
      }

      if(normalizeParameters)
      {
        p1.normalize_v();
        p2.normalize_v();
      }

      return false;
    }

    // Split the untrimmed geometry
    uncheckedSplit_v(v, p1, p2);

    // Split the trimming curves if necessary
    if(isTrimmed())
    {
      constexpr bool splitInU = false;
      splitTrimmingCurves(v, splitInU, p1.getTrimmingCurves(), p2.getTrimmingCurves());
    }

    if(normalizeParameters)
    {
      p1.normalize_v();
      p2.normalize_v();
    }

    return true;
  }

  /*!
   * \brief For a disk of radius r and center (u, v), split a NURBS surface into the portion inside/outside
   *
   * \param [in] u The x-coordinate of the center of the disk
   * \param [in] v The y-coordinate of the center of the disk
   * \param [in] r The radius of the disk 
   * \param [out] the_disk The NURBS surface inside the disk
   * \param [out] the_rest The NURBS surface outside the disk
   * \param [in] clipDisk If true, the returned disk is clipped to the disk boundary
   */
  void diskSplit(T u, T v, T r, NURBSPatch& the_disk, NURBSPatch& the_rest, bool clipDisk = true) const
  {
    bool isDiskInside = false;
    bool isDiskOutside = false;
    bool ignoreInteriorDisk = false;
    double disk_padding = clipDisk ? 0.0 : getParameterSpaceDiagonal();
    diskSplit(u, v, r, the_disk, the_rest, isDiskInside, isDiskOutside, ignoreInteriorDisk, disk_padding);
  }

  /*!
   * \brief Split a NURBS surface into two by cutting out a disk of radius r centered at (u, v)
   *
   * \param [in] u The x-coordinate of the center of the disk
   * \param [in] v The y-coordinate of the center of the disk
   * \param [in] r The radius of the disk 
   * \param [out] the_disk The NURBS surface inside the disk
   * \param [out] the_rest The NURBS surface outside the disk
   * \param [out] isDiskInside True if the disk is entirely inside the trimming curves
   * \param [out] isDiskOutside True if the disk is entirely outside the trimming curves
   * \param [in] ignoreInteriorDisk If true, don't perform subdivision if disk is entirely inside the trimming curves
   * \param [in] disk_padding How much space is retained around the disk in each direction after clipping
   * 
   * \note Function arguments suited for use in GWN evaluation
   */
  void diskSplit(T u,
                 T v,
                 T r,
                 NURBSPatch& the_disk,
                 NURBSPatch& the_rest,
                 bool& isDiskInside,
                 bool& isDiskOutside,
                 bool ignoreInteriorDisk,
                 double disk_padding) const
  {
    ParameterPointType uv_param({u, v});

    // Copy the control points and weights of the original patch, but not the trimming curves
    the_disk = NURBSPatch(m_controlPoints, m_weights, m_knotvec_u, m_knotvec_v);
    the_disk.markAsTrimmed();

    the_rest = *this;

    // the_rest needs trimming curves from the original patch, if any
    if(!isTrimmed())
    {
      the_rest.makeTriviallyTrimmed();
    }

    // Intersect all trimming curves with a circle of radius r, centered at (u, v).
    //  Record each intersection with the circle, and split the trimming curve at each intersection
    primal::Sphere<T, 2> circle_obj(ParameterPointType {u, v}, r);
    TrimmingCurveVec split_trimming_curves;
    TrimmingCurveVec circle_trimming_curves;

    axom::Array<T> circle_params;
    circle_params.push_back(0.0);
    for(const auto& curve : the_rest.m_trimmingCurves)
    {
      axom::Array<T> curve_params;
      // Compute intersections for each subcurve individually to avoid
      //  a circular dependency with intersect.hpp
      {
        // Default parameters for intersection routine
        const double sq_tol = 1e-14;
        const double EPS = 1e-6;

        // Extract the Bezier curves of the NURBS curve, checking each for intersection
        axom::Array<T> knot_vals = curve.getKnots().getUniqueKnots();
        const auto beziers = curve.extractBezier();
        for(int i = 0; i < beziers.size(); ++i)
        {
          axom::Array<T> temp_curve_p;
          axom::Array<T> temp_circle_p;
          detail::intersect_circle_bezier(circle_obj,
                                          beziers[i],
                                          temp_circle_p,
                                          temp_curve_p,
                                          sq_tol,
                                          EPS,
                                          beziers[i].getOrder(),
                                          0.,
                                          1.);

          // If the number of recorded intersection points is too great (as defined by Bezout's theorem),
          //   then they can be assumed to be completely overlapping, and no intersections are recorded.
          if(temp_curve_p.size() > 6 * beziers[i].getOrder())
          {
            continue;
          }

          // Scale the intersection parameters back into the span of the NURBS curve
          for(int j = 0; j < temp_curve_p.size(); ++j)
          {
            circle_params.push_back(temp_circle_p[j]);
            curve_params.push_back(knot_vals[i] + temp_curve_p[j] * (knot_vals[i + 1] - knot_vals[i]));
          }
        }
      }

      // Split all trimming curves at the intersection points
      if(curve_params.size() > 0)
      {
        // Sorting this keeps the splitting logic simpler
        std::sort(curve_params.begin(), curve_params.end());

        TrimmingCurveType c1, c2(curve);
        for(const auto& param : curve_params)
        {
          if(param <= c2.getMinKnot() || param >= c2.getMaxKnot())
          {
            continue;
          }

          c2.split(param, c1, c2);
          split_trimming_curves.push_back(c1);
        }
        split_trimming_curves.push_back(c2);
      }
      else
      {
        split_trimming_curves.push_back(curve);
      }
    }
    circle_params.push_back(2 * M_PI);

    // Handle special cases where 0 intersections are recorded
    isDiskInside = isDiskOutside = false;
    if(circle_params.size() == 2)
    {
      // If the circle is entirely inside the trimming curves,
      //  the_disk is a complete disk
      if(isVisible(u, v))
      {
        isDiskInside = true;

        if(ignoreInteriorDisk)
        {
          the_disk.m_trimmingCurves.clear();
          return;
        }

        TrimmingCurveType c1 = TrimmingCurveType::make_circular_arc_nurbs(0.0, 2 * M_PI, u, v, r);

        the_disk.m_trimmingCurves.clear();
        the_disk.addTrimmingCurve(c1);

        // Clip the_disk according to the width of the disk and the padding parameter
        the_disk.uncheckedClip(u - r - disk_padding,
                               u + r + disk_padding,
                               v - r - disk_padding,
                               v + r + disk_padding);

        c1.reverseOrientation();
        the_rest.addTrimmingCurve(c1);

        return;
      }
      else
      {
        // If the circle is entirely outside the trimming curves,
        //  the_rest is unchanged and the_disk is empty
        isDiskOutside = true;
        the_disk.m_trimmingCurves.clear();

        return;
      }
    }

    // Sort the circle parameters
    std::sort(circle_params.begin(), circle_params.end());

    for(int i = 0; i < circle_params.size() - 1; ++i)
    {
      // Skip any duplicate parameters
      if(circle_params[i + 1] - circle_params[i] < 1e-10)
      {
        continue;
      }

      // Determine if the circle arc is kept by the original surface
      ParameterPointType mid_arc_point {
        u + r * std::cos(0.5 * (circle_params[i] + circle_params[i + 1])),
        v + r * std::sin(0.5 * (circle_params[i] + circle_params[i + 1]))};
      bool isArcVisible = isVisible(mid_arc_point[0], mid_arc_point[1]);

      if(isArcVisible)
      {
        circle_trimming_curves.push_back(
          TrimmingCurveType::make_circular_arc_nurbs(circle_params[i], circle_params[i + 1], u, v, r));
      }
    }

    // Clear the trimming curves from each patch.
    //  Let the_rest be the "big" patch, and the_disk be the "small" patch
    //  the_rest gets all trimming curves *outside* the circle, and the reverse of all circle curves
    //  the_disk gets all trimming curves *inside* the circle, and all circle curves

    the_rest.m_trimmingCurves.clear();
    the_disk.m_trimmingCurves.clear();

    for(const auto& curve : split_trimming_curves)
    {
      auto curve_midpoint = curve.evaluate(0.5 * (curve.getMinKnot() + curve.getMaxKnot()));
      bool isInDisk = circle_obj.computeSignedDistance(curve_midpoint) < 0;

      if(isInDisk)
      {
        the_disk.addTrimmingCurve(curve);
      }
      else
      {
        the_rest.addTrimmingCurve(curve);
      }
    }

    for(auto& curve : circle_trimming_curves)
    {
      the_disk.addTrimmingCurve(curve);
      curve.reverseOrientation();
      the_rest.addTrimmingCurve(curve);
    }

    // Clip the_disk according to the width of the disk and the padding parameter
    the_disk.uncheckedClip(u - r - disk_padding,
                           u + r + disk_padding,
                           v - r - disk_padding,
                           v + r + disk_padding);
  }
  //@}

  /*!
     * \brief Simple formatted print of a NURBS Patch instance
     *
     * \param os The output stream to write to
     * \return A reference to the modified ostream
     */
  std::ostream& print(std::ostream& os) const
  {
    auto patch_shape = m_controlPoints.shape();

    int deg_u = m_knotvec_u.getDegree();
    int deg_v = m_knotvec_v.getDegree();

    int nkts_u = m_knotvec_u.getNumKnots();
    int nkts_v = m_knotvec_v.getNumKnots();

    os << "{ degree (" << deg_u << ", " << deg_v << ") NURBS Patch, ";
    os << "control points [";
    for(int p = 0; p < patch_shape[0]; ++p)
    {
      for(int q = 0; q < patch_shape[1]; ++q)
      {
        os << m_controlPoints(p, q)
           << ((p < patch_shape[0] - 1 || q < patch_shape[1] - 1) ? "," : "]");
      }
    }

    if(isRational())
    {
      os << ", weights [";
      for(int p = 0; p < patch_shape[0]; ++p)
      {
        for(int q = 0; q < patch_shape[1]; ++q)
        {
          os << m_weights(p, q) << ((p < patch_shape[0] - 1 || q < patch_shape[1] - 1) ? "," : "]");
        }
      }
    }

    os << ", knot vector u [";
    for(int i = 0; i < nkts_u; ++i)
    {
      os << m_knotvec_u[i] << ((i < nkts_u - 1) ? "," : "]");
    }

    os << ", knot vector v [";
    for(int i = 0; i < nkts_v; ++i)
    {
      os << m_knotvec_v[i] << ((i < nkts_v - 1) ? "," : "]");
    }

    if(isTrimmed())
    {
      os << ", trimming curves [";
      for(int i = 0; i < m_trimmingCurves.size(); ++i)
      {
        os << m_trimmingCurves[i];
        if(i < m_trimmingCurves.size() - 1)
        {
          os << ", ";
        }
      }
      os << "]";
    }

    return os;
  }

private:
  /// \brief Private function to rescale trimming curves from (a, b) to (c, d) in x
  /// \warning Does not check that the resulting curves are valid
  void rescaleTrimmingCurves_u(T a, T b, T c, T d)
  {
    SLIC_ASSERT(a < b);
    SLIC_ASSERT(c < d);

    for(auto& curve : m_trimmingCurves)
    {
      for(int i = 0; i < curve.getNumControlPoints(); ++i)
      {
        curve[i][0] = c + (d - c) * (curve[i][0] - a) / (b - a);
      }
    }
  }

  /// \brief Private function to rescale trimming curves from (a, b) to (c, d) in y
  /// \warning Does not check that the resulting curves are valid
  void rescaleTrimmingCurves_v(T a, T b, T c, T d)
  {
    SLIC_ASSERT(a < b);
    SLIC_ASSERT(c < d);

    for(auto& curve : m_trimmingCurves)
    {
      for(int i = 0; i < curve.getNumControlPoints(); ++i)
      {
        curve[i][1] = c + (d - c) * (curve[i][1] - a) / (b - a);
      }
    }
  }

  /// \brief Clip the edges of the patch along the given knot values, if necessary
  ///  but do so *without* checking any of the existing trimming curves for efficiency
  /// \sa NURBSPatch::clip()
  void uncheckedClip(T min_u, T max_u, T min_v, T max_v)
  {
    SLIC_ASSERT(min_u < max_u);
    SLIC_ASSERT(min_v < max_v);
    NURBSPatch dummy_patch;

    if(min_u > getMinKnot_u()) this->uncheckedSplit_u(min_u, dummy_patch, *this);
    if(min_v > getMinKnot_v()) this->uncheckedSplit_v(min_v, dummy_patch, *this);
    if(max_u < getMaxKnot_u()) this->uncheckedSplit_u(max_u, *this, dummy_patch);
    if(max_v < getMaxKnot_v()) this->uncheckedSplit_v(max_v, *this, dummy_patch);
  }

  /// \brief Private function to split patch geometry at a given u isoline,
  ///   without checking the trimming curves or normalizing the resulting knot vectors
  /// \sa NURBSPatch::split_u()
  void uncheckedSplit_u(T u, NURBSPatch& p1, NURBSPatch& p2) const
  {
    SLIC_ASSERT_MSG(
      isValidParameter_u(u, 1e-5),
      axom::fmt::format("Requested u-parameter {} for subdivision is outside valid range ({},{})",
                        u,
                        getMinKnot_u(),
                        getMaxKnot_u(),
                        1e-5));

    const bool isRationalPatch = isRational();

    const int p = getDegree_u();
    const int nq = getNumControlPoints_v() - 1;

    p1 = *this;
    p2.m_isTrimmed = m_isTrimmed;

    // Will make the multiplicity of the knot at u equal to p
    const auto k = p1.insertKnot_u(u, p);
    auto nkts1 = p1.getNumKnots_u();
    auto npts1 = p1.getNumControlPoints_u();

    // Split the knot vector, add to the returned curves
    KnotVectorType k1, k2;
    p1.getKnots_u().splitBySpan(k, k1, k2);

    p1.m_knotvec_u = k1;
    p1.m_knotvec_v = m_knotvec_v;

    p2.m_knotvec_u = k2;
    p2.m_knotvec_v = m_knotvec_v;

    // Copy the control points
    p2.m_controlPoints.resize(nkts1 - k - 1, nq + 1);
    if(isRationalPatch)
    {
      p2.m_weights.resize(nkts1 - k - 1, nq + 1);
    }
    else
    {
      p2.m_weights.resize(0, 0);
    }

    for(int i = 0; i < p2.m_controlPoints.shape()[0]; ++i)
    {
      for(int j = 0; j < p2.m_controlPoints.shape()[1]; ++j)
      {
        p2.m_controlPoints(nkts1 - k - 2 - i, j) = p1(npts1 - 1 - i, j);
        if(isRationalPatch)
        {
          p2.m_weights(nkts1 - k - 2 - i, j) = p1.getWeight(npts1 - 1 - i, j);
        }
      }
    }

    // Assumes that the resizing is done on the *flattened* array
    p1.m_controlPoints.resize(k - p + 1, nq + 1);
    if(isRationalPatch)
    {
      p1.m_weights.resize(k - p + 1, nq + 1);
    }
    else
    {
      p1.m_weights.resize(0, 0);
    }
  }

  /// \brief Private function to split patch geometry at a given v isoline,
  ///   without checking the trimming curves or normalizing the resulting knot vectors
  /// \sa NURBSPatch::split_v()
  void uncheckedSplit_v(T v, NURBSPatch& p1, NURBSPatch& p2) const
  {
    SLIC_ASSERT_MSG(
      isValidParameter_v(v, 1e-5),
      axom::fmt::format("Requested v-parameter {} for subdivision is outside valid range ({},{})",
                        v,
                        getMinKnot_v(),
                        getMaxKnot_v(),
                        1e-5));

    const bool isRationalPatch = isRational();

    const int np = getNumControlPoints_u() - 1;
    const int q = getDegree_v();

    p1 = *this;
    p2.m_isTrimmed = m_isTrimmed;

    // Will make the multiplicity of the knot at v equal to q
    const auto k = p1.insertKnot_v(v, q);
    auto nkts1 = p1.getNumKnots_v();
    auto npts1 = p1.getNumControlPoints_v();

    // Split the knot vector, add to the returned curves
    KnotVectorType k1, k2;
    p1.getKnots_v().splitBySpan(k, k1, k2);

    p1.m_knotvec_u = m_knotvec_u;
    p1.m_knotvec_v = k1;

    p2.m_knotvec_u = m_knotvec_u;
    p2.m_knotvec_v = k2;

    // Copy the control points
    p2.m_controlPoints.resize(np + 1, nkts1 - k - 1);
    if(isRationalPatch)
    {
      p2.m_weights.resize(np + 1, nkts1 - k - 1);
    }
    else
    {
      p2.m_weights.resize(0, 0);
    }

    for(int i = 0; i < p2.m_controlPoints.shape()[0]; ++i)
    {
      for(int j = 0; j < p2.m_controlPoints.shape()[1]; ++j)
      {
        p2.m_controlPoints(i, nkts1 - k - 2 - j) = p1(i, npts1 - 1 - j);
        if(isRationalPatch)
        {
          p2.m_weights(i, nkts1 - k - 2 - j) = p1.getWeight(i, npts1 - 1 - j);
        }
      }
    }

    // Rearrange the control points and weights by their flat index
    //  so that the `resize` method takes the correct submatrix
    for(int i = 0; i < np + 1; ++i)
    {
      for(int j = 0; j < k - q + 1; ++j)
      {
        p1.m_controlPoints.flatIndex(j + i * (k - q + 1)) = p1.m_controlPoints(i, j);
        if(isRationalPatch)
        {
          p1.m_weights.flatIndex(j + i * (k - q + 1)) = p1.m_weights(i, j);
        }
      }
    }

    // Resize the 2D arrays
    p1.m_controlPoints.resize(np + 1, k - q + 1);
    if(isRationalPatch)
    {
      p1.m_weights.resize(np + 1, k - q + 1);
    }
    else
    {
      p1.m_weights.resize(0, 0);
    }
  }
  /// \brief Private function to split a patch's trimming curves along a u/v isoline
  void splitTrimmingCurves(T uv,
                           bool splitInU,
                           TrimmingCurveVec& outCurvesFirst,
                           TrimmingCurveVec& outCurvesSecond) const
  {
    // Store a ray that is used as the splitting line
    primal::Ray<T, 2> ray_obj(Point<T, 2> {getMaxKnot_u() + 1.0, uv}, Vector<T, 2> {-1.0, 0.0});
    if(splitInU)
    {
      ray_obj = primal::Ray<T, 2>(Point<T, 2> {uv, getMinKnot_v() - 1.0}, Vector<T, 2> {0.0, 1.0});
    }
    TrimmingCurveVec split_trimming_curves;
    TrimmingCurveVec ray_trimming_curves;

    axom::Array<T> ray_params;
    for(const auto& curve : m_trimmingCurves)
    {
      axom::Array<T> curve_params;

      // Compute intersections for each subcurve individually to avoid
      //  a circular dependency with intersect.hpp
      {
        // Default parameters for intersection routine
        const double sq_tol = 1e-14;
        const double EPS = 1e-6;

        // Extract the Bezier curves of the NURBS curve, and check each for intersection
        axom::Array<T> knot_vals = curve.getKnots().getUniqueKnots();
        const auto beziers = curve.extractBezier();
        for(int i = 0; i < beziers.size(); ++i)
        {
          axom::Array<T> temp_curve_p;
          axom::Array<T> temp_ray_p;

          // Perform an initial check to see if the curve is linear and completely overlaps the ray
          if(beziers[i].isLinear(sq_tol) &&
             axom::utilities::isNearlyEqual(beziers[i][0][splitInU ? 0 : 1], uv) &&
             axom::utilities::isNearlyEqual(beziers[i][beziers[i].getOrder()][splitInU ? 0 : 1], uv))
          {
            continue;
          }

          detail::intersect_ray_bezier(ray_obj,
                                       beziers[i],
                                       temp_ray_p,
                                       temp_curve_p,
                                       sq_tol,
                                       EPS,
                                       beziers[i].getOrder(),
                                       0.,
                                       1.,
                                       false);

          // Scale the intersection parameters back into the span of the NURBS curve
          for(int j = 0; j < temp_curve_p.size(); ++j)
          {
            ray_params.push_back(temp_ray_p[j]);
            curve_params.push_back(knot_vals[i] + temp_curve_p[j] * (knot_vals[i + 1] - knot_vals[i]));
          }
        }
      }

      // Split all trimming curves at the intersection points
      if(curve_params.size() > 0)
      {
        // Sorting this keeps the splitting logic simpler
        std::sort(curve_params.begin(), curve_params.end());

        TrimmingCurveType c1, c2(curve);
        for(const auto& param : curve_params)
        {
          if(param <= c2.getMinKnot() || param >= c2.getMaxKnot())
          {
            continue;
          }

          c2.split(param, c1, c2);
          split_trimming_curves.push_back(c1);
        }
        split_trimming_curves.push_back(c2);
      }
      else
      {
        split_trimming_curves.push_back(curve);
      }
    }

    if(ray_params.size() != 0)
    {
      // Sort the ray parameters
      std::sort(ray_params.begin(), ray_params.end());

      for(int i = 0; i < ray_params.size() - 1; ++i)
      {
        // Skip any duplicate parameters
        if(ray_params[i + 1] - ray_params[i] < 1e-10)
        {
          continue;
        }

        // Determine if the ray segment is kept by the original surface
        ParameterPointType mid_ray_point(ray_obj.at(0.5 * (ray_params[i] + ray_params[i + 1])));
        bool isSegmentVisible = isVisible(mid_ray_point[0], mid_ray_point[1]);

        if(isSegmentVisible)
        {
          ray_trimming_curves.push_back(
            TrimmingCurveType::make_linear_segment_nurbs(ray_obj.at(ray_params[i]),
                                                         ray_obj.at(ray_params[i + 1])));
        }
      }
    }

    // Clear the output vectors
    outCurvesFirst.clear();
    outCurvesSecond.clear();

    // For all of the resulting trimming curves,
    //   add them to the right or left depending on the side of the ray
    for(auto& curve : split_trimming_curves)
    {
      auto eval_pt = curve.evaluate(0.5 * (curve.getMinKnot() + curve.getMaxKnot()));
      if(axom::utilities::isNearlyEqual(eval_pt[splitInU ? 0 : 1], uv, 1e-10))
      {
        // Can only happen if the curve is co-linear with the ray,
        //  decide what to do with it based on the orientation of the curve
        if(curve[0][splitInU ? 1 : 0] < curve[curve.getNumControlPoints() - 1][splitInU ? 1 : 0])
        {
          splitInU ? outCurvesFirst.push_back(curve) : outCurvesSecond.push_back(curve);
        }
        else
        {
          splitInU ? outCurvesSecond.push_back(curve) : outCurvesFirst.push_back(curve);
        }
      }
      else if(eval_pt[splitInU ? 0 : 1] < uv)
      {
        outCurvesFirst.push_back(curve);
      }
      else
      {
        outCurvesSecond.push_back(curve);
      }
    }

    for(auto& line : ray_trimming_curves)
    {
      outCurvesFirst.push_back(line);
      line.reverseOrientation();
      outCurvesSecond.push_back(line);
    }
  }

  /*!
   * \brief Equality operator for NURBS patches
   * 
   * \param [in] lhs The left-hand side NURBS patch
   * \param [in] rhs The right-hand side NURBS patch
   * 
   * \return True if the two patches are equal, false otherwise
   */
  friend inline bool operator==(const NURBSPatch<T, NDIMS>& lhs, const NURBSPatch<T, NDIMS>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) && (lhs.m_weights == rhs.m_weights) &&
      (lhs.m_knotvec_u == rhs.m_knotvec_u) && (lhs.m_knotvec_v == rhs.m_knotvec_v) &&
      (lhs.m_isTrimmed == rhs.m_isTrimmed) && (lhs.m_trimmingCurves == rhs.m_trimmingCurves);
  }

  /*!
   * \brief Inequality operator for NURBS patches
   * 
   * \param [in] lhs The left-hand side NURBS patch
   * \param [in] rhs The right-hand side NURBS patch
   * 
   * \return True if the two patches are not equal, false otherwise
   */
  friend inline bool operator!=(const NURBSPatch<T, NDIMS>& lhs, const NURBSPatch<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

private:
  CoordsMat m_controlPoints;
  WeightsMat m_weights;
  KnotVectorType m_knotvec_u, m_knotvec_v;

  bool m_isTrimmed {false};
  TrimmingCurveVec m_trimmingCurves;
};

//------------------------------------------------------------------------------
/// Free functions related to NURBSPatch
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSPatch<T, NDIMS>& nPatch)
{
  nPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::NURBSPatch using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::NURBSPatch<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_NURBSPATCH_HPP_
