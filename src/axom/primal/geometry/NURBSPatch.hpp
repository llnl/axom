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

/*! \brief Overloaded output operator for NURBS Patches*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSPatch<T, NDIMS>& nPatch);

/*!
 * \class NURBSPatch
 *
 * \brief Represents a 3D NURBS patch defined by a 2D array of control points
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
 * Untrimmed NURBS patches are identified by an internal flag, but a nonempty
 *  trimming curve vector should be marked as trimmed.
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

  AXOM_STATIC_ASSERT_MSG((NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
                         "A NURBS Patch object may be defined in 1-, 2-, or 3-D");

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A NURBS Patch must be defined using an arithmetic type");

public:
  /*! 
   * \brief Default constructor for an empty (invalid) NURBS patch
   *
   * \note An empty NURBS patch is not valid
   */
  NURBSPatch()
  {
    m_controlPoints.resize(0, 0);
    m_knotvec_u = KnotVectorType();
    m_knotvec_v = KnotVectorType();

    makeNonrational();
    makeUntrimmed();
  }

  /*!
   * \brief Constructor for a simple NURBS surface that reserves space for
   *  the minimum (sensible) number of points for the given degrees
   *
   * Constructs an empty patch by default (no nodes/weights on either axis)
   * 
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre deg_u, deg_v greater than or equal to 0.
   */
  NURBSPatch(int deg_u, int deg_v)
  {
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(deg_u + 1, deg_v + 1);
    m_knotvec_u = KnotVectorType(deg_u + 1, deg_u);
    m_knotvec_v = KnotVectorType(deg_v + 1, deg_v);

    makeNonrational();
    makeUntrimmed();
  }

  /*!
   * \brief Constructor for an empty NURBS surface from its size
   *
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * 
   * \pre Requires npts_d > deg_d and deg_d >= 0 for d = u, v 
   */
  NURBSPatch(int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(npts_u > deg_u && npts_v > deg_v);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);

    makeNonrational();
    makeUntrimmed();
  }

  /*!
   * \brief Constructor for a NURBS surface from a Bezier surface
   *
   * \param [in] bezierPatch the Bezier patch to convert to a NURBS patch 
   */
  explicit NURBSPatch(const BezierPatch<T, NDIMS>& bezierPatch)
  {
    m_controlPoints = bezierPatch.getControlPoints();
    m_weights = bezierPatch.getWeights();

    int deg_u = bezierPatch.getOrder_u();
    int deg_v = bezierPatch.getOrder_v();

    m_knotvec_u = KnotVectorType(deg_u + 1, deg_u);
    m_knotvec_v = KnotVectorType(deg_v + 1, deg_v);

    makeUntrimmed();
  }

  /*!
   * \brief Constructor for a NURBS Patch from an array of coordinates and degrees
   *
   * \param [in] pts A 1D C-style array of npts_u*npts_v control points
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   * 
   * Elements of pts[k] are mapped to control nodes (p, q) lexicographically, i.e.
   * pts[k] = nodes[ k // (npts_u + 1), k % npts_v ]
   */
  NURBSPatch(const PointType* pts, int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(npts_u >= deg_u + 1 && npts_v >= deg_v + 1);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    makeNonrational();
    makeUntrimmed();

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D C-style array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   * 
   * Elements of pts[k] are mapped to control nodes (p, q) lexicographically, i.e.
   * pts[k] = nodes[ k // (npts_u + 1), k % npts_v ]
   */
  NURBSPatch(const PointType* pts, const T* weights, int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(weights != nullptr);
    SLIC_ASSERT(npts_u >= deg_u + 1 && npts_v >= deg_v + 1);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    m_weights.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_weights.flatIndex(t) = weights[t];
    }

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);

    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 1D arrays of coordinates and degrees
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   * 
   * Elements of pts[k] are mapped to control nodes (p, q) lexicographically, i.e.
   * pts[k] = nodes[ k // (npts_u + 1), k % npts_v ]
   */
  NURBSPatch(const CoordsVec& pts, int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(npts_u >= deg_u + 1 && npts_v >= deg_v + 1);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < pts.size(); ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    makeNonrational();
    makeUntrimmed();

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 1D arrays of coordinates and weights
   *
   * \param [in] pts A 1D axom::Array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D axom::Array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   * 
   * Elements of pts[k] are mapped to control nodes (p, q) lexicographically, i.e.
   * pts[k] = nodes[ k // (npts_u + 1), k % npts_v ]
   */
  NURBSPatch(const CoordsVec& pts, const WeightsVec& weights, int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(npts_u > deg_u && npts_v > deg_v);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < pts.size(); ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    m_weights.resize(npts_u, npts_v);
    for(int t = 0; t < weights.size(); ++t)
    {
      m_weights.flatIndex(t) = weights[t];
    }

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);

    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 2D arrays of coordinates and degrees
   *
   * \param [in] pts A 2D axom::Array of (npts_u, npts_v) control points
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   */
  NURBSPatch(const CoordsMat& pts, int deg_u, int deg_v) : m_controlPoints(pts)
  {
    const auto pts_shape = pts.shape();

    SLIC_ASSERT(pts_shape[0] >= deg_u + 1 && pts_shape[1] >= deg_v + 1);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    makeNonrational();
    makeUntrimmed();

    m_knotvec_u = KnotVectorType(pts_shape[0], deg_u);
    m_knotvec_v = KnotVectorType(pts_shape[1], deg_v);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 2D arrays of coordinates and weights
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D axom::Array of (ord_u+1, ord_v+1) positive weights
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   */
  NURBSPatch(const CoordsMat& pts, const WeightsMat& weights, int deg_u, int deg_v)
    : m_controlPoints(pts)
    , m_weights(weights)
  {
    const auto pts_shape = pts.shape();

    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);
    SLIC_ASSERT(pts_shape[0] >= deg_u + 1 && pts_shape[1] >= deg_v + 1);
    SLIC_ASSERT(pts_shape == weights.shape());

    m_knotvec_u = KnotVectorType(pts_shape[0], deg_u);
    m_knotvec_v = KnotVectorType(pts_shape[1], deg_v);

    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from C-style arrays of coordinates and knot vectors
   *
   * \param [in] pts A 1D C-style array of npts_u*npts_v control points
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] knots_u A 1D C-style array of npts_u + deg_u + 1 knots
   * \param [in] nkts_u The number of knots in the u direction
   * \param [in] knots_v A 1D C-style array of npts_v + deg_v + 1 knots
   * \param [in] nkts_v The number of knots in the v direction
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires valid pointers and knot vectors
   */
  NURBSPatch(const PointType* pts,
             int npts_u,
             int npts_v,
             const T* knots_u,
             int nkts_u,
             const T* knots_v,
             int nkts_v)
  {
    SLIC_ASSERT(pts != nullptr && knots_u != nullptr && knots_v != nullptr);
    SLIC_ASSERT(npts_u >= 0 && npts_v >= 0);
    SLIC_ASSERT(nkts_u >= 0 && nkts_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    makeNonrational();
    makeUntrimmed();

    m_knotvec_u = KnotVectorType(knots_u, nkts_u, nkts_u - npts_u - 1);
    m_knotvec_v = KnotVectorType(knots_v, nkts_v, nkts_v - npts_v - 1);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from C-style arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of npts_u*npts_v control points
   * \param [in] weights A 1D C-style array of npts_u*npts_v positive weights
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] knots_u A 1D C-style array of npts_u + deg_u + 1 knots
   * \param [in] nkts_u The number of knots in the u direction
   * \param [in] knots_v A 1D C-style array of npts_v + deg_v + 1 knots
   * \param [in] nkts_v The number of knots in the v direction
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
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
  {
    SLIC_ASSERT(pts != nullptr && weights != nullptr && knots_u != nullptr && knots_v != nullptr);
    SLIC_ASSERT(npts_u >= 0 && npts_v >= 0);
    SLIC_ASSERT(nkts_u >= 0 && nkts_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    m_weights.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_weights.flatIndex(t) = weights[t];
    }

    m_knotvec_u = KnotVectorType(knots_u, nkts_u, nkts_u - npts_u - 1);
    m_knotvec_v = KnotVectorType(knots_v, nkts_v, nkts_v - npts_v - 1);

    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates and knots
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             int npts_u,
             int npts_v,
             const axom::Array<T>& knots_u,
             const axom::Array<T>& knots_v)
  {
    SLIC_ASSERT(npts_u >= 0 && npts_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < pts.size(); ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    makeNonrational();
    makeUntrimmed();

    m_knotvec_u = KnotVectorType(knots_u, knots_u.size() - npts_u - 1);
    m_knotvec_v = KnotVectorType(knots_v, knots_v.size() - npts_v - 1);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates, weights, and knots
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] weights A 1D axom::Array of npts_u*npts_v positive weights
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             const WeightsVec& weights,
             int npts_u,
             int npts_v,
             const axom::Array<T>& knots_u,
             const axom::Array<T>& knots_v)
  {
    SLIC_ASSERT(npts_u >= 0 && npts_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < pts.size(); ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    m_weights.resize(npts_u, npts_v);
    for(int t = 0; t < weights.size(); ++t)
    {
      m_weights.flatIndex(t) = weights[t];
    }

    m_knotvec_u = KnotVectorType(knots_u, knots_u.size() - npts_u - 1);
    m_knotvec_v = KnotVectorType(knots_v, knots_v.size() - npts_v - 1);

    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates and KnotVectors
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] knotvec_u An KnotVector object for the first axis
   * \param [in] knotvec_v An KnotVector object for the second axis
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             int npts_u,
             int npts_v,
             const KnotVectorType& knotvec_u,
             const KnotVectorType& knotvec_v)
    : m_knotvec_u(knotvec_u)
    , m_knotvec_v(knotvec_v)
  {
    SLIC_ASSERT(npts_u >= 0 && npts_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < pts.size(); ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    makeNonrational();
    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array arrays of coordinates, weights, and KnotVectors
   *
   * \param [in] pts A 1D axom::Array of npts_u*npts_v control points
   * \param [in] weights A 1D axom::Array of npts_u*npts_v positive weights
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] knotvec_u An KnotVector object for the first axis
   * \param [in] knotvec_v An KnotVector object for the second axis
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsVec& pts,
             const WeightsVec& weights,
             int npts_u,
             int npts_v,
             const KnotVectorType& knotvec_u,
             const KnotVectorType& knotvec_v)
    : m_knotvec_u(knotvec_u)
    , m_knotvec_v(knotvec_v)
  {
    SLIC_ASSERT(npts_u >= 0 && npts_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < pts.size(); ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    m_weights.resize(npts_u, npts_v);
    for(int t = 0; t < weights.size(); ++t)
    {
      m_weights.flatIndex(t) = weights[t];
    }

    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 2D axom::Array array of coordinates and array of knots
   *
   * \param [in] pts A 2D axom::Array of (npts_u, npts_v) control points
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts, const axom::Array<T>& knots_u, const axom::Array<T>& knots_v)
    : m_controlPoints(pts)
  {
    auto pts_shape = pts.shape();

    SLIC_ASSERT(pts_shape[0] >= 0 && pts_shape[1] >= 0);

    makeNonrational();
    makeUntrimmed();

    m_knotvec_u = KnotVectorType(knots_u, knots_u.size() - pts_shape[0] - 1);
    m_knotvec_v = KnotVectorType(knots_v, knots_v.size() - pts_shape[1] - 1);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 2D axom::Array array of coordinates, weights, and array of knots
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D axom::Array of (ord_u+1, ord_v+1) positive weights
   * \param [in] knots_u An axom::Array of npts_u + deg_u + 1 knots
   * \param [in] knots_v An axom::Array of npts_v + deg_v + 1 knots
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts,
             const WeightsMat& weights,
             const axom::Array<T>& knots_u,
             const axom::Array<T>& knots_v)
    : m_controlPoints(pts)
    , m_weights(weights)
  {
    auto pts_shape = pts.shape();

    SLIC_ASSERT(pts_shape[0] >= 0 && pts_shape[1] >= 0);

    m_knotvec_u = KnotVectorType(knots_u, knots_u.size() - pts_shape[0] - 1);
    m_knotvec_v = KnotVectorType(knots_v, knots_v.size() - pts_shape[1] - 1);

    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 1D axom::Array array of coordinates and KnotVector objects
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] knotvec_u A KnotVector object for the first axis
   * \param [in] knotvec_v A KnotVector object for the second axis
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts, const KnotVectorType& knotvec_u, const KnotVectorType& knotvec_v)
    : m_controlPoints(pts)
    , m_knotvec_u(knotvec_u)
    , m_knotvec_v(knotvec_v)
  {
    makeNonrational();
    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Patch from 2D axom::Array array of coordinates, weights, and KnotVector objects
   *
   * \param [in] pts A 2D axom::Array of (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D axom::Array of (ord_u+1, ord_v+1) positive weights
   * \param [in] knotvec_u A KnotVector object for the first axis
   * \param [in] knotvec_v A KnotVector object for the second axis
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts_d > deg_d
   */
  NURBSPatch(const CoordsMat& pts,
             const WeightsMat& weights,
             const KnotVectorType& knotvec_u,
             const KnotVectorType& knotvec_v)
    : m_controlPoints(pts)
    , m_weights(weights)
    , m_knotvec_u(knotvec_u)
    , m_knotvec_v(knotvec_v)
  {
    makeUntrimmed();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Evaluate the untrimmed NURBS surface at a particular parameter value \a t
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
    SLIC_ASSERT(m_knotvec_u.isValidParameter(u));
    SLIC_ASSERT(m_knotvec_v.isValidParameter(v));

    u = axom::utilities::clampVal(u, m_knotvec_u[0], m_knotvec_u[m_knotvec_u.getNumKnots() - 1]);
    v = axom::utilities::clampVal(v, m_knotvec_v[0], m_knotvec_v[m_knotvec_v.getNumKnots() - 1]);

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
   * \brief Reset the degree and resize arrays of points (and weights)
   * 
   * \param [in] npts_u The target number of control points on the first axis
   * \param [in] npts_v The target number of control points on the second axis
   * \param [in] deg_u The target degree on the first axis
   * \param [in] deg_v The target degree on the second axis
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

  /// \brief Returns the degree of the NURBS Patch on the first axis
  int getDegree_u() const { return m_knotvec_u.getDegree(); }

  /// \brief Returns the degree of the NURBS Patch on the second axis
  int getDegree_v() const { return m_knotvec_v.getDegree(); }

  /// \brief Returns the order (degree + 1) of the NURBS Patch on the first axis
  int getOrder_u() const { return m_knotvec_u.getDegree() + 1; }

  /// \brief Returns the order of the NURBS Patch on the second axis
  int getOrder_v() const { return m_knotvec_v.getDegree() + 1; }

  /// \brief Return a copy of the KnotVector instance on the first axis
  KnotVectorType getKnots_u() const { return m_knotvec_u; }

  /// \brief Return an array of knot values on the first axis
  axom::Array<T> getKnotsArray_u() const { return m_knotvec_u.getArray(); }

  /// \brief Return a copy of the KnotVector instance on the second axis
  KnotVectorType getKnots_v() const { return m_knotvec_v; }

  /// \brief Return an array of knot values on the second axis
  axom::Array<T> getKnotsArray_v() const { return m_knotvec_v.getArray(); }

  /// \brief Returns the number of control points in the NURBS Patch on the first axis
  int getNumControlPoints_u() const { return static_cast<int>(m_controlPoints.shape()[0]); }

  /// \brief Returns the number of control points in the NURBS Patch on the second axis
  int getNumControlPoints_v() const { return static_cast<int>(m_controlPoints.shape()[1]); }

  /// \brief Return the length of the knot vector on the first axis
  int getNumKnots_u() const { return m_knotvec_u.getNumKnots(); }

  /// \brief Return the length of the knot vector on the second axis
  int getNumKnots_v() const { return m_knotvec_v.getNumKnots(); }

  /// Make trivially rational. If already rational, do nothing
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

  /// \brief Use array size as flag for rationality
  bool isRational() const { return !m_weights.empty(); }

  /// \brief Get array of trimming curvse
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
    m_trimmingCurves.insert(m_trimmingCurves.end(), curves.begin(), curves.end());
  }

  /// \brief Clear trimming curves, but DON'T mark as untrimmed
  void clearTrimmingCurves() { m_trimmingCurves.clear(); }

  /// \brief Get number of trimming curves
  int getNumTrimmingCurves() const { return m_trimmingCurves.size(); }

  /// \brief use array size as flag for trimmed-ness
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
    TrimmingCurveType curve;
    constexpr int num_points = 2;
    constexpr int degree = 1;
    curve.setParameters(num_points, degree);

    // Top
    curve[0] = ParameterPointType({max_u, max_v});
    curve[1] = ParameterPointType({min_u, max_v});
    addTrimmingCurve(curve);

    // Right
    curve[0] = ParameterPointType({max_u, min_v});
    curve[1] = ParameterPointType({max_u, max_v});
    addTrimmingCurve(curve);

    // Bottom
    curve[0] = ParameterPointType({min_u, min_v});
    curve[1] = ParameterPointType({max_u, min_v});
    addTrimmingCurve(curve);

    // Left
    curve[0] = ParameterPointType({min_u, max_v});
    curve[1] = ParameterPointType({min_u, min_v});
    addTrimmingCurve(curve);

    m_isTrimmed = true;
  }

  /*!
   * \brief Check if a parameter point is visible on the NURBS patch
   *
   * \param [in] uv The parameter point to check
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
      gwn += detail::nurbs_winding_number(uv, curve);
    }

    return std::lround(gwn) % 2 != 0;
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

  /// Retrieves the control point at index \a (idx_p, idx_q)
  PointType& operator()(int ui, int vi) { return m_controlPoints(ui, vi); }

  /// Retrieves the vector of control points at index \a idx
  const PointType& operator()(int ui, int vi) const { return m_controlPoints(ui, vi); }

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

  /// Returns a copy of the NURBS patch's control points
  CoordsMat getControlPoints() const { return m_controlPoints; }

  /// Returns a copy of the NURBS patch's weights
  WeightsMat getWeights() const { return m_weights; }

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
    }
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
   * \brief Returns an isocurve with a fixed value of u
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
    SLIC_ASSERT(m_knotvec_u.isValidParameter(u));
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
   * \brief Returns an isocurve with a fixed value of v
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
    SLIC_ASSERT(m_knotvec_v.isValidParameter(v));
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
   * \brief Evaluate the untrimmed surface and the first \a d derivatives at parameter \a u, \a v
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
    SLIC_ASSERT(m_knotvec_u.isValidParameter(u));
    SLIC_ASSERT(m_knotvec_v.isValidParameter(v));

    u = axom::utilities::clampVal(u, m_knotvec_u[0], m_knotvec_u[m_knotvec_u.getNumKnots() - 1]);
    v = axom::utilities::clampVal(v, m_knotvec_v[0], m_knotvec_v[m_knotvec_v.getNumKnots() - 1]);

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
   * \brief Evaluates all first derivatives of the untrimmed patch at (\a u, \a v)
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
   * \brief Evaluates all linear derivatives of the untrimmed patch at (\a u, \a v)
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
   * \brief Evaluates all second derivatives of the untrimmed patch at (\a u, \a v)
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
   * \brief Computes a tangent in u of the untrimmed patch at (\a u, \a v)
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
   * \brief Computes a tangent in v of the NURBS patch at (\a u, \a v)
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
   * \brief Computes the second derivative in u of the untrimmed patch at (\a u, \a v)
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
   * \brief Computes the second derivative in v of the untrimmed patch at (\a u, \a v)
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
   * \brief Computes the mixed second derivative in u and v of the untrimmed patch at (\a u, \a v)
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
   * \brief Computes the mixed second derivative in u and v of the untrimmed patch at (\a u, \a v)
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
   * \brief Computes the normal vector to the untrimmed patch at (\a u, \a v)
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
    SLIC_ASSERT(m_knotvec_u.isValidParameter(u));
    u = axom::utilities::clampVal(u, m_knotvec_u[0], m_knotvec_u[m_knotvec_u.getNumKnots() - 1]);

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
    SLIC_ASSERT(m_knotvec_v.isValidParameter(v));
    v = axom::utilities::clampVal(v, m_knotvec_v[0], m_knotvec_v[m_knotvec_v.getNumKnots() - 1]);

    SLIC_ASSERT(target_multiplicity > 0);

    const bool isRationalPatch = isRational();

    const int np = getNumControlPoints_u() - 1;
    const int p = getDegree_u();

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
    axom::Array<T, 2> alpha(p - s, r + 1);
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

  /*!
    * \brief Splits an untrimmed NURBS patch into four NURBS patches
    *
    * \param [in] u parameter value at which to bisect on the first axis
    * \param [in] v parameter value at which to bisect on the second axis
    * \param [out] p1 First output NURBS patch
    * \param [out] p2 Second output NURBS patch
    * \param [out] p3 Third output NURBS patch
    * \param [out] p4 Fourth output NURBS patch
    *
    *   v = 1
    *   ----------------------
    *   |         |          |
    *   |   p3    |    p4    |
    *   |         |          |
    *   --------(u,v)---------
    *   |         |          |
    *   |   p1    |    p2    |
    *   |         |          |
    *   ---------------------- u = 1
    *
    * \pre Parameter \a u and \a v must be *strictly interior* to the knot span
    * \pre The patch must be untrimmed
    */
  void split(T u, T v, NURBSPatch& p1, NURBSPatch& p2, NURBSPatch& p3, NURBSPatch& p4) const
  {
    SLIC_ASSERT(m_knotvec_u.isValidInteriorParameter(u));
    SLIC_ASSERT(m_knotvec_v.isValidInteriorParameter(v));
    SLIC_ASSERT_MSG(!isTrimmed(), "Splitting a trimmed patch is not yet supported");

    // Bisect the patch along the u direction
    split_u(u, p1, p2);

    // Temporarily store the result in each half and split again
    NURBSPatch p0(p1);
    p0.split_v(v, p1, p2);

    p0 = p2;
    p0.split_v(v, p3, p4);
  }

  /*!
   * \brief Split the untrimmed NURBS patch in two along the u direction
   *
   * \pre The patch must be untrimmed
   */
  void split_u(T u, NURBSPatch& p1, NURBSPatch& p2, bool normalize = false) const
  {
    SLIC_ASSERT(m_knotvec_u.isValidInteriorParameter(u));
    SLIC_ASSERT_MSG(!isTrimmed(), "Splitting a trimmed patch is not yet supported");

    const bool isRationalPatch = isRational();

    const int p = getDegree_u();
    const int nq = getNumControlPoints_v() - 1;

    p1 = *this;

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

    if(normalize)
    {
      p1.normalize();
      p2.normalize();
    }
  }

  /*!
   * \brief Split the untrimmed NURBS patch in two along the v direction
   *
   * \pre The patch must be untrimmed
   */
  void split_v(T v, NURBSPatch& p1, NURBSPatch& p2, bool normalize = false) const
  {
    SLIC_ASSERT(m_knotvec_v.isValidInteriorParameter(v));
    SLIC_ASSERT_MSG(!isTrimmed(), "Splitting a trimmed patch is not yet supported");

    const bool isRationalPatch = isRational();

    const int np = getNumControlPoints_u() - 1;
    const int q = getDegree_v();

    p1 = *this;

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

    if(normalize)
    {
      p1.normalize();
      p2.normalize();
    }
  }

  /*!
   * \brief Splits the untrimmed NURBS surface (at each internal knot) into several Bezier patches
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

  /// \brief Normalize the knot vectors to the span [0, 1]
  void normalize()
  {
    m_knotvec_u.normalize();
    m_knotvec_v.normalize();

    rescaleTrimmingCurves_u(0.0, 1.0);
    rescaleTrimmingCurves_v(0.0, 1.0);
  }

  /// \brief Normalize the knot vector in u to the span [0, 1]
  void normalize_u()
  {
    m_knotvec_u.normalize();
    rescaleTrimmingCurves_u(0.0, 1.0);
  }

  /// \brief Normalize the knot vector in v to the span [0, 1]
  void normalize_v()
  {
    m_knotvec_v.normalize();
    rescaleTrimmingCurves_v(0.0, 1.0);
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
    m_knotvec_u.rescale(a, b);
    m_knotvec_v.rescale(a, b);

    rescaleTrimmingCurves_u(a, b);
    rescaleTrimmingCurves_v(a, b);
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
    m_knotvec_u.rescale(a, b);

    rescaleTrimmingCurves_u(a, b);
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
    m_knotvec_v.rescale(a, b);

    rescaleTrimmingCurves_v(a, b);
  }

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
  bool isValidInteriorParameter(T t) const { return m_knotvec_u.isValidInteriorParameter(t); }

private:
  CoordsMat m_controlPoints;
  WeightsMat m_weights;
  KnotVectorType m_knotvec_u, m_knotvec_v;

  bool m_isTrimmed;
  TrimmingCurveVec m_trimmingCurves;

  void rescaleTrimmingCurves_u(T a, T b)
  {
    for(auto& curve : m_trimmingCurves)
    {
      for(int i = 0; i < curve.getNumControlPoints(); ++i)
      {
        curve[i][0] = a + (b - a) * curve[i][0];
      }
    }
  }

  void rescaleTrimmingCurves_v(T a, T b)
  {
    for(auto& curve : m_trimmingCurves)
    {
      for(int i = 0; i < curve.getNumControlPoints(); ++i)
      {
        curve[i][1] = a + (b - a) * curve[i][1];
      }
    }
  }
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
