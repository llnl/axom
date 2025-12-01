// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file BezierPatch.hpp
 *
 * \brief A BezierPatch primitive
 */

#ifndef AXOM_PRIMAL_BEZIERPATCH_HPP_
#define AXOM_PRIMAL_BEZIERPATCH_HPP_

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
class BezierPatch;

/*! \brief Overloaded output operator for Bezier Patches*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierPatch<T, NDIMS>& bPatch);

/*!
 * \class BezierPatch
 *
 * \brief Represents a 3D Bezier patch defined by a 2D array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * The order of a Bezier patch with (N+1)(M+1) control points is (N, M).
 * The patch is approximated by the control points,
 * parametrized from u=0 to u=1 and v=0 to v=1.
 * 
 * Contains a 2D array of positive weights to represent a rational Bezier patch.
 * Polynomial (nonrational) Bezier patches are identified by an empty weights array.
 * 
 * Algorithms for Rational Bezier curves derived from 
 * Gerald Farin, "Algorithms for rational Bezier curves"
 * Computer-Aided Design, Volume 15, Number 2, 1983,
 */
template <typename T, int NDIMS>
class BezierPatch
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using PlaneType = Plane<T, NDIMS>;

  using CoordsVec = axom::Array<PointType, 1>;
  using CoordsMat = axom::Array<PointType, 2>;
  using WeightsVec = axom::Array<T, 1>;
  using WeightsMat = axom::Array<T, 2>;

  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;
  using BezierCurveType = primal::BezierCurve<T, NDIMS>;

  AXOM_STATIC_ASSERT_MSG((NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
                         "A Bezier Patch object may be defined in 1-, 2-, or 3-D");

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A Bezier Patch must be defined using an arithmetic type");

public:
  ///@{
  /** 
   * @name Constructors for BezierPatch
   *
   * The constructors allow for flexible initialization of BezierPatch objects from:
   * - 1D or 2D Axom arrays of control points and weights,
   * - C-style arrays of control points and weights,
   * - Specified polynomial orders in each axis,
   * - Rational or polynomial (nonrational) patches, depending on the presence of weights.
   *
   * For 1D arrays, the mapping of control points and weights to the patch is lexicographical, i.e.
     \verbatim
      pts[0]               -> nodes[0, 0],     ..., pts[ord_v]       -> nodes[0, ord_v]
      pts[ord_v+1]         -> nodes[1, 0],     ..., pts[2*ord_v]     -> nodes[1, ord_v]
                                               ...
      pts[ord_u*(ord_v-1)] -> nodes[ord_u, 0], ..., pts[ord_u*ord_v] -> nodes[ord_u, ord_v]
     \endverbatim
   * 
   * The patch is parametrized from u=0 to u=1 and v=0 to v=1. 
   * 
   * Rational patches are identified by a non-empty weights array, 
   * and nonrational patches by an empty weights array.
   * All weights must be greater than 0 in a rational patch.
   */

  /**
   * \brief Constructor from an ArrayView over the control points and weights
   * 
   * \param [in] controlPoints ArrayView of control points (size: (ord_u+1, ord_v+1) or (0,0))
   * \param [in] weights ArrayView of weights (size: (ord_u+1, ord_v+1) or (0,0))
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * 
   * If \a controlPoints is empty, we still allocate space for (ord_u+1, ord_v+1) control points
   * \pre ord_u and ord_v must either both be at least 0 for a valid patch, or both must be -1 for an empty patch
   */
  BezierPatch(axom::ArrayView<const PointType, 2> controlPoints,
              axom::ArrayView<const T, 2> weights,
              int ord_u,
              int ord_v)
  {
    SLIC_ASSERT((ord_u == -1 && ord_v == -1) || (ord_u >= 0 && ord_u >= 0));
    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    SLIC_ASSERT(controlPoints.size() >= weights.size());

    // note: always allocate space for control points
    if(controlPoints.empty())
    {
      m_controlPoints.resize(sz_u, sz_v);
    }
    else
    {
      SLIC_ASSERT(controlPoints.data() != nullptr);
      SLIC_ASSERT(controlPoints.shape()[0] == sz_u && controlPoints.shape()[1] == sz_v);
      m_controlPoints = controlPoints;
    }

    // note: only allocate space for weights when they are supplied
    if(!weights.empty())
    {
      SLIC_ASSERT(weights.data() != nullptr);
      SLIC_ASSERT(weights.shape()[0] == sz_u && weights.shape()[1] == sz_v);
      m_weights = weights;
      SLIC_ASSERT(isValidRational());
    }
  }

  /// Constructor for polynomial BezierPatch from ArrayView of control points
  BezierPatch(axom::ArrayView<PointType, 2> controlPoints, int ord_u, int ord_v)
    : BezierPatch(axom::ArrayView<const PointType, 2>(controlPoints.data(), controlPoints.shape()),
                  axom::ArrayView<const T, 2>(nullptr, {0, 0}),
                  ord_u,
                  ord_v)
  { }

  /// Constructor for rational BezierPatch from ArrayView of control points
  BezierPatch(axom::ArrayView<PointType, 2> controlPoints,
              axom::ArrayView<T, 2> weights,
              int ord_u,
              int ord_v)
    : BezierPatch(axom::ArrayView<const PointType, 2>(controlPoints.data(), controlPoints.shape()),
                  axom::ArrayView<const T, 2>(weights.data(), weights.shape()),
                  ord_u,
                  ord_v)
  { }

  /*!
   * \brief Constructor for a polynomial (nonrational) Bezier Patch that reserves 
   *  space for the given order of the surface
   * 
   * \param [in] ord_u, ord_v The patch's polynomial orders
   * \pre ord_u, ord_v greater than or equal to -1.
   */
  BezierPatch(int ord_u = -1, int ord_v = -1)
    : BezierPatch(axom::ArrayView<const PointType, 2>(nullptr, {0, 0}),
                  axom::ArrayView<const T, 2>(nullptr, {0, 0}),
                  ord_u,
                  ord_v)
  { }

  /*!
   * \brief Constructor for a Bezier Patch from an array of coordinates
   *
   * \param [in] pts A 1D C-style array of (ord_u+1)*(ord_v+1) control points
   * \param [in] ord_u, ord_v The patch's polynomial orders
   * \pre order in both directions is greater than or equal to zero
   */
  BezierPatch(const PointType* pts, int ord_u, int ord_v)
    : BezierPatch(axom::ArrayView<const PointType, 2>(pts, {ord_u + 1, ord_v + 1}),
                  axom::ArrayView<const T, 2>(nullptr, {0, 0}),
                  ord_u,
                  ord_v)
  { }

  /*!
   * \brief Constructor for a rational Bezier Patch from arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D C-style array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] ord_u, ord_v The patch's polynomial orders, both greater than or equal to zero
   *  
   * If \a weights is the null pointer, creates a nonrational curve
   */
  BezierPatch(PointType* pts, T* weights, int ord_u, int ord_v)
    : BezierPatch(axom::ArrayView<const PointType, 2>(pts, {ord_u + 1, ord_v + 1}),
                  axom::ArrayView<const T, 2>(weights, {ord_u + 1, ord_v + 1}),
                  ord_u,
                  ord_v)
  { }

  /*!
   * \brief Constructor for a Bezier Patch from a 1D Axom array of coordinates
   *
   * \param [in] pts A 1D Axom array of (ord_u+1)*(ord_v+1) control points
   * \param [in] ord_u, ord_v The patch's polynomial orders, both greater than or equal to zero
   */
  BezierPatch(const CoordsVec& pts, int ord_u, int ord_v)
    : BezierPatch(axom::ArrayView<const PointType, 2>(pts.data(), {ord_u + 1, ord_v + 1}),
                  axom::ArrayView<const T, 2>(nullptr, {0, 0}),
                  ord_u,
                  ord_v)
  { }

  /*!
   * \brief Constructor for a Rational Bezier Patch from 1D Axom arrays of coordinates and weights
   *
   * \param [in] pts A 1D Axom array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D Axom array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] ord_u, ord_v The patch's polynomial orders, both greater than or equal to zero
   */
  BezierPatch(const CoordsVec& pts, const WeightsVec& weights, int ord_u, int ord_v)
    : BezierPatch(axom::ArrayView<const PointType, 2>(pts.data(), {ord_u + 1, ord_v + 1}),
                  axom::ArrayView<const T, 2>(weights.data(), {ord_u + 1, ord_v + 1}),
                  ord_u,
                  ord_v)
  { }

  /*!
   * \brief Constructor for a Bezier Patch from an Axom array of coordinates
   *
   * \param [in] pts A 2D Axom array with (ord_u+1, ord_v+1) control points
   * \param [in] ord_u, ord_v The patch's polynomial orders, both greater than or equal to zero
   */
  BezierPatch(const CoordsMat& pts, int ord_u, int ord_v)
    : BezierPatch(pts.view(), axom::ArrayView<const T, 2>(nullptr, {0, 0}), ord_u, ord_v)
  { }

  /*!
   * \brief Constructor for a rational Bezier Patch from a 2D Axom array  of weights and coordinates
   *
   * \param [in] pts A 2D Axom array with (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D Axom array with (ord_u+1, ord_v+1) weights
   * \param [in] ord_u, ord_v The patch's polynomial orders, both greater than or equal to zero
   */
  BezierPatch(const CoordsMat& pts, const WeightsMat& weights, int ord_u, int ord_v)
    : BezierPatch(pts.view(), weights.view(), ord_u, ord_v)
  { }

  ///@}

  ///@{
  /// \name Query/modify patch properties (order, rationality, ...)

  /*!
   * \brief Sets the order of Bezier patch
   *
   * \param [in] ord_u, ord_v The patch's polynomial orders, 
   * 
   * \pre ord_u and ord_v must both be -1 or both greater than or equal to zero
   *
   * \note Will only resize the arrays, and likely make the patch invalid
   */
  void setOrder(int ord_u, int ord_v)
  {
    SLIC_ASSERT((ord_u == -1 && ord_v == -1) || (ord_u >= 0 && ord_v >= 0));
    const int SZ_U = utilities::max(0, ord_u + 1);
    const int SZ_V = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(SZ_U, SZ_V);
    if(isRational())
    {
      m_weights.resize(SZ_U, SZ_V);
    }
  }

  /// Returns the order of the Bezier Patch on the first axis
  int getOrder_u() const { return static_cast<int>(m_controlPoints.shape()[0]) - 1; }

  /// Returns the order of the Bezier Patch on the second axis
  int getOrder_v() const { return static_cast<int>(m_controlPoints.shape()[1]) - 1; }

  /// Clears the list of control points, make nonrational
  void clear()
  {
    m_controlPoints.clear();
    makeNonrational();
  }

  /// Use array size as flag for rationality
  bool isRational() const { return !m_weights.empty(); }

  /// Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      const int ord_u = getOrder_u();
      const int ord_v = getOrder_v();

      m_weights.resize(ord_u + 1, ord_v + 1);
      m_weights.fill(1.0);
    }
  }

  /// Make nonrational by shrinking array of weights
  void makeNonrational() { m_weights.clear(); }

  ///@}

  ///@{
  /// \name Query/modify patch's geometry (control points, weights, bounding box, ...)

  /// Retrieves the control point at index \a (idx_p, idx_q)
  PointType& operator()(int ui, int vi) { return m_controlPoints(ui, vi); }

  /// Retrieves the vector of control points at index \a idx
  const PointType& operator()(int ui, int vi) const { return m_controlPoints(ui, vi); }

  /// Returns a reference to the Bezier patch's control points
  CoordsMat& getControlPoints() { return m_controlPoints; }

  /// Returns a const reference to the Bezier patch's control points
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
  };

  /// Returns a reference to the Bezier patch's weights
  WeightsMat& getWeights() { return m_weights; }

  /// Returns a const reference to the Bezier patch's weights
  const WeightsMat& getWeights() const { return m_weights; }

  /// Returns an axis-aligned bounding box containing the Bezier patch
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  /// Returns an oriented bounding box containing the Bezier patch
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(), static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Predicate to check if the Bezier patch is approximately planar (i.e. flat)
   *
   * This function checks if all control points of the BezierPatch
   * are approximately on the plane defined by its four corners
   *
   * \param [in] sq_tol Threshold for sum of squared distances
   * \param [in] EPS Threshold for nearness to zero
   * \return True if the object is planar up to tolerance \a sq_tol
   */
  bool isPlanar(double sq_tol = 1e-8, double EPS = 1e-8) const
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    if(ord_u <= 0 && ord_v <= 0)
    {
      return true;
    }
    if(ord_u == 1 && ord_v == 0)
    {
      return true;
    }
    if(ord_u == 0 && ord_v == 1)
    {
      return true;
    }

    // Check that the four corners aren't coplanar
    VectorType v1(m_controlPoints(0, 0), m_controlPoints(0, ord_v));
    VectorType v2(m_controlPoints(0, 0), m_controlPoints(ord_u, 0));
    VectorType v3(m_controlPoints(0, 0), m_controlPoints(ord_u, ord_v));
    if(!axom::utilities::isNearlyEqual(VectorType::scalar_triple_product(v1, v2, v3), 0.0, EPS))
    {
      return false;
    }

    // Find three points that produce a nonzero normal
    VectorType plane_normal = VectorType::cross_product(v1, v2);
    if(axom::utilities::isNearlyEqual(plane_normal.norm(), 0.0, EPS))
    {
      plane_normal = VectorType::cross_product(v1, v3);
    }
    if(axom::utilities::isNearlyEqual(plane_normal.norm(), 0.0, EPS))
    {
      plane_normal = VectorType::cross_product(v2, v3);
    }
    plane_normal = plane_normal.unitVector();

    // Check all control points for simplicity
    for(int p = 0; p <= ord_u; ++p)
    {
      for(int q = ((p == 0) ? 1 : 0); q <= ord_v; ++q)
      {
        const double signedDist = plane_normal.dot(m_controlPoints(p, q) - m_controlPoints(0, 0));

        if(signedDist * signedDist > sq_tol)
        {
          return false;
        }
      }
    }
    return true;
  }

  /*!
   * \brief Predicate to check if the patch can be approximated by a polygon
   *
   * This function checks if a BezierPatch lies in a plane
   *  and that the edges are linear up to tolerance `sq_tol`
   *
   * \param [in] tol Threshold for sum of squared distances
   * \param [in] EPS Threshold for nearness to zero
   * \return True if the patch is a planar polygon, up to the desired tolerances
   */
  bool isPolygonal(double sq_tol = 1e-8, double EPS = 1e-8) const
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    if(ord_u <= 0 && ord_v <= 0)
    {
      return true;
    }
    if(ord_u == 1 && ord_v == 0)
    {
      return true;
    }
    if(ord_u == 0 && ord_v == 1)
    {
      return true;
    }

    // Check if the patch is planar
    if(!isPlanar(sq_tol, EPS))
    {
      return false;
    }

    // Check if each bounding curve is linear
    if(!isocurve_u(0).isLinear(sq_tol))
    {
      return false;
    }
    if(!isocurve_v(0).isLinear(sq_tol))
    {
      return false;
    }
    if(!isocurve_u(1).isLinear(sq_tol))
    {
      return false;
    }
    if(!isocurve_v(1).isLinear(sq_tol))
    {
      return false;
    }

    return true;
  }

  /*!
   * \brief Predicate to check if the Bezier patch is approximately bilinear
   *
   * This function checks if the patch is (nearly) bilinear.
   * A necessary condition for a geometrically bilinear patch is that each line of
   *  control points in the net is approximately linear. 
   * A necessary condition for a parametrically bilinear patch is that the control
   *  points are coincident with the surface of the bilinear patch defined by
   *  its corners evaluated at uniform parameter values, 
   *  i.e. the control points are also equally spaced on the net.
   *
   * \param [in] sq_tol Threshold for absolute squared distances
   * \param [in] useStrictBilinear If true, require the patch be parametrically bilinear
   * \return True if patch is bilinear up to tolerance \a sq_tol
   */
  bool isBilinear(double sq_tol = 1e-8, bool useStrictBilinear = false) const
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    if(ord_u <= 1 && ord_v <= 1)
    {
      return true;
    }

    if(useStrictBilinear)
    {
      // Anonymous function to evaluate the bilinear patch defined by the corners
      auto bilinear_patch = [&](T u, T v) -> PointType {
        PointType val;
        for(int N = 0; N < NDIMS; ++N)
        {
          val[N] = axom::utilities::lerp(
            axom::utilities::lerp(m_controlPoints(0, 0)[N], m_controlPoints(0, ord_v)[N], v),
            axom::utilities::lerp(m_controlPoints(ord_u, 0)[N], m_controlPoints(ord_u, ord_v)[N], v),
            u);
        }
        return val;
      };

      for(int u = 0; u <= ord_u; ++u)
      {
        for(int v = 0; v <= ord_v; ++v)
        {
          // Don't need to check the corners
          if((u == 0 && v == 0) || (u == 0 && v == ord_v) || (u == ord_u && v == 0) ||
             (u == ord_u && v == ord_v))
          {
            continue;
          }

          // Evaluate where the control point would be if the patch *was* bilinear
          PointType bilinear_point =
            bilinear_patch(u / static_cast<T>(ord_u), v / static_cast<T>(ord_v));

          if(squared_distance(m_controlPoints(u, v), bilinear_point) > sq_tol)
          {
            return false;
          }
        }
      }
    }
    else
    {
      for(int p = 0; p <= ord_u; ++p)
      {
        Segment<T, 3> seg(m_controlPoints(p, 0), m_controlPoints(p, ord_v));
        for(int q = 1; q < ord_v; ++q)
        {
          if(squared_distance(m_controlPoints(p, q), seg) > sq_tol)
          {
            return false;
          }
        }
      }

      for(int q = 0; q <= ord_v; ++q)
      {
        Segment<T, 3> seg(m_controlPoints(0, q), m_controlPoints(ord_u, q));
        for(int p = 1; p < ord_u; ++p)
        {
          if(squared_distance(m_controlPoints(p, q), seg) > sq_tol)
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
  /// \name Query/modify patch parametrization

  /*!
   * \brief Reverses the order of one direction of the Bezier patch's control points and weights
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

  /// Reverses the order of the Bezier patch's control points and weights on the first axis
  void reverseOrientation_u()
  {
    const int ord_u = getOrder_u();
    const int mid_u = (ord_u + 1) / 2;

    const int ord_v = getOrder_v();

    for(int q = 0; q <= ord_v; ++q)
    {
      for(int i = 0; i < mid_u; ++i)
      {
        axom::utilities::swap(m_controlPoints(i, q), m_controlPoints(ord_u - i, q));
      }

      if(isRational())
      {
        for(int i = 0; i < mid_u; ++i)
        {
          axom::utilities::swap(m_weights(i, q), m_weights(ord_u - i, q));
        }
      }
    }
  }

  /// Reverses the order of the Bezier patch's control points and weights on the second axis
  void reverseOrientation_v()
  {
    const int ord_u = getOrder_u();

    const int ord_v = getOrder_v();
    const int mid_v = (ord_v + 1) / 2;

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int i = 0; i < mid_v; ++i)
      {
        axom::utilities::swap(m_controlPoints(p, i), m_controlPoints(p, ord_v - i));
      }

      if(isRational())
      {
        for(int i = 0; i < mid_v; ++i)
        {
          axom::utilities::swap(m_weights(p, i), m_weights(p, ord_v - i));
        }
      }
    }
  }

  /// Swap the axes such that s(u, v) becomes s(v, u)
  void swapAxes()
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    CoordsMat new_controlPoints(ord_v + 1, ord_u + 1);

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        new_controlPoints(q, p) = m_controlPoints(p, q);
      }
    }

    m_controlPoints = new_controlPoints;

    if(isRational())
    {
      WeightsMat new_weights(ord_v + 1, ord_u + 1);
      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          new_weights(q, p) = m_weights(p, q);
        }
      }

      m_weights = new_weights;
    }
  }

  ///@}

  ///@{
  /// \name Functions to evaluate patch and its derivatives and normals along points or lines

  /*!
   * \brief Evaluates a slice Bezier patch for a fixed parameter value of \a u or \a v
   *
   * \param [in] u parameter value at which to evaluate the first axis
   * \param [in] v parameter value at which to evaluate the second axis
   * \param [in] axis orientation of curve. 0 for fixed u, 1 for fixed v
   * \return p the value of the Bezier patch at (u, v)
   *
   * \note We typically evaluate the patch at \a u or \a v between 0 and 1
   */
  BezierCurveType isocurve(T uv, int axis) const
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
   */
  BezierCurveType isocurve_u(T u) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    BezierCurveType c(ord_v);
    axom::Array<T> dCarray(ord_u + 1);

    if(!isRational())
    {
      // If looking for an edge, don't need to use De Casteljau
      if(u == 0.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(0, q);
        }
        return c;
      }
      if(u == 1.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(ord_u, q);
        }
        return c;
      }

      // Otherwise, do De Casteljau along the u-axis
      for(int q = 0; q <= ord_v; ++q)
      {
        for(int i = 0; i < NDIMS; ++i)
        {
          for(int p = 0; p <= ord_u; ++p)
          {
            dCarray[p] = m_controlPoints(p, q)[i];
          }

          for(int p = 1; p <= ord_u; ++p)
          {
            const int end = ord_u - p;
            for(int k = 0; k <= end; ++k)
            {
              dCarray[k] = lerp(dCarray[k], dCarray[k + 1], u);
            }
          }

          c[q][i] = dCarray[0];
        }
      }

      return c;
    }
    else
    {
      c.makeRational();

      // If looking for an edge, don't need to use De Casteljau
      if(u == 0.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(0, q);
          c.setWeight(q, m_weights(0, q));
        }
        return c;
      }
      if(u == 1.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(ord_u, q);
          c.setWeight(q, m_weights(ord_u, q));
        }
        return c;
      }

      // Otherwise, do de Casteljau on the homogeneous control points

      // Store BezierPatch of projective weights, (wx, wy, wz)
      //  and BezierPatch of weights (w)
      BezierPatch<T, NDIMS> projective(ord_u, ord_v);
      BezierPatch<T, 1> weights(ord_u, ord_v);

      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          weights(p, q)[0] = m_weights(p, q);

          for(int i = 0; i < NDIMS; ++i)
          {
            projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
          }
        }
      }

      BezierCurve<T, NDIMS> P = projective.isocurve_u(u);
      BezierCurve<T, 1> W = weights.isocurve_u(u);

      for(int q = 0; q <= ord_v; ++q)
      {
        for(int i = 0; i < NDIMS; ++i)
        {
          c[q][i] = P[q][i] / W[q][0];
          c.setWeight(q, W[q][0]);
        }
      }

      return c;
    }
  }

  /*!
   * \brief Returns an isocurve with a fixed value of v
   *
   * \param [in] v Parameter value fixed in the isocurve
   * \return c The isocurve C(u) = S(u, v) for fixed v
   */
  BezierCurveType isocurve_v(T v) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    BezierCurveType c(ord_u);
    axom::Array<T> dCarray(ord_v + 1);

    if(!isRational())
    {
      // If looking for an edge, don't need to use De Casteljau
      if(v == 0.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, 0);
        }
        return c;
      }
      if(v == 1.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, ord_v);
        }
        return c;
      }

      // Otherwise, do De Casteljau along the u-axis
      for(int p = 0; p <= ord_u; ++p)
      {
        for(int i = 0; i < NDIMS; ++i)
        {
          for(int q = 0; q <= ord_v; ++q)
          {
            dCarray[q] = m_controlPoints(p, q)[i];
          }

          for(int q = 1; q <= ord_v; ++q)
          {
            const int end = ord_v - q;
            for(int k = 0; k <= end; ++k)
            {
              dCarray[k] = lerp(dCarray[k], dCarray[k + 1], v);
            }
          }

          c[p][i] = dCarray[0];
        }
      }

      return c;
    }
    else
    {
      c.makeRational();

      // If looking for an edge, don't need to use De Casteljau
      if(v == 0.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, 0);
          c.setWeight(p, m_weights(p, 0));
        }
        return c;
      }
      if(v == 1.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, ord_v);
          c.setWeight(p, m_weights(p, ord_v));
        }
        return c;
      }

      // Otherwise, do de Casteljau on the homogeneous control points

      // Store BezierPatch of projective weights, (wx, wy, wz)
      //  and BezierPatch of weights (w)
      BezierPatch<T, NDIMS> projective(ord_u, ord_v);
      BezierPatch<T, 1> weights(ord_u, ord_v);

      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          weights(p, q)[0] = m_weights(p, q);

          for(int i = 0; i < NDIMS; ++i)
          {
            projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
          }
        }
      }

      BezierCurve<T, NDIMS> P = projective.isocurve_v(v);
      BezierCurve<T, 1> W = weights.isocurve_v(v);

      for(int p = 0; p <= ord_u; ++p)
      {
        for(int i = 0; i < NDIMS; ++i)
        {
          c[p][i] = P[p][i] / W[p][0];
          c.setWeight(p, W[p][0]);
        }
      }

      return c;
    }
  }

  /*!
   * \brief Evaluates a Bezier patch at a particular parameter value (\a u, \a v)
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return p the value of the Bezier patch at (u, v)
   *
   * \note We typically evaluate the patch at \a u and \a v between 0 and 1
   */
  PointType evaluate(T u, T v) const
  {
    if(getOrder_u() >= getOrder_v())
    {
      return isocurve_u(u).evaluate(v);
    }
    else
    {
      return isocurve_v(v).evaluate(u);
    }
  }

  /*!
   * \brief Evaluates all first derivatives Bezier patch at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the Bezier patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   *
   * \note We typically evaluate the patch at \a u and \a v between 0 and 1
   */
  void evaluateFirstDerivatives(T u,
                                T v,
                                Point<T, NDIMS>& eval,
                                Vector<T, NDIMS>& Du,
                                Vector<T, NDIMS>& Dv) const
  {
    using axom::utilities::lerp;
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);

    if(!isRational())
    {
      // Do de Casteljau until we get a 2x2
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          for(int q = 0; q <= ord_v; ++q)
          {
            dCmat(p, q) = m_controlPoints(p, q)[i];
          }
        }

        // Store the size after each de Casteljau reduction is made
        int end_u = ord_u;
        int end_v = ord_v;

        // Do de Casteljau over the longer direction first
        if(ord_u >= ord_v)
        {
          // Stop 1 steps early in u
          while(end_u > 1)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= end_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }

          // Stop 1 steps early in v
          while(end_v > 1)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }
        }
        else
        {
          // Stop 1 steps early in v
          while(end_v > 1)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }

          // Stop 2 steps early in v over the three columns
          while(end_u > 1)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= end_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }
        }

        // By now, the first 2x2 submatrix of dCmat should have
        //  all the intermediate values needed

        // Compute first order derivatives
        // clang-format off
        if( end_u == 1 && end_v == 1 )
        {
          Du[i] = ord_u * ( (1 - v) * (dCmat(1, 0) - dCmat(0, 0)) + 
                                  v * (dCmat(1, 1) - dCmat(0, 1)) );
          Dv[i] = ord_v * ( (1 - u) * (dCmat(0, 1) - dCmat(0, 0)) + 
                                  u * (dCmat(1, 1) - dCmat(1, 0)) );
          eval[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) + (1 - u) * v * dCmat(0, 1) + u * v * dCmat(1, 1);
        }
        else if (end_u == 1 && end_v == 0)
        {
          Du[i] = ord_u * (dCmat(1, 0) - dCmat(0, 0));
          Dv[i] = 0.0;
          eval[i] = (1 - u) * dCmat(0, 0) + u * dCmat(1, 0);
        }
        else if (end_u == 0 && end_v == 1)
        {
          Du[i] = 0.0;
          Dv[i] = ord_v * (dCmat(0, 1) - dCmat(0, 0));
          eval[i] = (1 - v) * dCmat(0, 0) + v * dCmat(0, 1);
        }
        else
        {
          Du[i] = 0.0;
          Dv[i] = 0.0;
          eval[i] = dCmat(0, 0);
        }
        // clang-format on
      }
    }
    else
    {
      // Store BezierPatch of projective weights, (wx, wy, wz)
      //  and BezierPatch of weights (w)
      BezierPatch<T, NDIMS> projective(ord_u, ord_v);
      BezierPatch<T, 1> weights(ord_u, ord_v);

      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          weights(p, q)[0] = m_weights(p, q);

          for(int i = 0; i < NDIMS; ++i)
          {
            projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
          }
        }
      }

      Point<T, NDIMS> P;
      Vector<T, NDIMS> P_u, P_v, P_uu, P_vv, P_uv;

      Point<T, 1> W;
      Vector<T, 1> W_u, W_v, W_uu, W_vv, W_uv;

      projective.evaluateSecondDerivatives(u, v, P, P_u, P_v, P_uu, P_vv, P_uv);
      weights.evaluateSecondDerivatives(u, v, W, W_u, W_v, W_uu, W_vv, W_uv);

      for(int i = 0; i < NDIMS; ++i)
      {
        eval[i] = P[i] / W[0];
        Du[i] = (P_u[i] - eval[i] * W_u[0]) / W[0];
        Dv[i] = (P_v[i] - eval[i] * W_v[0]) / W[0];
      }
    }
  }

  /*!
   * \brief Evaluates all linear derivatives Bezier patch at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the Bezier patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   *
   * \note We typically evaluate the patch at \a u and \a v between 0 and 1
   */
  void evaluate_linear_derivatives(T u,
                                   T v,
                                   Point<T, NDIMS>& eval,
                                   Vector<T, NDIMS>& Du,
                                   Vector<T, NDIMS>& Dv,
                                   Vector<T, NDIMS>& DuDv) const
  {
    using axom::utilities::lerp;
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);

    if(!isRational())
    {
      // Do de Casteljau until we get a 2x2
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          for(int q = 0; q <= ord_v; ++q)
          {
            dCmat(p, q) = m_controlPoints(p, q)[i];
          }
        }

        // Store the size after each de Casteljau reduction is made
        int end_u = ord_u;
        int end_v = ord_v;

        // Do de Casteljau over the longer direction first
        if(ord_u >= ord_v)
        {
          // Stop 1 steps early in u
          while(end_u > 1)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= end_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }

          // Stop 1 steps early in v
          while(end_v > 1)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }
        }
        else
        {
          // Stop 1 steps early in v
          while(end_v > 1)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }

          // Stop 2 steps early in v over the three columns
          while(end_u > 1)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= end_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }
        }

        // By now, the first 2x2 submatrix of dCmat should have
        //  all the intermediate values needed

        // Compute first order derivatives
        // clang-format off
        if( end_u == 1 && end_v == 1 )
        {
          Du[i] = ord_u * ( (1 - v) * (dCmat(1, 0) - dCmat(0, 0)) + 
                                  v * (dCmat(1, 1) - dCmat(0, 1)) );
          Dv[i] = ord_v * ( (1 - u) * (dCmat(0, 1) - dCmat(0, 0)) + 
                                  u * (dCmat(1, 1) - dCmat(1, 0)) );
          DuDv[i] = ord_u * ord_v * (dCmat(1, 1) - dCmat(1, 0) - dCmat(0, 1) + dCmat(0, 0));
          eval[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) + (1 - u) * v * dCmat(0, 1) + u * v * dCmat(1, 1);
        }
        else if (end_u == 1 && end_v == 0)
        {
          Du[i] = ord_u * (dCmat(1, 0) - dCmat(0, 0));
          Dv[i] = 0.0;
          DuDv[i] = 0.0;
          eval[i] = (1 - u) * dCmat(0, 0) + u * dCmat(1, 0);
        }
        else if (end_u == 0 && end_v == 1)
        {
          Du[i] = 0.0;
          Dv[i] = ord_v * (dCmat(0, 1) - dCmat(0, 0));
          DuDv[i] = 0.0;
          eval[i] = (1 - v) * dCmat(0, 0) + v * dCmat(0, 1);
        }
        else
        {
          Du[i] = 0.0;
          Dv[i] = 0.0;
          DuDv[i] = 0.0;
          eval[i] = dCmat(0, 0);
        }
        // clang-format on
      }
    }
    else
    {
      // Store BezierPatch of projective weights, (wx, wy, wz)
      //  and BezierPatch of weights (w)
      BezierPatch<T, NDIMS> projective(ord_u, ord_v);
      BezierPatch<T, 1> weights(ord_u, ord_v);

      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          weights(p, q)[0] = m_weights(p, q);

          for(int i = 0; i < NDIMS; ++i)
          {
            projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
          }
        }
      }

      Point<T, NDIMS> P;
      Vector<T, NDIMS> P_u, P_v, P_uu, P_vv, P_uv;

      Point<T, 1> W;
      Vector<T, 1> W_u, W_v, W_uu, W_vv, W_uv;

      projective.evaluateSecondDerivatives(u, v, P, P_u, P_v, P_uu, P_vv, P_uv);
      weights.evaluateSecondDerivatives(u, v, W, W_u, W_v, W_uu, W_vv, W_uv);

      for(int i = 0; i < NDIMS; ++i)
      {
        eval[i] = P[i] / W[0];
        Du[i] = (P_u[i] - eval[i] * W_u[0]) / W[0];
        Dv[i] = (P_v[i] - eval[i] * W_v[0]) / W[0];
        DuDv[i] = (P_uv[i] - Du[i] * W_v[0] - Dv[i] * W_u[0] - eval[i] * W_uv[0]) / W[0];
      }
    }
  }

  /*!
   * \brief Evaluates all second derivatives Bezier patch at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \param [out] eval The point value of the Bezier patch at (u, v)
   * \param [out] Du The vector value of S_u(u, v)
   * \param [out] Dv The vector value of S_v(u, v)
   * \param [out] DuDu The vector value of S_uu(u, v)
   * \param [out] DvDv The vector value of S_vv(u, v)
   * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
   *
   * \note We typically evaluate the patch at \a u and \a v between 0 and 1
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
    using axom::utilities::lerp;
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);

    if(!isRational())
    {
      // Do de Casteljau until we get a 3x3
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          for(int q = 0; q <= ord_v; ++q)
          {
            dCmat(p, q) = m_controlPoints(p, q)[i];
          }
        }

        // Store the size after each de Casteljau reduction is made
        int end_u = ord_u;
        int end_v = ord_v;

        // Do de Casteljau over the longer direction first
        if(ord_u >= ord_v)
        {
          // Stop 2 steps early in u
          while(end_u > 2)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= end_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }

          // Stop 2 steps early in v
          while(end_v > 2)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }
        }
        else
        {
          // Stop 2 steps early in v
          while(end_v > 2)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }

          // Stop 2 steps early in v over the three columns
          while(end_u > 2)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= end_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }
        }

        // By now, the first 3x3 submatrix of dCmat should have
        //  all the intermediate values needed

        // clang-format off
        // Get second order derivatives
        if( ord_u >= 2 )
        {
          if (ord_v == 0)
          {
            DuDu[i] = (ord_u - 1) * ord_u * 
                                    (dCmat(2, 0) - 2 * dCmat(1, 0) + dCmat(0, 0));
          }
          else if (ord_v == 1)
          {
            DuDu[i] = (ord_u - 1) * ord_u * 
                        ( (1 - v) * (dCmat(2, 0) - 2 * dCmat(1, 0) + dCmat(0, 0)) + 
                                v * (dCmat(2, 1) - 2 * dCmat(1, 1) + dCmat(0, 1)) ); 
          }
          else
          {
            DuDu[i] = (ord_u - 1) * ord_u * 
              ( (1 - v) * (1 - v) * (dCmat(2, 0) - 2 * dCmat(1, 0) + dCmat(0, 0)) + 
                  2 * v * (1 - v) * (dCmat(2, 1) - 2 * dCmat(1, 1) + dCmat(0, 1)) + 
                            v * v * (dCmat(2, 2) - 2 * dCmat(1, 2) + dCmat(0, 2)) );
          }
        }
        else
        {
          DuDu[i] = 0.0;
        }

        if( ord_v >= 2 )
        {
          if( ord_u == 0 )
          {
            DvDv[i] = (ord_v - 1) * ord_v * 
                                    (dCmat(0, 2) - 2 * dCmat(0, 1) + dCmat(0, 0)); 
          }
          else if (ord_u == 1)
          {
            DvDv[i] = (ord_v - 1) * ord_v * 
                        ( (1 - u) * (dCmat(0, 2) - 2 * dCmat(0, 1) + dCmat(0, 0)) + 
                                u * (dCmat(1, 2) - 2 * dCmat(1, 1) + dCmat(1, 0)) );
          }
          else
          {
            DvDv[i] = (ord_v - 1) * ord_v * 
              ( (1 - u) * (1 - u) * (dCmat(0, 2) - 2 * dCmat(0, 1) + dCmat(0, 0)) + 
                  2 * u * (1 - u) * (dCmat(1, 2) - 2 * dCmat(1, 1) + dCmat(1, 0)) + 
                            u * u * (dCmat(2, 2) - 2 * dCmat(2, 1) + dCmat(2, 0)) ); 
          }
        }
        else
        {
          DvDv[i] = 0.0;
        }
        
        // Do one de Casteljau iteration in u
        while( end_u > 1 )
        {
          end_u -= 1;
          for(int k = 0; k <= end_u; ++k)
          {
            for(int q = 0; q <= end_v; ++q)
            {
              dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
            }
          }
        }

        // Do (at most) one de Casteljau iteration in v
        while(end_v > 1)
        {
          end_v -= 1;
          for(int k = 0; k <= end_v; ++k)
          {
            for(int p = 0; p <= end_u; ++p)
            {
              dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
            }
          }
        }

        // Compute first order derivatives
        // clang-format off
        if( end_u == 1 && end_v == 1 )
        {
          Du[i] = ord_u * ( (1 - v) * (dCmat(1, 0) - dCmat(0, 0)) + 
                                  v * (dCmat(1, 1) - dCmat(0, 1)) );
          Dv[i] = ord_v * ( (1 - u) * (dCmat(0, 1) - dCmat(0, 0)) + 
                                  u * (dCmat(1, 1) - dCmat(1, 0)) );
          DuDv[i] = ord_u * ord_v * (dCmat(1, 1) - dCmat(1, 0) - dCmat(0, 1) + dCmat(0, 0));
          eval[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) + (1 - u) * v * dCmat(0, 1) + u * v * dCmat(1, 1);
        }
        else if (end_u == 1 && end_v == 0)
        {
          Du[i] = ord_u * (dCmat(1, 0) - dCmat(0, 0));
          Dv[i] = 0.0;
          DuDv[i] = 0.0;
          eval[i] = (1 - u) * dCmat(0, 0) + u * dCmat(1, 0);
        }
        else if (end_u == 0 && end_v == 1)
        {
          Du[i] = 0.0;
          Dv[i] = ord_v * (dCmat(0, 1) - dCmat(0, 0));
          DuDv[i] = 0.0;
          eval[i] = (1 - v) * dCmat(0, 0) + v * dCmat(0, 1);
        }
        else
        {
          Du[i] = 0.0;
          Dv[i] = 0.0;
          DuDv[i] = 0.0;
          eval[i] = dCmat(0, 0);
        }
        // clang-format on
      }
    }
    else
    {
      // Store BezierPatch of projective weights, (wx, wy, wz)
      //  and BezierPatch of weights (w)
      BezierPatch<T, NDIMS> projective(ord_u, ord_v);
      BezierPatch<T, 1> weights(ord_u, ord_v);

      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          weights(p, q)[0] = m_weights(p, q);

          for(int i = 0; i < NDIMS; ++i)
          {
            projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
          }
        }
      }

      Point<T, NDIMS> P;
      Vector<T, NDIMS> P_u, P_v, P_uu, P_vv, P_uv;

      Point<T, 1> W;
      Vector<T, 1> W_u, W_v, W_uu, W_vv, W_uv;

      projective.evaluateSecondDerivatives(u, v, P, P_u, P_v, P_uu, P_vv, P_uv);
      weights.evaluateSecondDerivatives(u, v, W, W_u, W_v, W_uu, W_vv, W_uv);

      for(int i = 0; i < NDIMS; ++i)
      {
        eval[i] = P[i] / W[0];
        Du[i] = (P_u[i] - eval[i] * W_u[0]) / W[0];
        Dv[i] = (P_v[i] - eval[i] * W_v[0]) / W[0];
        DuDu[i] = (P_uu[i] - 2 * W_u[0] * Du[i] - eval[i] * W_uu[0]) / W[0];
        DvDv[i] = (P_vv[i] - 2 * W_v[0] * Dv[i] - eval[i] * W_vv[0]) / W[0];
        DuDv[i] = (P_uv[i] - Du[i] * W_v[0] - Dv[i] * W_u[0] - eval[i] * W_uv[0]) / W[0];
      }
    }
  }

  /*!
   * \brief Computes a tangent of a Bezier patch at a particular parameter value (\a u, \a v) along the u axis
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return vec a tangent vector in u of the Bezier patch at (u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType du(T u, T v) const { return isocurve_v(v).dt(u); }

  /*!
   * \brief Computes the second derivative of a Bezier patch at (\a u, \a v) along the u axis
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \return vec The vector value of S_uu(u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType dudu(T u, T v) const { return isocurve_v(v).dtdt(u); }

  /*!
   * \brief Computes a tangent of a Bezier patch at a particular parameter value (\a u, \a v) along the v axis
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return vec a tangent vector in v of the Bezier patch at (u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType dv(T u, T v) const { return isocurve_u(u).dt(v); }

  /*!
   * \brief Computes the second derivative of a Bezier patch at (\a u, \a v) along the v axis
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \return vec The vector value of S_vv(u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType dvdv(T u, T v) const { return isocurve_u(u).dtdt(v); }

  /*!
   * \brief Computes the mixed second derivative of a Bezier patch at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \return vec The vector value of S_uv(u, v) == S_vu(u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType dudv(T u, T v) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    // If the curve is nonrational, we can use standard de Casteljau
    if(!isRational())
    {
      Vector<T, NDIMS> val;

      axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);
      axom::Array<T, 2> dWmat(ord_u + 1, ord_v + 1);

      // Do de Casteljau until we get a 2x2
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          for(int q = 0; q <= ord_v; ++q)
          {
            dCmat(p, q) = m_controlPoints(p, q)[i];
          }
        }

        // Store the size after each de Casteljau reduction is made
        int end_u = ord_u;
        int end_v = ord_v;

        // Do de Casteljau over the longer direction first
        if(ord_u >= ord_v)
        {
          // Stop 1 steps early in u
          while(end_u > 1)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= ord_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }

          // Stop 1 steps early in u
          while(end_v > 1)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }
        }
        else
        {
          // Stop 1 step early in v
          while(end_v > 1)
          {
            end_v -= 1;
            for(int k = 0; k <= end_v; ++k)
            {
              for(int p = 0; p <= end_u; ++p)
              {
                dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
              }
            }
          }

          // Stop 1 step early in v over the two columns
          while(end_u > 1)
          {
            end_u -= 1;
            for(int k = 0; k <= end_u; ++k)
            {
              for(int q = 0; q <= ord_v; ++q)
              {
                dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
              }
            }
          }
        }

        // Do the evaluation, accounting for insufficient control points
        if(end_u == 0 || end_v == 0)
        {
          val[i] = 0.0;
        }
        else
        {
          // clang-format off
        val[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) +
                       (1 - u) * v * dCmat(0, 1) +       u * v * dCmat(1, 1);
          // clang-format on
        }
      }

      return val;
    }
    // If rational, construct the 4D homogeneous surface,
    //  which requires all first derivatives
    else
    {
      Vector<T, NDIMS> val;

      // Store BezierPatch of projective weights, (wx, wy, wz)
      //  and BezierPatch of weights (w)
      BezierPatch<T, NDIMS> projective(ord_u, ord_v);
      BezierPatch<T, 1> weights(ord_u, ord_v);

      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          weights(p, q)[0] = m_weights(p, q);

          for(int i = 0; i < NDIMS; ++i)
          {
            projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
          }
        }
      }

      Point<T, NDIMS> P;
      Vector<T, NDIMS> P_u, P_v, P_uv;

      Point<T, 1> W;
      Vector<T, 1> W_u, W_v, W_uv;

      projective.evaluate_linear_derivatives(u, v, P, P_u, P_v, P_uv);
      weights.evaluate_linear_derivatives(u, v, W, W_u, W_v, W_uv);

      // Store values used in each coordinate computation
      double weight_prod = 2 * W_u[0] * W_v[0] - W[0] * W_uv[0];
      double weight_cubed = W[0] * W[0] * W[0];
      for(int i = 0; i < NDIMS; ++i)
      {
        val[i] =
          W[0] * W[0] * P_uv[i] - W[0] * (P_u[i] * W_v[0] + P_v[i] * W_u[0]) + P[i] * weight_prod;
        val[i] /= weight_cubed;
      }

      return val;
    }
  }

  /*!
   * \brief Computes the mixed second derivative of a Bezier patch at (\a u, \a v)
   *
   * \param [in] u Parameter value at which to evaluate on the first axis
   * \param [in] v Parameter value at which to evaluate on the second axis
   * \return vec The vector value of S_uv(u, v) == S_vu(u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType dvdu(T u, T v) const { return dudv(u, v); }

  /*!
   * \brief Computes the normal vector of a Bezier patch at a particular parameter value (\a u, \a v)
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return vec the normal vector of the Bezier patch at (u, v)
   *
   * \note We typically find the normal of the patch at \a u and \a v between 0 and 1
   */
  VectorType normal(T u, T v) const
  {
    Point<T, NDIMS> eval;
    Vector<T, NDIMS> Du, Dv;
    evaluateFirstDerivatives(u, v, eval, Du, Dv);
    return VectorType::cross_product(Du, Dv);
  }

  ///@}

  ///@{
  /// \name Functions dealing with patch subdivision

  /*!
   * \brief Splits a Bezier patch into two Bezier patches
   *
   * \param [in] uv parameter value between 0 and 1 at which to bisect the patch
   * \param [in] axis orientation of split. 0 for fixed u, 1 for fixed v
   * \param [out] p1 First output Bezier patch
   * \param [out] p2 Second output Bezier patch
   *
   * \pre Parameter \a uv must be between 0 and 1
   */
  void split(T uv, int axis, BezierPatch& p1, BezierPatch& p2) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));

    if(axis == 0)
    {
      split_u(uv, p1, p2);
    }
    else
    {
      split_v(uv, p1, p2);
    }
  }

  /// Split the patch along a fixed value of u
  void split_u(T u, BezierPatch& p1, BezierPatch& p2) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    p1.setOrder(ord_u, ord_v);

    // Note: The second patch's control points are computed inline
    //       as we find the first patch's control points
    p2 = *this;

    if(isRational())
    {
      p1.makeRational();  // p2 already rational

      for(int q = 0; q <= ord_v; ++q)
      {
        // Do the rational de Casteljau algorithm
        p1(0, q) = p2(0, q);
        p1.setWeight(0, q, p2.getWeight(0, q));

        for(int p = 1; p <= ord_u; ++p)
        {
          const int end = ord_u - p;

          for(int k = 0; k <= end; ++k)
          {
            double temp_weight = lerp(p2.getWeight(k, q), p2.getWeight(k + 1, q), u);

            for(int i = 0; i < NDIMS; ++i)
            {
              p2(k, q)[i] =
                lerp(p2.getWeight(k, q) * p2(k, q)[i], p2.getWeight(k + 1, q) * p2(k + 1, q)[i], u) /
                temp_weight;
            }

            p2.setWeight(k, q, temp_weight);
          }

          p1(p, q) = p2(0, q);
          p1.setWeight(p, q, p2.getWeight(0, q));
        }
      }
    }
    else
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        p1(0, q) = m_controlPoints(0, q);

        for(int i = 0; i < NDIMS; ++i)
        {
          for(int p = 1; p <= ord_u; ++p)
          {
            const int end = ord_u - p;
            for(int k = 0; k <= end; ++k)
            {
              p2(k, q)[i] = lerp(p2(k, q)[i], p2(k + 1, q)[i], u);
            }

            p1(p, q)[i] = p2(0, q)[i];
          }
        }
      }
    }
  }

  void split_v(T v, BezierPatch& p1, BezierPatch& p2) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    p1.setOrder(ord_u, ord_v);

    // Note: The second patch's control points are computed inline
    //       as we find the first patch's control points
    p2 = *this;

    if(isRational())
    {
      p1.makeRational();  // p2 already rational

      for(int p = 0; p <= ord_u; ++p)
      {
        // Do the rational de Casteljau algorithm
        p1(p, 0) = p2(p, 0);
        p1.setWeight(p, 0, p2.getWeight(p, 0));

        for(int q = 1; q <= ord_v; ++q)
        {
          const int end = ord_v - q;

          for(int k = 0; k <= end; ++k)
          {
            double temp_weight = lerp(p2.getWeight(p, k), p2.getWeight(p, k + 1), v);

            for(int i = 0; i < NDIMS; ++i)
            {
              p2(p, k)[i] =
                lerp(p2.getWeight(p, k) * p2(p, k)[i], p2.getWeight(p, k + 1) * p2(p, k + 1)[i], v) /
                temp_weight;
            }

            p2.setWeight(p, k, temp_weight);
          }

          p1(p, q) = p2(p, 0);
          p1.setWeight(p, q, p2.getWeight(p, 0));
        }
      }
    }
    else
    {
      for(int p = 0; p <= ord_u; ++p)
      {
        p1(p, 0) = m_controlPoints(p, 0);

        for(int i = 0; i < NDIMS; ++i)
        {
          for(int q = 1; q <= ord_v; ++q)
          {
            const int end = ord_v - q;
            for(int k = 0; k <= end; ++k)
            {
              p2(p, k)[i] = lerp(p2(p, k)[i], p2(p, k + 1)[i], v);
            }

            p1(p, q)[i] = p2(p, 0)[i];
          }
        }
      }
    }
  }

  /*!
   * \brief Splits a Bezier patch into four Bezier patches
   *
   * \param [in] u parameter value between 0 and 1 at which to bisect on the first axis
   * \param [in] v parameter value between 0 and 1 at which to bisect on the second axis
   * \param [out] p1 First output Bezier patch
   * \param [out] p2 Second output Bezier patch
   * \param [out] p3 Third output Bezier patch
   * \param [out] p4 Fourth output Bezier patch
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
   * \pre Parameter \a u and \a v must be between 0 and 1
   */
  void split(T u, T v, BezierPatch& p1, BezierPatch& p2, BezierPatch& p3, BezierPatch& p4) const
  {
    // Bisect the patch along the u direction
    split_u(u, p1, p2);

    // Temporarily store the result in each half and split again
    BezierPatch p0(p1);
    p0.split_v(v, p1, p3);

    p0 = p2;
    p0.split_v(v, p2, p4);
  }

  ///@}

  /*!
   * \brief Simple formatted print of a Bezier Patch instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    os << "{ order (" << ord_u << ',' << ord_v << ") Bezier Patch ";

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        os << m_controlPoints(p, q) << ((p < ord_u || q < ord_v) ? "," : "");
      }
    }

    if(isRational())
    {
      os << ", weights [";
      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          os << m_weights(p, q) << ((p < ord_u || q < ord_v) ? "," : "");
        }
      }
    }
    os << "}";

    return os;
  }

  /// Checks equality of two Bezier Patches
  friend inline bool operator==(const BezierPatch<T, NDIMS>& lhs, const BezierPatch<T, NDIMS>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) && (lhs.m_weights == rhs.m_weights);
  }

  friend inline bool operator!=(const BezierPatch<T, NDIMS>& lhs, const BezierPatch<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

private:
  /// Check that the weights used are positive, and
  ///  that there is one for each control node
  bool isValidRational() const
  {
    if(!isRational())
    {
      return true;
    }

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    if(m_weights.shape()[0] != (ord_u + 1) || m_weights.shape()[1] != (ord_v + 1))
    {
      return false;
    }

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        if(m_weights(p, q) <= 0)
        {
          return false;
        }
      }
    }

    return true;
  }

private:
  CoordsMat m_controlPoints;
  WeightsMat m_weights;
};

//------------------------------------------------------------------------------
/// Free functions related to BezierPatch
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierPatch<T, NDIMS>& bPatch)
{
  bPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_BEZIERPATCH_HPP_
