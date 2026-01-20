// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number_2d_memoization.hpp
 *
 * \brief Consists of data structures that accelerate GWN queries through "memoization,"
 *         i.e. dynamically caching and reusing intermediate curve subdivisions.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_2D_MEMOIZATION_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_2D_MEMOIZATION_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/KnotVector.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include "axom/primal/operators/is_convex.hpp"

#include <vector>
#include <ostream>
#include <math.h>

#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \struct BezierCurveData
 *
 * \brief Stores BezierCurves and relevant data/flags
 */
template <typename T>
struct BezierCurveData
{
  BezierCurveData() = default;

  BezierCurveData(const BezierCurve<T, 2>& a_curve, bool knownConvex, double bbExpansionAmount = 0.0)
    : m_curve(a_curve)
  {
    m_isConvexControlPolygon =
      knownConvex ? true : is_convex(Polygon<T, 2>(m_curve.getControlPoints()));
    m_boundingBox = m_curve.boundingBox().expand(bbExpansionAmount);
  }

  const auto& getCurve() const { return m_curve; }
  auto isConvexControlPolygon() const { return m_isConvexControlPolygon; }
  auto getBoundingBox() const { return m_boundingBox; }

  friend bool operator==(const BezierCurveData<T>& lhs, const BezierCurveData<T>& rhs)
  {
    // isConvexControlPolygon will be equal if the curves are
    return (lhs.m_curve == rhs.m_curve) && (lhs.m_boundingBox == rhs.m_boundingBox) &&
      (lhs.m_isConvexControlPolygon == rhs.m_isConvexControlPolygon);
  }

  friend bool operator!=(const BezierCurveData<T>& lhs, const BezierCurveData<T>& rhs)
  {
    return !(lhs == rhs);
  }

private:
  BezierCurve<T, 2> m_curve;
  bool m_isConvexControlPolygon;
  BoundingBox<T, 2> m_boundingBox;
};

// Forward declare the templated classes and operator functions
template <typename T>
class NURBSCurveGWNCache;

/// \brief Overloaded output operator for Cached Curves
template <typename T>
std::ostream& operator<<(std::ostream& os, const NURBSCurveGWNCache<T>& nCurveCache);

/*!
 * \class NURBSCurveGWNCache
 *
 * \brief Represents a NURBS curve and associated data for GWN evaluation
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * Stores subdivision Bezier curves, bounding boxes for each, and flags that track
 *  if the control polygon is known to be convex
 * 
 * \pre Assumes a 2D NURBS curve
 */
template <typename T>
class NURBSCurveGWNCache
{
public:
  using NumericType = T;
  using PointType = typename NURBSCurve<T, 2>::PointType;
  using VectorType = typename NURBSCurve<T, 2>::VectorType;
  using BoundingBoxType = typename NURBSCurve<T, 2>::BoundingBoxType;

public:
  NURBSCurveGWNCache() = default;

  /// \brief Initialize the cache with the data for the original curve
  NURBSCurveGWNCache(const NURBSCurve<T, 2>& a_curve, double bbExpansionAmount = 0.0)
  {
    m_boundingBox = a_curve.boundingBox();
    m_degree = a_curve.getDegree();
    m_numControlPoints = a_curve.getNumControlPoints();
    m_numSpans = a_curve.getNumKnotSpans();

    m_bezierSubdivisionMaps.resize(m_numSpans);
    auto beziers = a_curve.extractBezier();

    for(int idx = 0; idx < m_numSpans; ++idx)
    {
      m_bezierSubdivisionMaps[idx][std::make_pair(0, 0)] =
        BezierCurveData<T>(beziers[idx], false, bbExpansionAmount);
    }

    m_initPoint = a_curve[0];
    m_endPoint = a_curve[m_numControlPoints - 1];
  }

  /// \brief Initialize the cache with the data for a single Bezier curve
  NURBSCurveGWNCache(const BezierCurve<T, 2>& a_curve, double bbExpansionAmount = 0.0)
  {
    m_boundingBox = a_curve.boundingBox();
    m_degree = a_curve.getOrder();
    m_numControlPoints = a_curve.getOrder() + 1;

    if(a_curve.getOrder() <= 0)
    {
      m_numSpans = 0;
      m_bezierSubdivisionMaps.clear();

      m_initPoint = m_endPoint = Point<T, 2> {0.0, 0.0};
    }
    else
    {
      m_numSpans = 1;
      m_bezierSubdivisionMaps.resize(1);

      m_initPoint = a_curve[0];
      m_endPoint = a_curve[m_degree];

      m_bezierSubdivisionMaps[0][std::make_pair(0, 0)] =
        BezierCurveData<T>(a_curve, false, bbExpansionAmount);
    }
  }

  /// \brief Query the map. If curve is not found, add it and its pair from subdivision
  const BezierCurveData<T>& getSubdivisionData(int idx,
                                               int refinementLevel,
                                               int refinementIndex,
                                               double bbExpansionAmount = 0.0) const
  {
    using Key = std::pair<int, int>;
    auto& level_map = m_bezierSubdivisionMaps[idx];
    const Key hash_key {refinementLevel, refinementIndex};

    // If already there, return it
    if(auto it = level_map.find(hash_key); it != level_map.end())
    {
      return it->second;
    }

    // Otherwise, create (refinementLevel, refinementIndex) and sibling via their parent
    const Key parent_key {refinementLevel - 1, refinementIndex / 2};
    auto parent_it = level_map.find(parent_key);
    SLIC_ASSERT(parent_it != level_map.end());

    const BezierCurveData<T>& supercurve_data = parent_it->second;
    BezierCurve<T, 2> sub1, sub2;
    supercurve_data.getCurve().split(0.5, sub1, sub2);

    // Make keys for the requested curve and its "sibling" in the heirarchy
    const int base = refinementIndex - (refinementIndex % 2);
    const Key key1 {refinementLevel, base};
    const Key key2 {refinementLevel, base + 1};

    // Emplace both and return value associated with hash_key
    auto [it1, ins1] =
      level_map.try_emplace(key1, sub1, supercurve_data.isConvexControlPolygon(), bbExpansionAmount);
    auto [it2, ins2] =
      level_map.try_emplace(key2, sub2, supercurve_data.isConvexControlPolygon(), bbExpansionAmount);
    return (hash_key == key1) ? it1->second : it2->second;
  }

  ///@{
  //! \name Functions that mirror functionality of NURBSCurve and BezierCurve so signatures match in GWN evaluation.
  //!
  //! By limiting access to these functions, we ensure memoized information is always accurate
  auto getNumKnotSpans() const { return m_numSpans; }
  auto boundingBox() const { return m_boundingBox; }
  auto getNumControlPoints() const { return m_numControlPoints; }
  auto getDegree() const { return m_degree; }

  const auto& getInitPoint() const { return m_initPoint; }
  const auto& getEndPoint() const { return m_endPoint; }
  //@}

  friend bool operator==(const NURBSCurveGWNCache<T>& lhs, const NURBSCurveGWNCache<T>& rhs)
  {
    // numControlPoints, degree, and numSpans will be equal if the subdivision maps are
    return (lhs.m_bezierSubdivisionMaps == rhs.m_bezierSubdivisionMaps) &&
      (lhs.m_boundingBox == rhs.m_boundingBox);
  }

  friend bool operator!=(const NURBSCurveGWNCache<T>& lhs, const NURBSCurveGWNCache<T>& rhs)
  {
    return !(lhs == rhs);
  }

  std::ostream& print(std::ostream& os) const
  {
    os << "{ NURBSCurveGWNCache object with " << m_numSpans << " extracted bezier curves: ";

    if(m_numSpans >= 1)
    {
      os << m_bezierSubdivisionMaps[0][std::make_pair(0, 0)].getCurve();
    }
    for(int i = 1; i < m_numSpans; ++i)
    {
      os << ", " << m_bezierSubdivisionMaps[i][std::make_pair(0, 0)].getCurve();
    }
    os << "}";

    return os;
  }

private:
  BoundingBox<T, 2> m_boundingBox;
  int m_numControlPoints;
  int m_degree;
  int m_numSpans;

  Point<T, 2> m_initPoint, m_endPoint;

  mutable axom::Array<std::map<std::pair<int, int>, BezierCurveData<T>>> m_bezierSubdivisionMaps;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const NURBSCurveGWNCache<T>& nCurveCache)
{
  nCurveCache.print(os);
  return os;
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_2D_MEMOIZATION_HPP_
