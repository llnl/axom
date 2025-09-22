// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file NURBSCurve.hpp
 *
 * \brief A NURBS curve primitive
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
    : curve(a_curve)
  {
    isConvexControlPolygon = knownConvex ? true : is_convex(Polygon<T, 2>(curve.getControlPoints()));
    boundingBox = curve.boundingBox().expand(bbExpansionAmount);
  }

  auto getCurve() const { return curve; }
  auto isConvexControlPolygon() const { return isConvexControlPolygon; }
  auto getBoundingBox() const { return boundingBox; }

private:
  BezierCurve<T, 2> curve;
  bool isConvexControlPolygon;
  BoundingBox<T, 2> boundingBox;
};

/*!
 * \class NURBSCurveGWNCache
 *
 * \brief Represents a NURBS curve and associated data for GWN evaluation
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * Stores subdivision NURBS curves, and a flag that maintains if the control polygon is convex
 * 
 * \pre Assumes a 2D NURBS curve
 */
template <typename T>
class NURBSCurveGWNCache
{
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

  /// \brief Query the map. If curve is not found, add it and it's pair from subdivision
  const BezierCurveData<T>& getSubdivisionData(int idx,
                                               int refinementLevel,
                                               int refinementIndex,
                                               double bbExpansionAmount = 0.0) const
  {
    auto hash_key = std::make_pair(refinementLevel, refinementIndex);

    if(m_bezierSubdivisionMaps[idx].find(hash_key) == m_bezierSubdivisionMaps[idx].end())
    {
      const BezierCurveData<T>& supercurve_data =
        m_bezierSubdivisionMaps[idx][std::make_pair(refinementLevel - 1, refinementIndex / 2)];

      BezierCurve<T, 2> sub1, sub2;
      supercurve_data.curve().split(0.5, sub1, sub2);

      // Make keys for the requested curve and its "sibling" in the heirarchy
      const auto key1 = std::make_pair(refinementLevel, refinementIndex - refinementIndex % 2);
      const auto key2 = std::make_pair(refinementLevel, refinementIndex - refinementIndex % 2 + 1);

      // Populate the map
      m_bezierSubdivisionMaps[idx][key1] =
        BezierCurveData<T>(sub1, supercurve_data.isConvexControlPolygon(), bbExpansionAmount);
      m_bezierSubdivisionMaps[idx][key2] =
        BezierCurveData<T>(sub2, supercurve_data.isConvexControlPolygon(), bbExpansionAmount);
    }

    return m_bezierSubdivisionMaps[idx][hash_key];
  }

  auto getNumKnotSpans() const { return m_numSpans; }
  auto boundingBox() const { return m_boundingBox; }
  auto getNumControlPoints() const { return m_numControlPoints; }
  auto getDegree() const { return m_degree; }

  auto getInitPoint() const { return m_initPoint; }
  auto getEndPoint() const { return m_endPoint; }

private:
  BoundingBox<T, 2> m_boundingBox;
  int m_numControlPoints;
  int m_degree;
  int m_numSpans;

  Point<T, 2> m_initPoint, m_endPoint;

  mutable axom::Array<std::map<std::pair<int, int>, BezierCurveData<T>>> m_bezierSubdivisionMaps;
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_2D_MEMOIZATION_HPP_
