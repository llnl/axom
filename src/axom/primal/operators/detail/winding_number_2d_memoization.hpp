// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include <unordered_map>
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
    m_maxControlPointChordDistanceSq = computeMaxControlPointChordDistanceSq(m_curve);
  }

  const auto& getCurve() const { return m_curve; }
  auto isConvexControlPolygon() const { return m_isConvexControlPolygon; }
  auto getBoundingBox() const { return m_boundingBox; }
  bool isLinear(double tol) const { return m_maxControlPointChordDistanceSq <= tol; }

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
  static double computeMaxControlPointChordDistanceSq(const BezierCurve<T, 2>& curve)
  {
    const int order = curve.getOrder();
    if(order <= 1)
    {
      return 0.0;
    }

    Segment<T, 2> chord(curve[0], curve[order]);
    double maxSqDist = 0.0;
    for(int p = 1; p < order; ++p)
    {
      maxSqDist = std::max(maxSqDist, squared_distance(curve[p], chord));
    }
    return maxSqDist;
  }

  BezierCurve<T, 2> m_curve;
  bool m_isConvexControlPolygon;
  BoundingBox<T, 2> m_boundingBox;
  double m_maxControlPointChordDistanceSq {0.0};
};

// Forward declare the templated classes and operator functions
template <typename T>
class NURBSCurveGWNCache;

struct PairHash
{
  using argument_type = std::pair<int, int>;
  using result_type = std::size_t;

  result_type operator()(const argument_type& key) const noexcept
  {
    const auto a = static_cast<std::uint32_t>(key.first);
    const auto b = static_cast<std::uint32_t>(key.second);
    return (static_cast<result_type>(a) << 32) ^ static_cast<result_type>(b);
  }
};

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
      m_bezierSubdivisionMaps[idx].reserve(32);
      m_bezierSubdivisionMaps[idx].try_emplace(std::make_pair(0, 0),
                                               beziers[idx],
                                               false,
                                               bbExpansionAmount);
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

      m_bezierSubdivisionMaps[0].reserve(32);
      m_bezierSubdivisionMaps[0].try_emplace(std::make_pair(0, 0), a_curve, false, bbExpansionAmount);
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

  mutable axom::Array<std::unordered_map<std::pair<int, int>, BezierCurveData<T>, PairHash>>
    m_bezierSubdivisionMaps;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const NURBSCurveGWNCache<T>& nCurveCache)
{
  nCurveCache.print(os);
  return os;
}

}  // namespace detail

/*!
 * \brief Manage an array of NURBSCurveGWNCache<double>
 */
class NURBSCurveCacheManager
{
  using NURBSCache = axom::primal::detail::NURBSCurveGWNCache<double>;
  using NURBSCacheArray = axom::Array<NURBSCache>;
  using NURBSCacheArrayView = axom::ArrayView<const NURBSCache>;

  using CurveArrayView = axom::ArrayView<const axom::primal::NURBSCurve<double, 2>>;

public:
  NURBSCurveCacheManager() = default;

  NURBSCurveCacheManager(CurveArrayView curves, double bbExpansionAmount = 0.0)
  {
    for(auto& curve : curves)
    {
      m_nurbs_caches.push_back(NURBSCache(curve, bbExpansionAmount));
    }
  }

  /// A view of the manager object.
  struct View
  {
    NURBSCacheArrayView m_view;

    /// Return the NURBSCacheArrayView.
    NURBSCacheArrayView caches() const { return m_view; }
  };

  /// Return a view of this manager to pass into a device function.
  View view() const { return View {m_nurbs_caches.view()}; }

  /// Return if the underlying array is empty
  bool empty() const { return m_nurbs_caches.empty(); }

private:
  NURBSCacheArray m_nurbs_caches;
};

template <typename ExecSpace>
struct nurbs_cache_2d_traits
{
  using type = NURBSCurveCacheManager;
};

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
/*!
 * \brief Manage per-thread arrays of NURBSCurveGWNCache<double>
 */
class NURBSCurveCacheManagerOMP
{
  using NURBSCache = axom::primal::detail::NURBSCurveGWNCache<double>;
  using NURBSCachePerThreadArray = axom::Array<axom::Array<NURBSCache>>;
  using NURBSCachePerThreadArrayView = axom::ArrayView<const axom::Array<NURBSCache>>;
  using NURBSCacheArrayView = axom::ArrayView<const NURBSCache>;

  using CurveArrayView = axom::ArrayView<const axom::primal::NURBSCurve<double, 2>>;

public:
  NURBSCurveCacheManagerOMP() = default;

  NURBSCurveCacheManagerOMP(CurveArrayView curves, double bbExpansionAmount = 0.0)
  {
    const int nt = omp_get_max_threads();
    m_nurbs_caches.resize(nt);
    auto nurbs_caches_view = m_nurbs_caches.view();

    // Make the first one
    nurbs_caches_view[0].resize(curves.size());
    axom::for_all<axom::OMP_EXEC>(
      curves.size(),
      AXOM_LAMBDA(axom::IndexType i) {
        nurbs_caches_view[0][i] = NURBSCache(curves[i], bbExpansionAmount);
      });

    // Copy the constructed cache to the other threads' copies (less work than construction)
    axom::for_all<axom::OMP_EXEC>(
      1,
      nt,
      AXOM_LAMBDA(axom::IndexType t) { nurbs_caches_view[t] = nurbs_caches_view[0]; });
  }

  /// A view of the manager object.
  struct View
  {
    NURBSCachePerThreadArrayView m_views;

    /// Return the NURBSCacheArrayView for the current OMP thread.
    NURBSCacheArrayView caches() const { return m_views[omp_get_thread_num()].view(); }
  };

  /// Return a view of this manager to pass into a device function.
  View view() const { return View {m_nurbs_caches.view()}; }

  /// Return if the underlying array is empty
  bool empty() const { return m_nurbs_caches.empty(); }

private:
  NURBSCachePerThreadArray m_nurbs_caches;
};

template <>
struct nurbs_cache_2d_traits<axom::OMP_EXEC>
{
  using type = NURBSCurveCacheManagerOMP;
};
#endif

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_2D_MEMOIZATION_HPP_
