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
#include <cstdint>
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
  const auto& getBoundingBox() const { return m_boundingBox; }
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
  using SubdivisionKey = std::uint64_t;
  using ChildSubdivisionData = std::pair<const BezierCurveData<T>*, const BezierCurveData<T>*>;
  using SubdivisionMap = std::unordered_map<SubdivisionKey, BezierCurveData<T>>;

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
    m_rootBezierData.resize(m_numSpans);
    auto beziers = a_curve.extractBezier();

    for(int idx = 0; idx < m_numSpans; ++idx)
    {
      m_bezierSubdivisionMaps[idx].reserve(32);
      m_rootBezierData[idx] = BezierCurveData<T>(beziers[idx], false, bbExpansionAmount);
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
      m_rootBezierData.clear();
      m_bezierSubdivisionMaps.clear();

      m_initPoint = m_endPoint = Point<T, 2> {0.0, 0.0};
    }
    else
    {
      m_numSpans = 1;
      m_rootBezierData.resize(1);
      m_bezierSubdivisionMaps.resize(1);

      m_initPoint = a_curve[0];
      m_endPoint = a_curve[m_degree];

      m_bezierSubdivisionMaps[0].reserve(32);
      m_rootBezierData[0] = BezierCurveData<T>(a_curve, false, bbExpansionAmount);
    }
  }

  /// \brief Query the map. If curve is not found, add it and its pair from subdivision
  const BezierCurveData<T>& getSubdivisionData(int idx,
                                               int refinementLevel,
                                               int refinementIndex,
                                               double bbExpansionAmount = 0.0) const
  {
    if(refinementLevel == 0)
    {
      return m_rootBezierData[idx];
    }

    auto& subdivision_map = m_bezierSubdivisionMaps[idx];
    const auto hash_key = makeSubdivisionKey(refinementLevel, refinementIndex);

    // If already there, return it
    if(auto it = subdivision_map.find(hash_key); it != subdivision_map.end())
    {
      return it->second;
    }

    // Otherwise, create (refinementLevel, refinementIndex) and sibling via their parent
    const int base = refinementIndex - (refinementIndex % 2);
    const BezierCurveData<T>* parent_data = nullptr;
    if(refinementLevel == 1)
    {
      parent_data = &m_rootBezierData[idx];
    }
    else
    {
      const auto parent_key = makeSubdivisionKey(refinementLevel - 1, refinementIndex / 2);
      auto parent_it = subdivision_map.find(parent_key);
      SLIC_ASSERT(parent_it != subdivision_map.end());
      parent_data = &parent_it->second;
    }

    const auto [child1, child2] = insertChildSubdivisionData(subdivision_map,
                                                             refinementLevel,
                                                             base,
                                                             *parent_data,
                                                             bbExpansionAmount);
    return (refinementIndex % 2 == 0) ? *child1 : *child2;
  }

  ChildSubdivisionData getSubdivisionChildren(int idx,
                                              int refinementLevel,
                                              int refinementIndex,
                                              const BezierCurveData<T>& supercurveData,
                                              double bbExpansionAmount = 0.0) const
  {
    auto& subdivision_map = m_bezierSubdivisionMaps[idx];
    const int child_level = refinementLevel + 1;
    const int child_base_index = 2 * refinementIndex;
    const auto key1 = makeSubdivisionKey(child_level, child_base_index);

    if(auto it1 = subdivision_map.find(key1); it1 != subdivision_map.end())
    {
      const auto key2 = makeSubdivisionKey(child_level, child_base_index + 1);
      auto it2 = subdivision_map.find(key2);
      SLIC_ASSERT(it2 != subdivision_map.end());
      return {&it1->second, &it2->second};
    }

    return insertChildSubdivisionData(subdivision_map,
                                      child_level,
                                      child_base_index,
                                      supercurveData,
                                      bbExpansionAmount);
  }

  ///@{
  //! \name Functions that mirror functionality of NURBSCurve and BezierCurve so signatures match in GWN evaluation.
  //!
  //! By limiting access to these functions, we ensure memoized information is always accurate
  auto getNumKnotSpans() const { return m_numSpans; }
  const auto& boundingBox() const { return m_boundingBox; }
  auto getNumControlPoints() const { return m_numControlPoints; }
  auto getDegree() const { return m_degree; }

  const auto& getInitPoint() const { return m_initPoint; }
  const auto& getEndPoint() const { return m_endPoint; }
  //@}

  friend bool operator==(const NURBSCurveGWNCache<T>& lhs, const NURBSCurveGWNCache<T>& rhs)
  {
    // numControlPoints, degree, and numSpans will be equal if the subdivision maps are
    return (lhs.m_rootBezierData == rhs.m_rootBezierData) &&
      (lhs.m_bezierSubdivisionMaps == rhs.m_bezierSubdivisionMaps) &&
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
      os << m_rootBezierData[0].getCurve();
    }
    for(int i = 1; i < m_numSpans; ++i)
    {
      os << ", " << m_rootBezierData[i].getCurve();
    }
    os << "}";

    return os;
  }

private:
  static SubdivisionKey makeSubdivisionKey(int refinementLevel, int refinementIndex)
  {
    return (static_cast<SubdivisionKey>(static_cast<std::uint32_t>(refinementLevel)) << 32) |
      static_cast<SubdivisionKey>(static_cast<std::uint32_t>(refinementIndex));
  }

  static ChildSubdivisionData insertChildSubdivisionData(SubdivisionMap& subdivision_map,
                                                         int refinementLevel,
                                                         int refinementIndex,
                                                         const BezierCurveData<T>& supercurveData,
                                                         double bbExpansionAmount)
  {
    BezierCurve<T, 2> sub1, sub2;
    supercurveData.getCurve().split(0.5, sub1, sub2);

    const auto key1 = makeSubdivisionKey(refinementLevel, refinementIndex);
    const auto key2 = makeSubdivisionKey(refinementLevel, refinementIndex + 1);

    auto it1 = subdivision_map
                 .try_emplace(key1, sub1, supercurveData.isConvexControlPolygon(), bbExpansionAmount)
                 .first;
    auto it2 = subdivision_map
                 .try_emplace(key2, sub2, supercurveData.isConvexControlPolygon(), bbExpansionAmount)
                 .first;

    return {&it1->second, &it2->second};
  }

  BoundingBox<T, 2> m_boundingBox;
  int m_numControlPoints;
  int m_degree;
  int m_numSpans;

  Point<T, 2> m_initPoint, m_endPoint;

  axom::Array<BezierCurveData<T>> m_rootBezierData;
  mutable axom::Array<SubdivisionMap> m_bezierSubdivisionMaps;
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
