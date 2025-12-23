// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file CurvedPolygon.hpp
 *
 * \brief A polygon primitive whose edges are Bezier curves
 */

#ifndef AXOM_PRIMAL_CURVEDPOLYGON_HPP_
#define AXOM_PRIMAL_CURVEDPOLYGON_HPP_

#include "axom/slic.hpp"

#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

// For NURBSCurveGWNCache objects
#include "axom/primal/operators/detail/winding_number_3d_memoization.hpp"

#include <vector>
#include <ostream>

namespace axom
{
namespace primal
{
namespace internal
{
///@{
/// \name Type traits for a compile-time flag for the cached object
template <typename U>
struct has_cached_data : std::false_type
{ };

template <typename T>
struct has_cached_data<detail::NURBSCurveGWNCache<T>> : std::true_type
{ };
///@}
}  // namespace internal

// Forward declare the templated classes and operator functions
template <typename CurveType>
class CurvedPolygon;

/*! \brief Overloaded output operator for polygons */
template <typename CurveType>
std::ostream& operator<<(std::ostream& os, const CurvedPolygon<CurveType>& poly);

/*!
 * \class CurvedPolygon
 *
 * \brief Represents a polygon with generic curves for edges
 * \tparam NDIMS the number of dimensions
 * \note The component curves should be ordered in a counter clockwise
 *       orientation with respect to the polygon's normal vector
 */
template <typename CurveType>
class CurvedPolygon
{
public:
  using NumericType = typename CurveType::NumericType;
  using PointType = typename CurveType::PointType;
  using VectorType = typename CurveType::VectorType;
  using BoundingBoxType = typename CurveType::BoundingBoxType;

public:
  /// Default constructor for an empty polygon
  CurvedPolygon() = default;

  /*!
   * \brief Constructor for an empty CurvedPolygon that reserves space for
   *  the given number of Edges
   *
   * \param [in] numExpectedEdges number of edges for which to reserve space
   * \pre numExpectedEdges is at least 1
   *
   */
  explicit CurvedPolygon(int nEdges)
  {
    SLIC_ASSERT(nEdges >= 1);
    m_edges.reserve(nEdges);
    m_edges.resize(nEdges);
  }

  /*!
   * \brief Copy constructor for another curve type
   */
  template <typename OtherCurveType>
  CurvedPolygon(const CurvedPolygon<OtherCurveType>& other_poly)
  {
    m_edges.resize(other_poly.numEdges());
    for(int e = 0; e < other_poly.numEdges(); ++e)
    {
      m_edges[e] = CurveType(other_poly[e]);
    }
  }

  ///! \brief Constructor for CurvedPolygon from an ArrayView of curves
  CurvedPolygon(axom::ArrayView<const CurveType> curves) : m_edges(curves) { }

  ///! \brief Constructor from a c-style array of \a nEdges curves
  CurvedPolygon(const CurveType* curves, int nEdges)
    : CurvedPolygon(axom::ArrayView<const CurveType>(curves, nEdges))
  { }

  ///! \brief Constructor from a Axom array of curves
  CurvedPolygon(const axom::Array<CurveType>& curves) : CurvedPolygon(curves.view()) { }

  /// Clears the list of edges
  void clear() { m_edges.clear(); }

  /// Returns true if the polygon has no edges
  bool empty() const { return m_edges.empty(); }

  /// \name Operations on edges
  /// @{

  /// Return the number of edges in the polygon
  int numEdges() const { return m_edges.size(); }

  void setNumEdges(int ngon)
  {
    SLIC_ASSERT(ngon >= 0);
    m_edges.resize(ngon);
  }

  /// Appends a curve to the list of edges
  void addEdge(const CurveType& c1) { m_edges.push_back(c1); }

  /// Consumes then appends a curve to the list of edges
  void addEdge(CurveType&& c1) { m_edges.push_back(std::move(c1)); }

  /// Splits an edge "in place"
  void splitEdge(int idx, NumericType t)
  {
    SLIC_ASSERT(idx >= 0 && idx < static_cast<int>(m_edges.size()));
    AXOM_STATIC_ASSERT_MSG(!internal::has_cached_data<CurveType>::value,
                           "splitEdge cannot be called on objects with associated cache data");

    m_edges.reserve(numEdges() + 1);
    m_edges.insert(m_edges.begin() + idx + 1, 1, m_edges[idx]);
    auto& csplit = m_edges[idx];
    csplit.split(t, m_edges[idx], m_edges[idx + 1]);
  }

  /// Returns a reference to the underlying array of curved edges
  axom::Array<CurveType>& getEdges() { return m_edges; }

  /// Returns a const reference to the underlying array of curved edges
  const axom::Array<CurveType>& getEdges() const { return m_edges; }

  /// @}

  /*! Retrieves the curve at index idx */
  CurveType& operator[](int idx) { return m_edges[idx]; }
  /*! Retrieves the curve at index idx */
  const CurveType& operator[](int idx) const { return m_edges[idx]; }

  /// Tests equality of two CurvedPolygons
  friend inline bool operator==(const CurvedPolygon<CurveType>& lhs,
                                const CurvedPolygon<CurveType>& rhs)
  {
    return lhs.m_edges == rhs.m_edges;
  }

  /// Tests inequality of two CurvedPolygons
  friend inline bool operator!=(const CurvedPolygon<CurveType>& lhs,
                                const CurvedPolygon<CurveType>& rhs)
  {
    return !(lhs == rhs);
  }

  /*!
   * \brief Simple formatted print of a CurvedPolygon instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int sz = numEdges();

    os << "{" << sz << "-sided curved polygon:";
    for(int i = 0; i < sz - 1; ++i)
    {
      os << m_edges[i] << ",";
    }
    if(sz >= 2)
    {
      os << m_edges[sz - 1];
    }
    os << "}";

    return os;
  }

  /*!
   * \brief Check closedness of a CurvedPolygon
   *
   * A CurvedPolygon is closed when the endpoint of each edge coincides with startpoint of next edge
   * \return \a true, if the polygon is closed, \a false otherwise
   */
  bool isClosed(double tol = 1e-5) const
  {
    using axom::utilities::isNearlyEqual;

    const double sq_tol = tol * tol;
    const int nEdges = numEdges();

    // initial basic check: no edges, or one edge or linear or quadratic order cannot be closed
    if(nEdges < 1 || (nEdges == 1 && m_edges[0].getNumControlPoints() <= 2))
    {
      return false;
    }

    // foreach edge: check last vertex of previous edge against first vertex of current edge
    for(int cur = 0, prev = nEdges - 1; cur < nEdges; prev = cur++)
    {
      const auto& lastPrev = m_edges[prev].getEndPoint();
      const auto& firstCur = m_edges[cur].getInitPoint();
      if(!isNearlyEqual(squared_distance(lastPrev, firstCur), 0., sq_tol))
      {
        return false;
      }
    }

    return true;
  }

  /// \brief Reverses orientation of a CurvedPolygon
  void reverseOrientation()
  {
    AXOM_STATIC_ASSERT_MSG(
      !internal::has_cached_data<CurveType>::value,
      "reverseOrientation cannot be called on objects with associated cache data");

    const int ngon = numEdges();
    const int mid = ngon >> 1;
    const bool isOdd = (ngon & 1) != 0;

    // Swap matching left/right cases, unmatched center is dealt with below
    for(int i = 0; i < mid; ++i)
    {
      const int left = i;
      const int right = ngon - i - 1;
      m_edges[left].reverseOrientation();
      m_edges[right].reverseOrientation();
      axom::utilities::swap(m_edges[left], m_edges[right]);
    }

    // Handle unmatched center curve, if necessary
    if(isOdd)
    {
      m_edges[mid].reverseOrientation();
    }
  }

  /// Returns an axis-aligned bounding box containing the CurvedPolygon
  BoundingBoxType boundingBox() const
  {
    BoundingBoxType bbox;
    for(const auto& cp : m_edges)
    {
      bbox.addBox(cp.boundingBox());
    }
    return bbox;
  }

private:
  axom::Array<CurveType> m_edges;
};

//------------------------------------------------------------------------------
/// Free functions implementing Polygon's operators
//------------------------------------------------------------------------------
template <typename CurveType>
std::ostream& operator<<(std::ostream& os, const CurvedPolygon<CurveType>& poly)
{
  poly.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CURVEDPOLYGON_HPP_
