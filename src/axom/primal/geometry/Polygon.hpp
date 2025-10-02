// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Polygon.hpp
 *
 * \brief A Polygon primitive for primal
 */

#ifndef AXOM_PRIMAL_POLYGON_HPP_
#define AXOM_PRIMAL_POLYGON_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/StaticArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
namespace
{
static constexpr int DEFAULT_MAX_NUM_VERTICES = 8;
} /* end anonymous namespace */

/**
 * The polygon can have a dynamic or static array to store vertices
 */
enum class PolygonArray
{
  Dynamic,
  Static
};

// Forward declare the templated classes and operator functions
template <typename T, int NDIMS, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
class Polygon;

/// \brief Overloaded output operator for polygons
template <typename T, int NDIMS, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
std::ostream& operator<<(std::ostream& os, const Polygon<T, NDIMS, ARRAY_TYPE, MAX_VERTS>& poly);

/*!
 * \class Polygon
 *
 * \brief Represents a polygon defined by an array of points
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 * \tparam PolygonArray the array type of the polygon can be dynamic or
 *         static (default is dynamic). The static array type is for use in
 *         a device kernel.
 * \tparam MAX_VERTS the max number of vertices to preallocate space for
 *         a static array (default max is 8 vertices). MAX_VERTS is unused
 *         if array type is dynamic.
 * \note The polygon vertices should be ordered in a counter clockwise
 *       orientation with respect to the polygon's desired normal vector
 */
template <typename T, int NDIMS, PolygonArray ARRAY_TYPE = PolygonArray::Dynamic, int MAX_VERTS = DEFAULT_MAX_NUM_VERTICES>
class Polygon
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  // axom::Array for dynamic array type, StaticArray for static array type
  using ArrayType = std::conditional_t<ARRAY_TYPE == PolygonArray::Dynamic,
                                       axom::Array<PointType>,
                                       axom::StaticArray<PointType, MAX_VERTS>>;

public:
  /// Default constructor for an empty polygon (dynamic array specialization).
  /// Specializations are necessary to remove __host__ __device__ warning for
  /// axom::Array usage with the dynamic array type.
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Dynamic, int> = 0>
  Polygon()
  { }

  /// Default constructor for an empty polygon (static array specialization)
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Static, int> = 0>
  AXOM_HOST_DEVICE Polygon()
  { }

  /*!
   * \brief Destructor for Polygon. Suppress CUDA warnings for
   *        dynamic axom::Array.
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  ~Polygon() { m_vertices.clear(); }

  /*!
   * \brief Copy assignment operator for Polygon (static array specialization).
   *        Specializations are necessary to remove warnings.
   */
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Static, int> = 0>
  AXOM_HOST_DEVICE Polygon& operator=(const Polygon& other)
  {
    if(this == &other)
    {
      return *this;
    }

    m_vertices = other.m_vertices;
    return *this;
  }

  /*!
   * \brief Copy assignment operator for Polygon.
   *        (dynamic array specialization)
   */
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Dynamic, int> = 0>
  Polygon& operator=(const Polygon& other)
  {
    if(this == &other)
    {
      return *this;
    }

    m_vertices = other.m_vertices;
    return *this;
  }

  /// Copy constructor for Polygon. Specializations are necessary to
  /// remove __host__ __device__ warning for axom::Array usage with
  /// the dynamic array type.
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Dynamic, int> = 0>
  Polygon(const Polygon& other) : m_vertices(other.m_vertices)
  { }

  /// Copy constructor for Polygon (static array specialization)
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Static, int> = 0>
  AXOM_HOST_DEVICE Polygon(const Polygon& other) : m_vertices(other.m_vertices)
  { }

  /*!
   * \brief Constructor for an empty polygon that reserves space for
   *  the given number of vertices. Only available for dynamic arrays.
   *
   * \param [in] numExpectedVerts number of vertices for which to reserve space
   * \pre numExpectedVerts is not negative
   *
   */
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Dynamic, int> = 0>
  explicit Polygon(int numExpectedVerts)
  {
    SLIC_ASSERT(numExpectedVerts >= 0);
    m_vertices.reserve(numExpectedVerts);
  }

  /// \brief Constructor for a polygon with the given vertices
  explicit Polygon(const axom::Array<PointType>& vertices)
  {
    for(const auto& vertex : vertices)
    {
      m_vertices.push_back(vertex);
    }
  }

  // \brief Constructor for a polygon with an initializer list of Points
  //        (dynamic array specialization)
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Dynamic, int> = 0>
  explicit Polygon(std::initializer_list<PointType> vertices)
  {
    m_vertices = vertices;
  }

  // \brief Constructor for a polygon with an initializer list of Points
  //        (static array specialization)
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Static, int> = 0>
  AXOM_HOST_DEVICE explicit Polygon(std::initializer_list<PointType> vertices)
  {
    SLIC_ASSERT(static_cast<int>(vertices.size()) <= MAX_VERTS);

    for(const auto& vertex : vertices)
    {
      m_vertices.push_back(vertex);
    }
  }

  /// Return the number of vertices in the polygon
  AXOM_HOST_DEVICE
  int numVertices() const { return static_cast<int>(m_vertices.size()); }

  /*!
   * \brief Appends a vertex to the list of vertices
   *
   * \param [in] pt the point to be appended to the list of vertices
   *
   * \note If the array type is static and the list of vertices is full,
   *       addVertex will not modify the list of vertices.
   *
   * \sa axom::StaticArray::push_back() for behavior when array type is static
   *     and the list of vertices is full.
   */
  /// @{
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Static, int> = 0>
  AXOM_HOST_DEVICE void addVertex(const PointType& pt)
  {
    m_vertices.push_back(pt);
  }

  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Dynamic, int> = 0>
  AXOM_HOST_DEVICE void addVertex(const PointType& pt)
  {
#ifdef AXOM_DEVICE_CODE
    m_vertices.push_back_device(pt);
#else
    m_vertices.push_back(pt);
#endif
  }
  /// @}

  /// Clears the list of vertices (dynamic array specialization).
  /// Specializations are necessary to remove __host__ __device__ warning for
  /// axom::Array usage with the dynamic array type.
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Dynamic, int> = 0>
  void clear()
  {
    m_vertices.clear();
  }

  /// Clears the list of vertices (static array specialization)
  template <PolygonArray P_ARRAY_TYPE = ARRAY_TYPE,
            std::enable_if_t<P_ARRAY_TYPE == PolygonArray::Static, int> = 0>
  AXOM_HOST_DEVICE void clear()
  {
    m_vertices.clear();
  }

  /// Retrieves the vertex at index idx
  AXOM_HOST_DEVICE
  PointType& operator[](int idx) { return m_vertices[idx]; }

  /// Retrieves the vertex at index idx
  AXOM_HOST_DEVICE
  const PointType& operator[](int idx) const { return m_vertices[idx]; }

  /*! 
   * \brief Robustly returns the normal of the polygon (not normalized)
   * \pre This function is only valid when NDIMS = 3 and the polygon is valid
   * \return normal polygon normal vector
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, VectorType>::type normal() const
  {
    SLIC_ASSERT(isValid());
    const int nverts = numVertices();

    VectorType v0(m_vertices[nverts - 1]), v1(m_vertices[0]);

    VectorType normal = VectorType::cross_product(v0, v1);

    // Iterate over each pair of vertices
    for(int i = 1; i < nverts; ++i)
    {
      v0 = v1;
      v1 = VectorType(m_vertices[i]);
      normal += VectorType::cross_product(v0, v1);
    }

    return normal;
  }

  /// \brief Reverses orientation of the polygon in-place
  AXOM_HOST_DEVICE
  void reverseOrientation()
  {
    const int nverts = numVertices();
    const int mid = nverts >> 1;

    // Swap leftmost/rightmost vertices, midpoint unchanged
    for(int i = 0; i < mid; ++i)
    {
      const int left = i;
      const int right = nverts - i - 1;
      axom::utilities::swap(m_vertices[left], m_vertices[right]);
    }
  }

  /*!
   * \brief Computes the average of the polygon's vertex positions
   *
   * \return A point at the mean of the polygon's vertices
   * \pre  polygon.isValid() is true
   */
  AXOM_HOST_DEVICE
  PointType vertexMean() const
  {
    SLIC_ASSERT(isValid());

    PointType sum;

    const int sz = numVertices();
    for(int i = 0; i < sz; ++i)
    {
      sum.array() += m_vertices[i].array();
    }
    sum.array() /= sz;

    return sum;
  }

  /**
   * \brief Returns the area of the polygon (3D specialization)
   *
   * The algorithm sums up contributions from triangles defined by each edge 
   * of the polygon and a local origin (at the Polygon's \a vertexMean()).
   * The algorithm assumes a planar polygon embedded in 3D, however it will return 
   * a value for non-planar polygons.
   *
   * The area is always non-negative since 3D polygons do not have a unique orientation.
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, double>::type area() const
  {
    const int nVerts = numVertices();

    // check for early return
    if(nVerts < 3)
    {
      return 0.0;
    }

    // Add up areas of triangles connecting polygon edges the vertex average
    VectorType sum;
    const auto O = vertexMean();  // 'O' for (local) origin
    for(int curr = 0, prev = nVerts - 1; curr < nVerts; prev = curr++)
    {
      sum += VectorType::cross_product(m_vertices[prev] - O, m_vertices[curr] - O);
    }

    return 0.5 * axom::utilities::abs(sum.norm());
  }

  /**
   * \brief Returns the area of a planar polygon (2D specialization)
   *
   * The area is always non-negative.
   * \sa signedArea()
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type area() const
  {
    return axom::utilities::abs(signedArea());
  }

  /**
   * \brief Returns the signed area of a 2D polygon
   *
   * The signed area accounts for the orientation of the polygon.
   * It is positive when the vertices are oriented counter-clockwise
   * and negative when the vertices are oriented clockwise.
   * \note Signed area is only defined in 2D.
   * \sa area()
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type signedArea() const
  {
    const int nVerts = numVertices();
    double sum = 0.;

    // check for early return
    if(nVerts < 3)
    {
      return sum;
    }

    // use shoelace algorithm
    for(int curr = 0, prev = nVerts - 1; curr < nVerts; prev = curr++)
    {
      sum += m_vertices[prev][0] * m_vertices[curr][1]  //
        - m_vertices[curr][0] * m_vertices[prev][1];
    }

    return 0.5 * sum;
  }

  /*!
   * \brief Simple formatted print of a polygon instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int sz = numVertices();

    os << "{" << sz << "-gon:";
    for(int i = 0; i < sz - 1; ++i)
    {
      os << m_vertices[i] << ",";
    }
    if(sz >= 2)
    {
      os << m_vertices[sz - 1];
    }
    os << "}";

    return os;
  }

  /*!
   * \brief Simple check for validity of a polygon
   *
   * Initial check is that the polygon has three or more vertices
   * \return True, if the polygon is valid, False otherwise
   */
  AXOM_HOST_DEVICE
  bool isValid() const { return m_vertices.size() >= 3; }

private:
  ArrayType m_vertices {};
};

//------------------------------------------------------------------------------
/// Free functions implementing Polygon's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
std::ostream& operator<<(std::ostream& os, const Polygon<T, NDIMS, ARRAY_TYPE, MAX_VERTS>& poly)
{
  poly.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Polygon using fmt
template <typename T, int NDIMS, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
struct axom::fmt::formatter<axom::primal::Polygon<T, NDIMS, ARRAY_TYPE, MAX_VERTS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_POLYGON_HPP_
