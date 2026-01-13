// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_SEGMENT_HPP_
#define AXOM_PRIMAL_SEGMENT_HPP_

#include "axom/core/Macros.hpp"  // for Axom macros
#include "axom/slic.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int DIM>
class Segment;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Equality comparison operator for Segment
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE bool operator==(const Segment<T, NDIMS>& lhs, const Segment<T, NDIMS>& rhs);

/*!
 * \brief Inequality comparison operator for Segment
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE bool operator!=(const Segment<T, NDIMS>& lhs, const Segment<T, NDIMS>& rhs);

/*!
 * \brief Overloaded output operator for Segment
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Segment<T, NDIMS>& seg);

///@}

/*!
 * \class
 *
 * \brief Represents a directed straight line connecting two points
 *  \f$ A,B \in \mathcal{R}^n \f$ where point A is referred to as the source
 *  point and point B is the target.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Segment
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  enum
  {
    NUM_SEG_VERTS = 2
  };

public:
  Segment() = default;

  /*!
   * \brief Creates a segment instance from point A to point B.
   * \param A user-supplied source point
   * \param B user-supplied target point
   */
  AXOM_HOST_DEVICE Segment(const PointType& A, const PointType& B) : m_source(A), m_target(B) { };

  /*!
   * \brief Returns the source point of the segment.
   * \return s the source point of the segment.
   */
  AXOM_HOST_DEVICE const PointType& source() const { return m_source; };

  /*!
   * \brief Returns the target point of the segment.
   * \return t the target point of the segment.
   */
  AXOM_HOST_DEVICE const PointType& target() const { return m_target; };

  /*!
   * \brief Return the number of vertices in a Segment.
   *
   * \return The number of vertices in a Segment.
   */
  AXOM_HOST_DEVICE static constexpr int numVertices() { return NUM_SEG_VERTS; }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0 or 1
   */
  AXOM_HOST_DEVICE
  PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_SEG_VERTS);
    return idx == 0 ? m_source : m_target;
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0 or 1
   */
  AXOM_HOST_DEVICE
  const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_SEG_VERTS);
    return idx == 0 ? m_source : m_target;
  }

  /*!
   * \brief Returns a point \f$ (1 - t)A + tB \f$
   * \param [in] t user-supplied parameter for blending A and B
   * \post Return point P will always be collinear with the segment's
   *       start point A and end point B.
   * \post If \f$ t = 0, \f$ the return point \f$ P = A. \f$
   * \post If \f$ t = 1, \f$ the return point \f$ P = B. \f$
   */
  AXOM_HOST_DEVICE PointType at(const T& t) const { return PointType::lerp(m_source, m_target, t); }

  /*!
   * \brief Returns the length of the segment
   */
  AXOM_HOST_DEVICE double length() const { return VectorType(m_source, m_target).norm(); }

  /*!
   * \brief Returns a vector normal to the segment
   *
   * \note Only available in 2D
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, VectorType>::type normal() const
  {
    return VectorType {m_target[1] - m_source[1], m_source[0] - m_target[0]};
  }

  /*!
   * \brief Equality comparison operator for segments
   */
  AXOM_HOST_DEVICE
  friend inline bool operator==(const Segment& lhs, const Segment& rhs)
  {
    return lhs.m_source == rhs.m_source && lhs.m_target == rhs.m_target;
  }

  /*!
   * \brief Inequality operator for segments
   */
  AXOM_HOST_DEVICE
  friend inline bool operator!=(const Segment& lhs, const Segment& rhs) { return !(lhs == rhs); }

  /*!
   * \brief Simple formatted print of a segment instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{source:" << m_source << "; target:" << m_target << "}";

    return os;
  }

private:
  PointType m_source;
  PointType m_target;
};

}  // namespace primal
}  // namespace axom

//------------------------------------------------------------------------------
//  Segment Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
/// Free functions implementing Segments's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Segment<T, NDIMS>& seg)
{
  seg.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Segment using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Segment<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_SEGMENT_HPP_
