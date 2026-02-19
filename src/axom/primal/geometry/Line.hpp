// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_LINE_HPP_
#define AXOM_PRIMAL_LINE_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/slic/interface/slic.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class Line;

/*!
 * \brief Overloaded output operator for lines
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Line<T, NDIMS>& line);

/*!
 * \class Line
 *
 * \brief Represents a line, \f$ L(t) \in \mathcal{R}^d \f$ , defined by an
 *  origin point, \f$ P \f$ and a normalized direction vector, \f$ \vec{d} \f$,
 *  \f$ \ni L(t)= P + t\vec{d} \forall t \in \mathcal{R} \f$
 * 
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Line
{
public:
  using CoordType = T;
  using PointType = Point<T, NDIMS>;
  using SegmentType = Segment<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

public:
  // Disable the default constructor
  Line() = delete;

  /*!
     * \brief Constructs a line object with the given origin and direction.
     * \param [in] origin the origin of the line.
     * \param [in] direction the direction of the line.
     * \pre direction.squared_norm()!= 0.0
     */
  AXOM_HOST_DEVICE
  Line(const PointType& origin, const VectorType& direction);

  /*!
    * \brief Constructs a line object from a directed segment.
    * \params [in] S user-supplied segment 
    */
  explicit Line(const SegmentType& S);

  /*!
    * \brief Returns the origin of this Line instance.
    * \return origin a point instance corresponding to the origin of the line.
    */
  AXOM_HOST_DEVICE
  const PointType& origin() const { return m_origin; };

  /*!
    * \brief Returns a point along the line by evaluating \f$ L(t) \f$
    * \param [in] t user-supplied value for L(t).
    * \return p a point along the line.
    */
  AXOM_HOST_DEVICE
  PointType at(const T& t) const;

  /*!
    * \brief Returns the direction vector of this Line instance.
    * \return direction the direction vector of the line.
    * \post direction.norm()==1
    */
  AXOM_HOST_DEVICE
  const VectorType& direction() const { return m_direction; };

  /*!
    * \brief Simple formatted print of a line instance
    * \param os The output stream to write to
    * \return A reference to the modified ostream
    */
  std::ostream& print(std::ostream& os) const
  {
    os << "{origin:" << m_origin << "; direction:" << m_direction << "}";

    return os;
  }

private:
  PointType m_origin;
  VectorType m_direction;
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
// Line Implementation
//------------------------------------------------------------------------------

namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Line<T, NDIMS>::Line(const PointType& origin, const VectorType& direction)
  : m_origin(origin)
  , m_direction(direction.unitVector())
{
  SLIC_ASSERT(m_direction.squared_norm() != 0.0);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Line<T, NDIMS>::Line(const SegmentType& S)
  : m_origin(S.source())
  , m_direction(VectorType(S.source(), S.target()).unitVector())
{
  SLIC_ASSERT(m_direction.squared_norm() != 0.0);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> Line<T, NDIMS>::at(const T& t) const
{
  PointType p;
  for(int i = 0; i < NDIMS; ++i)
  {
    p[i] = m_origin[i] + t * m_direction[i];
  }
  return (p);
}

//------------------------------------------------------------------------------
/// Free functions implementing Line's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Line<T, NDIMS>& line)
{
  line.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Line using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Line<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_LINE_HPP_