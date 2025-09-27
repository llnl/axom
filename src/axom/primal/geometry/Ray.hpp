// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_RAY_HPP_
#define AXOM_PRIMAL_RAY_HPP_

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
template <typename T, int DIM>
class Ray;

/*!
 * \brief Overloaded output operator for rays
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Ray<T, NDIMS>& ray);

/*!
 * \class
 *
 * \brief Represents a ray, \f$ R(t) \in \mathcal{R}^d \f$ , defined by an
 *  origin point, \f$ P \f$ and a normalized direction vector, \f$ \vec{d} \f$,
 *  \f$ \ni R(t)= P + t\vec{d} \forall t \ge 0 \f$
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Ray
{
public:
  using CoordType = T;
  using PointType = Point<T, NDIMS>;
  using SegmentType = Segment<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

public:
  /// Disable the default constructor
  Ray() = delete;

  /*!
   * \brief Constructs a ray object with the given origin and direction.
   * \param [in] origin the origin of the ray.
   * \param [in] direction the direction of the ray.
   * \pre direction.squared_norm()!= 0.0
   */
  AXOM_HOST_DEVICE
  Ray(const PointType& origin, const VectorType& direction);

  /*!
   * \brief Constructs a ray object from a directed segment.
   * \param [in] S user-supplied segment.
   */
  AXOM_HOST_DEVICE
  explicit Ray(const SegmentType& S);

  /*!
   * \brief Returns the point of origin of this Ray instance.
   * \return origin a point instance corresponding to the origin of the ray.
   */
  AXOM_HOST_DEVICE
  const PointType& origin() const { return m_origin; };

  /*!
   * \brief Returns a point along the ray by evaluating \f$ R(t) \f$
   * \param [in] t user-supplied value for R(t).
   * \return p a point along the ray.
   * \pre \f$ t \ge 0 \f$
   */
  AXOM_HOST_DEVICE
  PointType at(const T& t) const;

  /*!
   * \brief Returns the direction vector of this Ray instance.
   * \return direction the direction vector of the ray.
   * \post direction.norm()==1
   */
  AXOM_HOST_DEVICE
  const VectorType& direction() const { return m_direction; };

  /*!
   * \brief Simple formatted print of a ray instance
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
// Ray Implementation
//------------------------------------------------------------------------------

namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Ray<T, NDIMS>::Ray(const PointType& origin, const VectorType& direction)
  : m_origin(origin)
  , m_direction(direction.unitVector())
{
  SLIC_ASSERT(m_direction.squared_norm() != 0.0);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Ray<T, NDIMS>::Ray(const SegmentType& S)
  : m_origin(S.source())
  , m_direction(VectorType(S.source(), S.target()).unitVector())
{
  SLIC_ASSERT(m_direction.squared_norm() != 0.0);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> Ray<T, NDIMS>::at(const T& t) const
{
  PointType p;
  for(int i = 0; i < NDIMS; ++i)
  {
    p[i] = m_origin[i] + t * m_direction[i];
  }
  return (p);
}

//------------------------------------------------------------------------------
/// Free functions implementing Ray's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Ray<T, NDIMS>& ray)
{
  ray.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Ray using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Ray<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_RAY_HPP_
