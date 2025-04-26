// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_CONE_HPP_
#define AXOM_PRIMAL_CONE_HPP_

#include "axom/core.hpp"

#include "axom/primal/constants.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/fmt.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
/*!
  \class Cone

  \brief Represents a cone defined by a base radius,
  a top radius and the length from base to top along its axis.
  \tparam T the coordinate type, e.g., double, float, etc.
  \tparam NDIMS the number of spatial dimensions
*/
template <typename T, int NDIMS>
class Cone
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

public:
  /// \brief Default constructor constructs a degenerate cone.
  AXOM_HOST_DEVICE Cone()
    : m_baseRadius(0.0)
    , m_topRadius(0.0)
    , m_length(0.0)
    , m_base(0.0, NDIMS)
    , m_direction(0.0, NDIMS)
    { }

  /*!
    \brief Construct a cone with a base at the origin,
    oriented along the first axis.
    \param [in] baseRadius
    \param [in] topRadius
    \param [in] length
  */
  AXOM_HOST_DEVICE Cone(T baseRadius, T topRadius, T length)
    : m_baseRadius(baseRadius)
    , m_topRadius(topRadius)
    , m_length(length)
    , m_base(0.0, NDIMS)
    , m_direction(0.0, NDIMS)
    {
      m_direction[0] = 1.0;
      assertValid();
    }

  /*!
    \brief Construct a cone rotated to a given direction
    and translated.
    \param [in] baseRadius
    \param [in] topRadius
    \param [in] length
    \param [in] base Coordinates of the base
    \param [in] direction Direction of axis, from base to top.
   */
  AXOM_HOST_DEVICE Cone(T baseRadius, T topRadius, T length, const PointType& base, const VectorType& direction)
    : m_baseRadius(baseRadius)
    , m_topRadius(topRadius)
    , m_length(length)
    , m_base(base)
    , m_direction(direction.unitVector())
    {
      assertValid();
    }

  /*!
   * \brief Return the base coordinates.
   */
  AXOM_HOST_DEVICE PointType& getBase()
  {
    return m_base;
  }

  /*!
    \brief Return the base coordinates.
  */
  AXOM_HOST_DEVICE const PointType& getBase() const
  {
    return m_base;
  }

  /*!
    \brief Return the base radius.
  */
  AXOM_HOST_DEVICE T getBaseRadius()
  {
    return m_baseRadius;
  }

  /*!
    \brief Return the top radius.
  */
  AXOM_HOST_DEVICE T getTopRadius()
  {
    return m_topRadius;
  }

  /*!
    \brief Return the axis direction.
  */
  AXOM_HOST_DEVICE const VectorType& getDirection() const
  {
    return m_direction;
  }

  /*!
    \brief Simple formatted print of a cone instance
   \param os The output stream to write to
   \return A reference to the modified ostream
  */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_baseRadius << '-' << m_topRadius << "x" << m_length << " at " << m_base << " along " << m_direction << "}";

    return os;
  }

  /*!
   \brief Returns the volume of the cone

   Volume is only defined when NDIMS == 3.
  */
  AXOM_HOST_DEVICE
  template <int TDIM = NDIMS>
  typename std::enable_if<TDIM == 3, T>::type volume() const
  {
    T vol =
      (m_baseRadius*m_baseRadius + m_baseRadius*m_topRadius + m_topRadius*m_topRadius)
      * 1/3.0 * M_PI * m_length;
    return vol;
  }

private:
  T m_baseRadius;
  T m_topRadius;
  T m_length;
  PointType m_base;
  VectorType m_direction;

  AXOM_HOST_DEVICE void assertValid() const
  {
    SLIC_ASSERT(m_baseRadius >= 0.0);
    SLIC_ASSERT(m_topRadius >= 0.0);
    SLIC_ASSERT(m_length >= 0.0);
  }
};

//------------------------------------------------------------------------------
/// Free functions implementing Cone's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Cone<T, NDIMS>& Cone)
{
  Cone.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Cone using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Cone<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_CONE_HPP_
