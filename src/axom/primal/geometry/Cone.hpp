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
 * \class Cone
 *
 * \brief Represents a cone defined by a base radius,
 * a top radius, the length from base to top along its axis,
 * the orientation of the axis and a translation.
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 *
 * A negative length is allowed and results in a negative volume.
 * Radii must be non-negative.
 *
 * A cylinder can be represented using equal base and top radii.
 */
template <typename T, int NDIMS>
class Cone
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

public:
  /*!
   * \brief Default constructor constructs a degenerate cone,
   * oriented along the first axis.
   */
  AXOM_HOST_DEVICE Cone()
    : m_baseZ(0.0)
    , m_baseRad(0.0)
    , m_topZ(0.0)
    , m_topRad(0.0)
    , m_direction(0.0, NDIMS)
    , m_origin(0.0, NDIMS)
  {
    m_direction[0] = 1.0;
  }

  /*!
   * \brief Construct a cone with a base at the origin,
   * oriented along the first axis.
   * \param [in] baseRadius
   * \param [in] topRadius
   * \param [in] length Negative value is allowed and leads
   *   to negative volume.
   */
  AXOM_HOST_DEVICE Cone(T baseRadius, T topRadius, T length)
    : m_baseZ(0.0)
    , m_baseRad(baseRadius)
    , m_topZ(length)
    , m_topRad(topRadius)
    , m_direction(0.0, NDIMS)
  {
    m_direction[0] = 1.0;
    assertValid();
  }

  /*!
   * \brief Construct a cone at an arbitrary position and orientation.
   * \param [in] baseZ
   * \param [in] baseRadius
   * \param [in] topZ
   * \param [in] topRadius
   * \param [in] direction Direction of axis, from base to top.
   * \param [in] origin Coordinates of the z=0 point
   *
   * topZ < baseZ is allowed and leads to negative length and volume.
   */
  AXOM_HOST_DEVICE Cone(T baseZ,
                        T baseRadius,
                        T topZ,
                        T topRadius,
                        const VectorType& direction,
                        const PointType& origin)
    : m_baseZ(baseZ)
    , m_baseRad(baseRadius)
    , m_topZ(topZ)
    , m_topRad(topRadius)
    , m_direction(direction.unitVector())
    , m_origin(origin)
  {
    assertValid();
  }

  //! \brief Return the z-coordinate of the base.
  AXOM_HOST_DEVICE T getBaseZ() const { return m_baseZ; }

  //! \brief Return the radius at the base.
  AXOM_HOST_DEVICE T getBaseRadius() const { return m_baseRad; }

  //! \brief Return the z-coordinate of the top.
  AXOM_HOST_DEVICE T getTopZ() const { return m_topZ; }

  //! \brief Return the radius at the top.
  AXOM_HOST_DEVICE T getTopRadius() const { return m_topRad; }

  //! \brief Return algebraic length of the cone.
  AXOM_HOST_DEVICE T getLength() const { return m_topZ - m_baseZ; }

  //! \brief Return the coordinates at the (z=0, r=0) point
  AXOM_HOST_DEVICE const PointType& getOrigin() const { return m_origin; }

  //! \brief Return the axis direction.
  AXOM_HOST_DEVICE const VectorType& getDirection() const { return m_direction; }

  //! \brief Return the interpolated/extrapolated radius at a given z.
  AXOM_HOST_DEVICE double getRadiusAt(double z) const
  {
    return m_baseRad + (m_topRad - m_baseRad)/(m_topZ - m_baseZ) * (z - m_baseZ);
  }

  /*!
   * \brief Simple formatted print of a cone instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "Cone{ base(" << m_baseZ << ',' << m_baseRad << "), top(" << m_topZ << ',' << m_topRad << ", axis at " << m_origin << " along " << m_direction << '}';

    return os;
  }

  /*!
   * \brief Returns the algebraic volume of the cone
   *
   * The volume returned is positive when the top-z coordinate is larger
   * than the base-z coordinate.  Otherwise, it's non-positive.
   *
   * Volume is only defined when NDIMS == 3.
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<TDIM == 3, T>::type volume() const
  {
    T vol = (m_baseRad * m_baseRad + m_baseRad * m_topRad + m_topRad * m_topRad) / 3.0 * M_PI * (m_topZ - m_baseZ);
    return vol;
  }

private:
  T m_baseZ;
  T m_baseRad;
  T m_topZ;
  T m_topRad;
  VectorType m_direction;
  PointType m_origin;

  AXOM_HOST_DEVICE void assertValid() const
  {
    SLIC_ASSERT(m_baseRad >= 0.0);
    SLIC_ASSERT(m_topRad >= 0.0);
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
