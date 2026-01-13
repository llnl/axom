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

#include "axom/core/NumericLimits.hpp"
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
 * a top radius, the length, the orientation and the
 * location of the base center.
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 *
 * Negative length is allowed and leads to negative volume.
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
  static_assert(NDIMS >= 1);

public:
  /*!
   * \brief Default constructor constructs a cone with lengh 1,
   * base radius 1 and top radius 0, at the origin,
   * oriented along the first spatial direction.
   */
  AXOM_HOST_DEVICE Cone()
    : m_baseRad(1.0)
    , m_topRad(0.0)
    , m_length(1.0)
    , m_direction(0.0, NDIMS)
    , m_baseCenter(0.0, NDIMS)
  {
    m_direction[0] = 1.0;
  }

  /*!
   * \brief Construct a cone with a base center at the origin,
   * oriented along the first axis.
   * \param [in] baseRadius base radius
   * \param [in] topRadius top radius
   * \param [in] length Negative value is allowed and leads
   *   to negative volume.
   */
  AXOM_HOST_DEVICE Cone(T baseRadius, T topRadius, T length)
    : m_baseRad(baseRadius)
    , m_topRad(topRadius)
    , m_length(length)
    , m_direction(0.0, NDIMS)
    , m_baseCenter(0.0, NDIMS)
  {
    m_direction[0] = 1.0;
    assertValid();
  }

  /*!
   * \brief Construct a cone at an arbitrary position and orientation.
   * \param [in] baseRadius base radius
   * \param [in] topRadius top radius
   * \param [in] length length
   * \param [in] direction Direction of axis, from base to top.
   * \param [in] baseCenter Coordinates of the base center.
   */
  AXOM_HOST_DEVICE Cone(T baseRadius,
                        T topRadius,
                        T length,
                        const VectorType& direction,
                        const PointType& baseCenter)
    : m_baseRad(baseRadius)
    , m_topRad(topRadius)
    , m_length(length)
    , m_direction(direction.unitVector())
    , m_baseCenter(baseCenter)
  {
    assertValid();
  }

  //! \brief Return the radius at the base.
  AXOM_HOST_DEVICE T getBaseRadius() const { return m_baseRad; }

  //! \brief Return the radius at the top.
  AXOM_HOST_DEVICE T getTopRadius() const { return m_topRad; }

  //! \brief Return the length from base to top.
  AXOM_HOST_DEVICE T getLength() const { return m_length; }

  //! \brief Return the coordinates of the base center.
  AXOM_HOST_DEVICE const PointType& getBaseCenter() const { return m_baseCenter; }

  //! \brief Return the axis direction.
  AXOM_HOST_DEVICE const VectorType& getDirection() const { return m_direction; }

  /*!
   * \brief Return the interpolated/extrapolated radius at a given
   * distance from the base in the direction of \c getDirection().
   */
  AXOM_HOST_DEVICE double getRadiusAt(double z) const
  {
    if(std::abs(m_length) < axom::numeric_limits<double>::min())
    {
      return numeric_limits<T>::quiet_NaN();
    }
    return m_baseRad + (m_topRad - m_baseRad) / m_length * z;
  }

  /*!
   * \brief Simple formatted print of a cone instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "Cone {base radius" << m_baseRad << ", top radius " << m_topRad << ", axis at "
       << m_baseCenter << " along " << m_direction << '}';

    return os;
  }

  /*!
   * \brief Returns the signed volume of the cone,
   * which is negative if the length is negative.
   *
   * Volume is only defined when NDIMS == 3.
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, T>::type volume() const
  {
    T vol =
      (m_baseRad * m_baseRad + m_baseRad * m_topRad + m_topRad * m_topRad) / 3.0 * M_PI * m_length;
    return vol;
  }

private:
  T m_baseRad;
  T m_topRad;
  T m_length;
  VectorType m_direction;
  PointType m_baseCenter;

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
