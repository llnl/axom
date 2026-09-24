// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_ELLIPSOID_HPP_
#define AXOM_PRIMAL_ELLIPSOID_HPP_

#include "axom/core/Macros.hpp"

#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/matvecops.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/fmt.hpp"

#include <math.h>

namespace axom
{
namespace primal
{
/// \name Forward Declared Overloaded Operators
/// @{

template <typename T, int NDIMS>
class Ellipsoid;

template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Ellipsoid<T, NDIMS>& s);

/// @}

/*!
 * \class Ellipsoid
 *
 * \brief Defines an oriented Ellipsoid in 2-D (i.e., an ellipse) or 3-D 
 *  given by its center, \f$ \mathcal{X} \f$, and three points, 
 *  \f$ \mathcal{P_1, P_2, P_3} \f$.  These three points are the end points 
 *  of the axes of the ellipsoid and should be chosen so that the segments
 *  \f$ \mathcal{XP_1, XP_2, XP_3} \f$ are mutually perpendicular.  If this 
 *  is not the case, \f$ \mathcal{P_2} \f$ and \f$ \mathcal{P_3} \f$ will 
 *  be adjusted so the axes are mutually perpendicular.  The Ellipsoid
 *  object provides associated operations on an ellipsoid, such as signed 
 *  distance and orientation.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Ellipsoid
{
public:
  using PointType = primal::Point<T, NDIMS>;

  static_assert((NDIMS == 2 || NDIMS == 3), "An Ellipsoid object may be defined in 2-D or 3-D");

private:
  using VectorType = primal::Vector<T, NDIMS>;

public:
  /// \name Constructors
  /// @{

  /*!
   * \brief Constructs a Ellipsoid centered at origin with the given radius
   * in all dimensions: a circle or sphere.
   * \param [in] radius the radius of the Ellipsoid (optional).
   * \note If a radius is not supplied, the default radius is 1.0.
   */
  AXOM_HOST_DEVICE
  explicit Ellipsoid(T radius = 1.0)
    : m_center(PointType::zero())
    , m_P1(PointType::zero())
    , m_P2(PointType::zero())
    , m_P3(PointType::zero())
  {
    m_P1[0] = radius;
    m_P2[1] = radius;
    if constexpr (NDIMS == 3)
    {
      m_P3[2] = radius;
    }
  }

  /*!
   * \brief Constructs a Ellipsoid with the given center and radius
   * in all dimensions: a circle or sphere.
   *
   * \param [in] center user-supplied center.
   * \param [in] radius the radius of the Ellipsoid (optional).
   *
   * \note If a radius is not supplied, the default radius is 1.0.
   */
  AXOM_HOST_DEVICE
  explicit Ellipsoid(const PointType& center, T radius = 1.0)
    : m_center(center)
    , m_P1(PointType::zero())
    , m_P2(PointType::zero())
    , m_P3(PointType::zero())
  {
    m_P1[0] = radius;
    m_P1 = m_P1 + m_center;
    m_P2[1] = radius;
    m_P2 = m_P2 + m_center;
    if constexpr(NDIMS == 3)
    {
      m_P3[2] = radius;
      m_P3 = m_P3 + m_center;
    }
  }

  /*!
   * \brief Constructs a Ellipsoid with the given center and axis endpoints,
   * in two dimensions: an ellipse.
   *
   * \param [in] center user-supplied center.
   * \param [in] P1 the endpoint of the first semimajor axis
   * \param [in] P2 the endpoint of the second semimajor axis
   */
  AXOM_HOST_DEVICE
  explicit Ellipsoid(const PointType& center,
                     const PointType& P1,
                     const PointType& P2,
                     typename std::enable_if<TDIM == 2, bool>::type dummy = false)
    : m_center(center)
    , m_P1(P1)
    , m_P2(P2)
  { 
      orthogonalize();
  }

  /*!
   * \brief Constructs a Ellipsoid with the given center and axis endpoints,
   * in three dimensions: an ellipsoid.
   *
   * \param [in] center user-supplied center.
   * \param [in] P1 the endpoint of the first semimajor axis
   * \param [in] P2 the endpoint of the second semimajor axis
   * \param [in] P3 the endpoint of the last semimajor axis
   */
  AXOM_HOST_DEVICE
  explicit Ellipsoid(const PointType& center,
                     const PointType& P1,
                     const PointType& P2,
                     const PointType& P3,
                     typename std::enable_if<TDIM == 3, bool>::type dummy = false)
    : m_center(center)
    , m_P1(P1)
    , m_P2(P2)
    , m_P3(P3)
  {
    orthogonalize();
  }

  /*!
   * \brief If necessary, adjusts ellipsoid semi-axes to be orthogonal.
   */
  AXOM_HOST_DEVICE
  void orthogonalize() 
  {
    VectorType A(m_center, m_P1);
    VectorType B(m_center, m_P2);

    // If A dot B > tolerance, adjust B so that A dot B < tolerance

    if constexpr (NDIMS == 3)
    {
      VectorType C(m_center, m_P3);

      // VectorType maybe_C = A cross B
      // If necessary, negate maybe_C
      // If unit(C) dot unit(maybe_C) < 1 - tolerance, set C = maybe_C
    }
  }

  /*!
   * \brief Constructs a Ellipsoid with the given center and radius.
   *
   * \param [in] center user-supplied center.
   * \param [in] radius the radius of the Ellipsoid (optional).
   *
   * \note If a radius is not supplied, the default radius is 1.0.
   *
   * \pre center != nullptr
   */
  //AXOM_HOST_DEVICE
  //explicit Ellipsoid(const T* center, T radius = 1.0);

  /// @}

  /*!
   * \brief Returns the radius of the Ellipsoid.
   * \return r the radius of the Ellipsoid.
   */
  //AXOM_HOST_DEVICE
  //inline T getRadius() const { return m_radius; };

  /*!
   * \brief Returns the center of the Ellipsoid.
   *
   * \return c pointer to array that holds the center of the Ellipsoid.
   * \note c points to an array that is NDIMS long.
   * \post c != nullptr
   */
  AXOM_HOST_DEVICE
  inline const PointType& getCenter() const { return m_center; };

  AXOM_HOST_DEVICE
  inline void getRadii(T& a, T& b, T& c) const
  {
    a = VectorType(m_center, m_P1).norm();
    b = VectorType(m_center, m_P2).norm();
    if constexpr(NDIMS == 3)
    {
      c = VectorType(m_center, m_P3).norm();
    }
  }

  /*!
   * \brief Returns the volume of the Ellipsoid.
   */
  AXOM_HOST_DEVICE
  inline T getVolume() const 
  {
    T a, b, c;
    getRadii(a, b, c);
    return 4.0 / 3 * M_PI * a * b * c; 
  };

  /*!
   * \brief Computes the signed distance of a point to the Ellipsoid's boundary.
   *
   * \param [in] q The test point
   * \return d the computed signed distance of the point \a q to the Ellipsoid.
   *
   * \note The signed distance of a point \a q is:
   *  <ul>
   *   <li> negative inside the Ellipsoid </li>
   *   <li> positive outside the Ellipsoid </li>
   *   <li> zero on the boundary </li>
   *  </ul>
   */
  /*AXOM_HOST_DEVICE inline T computeSignedDistance(const PointType& q) const
  {
    return VectorType(m_center, q).norm() - m_radius;
  }*/

  /*!
   * \brief Computes the orientation of a point with respect to the Ellipsoid.
   *
   * \param [in] q The test point
   * \param [in] TOL user-supplied tolerance. Optional. Default is 1.e-9.
   * \return orient the orientation of \a q with respect to the Ellipsoid.
   *
   *  \note This method returns one of the following values:
   *   <ul>
   *    <li> <b>ON_BOUNDARY</b>      : if \a q is on the Ellipsoid's boundary </li>
   *    <li> <b>ON_POSITIVE_SIDE</b> : if \a q is outside the Ellipsoid </li>
   *    <li> <b>ON_NEGATIVE_SIDE</b> : if \a q is inside the Ellipsoid </li>
   *  </ul>
   *
   * \see OrientationResult for the list of possible return values.
   *
   */
  AXOM_HOST_DEVICE
  inline int getOrientation(const PointType& q, double TOL = 1.e-9) const
  {
    T a, b, c;
    getRadii(a, b, c);

    T testval = (q[0] * q[0]) / (a * a) + (q[1] * q[1]) / (b * b);
    if constexpr (NDIMS == 3)
    {
      testval += (q[2] * q[2]) / (c * c);
    }

    if(axom::utilities::isNearlyEqual(testval, T {1}, TOL))
    {
      return primal::ON_BOUNDARY;
    }

    return (signed_distance < T {1}) ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
  }

  /*!
   * \brief Tests if this Ellipsoid instance intersects with another Ellipsoid.
   *
   * \param [in] Ellipsoid the Ellipsoid object to check for intersection
   * \param [in] TOL tolerance for intersection test. Optional. If not specified
   *  the default tolerance is set to 1.e-9.
   *
   * \return status true if the Ellipsoid intersects, false otherwise.
   */
  /*AXOM_HOST_DEVICE
  inline bool intersectsWith(const Ellipsoid<T, NDIMS>& Ellipsoid, double TOL = 1.e-9) const;*/

  /*!
   * \brief Prints the Ellipsoid information in the given output stream.
   * \param [in,out] os the output stream to write to.
   * \note This method is primarily used for debugging.
   * \return s the modified output stream object.
   */
  std::ostream& print(std::ostream& os) const;

private:
  PointType m_center; /*!< Ellipsoid center */
  PointType m_P1, m_P2, m_P3;         /*!< Ellipsoid axis end-points */
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
// Ellipsoid Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Ellipsoid<T, NDIMS>::Ellipsoid(const T* center, T radius) : m_radius(radius)
{
  SLIC_ASSERT(center != nullptr);
  for(int i = 0; i < NDIMS; ++i)
  {
    m_center[i] = center[i];
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& Ellipsoid<T, NDIMS>::print(std::ostream& os) const
{
  os << "{center: " << m_center << ", radius: " << m_radius << "}";
  return os;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline bool Ellipsoid<T, NDIMS>::intersectsWith(const Ellipsoid<T, NDIMS>& Ellipsoid,
                                                              double TOL) const
{
  const T distance_squared = VectorType(Ellipsoid.getCenter(), m_center).squared_norm();
  const T sum_of_radii = m_radius + Ellipsoid.getRadius();
  const T sum_of_radii_2 = sum_of_radii * sum_of_radii;

  return (distance_squared < sum_of_radii_2 ||
          utilities::isNearlyEqual(distance_squared, sum_of_radii_2, TOL));
}

//------------------------------------------------------------------------------
//  implementation of free functions
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Ellipsoid<T, NDIMS>& s)
{
  return (s.print(os));
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Ellipsoid using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Ellipsoid<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_ELLIPSOID_HPP_
