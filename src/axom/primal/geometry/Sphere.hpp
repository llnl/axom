// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_SPHERE_HPP_
#define AXOM_PRIMAL_SPHERE_HPP_

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
class Sphere;

template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Sphere<T, NDIMS>& s);

/// @}

/*!
 * \class Sphere
 *
 * \brief Defines an oriented Sphere in 2-D (i.e., a circle) or 3-D given by
 *  its center, \f$ \mathcal{X} \f$ and radius \f$ \mathcal{R} \f$. The Sphere
 *  object provides associated operations on a sphere, such as, signed distance
 *  and orientation.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Sphere
{
public:
  using PointType = primal::Point<T, NDIMS>;

  static_assert((NDIMS == 2 || NDIMS == 3), "A Sphere object may be defined in 2-D or 3-D");

private:
  using VectorType = primal::Vector<T, NDIMS>;

public:
  /// \name Constructors
  /// @{

  /*!
   * \brief Constructs a Sphere centered at origin with the given radius.
   * \param [in] radius the radius of the Sphere (optional).
   * \note If a radius is not supplied, the default radius is 1.0.
   */
  AXOM_HOST_DEVICE
  explicit Sphere(T radius = 1.0) : m_center(PointType::zero()), m_radius(radius) { }

  /*!
   * \brief Constructs a Sphere with the given center and radius.
   *
   * \param [in] center user-supplied center.
   * \param [in] radius the radius of the Sphere (optional).
   *
   * \note If a radius is not supplied, the default radius is 1.0.
   */
  AXOM_HOST_DEVICE
  explicit Sphere(const PointType& center, T radius = 1.0) : m_center(center), m_radius(radius) { }

  /*!
   * \brief Constructs a Sphere with the given center and radius.
   *
   * \param [in] center user-supplied center.
   * \param [in] radius the radius of the Sphere (optional).
   *
   * \note If a radius is not supplied, the default radius is 1.0.
   *
   * \pre center != nullptr
   */
  AXOM_HOST_DEVICE
  explicit Sphere(const T* center, T radius = 1.0);

  /// @}

  /*!
   * \brief Returns the radius of the Sphere.
   * \return r the radius of the Sphere.
   */
  AXOM_HOST_DEVICE
  inline T getRadius() const { return m_radius; };

  /*!
   * \brief Returns the center of the Sphere.
   *
   * \return c pointer to array that holds the center of the Sphere.
   * \note c points to an array that is NDIMS long.
   * \post c != nullptr
   */
  AXOM_HOST_DEVICE
  inline const PointType& getCenter() const { return m_center; };

  /*!
   * \brief Returns the volume of the Sphere.
   */
  AXOM_HOST_DEVICE
  inline T getVolume() const { return 4.0 / 3 * M_PI * m_radius * m_radius * m_radius; };

  /*!
   * \brief Computes the signed distance of a point to the Sphere's boundary.
   *
   * \param [in] q The test point
   * \return d the computed signed distance of the point \a q to the sphere.
   *
   * \note The signed distance of a point \a q is:
   *  <ul>
   *   <li> negative inside the sphere </li>
   *   <li> positive outside the sphere </li>
   *   <li> zero on the boundary </li>
   *  </ul>
   */
  AXOM_HOST_DEVICE inline T computeSignedDistance(const PointType& q) const
  {
    return VectorType(m_center, q).norm() - m_radius;
  }

  /*!
   * \brief Computes the orientation of a point with respect to the Sphere.
   *
   * \param [in] q The test point
   * \param [in] TOL user-supplied tolerance. Optional. Default is 1.e-9.
   * \return orient the orientation of \a q with respect to the sphere.
   *
   *  \note This method returns one of the following values:
   *   <ul>
   *    <li> <b>ON_BOUNDARY</b>      : if \a q is on the sphere's boundary </li>
   *    <li> <b>ON_POSITIVE_SIDE</b> : if \a q is outside the sphere </li>
   *    <li> <b>ON_NEGATIVE_SIDE</b> : if \a q is inside the sphere </li>
   *  </ul>
   *
   * \see OrientationResult for the list of possible return values.
   *
   */
  AXOM_HOST_DEVICE
  inline int getOrientation(const PointType& q, double TOL = 1.e-9) const
  {
    const T signed_distance = this->computeSignedDistance(q);

    if(axom::utilities::isNearlyEqual(signed_distance, 0., TOL))
    {
      return primal::ON_BOUNDARY;
    }

    return (signed_distance < T {0}) ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
  }

  /*!
   * \brief Tests if this sphere instance intersects with another sphere.
   *
   * \param [in] sphere the sphere object to check for intersection
   * \param [in] TOL tolerance for intersection test. Optional. If not specified
   *  the default tolerance is set to 1.e-9.
   *
   * \return status true if the sphere intersects, false otherwise.
   */
  AXOM_HOST_DEVICE
  inline bool intersectsWith(const Sphere<T, NDIMS>& sphere, double TOL = 1.e-9) const;

  /*!
   * \brief Tests if this sphere completely contains another sphere.
   *
   * \param [in] other The sphere object to check for containment
   * \param [in] margin Amount to add to other's radius before comparing.
   *
   * Note: a sphere does contain itself.
   *
   * \return true if this sphere contains the other, false otherwise.
   */
  AXOM_HOST_DEVICE
  inline bool contains(const Sphere<T, NDIMS>& other, double margin = 0.0) const;

  /*!
   * \brief Prints the Sphere information in the given output stream.
   * \param [in,out] os the output stream to write to.
   * \note This method is primarily used for debugging.
   * \return s the modified output stream object.
   */
  std::ostream& print(std::ostream& os) const;

private:
  PointType m_center; /*!< sphere center */
  T m_radius;         /*!< sphere radius */
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
// Sphere Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Sphere<T, NDIMS>::Sphere(const T* center, T radius) : m_radius(radius)
{
  SLIC_ASSERT(center != nullptr);
  for(int i = 0; i < NDIMS; ++i)
  {
    m_center[i] = center[i];
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& Sphere<T, NDIMS>::print(std::ostream& os) const
{
  os << "{center: " << m_center << ", radius: " << m_radius << "}";
  return os;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline bool Sphere<T, NDIMS>::intersectsWith(const Sphere<T, NDIMS>& sphere,
                                                              double TOL) const
{
  const T distance_squared = VectorType(sphere.getCenter(), m_center).squared_norm();
  const T sum_of_radii = m_radius + sphere.getRadius();
  const T sum_of_radii_2 = sum_of_radii * sum_of_radii;

  return (distance_squared < sum_of_radii_2 ||
          utilities::isNearlyEqual(distance_squared, sum_of_radii_2, TOL));
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline bool Sphere<T, NDIMS>::contains(const Sphere<T, NDIMS>& other,
                                                        double margin) const
{
  const T center_sep = VectorType(other.getCenter(), m_center).norm();
  return (m_radius >= center_sep + other.getRadius() + margin);
}

//------------------------------------------------------------------------------
//  implementation of free functions
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Sphere<T, NDIMS>& s)
{
  return (s.print(os));
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Sphere using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Sphere<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_SPHERE_HPP_
