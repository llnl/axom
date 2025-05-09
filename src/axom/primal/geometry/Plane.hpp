// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_PLANE_HPP_
#define AXOM_PRIMAL_PLANE_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/numerics/matvecops.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Point.hpp"

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace primal
{
/// \name Forward Declared Overloaded Operators
/// @{

template <typename T, int NDIMS>
class Plane;

template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Plane<T, NDIMS>& p);

/// @}

/// \name Forward-declared helper functions
/// @{

/*!
 * \brief Constructs a Plane in 2D that goes through the specified points.
 *
 * \param [in] x1 coordinates of the first point.
 * \param [in] x2 coordinates of the second point.
 *
 * \return plane the plane/line defined by the two points
 *
 * \pre x1 should not be equal to x2
 */
template <typename T>
AXOM_HOST_DEVICE Plane<T, 2> make_plane(const Point<T, 2>& x1, const Point<T, 2>& x2);

/*!
 * \brief Constructs a Plane in 3D that goes through the specified points.
 *
 * \param [in] x1 coordinates of the first point.
 * \param [in] x2 coordinates of the second point.
 * \param [in] x3 coordinates of the third point.
 *
 * \return plane the plane defined by the three points
 *
 * \pre The user-supplied points, x1, x2, x3 should not be collinear.
 */
template <typename T>
AXOM_HOST_DEVICE Plane<T, 3> make_plane(const Point<T, 3>& x1,
                                        const Point<T, 3>& x2,
                                        const Point<T, 3>& x3);

/// @}

/*!
 * \class Plane
 *
 * \brief Defines an oriented plane within a 2-D or 3-D Euclidean space
 *  and provides associated operators, such as, projection, signed distance,
 *  orientation, etc.
 *
 *  The Plane object defines an oriented plane in Hessian Normal form:
 *
 *           \f$ \mathcal{N} \cdot x - d = 0 \f$
 *
 *  where, \f$ \mathcal{N} \f$ is a normal and \f$ d \f$ is the offset of
 *  the plane (in the direction of the specified normal) to the origin.
 *
 *  The Plane class defines a co-dimension-one plane that splits the ambient
 *  space into two halfspaces, such that, the signed distance for any point
 *  is: (a) negative below the plance, (b) positive above the plane and (c)
 *  zero on the plane.
 *
 *  A Plane object may be constructed in three ways:
 *  1. Specifying a normal and a point that the plane passes through
 *  2. Specifying a normal and an offset from the origin, or
 *  3. Specifying the points that the plane passes through. A minimum of two
 *     points are required to define a plane in 2D, or, three points for 3D.
 *
 * \tparam T the underlying data type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * \pre NDIMS==2 || NDIMS==3
 */
template <typename T, int NDIMS>
class Plane
{
public:
  using VectorType = primal::Vector<T, NDIMS>;
  using PointType = primal::Point<T, NDIMS>;

  static_assert(NDIMS == 2 || NDIMS == 3, "A plane object may be defined in 2-D or 3-D");

public:
  /// \name Constructors
  /// @{

  /*!
   * \brief Constructs a Plane with a zero normal vector and zero offset
   *
   * \post isValid() == false
   */
  Plane() = default;

  /*!
   * \brief Constructs a Plane with a given normal, \f$ \mathcal{N} \f$, that
   *  passes through the specified point, \f$  x \f$
   *
   * \param [in] normal the supplied plane normal, \f$ \mathcal{N} \f$
   * \param [in] x a specified point that the plane passes through, \f$ x \f$
   * \param [in] normalize if true, normalize the given vector (default: true)
   */
  AXOM_HOST_DEVICE Plane(const VectorType& normal, const PointType& x, bool normalize = true);

  /*!
   * \brief Constructs a plane with a specified normal, \f$ \mathcal{N} \f$,
   *  located at the given offset from origin.
   *
   * \param [in] normal the supplied plane normal, \f$ \mathcal{N} \f$
   * \param [in] offset the plane offset from origin
   * \param [in] normalize if true, normalize the given vector (default: true)
   */
  AXOM_HOST_DEVICE Plane(const VectorType& normal, T offset, bool normalize = true);

  /// @}

  /*!
   * \brief Returns the dimension of the Plane
   * \return dim the dimension of the Plane.
   * \post (dim==2) || (dim==3)
   */
  AXOM_HOST_DEVICE static constexpr int getDimension() { return NDIMS; };

  /*!
   * \brief Returns a const pointer to the plane's normal vector.
   * \return N the normal vector.
   */
  AXOM_HOST_DEVICE const VectorType& getNormal() const { return m_normal; }

  /*!
   * \brief Returns the offset of the plane from origin.
   * \return offSet the offset of the plane from origin.
   */
  AXOM_HOST_DEVICE T getOffset() const { return m_offset; }

  /*!
   * \brief Computes the signed distance of the specified point to this Plane.
   *
   * \param [in] x buffer consisting of the query point coordinates.
   * \return d the signed distance of the point x to the plane.
   */
  AXOM_HOST_DEVICE T signedDistance(const PointType& x) const
  {
    return m_normal.dot(VectorType(x)) - m_offset;
  }

  /*!
   * \brief Computes the projection of a given point, x, onto this Plane.
   *
   * \param [in] x buffer consisting of the coordinates of the point to project.
   * \return projx the coordinates of the projected point.
   *
   * \post this->getOrientedSide( projx ) == ON_BOUNDARY
   */
  AXOM_HOST_DEVICE PointType projectPoint(const PointType& x) const;

  /*!
   * \brief Computes the reflection of a given point, x, across this Plane.
   *
   * \param [in] x buffer consisting of the coordinates of the point to project.
   * \return refx the coordinates of the reflected point.
   *
   * \post The reflected point will be on the Plane if x is on the Plane,
   *       otherwise it will be on the opposite side of the Plane from x.
   */
  AXOM_HOST_DEVICE PointType reflectPoint(const PointType& x) const;

  /*!
   * \brief Flips the orientation of the plane.
   */
  AXOM_HOST_DEVICE void flip();

  /*!
   * \brief Computes the orientation of the point with respect to this Plane.
   *
   * \param [in] x buffer consisting of the query point coordinates.
   * \param [in] TOL user-supplied tolerance. Optional. Default is 1.e-9.
   * \return orientation the orientation of the point
   *
   * \note This method returns one of the following values:
   *  <ul>
   *   <li> <b>ON_BOUNDARY</b>      : if `x` lies on the plane  </li>
   *   <li> <b>ON_POSITIVE_SIDE</b> : if `x` is above the plane </li>
   *   <li> <b>ON_NEGATIVE_SIDE</b> : if `x` is below the plane </li>
   *  </ul>
   *
   * \see OrientationResult
   */
  AXOM_HOST_DEVICE int getOrientation(const PointType& x, double TOL = 1.e-9) const;

  /*!
   * \brief Simple check for validity of a Plane
   *
   * Check that the normal is not the zero vector.
   *
   * \param [in] TOL user-supplied tolerance. Optional. Default is 1.0e-50.
   *
   * \return True, if the Plane is valid, False otherwise
   */
  AXOM_HOST_DEVICE bool isValid(double TOL = PRIMAL_TINY) const;

  /*!
   * \brief Prints the Plane information in the given output stream.
   * \param [in] os The output stream to write to.
   * \note This method is primarily used for debugging.
   * \return s the modified output stream object.
   */
  std::ostream& print(std::ostream& os) const;

private:
  /*!
   * \brief Sets the normal of this Plane instance.
   *
   * \param [in] normal pointer to a buffer consisting of the normal
   * \note The supplied buffer must be at least NDIMS long.
   * \param [in] normalize if true, normalize the given vector
   */
  AXOM_HOST_DEVICE void setNormal(const VectorType& normal, bool normalize = true);

  VectorType m_normal;              /*!< plane normal vector (defaults to the zero vector) */
  T m_offset {static_cast<T>(0.0)}; /*!< offset from origin */
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
//  IMPLEMENTATION
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Plane<T, NDIMS>::Plane(const VectorType& normal, const PointType& x, bool normalize)
{
  SLIC_ASSERT_MSG(!normal.is_zero(), "Normal vector of a plane should be non-zero");

  this->setNormal(normal, normalize);

  m_offset = m_normal.dot(VectorType(x.array()));
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Plane<T, NDIMS>::Plane(const VectorType& normal, T offset, bool normalize)
  : m_offset(offset)
{
  SLIC_ASSERT_MSG(!normal.is_zero(), "Normal vector of a plane should be non-zero");

  this->setNormal(normal, normalize);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE typename Plane<T, NDIMS>::PointType Plane<T, NDIMS>::projectPoint(const PointType& x) const
{
  const T signed_distance = this->signedDistance(x);
  return x - signed_distance * m_normal;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE typename Plane<T, NDIMS>::PointType Plane<T, NDIMS>::reflectPoint(const PointType& x) const
{
  const T signed_distance = this->signedDistance(x);
  return x - static_cast<T>(2.0) * signed_distance * m_normal;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE void Plane<T, NDIMS>::flip()
{
  m_normal *= static_cast<T>(-1.0);
  m_offset *= static_cast<T>(-1.0);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE int Plane<T, NDIMS>::getOrientation(const PointType& x, double TOL) const
{
  const T signed_distance = this->signedDistance(x);

  if(utilities::isNearlyEqual(signed_distance, static_cast<T>(0.0), static_cast<T>(TOL)))
  {
    return primal::ON_BOUNDARY;
  }
  return signed_distance < static_cast<T>(0.0) ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
}

template <typename T, int NDIMS>
AXOM_HOST_DEVICE bool Plane<T, NDIMS>::isValid(double TOL) const
{
  for(int i = 0; i < NDIMS; ++i)
  {
    if(!utilities::isNearlyEqual(m_normal[i], static_cast<T>(0.0), static_cast<T>(TOL)))
    {
      return true;
    }
  }

  return false;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& Plane<T, NDIMS>::print(std::ostream& os) const
{
  os << "normal: [ ";
  for(int i = 0; i < NDIMS; ++i)
  {
    os << m_normal[i] << " ";
  }
  os << "] offset: " << m_offset;
  return (os);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE void Plane<T, NDIMS>::setNormal(const VectorType& normal, bool normalize)
{
  if(normalize)
  {
    m_normal = normal.unitVector();
  }
  else
  {
    m_normal = normal;
  }
}

//------------------------------------------------------------------------------
//  implementation of free functions
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Plane<T, NDIMS>& p)
{
  return (p.print(os));
}

template <typename T>
AXOM_HOST_DEVICE Plane<T, 2> make_plane(const Point<T, 2>& x1, const Point<T, 2>& x2)
{
  Vector<T, 2> normal;
  normal[0] = x1[1] - x2[1];  // -dy
  normal[1] = x2[0] - x1[0];  //  dx

  // check for degenerate line or triangle
  SLIC_CHECK_MSG(!normal.is_zero(), "Supplied points form a degenerate line");

  return Plane<T, 2>(normal, x1);
}

//------------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE Plane<T, 3> make_plane(const Point<T, 3>& x1,
                                        const Point<T, 3>& x2,
                                        const Point<T, 3>& x3)
{
  Vector<T, 3> r1(x1, x2);
  Vector<T, 3> r2(x1, x3);

  Vector<T, 3> normal = Vector<T, 3>::cross_product(r1, r2);

  // check for degenerate line or triangle
  SLIC_CHECK_MSG(!normal.is_zero(), "Supplied points form a degenerate triangle");

  return Plane<T, 3>(normal, x1);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_PLANE_HPP_
