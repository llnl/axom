// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_POINT_HPP_
#define AXOM_PRIMAL_POINT_HPP_

#include "axom/core/NumericArray.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/matvecops.hpp"
#include "axom/slic/interface/slic.hpp"

// C/C++ includes
#include <cstring>
#include <ostream>
#include <type_traits>
#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class Point;

namespace detail
{
/// \name Type trait for Point class
///@{
template <typename T>
struct is_point : std::false_type
{ };

template <typename T, int NDIMS>
struct is_point<Point<T, NDIMS>> : std::true_type
{ };

template <typename T>
constexpr bool is_point_v = is_point<T>::value;
///@}

}  // namespace detail

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Equality comparison operator for points
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE bool operator==(const Point<T, NDIMS>& lhs, const Point<T, NDIMS>& rhs);

/*!
 * \brief Inequality comparison operator for points
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE bool operator!=(const Point<T, NDIMS>& lhs, const Point<T, NDIMS>& rhs);

/*!
 * \brief Overloaded output operator for points
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Point<T, NDIMS>& pt);

///@}

/*!
 * \accelerated
 * \class Point
 *
 * \brief The point class represents a point, \f$ p \in \mathcal{R}^d \f$ . It
 *  provides access methods to set and query the point coordinates.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Point
{
public:
  enum
  {
    DIMENSION = NDIMS,
    NBYTES = NDIMS * sizeof(T)
  };

  using NumericArray = axom::NumericArray<T, NDIMS>;
  using PointType = Point<T, NDIMS>;
  using CoordType = T;

public:
  /*!
   * \brief Fill the first sz coordinates with val and zeros the rest
   * \param [in] val The value to set the coordinates to.  Defaults to zero
   * \param [in] sz The number of coordinates to set to val.
   * The rest will be set to zero.  Defaults is NDIMS.
   * If sz is greater than NDIMS, we set all coordinates to val
   */
  AXOM_HOST_DEVICE
  explicit Point(T val = T(), int sz = NDIMS) : m_components(val, sz) { }

  /*!
   * \brief Constructor from a numeric array
   * \param [in] arr The numeric array to copy from
   */
  AXOM_HOST_DEVICE
  explicit Point(const NumericArray& arr) : m_components(arr) { }

  /*!
   * \brief Creates a point from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz num values to copy from the vals array. Defaults to NDIMS.
   * \note If sz is greater than NDIMS, we only take the first NDIMS values.
   */
  AXOM_HOST_DEVICE
  explicit Point(const T* vals, int sz = NDIMS) : m_components(vals, sz) { }

  /*!
   * \brief Creates a point from an initializer list
   * \param [in] values an initializer list containing the values of the
   * point. If the size is not the same as the size of this point, this
   * behaves the same way as the constructor which takes a pointer and size.
   */
  AXOM_HOST_DEVICE
  Point(std::initializer_list<T> values) : Point {values.begin(), static_cast<int>(values.size())}
  { }

  /*!
   * \brief Returns the dimension of this point instance.
   * \return d the dimension of the point.
   * \post d >= 1.
   */
  AXOM_HOST_DEVICE
  static constexpr int dimension() { return NDIMS; };

  /// \name Overloaded [] operator methods
  ///@{

  /*!
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return p[i] the value at the given component index.
   * \pre (i >= 0) && (i < ndims)
   */
  AXOM_HOST_DEVICE
  const T& operator[](int i) const { return m_components[i]; }

  AXOM_HOST_DEVICE
  T& operator[](int i) { return m_components[i]; }

  ///@}

  /// \name Raw data member access methods
  ///@{

  /*!
   * \brief Returns a pointer to the underlying data.
   */
  AXOM_HOST_DEVICE const T* data() const { return m_components.data(); }

  AXOM_HOST_DEVICE T* data() { return m_components.data(); }

  ///@}

  /*!
   * \brief Returns a reference to the underlying NumericArray.
   */
  AXOM_HOST_DEVICE
  const NumericArray& array() const { return m_components; }

  AXOM_HOST_DEVICE
  NumericArray& array() { return m_components; }

  /*!
   * \brief Output the point's coordinates to the array
   * \param arr The array that we are outputting to.
   * \pre The user needs to make sure that the array has been allocated
   * and has sufficient space for NDIMS coordinates.
   */
  AXOM_HOST_DEVICE
  void to_array(T* arr) const { m_components.to_array(arr); }

  /*!
   * \brief Equality comparison operator for points
   */
  AXOM_HOST_DEVICE
  friend inline bool operator==(const Point& lhs, const Point& rhs)
  {
    return lhs.m_components == rhs.m_components;
  }

  /*!
   * \brief Inequality operator for points
   */
  AXOM_HOST_DEVICE
  friend inline bool operator!=(const Point& lhs, const Point& rhs) { return !(lhs == rhs); }

  /*!
   * \brief Simple formatted print of a point instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const;

  /*!
   * \brief Utility function to constructs a Point with the given coordinates.
   * \param [in] x the x--coordinate of the point.
   * \param [in] y the y--coordinate of the point.
   * \param [in] z the z--coordinate of the point. Default is 0.0.
   * \return p a Point instance with the given coordinates.
   */
  AXOM_HOST_DEVICE
  static Point make_point(const T& x, const T& y, const T& z = 0.0);

  /*!
   * \brief Returns the midpoint between point A and point B.
   * \param [in] A user-supplied point
   * \param [in] B user-supplied point
   * \return p point at the midpoint A and B.
   */
  AXOM_HOST_DEVICE
  static Point midpoint(const Point& A, const Point& B);

  /*!
   * \brief Performs linear interpolation/extrapolation between two points.
   *
   * Given points A, B and a weight \f$ \alpha \f$, this method computes
   * point \f$ P\f$ defined by \f$ P = (1-\alpha) \cdot A + \alpha \cdot B\f$
   *
   * \param [in] A user-supplied point
   * \param [in] B user-supplied point
   * \param [in] alpha user-supplied scalar weight \f$ \alpha\f$
   *
   * \return P the computed point.
   *
   * \note \f$ \forall\alpha \in [0,1] \f$ this method linearly interpolates
   *  between the two points A, B.
   * \note \f$ \forall\alpha \not\in [0,1] \f$ the method extrapolates.
   *
   * \post \f$ P==A\f$ when \f$ \alpha=0.0\f$
   * \post \f$ P==B\f$ when \f$ \alpha=1.0\f$
   * \post The return point, P, and the user-supplied points A, B are collinear.
   */
  AXOM_HOST_DEVICE
  static Point lerp(const Point& A, const Point& B, T alpha);

  /*!
   * \brief Helper function to return a point whose coordinates are all 0
   */
  AXOM_HOST_DEVICE
  static Point zero() { return Point(); }

  /*!
   * \brief Helper function to return a point whose coordinates are all 1
   * This is equivalent to using the single value constructor: Point(1)
   * (with the appropriate casting) and is only valid for Points with
   * a numerical type (i.e. where static_cast<T>(1) is valid.
   */
  AXOM_HOST_DEVICE
  static Point ones() { return Point(static_cast<T>(1)); }

private:
  NumericArray m_components;
};

/// \name Pre-defined point types
/// @{

using Point2D = Point<double, 2>;
using Point3D = Point<double, 3>;

/// @}

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  Point implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> Point<T, NDIMS>::make_point(const T& x, const T& y, const T& z)
{
  T tmp_array[3] = {x, y, z};
  return Point(tmp_array, NDIMS);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> Point<T, NDIMS>::midpoint(const Point<T, NDIMS>& A,
                                                                  const Point<T, NDIMS>& B)
{
  Point<T, NDIMS> mid_point;

  for(int i = 0; i < NDIMS; ++i)
  {
    mid_point[i] = static_cast<T>(0.5 * (A[i] + B[i]));
  }

  return mid_point;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> Point<T, NDIMS>::lerp(const Point<T, NDIMS>& A,
                                                              const Point<T, NDIMS>& B,
                                                              T alpha)
{
  PointType res;
  const T beta = 1. - alpha;
  for(int i = 0; i < NDIMS; ++i)
  {
    res[i] = beta * A[i] + alpha * B[i];
  }
  return res;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& Point<T, NDIMS>::print(std::ostream& os) const
{
  os << "(";
  for(int dim = 0; dim < NDIMS - 1; ++dim)
  {
    os << static_cast<typename NonChar<T>::type>(m_components[dim]) << ",";
  }
  os << static_cast<typename NonChar<T>::type>(m_components[NDIMS - 1]) << ")";

  return os;
}

//------------------------------------------------------------------------------
/// Free functions implementing Point's operators
//------------------------------------------------------------------------------

template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Point<T, NDIMS>& pt)
{
  pt.print(os);
  return os;
}

/*!
 * \brief Transform a point using a transformation matrix.
 *
 * \param pt The point to be transformed.
 * \param transform The transformation matrix.
 *
 * \return The transformed point.
 *
 * \note We allow a matrix 1 rank larger than NDIMS so we can perform translations.
 */
template <typename T, int NDIMS>
Point<T, NDIMS> transform_point(
  const Point<T, NDIMS>& pt,
  const axom::numerics::Matrix<T>& transform = axom::numerics::Matrix<T>::identity(NDIMS))
{
  const int nr = transform.getNumRows();
  SLIC_ASSERT(nr == transform.getNumColumns());
  SLIC_ASSERT(nr == NDIMS || nr == (NDIMS + 1));
  SLIC_ASSERT(nr > 0 && nr <= 4);

  // Make a column vector to hold the point.
  T vec[4] = {0, 0, 0, 0};
  for(int row = 0; row < NDIMS; row++)
  {
    vec[row] = pt[row];
  }
  if(nr > NDIMS)
  {
    vec[nr - 1] = T {1};
  }

  // Transform the point.
  T transformedPt[4] = {0, 0, 0, 0};
  axom::numerics::matrix_vector_multiply(transform, vec, transformedPt);

  // Return the transformed point.
  return Point<T, NDIMS>(transformedPt);
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Point using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Point<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_POINT_HPP_
