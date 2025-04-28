// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file InletVector.hpp
 *
 * \brief This file contains the class definition of Inlet's InletVector class,
 * which wraps Primal's Vector class
 *******************************************************************************
 */

#ifndef INLET_INLETVECTOR_HPP
#define INLET_INLETVECTOR_HPP

#include "axom/primal/geometry/Vector.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \brief A wrapper over Primal's Vector3D that also includes dimension
 * information
 * 
 * Vector3D is a statically-sized (stack-allocated) three-element vector.
 * To represent two-element vectors using this type, additional dimension
 * information is required
 * 
 * \note To use vector operations on this class, perform them using the \p vec
 * member directly
 *******************************************************************************
 */
struct InletVector
{
  primal::Vector3D vec;
  int dim = 3;

  /*!
   *******************************************************************************
   * \brief Constructs an empty vector (size defaults to 3)
   *******************************************************************************
   */
  InletVector() = default;

  /*!
   *******************************************************************************
   * \brief Constructs a vector with an initializer list
   * 
   * \param [in] values The vector components to construct with
   *******************************************************************************
   */
  InletVector(std::initializer_list<double> values)
    : vec(values)
    , dim(static_cast<int>(values.size()))
  { }

  /*!
   *******************************************************************************
   * \brief Constructs a vector with an existing Primal vector and a dimension
   * 
   * \param [in] v The existing Primal vector
   * \param [in] d The dimension of the vector
   *******************************************************************************
   */
  explicit InletVector(primal::Vector3D&& v, int d = 3) : vec(std::move(v)), dim(d) { }

  /*!
   *******************************************************************************
   * \brief Constructs a vector with a pointer and a dimension
   * 
   * \param [in] values The pointer to the vector data
   * \param [in] d The dimension of the vector (length of the data)
   * 
   * \note Data is copied from the pointer - lifetime of the constructed InletVector
   * is not dependent on the lifetime of the pointer.
   *******************************************************************************
   */
  explicit InletVector(const double* values, int d = 3) : vec(values, d), dim(d) { }

  /*!
   *******************************************************************************
   * \brief Retrieves an element of the vector
   * 
   * \param [in] i The index of the element to retrieve (zero-indexed)
   *******************************************************************************
   */
  double operator[](const int i) const { return vec[i]; }
  /// \overload
  double& operator[](const int i) { return vec[i]; }

  /*!
   *******************************************************************************
   * \brief Retrieves the underlying Primal vector
   *******************************************************************************
   */
  operator axom::primal::Vector3D&() { return vec; }
  /// \overload
  operator const axom::primal::Vector3D&() const { return vec; }
};

/*!
 *******************************************************************************
 * \brief Compares two InletVectors
 * 
 * \param [in] u The comparison LHS
 * \param [in] v The comparison RHS
 * 
 * \return Whether the two vectors are equivalent
 *******************************************************************************
 */
inline bool operator==(const InletVector& u, const InletVector& v)
{
  return (u.vec == v.vec) && (u.dim == v.dim);
}

/*!
 *******************************************************************************
 * \brief "Prints" the vector to a stream
 * 
 * \param [inout] os The stream to insert into
 * \param [in] v The vector to insert into the stream
 *******************************************************************************
 */
inline std::ostream& operator<<(std::ostream& os, const InletVector& v)
{
  os << "<";
  for(int i = 0; i < v.dim - 1; i++)
  {
    os << v[i] << ",";
  }
  os << v[v.dim - 1] << ">";
  return os;
}

}  // end namespace inlet
}  // end namespace axom

/// Overload to format an inlet::InletVector using fmt
template <>
struct axom::fmt::formatter<axom::inlet::InletVector> : ostream_formatter
{ };

#endif  // INLET_INLETVECTOR_HPP
