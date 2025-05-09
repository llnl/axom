// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_ZIP_POINT_HPP_
#define AXOM_PRIMAL_ZIP_POINT_HPP_

#include "axom/config.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/utils/ZipIndexable.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Implements ZipIndexable for a primal::Point instantiation
 */
template <typename T, int NDIMS>
struct ZipBase<Point<T, NDIMS>>
{
  using GeomType = Point<T, NDIMS>;

  static constexpr bool Exists = true;

  /// Default constructor for a ZipBase of primal::Point
  ZipBase()
  {
    for(int d = 0; d < NDIMS; ++d)
    {
      pts_arrays[d] = nullptr;
    }
  }

  /*!
   * \brief Creates a ZipIndexable over a set of arrays.
   * \param [in] arrays the arrays storing coordinate data for each dimension
   *
   * \pre Size >= NDIMS
   */
  template <size_t Size>
  ZipBase(const T* const (&arrays)[Size])
  {
    AXOM_STATIC_ASSERT_MSG(Size >= NDIMS, "Must provide at least NDIMS arrays");
    for(int d = 0; d < NDIMS; ++d)
    {
      pts_arrays[d] = arrays[d];
    }
  }

  /*!
   * \brief Returns the Point at an index i.
   * \param [in] i the index to access
   */
  AXOM_HOST_DEVICE GeomType operator[](int i) const
  {
    StackArray<T, NDIMS> pt_data;
    for(int d = 0; d < NDIMS; ++d)
    {
      pt_data[d] = pts_arrays[d][i];
    }
    return GeomType(pt_data);
  }

private:
  const T* pts_arrays[NDIMS];
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ZIP_POINT_HPP_
