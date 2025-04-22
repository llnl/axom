// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_ZIP_RAY_HPP_
#define AXOM_PRIMAL_ZIP_RAY_HPP_

#include "axom/config.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/utils/ZipIndexable.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Implements ZipIndexable for a primal::Ray instantiation
 */
template <typename T, int NDIMS>
struct ZipBase<Ray<T, NDIMS>>
{
  using GeomType = Ray<T, NDIMS>;

  static constexpr bool Exists = true;

  /// Default constructor for a ZipBase of primal::Ray
  ZipBase()
  {
    for(int d = 0; d < NDIMS; ++d)
    {
      ray_origs[d] = nullptr;
      ray_dirs[d] = nullptr;
    }
  }

  /*!
   * \brief Creates a ZipIndexable from a set of arrays
   * \param [in] orig_arrays the arrays for each dimension storing the origin
   *  of the ray
   * \param [in] dir_arrays the arrays storing for each dimension the direction
   *  of the ray
   *
   * \pre Size1 >= NDIMS
   * \pre Size2 >= NDIMS
   */
  template <size_t Size1, size_t Size2>
  ZipBase(const T* const (&orig_arrays)[Size1], const T* const (&dir_arrays)[Size2])
  {
    AXOM_STATIC_ASSERT_MSG(Size1 >= NDIMS, "Must provide at least NDIMS arrays");
    AXOM_STATIC_ASSERT_MSG(Size2 >= NDIMS, "Must provide at least NDIMS arrays");
    for(int d = 0; d < NDIMS; ++d)
    {
      ray_origs[d] = orig_arrays[d];
      ray_dirs[d] = dir_arrays[d];
    }
  }

  /*!
   * \brief Returns the Ray at an index i.
   * \param [in] i the index to access
   */
  AXOM_HOST_DEVICE GeomType operator[](int i) const
  {
    using PointType = typename GeomType::PointType;
    using VectorType = typename GeomType::VectorType;
    StackArray<T, NDIMS> orig_data, dir_data;
    for(int d = 0; d < NDIMS; ++d)
    {
      orig_data[d] = ray_origs[d][i];
      dir_data[d] = ray_dirs[d][i];
    }
    return GeomType(PointType(orig_data), VectorType(dir_data));
  }

private:
  const T* ray_origs[NDIMS];
  const T* ray_dirs[NDIMS];
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ZIP_RAY_HPP_
