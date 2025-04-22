// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_ZIP_BOUNDINGBOX_HPP_
#define AXOM_PRIMAL_ZIP_BOUNDINGBOX_HPP_

#include "axom/config.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/utils/ZipIndexable.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Implements ZipIndexable for a primal::BoundingBox instantiation
 */
template <typename T, int NDIMS>
struct ZipBase<BoundingBox<T, NDIMS>>
{
  using GeomType = BoundingBox<T, NDIMS>;

  static constexpr bool Exists = true;

  /// Default constructor for a ZipBase of primal::BoundingBox
  ZipBase()
  {
    for(int d = 0; d < NDIMS; ++d)
    {
      bb_min_arrays[d] = nullptr;
      bb_max_arrays[d] = nullptr;
    }
  }

  /*!
   * \brief Creates a ZipIndexable from a set of arrays
   * \param [in] min_arrays the arrays storing the min coordinate for each dimension
   * \param [in] max_arrays the arrays storing the max coordinate for each dimension
   *
   * \pre Size1 >= NDIMS
   * \pre Size2 >= NDIMS
   */
  template <size_t Size1, size_t Size2>
  ZipBase(const T* const (&min_arrays)[Size1], const T* const (&max_arrays)[Size2])
  {
    AXOM_STATIC_ASSERT_MSG(Size1 >= NDIMS, "Must provide at least NDIMS arrays");
    AXOM_STATIC_ASSERT_MSG(Size2 >= NDIMS, "Must provide at least NDIMS arrays");
    for(int d = 0; d < NDIMS; ++d)
    {
      bb_min_arrays[d] = min_arrays[d];
      bb_max_arrays[d] = max_arrays[d];
    }
  }

  /*!
   * \brief Returns the BoundingBox at an index i.
   * \param [in] i the index to access
   */
  AXOM_HOST_DEVICE GeomType operator[](int i) const
  {
    using PointType = typename GeomType::PointType;
    StackArray<T, NDIMS> min_data, max_data;
    for(int d = 0; d < NDIMS; ++d)
    {
      min_data[d] = bb_min_arrays[d][i];
      max_data[d] = bb_max_arrays[d][i];
    }
    return GeomType(PointType(min_data), PointType(max_data));
  }

private:
  const T* bb_min_arrays[NDIMS];
  const T* bb_max_arrays[NDIMS];
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ZIP_BOUNDINGBOX_HPP_
