// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_ZIP_INDEXABLE_HPP_
#define AXOM_PRIMAL_ZIP_INDEXABLE_HPP_

#include "axom/config.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
template <typename T>
class ZipBase
{
  static constexpr bool Exists = false;
};

}  // namespace detail

/*!
 * \class ZipIndexable
 *
 * \brief Utility class which allows iterating over Primal objects stored in
 *  memory as a structure of arrays.
 *
 *  \tparam GeomType the Primal instance type
 */
template <typename GeomType>
class ZipIndexable : detail::ZipBase<GeomType>
{
public:
  using CoordType = typename GeomType::CoordType;

  static_assert(detail::ZipBase<GeomType>::Exists,
                "A ZipBase specialization is not defined for this geometry "
                "type. Check that you are passing in a supported Primal type "
                "as a template parameter.");

  /**
   * \brief Inherited constructor from detail instance.
   */
  using detail::ZipBase<GeomType>::ZipBase;

  /**
   * \brief Inherited implementation of array index operator.
   */
  using detail::ZipBase<GeomType>::operator[];
};

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ZIP_INDEXABLE_HPP
