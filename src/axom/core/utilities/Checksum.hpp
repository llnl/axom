// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_UTILITIES_CHECKSUM_HPP_
#define AXOM_UTILITIES_CHECKSUM_HPP_

#include <cmath>

#include "axom/config.hpp"  // for compile-time definitions
#include "axom/core/ArrayView.hpp"

namespace axom
{
namespace utilities
{
// Checksum types
using CheckSum = long double;
using ScaleFactor = double;

/*!
 * \brief Calculate and return checksum for data arrays.
 *
 * \param view The view that contains the data.
 *
 * \note Adapted from RAJAPerf at https://github.com/LLNL/RAJAPerf/blob/cda42470851fff2b7c8e6a9b5b11ab83f33a5a07/src/common/DataUtils.cpp#L598-L622
 *
 * \return A CheckSum value for the array view.
 */
template <typename DataGetter>
inline CheckSum calculateChecksum(DataGetter data, axom::IndexType len)
{
  CheckSum tchk = 0.0;
  CheckSum ckahan = 0.0;
  for (axom::IndexType j = 0; j < len; ++j) {
    const auto value = data(j);
    CheckSum x = (std::abs(std::sin(j+1.0))+0.5) * value;
    CheckSum y = x - ckahan;
    volatile CheckSum t = tchk + y;
    volatile CheckSum z = t - tchk;
    ckahan = z - y;
    tchk = t;
  }
  return tchk;
}

/*!
 * \brief Calculate and return checksum.
 *
 * \param value The value we want to checksum.
 *
 * \return A CheckSum value for the array view.
 */
template <typename T>
inline CheckSum checksum(T value, const ScaleFactor scaleFactor = ScaleFactor{1})
{
  return calculateChecksum([=](axom::IndexType) { return static_cast<CheckSum>(value); }, 1) * scaleFactor;
}

/*!
 * \brief Calculate and return checksum for an array view.
 *
 * \param view The view that contains the data we want to checksum.
 *
 * \return A CheckSum value for the array view.
 */
template <typename T>
inline CheckSum checksum(axom::ArrayView<T> view, const ScaleFactor scaleFactor = ScaleFactor{1})
{
  return calculateChecksum([=](axom::IndexType i) { return static_cast<CheckSum>(view[i]); }, view.size()) * scaleFactor;
}

}  // namespace utilities
}  // namespace axom

#endif  // AXOM_UTILITIES_CHECKSUM_HPP_
