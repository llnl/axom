// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_UTILITIES_HPP_
#define AXOM_BUMP_UTILITIES_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

#include <cstdint>

namespace axom
{
namespace bump
{
namespace utilities
{
//------------------------------------------------------------------------------
/*!
 * \brief This class and its specializations provide a type trait that lets us
 *        determine the type that should be used to accumulate values when we
 *        do floating point math.
 *
 * \note this belongs in algorithm utilities, maybe core.
 */
template <typename T>
struct accumulation_traits
{
  using type = float;
};

template <>
struct accumulation_traits<double>
{
  using type = double;
};

template <>
struct accumulation_traits<long>
{
  using type = double;
};

template <>
struct accumulation_traits<unsigned long>
{
  using type = double;
};

//------------------------------------------------------------------------------
/*!
 * \brief Fill an ArrayView with a value.
 *
 * \tparam ExecSpace The execution space where the fill will be done.
 * \tparam T The data type of the values in the ArrayView.
 *
 * \param view The ArrayView being filled.
 * \param fillValue The value to be used for filling the ArrayView.
 */
template <typename ExecSpace, typename T>
void fill(axom::ArrayView<T> view, T fillValue)
{
  axom::for_all<ExecSpace>(
    view.size(),
    AXOM_LAMBDA(axom::IndexType index) { view[index] = fillValue; });
}

//------------------------------------------------------------------------------
/*!
 * \brief Use binary search to find the index of the \a value in the supplied
 *        sorted view.
 *
 * \tparam T The type of elements we're handling.
 * \tparam ContainerT A container or view of T.
 *
 * \param[in] value The search value.
 * \param[in] view A view or container that contains the sorted search data values.
 *
 * \return The index where value was located in view or -1 if not found.
 */
template <typename T, typename ContainerT>
AXOM_HOST_DEVICE std::int32_t bsearch(T value, const ContainerT &view)
{
  std::int32_t index = -1;
  std::int32_t left = 0;
  std::int32_t right = view.size() - 1;
  while(left <= right)
  {
    std::int32_t m = (left + right) / 2;
    if(view[m] < value)
      left = m + 1;
    else if(view[m] > value)
      right = m - 1;
    else
    {
      index = m;
      break;
    }
  }

  return index;
}

//------------------------------------------------------------------------------
/*!
 * \brief Hash a stream of bytes into a uint64_t hash value.
 *
 * \param[in] data The bytes to be hashed.
 * \param[in] length The number of bytes to be hashed.
 *
 * \return A uint64_t hash for the byte stream.
 *
 * \note The byte stream is hashed using a Jenkins-hash algorithm forwards and
 *       backwards and the two results are merged into a uint64_t. The length is
 *       also part of the hash to guard against a lot of repeated values in the
 *       byte stream hashing to the same thing.
 *
 * \note We make this function inline since it is not a template and we want to
 *       use it in both host and device code.
 */
AXOM_HOST_DEVICE
inline std::uint64_t hash_bytes(const std::uint8_t *data, std::uint32_t length)
{
  std::uint32_t hash = 0;

  // Build the length into the hash.
  const auto ldata = reinterpret_cast<const std::uint8_t *>(&length);
  for(int e = 0; e < 4; e++)
  {
    hash += ldata[e];
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  std::uint32_t hashr = hash;
  for(std::uint32_t i = 0; i < length; i++)
  {
    hash += data[i];
    hash += hash << 10;
    hash ^= hash >> 6;

    hashr += data[length - 1 - i];
    hashr += hashr << 10;
    hashr ^= hashr >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;

  hashr += hashr << 3;
  hashr ^= hashr >> 11;
  hashr += hashr << 15;

  return (static_cast<std::uint64_t>(hash) << 32) | hashr;
}

//------------------------------------------------------------------------------
/*!
 * \brief Base template for computing a shape's area or volume.
 */
template <int NDIMS>
struct ComputeShapeAmount
{ };

/*!
 * \brief 2D specialization for shapes to compute area.
 */
template <>
struct ComputeShapeAmount<2>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.area();
  }
};

/*!
 * \brief 3D specialization for shapes to compute volume.
 */
template <>
struct ComputeShapeAmount<3>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.volume();
  }
};

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom

#endif
