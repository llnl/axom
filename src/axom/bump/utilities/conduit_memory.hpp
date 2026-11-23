// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/bump/utilities/conduit_array_view.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/execution/synchronize.hpp"
#include "axom/slic.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"
#include "axom/export/bump.h"

#include <conduit/conduit.hpp>

#include <string>
#include <type_traits>

namespace axom
{
namespace bump
{
namespace utilities
{
//------------------------------------------------------------------------------

/*!
 * \brief Make an axom::ArrayView from a Conduit node.
 *
 * \tparam T The type for the array view elements.
 *
 * \param n The conduit node for which we want an array view.
 *
 * \return An axom::ArrayView that wraps the data in the Conduit node.
 */
/// @{
template <typename T>
inline axom::ArrayView<T> make_array_view(conduit::Node& n)
{
  return detail::make_conduit_array_view<T>(n);
}

template <typename T>
inline axom::ArrayView<T> make_array_view(const conduit::Node& n)
{
  return detail::make_conduit_array_view<T>(n);
}
/// @}

//------------------------------------------------------------------------------
/*!
 * \brief Returns whether the data pointer in the Conduit node appears to be
 *        allocated on a device.
 *
 * \param n The node whose data we're checking for device allocation.
 *
 * \return True if the data looks device-allocated; false otherwise.
 */
bool isDeviceAllocated(const conduit::Node& n);

//------------------------------------------------------------------------------
/*!
 * \brief Copies a Conduit tree in the \a src node to a new Conduit \a dest node,
 *        compacting the source into a single destination allocation in the
 *        appropriate memory space for the execution space.
 *
 * \tparam ExecSpace The destination execution space (e.g. axom::SEQ_EXEC).
 *
 * \param dest The conduit node that will receive the copied data.
 * \param src The source data to be copied.
 * \param destAllocatorID The allocator for the destination. It defaults to the allocator for ExecSpace.
 *
 * \note The destination root node keeps the supplied allocator, and the compacted
 *       data buffer is allocated in that memory space. However, Conduit's
 *       `compact_to()` currently rebuilds nested destination child nodes using the
 *       source allocator id for their `Node::allocator()` metadata. Current `bump`
 *       code relies on the actual data pointers rather than nested `Node::allocator()`
 *       values, but this caveat may matter if future code starts allocating new
 *       children under copied nested subtrees.
 */
template <typename ExecSpace>
void copy(conduit::Node& dest,
          const conduit::Node& src,
          int destAllocatorID = axom::execution_space<ExecSpace>::allocatorID())
{
  // No explicit dest.reset() is needed here:
  // - set_allocator() resets if the allocator changes, and
  // - compact_to() resets the destination before rebuilding the compact tree.
  dest.set_allocator(axom::sidre::ConduitMemory::axomAllocIdToConduit(destAllocatorID));
  src.compact_to(dest);
  axom::synchronize<ExecSpace>();
}

//------------------------------------------------------------------------------
/*!
 * \brief Reads a Conduit value as a specific value type, using Node::to_* methods
 *        to allow some flexibility in the types.
 *
 * \tparam T The desired type
 *
 * \param n The Node to read
 *
 * \return The Node's value converted to type T.
 */
template <typename T>
inline T readValue(const conduit::Node& n)
{
  if constexpr(std::is_same_v<T, std::string>)
  {
    return n.as_string();
  }
  else if constexpr(std::is_same_v<T, conduit::index_t>)
  {
    return n.to_index_t();
  }
  else if constexpr(std::is_same_v<T, float>)
  {
    return n.to_float();
  }
  else if constexpr(std::is_same_v<T, double>)
  {
    return n.to_double();
  }
  else if constexpr(std::is_integral_v<T> && std::is_signed_v<T>)
  {
    return static_cast<T>(n.to_int64());
  }
  else if constexpr(std::is_integral_v<T> && std::is_unsigned_v<T>)
  {
    return static_cast<T>(n.to_uint64());
  }
  else
  {
    static_assert(!std::is_same_v<T, T>, "readValue does not support this type.");
  }
}

//------------------------------------------------------------------------------
/*!
 * \brief Read a Conduit leaf-node value on the host, copying just that leaf to
 *        host first when its data are device-resident.
 *
 * \tparam T The host value type to return.
 *
 * \param n The Conduit leaf node to read.
 *
 * \return The requested host value.
 *
 * \note This helper is only valid for leaf nodes. It raises `SLIC_ERROR` for
 *       empty, object, or list nodes.
 *
 * \note Non-string leaves must contain exactly one element.
 *
 * \note Numeric conversions use Conduit's `to_*()` helpers so callers can read
 *       compatible scalar types without requiring an exact stored type match.
 */
template <typename T>
T hostRead(const conduit::Node& n)
{
  const auto& dtype = n.dtype();
  const auto dtypeID = dtype.id();
  SLIC_ERROR_IF(dtypeID == conduit::DataType::EMPTY_ID ||
                  dtypeID == conduit::DataType::OBJECT_ID ||
                  dtypeID == conduit::DataType::LIST_ID,
                "hostRead only supports Conduit leaf nodes.");
  SLIC_ERROR_IF(!dtype.is_string() && dtype.number_of_elements() != 1,
                "hostRead only supports scalar leaves and string leaves.");

  // Use the actual data pointer to decide whether a host copy is needed.
  // This avoids relying on nested Node::allocator() metadata, which Conduit can
  // currently rebuild incorrectly after compact_to().
  if(isDeviceAllocated(n))
  {
    conduit::Node hostNode;
    copy<axom::SEQ_EXEC>(hostNode, n);
    return readValue<T>(hostNode);
  }

  return readValue<T>(n);
}

//------------------------------------------------------------------------------
/*!
 * \brief Read a Conduit leaf-node array value on the host, copying just that
 *        leaf to host first when its data are device-resident.
 *
 * \tparam T The host array element type to return.
 *
 * \param n The Conduit leaf node to read.
 *
 * \return An `axom::Array<T>` containing the node values.
 *
 * \note This helper is only valid for leaf nodes. It raises `SLIC_ERROR` for
 *       empty, object, or list nodes.
 */
template <typename T>
axom::Array<T> hostReadArray(const conduit::Node& n)
{
  const auto& dtype = n.dtype();
  const auto dtypeID = dtype.id();

  SLIC_ERROR_IF(dtypeID == conduit::DataType::EMPTY_ID ||
                dtypeID == conduit::DataType::OBJECT_ID ||
                dtypeID == conduit::DataType::LIST_ID,
                "hostReadArray only supports Conduit leaf nodes.");

  axom::Array<T> values(dtype.number_of_elements(), dtype.number_of_elements());

  // Decide whether a host copy is needed.
  if(isDeviceAllocated(n))
  {
    // Copy to the host.
    conduit::Node hostNode;
    copy<axom::SEQ_EXEC>(hostNode, n);

    const conduit::DataAccessor<T> accessor = hostNode.value();
    for(conduit::index_t i = 0; i < accessor.number_of_elements(); i++)
    {
      values[i] = accessor[i];
    }
  }
  else
  {
    const conduit::DataAccessor<T> accessor = n.value();
    for(conduit::index_t i = 0; i < accessor.number_of_elements(); i++)
    {
      values[i] = accessor[i];
    }
  }

  return values;
}

//------------------------------------------------------------------------------
/*!
 * \brief Read a Conduit leaf node as an `int` on the host.
 *
 * \param n The Conduit leaf node to read.
 *
 * \return The node value converted to `int`.
 */
inline int hostReadInt(const conduit::Node& n) { return hostRead<int>(n); }

/*!
 * \brief Read a Conduit leaf node as a `float` on the host.
 *
 * \param n The Conduit leaf node to read.
 *
 * \return The node value converted to `float`.
 */
inline double hostReadFloat(const conduit::Node& n) { return hostRead<float>(n); }

/*!
 * \brief Read a Conduit leaf node as a `double` on the host.
 *
 * \param n The Conduit leaf node to read.
 *
 * \return The node value converted to `double`.
 */
inline double hostReadDouble(const conduit::Node& n) { return hostRead<double>(n); }

/*!
 * \brief Read a Conduit leaf node as a `std::string` on the host.
 *
 * \param n The Conduit leaf node to read.
 *
 * \return The node value converted to `std::string`.
 */
inline std::string hostReadString(const conduit::Node& n)
{
  return hostRead<std::string>(n);
}

//------------------------------------------------------------------------------
/*!
 * \brief Fill an array-like object from a Conduit node.
 *
 * \tparam ArrayType The array type being filled. It must supply size(), operator[].
 *
 * \param n The node that contains the data.
 * \param key The name of the node that contains the data in \a n.
 * \param[out] arr The array being filled.
 *
 * \return true if \a key exists in \a n, in which case \a arr is filled.
 *
 * \note Fills the first arr.size() entries. Raises `SLIC_ERROR` if the node has
 *       fewer values than that. Extra values in the node are ignored.
 *
 * \note If the source leaf data are device-resident, they are first copied to a
 *       temporary host node before the destination array is filled.
 */
template <typename ArrayType>
bool fillFromNode(const conduit::Node& n, const std::string& key, ArrayType& arr)
{
  using T = typename ArrayType::value_type;

  bool found = false;
  if((found = n.has_path(key)) == true)
  {
    const auto fill = [&](const conduit::Node& valuesNode) {
      const auto& dtype = valuesNode.dtype();
      const auto dtypeID = dtype.id();

      SLIC_ERROR_IF(dtypeID == conduit::DataType::EMPTY_ID ||
                      dtypeID == conduit::DataType::OBJECT_ID ||
                      dtypeID == conduit::DataType::LIST_ID,
                    "Node '" << key << "' must be a Conduit leaf node.");

      const conduit::DataAccessor<T> accessor = valuesNode.value();
      const auto numValues = static_cast<axom::IndexType>(accessor.number_of_elements());
      const auto numRequired = static_cast<axom::IndexType>(arr.size());

      SLIC_ERROR_IF(numValues < numRequired,
                    "Node '" << key << "' has " << numValues << " values, but "
                             << numRequired << " are required.");

      for(axom::IndexType i = 0; i < numRequired; i++)
      {
        arr[i] = accessor[i];
      }
    };

    const conduit::Node& valuesNode = n.fetch_existing(key);
    if(isDeviceAllocated(valuesNode))
    {
      conduit::Node hostNode;
      copy<axom::SEQ_EXEC>(hostNode, valuesNode);
      fill(hostNode);
    }
    else
    {
      fill(valuesNode);
    }
  }
  return found;
}

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom
