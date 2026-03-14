// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_CONDUIT_MEMORY_HPP_
#define AXOM_BUMP_CONDUIT_MEMORY_HPP_

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
inline axom::ArrayView<T> make_array_view(conduit::Node &n)
{
  SLIC_ASSERT_MSG(cpp2conduit<T>::id == n.dtype().id(),
                  axom::fmt::format("Cannot create ArrayView<{}> for Conduit {} data.",
                                    cpp2conduit<T>::name,
                                    n.dtype().name()));
  return axom::ArrayView<T>(static_cast<T *>(n.data_ptr()), n.dtype().number_of_elements());
}

template <typename T>
inline axom::ArrayView<T> make_array_view(const conduit::Node &n)
{
  SLIC_ASSERT_MSG(cpp2conduit<T>::id == n.dtype().id(),
                  axom::fmt::format("Cannot create ArrayView<{}> for Conduit {} data.",
                                    cpp2conduit<T>::name,
                                    n.dtype().name()));
  return axom::ArrayView<T>(static_cast<T *>(const_cast<void *>(n.data_ptr())),
                            n.dtype().number_of_elements());
}
/// @}

//------------------------------------------------------------------------------
namespace internal
{
/*!
 * \brief Copies a Conduit tree in the \a src node to a new Conduit \a dest node,
 *        making sure to allocate array data in the appropriate memory space for
 *        the execution space.
 *
 * \tparam ExecSpace The destination execution space (e.g. axom::SEQ_EXEC).
 *
 * \param dest The conduit node that will receive the copied data.
 * \param src The source data to be copied.
 * \param destAllocatorID The allocator for the destination. It defaults to the allocator for ExecSpace.
 */
template <typename ExecSpace>
void copyImpl(conduit::Node &dest,
              const conduit::Node &src,
              int destAllocatorID,
              bool destAllocatorForDevice)
{
  dest.reset();
  if(src.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < src.number_of_children(); i++)
    {
      copyImpl<ExecSpace>(dest[src[i].name()], src[i], destAllocatorID, destAllocatorForDevice);
    }
  }
  else
  {
    const int srcAllocatorID = axom::getAllocatorIDFromPointer(src.data_ptr());
    const bool srcDataOnDevice =
      (srcAllocatorID == INVALID_ALLOCATOR_ID) ? false : isDeviceAllocator(srcAllocatorID);
    const bool deviceInvolved = srcDataOnDevice || destAllocatorForDevice;

    if(deviceInvolved || (!src.dtype().is_string() && src.dtype().number_of_elements() > 1))
    {
      // Allocate the node's memory in the right place.
      dest.reset();
      dest.set_allocator(axom::sidre::ConduitMemory::axomAllocIdToConduit(destAllocatorID));
      dest.set(conduit::DataType(src.dtype().id(), src.dtype().number_of_elements()));

      // Copy the data to the destination node. Axom uses Umpire to manage that.
      if(src.is_compact())
        axom::copy(dest.data_ptr(), src.data_ptr(), src.dtype().bytes_compact());
      else
      {
        // NOTE: This assumes that src is on the host.
        conduit::Node tmp;
        src.compact_to(tmp);
        axom::copy(dest.data_ptr(), tmp.data_ptr(), tmp.dtype().bytes_compact());
      }
      axom::synchronize<ExecSpace>();
    }
    else
    {
      // The data fits in the node or is a string. It's on the host.
      dest.set(src);
    }
  }
}

}  // end namespace internal

/*!
 * \brief Copies a Conduit tree in the \a src node to a new Conduit \a dest node,
 *        making sure to allocate array data in the appropriate memory space for
 *        the execution space.
 *
 * \tparam ExecSpace The destination execution space (e.g. axom::SEQ_EXEC).
 *
 * \param dest The conduit node that will receive the copied data.
 * \param src The source data to be copied.
 * \param destAllocatorID The allocator for the destination. It defaults to the allocator for ExecSpace.
 */
template <typename ExecSpace>
void copy(conduit::Node &dest,
          const conduit::Node &src,
          int destAllocatorID = axom::execution_space<ExecSpace>::allocatorID())
{
  const bool destAllocatorForDevice = isDeviceAllocator(destAllocatorID);
  internal::copyImpl<ExecSpace>(dest, src, destAllocatorID, destAllocatorForDevice);
}

//------------------------------------------------------------------------------
/*!
 * \brief Fill an array with int values from a Conduit node.
 *
 * \tparam ArrayType The array type being filled. It must supply size(), operator[].
 *
 * \param n The node that contains the data.
 * \param key The name of the node that contains the data in \a n.
 * \param[out] arr The array being filled.
 * \param moveToHost Sometimes data are on device and need to be moved to host first.
 */
template <typename ArrayType>
bool fillFromNode(const conduit::Node &n, const std::string &key, ArrayType &arr, bool moveToHost = false)
{
  bool found = false;
  if((found = n.has_path(key)) == true)
  {
    if(moveToHost)
    {
      // Make sure data are on host.
      conduit::Node hostNode;
      copy<axom::SEQ_EXEC>(hostNode, n.fetch_existing(key));

      const auto acc = hostNode.as_int_accessor();
      for(int i = 0; i < arr.size(); i++)
      {
        arr[i] = acc[i];
      }
    }
    else
    {
      const auto acc = n.fetch_existing(key).as_int_accessor();
      for(int i = 0; i < arr.size(); i++)
      {
        arr[i] = acc[i];
      }
    }
  }
  return found;
}

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom

#endif
