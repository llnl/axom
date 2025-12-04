// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_CONDUIT_MEMORY_HPP_
#define AXOM_BUMP_CONDUIT_MEMORY_HPP_

#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
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
/*!
 * \brief This class registers a Conduit allocator that can make Conduit allocate
 *        through Axom's allocate/deallocate functions using a specific allocator.
 *        This permits Conduit to allocate through Axom's UMPIRE logic.
 *
 * \tparam ExecSpace The execution space.
 */
template <typename ExecSpace>
class ConduitAllocateThroughAxom
{
public:
  /*!
   * \brief Get the Conduit allocator ID for this ExecSpace.
   *
   * \return The Conduit allocator ID for this ExecSpace.
   */
  static conduit::index_t getConduitAllocatorID()
  {
    constexpr conduit::index_t NoAllocator = -1;
    static conduit::index_t conduitAllocatorID = NoAllocator;
    if(conduitAllocatorID == NoAllocator)
    {
      conduitAllocatorID = conduit::utils::register_allocator(internal_allocate, internal_free);
    }
    return conduitAllocatorID;
  }

private:
  /*!
   * \brief A function we register with Conduit to allocate memory.
   *
   * \param items The number of items to allocate.
   * \param item_size The size of each item in bytes.
   *
   * \brief A block of newly allocated memory large enough for the requested items.
   */
  static void *internal_allocate(size_t items, size_t item_size)
  {
    int axomAllocatorID;
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_GPU)
    constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

    axomAllocatorID = on_device ? axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
                                : axom::execution_space<ExecSpace>::allocatorID();
#else
    axomAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
#endif

    void *ptr = static_cast<void *>(axom::allocate<std::uint8_t>(items * item_size, axomAllocatorID));
    //std::cout << axom::execution_space<ExecSpace>::name()
    //  << ": Allocated for Conduit via axom: items=" << items
    //  << ", item_size=" << item_size << ", ptr=" << ptr << std::endl;
    return ptr;
  }

  /*!
   * \brief A deallocation function we register with Conduit.
   */
  static void internal_free(void *ptr)
  {
    //std::cout << axom::execution_space<ExecSpace>::name()
    //  << ": Dellocating for Conduit via axom: ptr=" << ptr << std::endl;
    axom::deallocate(ptr);
  }
};

//------------------------------------------------------------------------------
/*!
 * \brief Copies a Conduit tree in the \a src node to a new Conduit \a dest node,
 *        making sure to allocate array data in the appropriate memory space for
 *        the execution space.
 *
 * \tparam ExecSpace The destination execution space (e.g. axom::SEQ_EXEC).
 *
 * \param dest The conduit node that will receive the copied data.
 * \param src The source data to be copied.
 */
template <typename ExecSpace>
void copy(conduit::Node &dest, const conduit::Node &src)
{
  ConduitAllocateThroughAxom<ExecSpace> c2a;
  dest.reset();
  if(src.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < src.number_of_children(); i++)
    {
      copy<ExecSpace>(dest[src[i].name()], src[i]);
    }
  }
  else
  {
    const int allocatorID = axom::getAllocatorIDFromPointer(src.data_ptr());
    bool deviceAllocated =
      (allocatorID == INVALID_ALLOCATOR_ID) ? false : isDeviceAllocator(allocatorID);
    if(deviceAllocated || (!src.dtype().is_string() && src.dtype().number_of_elements() > 1))
    {
      // Allocate the node's memory in the right place.
      dest.reset();
      dest.set_allocator(c2a.getConduitAllocatorID());
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
    }
    else
    {
      // The data fits in the node or is a string. It's on the host.
      dest.set(src);
    }
  }
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
