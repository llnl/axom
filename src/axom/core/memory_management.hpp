// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MEMORYMANAGEMENT_HPP_
#define AXOM_MEMORYMANAGEMENT_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"

// Umpire includes
#ifdef AXOM_USE_UMPIRE
  #include "umpire/config.hpp"
  #include "umpire/ResourceManager.hpp"
  #include "umpire/op/MemoryOperationRegistry.hpp"
  #include "umpire/resource/MemoryResourceTypes.hpp"
  #include "umpire/strategy/QuickPool.hpp"
#else
  #include <cstring>  // for std::memcpy
  #include <cstdlib>  // for std::malloc, std::realloc, std::free
#endif

#include <iostream>
#include <type_traits>

namespace axom
{
// To co-exist with Umpire allocator ids, use negative values here.
constexpr int INVALID_ALLOCATOR_ID = -1;  //!< Place holder for no/unknown allocator
constexpr int MALLOC_ALLOCATOR_ID = -3;   //!< Refers to MemorySpace::Malloc

// _memory_space_start
/*!
 * \brief Memory spaces supported by Array-like types
 *
 * This abstraction is not implemented using Umpire's MemoryResourceType enum
 * in order to also include
 * - a "Malloc" option that uses malloc and free.
 * - a "Dynamic" option as a default template parameter
 *   for use in Array-like types (see axom::Array).  If using
 *   Umpire, "Dynamic" refers to the default Umpire allocator.
 *   If not using Umpire, "Dynamic" falls back on "Malloc".
 *   (See axom::setDefaultAllocator() and axom::getDefaultAllocator().)
 */
enum class MemorySpace
{
  Malloc,   //!< Host memory using malloc, free and realloc
  Dynamic,  //!< Refers to Umpire's current default allocator
#ifdef AXOM_USE_UMPIRE
  Host,     //!< Umpire's host memory space
  Device,   //!< Umpire's device memory space
  Unified,  //!< Umpire's unified memory space
  Pinned,   //!< Umpire's pinned memory space
  Constant  //!< Umpire's constant memory space
#endif
};
// _memory_space_end

// _memory_management_routines_start
/// \name Memory Management Routines
/// @{

#ifdef AXOM_USE_UMPIRE

/*!
 * \brief Returns the ID of the predefined allocator for a given resource.
 * \param [in] resource_type the Umpire resource type
 * \return ID the id of the predefined umpire allocator.
 */
inline int getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType resource_type)
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator alloc = rm.getAllocator(resource_type);
  return alloc.getId();
}

/*!
 * \brief Sets the default memory allocator to use.
 * \param [in] resource_type the Umpire resource type
 */
inline void setDefaultAllocator(umpire::resource::MemoryResourceType resource_type)
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(resource_type);
  rm.setDefaultAllocator(allocator);
}

#endif

/*!
 * \brief Sets the default memory allocator to use.
 * \param [in] allocId the Umpire allocator id
 * 
 * \note This function has no effect when Axom is not compiled with Umpire.
 */
inline void setDefaultAllocator(int allocId)
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(allocId);
  rm.setDefaultAllocator(allocator);
#else
  AXOM_UNUSED_VAR(allocId);
#endif
}

/*!
 * \brief Returns the ID of the current default Umpire allocator
 * or MALLOC_ALLOCATOR_ID if Umpire is not used.
 *
 * \return ID of the current Umpire default allocator or MALLOC_ALLOCATOR_ID.
 */
inline int getDefaultAllocatorID()
{
#ifdef AXOM_USE_UMPIRE
  return umpire::ResourceManager::getInstance().getDefaultAllocator().getId();
#else
  return MALLOC_ALLOCATOR_ID;
#endif
}

/*!
 * \brief Get the allocator id from which data has been allocated.
 * \return Allocator id.  If Umpire doesn't have an allocator for the
 * pointer, or if Axom wasn't configured with Umpire, assume the
 * non-null pointers are from a malloc and return \c
 * axom::MALLOC_ALLOCATOR_ID.
 *
 * \pre ptr has a valid pointer value.
 */
inline int getAllocatorIDFromPointer(const void* ptr)
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  if(rm.hasAllocator(const_cast<void*>(ptr)))
  {
    umpire::Allocator allocator = rm.getAllocator(const_cast<void*>(ptr));
    return allocator.getId();
  }
#endif
  return ptr == nullptr ? INVALID_ALLOCATOR_ID : MALLOC_ALLOCATOR_ID;
}

/*!
 * \brief Determines whether an allocator id is for shared memory.
 *
 * \param allocID An allocator id.
 *
 * \return True if the allocator id is for shared memory; false otherwise.
 */
bool isSharedMemoryAllocator(int allocID);

/*!
 * \brief Get the allocator ID for Axom's shared memory allocator.
 *
 * \return The allocator ID for Axom's shared memory allocator (if Axom is using Umpire),
 *         or INVALID_ALLOCATOR_ID otherwise.
 */
int getSharedMemoryAllocatorID();

/*!
 * \brief Allocates a chunk of memory of type T.
 *
 * \param [in] n the number of elements to allocate.
 * \param [in] allocID the Umpire allocator to use (optional)
 *
 * \tparam T the type of pointer returned.
 *
 * \note By default allocate() will use the current default allocator. The
 *  caller may explicitly specify a different allocator to use by supplying the
 *  second, optional argument, or change the default allocator by calling
 *  axom::setDefaultAllocator().
 *
 * \return p pointer to the new allocation or a nullptr if allocation failed.
 */
template <typename T>
inline T* allocate(std::size_t n, int allocID = getDefaultAllocatorID()) noexcept;

/*!
 * \brief Frees the chunk of memory pointed to by the supplied pointer, p.
 * \param [in/out] p a pointer to memory allocated with allocate/reallocate or a
 * nullptr.
 * \post p == nullptr
 */
template <typename T>
inline void deallocate(T*& p) noexcept;

/*!
 * \brief Reallocates the chunk of memory pointed to by the supplied pointer.
 *
 * \param [in] p pointer to memory allocated with allocate/reallocate, or a
 * nullptr.
 * \param [in] n the number of elements to allocate.
 * \param [in] allocID the ID of the allocator to use if pointer is null
 * (optional)
 *
 * \tparam T the type pointer p points to.
 *
 * \return p pointer to the new allocation or a nullptr if allocation failed.
 *
 * \note When n == 0, this function returns a valid pointer (of size 0) in the
 * current allocator's memory space. This follows the semantics of
 * Umpire's reallocate function.
 * \note When p is a null pointer, allocID is used to allocate the data.
 * Otherwise, it is unused.
 */
template <typename T>
inline T* reallocate(T* p, std::size_t n, int allocID = getDefaultAllocatorID()) noexcept;

/*!
 * \brief Copies memory from the source to the destination.
 *
 * \param [in/out] dst the destination to copy to.
 * \param [in] src the source to copy from.
 * \param [in] numbytes the number of bytes to copy.
 *
 * \note When using Umpire if either src or dst is not registered with the
 *  ResourceManager then the default host allocation strategy is assumed for
 *  that pointer.
 */
inline void copy(void* dst, const void* src, std::size_t numbytes) noexcept;

/*!
 * \brief Fills memory with a value.
 *
 * \param [in/out] dst the destination to copy to.
 * \param [in] n the number of items to copy.
 * \param [in] The value to copy. It must be trivially copyable for use with GPU.
 *
 * \note When using Umpire if dst is not registered with the
 *  ResourceManager then the default host allocation strategy is assumed for
 *  that pointer.
 */
template <typename T>
inline void fill(void* dst, std::size_t n, const T& value) noexcept;

/// @}
// _memory_management_routines_end

/*!
 * \brief Wrapper type representing an Umpire allocator ID.
 *
 *  This type is intended for use in function and constructor arguments, in
 *  order to avoid ambiguities in overload resolution.
 */
struct Allocator
{
public:
  explicit Allocator(int alloc_id = axom::getDefaultAllocatorID()) : m_id {alloc_id} { }

  /// \brief Returns the allocator ID.
  int getID() const { return m_id; }

private:
  int m_id;
};

//------------------------------------------------------------------------------
//                        IMPLEMENTATION
//------------------------------------------------------------------------------

template <typename T>
inline T* allocate(std::size_t n, int allocID) noexcept
{
  const std::size_t numbytes = n * sizeof(T);

#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  if(rm.isAllocator(allocID))
  {
    umpire::Allocator allocator = rm.getAllocator(allocID);
    return static_cast<T*>(allocator.allocate(numbytes));
  }
#endif

  if(allocID == MALLOC_ALLOCATOR_ID)
  {
    return static_cast<T*>(std::malloc(numbytes));
  }

  std::cerr << "Unrecognized allocator id " << allocID << std::endl;
  axom::utilities::processAbort();

  return nullptr;  // Silence warning.
}
//------------------------------------------------------------------------------
template <typename T>
inline void deallocate(T*& pointer) noexcept
{
  if(pointer == nullptr)
  {
    return;
  }

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  if(rm.hasAllocator(pointer))
  {
    rm.deallocate(pointer);
    pointer = nullptr;
    return;
  }

#endif

  std::free(pointer);
  pointer = nullptr;
}

//------------------------------------------------------------------------------
template <typename T>
inline T* reallocate(T* pointer, std::size_t n, int allocID) noexcept
{
  assert(allocID != INVALID_ALLOCATOR_ID);

  const std::size_t numbytes = n * sizeof(T);

#if defined(AXOM_USE_UMPIRE)

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  if(rm.isAllocator(allocID))
  {
    if(pointer == nullptr)
    {
      pointer = axom::allocate<T>(n, allocID);
    }
    else
    {
      if(rm.hasAllocator(pointer))
      {
        pointer = static_cast<T*>(rm.reallocate(pointer, numbytes));
      }
      else
      {
        /*
         * Reallocate from non-Umpire to Umpire, manually, using
         * allocate, copy and deallocate.  Because we don't know the
         * current size, we first do a (extra) reallocate within the
         * current space just so we have the size for the copy.
         * Is there a better way?
         */
        auto tmpPointer = std::realloc(pointer, numbytes);
        pointer = axom::allocate<T>(n, allocID);
        copy(pointer, tmpPointer, numbytes);
        deallocate(tmpPointer);
      }
    }
    return pointer;
  }

#else

  if(allocID == MALLOC_ALLOCATOR_ID)
  {
    pointer = static_cast<T*>(std::realloc(pointer, numbytes));
  }
  else
  {
    std::cerr << "*** Unrecognized allocator id "
              << allocID << ".  Axom was NOT built with Umpire, so the only valid allocator id is MALLOC_ALLOCATOR_ID ("
              << MALLOC_ALLOCATOR_ID << ")." << std::endl;
    axom::utilities::processAbort();
  }

  // Consistently handle realloc(0) for std::realloc to match Umpire's behavior
  if(n == 0 && pointer == nullptr)
  {
    pointer = axom::allocate<T>(0);
  }

#endif

  return pointer;
}

//------------------------------------------------------------------------------
inline void copy(void* dst, const void* src, std::size_t numbytes) noexcept
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::op::MemoryOperationRegistry& op_registry =
    umpire::op::MemoryOperationRegistry::getInstance();

  auto dstStrategy = rm.getAllocator("HOST").getAllocationStrategy();
  auto srcStrategy = dstStrategy;

  using AllocationRecord = umpire::util::AllocationRecord;
  AllocationRecord* dstRecord = nullptr;
  AllocationRecord* srcRecord = nullptr;

  if(rm.hasAllocator(dst))
  {
    dstRecord = const_cast<AllocationRecord*>(rm.findAllocationRecord(dst));
    dstStrategy = dstRecord->strategy;
  }

  if(rm.hasAllocator(const_cast<void*>(src)))
  {
    srcRecord = const_cast<AllocationRecord*>(rm.findAllocationRecord(const_cast<void*>(src)));
    srcStrategy = srcRecord->strategy;
  }

  auto op = op_registry.find("COPY", srcStrategy, dstStrategy);
  op->transform(const_cast<void*>(src), &dst, srcRecord, dstRecord, numbytes);
#else
  std::memcpy(dst, src, numbytes);
#endif
}

//------------------------------------------------------------------------------
template <typename T>
inline void fill(void* dst, std::size_t n, const T& value) noexcept
{
  bool doHostFill = true;
#ifdef AXOM_USE_UMPIRE
  // Since data might be copied to GPU, it needs to be trivially copyable.
  static_assert(std::is_trivially_copyable<T>::value, "value must be trivially copyable.");
  auto& rm = umpire::ResourceManager::getInstance();

  if(rm.hasAllocator(dst))
  {
    auto alloc = rm.getAllocator(dst);
    if((alloc.getPlatform() != umpire::Platform::host))
    {
      doHostFill = false;

      // Device memory: fill on host, then copy to device
      const auto num_bytes = n * sizeof(T);
      T* src = allocate<T>(num_bytes, rm.getDefaultAllocator().getId());
      for(std::size_t i = 0; i < n; ++i)
      {
        src[i] = value;
      }
      rm.copy(dst, src, num_bytes);
      deallocate<T>(src);
    }
  }
#endif
  if(doHostFill)
  {
    T* typed_dst = static_cast<T*>(dst);
    for(std::size_t i = 0; i < n; ++i)
    {
      typed_dst[i] = value;
    }
  }
}

namespace detail
{
/// \brief Translates between the MemorySpace enum and Umpire allocator IDs
template <MemorySpace SPACE>
inline int getAllocatorID();

template <>
inline int getAllocatorID<MemorySpace::Dynamic>()
{
  /*
    With Umpire enabled, this returns the current default Umpire id.
    Without Umpire, it returns MALLOC_ALLOCATOR_ID.
  */
  return axom::getDefaultAllocatorID();
}

template <>
inline int getAllocatorID<MemorySpace::Malloc>()
{
  return axom::MALLOC_ALLOCATOR_ID;
}

/**
 * @brief Return the Axom MemorySpace for the given Axom allocator id.
 *
 * For Umpire allocator ids, the MemorySpace is the corresponding Axom
 * memory space.  For MALLOC_ALLOCATOR_ID, the MemorySpace is
 * MemorySpace::Malloc.  Other values have no corresponding MemorySpace
 * and will cause an abort.
 */
inline MemorySpace getAllocatorSpace(int allocatorId)
{
#ifdef AXOM_USE_UMPIRE
  using ump_res_type = typename umpire::MemoryResourceTraits::resource_type;

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  if(rm.isAllocator(allocatorId))
  {
    auto umpResType = rm.getAllocator(allocatorId).getAllocationStrategy()->getTraits().resource;
    switch(umpResType)
    {
    case ump_res_type::host:
      return MemorySpace::Host;
    case ump_res_type::device:
      return MemorySpace::Device;
    case ump_res_type::device_const:
      return MemorySpace::Constant;
    case ump_res_type::pinned:
      return MemorySpace::Pinned;
    case ump_res_type::um:
      return MemorySpace::Unified;
    default:
      return MemorySpace::Dynamic;
    }
  }
#endif
  if(allocatorId == MALLOC_ALLOCATOR_ID)
  {
    return MemorySpace::Malloc;
  }

  std::cerr << "*** Unrecognized allocator id " << allocatorId << "." << std::endl;
  axom::utilities::processAbort();

  return MemorySpace::Malloc;  // Silence warning.
}

#ifdef AXOM_USE_UMPIRE

template <>
inline int getAllocatorID<MemorySpace::Host>()
{
  return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
}

template <>
inline int getAllocatorID<MemorySpace::Device>()
{
  return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
}

template <>
inline int getAllocatorID<MemorySpace::Unified>()
{
  return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
}

template <>
inline int getAllocatorID<MemorySpace::Pinned>()
{
  return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Pinned);
}

template <>
inline int getAllocatorID<MemorySpace::Constant>()
{
  return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Constant);
}

#endif

}  // namespace detail

/*!
 * \brief Determines whether an allocator id is on device.
 *
 * \param allocator_id An allocator id.
 *
 * \return True if the allocator id is for a device; false otherwise.
 */
#if defined(AXOM_USE_UMPIRE)
inline bool isDeviceAllocator(int allocator_id)
{
  return axom::detail::getAllocatorSpace(allocator_id) == axom::MemorySpace::Device;
}
#else
inline bool isDeviceAllocator(int AXOM_UNUSED_PARAM(allocator_id)) { return false; }
#endif

}  // namespace axom

#endif /* AXOM_MEMORYMANAGEMENT_HPP_ */
