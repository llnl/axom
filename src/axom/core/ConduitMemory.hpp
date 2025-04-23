// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ConduitMemory.hpp
 *
 * \brief   Call-backs for using Axom memory mamagement in Conduit.
 *
 ******************************************************************************
 */

#ifndef SIDRE_CONDUITMEMORY_HPP_
#define SIDRE_CONDUITMEMORY_HPP_

// Standard C++ headers
#include <string>
#include <set>
#include <memory>

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "conduit_node.hpp"
#include "conduit_utils.hpp"

namespace axom
{

/*!
  @brief Object to do Conduit memory operations through Axom.

  This class has no public constructor.  Use instanceForAxomId(int
  axomAllocId) to access the instance for a specific Axom allocator
  id.  The construction registers the appropriate callbacks with
  Conduit, including the required memset and memcopy callbacks.

  Allocator ids have a 1-to-1 relationship with allocators.

  Axom's allocator is an extension of the Umpire allocator
  when Umpire is used.  Conduit's allocator is opaque, but when used
  by this class, it is associated with an Axom allocator (which is
  an Umpire allocator).

  Examples for setting Conduit allocator ids when you have Axom
  allocator ids:

  @code{.cpp}
    void foo(conduit::Node& n, int axomAllocId) {
      n.set_allocator(axomAllocIdToConduit(axomAllocId));
    }

    void bar(conduit::Node& n, int axomAllocId) {
      const auto& instance = getInstance(axomAllocId);
      assert(instance.axomId() == axomAllocId);
      n.set_allocator(instance.conduitId());
    }
  @endcode
*/
struct ConduitMemory
{
  //!@brief Return the Axom allocator id.
  int axomId() const { return m_axomId; }

  //!@brief Return the Conduit allocator id coresponding to axomId().
  conduit::index_t conduitId() const { return m_conduitId; }

  /*!
    @brief Convert an Axom allocator id to Conduit, registering
    a new Conduit allocator if needed.
  */
  static conduit::index_t axomAllocIdToConduit(int axomAllocId)
  {
    return instanceForAxomId(axomAllocId).conduitId();
  }

  /*!
    @brief Convert a Conduit allocator id to Axom.

    The allocator must have been registered by a prior
    instanceForAxomId() call.
  */
  static int conduitAllocIdToAxom(conduit::index_t conduitAllocId)
  {
    return instanceForConduitId(conduitAllocId).axomId();
  }

  //!@brief Return the instance for the given Axom allocator id.
  static const ConduitMemory& instanceForAxomId(int axomAllocId);

  //!@brief Return the instance for the given Conduit allocator id.
  static const ConduitMemory& instanceForConduitId(conduit::index_t conduitAllocId);

  ~ConduitMemory() { }

private:
  //!@brief Mapping from Axom allocator to an instance.
  static std::map<int, std::shared_ptr<ConduitMemory>> s_axomToInstance;

  //!@brief Mapping from Conduit allocator to an instance.
  static std::map<conduit::index_t, std::shared_ptr<ConduitMemory>> s_conduitToInstance;

  //!@brief Axom's allocator id.
  int m_axomId;
  //!@brief Conduit's allocator id equivalent to m_axomId.
  conduit::index_t m_conduitId;

#if(CONDUIT_VERSION_MAJOR >= 0 && CONDUIT_VERSION_MINOR >= 9 && CONDUIT_VERSION_PATCH > 4) || \
  (CONDUIT_VERSION_MAJOR >= 0 && CONDUIT_VERSION_MINOR >= 10) || CONDUIT_VERSION_MAJOR >= 1
  #define AXOM_CONDUIT_USES_STD_FUNCTION 1
#endif

#if defined(AXOM_CONDUIT_USES_STD_FUNCTION)
  using AllocatorCallback = std::function<void*(size_t, size_t)>;
  using DeallocCallback = std::function<void(void*)>;
#else
  typedef void* (*AllocatorCallback)(size_t, size_t);
  typedef void (*DeallocCallback)(void*);
#endif

  AllocatorCallback m_allocCallback;
  DeallocCallback m_deallocCallback;

  static void staticDeallocator(void* ptr)
  {
    char* cPtr = (char*)(ptr);
    axom::deallocate<char>(cPtr);
  }

  ConduitMemory() = delete;

  /*!
    @brief Constructor creates allocator/deallocator function and registers
    them with Conduit.
  */
  ConduitMemory(int axomAllocId) : m_axomId(axomAllocId)
  {
    using conduit::utils::register_allocator;
    auto deallocator = [](void* ptr) -> void {
      char* cPtr = (char*)(ptr);
      axom::deallocate<char>(cPtr);
    };
    m_deallocCallback = deallocator;
#if defined(AXOM_USE_CONDUIT_STD_FUNCTION)
    m_allocCallback = [=](size_t itemCount, size_t itemByteSize) -> void* {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, axomAllocId);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
#else
    /*
      Note: conduit-0.9.4 allows the callbacks as std::function types.
      Once we are there, we can use a single allocator, eliminating
      the need for these if-else blocks.
    */
    if(axomAllocId == MALLOC_ALLOCATOR_ID)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, MALLOC_ALLOCATOR_ID);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 0)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 0);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 1)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 1);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 2)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 2);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 3)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 3);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 4)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 4);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 5)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 5);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else
    {
      std::cerr << "*** Work-around for conduit::utils::register_allocator "
                   "needs case for "
                   "axomAllocId = "
                << std::to_string(axomAllocId) << ".  Please add it to ConduitMemory.hpp.";
      axom::utilities::processAbort();
    }
#endif
  }
};

} /* end namespace axom */

#endif  // AXOM_USE_CONDUIT
