// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ConduitMemory.hpp
 *
 * \brief   Call-backs for using Axom memory management in Conduit.
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
namespace sidre
{

/*!
 * @brief Object to do Conduit memory operations through Axom.
 *
 * This class has no public constructor.  Use instanceForAxomId(int
 * axomAllocId) to access the instance for a specific Axom allocator
 * id.  The construction registers the appropriate callbacks with
 * Conduit, including the required memset and memcopy callbacks.
 *
 * Allocator ids have a 1-to-1 relationship with allocators.
 *
 * Axom's allocator is an extension of the Umpire allocator
 * when Umpire is used.  Conduit's allocator is opaque, but when used
 * by this class, it is associated with an Axom allocator (which is
 * an Umpire allocator).
 *
 * Examples for setting Conduit allocator ids when you have Axom
 * allocator ids:
 *
 * @code{.cpp}
 *   void foo(conduit::Node& n, int axomAllocId) {
 *     n.set_allocator(axomAllocIdToConduit(axomAllocId));
 *   }
 *
 *   void bar(conduit::Node& n, int axomAllocId) {
 *     const auto& instance = getInstance(axomAllocId);
 *     assert(instance.axomId() == axomAllocId);
 *     n.set_allocator(instance.conduitId());
 *   }
 * @endcode
 */
struct ConduitMemory
{
  //!@brief Return the Axom allocator id.
  int axomId() const { return m_axomId; }

  //!@brief Return the Conduit allocator id coresponding to axomId().
  conduit::index_t conduitId() const { return m_conduitId; }

  /*!
   * @brief Convert an Axom allocator id to Conduit, registering
   * a new Conduit allocator if needed.
   */
  static conduit::index_t axomAllocIdToConduit(int axomAllocId)
  {
    return instanceForAxomId(axomAllocId).conduitId();
  }

  /*!
   * @brief Convert a Conduit allocator id to Axom.
   *
   * The allocator must have been registered by a prior
   * instanceForAxomId() call.
   */
  static int conduitAllocIdToAxom(conduit::index_t conduitAllocId)
  {
    return instanceForConduitId(conduitAllocId).axomId();
  }

  /*!
   * @brief Return the instance for the given Axom allocator id.
   *
   * This method IS NOT thread safe for new values of @c axomAllocId.
   */
  static const ConduitMemory& instanceForAxomId(int axomAllocId);

  /*!
   * @brief Return the instance for the given Conduit allocator id.
   *
   * If @c conduitAllocId doesn't correspond to an Axom allocator,
   * an object corresponding to axom::INVALID_ALLOCATOR_ID will be returned.
   *
   * This method IS thread safe.
   */
  static const ConduitMemory& instanceForConduitId(conduit::index_t conduitAllocId);

  //!@brief Return the default conduit allocator id.
  static conduit::index_t defaultConduitId() { return s_defaultConduitId; }

  ~ConduitMemory() { }

private:
  //!@brief Mapping from Axom allocator to an instance.
  static std::map<int, std::shared_ptr<ConduitMemory>> s_axomToInstance;

  //!@brief Mapping from Conduit allocator to an instance.
  static std::map<conduit::index_t, std::shared_ptr<ConduitMemory>> s_conduitToInstance;

  //!@brief Conduit's default allocator id.
  static const conduit::index_t s_defaultConduitId;

  //!@brief Axom's allocator id.
  int m_axomId;
  //!@brief Conduit's allocator id equivalent to m_axomId.
  conduit::index_t m_conduitId;

#if (CONDUIT_VERSION_MAJOR >= 0 && CONDUIT_VERSION_MINOR >= 9 && CONDUIT_VERSION_PATCH > 4) || \
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
   * @brief Constructor creates allocator/deallocator function and registers
   * them with Conduit.
   */
  explicit ConduitMemory(int axomAllocId) : m_axomId(axomAllocId) { privateRegisterAllocator(); }

  void privateRegisterAllocator();
};

} /* end namespace sidre */
} /* end namespace axom */

#endif  // AXOM_USE_CONDUIT
