// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file HandleMesh.cpp
 *
 * \brief Slam example where the ElementType of a set is a handle
 * that wraps a mesh index.
 */

#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <ostream>

namespace slam = axom::slam;

namespace
{
/**
 * \class Handle
 * \brief A handle is a wrapper around a variable.
 *
 * This class is a stub for some mesh data structures
 * that use typed index handles rather than (untyped) indices
 * to refer to mesh elements.
 */
template <typename T>
struct Handle
{
  Handle() : mID(T()) { }
  explicit Handle(T id) : mID(id) { }
  Handle(const Handle& h) : mID(h.mID) { }
  Handle& operator=(const Handle& h) = default;

  bool operator==(const Handle& h) const { return mID == h.mID; }

  static Handle make_handle(T id) { return Handle(id); }

  void print(std::ostream& os) const { os << "Handle(" << mID << ")"; }

  // Simple function to call on a Handle
  T twiceIndex() const { return 2 * mID; }

  T mID;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Handle<T>& h)
{
  h.print(os);
  return os;
}

}  // end anonymous namespace

int main(int, char**)
{
  axom::slic::SimpleLogger logger;

  using PosType = slam::DefaultPositionType;
  using HandleType = Handle<PosType>;

  using HandleSet = slam::VectorIndirectionSet<PosType, HandleType>;

  const int sz = 5;
  HandleSet::IndirectionBufferType vecHandle(sz);

  // Create a set of handles
  HandleSet hSet = HandleSet::SetBuilder()  //
                     .size(sz)              //
                     .data(&vecHandle);

  // Add handles with (somewhat) arbitrary IDs to the set
  for(auto i : hSet.positions())
  {
    hSet[i] = HandleType::make_handle(2 * sz - i);
  }

  // Iterate over the set
  SLIC_INFO("Iterating a set of Handle elements: ");
  for(auto i : hSet.positions())
  {
    auto it = hSet.begin() + i;

    SLIC_INFO("  " << i << ": " << hSet[i] << " -- double of index is: " << it->twiceIndex());
  }

  // Check equality
  SLIC_INFO("Checking equality of Handle elements: ");
  SLIC_INFO("  hSet[0] == hSet[0] ? " << (hSet[0] == hSet[0] ? "yes" : "no"));
  SLIC_INFO("  hSet[0] == hSet[1] ? " << (hSet[0] == hSet[1] ? "yes" : "no"));

  return 0;
}
