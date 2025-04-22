// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file NullSet.hpp
 *
 * \brief Sentinel set type indicating an empty set.
 *
 */

#ifndef SLAM_NULL_SET_H_
#define SLAM_NULL_SET_H_

#include "axom/slic.hpp"
#include "axom/slam/Set.hpp"

namespace axom
{
namespace slam
{
/**
 * \class NullSet
 *
 * \brief An indexed set (a tuple) of entities in a simulation
 */
template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
class NullSet : public Set<PosType, ElemType>
{
public:
  using ParentSet = Set<PosType, ElemType>;
  using PositionType = typename ParentSet::PositionType;
  using ElementType = typename ParentSet::ElementType;

public:
  NullSet() { }

  AXOM_HOST_DEVICE inline PositionType size() const { return PositionType(); }

  inline ElementType at(PositionType pos) const
  {
    verifyPosition(pos);
    return PositionType();
  }

  inline ElementType operator[](PositionType pos) const { return at(pos); }

  inline bool isSubset() const { return false; }
  const ParentSet* parentSet() const { return this; }

  bool isValid(bool AXOM_UNUSED_PARAM(verboseOutput) = false) const { return true; }

  AXOM_HOST_DEVICE bool empty() const { return true; }

  // TODO: Do we need to add iterator stubs here to satisfy some interface?
  //       The result will be invalid, but it may be useful to get the code
  //       to compile, or avoid special logic in the code...
  // iterator begin();
  // iterator end();
  // iterator_pair range();

private:
  void verifyPosition(PositionType AXOM_DEBUG_PARAM(pos)) const
  {
    SLIC_ASSERT_MSG(false,
                    "Subscripting on NullSet is never valid."
                      << "\n\tAttempted to access item at index " << pos << ".");
  }
};

#if 0
/**
 * \brief NullSets are always equal
 * \note Two sets of different types are (currently) considered to be unequal
 */
inline bool operator==(NullSet const&, NullSet const&)
{
  SLIC_WARNING("operator==(NullSet,NullSet)"); return true;
}
/**
 * \brief NullSets are always equal
 * \note Two sets of different types are (currently) considered to be unequal
 */
inline bool operator!=(NullSet const&, NullSet const&)
{
  SLIC_WARNING("operator!=(NullSet,NullSet)"); return false;
}
#endif

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_NULL_SET_H_
