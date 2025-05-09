// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Set.hpp
 *
 * \brief Basic API for a Set of entities in a simulation
 *
 */

#ifndef SLAM_SET_H_
#define SLAM_SET_H_

#include <cstddef>
#include <vector>
#include <type_traits>  // for std::common_type

#include "axom/core/utilities/Utilities.hpp"
#include "axom/slam/Utilities.hpp"

namespace axom
{
namespace slam
{
/**
 * \class Set
 *
 * \brief Abstract base class for a Set of entities in a simulation
 *
 * This class defines the minimal required API for a slam Set,
 * a container class for a set of entities in a simulation.
 * Each entity has an index.
 *
 * Examples of sets include:
 * <ol>
 *  <li> Mesh elements: vertices, edges, faces, cells
 *  <li> Subzonal elements: sides, corners, finite element degrees of freedom
 *  <li> Boundary elements: external surfaces and springs
 *  <li> Elements of a space partition: e.g. Domains in a block structured mesh,
 *       leaf nodes of an octree/kd-tree
 *  <li> AMR bricks / tiles
 *  <li> Thread ids, MPI ranks, warps, thread groups, etc...
 *  <li> particles
 *  <li> boundary conditions
 *  <li> ...
 * </ol>
 *
 * Examples of subsets include:
 * <ol>
 *  <li> Regions
 *  <li> Ghost cells -- send, receive
 *  <li> Boundary cells -- external surface
 * </ol>
 *
 * Note: Elements of a set do not necessarily need explicit indices.
 * E.g. if we have a contiguous range of elements (or slices of contiguous
 * ranges),
 * they can be implicitly encoded.
 *
 * Thus, we can have
 * <ol>
 *  <li> Implicit indexes -- all we need here is a size operator
 *  <li> Sliced indices -- here we need the dimension and the striding
 *  <li> Explicit indices -- for a subset, we need the indices with respect to
 *       some other indexing scheme
 * </ol>
 *
 * The interface is for constant access to the elements.
 */
template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
class Set
{
public:
  using PositionType = PosType;
  using ElementType = ElemType;

public:
  // Set () {}
  virtual ~Set() = default;

  /**
   * \brief Random access to the entities of the set
   * \param The index of the desired element
   * \return The value of the element at the given position
   * \pre The position must be less than the number of elements in the set (
   * size() )
   * \note Concrete realizations of Set also support subscript operator --
   * operator[].
   * \note How are we planning to handle indexes that are out or range
   *(accidentally)?
   *       Are we planning to handle indexes that are intentionally out of range
   *       (e.g. to indicate a problem, or a missing element etc..)?
   */
  virtual ElementType at(PositionType) const = 0;

  /**
   * \brief Get the number of entities in the set
   * \return The number of entities in the set.
   */
  AXOM_HOST_DEVICE virtual PositionType size() const = 0;

  /**
   * \brief Determines if the Set is a Subset of another set.
   * \return true if the set is a subset of another set, otherwise false.
   */
  virtual bool isSubset() const = 0;

  /**
   * \brief Checks whether the set is valid.
   * \return true if the underlying indices are valid, false otherwise.
   */
  virtual bool isValid(bool verboseOutput = false) const = 0;

  /**
   * \brief Checks if there are any elements in the set -- equivalent to:
   * set.size() == 0
   */
  AXOM_HOST_DEVICE virtual bool empty() const = 0;

#if 0
  /**
   * \brief Returns true if the set contains the given element.
   *
   * Alternatively, we can return the position in the set containing the
   * element,
   * with some value for not containing the element
   */
  virtual bool          contains(const SetElement & elt) const = 0;

  //Possible other useful functions
  void                  reset(size_type) { throw NotImplementedException(); }
#endif

private:
  /**
   * \brief Utility function to verify that the given SetPosition is in a valid
   * range.
   */
  virtual void verifyPosition(PositionType) const = 0;
};

/**
 * \brief General equality operator for two sets.
 * \details Two sets are considered equal if they have the same number of
 * elements and their ordered indices agree.
 */
template <typename P1, typename E1, typename P2, typename E2>
inline bool operator==(const Set<P1, E1>& set1, const Set<P2, E2>& set2)
{
  using PosType = typename std::common_type<P1, P2>::type;
  using ElemType = typename std::common_type<E1, E2>::type;

  PosType const numElts = set1.size();

  // Sets are different if they have a different size
  if(set2.size() != numElts)
  {
    return false;
  }

  // Otherwise, compare the indices element wise
  for(PosType pos = PosType(); pos < numElts; ++pos)
  {
    auto&& e1 = static_cast<ElemType&&>(set1.at(static_cast<P1>(pos)));
    auto&& e2 = static_cast<ElemType&&>(set2.at(static_cast<P2>(pos)));
    if(e1 != e2)
    {
      return false;
    }
  }
  return true;
}
/**
 * \brief Set inequality operator
 */
template <typename P1, typename E1, typename P2, typename E2>
inline bool operator!=(const Set<P1, E1>& set1, const Set<P2, E2>& set2)
{
  return !(set1 == set2);
}

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_SET_H_
