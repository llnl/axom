// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Relation.hpp
 *
 * \brief Basic API for a topological relation between two sets
 *
 */

#ifndef SLAM_RELATION_HPP_
#define SLAM_RELATION_HPP_

#include <vector>

#include "axom/slam/Set.hpp"
#include "axom/slam/NullSet.hpp"

namespace axom
{
namespace slam
{
template <typename PosType, typename ElemType>
class NullSet;

template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
class Relation
{
public:
  using SetPosition = typename Set<PosType, ElemType>::PositionType;
  using SetElement = typename Set<PosType, ElemType>::ElementType;

  using RelationVec = std::vector<SetPosition>;
  using RelationVecIterator = typename RelationVec::iterator;
  using RelationVecIteratorPair = std::pair<RelationVecIterator, RelationVecIterator>;
  using RelationVecConstIterator = typename RelationVec::const_iterator;
  using RelationVecConstIteratorPair = std::pair<RelationVecConstIterator, RelationVecConstIterator>;

  static NullSet<PosType, ElemType> s_nullSet;

public:
  virtual ~Relation() { }

  virtual RelationVecConstIterator begin(SetPosition fromSetIndex) const = 0;

  virtual RelationVecConstIterator end(SetPosition fromSetIndex) const = 0;

  virtual RelationVecConstIteratorPair range(SetPosition fromSetIndex) const = 0;

  virtual SetPosition size(SetPosition fromSetIndex) const = 0;

  virtual bool isValid(bool verboseOutput = false) const = 0;

#if 0
  // Go through the relation's data and ensure that no entity from the ToSet
  // is mapped to more than once by an element of the FromSet.
  // This operation will compact the relation
  // We can optionally sort the entities as well...
  // Not yet implemented
  virtual void        makeUnique() = 0;

  // Accessors to the underlying sets -- allows setting and getting the sets
  virtual Set* fromSet()       = 0;
  virtual const Set* fromSet() const = 0;
  virtual Set* toSet()         = 0;
  virtual const Set* toSet()   const = 0;

  // This function differs for each concrete relation type..
  void bindRelationData(..args..) = 0;

  // Accessors to the underlying relation data
  // -- differs depending on the implementation (e.g. structured vs.
  // unstructured
  ArrType* relationData() = 0;
#endif
};

/**
 * \brief Definition of static instance of nullSet for all relations
 * \note Should this be a singleton or a global object?  Should the scope be
 *  public?
 */
template <typename PosType, typename ElemType>
NullSet<PosType, ElemType> Relation<PosType, ElemType>::s_nullSet;

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_RELATION_HPP_
