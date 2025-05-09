// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_MAPPED_RELATION_SET_H_
#define SLAM_MAPPED_RELATION_SET_H_

#include "axom/slam/RangeSet.hpp"
#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/policies/BivariateSetInterfacePolicies.hpp"

namespace axom
{
namespace slam
{
/**
 * \class RelationSet
 *
 * \brief Models a Set whose elements are derived from a relation, one element
 *        per fromSet and toSet pair in the relation.
 *
 *  RelationSet models a subset of the Cartesian product of two sets. Users
 *  should refer to the BivariateSet documentation for descriptions of the
 *  different indexing names (SparseIndex, DenseIndex, FlatIndex).
 *
 * \tparam  Relation  The Relation type that this set uses.
 *
 * \see   BivariateSet
 */

template <typename Relation,
          typename SetType1 = typename Relation::FromSetType,
          typename SetType2 = typename Relation::ToSetType,
          typename InterfaceType = policies::VirtualInterface>
class RelationSet final : public policies::BivariateSetInterface<InterfaceType, SetType1, SetType2>
{
public:
  using FirstSetType = SetType1;
  using SecondSetType = SetType2;

  using RelationType = Relation;

private:
  using BaseType = policies::BivariateSetInterface<InterfaceType, SetType1, SetType2>;
  using RangeSetType = typename BaseType::RangeSetType;
  using BaseSubsetType = typename BaseType::SubsetType;

public:
  using PositionType = typename RelationType::SetPosition;
  using ElementType = typename RelationType::SetElement;

  using RelationSubset = typename RelationType::RelationSubset;
  using SubsetType =
    std::conditional_t<std::is_same<void, BaseSubsetType>::value, RelationSubset, BaseSubsetType>;

  using BaseType::INVALID_POS;

  using IteratorType = BivariateSetIterator<RelationSet>;

public:
  using ConcreteSet = RelationSet<Relation, SetType1, SetType2, policies::ConcreteInterface>;
  using VirtualSet = RelationSet<Relation, SetType1, SetType2, policies::VirtualInterface>;

  using OtherSet =
    std::conditional_t<std::is_same<InterfaceType, policies::VirtualInterface>::value, ConcreteSet, VirtualSet>;

  RelationSet(const OtherSet& other)
    : BaseType(other.getFirstSet(), other.getSecondSet())
    , m_relation(other.getRelation())
  { }

public:
  RelationSet() = default;

  /**
   * \brief Constructor taking in the relation this BivariateSet is based on.
   * \pre relation pointer must not be a null pointer
   */
  RelationSet(RelationType* relation)
    : BaseType(relation ? relation->fromSet() : policies::EmptySetTraits<FirstSetType>::emptySet(),
               relation ? relation->toSet() : policies::EmptySetTraits<SecondSetType>::emptySet())
    , m_relation(relation)
  {
    SLIC_ASSERT(relation != nullptr);
  }

  /**
   * \brief Searches for the SparseIndex of the element given its DenseIndex.
   * \detail If the element (i,j) is the k<sup>th</sup> non-zero in the row,
   *         then `findElementIndex(i,j)` returns `k`. If `element(i,j)` does
   *         not exist (such as the case of a zero in a sparse matrix), then
   *         `INVALID_POS` is returned.
   *
   * \warning This function can be slow, since a linear search is performed on
   *          the row each time.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The DenseIndex of the given element, or INVALID_POS if such
   *          element is missing from the set.
   * \pre   0 <= pos1 <= set1.size() && 0 <= pos2 <= size2.size()
   */

  PositionType findElementIndex(PositionType pos1, PositionType pos2) const
  {
    RelationSubset ls = (*m_relation)[pos1];
    for(PositionType i = 0; i < ls.size(); i++)
    {
      if(ls[i] == pos2)
      {
        return i;
      }
    }
    return BaseType::INVALID_POS;
  }

  /**
   * \brief Search for the FlatIndex of the element given its DenseIndex.
   * \warning This function can be slow, since a linear search is performed on
   *          the row each time.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The element's FlatIndex
   * \pre   0 <= pos1 <= set1.size() && 0 <= pos2 <= size2.size()
   */
  AXOM_HOST_DEVICE PositionType findElementFlatIndex(PositionType s1, PositionType s2) const
  {
    RelationSubset ls = (*m_relation)[s1];
    for(PositionType i = 0; i < ls.size(); i++)
    {
      if(ls[i] == s2)
      {
        return ls.offset() + i;
      }
    }
    return BaseType::INVALID_POS;
  }

  /**
   * \brief Given the from-set index pos1, return the FlatIndex of the first
   *        existing to-set element in the relation pair, or `INVALID_POS` if
   *        this row contains no elements.
   *
   * \param pos1  Index into the from-set.
   * \param pos2  Index into the to-set.
   *
   * \return  The FlatIndex of the first existing to-set element.
   */
  PositionType findElementFlatIndex(PositionType pos1) const
  {
    RelationSubset ls = (*m_relation)[pos1];

    if(ls.size() > 0)
    {
      return ls.offset();
    }

    return BaseType::INVALID_POS;
  }

  /**
   * \brief Given the flat index, return the associated to-set index in the
   *        relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos2  The to-set index.
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE PositionType flatToSecondIndex(PositionType flatIndex) const
  {
    if(flatIndex < 0 || flatIndex > size())
    {
      SLIC_ASSERT("Flat index out of bounds of the relation set.");
    }
    return m_relation->relationData()[flatIndex];
  }

  /**
   * \brief Given the flat index, return the associated from-set index in the
   *        relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos1  The from-set index.
   */
  AXOM_HOST_DEVICE PositionType flatToFirstIndex(PositionType flatIndex) const
  {
    if(flatIndex < 0 || flatIndex > size())
    {
      SLIC_ASSERT("Flat index out of bounds of the relation set.");
    }
    return m_relation->firstIndex(flatIndex);
  }

  AXOM_HOST_DEVICE RangeSetType elementRangeSet(PositionType pos1) const
  {
    return
      typename RangeSetType::SetBuilder().size(m_relation->size(pos1)).offset(m_relation->offset(pos1));
  }

  /**
   * \brief A set of elements with the given first set index.
   *
   * \param s1  The first set index.
   * \return  An OrderedSet containing the elements in the row.
   * \pre  0 <= pos1 <= set1.size()
   */
  SubsetType getElements(PositionType s1) const { return (*m_relation)[s1]; }

  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE ElementType at(PositionType pos) const
  {
#ifndef AXOM_DEVICE_CODE
    RelationSet::verifyPosition(pos);
#endif
    return m_relation->relationData()[pos];
  }

  /** \brief Returns the relation pointer   */
  RelationType* getRelation() const { return m_relation; }

  RelationType* getRelation() { return m_relation; }

  /** \brief Return the size of the relation   */
  PositionType totalSize() const { return PositionType(m_relation->relationData().size()); }

  /**
   * \brief Return the size of a row, which is the number of to-set
   *        elements associated with the given from-set index.
   *
   * \param pos The from-set position.
   */
  PositionType size(PositionType pos) const { return m_relation->size(pos); }

  /*!
   * \brief Return an iterator to the first pair of set elements in the
   *  relation.
   */
  IteratorType begin() const { return IteratorType(this, 0); }

  /*!
   * \brief Return an iterator to one past the last pair of set elements in the
   *  relation.
   */
  IteratorType end() const { return IteratorType(this, totalSize()); }

  bool isValid(bool verboseOutput = false) const
  {
    if(m_relation == nullptr)
    {
      if(verboseOutput)
      {
        std::cout << "\n*** RelationSet is not valid:\n"
                  << "\t* Relation pointer should not be null.\n"
                  << std::endl;
      }
      return false;
    }
    return m_relation->isValid(verboseOutput);
  }

public:
  //hiding size() from the Set base class, replaced with totalSize().
  //but still implemented due to the function being virtual
  //(and can be called from base ptr)
  // KW -- made this public to use from BivariateMap
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE PositionType size() const
  {
    return PositionType(m_relation->relationData().size());
  }

private:
  //range check only
  bool isValidIndex(PositionType s1, PositionType s2) const
  {
    return s1 >= 0 && s1 < m_relation->fromSet()->size() && s2 >= 0 && s2 < m_relation->size(s1);
  }

  void verifyPosition(PositionType AXOM_DEBUG_PARAM(sPos)) const
  {
    SLIC_ASSERT_MSG(sPos >= 0 && sPos < size(),
                    "SLAM::RelationSet -- requested out-of-range element at position "
                      << sPos << ", but set only has " << size() << " elements.");
  }

  void verifyPosition(PositionType AXOM_DEBUG_PARAM(s1), PositionType AXOM_DEBUG_PARAM(s2)) const
  {
    SLIC_ASSERT_MSG(isValidIndex(s1, s2),
                    "SLAM::RelationSet -- requested out-of-range element at position ("
                      << s1 << "," << s2 << "), but set only has " << this->firstSetSize() << "x"
                      << this->secondSetSize() << " elements.");
  }

private:
  RelationType* m_relation;  //the relation that this set is based off of
};

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_MAPPED_RELATION_SET_H_
