// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file DynamicConstantRelation.hpp
 *
 * \brief API for a topological relation between two sets in which entities
 * from the first set can be related to a constant number of entities from
 * the second set. For example, in a triangle mesh, each triangle is
 * incident to three vertices.
 *
 * This relation is dynamic; the related entities can change at runtime.
 */

#ifndef SLAM_DYNAMIC_CONSTANT_RELATION_HPP_
#define SLAM_DYNAMIC_CONSTANT_RELATION_HPP_

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/Relation.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include "axom/fmt.hpp"

#include <vector>

namespace axom
{
namespace slam
{
/**
 * \class DynamicConstantRelation
 * \brief  A relation class with constant cardinality that supports
 * adding, removing and modifying set relations.
 *
 * A DynamicConstantRelation encodes the relation between two sets,
 * A FromSet and a ToSet, where the cardinality of the relation from
 * each element of the FromSet to the ToSet is fixed to a constant value.
 * For example, each triangle in the triangle set of a triangle mesh
 * has three incident vertices from set of vertices.
 *
 * The relation from an element of the FromSet to an element
 * of the ToSet is considered to be valid if its entry in the FromSet
 * is valid and at least one of its relation entities in the ToSet is valid
 * (i.e. not equal to INVALID_INDEX).
 *
 * \note The current implementation fixes the value of INVALID_INDEX.
 * A future update will allow users to set the value of INVALID_INDEX to a
 * more convenient value, when necessary.
 */
template <typename PosType,   //= slam::DefaultPositionType,
          typename ElemType,  // = slam::DefaultElementType,
          typename CardinalityPolicy>
class DynamicConstantRelation : public /*Relation,*/ CardinalityPolicy
{
public:
  enum
  {
    INVALID_INDEX = ~0  ///< value to mark indices of deleted elements
  };

  using SetPosition = PosType;
  using SetElement = ElemType;
  using RelationVec = std::vector<SetPosition>;

  using FromSetType = DynamicSet<PosType, ElemType>;
  using ToSetType = DynamicSet<PosType, ElemType>;

  using BeginsSizePolicy = typename CardinalityPolicy::RelationalOperatorSizeType;

  using STLIndirection = policies::STLVectorIndirection<SetPosition, SetElement>;
  using RelationSubset = OrderedSet<SetPosition,
                                    SetElement,
                                    BeginsSizePolicy,
                                    policies::RuntimeOffset<SetPosition>,
                                    policies::StrideOne<SetPosition>,
                                    STLIndirection>;

  // types for iterator
  using RelationIterator = typename RelationSubset::iterator;
  using RelationIteratorPair = typename RelationSubset::iterator_pair;

  using RelationConstIterator = typename RelationSubset::const_iterator;
  using RelationConstIteratorPair = typename RelationSubset::const_iterator_pair;

public:
  /**
   * \brief Default constructor with empty set for toSet and fromSet
   */
  DynamicConstantRelation()
    : m_fromSet(policies::EmptySetTraits<FromSetType>::emptySet())
    , m_toSet(policies::EmptySetTraits<ToSetType>::emptySet())
  { }

  /**
   * \brief Construct a DynamicConstantRelation from the given \a fromSet
   * to \a toSet
   */
  DynamicConstantRelation(FromSetType* fromSet, ToSetType* toSet)
    : CardinalityPolicy(policies::EmptySetTraits<FromSetType>::isEmpty(fromSet) ? 0 : fromSet->size())
    , m_fromSet(fromSet)
    , m_toSet(toSet)
    , m_currentFromSize(fromSet == nullptr ? 0 : m_fromSet->size())
  {
    updateSizes();
  };

public:
  /// \name DynamicConstantRelation iterator interface
  /// @{

  /**
   * \brief Returns a begin iterator to the set of entities in the ToSet
   * that are related to the element with index \a fromSetInd in the FromSet
   *
   * \param fromSetInd The index of the element in the FromSet
   * \return A begin iterator to the set of related elements in ToSet
   */
  RelationIterator begin(SetPosition fromSetInd)
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].begin();
  }

  /**
   * \brief Returns a begin const iterator to the set of entities in the ToSet
   * that are related to the element with index \a fromSetInd in the FromSet
   *
   * \param fromSetInd The index of the element in the FromSet
   * \return A const begin iterator to the set of related elements in ToSet
   */
  RelationConstIterator begin(SetPosition fromSetInd) const
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].begin();
  }

  /**
   * \brief Returns an end iterator to the set of entities in the ToSet
   * that are related to the element with index \a fromSetInd in the FromSet
   *
   * \param fromSetInd The index of the element in the FromSet
   * \return An end iterator to the set of related elements in ToSet
   */
  RelationIterator end(SetPosition fromSetInd)
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].end();
  }

  /**
   * \brief Returns a end const iterator to the set of entities in the ToSet
   * that are related to the element with index \a fromSetInd in the FromSet
   *
   * \param fromSetInd The index of the element in the FromSet
   * \return A const end iterator to the set of related elements in ToSet
   */
  RelationConstIterator end(SetPosition fromSetInd) const
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].end();
  }

  /**
   * \brief Returns an iterator range to the set of entities in the ToSet
   * that are related to the element with index \a fromSetInd in the FromSet
   *
   * \param fromSetInd The index of the element in the FromSet
   * \return An iterator range (begin/end pair) to the set of related
   * elements in ToSet
   */
  RelationIteratorPair range(SetPosition fromSetInd) { return (*this)[fromSetInd].range(); }

  /**
   * \brief Returns a const iterator range to the set of entities in the ToSet
   * that are related to the element with index \a fromSetInd in the FromSet
   *
   * \param fromSetInd The index of the element in the FromSet
   * \return A const iterator range (begin/end pair) to the set of related
   * elements in ToSet
   */
  RelationConstIteratorPair range(SetPosition fromSetInd) const
  {
    return (*this)[fromSetInd].range();
  }

  /// @}

public:
  /// \name DynamicConstantRelation per-element relation access functions
  /// @{
  ///

  /**
   * \brief Returns the const set of entities in the ToSet related to the
   * element with index \a fromSetIndex in the FromSet
   * \param fromSetIndex The index of an element in the FromSet
   */
  RelationSubset const at(SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return operator[](fromSetIndex);
  }

  RelationSubset at(SetPosition fromSetIndex)
  {
    verifyPosition(fromSetIndex);
    return operator[](fromSetIndex);
  }

  /**
   * \brief Returns the const set of entities in the ToSet related to the
   * element with index \a fromSetIndex in the FromSet
   * \param fromSetIndex The index of an element in the FromSet
   * \note This function does not modify the size of the relation. 
   * Use updateSizes(), insert(), or modify() functions to change sizes
   */
  RelationSubset const operator[](SetPosition fromSetIndex) const
  {
    // NOTE: Need to const_cast the pointer to the vector
    // since SetBuilder, and the IndirectionPolicy don't
    // currently support const buffers
    // TODO: Fix this!

    verifyPosition(fromSetIndex);
    using SetBuilder = typename RelationSubset::SetBuilder;
    return SetBuilder()
      .size(relationCardinality())
      //.offset( CardinalityPolicy::offset( fromSetIndex) )
      .offset(fromSetIndex * relationCardinality())
      .data(const_cast<RelationVec*>(&m_relationsVec));
  }

  RelationSubset operator[](SetPosition fromSetIndex)
  {
    verifyPosition(fromSetIndex);
    using SetBuilder = typename RelationSubset::SetBuilder;
    return SetBuilder()
      .size(relationCardinality())
      //.offset( CardinalityPolicy::offset( fromSetIndex) )
      .offset(fromSetIndex * relationCardinality())
      .data(&m_relationsVec);
  }

  /**
   * \brief Returns the cardinality of the set of entities in the ToSet
   * related to the element with index \a fromSetIndex in the FromSet
   * \param fromSetIndex The index of an element in the FromSet
   */
  SetPosition size(SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return relationCardinality();
  }

  /// @}

  /// \brief Returns the cardinality of the FromSet
  inline SetPosition size() const { return m_currentFromSize; }

public:
  /// \name DynamicConstantRelation validity check functions
  /// @{
  ///

  /**
   * \brief Returns the number of valid entries in the FromSet
   *
   * An element of the FromSet is considered valid with respect to a
   * DynamicConstantRelation when it is valid in the FromSet and when
   * its relation set is not marked as invalid.
   * \sa isValidEntry()
   */
  SetPosition numberOfValidEntries() const
  {
    SetPosition nvalid = 0;
    const int N = size();
    for(int i = 0; i < N; ++i)
    {
      nvalid += isValidEntry(i);
    }
    return nvalid;
  }

  /**
   * \brief return if an entry is valid or not.
   * \details an entry is considered valid if it is valid in the from set and has at least one valid value
   */
  bool isValidEntry(SetPosition idx) const
  {
    if(m_fromSet->isValidEntry(idx))
    {
      const auto SZ = relationCardinality();
      const auto beg_idx = idx * SZ;
      for(auto i = beg_idx; i < (beg_idx + SZ); ++i)
      {
        if(m_relationsVec[i] != INVALID_INDEX)
        {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * \brief Predicate to check if the DynamicConstantRelation instance is valid
   */
  bool isValid(bool verboseOutput = false) const;

  /// @}

public:
  /// \name DynamicConstantRelation functions that modify the relation
  /// @{
  ///

  /**
   * \brief Inserts a new entry into the relation at the first invalid index
   * \param fromSetIndex The index of the element in the FromSet
   * \param toSetIndex The index of the element in the ToSet
   * to associate with \a fromSetIndex
   */
  void insert(SetPosition fromSetIndex, SetPosition toSetIndex)
  {
    expandSizeIfNeeded(fromSetIndex + 1);
    verifyPosition(fromSetIndex);

    //find the first invalid place to put it
    const auto SZ = relationCardinality();
    for(int i = 0; i < SZ; ++i)
    {
      const auto idx = SZ * fromSetIndex + i;
      if(m_relationsVec[idx] == INVALID_INDEX)
      {
        m_relationsVec[idx] = toSetIndex;
        return;
      }
    }

    //The entry was not inserted
    SLIC_WARNING("Relation from " << fromSetIndex << " to " << toSetIndex
                                  << " was not inserted because the entry is full.");
  }

  /**
   * \brief Function to modify the value at offset \a offset of the
   * FromSet index \fromSetIndex to the value \a toSetIndex
   *
   * \note This is a temporary function until operator[]
   *  allows us to modify values.
   *
   * This should be replaced with operator[] which returns a non-const
   * RelationSubset so users can more naturally update the relation.
   * E.g. relation[fromSetIndex][offset] = toSetIndex;
   */
  void modify(SetPosition fromSetIndex, SetPosition offset, SetPosition toSetIndex)
  {
    expandSizeIfNeeded(fromSetIndex + 1);
    verifyPosition(fromSetIndex);
    m_relationsVec[relationCardinality() * fromSetIndex + offset] = toSetIndex;
  }

  /// \brief Mark all values in entry \a fromSetIndex as invalid.
  void remove(SetPosition fromSetIndex)
  {
    if(!isValidEntry(fromSetIndex))
    {
      return;
    }

    const auto SZ = relationCardinality();
    for(int i = 0; i < SZ; ++i)
    {
      m_relationsVec[SZ * fromSetIndex + i] = INVALID_INDEX;
    }
  }

  void updateSizes()
  {
    m_currentFromSize = m_fromSet->size();
    m_relationsVec.resize(m_currentFromSize * relationCardinality(), INVALID_INDEX);
  }

  /// @}

public:
  /** \brief Direct access to the relation data  */
  RelationVec& data() { return m_relationsVec; }

  /** \brief Direct const access to the relation data  */
  const RelationVec& data() const { return m_relationsVec; }

private:
  inline constexpr SetPosition relationCardinality() const
  {
    return CardinalityPolicy::size(SetPosition());
  }

  /**
   * \brief Helper function to expand the relation data storage
   * \param s The requested size
   */
  void expandSizeIfNeeded(SetPosition s)
  {
    if(s > m_currentFromSize)
    {
      m_currentFromSize = m_fromSet->size();
      m_relationsVec.resize(m_currentFromSize * relationCardinality(), INVALID_INDEX);
    }

    SLIC_ASSERT_MSG(s <= m_currentFromSize,
                    fmt::format("Expanded size {} is larger than relation's 'from' set of {}",
                                s,
                                m_fromSet->size()));
  }

  /**
   * \brief Debug check that an index in the FromSet is not out-of-range
   * \param fromSetIndex An (alleged) index in the FromSet
   */
  inline void verifyPosition(SetPosition AXOM_DEBUG_PARAM(fromSetIndex)) const
  {
    SLIC_ASSERT_MSG(fromSetIndex >= 0 && fromSetIndex < m_currentFromSize,
                    fmt::format("Index {} out of range [0,{})", fromSetIndex, m_currentFromSize));
  }

private:
  FromSetType* m_fromSet;
  ToSetType* m_toSet;

  RelationVec m_relationsVec;
  IndexType m_currentFromSize {0};
};

/* Checks whether the relation is valid.  */
template <typename PosType, typename ElemType, typename CardinalityPolicy>
bool DynamicConstantRelation<PosType, ElemType, CardinalityPolicy>::isValid(bool verboseOutput) const
{
  fmt::memory_buffer out;

  bool setsAreValid = true;
  bool relationdataIsValid = true;

  // Check if the sets are valid
  const bool isFromSetNull = (m_fromSet == nullptr);
  const bool isToSetNull = (m_toSet == nullptr);

  if(isFromSetNull || isToSetNull)
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\t Static relations require both the fromSet and toSet "
                     "to be non-null:"
                     "\t-- fromSet was {} null"
                     "\t-- toSet was {} null",
                     isFromSetNull ? "" : " not ",
                     isToSetNull ? "" : " not ");
    }

    setsAreValid = false;
  }

  // Check the sizes of fromSet matches relationVec
  if(setsAreValid)
  {
    if(m_fromSet->size() != m_currentFromSize)
    {
      if(verboseOutput)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\t Internal size does not match fromSet size:"
                       "\t-- fromSet size is {}"
                       "\t-- internal size is {}",
                       m_fromSet->size(),
                       m_currentFromSize);
      }
      setsAreValid = false;
    }

    if(m_fromSet->size() * relationCardinality() != (int)m_relationsVec.size())
    {
      if(verboseOutput)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\t Size of relationVec does not match toSet size:"
                       "\t-- fromSet size is {}"
                       "\t-- m_relationsVec size is {}",
                       m_fromSet->size(),
                       m_relationsVec.size());
      }
      setsAreValid = false;
    }
  }

  // Check if the relation data is valid
  if(setsAreValid)
  {
    // Check that invalid set entry points to invalid relation.
    // Note: the reverse can be valid. ie. valid set entry may have invalid relation entry.
    for(auto pos : m_fromSet->positions())
    {
      if(m_fromSet->at(pos) == FromSetType::INVALID_ENTRY && isValidEntry(pos))
      {
        if(verboseOutput)
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\t* invalid entries in fromSet; has a valid relation at index "
                         "{}, but element not in from set. Values: {}",
                         pos,
                         (*this)[pos]);
        }
        relationdataIsValid = false;
      }
    }

    // Check that all relation indices are in range for m_toSet
    for(auto from_idx : m_fromSet->positions())
    {
      if(m_fromSet->isValidEntry(from_idx))
      {
        for(auto idx = 0; idx < relationCardinality(); ++idx)
        {
          const auto pos = from_idx * relationCardinality() + idx;
          const auto val = m_relationsVec[pos];
          if(val != INVALID_INDEX && !m_toSet->isValidEntry(val))
          {
            if(verboseOutput)
            {
              fmt::format_to(std::back_inserter(out),
                             "\n\t* Relation index out of range or invalid:"
                             "\n\t-- position {} ({}-{}) with value {} needs "
                             "to be in range [0,{}) and index a valid entry",
                             pos,
                             from_idx,
                             idx,
                             val,
                             m_toSet->size());
            }
            relationdataIsValid = false;
          }
        }
      }
    }
  }

  // We are done.  Output the messages if applicable and return
  bool bValid = setsAreValid && relationdataIsValid;

  if(verboseOutput && !bValid)
  {
    SLIC_INFO(fmt::to_string(out));
  }

  return bValid;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_DYNAMIC_CONSTANT_RELATION_HPP_
