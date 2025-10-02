// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file DynamicVariableRelation.hpp
 *
 * \brief API for a topological relation between two sets in which entities from
 * the first set can be related to an arbitrary number of entities from the
 * second set. This relation is dynamic; the related entities can change
 * at runtime.
 */

#ifndef SLAM_DYNAMIC_VARIABLE_RELATION_HPP_
#define SLAM_DYNAMIC_VARIABLE_RELATION_HPP_

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/policies/PolicyTraits.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/Relation.hpp"

#include <vector>
#include <sstream>
#include <iterator>

namespace axom
{
namespace slam
{
template <typename FirstSetType = slam::Set<>, typename SecondSetType = slam::Set<>>
class DynamicVariableRelation
  : public Relation<typename FirstSetType::PositionType, typename FirstSetType::ElementType>
{
public:
  using FromSetType = FirstSetType;
  using ToSetType = SecondSetType;

  using SetPosition = typename FirstSetType::PositionType;
  using SetElement = typename FirstSetType::ElementType;

  using RelationVec = std::vector<SetPosition>;
  using RelationVecIterator = typename RelationVec::iterator;
  using RelationVecIteratorPair = std::pair<RelationVecIterator, RelationVecIterator>;
  using RelationVecConstIterator = typename RelationVec::const_iterator;
  using RelationVecConstIteratorPair = std::pair<RelationVecConstIterator, RelationVecConstIterator>;
  using RelationsContainer = std::vector<RelationVec>;
  using RelationsContainerCIt = typename RelationsContainer::const_iterator;
  using RelationsContainerIt = typename RelationsContainer::iterator;

public:
  DynamicVariableRelation(FirstSetType* fromSet = policies::EmptySetTraits<FirstSetType>::emptySet(),
                          SecondSetType* toSet = policies::EmptySetTraits<SecondSetType>::emptySet())
    : m_fromSet(fromSet)
    , m_toSet(toSet)
  {
    if(m_fromSet)
    {
      m_relationsVec.resize(m_fromSet->size());
    }
  }

  ~DynamicVariableRelation() { }

public:
  /// \name DynamicVariableRelation iterator interface
  /// @{
  RelationVecConstIterator begin(SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).begin();
  }

  RelationVecConstIterator end(SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).end();
  }

  RelationVecConstIteratorPair range(SetPosition fromSetIndex) const
  {
    return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
  }
  /// @}

  RelationVec const& operator[](SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return m_relationsVec[fromSetIndex];
  }

  SetPosition size(SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).size();
  }

  SetPosition totalSize() const
  {
    SetPosition sz = 0;
    for(auto& vec : m_relationsVec)
    {
      sz += vec.size();
    }
    return sz;
  }

  bool hasFromSet() const { return !policies::EmptySetTraits<FromSetType>::isEmpty(m_fromSet); }
  FromSetType* fromSet() { return m_fromSet; }
  const FromSetType* fromSet() const { return m_fromSet; }

  bool hasToSet() const { return !policies::EmptySetTraits<ToSetType>::isEmpty(m_toSet); }
  ToSetType* toSet() { return m_toSet; }
  const ToSetType* toSet() const { return m_toSet; }

  SetPosition fromSetSize() const { return m_relationsVec.size(); }

  bool isValid(bool verboseOutput = false) const;

public:  // Modifying functions
  void insert(SetPosition fromSetIndex, SetPosition toSetIndex)
  {
    verifyPosition(fromSetIndex);
    m_relationsVec[fromSetIndex].push_back(toSetIndex);
  }

  RelationVec& operator[](SetPosition fromSetIndex)
  {
    verifyPosition(fromSetIndex);
    return m_relationsVec[fromSetIndex];
  }

public:
  /**
   * \name DirectDataAccess
   * \brief Accessor functions to get the underlying relation data for each
   *  element

   * \note We will have to figure out a good way
   * to limit this access to situations where it makes sense.
   */

  /// \{

  /**
   * \brief Access the set of positions in the 'toSet'
   * associated with the given position in 'fromSet'
   *
   * \param fromSetPos The position within the 'fromSet'
   * whose relation data (in the 'toSet') we are requesting
   */
  RelationVec& data(SetPosition fromSetPos)
  {
    verifyPosition(fromSetPos);
    return m_relationsVec[fromSetPos];
  }

  /**
   * \brief Access the set of positions in the 'toSet'
   * associated with the given position in 'fromSet'
   *
   * \param fromSetPos The position within the 'fromSet'
   * whose relation data (in the 'toSet') we are requesting
   */
  const RelationVec& data(SetPosition fromSetPos) const
  {
    verifyPosition(fromSetPos);
    return m_relationsVec[fromSetPos];
  }

  /// \}

private:
  inline void verifyPosition(SetPosition AXOM_DEBUG_PARAM(fromSetIndex)) const
  {
    SLIC_ASSERT_MSG(fromSetIndex >= 0 && fromSetIndex < static_cast<SetPosition>(m_fromSet->size()),
                    "Index " << fromSetIndex << " out of range [0," << m_fromSet->size() << ")");
  }

  inline RelationVec& fromSetRelationsVec(SetPosition fromSetIndex)
  {
    return m_relationsVec[fromSetIndex];
  }
  inline RelationVec const& fromSetRelationsVec(SetPosition fromSetIndex) const
  {
    return m_relationsVec[fromSetIndex];
  }

private:
  FromSetType* m_fromSet;
  ToSetType* m_toSet;

  RelationsContainer m_relationsVec;
};

template <typename FirstSetType, typename SecondSetType>
bool DynamicVariableRelation<FirstSetType, SecondSetType>::isValid(bool verboseOutput) const
{
  bool bValid = true;

  std::stringstream sstr;

  if(!hasFromSet() || !hasToSet())
  {
    if(!m_relationsVec.empty())
    {
      if(verboseOutput)
      {
        sstr << "\n\t* relations vector was not empty "
             << " -- fromSet was " << (!hasFromSet() ? "" : " not ") << "null"
             << " , toSet was " << (!hasToSet() ? "" : " not ") << "null";
      }

      bValid = false;
    }
  }
  else
  {
    if(verboseOutput)
    {
      sstr << "\n\t* Neither set was null";
    }

    // Check that the the relations vector has the right size
    // (should be same as fromSet's size() )
    if(static_cast<SetPosition>(m_relationsVec.size()) != m_fromSet->size())
    {
      if(verboseOutput)
      {
        sstr << "\n\t* relations vector has the wrong size."
             << "\n\t-- from set size is: " << m_fromSet->size()
             << "\n\t-- expected relation size: " << m_fromSet->size()
             << "\n\t-- actual size: " << m_relationsVec.size();
      }
      bValid = false;
    }

    // Check that all elements of the relations vector point to
    // valid  set elements in the toSet
    for(SetPosition fromIdx = 0; fromIdx < m_fromSet->size(); ++fromIdx)
    {
      SetPosition idx = fromIdx;
      for(RelationVecConstIterator rIt = begin(idx), rEnd = end(idx); rIt < rEnd; ++rIt)
      {
        if(*rIt >= m_toSet->size())
        {
          if(verboseOutput)
          {
            sstr << "\n\t* relation for element " << m_fromSet->at(fromIdx)
                 << " of fromSet had an out-of-range element.-- value "
                 << std::distance(begin(idx), rIt) << " was " << *rIt
                 << ". Max possible value should be " << m_toSet->size() << ".";
          }
          bValid = false;
        }
      }
    }
  }

  if(verboseOutput)
  {
    std::stringstream sstr2;
    sstr2 << "\n*** Detailed results of isValid on the relation.\n";
    if(bValid)
    {
      sstr2 << "(dynamic,variable) Relation was valid." << std::endl;
    }
    else
    {
      sstr2 << "Relation was NOT valid.\n" << sstr.str() << std::endl;
    }

    if(m_fromSet)
    {
      sstr2 << "\n** fromSet has size " << m_fromSet->size() << ": ";
    }
    if(m_toSet)
    {
      sstr2 << "\n** toSet has size " << m_toSet->size() << ": ";
    }

    if(m_relationsVec.empty())
    {
      sstr2 << "\n** relations vec is empty:";
    }
    else
    {
      SetPosition overallCount = 0;
      sstr2 << "\n** relations vec elements:";

      for(SetPosition fromIdx = 0; fromIdx < m_fromSet->size(); ++fromIdx)
      {
        SetPosition idx = fromIdx;
        sstr2 << "\n\t" << m_fromSet->at(fromIdx) << " (" << size(idx) << "):\t";
        std::copy(begin(idx), end(idx), std::ostream_iterator<SetPosition>(sstr2, " "));
        overallCount += size(idx);
      }
      sstr2 << "\n\n\tOverall size of relation" << overallCount << std::endl;

      SLIC_INFO(sstr2.str());
    }
  }

  return bValid;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_DYNAMIC_VARIABLE_RELATION_HPP_
