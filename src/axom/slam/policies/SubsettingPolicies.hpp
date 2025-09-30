// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file SubsettingPolicies.hpp
 *
 * \brief Subsetting policies for SLAM
 *
 * Subsetting policies encompass the type and availability of a set's parent
 * A valid subset policy must support the following interface:
 *   * [required]
 *   * isSubset(): bool -- returns whether the set is a subset of another set
 *   * parentSet() : ParentSetType -- returns a pointer to the parent set.
 *                                     nullptr when isSubset() is false
 *   * isValid() : bool -- indicates whether the Subsetting policy of the set is valid
 *   * [optional]
 *   * operator(): IntType -- alternate accessor for indirection
 */

#ifndef SLAM_POLICIES_SUBSET_H_
#define SLAM_POLICIES_SUBSET_H_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

#include "axom/slam/NullSet.hpp"
#include "axom/export/slam.h"

#include <set>

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \name OrderedSet_Subsetting_Policies
 * \brief A few default policies for the subsetting of an OrderedSet
 */

/// \{

struct NoSubset
{
  AXOM_SLAM_EXPORT static const NullSet<> s_nullSet;
  using ParentSetType = const Set<>;

  AXOM_HOST_DEVICE NoSubset() { }

  // This empty .ctor is here to satisfy the SubsettingPolicy API
  NoSubset(ParentSetType*) { }

  /**
   * \brief Checks whether the set containing this policy class is a subset
   */
  bool isSubset() const { return false; }
  const ParentSetType* parentSet() const { return nullptr; }

  template <typename OrderedSetIt>
  bool isValid(OrderedSetIt, OrderedSetIt, bool) const
  {
    return true;
  }
};

struct VirtualParentSubset
{
  AXOM_SLAM_EXPORT static NullSet<> s_nullSet;

  using ParentSetType = Set<>;

  VirtualParentSubset(ParentSetType* parSet = &s_nullSet) : m_parentSet(parSet) { }

  /**
   * \brief Checks whether the set containing this policy class is a subset
   */
  bool isSubset() const { return *m_parentSet != s_nullSet; }
  const ParentSetType* parentSet() const { return m_parentSet; }
  ParentSetType*& parentSet() { return m_parentSet; }

  template <typename OrderedSetIt>
  bool isValid(OrderedSetIt beg, OrderedSetIt end, bool verboseOutput = false) const
  {
    // We allow parent sets to be null (i.e. the subset feature is deactivated)
    if(!isSubset() || m_parentSet == nullptr)
    {
      return true;
    }

    // Next, check if child is empty -- null set is a subset of all sets
    bool childIsEmpty = (beg == end);
    if(childIsEmpty)
    {
      return true;
    }

    // Next, since child has at least one element, the parent cannot be empty
    if(verboseOutput)
    {
      bool bValid = (m_parentSet->size() > 0);
      AXOM_UNUSED_VAR(bValid);
      SLIC_CHECK_MSG(bValid,
                     "VirtualParentSubset -- if we are a subset and input set is "
                       << "non-empty then parent set must be non-empty");
    }

    // At this point, parent and child are both non-null
    //  -- ensure that all elts of child are in parent
    using ElType = typename OrderedSetIt::value_type;
    std::set<ElType> pSet;
    for(auto pos = 0; pos < m_parentSet->size(); ++pos)
    {
      pSet.insert(m_parentSet->at(pos));
    }
    for(; beg != end; ++beg)
    {
      if(pSet.find(*beg) == pSet.end())
      {
        // For now, we only warn about the first failure '
        // Note: in the future, we might want to show all problems.
        SLIC_CHECK_MSG(verboseOutput,
                       "VirtualParentSubset :: parent set does not contain element "
                         << *beg << " so child cannot be a subset of parent");
        return false;
      }
    }
    return true;
  }

private:
  ParentSetType* m_parentSet;
};

template <typename TheParentSetType>
struct ConcreteParentSubset
{
  using ParentSetType = const TheParentSetType;

  ConcreteParentSubset(ParentSetType* parSet = nullptr) : m_parentSet(parSet) { }

  /**
   * \brief Checks whether the set containing this policy class is a subset
   */
  bool isSubset() const { return m_parentSet != nullptr; }

  ParentSetType* const& parentSet() const { return m_parentSet; }
  ParentSetType*& parentSet() { return m_parentSet; }

  template <typename OrderedSetIt>
  bool isValid(OrderedSetIt beg, OrderedSetIt end, bool verboseOutput = false) const
  {
    // We allow parent sets to be null (i.e. the subset feature is deactivated)
    if(!isSubset())
    {
      return true;
    }

    // Next, check if child is empty -- null set is a subset of all sets
    bool childIsEmpty = (beg == end);
    if(childIsEmpty)
    {
      return true;
    }

    // Next, since child has at least one element, the parent cannot be empty
    if(verboseOutput)
    {
      bool bValid = (m_parentSet->size() > 0);
      AXOM_UNUSED_VAR(bValid);
      SLIC_CHECK_MSG(
        bValid,
        "VirtualParentSubset -- if input set is non-empty " << " parent set must be non-empty");
    }

    // At this point, parent and child are both non-null
    std::set<typename ParentSetType::ElementType> pSet;
    for(auto pos = 0; pos < m_parentSet->size(); ++pos)
    {
      pSet.insert((*m_parentSet)[pos]);
    }
    for(; beg != end; ++beg)
    {
      if(pSet.find(*beg) == pSet.end())
      {
        // For now, we only warn about the first failure '
        // Note: in the future, we might want to show all problems.
        SLIC_CHECK_MSG(verboseOutput,
                       "ConcreteParentSubset :: parent set does not contain element "
                         << *beg << " so child cannot be a subset of parent");
        return false;
      }
    }
    return true;
  }

private:
  ParentSetType* m_parentSet;
};

/// \}

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_SUBSET_H_
