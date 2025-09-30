// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_BivarSetIfacePolicies_HPP
#define SLAM_BivarSetIfacePolicies_HPP

#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/policies/InterfacePolicies.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
namespace detail
{
/*!
 * \brief The base class of concrete-interface bivariate sets.
 *
 *  This is effectively slam::BivariateSet without the virtual functions.
 */
template <typename Set1 = slam::Set<>, typename Set2 = slam::Set<>>
struct ConcreteBivariateSet
{
  using FirstSetType = Set1;
  using SecondSetType = Set2;

  using PositionType = typename FirstSetType::PositionType;
  using ElementType = typename FirstSetType::ElementType;

  using SubsetType = void;
  using RangeSetType = typename RangeSet<PositionType, ElementType>::ConcreteSet;

  static constexpr PositionType INVALID_POS = PositionType(-1);

  /**
   * \brief Constructor taking pointers to the two sets that defines the range
   *        of the indices of the BivariateSet.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */
  ConcreteBivariateSet(const Set1* set1 = policies::EmptySetTraits<Set1>::emptySet(),
                       const Set2* set2 = policies::EmptySetTraits<Set2>::emptySet())
    : m_set1(set1)
    , m_set2(set2)
  { }

  /** \brief Size of the first set.   */
  AXOM_HOST_DEVICE inline PositionType firstSetSize() const
  {
    return getSize<FirstSetType>(m_set1);
  }
  /** \brief Size of the second set.   */
  AXOM_HOST_DEVICE inline PositionType secondSetSize() const
  {
    return getSize<SecondSetType>(m_set2);
  }

  /** \brief Returns pointer to the first set.   */
  const FirstSetType* getFirstSet() const { return m_set1; }
  /** \brief Returns pointer to the second set.   */
  const SecondSetType* getSecondSet() const { return m_set2; }

  bool isValid(bool verboseOutput = false) const
  {
    if(getFirstSet() == nullptr || getSecondSet() == nullptr)
    {
      if(verboseOutput)
      {
        SLIC_INFO("BivariateSet is not valid: " << " Set pointers should not be null.");
      }
      return false;
    }
    return getFirstSet()->isValid(verboseOutput) && getSecondSet()->isValid(verboseOutput);
  }

private:
  template <typename SetType>
  typename std::enable_if<std::is_abstract<SetType>::value, PositionType>::type getSize(
    const SetType* s) const
  {
    SLIC_ASSERT_MSG(s != nullptr, "nullptr in BivariateSet::getSize()");
    return s->size();
  }

  template <typename SetType>
  AXOM_HOST_DEVICE typename std::enable_if<!std::is_abstract<SetType>::value, PositionType>::type
  getSize(const SetType* s) const
  {
    SLIC_ASSERT_MSG(s != nullptr, "nullptr in BivariateSet::getSize()");
    return static_cast<SetType>(*s).size();
  }

protected:
  const FirstSetType* m_set1;
  const SecondSetType* m_set2;
};

/*!
 * \brief Templated helper class to select the correct base class to inherit
 *  for instances of RelationSets or ProductSets.
 *
 *  This class helps select a virtual or non-virtual interface depending on
 *  the interface type:
 *  * In the ConcreteInterface case, the base class ConcreteBivariateSet
 *    is used to avoid virtual function calls in the derived instance.
 *  * In the VirtualInterface case, the base class is BivariateSet, which
 *    contains a common virtual interface between Relation/ProductSets.
 *
 * \tparam InterfacePolicy the interface type to select (virtual or concrete)
 *
 * \see InterfacePolicies.hpp
 */
template <typename InterfacePolicy, typename Set1, typename Set2>
struct BSetInterfaceSelector;

template <typename InterfacePolicy, typename Set1, typename Set2>
struct BSetInterfaceSelector
{
  static_assert(std::is_same<InterfacePolicy, ConcreteInterface>::value,
                "InterfacePolicy must be one of policies::ConcreteInterface or "
                "policies::VirtualInterface.");
  using Type = ConcreteBivariateSet<Set1, Set2>;
};

template <typename Set1, typename Set2>
struct BSetInterfaceSelector<VirtualInterface, Set1, Set2>
{
  using Type = BivariateSet<Set1, Set2>;
};
}  // namespace detail

template <typename InterfacePolicy, typename FromSet, typename ToSet>
using BivariateSetInterface =
  typename detail::BSetInterfaceSelector<InterfacePolicy, FromSet, ToSet>::Type;

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_BivarSetIfacePolicies_HPP
