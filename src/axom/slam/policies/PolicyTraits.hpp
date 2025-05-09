// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_POLICY_TRAITS_H_
#define SLAM_POLICY_TRAITS_H_

/**
 * \file PolicyTraits.hpp
 *
 * A collection of utility policy and traits classes for Slam.
 */

#include "axom/slam/Set.hpp"
#include "axom/slam/NullSet.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

#include <utility>

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \brief Definition of a type trait to adapt a StridePolicy into a SizePolicy
 */
template <typename StridePolicyType, typename IntType, int VAL = 1>
struct StrideToSize
{
  using SizeType = CompileTimeSize<IntType, VAL>;
};

/**
 * \brief Specialization of StrideToSize trait for a RuntimeStride
 */
template <typename IntType>
struct StrideToSize<RuntimeStride<IntType>, IntType>
{
  using SizeType = RuntimeSize<IntType>;
};

/**
 * \brief Specialization of StrideToSize trait for a CompileTimeStride
 */
template <typename IntType, int VAL>
struct StrideToSize<CompileTimeStride<IntType, IntType(VAL)>, IntType, VAL>
{
  using SizeType = CompileTimeSize<IntType, IntType(VAL)>;
};

/**
 * \brief Type traits for null sets.
 *
 * The null pointer for most sets is nullptr
 */
template <typename SetType, typename P = typename SetType::PositionType, typename E = typename SetType::ElementType>
struct EmptySetTraits
{
  using EmptySetType = SetType;

  AXOM_HOST_DEVICE static EmptySetType* emptySet() { return nullptr; }

  AXOM_HOST_DEVICE static bool isEmpty(const EmptySetType* set)
  {
    return (set == emptySet() || set->empty());
  }
};

/**
 * \brief Specialization of NullSetTraits for the base class Set
 *
 * The null pointer is of type NullSet
 */
template <typename P, typename E>
struct EmptySetTraits<slam::Set<P, E>>
{
  using EmptySetType = slam::Set<P, E>;

  AXOM_HOST_DEVICE static EmptySetType* emptySet()
  {
#ifndef AXOM_DEVICE_CODE
    static slam::NullSet<P, E> s_nullSet {};
    return &s_nullSet;
#else
    return nullptr;
#endif
  }

  AXOM_HOST_DEVICE static bool isEmpty(const EmptySetType* set)
  {
    return set == nullptr || set->empty();
  }
};

}  // end namespace policies

namespace traits
{
// Implementation of void_t (from C++17) with bug fix for
// earlier versions of gcc. Credit: https://stackoverflow.com/a/35754473
namespace void_details
{
template <class...>
struct make_void
{
  using type = void;
};
}  // namespace void_details

template <class... T>
using void_t = typename void_details::make_void<T...>::type;

///\name has_relation_ptr traits class
///@{

template <class T, class = void>
struct has_relation_ptr : std::false_type
{ };

template <class T>
struct has_relation_ptr<T, void_t<decltype(std::declval<T>().getRelation())>> : std::true_type
{ };

///@}

///\name indices_use_indirection traits class for BivariateSetTypes
///@{

template <class T, class = void>
struct indices_use_indirection : std::true_type
{ };

template <class T>
struct indices_use_indirection<T, void_t<typename T::ProductSetType>> : std::false_type
{ };

///@}

}  // end namespace traits

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICY_TRAITS_H_
