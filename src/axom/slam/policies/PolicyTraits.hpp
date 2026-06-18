// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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

#include <type_traits>
#include <utility>

namespace axom::slam::policies
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

}  // end namespace axom::slam::policies

namespace axom::slam::traits
{
///\name has_relation_ptr traits class
///@{

template <class T, class = void>
struct has_relation_ptr : std::false_type
{ };

template <class T>
struct has_relation_ptr<T, std::void_t<decltype(std::declval<T>().getRelation())>> : std::true_type
{ };

///@}

///\name indices_use_indirection traits class for BivariateSetTypes
///@{

template <class T, class = void>
struct indices_use_indirection : std::true_type
{ };

template <class T>
struct indices_use_indirection<T, std::void_t<typename T::ProductSetType>> : std::false_type
{ };

///@}

}  // end namespace axom::slam::traits

#endif  // SLAM_POLICY_TRAITS_H_
