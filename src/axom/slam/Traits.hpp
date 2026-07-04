// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Traits.hpp
 *
 * \brief Named, compile-time trait predicates for slam's core concepts.
 *
 * These variable templates give slam's implicit, duck-typed policy and container contracts
 * a name, documentation, and enforcement within C++17, and will be converted to C++20 concepts.
 *
 * Everything in this header is compile-time, so the predicates are usable in
 * host and kernel contexts.
 * 
 * The predicates categorize a type by the interface it exposes, using the
 * distinguishing member typedefs each slam concept carries:
 *   - a Set exposes \c PositionType and \c ElementType and a zero-argument \c size();
 *   - a Relation exposes \c FromSetType and \c ToSetType;
 *   - a BivariateSet exposes \c FirstSetType and \c SecondSetType;
 *   - a Map exposes \c SetType.
 * Because a BivariateSet also carries \c PositionType, \c ElementType and \c size(),
 * \c is_set_like_v excludes bivariate sets explicitly.
 */

#ifndef SLAM_TRAITS_H_
#define SLAM_TRAITS_H_

#include <type_traits>
#include <utility>

namespace axom::slam
{
namespace detail
{

//---------------------------------------------------------------------------
// Member-typedef detectors (void_t idiom).
//---------------------------------------------------------------------------

template <typename T, typename = void>
struct has_position_element : std::false_type
{ };
template <typename T>
struct has_position_element<T, std::void_t<typename T::PositionType, typename T::ElementType>>
  : std::true_type
{ };

template <typename T, typename = void>
struct has_from_to_set : std::false_type
{ };
template <typename T>
struct has_from_to_set<T, std::void_t<typename T::FromSetType, typename T::ToSetType>> : std::true_type
{ };

template <typename T, typename = void>
struct has_first_second_set : std::false_type
{ };
template <typename T>
struct has_first_second_set<T, std::void_t<typename T::FirstSetType, typename T::SecondSetType>>
  : std::true_type
{ };

template <typename T, typename = void>
struct has_set_type : std::false_type
{ };
template <typename T>
struct has_set_type<T, std::void_t<typename T::SetType>> : std::true_type
{ };

template <typename T, typename = void>
struct has_indirection_result : std::false_type
{ };
template <typename T>
struct has_indirection_result<T, std::void_t<typename T::IndirectionResult>> : std::true_type
{ };

//---------------------------------------------------------------------------
// Member-function detectors (const-callable).
//---------------------------------------------------------------------------

template <typename T, typename = void>
struct has_size : std::false_type
{ };
template <typename T>
struct has_size<T, std::void_t<decltype(std::declval<const T&>().size())>> : std::true_type
{ };

template <typename T, typename = void>
struct has_begin_end : std::false_type
{ };
template <typename T>
struct has_begin_end<
  T,
  std::void_t<decltype(std::declval<const T&>().begin()), decltype(std::declval<const T&>().end())>>
  : std::true_type
{ };

// The value policies expose a zero argument value(); sets and relations do not.
template <typename T, typename = void>
struct has_value : std::false_type
{ };
template <typename T>
struct has_value<T, std::void_t<decltype(std::declval<const T&>().value())>> : std::true_type
{ };

template <typename T, typename = void>
struct has_stride_accessor : std::false_type
{ };
template <typename T>
struct has_stride_accessor<T, std::void_t<decltype(std::declval<const T&>().stride())>>
  : std::true_type
{ };

template <typename T, typename = void>
struct has_offset_accessor : std::false_type
{ };
template <typename T>
struct has_offset_accessor<T, std::void_t<decltype(std::declval<const T&>().offset())>>
  : std::true_type
{ };

// Checks that M is a Map whose SetType is S
template <typename M, typename S, bool = has_set_type<M>::value>
struct map_over : std::false_type
{ };
template <typename M, typename S>
struct map_over<M, S, true> : std::is_same<typename M::SetType, S>
{ };

}  // namespace detail

//---------------------------------------------------------------------------
// Structural concepts
//---------------------------------------------------------------------------

/*!
 * \brief True if \a T models a bivariate set (i.e. has \c FirstSetType and \c SecondSetType).
 *  Examples: \c ProductSet and \c RelationSet.
 */
template <typename T>
inline constexpr bool is_bivariate_set_like_v = detail::has_first_second_set<T>::value;

/*!
 * \brief True if \a T models a (univariate) set (i.e. has \c PositionType and \c ElementType typedefs
 *  and a zero argument \c size(). Bivariate sets are excluded even though they satisfy those members.
 */
template <typename T>
inline constexpr bool is_set_like_v = detail::has_position_element<T>::value &&
  detail::has_size<T>::value && !is_bivariate_set_like_v<T>;

/*!
 * \brief True if \a T is a set that also supports ordered iteration via
 *  const \c begin() / \c end() (e.g. \c OrderedSet, \c RangeSet).
 */
template <typename T>
inline constexpr bool is_ordered_set_like_v = is_set_like_v<T> && detail::has_begin_end<T>::value;

/*!
 * \brief True if \a T models a relation between two sets
 *  (i.e. is has   \c FromSetType and \c ToSetType )
 *  Examples: \c StaticRelation and the \c VariableRelation / \c ConstantRelation aliases.
 */
template <typename T>
inline constexpr bool is_relation_like_v = detail::has_from_to_set<T>::value;

/*!
 * \brief True if \a T models a map over a set (i.e. it has a \c SetType typedef).
 *  Examples: \c Map / \c BivariateMap.
 */
template <typename T>
inline constexpr bool is_map_like_v = detail::has_set_type<T>::value;

/*!
 * \brief True if \a M is a map whose \c SetType is \a S.
 *  SFINAE-safe for non-map \a M (evaluates to false rather than erroring).
 */
template <typename M, typename S>
inline constexpr bool is_map_over_v = detail::map_over<M, S>::value;

//---------------------------------------------------------------------------
// Policy concepts
//---------------------------------------------------------------------------

/*!
 * \brief True if \a T is a value policy
 *
 *  It exposes a zero argument \c value() accessor 
 *  ( \c RuntimeValue / \c CompileTimeValue and their leaf policies).
 *  Sets also expose \c size() but not \c value(), which keeps the two apart.
 */
template <typename T>
inline constexpr bool is_value_policy_v = detail::has_value<T>::value;

/*!
 * \brief True if \a T is a size policy
 *
 *  A value policy that also exposes the named \c size() accessor.
 */
template <typename T>
inline constexpr bool is_size_policy_v = detail::has_value<T>::value && detail::has_size<T>::value;

/*!
 * \brief True if \a T is a stride policy
 * 
 *  A value policy that also exposes the named \c stride() accessor.
 */
template <typename T>
inline constexpr bool is_stride_policy_v =
  detail::has_value<T>::value && detail::has_stride_accessor<T>::value;

/*!
 * \brief True if \a T is an offset policy
 * 
 *  A value policy that also exposes the named \c offset() accessor.
 */
template <typename T>
inline constexpr bool is_offset_policy_v =
  detail::has_value<T>::value && detail::has_offset_accessor<T>::value;

/*!
 * \brief True if \a T is an indirection policy (i.e. has an \c IndirectionResult typedef)
 * 
 *  Examples: \c NoIndirection, \c ArrayIndirection, \c ArrayViewIndirection,
 *  \c STLVectorIndirection, \c CArrayIndirection.
 */
template <typename T>
inline constexpr bool is_indirection_policy_v = detail::has_indirection_result<T>::value;

//---------------------------------------------------------------------------
// Index-space concepts for handle-like elements
//---------------------------------------------------------------------------

/// \brief True if \a T can serve as a position in a slam index space
template <typename T>
inline constexpr bool is_position_like_v = std::is_integral_v<std::remove_cv_t<T>>;

/*!
 * \brief True if \a T can serve as a slam element handle
 *
 *  It must be trivially copyable so it survives capture-by-value into device kernels.
 *  A type with a user-provided copy constructor fails this by design.
 */
template <typename T>
inline constexpr bool is_handle_like_v = std::is_trivially_copyable_v<T>;

}  // namespace axom::slam

#endif  // SLAM_TRAITS_H_
