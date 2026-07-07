// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Aliases.hpp
 *
 * \brief Convenience aliases for Slam's most common set and relation configurations.
 *
 * Slam's sets, relations and maps are assembled from policy template parameters
 * (cardinality, stride, indirection, offset, size, subsetting, interface).
 * This lets users customize data structure to their required data model,
 * avoiding unnecessary computational and memory overhead.
 * See the Slam user guide for how to choose policies.
 *
 * A handful of \e relations are commonly used, warranting named shorthands.
 * The aliases in this file for the  common relation types, as well as some set types
 * use `axom::Array` for storage the object manages itself and `axom::ArrayView` for
 * storage managed elsewhere:
 *
 *   - \c RangeSet<P,E>                -- a contiguous range of positions (from RangeSet.hpp)
 *   - \c ArraySet<P,E>                -- a set indexed through an \c axom::Array it manages
 *   - \c ArrayViewSet<P,E>            -- a set indexed through an \c axom::ArrayView managed elsewhere
 *   - \c VariableRelation<F,T>        -- a static relation with variable cardinality, managing its buffers
 *   - \c VariableRelationView<F,T>    -- as above, viewing buffers managed elsewhere
 *   - \c ConstantRelation<F,T,N>      -- a static relation with compile-time cardinality N, managing its buffer
 *   - \c ConstantRelationView<F,T,N>  -- as above, viewing a buffer managed elsewhere
 *   - \c RuntimeConstantRelation<F,T> -- a static relation with runtime-constant cardinality, managing its buffer
 *   - \c RuntimeConstantRelationView<F,T> -- as above, viewing a buffer managed elsewhere
 *
 * These aliases are a convenience for the commonly used configurations.
 * Use the underlying policies directly whenever you need other configuration,
 * e.g. -- a specific stride, offset, size, subsetting or interface, a specialized cardinality such as
 * \c policies::MappedVariableCardinality, or a different buffer such as \c std::vector
 * (\c policies::STLVectorIndirection) or a raw pointer (\c policies::CArrayIndirection).
 *
 * \note This file does not currently alias map types.  A \c Map or \c BivariateMap already
 *  defaults to the common case of `axom::Array` indirection with stride one,
 *  while cases that do vary the map (a view into an externally-managed buffer,
 *  or a runtime stride) also tend to vary in ways a small fixed set of names cannot capture cleanly.
 *
 * \note "Manages its buffer" vs. "views a buffer" refers to lifetime, not to whether
 *  the data is logically owned by the object. An \c axom::Array indirection allocates
 *  and frees the buffer as part of the Slam object's lifetime, while an \c axom::ArrayView
 *  indirection refers to a buffer whose lifetime is managed elsewhere and must outlive the Slam object.
 */

#ifndef SLAM_ALIASES_H_
#define SLAM_ALIASES_H_

#include "axom/slam/RangeSet.hpp"
#include "axom/slam/IndirectionSet.hpp"
#include "axom/slam/StaticRelation.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/BivariateMap.hpp"

#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

namespace axom
{
namespace slam
{
/*!
 * \brief A set whose \a ElemType elements are read from an \c axom::Array.
 *
 * The \c axom::Array is managed outside the set and must outlive it.
 * \sa ArrayViewSet for the device-capturable form
 */
template <typename PosType = DefaultPositionType, typename ElemType = DefaultElementType>
using ArraySet = ArrayIndirectionSet<PosType, ElemType>;

/*!
 * \brief A set whose \a ElemType elements are read through an \c axom::ArrayView.
 *
 * The view is held by value and refers to a buffer managed elsewhere (which must outlive the set);
 * because \c axom::ArrayView is trivially copyable, it can be used for device capture.
 */
template <typename PosType = DefaultPositionType, typename ElemType = DefaultElementType>
using ArrayViewSet = ArrayViewIndirectionSet<PosType, ElemType>;

/*!
 * \brief A static relation from \a FromSet to \a ToSet with per-element (variable) cardinality,
 *  reading offsets and indices from \c axom::Array buffers managed elsewhere.
 *
 * This is the type produced by the \c axom::Array overload of \c make_variable_relation.
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename ToSet::PositionType>
using VariableRelation =
  StaticRelation<PosType,
                 ElemType,
                 policies::VariableCardinality<PosType, policies::ArrayIndirection<PosType, PosType>>,
                 policies::ArrayIndirection<PosType, ElemType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation view from \a FromSet to \a ToSet with per-element
 *  (variable) cardinality, binding offsets and indices through \c axom::ArrayView.
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename ToSet::PositionType>
using VariableRelationView =
  StaticRelation<PosType,
                 ElemType,
                 policies::VariableCardinality<PosType, policies::ArrayViewIndirection<PosType, PosType>>,
                 policies::ArrayViewIndirection<PosType, ElemType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation from \a FromSet to \a ToSet with fixed cardinality \a N
 *  (each from-element maps to exactly N to-elements), reading indices from an \c axom::Array buffer managed elsewhere.
 */
template <typename FromSet,
          typename ToSet,
          int N,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename ToSet::PositionType>
using ConstantRelation =
  StaticRelation<PosType,
                 ElemType,
                 policies::ConstantCardinality<PosType, policies::CompileTimeStride<PosType, N>>,
                 policies::ArrayIndirection<PosType, ElemType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation view from \a FromSet to \a ToSet with fixed cardinality \a N,
 *  binding indices through \c axom::ArrayView.
 */
template <typename FromSet,
          typename ToSet,
          int N,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename ToSet::PositionType>
using ConstantRelationView =
  StaticRelation<PosType,
                 ElemType,
                 policies::ConstantCardinality<PosType, policies::CompileTimeStride<PosType, N>>,
                 policies::ArrayViewIndirection<PosType, ElemType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation from \a FromSet to \a ToSet with runtime constant cardinality,
 *  reading indices from an \c axom::Array buffer managed elsewhere.
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename ToSet::PositionType>
using RuntimeConstantRelation =
  StaticRelation<PosType,
                 ElemType,
                 policies::ConstantCardinality<PosType, policies::RuntimeStride<PosType>>,
                 policies::ArrayIndirection<PosType, ElemType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation view from \a FromSet to \a ToSet with runtime constant cardinality,
 *  binding indices through \c axom::ArrayView.
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename ToSet::PositionType>
using RuntimeConstantRelationView =
  StaticRelation<PosType,
                 ElemType,
                 policies::ConstantCardinality<PosType, policies::RuntimeStride<PosType>>,
                 policies::ArrayViewIndirection<PosType, ElemType>,
                 FromSet,
                 ToSet>;

}  // namespace slam
}  // namespace axom

#endif  // SLAM_ALIASES_H_
