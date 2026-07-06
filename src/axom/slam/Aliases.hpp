// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Aliases.hpp
 *
 * \brief Blessed, user-facing aliases for Slam's common set/relation/map types.
 *
 * Slam's containers are assembled from policy template parameters, which can be verbose.
 * This header provides aliases for common configurations using Axom-native storage:
 *
 *   - \c RangeSet<P,E>                -- a contiguous range of positions (from RangeSet.hpp)
 *   - \c ArraySet<P,E>                -- a set that binds an external \c axom::Array
 *   - \c ArrayViewSet<P,E>            -- a set that binds an \c axom::ArrayView by value
 *   - \c VariableRelation<F,T>        -- an Array-backed static relation with variable cardinality
 *   - \c VariableRelationView<F,T>    -- an ArrayView-backed static relation with variable cardinality
 *   - \c ConstantRelation<F,T,N>      -- an Array-backed static relation with compile-time cardinality N
 *   - \c ConstantRelationView<F,T,N>  -- an ArrayView-backed static relation with compile-time cardinality N
 *   - \c RuntimeConstantRelation<F,T> -- an Array-backed static relation with runtime constant cardinality
 *   - \c RuntimeConstantRelationView<F,T>
 *       -- an ArrayView-backed static relation with runtime constant cardinality
 *   - \c ArrayMap<S,T,STRIDE>         -- an Array-backed map from set S to T
 *   - \c ArrayViewMap<S,T,STRIDE>     -- an ArrayView-backed map from set S to T
 *   - \c RuntimeArrayMap<S,T>         -- an Array-backed map with runtime stride
 *   - \c RuntimeArrayViewMap<S,T>     -- an ArrayView-backed map with runtime stride
 *   - \c ArrayBivariateMap<BSet,T>     -- an Array-backed bivariate map
 *   - \c ArrayViewBivariateMap<BSet,T> -- an ArrayView-backed bivariate map
 *   - \c RuntimeArrayBivariateMap<BSet,T>
 *       -- an Array-backed bivariate map with runtime stride
 *   - \c RuntimeArrayViewBivariateMap<BSet,T>
 *       -- an ArrayView-backed bivariate map with runtime stride
 *
 * Prefer these aliases in user-facing and component-facing Axom code.
 * They keep the set/relation/map concepts visible while hiding the common policy stack.
 * Use the full policies when adapting legacy storage, such as
 * \c std::vector via \c policies::STLVectorIndirection, 
 * raw pointers via \c policies::CArrayIndirection, 
 * custom third-party buffers, or less common combinations of 
 * stride, offset, size, subsetting, indirection or interface policies.
 *
 * For maps, the ownership distinction follows Axom's container vocabulary:
 * \c ArrayMap and \c ArrayBivariateMap own an \c axom::Array buffer, 
 * while \c ArrayViewMap and \c ArrayViewBivariateMap store an \c axom::ArrayView by value
 * and view storage owned elsewhere. Owning maps are convenient host-side storage objects
 * while view maps are the preferred lightweight form to pass into kernels and generic algorithms.
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
 * \brief A set whose \a ElemType elements are stored in an external \c axom::Array.
 *
 * This is the canonical Axom indirection set for binding an \c axom::Array buffer.
 * The array object is managed outside the set and must outlive it.
 * \sa ArrayViewSet for the non-owning, device-capturable form
 */
template <typename PosType = DefaultPositionType, typename ElemType = DefaultElementType>
using ArraySet = ArrayIndirectionSet<PosType, ElemType>;

/*!
 * \brief A set whose \a ElemType elements are viewed through an \c axom::ArrayView.
 *
 * The view is stored by value, so this is the preferred indirection-set shape
 * for device-capturable views over storage owned elsewhere.
 */
template <typename PosType = DefaultPositionType, typename ElemType = DefaultElementType>
using ArrayViewSet = ArrayViewIndirectionSet<PosType, ElemType>;

/*!
 * \brief A static relation from \a FromSet to \a ToSet with per-element (variable) cardinality,
 *  binding offsets and indices in external \c axom::Array buffers.
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
 *  (each from-element maps to exactly N to-elements), binding indices in an external \c axom::Array buffer.
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
 *  binding indices in an external \c axom::Array buffer.
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

/*!
 * \brief A map from the set \a S to values of type \a T, with compile-time \a STRIDE
 *  components per element (default 1), owning its values in an \c axom::Array buffer.
 */
template <typename S, typename T = DefaultElementType, int STRIDE = 1>
using ArrayMap = Map<T,
                     S,
                     policies::ArrayIndirection<typename S::PositionType, T>,
                     policies::CompileTimeStride<typename S::PositionType, STRIDE>>;

/*!
 * \brief A non-owning map view from the set \a S to values of type \a T, with compile-time
 *  \a STRIDE components per element (default 1).
 */
template <typename S, typename T = DefaultElementType, int STRIDE = 1>
using ArrayViewMap = Map<T,
                         S,
                         policies::ArrayViewIndirection<typename S::PositionType, T>,
                         policies::CompileTimeStride<typename S::PositionType, STRIDE>>;

/*!
 * \brief A map from the set \a S to values of type \a T with runtime stride,
 *  owning its values in an \c axom::Array buffer.
 */
template <typename S, typename T = DefaultElementType>
using RuntimeArrayMap = Map<T,
                            S,
                            policies::ArrayIndirection<typename S::PositionType, T>,
                            policies::RuntimeStride<typename S::PositionType>>;

/*!
 * \brief A non-owning map view from the set \a S to values of type \a T with runtime stride.
 */
template <typename S, typename T = DefaultElementType>
using RuntimeArrayViewMap = Map<T,
                                S,
                                policies::ArrayViewIndirection<typename S::PositionType, T>,
                                policies::RuntimeStride<typename S::PositionType>>;

/*!
 * \brief A bivariate map of \a T over the bivariate set \a BSet
 *  (e.g. a \c ProductSet or \c RelationSet), owning its values in an \c axom::Array buffer.
 */
template <typename BSet, typename T = DefaultElementType, int STRIDE = 1>
using ArrayBivariateMap =
  BivariateMap<T,
               BSet,
               policies::ArrayIndirection<typename BSet::PositionType, T>,
               policies::CompileTimeStride<typename BSet::PositionType, STRIDE>>;

/*!
 * \brief A non-owning bivariate map view of \a T over the bivariate set \a BSet.
 */
template <typename BSet, typename T = DefaultElementType, int STRIDE = 1>
using ArrayViewBivariateMap =
  BivariateMap<T,
               BSet,
               policies::ArrayViewIndirection<typename BSet::PositionType, T>,
               policies::CompileTimeStride<typename BSet::PositionType, STRIDE>>;

/*!
 * \brief A bivariate map of \a T over the bivariate set \a BSet with runtime stride,
 *  owning its values in an \c axom::Array buffer.
 */
template <typename BSet, typename T = DefaultElementType>
using RuntimeArrayBivariateMap =
  BivariateMap<T,
               BSet,
               policies::ArrayIndirection<typename BSet::PositionType, T>,
               policies::RuntimeStride<typename BSet::PositionType>>;

/*!
 * \brief A non-owning bivariate map view of \a T over the bivariate set \a BSet with runtime stride.
 */
template <typename BSet, typename T = DefaultElementType>
using RuntimeArrayViewBivariateMap =
  BivariateMap<T,
               BSet,
               policies::ArrayViewIndirection<typename BSet::PositionType, T>,
               policies::RuntimeStride<typename BSet::PositionType>>;

}  // namespace slam
}  // namespace axom

#endif  // SLAM_ALIASES_H_
