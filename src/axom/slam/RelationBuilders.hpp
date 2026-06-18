// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file RelationBuilders.hpp
 *
 * \brief Free-function "make" helpers that construct SLAM relations while
 *  deducing the from/to set types and the policy stack.
 *
 * A static, variable-cardinality (CSR-style) relation is configured by
 * a cardinality policy, an indirection policy, and the from/to set types,
 * then built through a chained RelationBuilder over begins/indices SetBuilders:
 *
 * \code
 *   using Rel = slam::StaticRelation<P, E,
 *                 policies::VariableCardinality<P, STLIndirection>,
 *                 STLIndirection, FromSet, ToSet>;
 *   Rel r(Rel::RelationBuilder()
 *           .fromSet(&from).toSet(&to)
 *           .begins (Rel::RelationBuilder::BeginsSetBuilder ().size(off.size()).data(&off))
 *           .indices(Rel::RelationBuilder::IndicesSetBuilder().size(idx.size()).data(&idx)));
 * \endcode
 *
 * make_variable_relation collapses that to one call, deducing FromSet and ToSet from the set pointers:
 *
 * \code
 *   auto r = slam::make_variable_relation(&from, &to, offsets, indices);
 * \endcode
 *
 * \note See SetBuilders.hpp for why these are free functions rather than class-template-argument
 *  deduction guides (the builder argument is a non-deduced nested-name context).
 */

#ifndef SLAM_RELATION_BUILDERS_H_
#define SLAM_RELATION_BUILDERS_H_

#include "axom/slam/StaticRelation.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"

#include "axom/core/ArrayView.hpp"

#include <vector>

namespace axom::slam
{
/// \name Relation construction helpers
/// \{

/*!
 * \brief Make a static, variable-cardinality (CSR) relation 
 *  from \a fromSet to \a toSet, backed by std::vector storage for its begins and indices.
 *
 * The from/to set types are deduced from the pointers. 
   The begins offsets and flat indices are taken from \a begins and \a indices
 * (which must outlive the relation). The relation uses STL-vector indirection.
 *
 * \param fromSet pointer to the from-set (must outlive the relation)
 * \param toSet   pointer to the to-set (must outlive the relation)
 * \param begins  the per-from-element begin offsets (size == fromSet->size()+1)
 * \param indices the flat to-set indices
 * \return a StaticRelation with VariableCardinality and STLVector indirection
 *
 * \pre begins.size() == fromSet->size() + 1
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename FromSet::ElementType>
auto make_variable_relation(FromSet* fromSet,
                            ToSet* toSet,
                            std::vector<PosType>& begins,
                            std::vector<ElemType>& indices)
{
  using BeginsIndirection = policies::STLVectorIndirection<PosType, PosType>;
  using IndicesIndirection = policies::STLVectorIndirection<PosType, ElemType>;
  using Cardinality = policies::VariableCardinality<PosType, BeginsIndirection>;
  using RelationType =
    StaticRelation<PosType, ElemType, Cardinality, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(
        typename Builder::BeginsSetBuilder().size(static_cast<PosType>(begins.size())).data(&begins))
      .indices(
        typename Builder::IndicesSetBuilder().size(static_cast<PosType>(indices.size())).data(&indices)));
}

/*!
 * \brief Make a static, variable-cardinality (CSR) relation backed by C array storage.
 *
 * \param fromSet pointer to the from-set (must outlive the relation)
 * \param toSet   pointer to the to-set (must outlive the relation)
 * \param begins  pointer to begin offsets (size == fromSet->size()+1; must outlive the relation)
 * \param beginsSize number of begin offsets
 * \param indices pointer to flat indices (must outlive the relation)
 * \param indicesSize number of indices
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename FromSet::ElementType>
auto make_variable_relation(FromSet* fromSet,
                            ToSet* toSet,
                            PosType* begins,
                            PosType beginsSize,
                            ElemType* indices,
                            PosType indicesSize)
{
  using BeginsIndirection = policies::CArrayIndirection<PosType, PosType>;
  using IndicesIndirection = policies::CArrayIndirection<PosType, ElemType>;
  using Cardinality = policies::VariableCardinality<PosType, BeginsIndirection>;
  using RelationType =
    StaticRelation<PosType, ElemType, Cardinality, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(typename Builder::BeginsSetBuilder().size(beginsSize).data(begins))
      .indices(typename Builder::IndicesSetBuilder().size(indicesSize).data(indices)));
}

/*!
 * \brief Make a static, variable-cardinality (CSR) relation backed by ArrayView storage.
 *
 * \param fromSet pointer to the from-set (must outlive the relation)
 * \param toSet   pointer to the to-set (must outlive the relation)
 * \param begins  array view of begin offsets (size == fromSet->size()+1)
 * \param indices array view of flat indices
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename FromSet::ElementType>
auto make_variable_relation(FromSet* fromSet,
                            ToSet* toSet,
                            axom::ArrayView<PosType> begins,
                            axom::ArrayView<ElemType> indices)
{
  using BeginsIndirection = policies::ArrayViewIndirection<PosType, PosType>;
  using IndicesIndirection = policies::ArrayViewIndirection<PosType, ElemType>;
  using Cardinality = policies::VariableCardinality<PosType, BeginsIndirection>;
  using RelationType =
    StaticRelation<PosType, ElemType, Cardinality, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  return RelationType(
    Builder()
      .fromSet(fromSet)
      .toSet(toSet)
      .begins(
        typename Builder::BeginsSetBuilder().size(static_cast<PosType>(begins.size())).data(begins))
      .indices(
        typename Builder::IndicesSetBuilder().size(static_cast<PosType>(indices.size())).data(indices)));
}

/*!
 * \brief Make a static, constant-cardinality relation with a runtime stride, backed by std::vector indices.
 *
 * \param fromSet pointer to the from-set (must outlive the relation)
 * \param toSet   pointer to the to-set (must outlive the relation)
 * \param stride  number of to-set elements per from-set element
 * \param indices flat indices (size == fromSet->size() * stride; must outlive the relation)
 */
template <typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename FromSet::ElementType>
auto make_constant_relation(FromSet* fromSet, ToSet* toSet, PosType stride, std::vector<ElemType>& indices)
{
  using IndicesIndirection = policies::STLVectorIndirection<PosType, ElemType>;
  using CTy = policies::ConstantCardinality<PosType, policies::RuntimeStride<PosType>>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  auto begins_builder = typename Builder::BeginsSetBuilder().stride(stride);
  return RelationType(Builder()
                        .fromSet(fromSet)
                        .toSet(toSet)
                        .begins(begins_builder)
                        .indices(typename Builder::IndicesSetBuilder()
                                   .size(static_cast<PosType>(indices.size()))
                                   .data(&indices)));
}

/*!
 * \brief Make a static, constant-cardinality relation with a compile-time stride, backed by C array indices.
 *
 * \tparam STRIDE number of to-set elements per from-set element
 */
template <int STRIDE,
          typename FromSet,
          typename ToSet,
          typename PosType = typename FromSet::PositionType,
          typename ElemType = typename FromSet::ElementType>
auto make_constant_relation_ct(FromSet* fromSet, ToSet* toSet, ElemType* indices, PosType indicesSize)
{
  using IndicesIndirection = policies::CArrayIndirection<PosType, ElemType>;
  using StridePolicy = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using CTy = policies::ConstantCardinality<PosType, StridePolicy>;
  using RelationType = StaticRelation<PosType, ElemType, CTy, IndicesIndirection, FromSet, ToSet>;
  using Builder = typename RelationType::RelationBuilder;

  return RelationType(Builder().fromSet(fromSet).toSet(toSet).indices(
    typename Builder::IndicesSetBuilder().size(indicesSize).data(indices)));
}

/// \}

}  // end namespace axom::slam

#endif  // SLAM_RELATION_BUILDERS_H_
