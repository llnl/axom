// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file CardinalityPolicies.hpp
 *
 * \brief Cardinality policies for Slam
 *
 * Cardinality policies are meant to represent the cardinality of a relation
 * with respect to an element of a an OrderedSet, i.e., the number of elements
 * of a FromSet to which each element of a ToSet maps.
 *
 * This file implements two concrete cardinality policies:
 * * ConstantCardinality, in which every member of the FromSet maps to a fixed
 *   number of entries in the ToSet
 * * VariableCardinality, in which members of the FromSet map to an arbitrary
 *   number of entries in the ToSet
 *
 * A valid cardinality policy must support the following interface:
 *  * RelationalOperatorSizeType
 *    -- A public type that indicates the SizePolicy for the each entry in
 *        the Cardinality relation \see SizePolicies.hpp
 *  * size(ElementType idx) : const ElementType
 *    -- returns the cardinality of the relation for element with index
 *        idx of the FromSet
 *  * offset(ElementType idx) : const ElementType
 *     -- returns the offset to the first element of the ToSet for element
 *        with index idx of the FromSet
 *  * firstIndex(ElementType offset) : const ElementType
 *     -- returns the element index in the FromSet given an offset into the
 *        relation
 *  * totalSize(): int
 *     -- returns the total number of elements in this relation.
 *        That is, the sum of size(idx) for each element (with index idx)
 *        of the from set.
 *  * isValid(): bool
 *     -- indicates whether the CardinalityPolicy instance is valid
 *
 */

#ifndef SLAM_POLICIES_CARDINALITY_H_
#define SLAM_POLICIES_CARDINALITY_H_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include "axom/slam/OrderedSet.hpp"  // Note: Not a circular dependency since
// CardinalityPolicies are for relations

namespace axom
{
namespace slam
{
namespace policies
{
/*!
 * \class ConstantCardinality
 * \brief Represents a mapping between two sets, where each element in the
 *  first set maps to a fixed number of elements in the second set
 *
 * \tparam ElementType the index data type
 * \tparam StridePolicy policy for number of elements being mapped
 */
template <typename ElementType = int, typename StridePolicy = RuntimeStride<ElementType>>
struct ConstantCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StridePolicy;
  using BeginsIndirectionPolicy = NoIndirection<ElementType, ElementType>;

  // runtime size (fromSet.size()), striding from template parameter, no offset
  using BeginsSet =
    OrderedSet<ElementType, ElementType, BeginsSizePolicy, BeginsOffsetPolicy, BeginsStridePolicy>;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType =
    typename StrideToSize<BeginsStridePolicy, ElementType, BeginsStridePolicy::DEFAULT_VALUE>::SizeType;

  using IndirectionPtrType = typename BeginsIndirectionPolicy::IndirectionPtrType;

  ConstantCardinality() : m_begins() { }
  ConstantCardinality(BeginsSet begins) : m_begins(begins) { }
  ConstantCardinality(ElementType fromSetSize) { m_begins = BeginsSet(fromSetSize); }

  ConstantCardinality(ElementType fromSetSize, typename BeginsSet::SetBuilder& builder)
  {
    // needs a size and a stride (when runtime)
    builder.size(fromSetSize);
    m_begins = builder;
  }

  AXOM_HOST_DEVICE ElementType size(ElementType AXOM_UNUSED_PARAM(fromPos)) const
  {
    return m_begins.stride();
  }

  AXOM_HOST_DEVICE ElementType offset(ElementType fromPos) const { return m_begins[fromPos]; }

  AXOM_HOST_DEVICE ElementType firstIndex(ElementType offset) const
  {
    return offset / m_begins.stride();
  }

  IndirectionPtrType offsetData() { return m_begins.ptr(); }

  const IndirectionPtrType offsetData() const { return m_begins.ptr(); }

  void bindBeginOffsets(ElementType fromSetSize, ElementType stride)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize).stride(stride);
  }

  ElementType totalSize() const { return m_begins.stride() * m_begins.size(); }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool AXOM_UNUSED_PARAM(verboseOutput) = false) const
  {
    return m_begins.size() == fromSet->size();
  }

  BeginsSet m_begins;
};

/*!
 * \class VariableCardinality
 * \brief Represents a mapping between two sets, where each element in the
 *  first set maps to an arbitrary number of elements in the second set.
 *
 * \tparam ElementType the index data type
 * \tparam IndirectionPolicy the policy to use for storing offsets and indices
 */
template <typename ElementType = int,
          typename IndirectionPolicy = STLVectorIndirection<ElementType, ElementType>>
struct VariableCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StrideOne<ElementType>;
  using BeginsIndirectionPolicy = IndirectionPolicy;

  // runtime size (fromSet.size()), striding from template parameter, no offset
  using BeginsSet =
    OrderedSet<ElementType, ElementType, BeginsSizePolicy, BeginsOffsetPolicy, BeginsStridePolicy, IndirectionPolicy>;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType = BeginsSizePolicy;

  using IndirectionBufferType = typename IndirectionPolicy::IndirectionBufferType;
  using IndirectionPtrType = typename IndirectionPolicy::IndirectionPtrType;

  VariableCardinality() : m_begins() { }
  VariableCardinality(BeginsSet begins) : m_begins(begins) { }
  VariableCardinality(ElementType fromSetSize, typename BeginsSet::SetBuilder& builder)
  {
    builder.size(fromSetSize + 1);
    m_begins = builder;
  }

  void bindBeginOffsets(ElementType fromSetSize, IndirectionPtrType data)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize + 1).data(data);
  }

  AXOM_HOST_DEVICE ElementType size(ElementType fromPos) const
  {
    return offset(fromPos + 1) - offset(fromPos);
  }

  AXOM_HOST_DEVICE ElementType offset(ElementType fromPos) const { return m_begins[fromPos]; }

  AXOM_HOST_DEVICE ElementType firstIndex(ElementType relationOffset) const
  {
    for(ElementType firstIdx = 0; firstIdx < m_begins.size() - 1; firstIdx++)
    {
      if(offset(firstIdx + 1) > relationOffset)
      {
        return firstIdx;
      }
    }
    return -1;
  }

  IndirectionPtrType offsetData() { return m_begins.data(); }

  const IndirectionPtrType offsetData() const { return m_begins.data(); }

  ElementType totalSize() const
  {
    return m_begins.empty() ? ElementType() : offset(m_begins.size() - 1);
  }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool verboseOutput = false) const
  {
    return m_begins.size() == (fromSet->size() + 1) &&
      static_cast<IndirectionPolicy>(m_begins).isValid(m_begins.size(),
                                                       m_begins.offset(),
                                                       m_begins.stride(),
                                                       verboseOutput);
  }

  BeginsSet m_begins;
};

/*!
 * \class MappedVariableCardinality
 * \brief Represents a mapping between two sets, where each element in the
 *  first set maps to an arbitrary number of elements in the second set.
 *
 *  MappedVariableCardinality extends VariableCardinality to map "flat" indices
 *  in the associated RelationSet to first set indices.
 *
 * \tparam ElementType the index data type
 * \tparam IndirectionPolicy the policy to use for storing offsets and indices
 */
template <typename ElementType = int, typename IndirectionPolicy = ArrayIndirection<ElementType, ElementType>>
struct MappedVariableCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StrideOne<ElementType>;
  using BeginsIndirectionPolicy = IndirectionPolicy;

  // runtime size (fromSet.size()), striding from template parameter, no offset
  using IndexSet =
    OrderedSet<ElementType, ElementType, BeginsSizePolicy, BeginsOffsetPolicy, BeginsStridePolicy, IndirectionPolicy>;
  using BeginsSet = IndexSet;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType = BeginsSizePolicy;

  using IndirectionBufferType = typename IndirectionPolicy::IndirectionBufferType;
  using IndirectionPtrType = typename IndirectionPolicy::IndirectionPtrType;

  MappedVariableCardinality() : m_begins() { }
  MappedVariableCardinality(BeginsSet begins) : m_begins(begins) { }
  MappedVariableCardinality(ElementType fromSetSize, typename BeginsSet::SetBuilder& builder)
  {
    builder.size(fromSetSize + 1);
    m_begins = builder;
  }

  void bindBeginOffsets(ElementType fromSetSize, IndirectionPtrType data)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize + 1).data(data);
  }

  void bindFirstIndices(ElementType relationSize, IndirectionPtrType data, bool fillIndices = true)
  {
    m_firstIndexes = typename IndexSet::SetBuilder().size(relationSize + 1).data(data);
    if(fillIndices)
    {
      // Construct the flat-to-first mapping.
      for(int fromIdx = 0; fromIdx < m_begins.size() - 1; fromIdx++)
      {
        int beginIdx = offset(fromIdx);
        for(int slotIdx = 0; slotIdx < size(fromIdx); slotIdx++)
        {
          m_firstIndexes[slotIdx + beginIdx] = fromIdx;
        }
      }
    }
  }

  AXOM_HOST_DEVICE ElementType size(ElementType fromPos) const
  {
    return offset(fromPos + 1) - offset(fromPos);
  }

  AXOM_HOST_DEVICE ElementType offset(ElementType fromPos) const { return m_begins[fromPos]; }

  AXOM_HOST_DEVICE ElementType firstIndex(ElementType offset) const
  {
    return m_firstIndexes[offset];
  }

  IndirectionPtrType offsetData() { return m_begins.data(); }

  const IndirectionPtrType offsetData() const { return m_begins.data(); }

  ElementType totalSize() const
  {
    return m_begins.empty() ? ElementType() : offset(m_begins.size() - 1);
  }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool verboseOutput = false) const
  {
    return m_begins.size() == (fromSet->size() + 1) &&
      static_cast<IndirectionPolicy>(m_begins).isValid(m_begins.size(),
                                                       m_begins.offset(),
                                                       m_begins.stride(),
                                                       verboseOutput);
  }

  IndexSet m_firstIndexes;
  BeginsSet m_begins;
};

}  // end namespace policies

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_CARDINALITY_H_
