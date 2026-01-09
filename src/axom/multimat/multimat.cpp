// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file multimat.cpp
 *
 * \brief Implementation of the MultiMat class
 */
#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/atomics.hpp"
#include "axom/core/execution/scans.hpp"

#include "axom/multimat/multimat.hpp"
#include "axom/slic.hpp"

#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>

#include <cassert>

using namespace axom::multimat;

namespace
{
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
using GPU_Exec = axom::CUDA_EXEC<256>;
#elif defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_HIP)
using GPU_Exec = axom::HIP_EXEC<256>;
#else
using GPU_Exec = axom::SEQ_EXEC;
#endif

bool AllocatorOnDevice(int allocatorId)
{
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  axom::MemorySpace space = axom::detail::getAllocatorSpace(allocatorId);
  if(space == axom::MemorySpace::Device || space == axom::MemorySpace::Unified)
  {
    return true;
  }
#endif
  AXOM_UNUSED_VAR(allocatorId);
  return false;
}

template <typename Lambda>
void ExecLambdaForMemory(int length, int allocatorId, Lambda&& function)
{
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  if(AllocatorOnDevice(allocatorId))
  {
    axom::for_all<GPU_Exec>(length, function);
    return;
  }
#else
  AXOM_UNUSED_VAR(allocatorId);
#endif
  axom::for_all<axom::SEQ_EXEC>(length, function);
}
}  // namespace

template <>
axom::Array<unsigned char>& MultiMat::FieldBacking::getArray<unsigned char>()
{
  return m_ucharData;
}
template <>
axom::Array<int>& MultiMat::FieldBacking::getArray<int>()
{
  return m_intData;
}
template <>
axom::Array<float>& MultiMat::FieldBacking::getArray<float>()
{
  return m_floatData;
}
template <>
axom::Array<double>& MultiMat::FieldBacking::getArray<double>()
{
  return m_dblData;
}

template <>
void MultiMat::FieldBacking::setArrayView<unsigned char>(axom::ArrayView<unsigned char> view)
{
  m_ucharView = view;
}
template <>
void MultiMat::FieldBacking::setArrayView<int>(axom::ArrayView<int> view)
{
  m_intView = view;
}
template <>
void MultiMat::FieldBacking::setArrayView<float>(axom::ArrayView<float> view)
{
  m_floatView = view;
}
template <>
void MultiMat::FieldBacking::setArrayView<double>(axom::ArrayView<double> view)
{
  m_dblView = view;
}

template <>
axom::ArrayView<unsigned char> MultiMat::FieldBacking::getArrayView<unsigned char>()
{
  return m_isOwned ? m_ucharData.view() : m_ucharView;
}
template <>
axom::ArrayView<int> MultiMat::FieldBacking::getArrayView<int>()
{
  return m_isOwned ? m_intData.view() : m_intView;
}
template <>
axom::ArrayView<float> MultiMat::FieldBacking::getArrayView<float>()
{
  return m_isOwned ? m_floatData.view() : m_floatView;
}
template <>
axom::ArrayView<double> MultiMat::FieldBacking::getArrayView<double>()
{
  return m_isOwned ? m_dblData.view() : m_dblView;
}

MultiMat::MultiMat(DataLayout AXOM_UNUSED_PARAM(d), SparsityLayout AXOM_UNUSED_PARAM(s))
  : m_slamAllocatorId(axom::getDefaultAllocatorID())
  , m_fieldAllocatorId(axom::getDefaultAllocatorID())
  , m_ncells(0)
  , m_nmats(0)
  , m_sets(2, 2, m_slamAllocatorId)
  , m_staticRelations(2, 2, m_slamAllocatorId)
  , m_dynamicRelations(2, 2, m_slamAllocatorId)
  , m_sparseBivarSet(2, 2, m_slamAllocatorId)
  , m_denseBivarSet(2, 2, m_slamAllocatorId)
  , m_dynamic_mode(false)
{ }

MultiMat::IndBufferType& MultiMat::relBeginVec(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? m_cellMatRel_beginsVec : m_matCellRel_beginsVec;
}

MultiMat::IndBufferType& MultiMat::relIndVec(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? m_cellMatRel_indicesVec : m_matCellRel_indicesVec;
}

MultiMat::IndBufferType& MultiMat::relFirstIndVec(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? m_cellMatRel_firstIndicesVec
                                          : m_matCellRel_firstIndicesVec;
}

MultiMat::StaticVariableRelationType& MultiMat::relStatic(DataLayout layout)
{
  return m_staticRelations[(int)layout];
}

MultiMat::DynamicVariableRelationType& MultiMat::relDynamic(DataLayout layout)
{
  return m_dynamicRelations[(int)layout];
}

MultiMat::RangeSetType& MultiMat::relDominantSet(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? getCellSet() : getMatSet();
}

MultiMat::RangeSetType& MultiMat::relSecondarySet(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? getMatSet() : getCellSet();
}

MultiMat::RelationSetType& MultiMat::relSparseSet(DataLayout layout)
{
  return m_sparseBivarSet[(int)layout];
}

const MultiMat::RelationSetType& MultiMat::relSparseSet(DataLayout layout) const
{
  return m_sparseBivarSet[(int)layout];
}

MultiMat::ProductSetType& MultiMat::relDenseSet(DataLayout layout)
{
  return m_denseBivarSet[(int)layout];
}

const MultiMat::ProductSetType& MultiMat::relDenseSet(DataLayout layout) const
{
  return m_denseBivarSet[(int)layout];
}

bool MultiMat::hasValidStaticRelation(DataLayout layout) const
{
  const StaticVariableRelationType& rel = m_staticRelations[(int)layout];

  return rel.hasFromSet() && rel.hasToSet();
}

bool MultiMat::hasValidDynamicRelation(DataLayout layout) const
{
  const DynamicVariableRelationType& rel = m_dynamicRelations[(int)layout];

  return rel.hasFromSet() && rel.hasToSet();
}

// Copy constructor
MultiMat::MultiMat(const MultiMat& other)
  : m_slamAllocatorId(other.m_slamAllocatorId)
  , m_fieldAllocatorId(other.m_fieldAllocatorId)
  , m_ncells(other.m_ncells)
  , m_nmats(other.m_nmats)
  , m_sets(other.m_sets)
  , m_cellMatRel_beginsVec(other.m_cellMatRel_beginsVec)
  , m_cellMatRel_indicesVec(other.m_cellMatRel_indicesVec)
  , m_cellMatRel_firstIndicesVec(other.m_cellMatRel_firstIndicesVec)
  , m_matCellRel_beginsVec(other.m_matCellRel_beginsVec)
  , m_matCellRel_indicesVec(other.m_matCellRel_indicesVec)
  , m_matCellRel_firstIndicesVec(other.m_matCellRel_firstIndicesVec)
  , m_staticRelations(other.m_staticRelations)
  , m_dynamicRelations(other.m_dynamicRelations)
  , m_sparseBivarSet(other.m_sparseBivarSet)
  , m_denseBivarSet(other.m_denseBivarSet)
  , m_flatCellToMatIndexMap(other.m_flatCellToMatIndexMap)
  , m_flatMatToCellIndexMap(other.m_flatMatToCellIndexMap)
  , m_fieldNameVec(other.m_fieldNameVec)
  , m_fieldMappingVec(other.m_fieldMappingVec)
  , m_dataTypeVec(other.m_dataTypeVec)
  , m_fieldDataLayoutVec(other.m_fieldDataLayoutVec)
  , m_fieldSparsityLayoutVec(other.m_fieldSparsityLayoutVec)
  , m_fieldStrideVec(other.m_fieldStrideVec)
  , m_dynamic_mode(false)
{
  if(other.hasValidStaticRelation(DataLayout::CELL_DOM))
  {
    StaticVariableRelationType& cellMatRel = relStatic(DataLayout::CELL_DOM);

    cellMatRel = StaticVariableRelationType(&getCellSet(), &getMatSet());
    cellMatRel.bindBeginOffsets(getCellSet().size(), m_cellMatRel_beginsVec.view());
    cellMatRel.bindIndices(m_cellMatRel_indicesVec.size(), m_cellMatRel_indicesVec.view());
    cellMatRel.bindFirstIndices(m_cellMatRel_firstIndicesVec.size(),
                                m_cellMatRel_firstIndicesVec.view(),
                                false);
    relSparseSet(DataLayout::CELL_DOM) = RelationSetType(&cellMatRel);
    relDenseSet(DataLayout::CELL_DOM) = ProductSetType(&getCellSet(), &getMatSet());
  }
  if(other.hasValidStaticRelation(DataLayout::MAT_DOM))
  {
    StaticVariableRelationType& matCellRel = relStatic(DataLayout::MAT_DOM);

    matCellRel = StaticVariableRelationType(&getMatSet(), &getCellSet());
    matCellRel.bindBeginOffsets(getMatSet().size(), m_matCellRel_beginsVec.view());
    matCellRel.bindFirstIndices(m_matCellRel_firstIndicesVec.size(),
                                m_matCellRel_firstIndicesVec.view(),
                                false);
    matCellRel.bindIndices(m_matCellRel_indicesVec.size(), m_matCellRel_indicesVec.view());
    relSparseSet(DataLayout::MAT_DOM) = RelationSetType(&matCellRel);
    relDenseSet(DataLayout::MAT_DOM) = ProductSetType(&getMatSet(), &getCellSet());
  }
  for(size_t idx = 0; idx < other.m_fieldBackingVec.size(); idx++)
  {
    if(other.m_fieldBackingVec[idx] == nullptr)
    {
      m_fieldBackingVec.push_back(nullptr);
    }
    else if(other.m_fieldBackingVec[idx]->isOwned())
    {
      m_fieldBackingVec.emplace_back(new FieldBacking(*(other.m_fieldBackingVec[idx])));
    }
    else
    {
      m_fieldBackingVec.emplace_back(new FieldBacking);
      switch(other.m_dataTypeVec[idx])
      {
      case DataTypeSupported::TypeFloat:
        *(m_fieldBackingVec[idx]) =
          FieldBacking(other.m_fieldBackingVec[idx]->getArrayView<float>(), true, m_fieldAllocatorId);
        break;
      case DataTypeSupported::TypeDouble:
        *(m_fieldBackingVec[idx]) =
          FieldBacking(other.m_fieldBackingVec[idx]->getArrayView<double>(), true, m_fieldAllocatorId);
        break;
      case DataTypeSupported::TypeInt:
        *(m_fieldBackingVec[idx]) =
          FieldBacking(other.m_fieldBackingVec[idx]->getArrayView<int>(), true, m_fieldAllocatorId);
        break;
      case DataTypeSupported::TypeUnsignChar:
        *(m_fieldBackingVec[idx]) =
          FieldBacking(other.m_fieldBackingVec[idx]->getArrayView<unsigned char>(),
                       true,
                       m_fieldAllocatorId);
        break;
      case DataTypeSupported::TypeUnknown:
      default:
        SLIC_ERROR("Multimat: Unknown field type for field \"" << other.m_fieldNameVec[idx] << "\"");
        break;
      }
    }
  }
}

void MultiMat::setNumberOfMaterials(int n)
{
  SLIC_ASSERT(n > 0);
  m_nmats = n;

  getMatSet() = RangeSetType(0, m_nmats);
  SLIC_ASSERT(getMatSet().isValid());
}

void MultiMat::setNumberOfCells(int c)
{
  SLIC_ASSERT(c > 0);
  m_ncells = c;

  getCellSet() = RangeSetType(0, m_ncells);
  SLIC_ASSERT(getCellSet().isValid());
}

void MultiMat::setAllocatorID(int alloc_id)
{
  setSlamAllocatorID(alloc_id);
  setFieldAllocatorID(alloc_id);
}

void MultiMat::setSlamAllocatorID(int alloc_id)
{
  const bool hasCellDomRelation = hasValidStaticRelation(DataLayout::CELL_DOM);
  const bool hasMatDomRelation = hasValidStaticRelation(DataLayout::MAT_DOM);

  m_slamAllocatorId = alloc_id;
  m_sets = axom::Array<RangeSetType>(m_sets, m_slamAllocatorId);

  m_cellMatRel_beginsVec = IndBufferType(m_cellMatRel_beginsVec, m_slamAllocatorId);
  m_cellMatRel_indicesVec = IndBufferType(m_cellMatRel_indicesVec, m_slamAllocatorId);
  m_cellMatRel_firstIndicesVec = IndBufferType(m_cellMatRel_firstIndicesVec, m_slamAllocatorId);

  m_matCellRel_beginsVec = IndBufferType(m_matCellRel_beginsVec, m_slamAllocatorId);
  m_matCellRel_indicesVec = IndBufferType(m_matCellRel_indicesVec, m_slamAllocatorId);
  m_matCellRel_firstIndicesVec = IndBufferType(m_matCellRel_firstIndicesVec, m_slamAllocatorId);

  m_staticRelations = axom::Array<StaticVariableRelationType>(m_staticRelations, m_slamAllocatorId);
  m_dynamicRelations =
    axom::Array<DynamicVariableRelationType>(m_dynamicRelations, m_slamAllocatorId);
  m_sparseBivarSet = axom::Array<RelationSetType>(m_sparseBivarSet, m_slamAllocatorId);
  m_denseBivarSet = axom::Array<ProductSetType>(m_denseBivarSet, m_slamAllocatorId);

  m_flatCellToMatIndexMap = IndBufferType(m_flatCellToMatIndexMap, m_slamAllocatorId);
  m_flatMatToCellIndexMap = IndBufferType(m_flatMatToCellIndexMap, m_slamAllocatorId);

  if(hasCellDomRelation)
  {
    StaticVariableRelationType& cellMatRel = relStatic(DataLayout::CELL_DOM);

    cellMatRel = StaticVariableRelationType(&getCellSet(), &getMatSet());
    cellMatRel.bindBeginOffsets(getCellSet().size(), m_cellMatRel_beginsVec.view());
    cellMatRel.bindIndices(m_cellMatRel_indicesVec.size(), m_cellMatRel_indicesVec.view());
    cellMatRel.bindFirstIndices(m_cellMatRel_firstIndicesVec.size(),
                                m_cellMatRel_firstIndicesVec.view(),
                                false);
    relSparseSet(DataLayout::CELL_DOM) = RelationSetType(&cellMatRel);
    relDenseSet(DataLayout::CELL_DOM) = ProductSetType(&getCellSet(), &getMatSet());
  }
  if(hasMatDomRelation)
  {
    StaticVariableRelationType& matCellRel = relStatic(DataLayout::MAT_DOM);

    matCellRel = StaticVariableRelationType(&getMatSet(), &getCellSet());
    matCellRel.bindBeginOffsets(getMatSet().size(), m_matCellRel_beginsVec.view());
    matCellRel.bindIndices(m_matCellRel_indicesVec.size(), m_matCellRel_indicesVec.view());
    matCellRel.bindFirstIndices(m_matCellRel_firstIndicesVec.size(),
                                m_matCellRel_firstIndicesVec.view(),
                                false);
    relSparseSet(DataLayout::MAT_DOM) = RelationSetType(&matCellRel);
    relDenseSet(DataLayout::MAT_DOM) = ProductSetType(&getMatSet(), &getCellSet());
  }
}

void MultiMat::setFieldAllocatorID(int alloc_id)
{
  m_fieldAllocatorId = alloc_id;
  for(size_t idx = 0; idx < m_fieldBackingVec.size(); idx++)
  {
    if(m_fieldBackingVec[idx] != nullptr)
    {
      m_fieldBackingVec[idx]->moveSpaces(m_fieldAllocatorId);
    }
  }
}

void MultiMat::setCellMatRel(const std::vector<bool>& vecarr, DataLayout layout)
{
  //Setup the SLAM cell to mat relation
  //This step is necessary if the volfrac field is sparse

  SLIC_ASSERT(vecarr.size() == m_ncells * m_nmats);  //Check it's dense
  SLIC_ASSERT(!hasValidStaticRelation(layout));

  RangeSetType& set1 = relDominantSet(layout);
  RangeSetType& set2 = relSecondarySet(layout);

  axom::Array<SetPosType> counts(set1.size()), indices;

  for(SetPosType i = 0; i < set1.size(); ++i)
  {
    SetPosType cardinality = 0;
    for(SetPosType j = 0; j < set2.size(); ++j)
    {
      if(vecarr[i * set2.size() + j])
      {
        cardinality++;
        indices.push_back(j);
      }
    }
    counts[i] = cardinality;
  }

  if(AllocatorOnDevice(m_slamAllocatorId))
  {
    // Copy counts and indices to GPU memory.
    counts = IndBufferType(counts, m_slamAllocatorId);
    indices = IndBufferType(indices, m_slamAllocatorId);
  }

  setCellMatRel(counts, indices, layout);
}

/*!
 * \brief Given the cardinality array from a cell-material relation, generates
 *  arrays to be bound to the associated Slam StaticRelation type.
 *
 * \param [in] counts The number of second-set elements associated with the
 *                    first-set element
 * \param [out] begins For each first-set element, the first offset for a range of
 *                     associated second-set elements
 * \param [out] firstIndices For each flat index, the associated first-set element
 *
 * \pre Requires that Multimat is built with support for RAJA and Umpire
 */
template <typename ExecSpace, typename IndexType>
void ScanRelationOffsetsRAJA(const axom::ArrayView<const IndexType> counts,
                             const axom::ArrayView<IndexType> begins,
                             const axom::ArrayView<IndexType> firstIndices)
{
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

  // The begins array has one more element than the counts array. The last
  // element is the size of the relation, so we can just inclusive scan with an
  // offset to the begins array.
  axom::inclusive_scan<ExecSpace>(counts,
                                  axom::ArrayView<IndexType>(begins.data() + 1, counts.size()));

  // Generate the first indices array on the GPU.
  axom::for_all<ExecSpace>(
    counts.size(),
    AXOM_LAMBDA(int firstIndex) {
      for(int flatIndex = begins[firstIndex]; flatIndex < begins[firstIndex + 1]; flatIndex++)
      {
        firstIndices[flatIndex] = firstIndex;
      }
    });
#else
  AXOM_UNUSED_VAR(counts);
  AXOM_UNUSED_VAR(begins);
  AXOM_UNUSED_VAR(firstIndices);
  SLIC_ASSERT_MSG(false, "Calling ScanRelationOffsetsRAJA requires support for RAJA and Umpire.");
#endif
}

void MultiMat::setCellMatRel(axom::ArrayView<const SetPosType> cardinality,
                             axom::ArrayView<const SetPosType> indices,
                             DataLayout layout)
{
  StaticVariableRelationType& Rel_ptr = relStatic(layout);
  IndBufferType& Rel_beginsVec = relBeginVec(layout);
  IndBufferType& Rel_indicesVec = relIndVec(layout);
  IndBufferType& Rel_firstIndicesVec = relFirstIndVec(layout);

  // Check that we haven't already set up the cell-material relation.
  // TODO: should we allow resetting the relation?
  SLIC_ASSERT(!hasValidStaticRelation(layout));

  RangeSetType& set1 = relDominantSet(layout);
  RangeSetType& set2 = relSecondarySet(layout);

  SLIC_ASSERT(set1.size() == cardinality.size());

  // Offsets vector is just an inclusive scan of the cardinality
  Rel_beginsVec.resize(set1.size() + 1);
  Rel_firstIndicesVec.resize(indices.size());
  if(AllocatorOnDevice(m_slamAllocatorId))
  {
    auto begins = Rel_beginsVec.view();
    auto firstIndices = Rel_firstIndicesVec.view();
    ScanRelationOffsetsRAJA<GPU_Exec>(cardinality, begins, firstIndices);
  }
  else
  {
    Rel_beginsVec[0] = 0;
    for(SetPosType maj_idx = 1; maj_idx < cardinality.size() + 1; maj_idx++)
    {
      Rel_beginsVec[maj_idx] = Rel_beginsVec[maj_idx - 1] + cardinality[maj_idx - 1];
    }
  }

  // Just copy over the indices directly
  Rel_indicesVec = axom::Array<SetPosType>(indices, m_slamAllocatorId);

  Rel_ptr = StaticVariableRelationType(&set1, &set2);
  Rel_ptr.bindBeginOffsets(set1.size(), Rel_beginsVec.view());
  Rel_ptr.bindIndices(Rel_indicesVec.size(), Rel_indicesVec.view());

  // If we are on the CPU, let the MappedVariableRelation policy construct the
  // first indices array.
  bool constructFirstIndices = !AllocatorOnDevice(m_slamAllocatorId);
  Rel_ptr.bindFirstIndices(Rel_firstIndicesVec.size(),
                           Rel_firstIndicesVec.view(),
                           constructFirstIndices);

  SLIC_ASSERT(relBeginVec(layout).getAllocatorID() == m_slamAllocatorId);
  SLIC_ASSERT(relIndVec(layout).getAllocatorID() == m_slamAllocatorId);
  SLIC_ASSERT(relFirstIndVec(layout).getAllocatorID() == m_slamAllocatorId);

  SLIC_ASSERT(Rel_ptr.isValid());

  //Set-up both dense and sparse BivariateSets.
  relSparseSet(layout) = RelationSetType(&Rel_ptr);
  relDenseSet(layout) = ProductSetType(&set1, &set2);

  //Create a field for VolFrac as the 0th field
  m_fieldNameVec.push_back("Volfrac");
  m_fieldBackingVec.emplace_back(new FieldBacking);
  m_fieldMappingVec.push_back(FieldMapping::PER_CELL_MAT);
  m_dataTypeVec.push_back(DataTypeSupported::TypeDouble);
  m_fieldDataLayoutVec.push_back(DataLayout::CELL_DOM);
  m_fieldSparsityLayoutVec.push_back(SparsityLayout::SPARSE);
  m_fieldStrideVec.push_back(1);

  SLIC_ASSERT(m_fieldNameVec.size() == 1);
  SLIC_ASSERT(m_fieldMappingVec.size() == 1);
  SLIC_ASSERT(m_dataTypeVec.size() == 1);
  SLIC_ASSERT(m_fieldDataLayoutVec.size() == 1);
  SLIC_ASSERT(m_fieldSparsityLayoutVec.size() == 1);
  SLIC_ASSERT(m_fieldStrideVec.size() == 1);
}

void MultiMat::removeField(const std::string& field_name)
{
  int field_index = getFieldIdx(field_name);

  if(field_index == 0)
  {
    // This is the volume fractions array. Issue a warning to the user.
    SLIC_WARNING("Multimat Error: cannot remove volume fractions array.");
  }
  else if(field_index < 0)
  {
    // Field does not exist.
    SLIC_WARNING("Multimat Error: field with name \"" << field_name << "\" does not exist.");
  }
  else
  {
    m_fieldNameVec.erase(m_fieldNameVec.begin() + field_index);
    m_fieldMappingVec.erase(m_fieldMappingVec.begin() + field_index);
    m_fieldBackingVec.erase(m_fieldBackingVec.begin() + field_index);
    m_fieldDataLayoutVec.erase(m_fieldDataLayoutVec.begin() + field_index);
    m_fieldSparsityLayoutVec.erase(m_fieldSparsityLayoutVec.begin() + field_index);
    m_fieldStrideVec.erase(m_fieldStrideVec.begin() + field_index);
    m_dataTypeVec.erase(m_dataTypeVec.begin() + field_index);
  }

  SLIC_ASSERT(m_fieldNameVec.size() == m_dataTypeVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldMappingVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldDataLayoutVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldSparsityLayoutVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldStrideVec.size());
}

int MultiMat::setVolfracField(axom::ArrayView<const double> arr,
                              DataLayout layout,
                              SparsityLayout sparsity)
{
  // m_mapVec[0] should already be a volfrac map. Copy over the array to the
  // backing entry.
  SLIC_ASSERT(m_fieldNameVec[0] == "Volfrac");
  SLIC_ASSERT(m_fieldStrideVec[0] == 1);
  SLIC_ASSERT(arr.size() == get_mapped_biSet(layout, sparsity)->size());

  *(m_fieldBackingVec[0]) = FieldBacking(arr, true, m_fieldAllocatorId);
  m_fieldDataLayoutVec[0] = layout;
  m_fieldSparsityLayoutVec[0] = sparsity;

  return 0;
}

MultiMat::Field2D<double> MultiMat::getVolfracField() { return get2dFieldImpl<double>(0); }

MultiMat::Field2D<const double> MultiMat::getVolfracField() const
{
  return get2dFieldImpl<const double>(0);
}

int MultiMat::getFieldIdx(const std::string& field_name) const
{
  for(unsigned int i = 0; i < m_fieldNameVec.size(); i++)
  {
    if(m_fieldNameVec[i] == field_name)
    {
      return i;
    }
  }

  return -1;
}

std::string MultiMat::getFieldName(int field_idx) const
{
  if(field_idx < 0 || field_idx >= static_cast<int>(m_fieldNameVec.size()))
  {
    return "";
  }
  return m_fieldNameVec[field_idx];
}

MultiMat::IdSet MultiMat::getMatInCell(int c)
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::CELL_DOM));

  return relStatic(DataLayout::CELL_DOM)[c];
}

MultiMat::IdSet MultiMat::getCellContainingMat(int m)
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::MAT_DOM));

  return relStatic(DataLayout::MAT_DOM)[m];
}

MultiMat::IndexSet MultiMat::getSubfieldIndexingSet(int idx,
                                                    DataLayout layout,
                                                    SparsityLayout sparsity) const
{
  if(layout == DataLayout::CELL_DOM)
  {
    return getIndexingSetOfCell(idx, sparsity);
  }
  else
  {
    return getIndexingSetOfMat(idx, sparsity);
  }
}

MultiMat::IndexSet MultiMat::getIndexingSetOfCell(int c, SparsityLayout sparsity) const
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::CELL_DOM));
  SLIC_ASSERT(0 <= c && c < (int)m_ncells);

  return get_mapped_biSet(DataLayout::CELL_DOM, sparsity)->elementRangeSet(c);
}

MultiMat::IndexSet MultiMat::getIndexingSetOfMat(int m, SparsityLayout sparsity) const
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::MAT_DOM));
  SLIC_ASSERT(0 <= m && m < (int)m_nmats);

  return get_mapped_biSet(DataLayout::MAT_DOM, sparsity)->elementRangeSet(m);
}

void MultiMat::convertToDynamic()
{
  if(m_dynamic_mode)
  {
    return;
  }
  // If we weren't in dynamic mode before, we shouldn't already have dynamic
  // relations
  SLIC_ASSERT(!hasValidDynamicRelation(DataLayout::CELL_DOM));
  SLIC_ASSERT(!hasValidDynamicRelation(DataLayout::MAT_DOM));

  // Save what the current layout is for later
  //m_static_layout = Layout { m_dataLayout, m_sparsityLayout };  //old version with single layout for all MM
  m_layout_when_static.resize(m_fieldDataLayoutVec.size());
  const int SZ = static_cast<int>(m_fieldDataLayoutVec.size());
  for(auto i = 0; i < SZ; i++)
  {
    m_layout_when_static[i].data_layout = m_fieldDataLayoutVec[i];
    m_layout_when_static[i].sparsity_layout = m_fieldSparsityLayoutVec[i];

    // For now, handle dynamic by changing maps to dense,
    // and relation to DynamicRelation
    // if(isSparse()) convertLayoutToDense(); //old version with single layout for all MM
    if(m_fieldSparsityLayoutVec[i] == SparsityLayout::SPARSE)
    {
      convertFieldToDense(i);
    }
  }

  //create the dynamic relation
  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidStaticRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }
    StaticVariableRelationType& rel = relStatic(layout);

    RangeSetType* set1 = &relDominantSet(layout);
    RangeSetType* set2 = &relSecondarySet(layout);

    DynamicVariableRelationType relDyn(set1, set2);
    for(int i = 0; i < rel.fromSetSize(); i++)
    {
      auto&& rel_vec = rel[i];
      for(int j = 0; j < rel_vec.size(); j++)
      {
        relDyn.insert(i, rel_vec[j]);
      }
    }

    relDynamic(layout) = relDyn;

    SLIC_ASSERT(relDynamic(layout).isValid());
    SLIC_ASSERT(relDynamic(layout).totalSize() == relStatic(layout).totalSize());
  }

  relStatic(DataLayout::CELL_DOM) = StaticVariableRelationType {};
  relSparseSet(DataLayout::CELL_DOM) = RelationSetType {};
  relStatic(DataLayout::MAT_DOM) = StaticVariableRelationType {};
  relSparseSet(DataLayout::MAT_DOM) = RelationSetType {};

  m_dynamic_mode = true;
}

void MultiMat::convertToStatic()
{
  if(!m_dynamic_mode)
  {
    return;
  }

  // Change dynamicRelation back to staticRelation
  // change the layout to previously stored static layout

  SLIC_ASSERT(!hasValidStaticRelation(DataLayout::CELL_DOM));
  SLIC_ASSERT(!hasValidStaticRelation(DataLayout::MAT_DOM));

  //Create the static relations
  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidDynamicRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }

    DynamicVariableRelationType& relDyn = relDynamic(layout);

    RangeSetType& set1 = relDominantSet(layout);
    RangeSetType& set2 = relSecondarySet(layout);

    IndBufferType& rel_beginvec = relBeginVec(layout);
    IndBufferType& rel_indicesVec = relIndVec(layout);
    IndBufferType& rel_firstIndicesVec = relFirstIndVec(layout);

    SLIC_ASSERT(SetPosType(rel_beginvec.size()) == set1.size() + 1);
    int rel_data_size = 0;
    for(int i = 0; i < relDyn.fromSetSize(); i++)
    {
      auto& rel_vec = relDyn[i];
      rel_beginvec[i] = rel_data_size;
      rel_data_size += static_cast<int>(rel_vec.size());
    }
    rel_beginvec.back() = rel_data_size;
    rel_indicesVec.resize(rel_data_size);
    rel_firstIndicesVec.resize(rel_data_size);
    int idx = 0;
    for(int i = 0; i < relDyn.fromSetSize(); i++)
    {
      auto& rel_vec = relDyn[i];
      for(unsigned int j = 0; j < rel_vec.size(); j++)
      {
        rel_indicesVec[idx++] = rel_vec[j];
      }
    }
    SLIC_ASSERT(idx == relDyn.totalSize());

    StaticVariableRelationType& rel = relStatic(layout);

    rel = StaticVariableRelationType(&set1, &set2);
    rel.bindBeginOffsets(set1.size(), rel_beginvec.view());
    rel.bindIndices(rel_indicesVec.size(), rel_indicesVec.view());
    rel.bindFirstIndices(rel_firstIndicesVec.size(), rel_firstIndicesVec.view());

    SLIC_ASSERT(rel.isValid());

    RelationSetType& nzSet = relSparseSet(layout);
    SLIC_ASSERT(nzSet.getRelation() == nullptr);
    nzSet = RelationSetType(&rel);
    SLIC_ASSERT(nzSet.isValid());
  }

  m_dynamic_mode = false;

  //Change each field to their corresponding sparsity
  //if (m_static_layout.sparsity_layout == SparsityLayout::SPARSE) convertLayoutToSparse();
  const int SZ = static_cast<int>(m_layout_when_static.size());
  for(auto i = 0; i < SZ; i++)
  {
    if(m_layout_when_static[i].sparsity_layout == SparsityLayout::SPARSE)
    {
      convertFieldToSparse(i);
    }
  }
  m_layout_when_static.resize(0);

  // Unset the dynamic relations
  relDynamic(DataLayout::CELL_DOM) = DynamicVariableRelationType {};
  relDynamic(DataLayout::MAT_DOM) = DynamicVariableRelationType {};
}

bool MultiMat::addEntry(int cell_id, int mat_id)
{
  SLIC_ASSERT(m_dynamic_mode);

  //right now this is implemented by converting the data layout to dense so that
  // adding an entry only involve changing the relation data.

  bool searched = false;
  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidDynamicRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }

    DynamicVariableRelationType& relDyn = relDynamic(layout);

    std::pair<int, int> idx;
    if(layout == DataLayout::CELL_DOM)
    {
      idx = {cell_id, mat_id};
    }
    else
    {
      idx = {mat_id, cell_id};
    }

    auto& rel_vec = relDyn.data(idx.first);
    if(!searched)
    {
      // Search through the first valid relation we encounter to check if the
      // entry already exists
      auto&& found_iter = std::find(rel_vec.begin(), rel_vec.end(), idx.second);
      if(found_iter != rel_vec.end())
      {
        SLIC_ASSERT_MSG(false, "MultiMat::addEntry() -- entry already exists.");
        return false;
      }
    }

    relDyn.insert(idx.first, idx.second);
  }

  return true;
}

bool MultiMat::removeEntry(int cell_id, int mat_id)
{
  SLIC_ASSERT(m_dynamic_mode);

  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidDynamicRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }

    DynamicVariableRelationType& relDyn = relDynamic(layout);

    std::pair<int, int> idx;
    if(layout == DataLayout::CELL_DOM)
    {
      idx = {cell_id, mat_id};
    }
    else
    {
      idx = {mat_id, cell_id};
    }

    auto& rel_vec = relDyn.data(idx.first);
    auto&& found_iter = std::find(rel_vec.begin(), rel_vec.end(), idx.second);
    if(found_iter == rel_vec.end())
    {
      SLIC_ASSERT_MSG(false, "MultiMat::removeEntry() -- entry not found");
      return false;
    }

    rel_vec.erase(found_iter);
  }

  //TODO make all map value zero?

  return true;
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
/*!
 * \brief Transposes a sparse relation between cells and materials, generating
 *  the complementary data layout. Also generates transposition maps for flat
 *  indices, enabling quick transpositions of fields between cell-dominant and
 *  material-dominant layouts.
 *
 *  For example, a cell-dominant relation would be transposed into a material-
 *  dominant relation, and vice versa.
 *
 *  This function may run on the GPU if an appropriate execution space is
 *  provided.
 *
 * \param [in] oldRelationSet the sparse relation to transpose
 * \param [out] beginOffsets the new relation's begin offsets
 * \param [out] secondIndexes the new relation's second-set indices map
 * \param [out] firstIndexes the new relation's first-set indices map
 * \param [out] flatOldToNew mapping of flat indices from old to new relation
 * \param [out] flatNewToOld mapping of flat indices from old to new relation
 *
 * \tparam ExecSpace the execution space to run the transposition in
 *
 * \pre Requires that Multimat is built with support for RAJA and Umpire
 */
template <typename ExecSpace, typename IndexType>
void TransposeRelationImplRAJA(const MultiMat::RelationSetType* oldRelationSet,
                               axom::Array<IndexType>& beginOffsets,
                               axom::Array<IndexType>& secondIndexes,
                               axom::Array<IndexType>& firstIndexes,
                               axom::Array<IndexType>& flatOldToNew,
                               axom::Array<IndexType>& flatNewToOld)
{
  const auto countsView = beginOffsets.view();
  const auto firstIdxView = firstIndexes.view();
  const auto flatNewToOldView = flatNewToOld.view();

  // Count the number of entries for each second set index
  axom::for_all<ExecSpace>(
    oldRelationSet->totalSize(),
    AXOM_LAMBDA(int flatIndex) {
      int secondIdx = oldRelationSet->flatToSecondIndex(flatIndex);
      axom::atomicAdd<ExecSpace>(&countsView[secondIdx], IndexType {1});

      // We create the first indices array and flatNewToOld maps here.
      firstIdxView[flatIndex] = secondIdx;
      flatNewToOldView[flatIndex] = flatIndex;
    });

  // Scan to get the offsets array.
  axom::exclusive_scan_inplace<ExecSpace>(beginOffsets);

  // Stable sort to create the first-set indices array and new-to-old index
  // mapping for the new relation.
  axom::stable_sort_pairs<ExecSpace>(firstIndexes, flatNewToOld);

  const auto secondIdxView = secondIndexes.view();
  const auto flatOldToNewView = flatOldToNew.view();

  // With the new-to-old map, we can now fill in the second-set indices array
  // and the old-to-new map.
  axom::for_all<ExecSpace>(
    oldRelationSet->totalSize(),
    AXOM_LAMBDA(int newFlatIndex) {
      int oldFlatIndex = flatNewToOldView[newFlatIndex];

      int firstIdx = oldRelationSet->flatToFirstIndex(oldFlatIndex);
      secondIdxView[newFlatIndex] = firstIdx;
      flatOldToNewView[oldFlatIndex] = newFlatIndex;
    });
}
#endif  // defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

/*!
 * \brief Transposes a sparse relation between cells and materials, generating
 *  the complementary data layout. Also generates transposition maps for flat
 *  indices, enabling quick transpositions of fields between cell-dominant and
 *  material-dominant layouts.
 *
 *  For example, a cell-dominant relation would be transposed into a material-
 *  dominant relation, and vice versa.
 *
 * \param [in] oldRelationSet the sparse relation to transpose
 * \param [out] beginOffsets the new relation's begin offsets
 * \param [out] secondIndexes the new relation's second-set indices map
 * \param [out] firstIndexes the new relation's first-set indices map
 * \param [out] flatOldToNew mapping of flat indices from old to new relation
 * \param [out] flatNewToOld mapping of flat indices from old to new relation
 */
template <typename ExecSpace, typename IndexType>
void TransposeRelationImpl(const MultiMat::RelationSetType* oldRelationSet,
                           axom::Array<IndexType>& beginOffsets,
                           axom::Array<IndexType>& secondIndexes,
                           axom::Array<IndexType>& firstIndexes,
                           axom::Array<IndexType>& flatOldToNew,
                           axom::Array<IndexType>& flatNewToOld,
                           int slamAllocatorID)
{
  using IndBufferType = axom::Array<IndexType>;

  int numCols = oldRelationSet->secondSetSize() + 1;
  int relationSize = oldRelationSet->totalSize();
  beginOffsets = IndBufferType(numCols, numCols, slamAllocatorID);
  firstIndexes = IndBufferType(relationSize, relationSize, slamAllocatorID);
  secondIndexes = IndBufferType(relationSize, relationSize, slamAllocatorID);
  flatOldToNew = IndBufferType(relationSize, relationSize, slamAllocatorID);
  flatNewToOld = IndBufferType(relationSize, relationSize, slamAllocatorID);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  if(AllocatorOnDevice(slamAllocatorID))
  {
    TransposeRelationImplRAJA<GPU_Exec>(oldRelationSet,
                                        beginOffsets,
                                        secondIndexes,
                                        firstIndexes,
                                        flatOldToNew,
                                        flatNewToOld);
    return;
  }
#endif

  IndBufferType counts(numCols, numCols, slamAllocatorID);
  const auto countsView = counts.view();

  // Count the number of entries for each second set index
  axom::for_all<axom::SEQ_EXEC>(
    oldRelationSet->totalSize(),
    AXOM_LAMBDA(int flatIndex) {
      int secondIdx = oldRelationSet->flatToSecondIndex(flatIndex);
      countsView[secondIdx]++;
    });

  const auto beginView = beginOffsets.view();

  // Scan to get begin offsets array
  axom::for_all<axom::SEQ_EXEC>(
    oldRelationSet->secondSetSize(),
    AXOM_LAMBDA(int secondIdx) {
      beginView[secondIdx + 1] = beginView[secondIdx] + countsView[secondIdx];
    });

  const auto firstIdxView = firstIndexes.view();
  const auto secondIdxView = secondIndexes.view();

  const auto flatOldToNewView = flatOldToNew.view();
  const auto flatNewToOldView = flatNewToOld.view();

  // Fill first and second index arrays
  axom::for_all<axom::SEQ_EXEC>(
    oldRelationSet->totalSize(),
    AXOM_LAMBDA(int oldFlatIndex) {
      int firstIdx = oldRelationSet->flatToFirstIndex(oldFlatIndex);
      int secondIdx = oldRelationSet->flatToSecondIndex(oldFlatIndex);

      // Start from end of the second set element's range.
      int reverse_offset = beginView[secondIdx + 1];

      // Use countsView to keep track of indexes to place elements at.
      int newFlatIndex = reverse_offset - countsView[secondIdx];
      countsView[secondIdx]--;

      firstIdxView[newFlatIndex] = secondIdx;
      secondIdxView[newFlatIndex] = firstIdx;

      // Build transposition maps while we have both the old and the new flat index.
      flatOldToNewView[oldFlatIndex] = newFlatIndex;
      flatNewToOldView[newFlatIndex] = oldFlatIndex;
    });
}

void MultiMat::makeOtherRelation(DataLayout layout)
{
  DataLayout old_layout =
    (layout == DataLayout::CELL_DOM ? DataLayout::MAT_DOM : DataLayout::CELL_DOM);
  RelationSetType* oldRelationSet = &relSparseSet(old_layout);
  IndBufferType& newBeginVec = relBeginVec(layout);
  IndBufferType& newIndicesVec = relIndVec(layout);
  IndBufferType& newFirstIndicesVec = relFirstIndVec(layout);

  //construct the new transposed relation
  if(old_layout == DataLayout::CELL_DOM)
  {
    TransposeRelationImpl<axom::SEQ_EXEC>(oldRelationSet,
                                          newBeginVec,
                                          newIndicesVec,
                                          newFirstIndicesVec,
                                          m_flatCellToMatIndexMap,
                                          m_flatMatToCellIndexMap,
                                          m_slamAllocatorId);
  }
  else
  {
    TransposeRelationImpl<axom::SEQ_EXEC>(oldRelationSet,
                                          newBeginVec,
                                          newIndicesVec,
                                          newFirstIndicesVec,
                                          m_flatMatToCellIndexMap,
                                          m_flatCellToMatIndexMap,
                                          m_slamAllocatorId);
  }

  RangeSetType& newFirstSet = relDominantSet(layout);
  RangeSetType& newSecondSet = relSecondarySet(layout);

  SLIC_ASSERT(newBeginVec.size() == newFirstSet.size() + 1);

  StaticVariableRelationType& newRel = relStatic(layout);
  newRel = StaticVariableRelationType(&newFirstSet, &newSecondSet);
  newRel.bindBeginOffsets(newFirstSet.size(), newBeginVec.view());
  newRel.bindIndices(newIndicesVec.size(), newIndicesVec.view());
  newRel.bindFirstIndices(newFirstIndicesVec.size(), newFirstIndicesVec.view(), false);

  relSparseSet(layout) = RelationSetType(&newRel);
  relDenseSet(layout) = ProductSetType(&newFirstSet, &newSecondSet);
}

void MultiMat::convertLayoutToCellDominant()
{
  for(unsigned int i = 0; i < m_fieldNameVec.size(); ++i)
  {
    convertFieldToCellDom(i);
  }
}

void MultiMat::convertLayoutToMaterialDominant()
{
  for(unsigned int i = 0; i < m_fieldNameVec.size(); ++i)
  {
    convertFieldToMatDom(i);
  }
}

void MultiMat::convertLayoutToSparse()
{
  for(unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    convertFieldToSparse(map_i);
  }
}

void MultiMat::convertLayoutToDense()
{
  for(unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    convertFieldToDense(map_i);
  }
}

void MultiMat::convertLayout(DataLayout new_layout, SparsityLayout new_sparsity)
{
  //sparse/dense conversion
  if(new_sparsity == SparsityLayout::SPARSE)
  {
    convertLayoutToSparse();
  }
  else if(new_sparsity == SparsityLayout::DENSE)
  {
    convertLayoutToDense();
  }

  //cell/mat centric conversion
  if(new_layout == DataLayout::MAT_DOM)
  {
    convertLayoutToMaterialDominant();
  }
  else if(new_layout == DataLayout::CELL_DOM)
  {
    convertLayoutToCellDominant();
  }
}

void MultiMat::convertFieldLayout(int field_idx, SparsityLayout new_sparsity, DataLayout new_layout)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_fieldNameVec.size()));

  DataLayout field_data_layout = m_fieldDataLayoutVec[field_idx];
  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  if(new_layout == field_data_layout && new_sparsity == field_sparsity_layout)
  {
    return;
  }

  //sparse/dense conversion
  if(field_sparsity_layout == SparsityLayout::DENSE && new_sparsity == SparsityLayout::SPARSE)
  {
    convertFieldToSparse(field_idx);
  }
  else if(field_sparsity_layout == SparsityLayout::SPARSE && new_sparsity == SparsityLayout::DENSE)
  {
    convertFieldToDense(field_idx);
  }

  //cell/mat centric conversion
  if(field_data_layout == DataLayout::CELL_DOM && new_layout == DataLayout::MAT_DOM)
  {
    convertFieldToMatDom(field_idx);
  }
  else if(field_data_layout == DataLayout::MAT_DOM && new_layout == DataLayout::CELL_DOM)
  {
    convertFieldToCellDom(field_idx);
  }

  SLIC_ASSERT(false);
}

void MultiMat::convertFieldToSparse(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_fieldNameVec.size()));

  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  //no conversion needed unless the field is PER_CELL_MAT
  if(field_sparsity_layout == SparsityLayout::SPARSE ||
     m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT)
  {
    return;
  }

  if(!m_fieldBackingVec[field_idx]->isOwned())
  {
    SLIC_WARNING("Multimat: cannot convert unowned field \"" + m_fieldNameVec[field_idx] +
                 "\" to sparse layout. Skipping.");
    return;
  }

  if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeDouble)
  {
    convertToSparse_helper<double>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeFloat)
  {
    convertToSparse_helper<float>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeInt)
  {
    convertToSparse_helper<int>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeUnsignChar)
  {
    convertToSparse_helper<unsigned char>(field_idx);
  }
  else
  {
    SLIC_ASSERT(false);  //TODO
  }

  m_fieldSparsityLayoutVec[field_idx] = SparsityLayout::SPARSE;
}

void MultiMat::convertFieldToDense(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_fieldNameVec.size()));

  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  //no conversion needed unless the field is PER_CELL_MAT
  if(field_sparsity_layout == SparsityLayout::DENSE ||
     m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT)
  {
    return;
  }

  if(!m_fieldBackingVec[field_idx]->isOwned())
  {
    SLIC_WARNING("Multimat: cannot convert unowned field \"" + m_fieldNameVec[field_idx] +
                 "\" to dense layout. Skipping.");
    return;
  }

  if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeDouble)
  {
    convertToDense_helper<double>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeFloat)
  {
    convertToDense_helper<float>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeInt)
  {
    convertToDense_helper<int>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeUnsignChar)
  {
    convertToDense_helper<unsigned char>(field_idx);
  }
  else
  {
    SLIC_ASSERT(false);  //TODO
  }

  m_fieldSparsityLayoutVec[field_idx] = SparsityLayout::DENSE;
}

/*!
 * \brief Converts a dense-layout field to a sparse-layout field.
 *
 * \param [in] oldField the field to convert
 * \param [in] prodSet the associated relation set for the new field
 * \param [in] allocatorId allocator to use for the new array
 *
 * \return the associated sparse-layout data
 */
template <typename DataType>
axom::Array<DataType> ConvertToSparseImpl(const MultiMat::DenseField2D<DataType> oldField,
                                          const MultiMat::RelationSetType* relationSet,
                                          int allocatorId)
{
  int stride = oldField.stride();
  int sparseSize = relationSet->totalSize() * stride;
  axom::Array<DataType> sparseField(sparseSize, sparseSize, allocatorId);
  const auto sparseFieldView = sparseField.view();

  ExecLambdaForMemory(
    relationSet->totalSize() * stride,
    allocatorId,
    AXOM_LAMBDA(int index) {
      int flatIdx = index / stride;
      int comp = index % stride;

      auto firstIdx = relationSet->flatToFirstIndex(flatIdx);
      auto secondIdx = relationSet->flatToSecondIndex(flatIdx);

      sparseFieldView[index] = oldField(firstIdx, secondIdx, comp);
    });

  return sparseField;
}

template <typename DataType>
void MultiMat::convertToSparse_helper(int map_i)
{
  SLIC_ASSERT(m_fieldSparsityLayoutVec[map_i] != SparsityLayout::SPARSE);
  SLIC_ASSERT(m_fieldBackingVec[map_i]->isOwned());

  //Skip if no volume fraction array is set-up
  if(map_i == 0 && m_fieldBackingVec[0] == nullptr)
  {
    return;
  }

  const RelationSetType* rel_set = &relSparseSet(m_fieldDataLayoutVec[map_i]);

  DenseField2D<DataType> dense_field = getDense2dField<DataType>(m_fieldNameVec[map_i]);

  axom::Array<DataType> sparseFieldData =
    ConvertToSparseImpl(dense_field, rel_set, m_fieldAllocatorId);

  m_fieldBackingVec[map_i]->getArray<DataType>() = std::move(sparseFieldData);
}

/*!
 * \brief Converts a sparse-layout field to a dense-layout field.
 *
 * \param [in] oldField the field to convert
 * \param [in] prodSet the associated product set for the new field
 * \param [in] allocatorId allocator to use for the new array
 *
 * \return the associated dense-layout data
 */
template <typename DataType>
axom::Array<DataType> ConvertToDenseImpl(const MultiMat::SparseField2D<DataType> oldField,
                                         const MultiMat::ProductSetType* prodSet,
                                         int allocatorId)
{
  int stride = oldField.stride();
  int denseSize = prodSet->size() * stride;
  axom::Array<DataType> denseField(denseSize, denseSize, allocatorId);
  const auto denseFieldView = denseField.view();
  const auto* relationSetHost = oldField.set();

  ExecLambdaForMemory(
    relationSetHost->totalSize() * stride,
    allocatorId,
    AXOM_LAMBDA(int index) {
      int flatIdx = index / stride;
      int comp = index % stride;

      const auto* relationSet = oldField.set();
      auto firstIdx = relationSet->flatToFirstIndex(flatIdx);
      auto secondIdx = relationSet->flatToSecondIndex(flatIdx);

      int denseIdx = prodSet->findElementFlatIndex(firstIdx, secondIdx);

      denseFieldView[denseIdx * stride + comp] = oldField[index];
    });

  return denseField;
}

template <typename DataType>
void MultiMat::convertToDense_helper(int map_i)
{
  SLIC_ASSERT(m_fieldSparsityLayoutVec[map_i] != SparsityLayout::DENSE);
  SLIC_ASSERT(m_fieldBackingVec[map_i]->isOwned());

  //Skip if no volume fraction array is set-up
  if(map_i == 0 && m_fieldBackingVec[0] == nullptr)
  {
    return;
  }

  ProductSetType* prod_set = &relDenseSet(m_fieldDataLayoutVec[map_i]);

  SparseField2D<DataType> oldField = getSparse2dField<DataType>(m_fieldNameVec[map_i]);

  axom::Array<DataType> denseFieldData = ConvertToDenseImpl(oldField, prod_set, m_fieldAllocatorId);

  m_fieldBackingVec[map_i]->getArray<DataType>() = std::move(denseFieldData);
}

DataLayout MultiMat::getFieldDataLayout(int field_idx) const
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_fieldNameVec.size()));

  return m_fieldDataLayoutVec[field_idx];
}

SparsityLayout MultiMat::getFieldSparsityLayout(int field_idx) const
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_fieldNameVec.size()));

  return m_fieldSparsityLayoutVec[field_idx];
}

/*!
 * \brief Transposes a sparse field to a complementary data layout.
 *
 * \param [in] oldField the field to convert
 * \param [in] relationSet the associated relation set for the field
 * \param [in] flatTransposeMap the map of flat indices from the old relation
 *                              to the new relation
 * \param [in] allocatorId allocator to use for the new array
 *
 * \return the sparse field transposed into the complementary data layout
 */
template <typename DataType, typename IndexType>
axom::Array<DataType> TransposeSparseImpl(const MultiMat::SparseField2D<DataType> oldField,
                                          const MultiMat::RelationSetType* relationSet,
                                          const axom::Array<IndexType>& flatTransposeMap,
                                          int allocatorId)
{
  int stride = oldField.stride();
  int sparseSize = relationSet->totalSize() * stride;
  axom::Array<DataType> sparseField(sparseSize, sparseSize, allocatorId);

  const auto sparseFieldView = sparseField.view();
  const auto flatTransposeMapView = flatTransposeMap.view();

  ExecLambdaForMemory(
    relationSet->totalSize() * stride,
    allocatorId,
    AXOM_LAMBDA(int index) {
      int flatIdx = index / stride;
      int comp = index % stride;

      int newFlatIdx = flatTransposeMapView[flatIdx];

      sparseFieldView[newFlatIdx * stride + comp] = oldField[index];
    });

  return sparseField;
}

/*!
 * \brief Transposes a dense field to a complementary data layout.
 *
 * \param [in] oldField the field to convert
 * \param [in] relationSet the associated relation set for the field's layout
 * \param [in] allocatorId allocator to use for the new array
 *
 * \return the dense field transposed into the complementary data layout
 */
template <typename DataType>
axom::Array<DataType> TransposeDenseImpl(const MultiMat::DenseField2D<DataType> oldField,
                                         const MultiMat::RelationSetType* relationSet,
                                         int allocatorId)
{
  int stride = oldField.stride();

  const auto* firstSet = relationSet->getFirstSet();
  const auto* secondSet = relationSet->getSecondSet();
  int denseSize = firstSet->size() * secondSet->size() * stride;

  axom::Array<DataType> denseField(denseSize, denseSize, allocatorId);
  const auto denseFieldView = denseField.view();

  // Note: even though this is a dense field, we iterate over the relation set
  // in order to only copy over filled-in slots.
  ExecLambdaForMemory(
    relationSet->totalSize() * stride,
    allocatorId,
    AXOM_LAMBDA(int index) {
      int flatIdx = index / stride;
      int comp = index % stride;

      auto firstIdx = relationSet->flatToFirstIndex(flatIdx);
      auto secondIdx = relationSet->flatToSecondIndex(flatIdx);

      // Compute complementary dense index
      int denseIndex = secondIdx * firstSet->size() + firstIdx;
      denseIndex = denseIndex * stride + comp;

      denseFieldView[denseIndex] = oldField(firstIdx, secondIdx, comp);
    });

  return denseField;
}

template <typename DataType>
void MultiMat::transposeField_helper(int field_idx)
{
  SLIC_ASSERT(m_fieldMappingVec[field_idx] == FieldMapping::PER_CELL_MAT);

  //Skip if no volume fraction array is set-up
  if(field_idx == 0 && m_fieldBackingVec[0] == nullptr)
  {
    return;
  }

  DataLayout oldDataLayout = getFieldDataLayout(field_idx);
  DataLayout new_layout;
  if(oldDataLayout == DataLayout::CELL_DOM)
  {
    new_layout = DataLayout::MAT_DOM;
  }
  else
  {
    new_layout = DataLayout::CELL_DOM;
  }

  if(!hasValidStaticRelation(new_layout))
  {
    makeOtherRelation(new_layout);
  }

  axom::Array<DataType> arr_data;
  RelationSetType* fromRelSet = &relSparseSet(m_fieldDataLayoutVec[field_idx]);

  auto old_data_view = m_fieldBackingVec[field_idx]->getArrayView<DataType>();

  bool shallow_data_on_device = AllocatorOnDevice(old_data_view.getAllocatorID());
  bool field_allocator_on_device = AllocatorOnDevice(m_fieldAllocatorId);

  // An externally-managed array might have a memory space different from the
  // current allocator for internally-managed fields. We should create a
  // temporary copy of the field in order to correctly access the data on the
  // host or device.
  axom::Array<DataType> copied_old_data;
  if(shallow_data_on_device != field_allocator_on_device)
  {
    // Move data to expected field memory space.
    copied_old_data = axom::Array<DataType>(old_data_view, m_fieldAllocatorId);
    old_data_view = copied_old_data.view();
  }

  if(m_fieldSparsityLayoutVec[field_idx] == SparsityLayout::SPARSE)
  {
    typename RelationSetType::ConcreteSet rel_set = relSparseSet(oldDataLayout);
    SparseField2D<DataType> oldField(*this,
                                     rel_set,
                                     field_idx,
                                     old_data_view,
                                     m_fieldStrideVec[field_idx]);

    if(oldDataLayout == DataLayout::CELL_DOM)
    {
      arr_data =
        TransposeSparseImpl(oldField, fromRelSet, m_flatCellToMatIndexMap, m_fieldAllocatorId);
    }
    else
    {
      arr_data =
        TransposeSparseImpl(oldField, fromRelSet, m_flatMatToCellIndexMap, m_fieldAllocatorId);
    }
  }
  else  //dense
  {
    typename ProductSetType::ConcreteSet prod_set = relDenseSet(oldDataLayout);
    DenseField2D<DataType> oldField(*this,
                                    prod_set,
                                    field_idx,
                                    old_data_view,
                                    m_fieldStrideVec[field_idx]);

    arr_data = TransposeDenseImpl(oldField, fromRelSet, m_fieldAllocatorId);
  }

  SLIC_ASSERT(arr_data.getAllocatorID() == m_fieldAllocatorId);

  if(m_fieldBackingVec[field_idx]->isOwned())
  {
    m_fieldBackingVec[field_idx]->getArray<DataType>() = std::move(arr_data);
  }
  else
  {
    // We don't own the underlying buffer, just copy the data.
    const auto old_view = m_fieldBackingVec[field_idx]->getArrayView<DataType>();
    axom::copy(old_view.data(), arr_data.data(), old_view.size() * sizeof(DataType));
  }
  m_fieldDataLayoutVec[field_idx] = new_layout;
}

void MultiMat::transposeField(int field_idx)
{
  if(m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT)
  {
    return;
  }

  switch(m_dataTypeVec[field_idx])
  {
  case DataTypeSupported::TypeDouble:
    transposeField_helper<double>(field_idx);
    break;
  case DataTypeSupported::TypeFloat:
    transposeField_helper<float>(field_idx);
    break;
  case DataTypeSupported::TypeInt:
    transposeField_helper<int>(field_idx);
    break;
  case DataTypeSupported::TypeUnsignChar:
    transposeField_helper<unsigned char>(field_idx);
    break;
  default:
    SLIC_ASSERT(false);
  }
}

void MultiMat::convertFieldToMatDom(int field_idx)
{
  if(getFieldDataLayout(field_idx) == DataLayout::MAT_DOM)
  {
    return;
  }
  transposeField(field_idx);
}

void MultiMat::convertFieldToCellDom(int field_idx)
{
  if(getFieldDataLayout(field_idx) == DataLayout::CELL_DOM)
  {
    return;
  }
  transposeField(field_idx);
}

std::string MultiMat::getFieldMappingAsString(int field_i) const
{
  if(m_fieldMappingVec[field_i] == FieldMapping::PER_CELL)
  {
    return "Per-Cell";
  }
  else if(m_fieldMappingVec[field_i] == FieldMapping::PER_MAT)
  {
    return "Per-Material";
  }
  else if(m_fieldMappingVec[field_i] == FieldMapping::PER_CELL_MAT)
  {
    return "Per-Cell-Material";
  }
  else
  {
    SLIC_ASSERT(false);
  }
  return "";
}

std::string MultiMat::getFieldDataLayoutAsString(int field_i) const
{
  if(m_fieldDataLayoutVec[field_i] == DataLayout::CELL_DOM)
  {
    return "Cell-Dominant";
  }
  else if(m_fieldDataLayoutVec[field_i] == DataLayout::MAT_DOM)
  {
    return "Material-Dominant";
  }
  else
  {
    SLIC_ASSERT(false);
  }
  return "";
}

std::string MultiMat::getFieldSparsityLayoutAsString(int field_i) const
{
  if(m_fieldSparsityLayoutVec[field_i] == SparsityLayout::SPARSE)
  {
    return "Sparse";
  }
  else if(m_fieldSparsityLayoutVec[field_i] == SparsityLayout::DENSE)
  {
    return "Dense";
  }
  else
  {
    SLIC_ASSERT(false);
  }
  return "";
}

std::string MultiMat::getFieldDataTypeAsString(int field_i) const
{
  if(m_dataTypeVec[field_i] == DataTypeSupported::TypeUnknown)
  {
    return "Unknown";
  }
  else if(m_dataTypeVec[field_i] == DataTypeSupported::TypeInt)
  {
    return "int";
  }
  else if(m_dataTypeVec[field_i] == DataTypeSupported::TypeDouble)
  {
    return "double";
  }
  else if(m_dataTypeVec[field_i] == DataTypeSupported::TypeFloat)
  {
    return "float";
  }
  else if(m_dataTypeVec[field_i] == DataTypeSupported::TypeUnsignChar)
  {
    return "unsigned char";
  }
  else
  {
    SLIC_ASSERT(false);
  }
  return "";
}

void MultiMat::print() const
{
  std::stringstream sstr;

  sstr << "  Multimat Object Details:";
  sstr << "\nNumber of materials: " << m_nmats;
  sstr << "\nNumber of cells:     " << m_ncells;

  sstr << "\n\n Number of fields: " << m_fieldNameVec.size() << "\n";
  for(unsigned int i = 0; i < m_fieldNameVec.size(); i++)
  {
    sstr << "Field " << i << ": " << m_fieldNameVec[i].c_str();
    sstr << "  Mapping per ";
    switch(m_fieldMappingVec[i])
    {
    case FieldMapping::PER_CELL:
      sstr << "cell";
      break;
    case FieldMapping::PER_MAT:
      sstr << "material";
      break;
    case FieldMapping::PER_CELL_MAT:
      sstr << "cellXmaterial";
      sstr << "\n  Data layout: " << getFieldDataLayoutAsString(i);
      sstr << "\n  Sparsity layout: " << getFieldSparsityLayoutAsString(i);
      break;
    }
  }
  sstr << "\n\n";

  std::cout << sstr.str() << std::endl;
}

bool MultiMat::isValid(bool verboseOutput) const
{
  bool bValid = true;
  std::stringstream errStr;

  if(getCellSet().size() > 0 && getMatSet().size() > 0)
  {
    //It's a non-empty MultiMat object.
    //Check Volfrac exists
    if(m_fieldBackingVec[0] == nullptr)
    {
      errStr << "\n\t*No Volfrac field added.";
      bValid = false;
    }
    else
    {
      // TODO: fix constness of multimat
      auto volfrac_map = getVolfracField();
      auto volfrac_sparsity = m_fieldSparsityLayoutVec[0];
      auto volfrac_layout = m_fieldDataLayoutVec[0];

      //Check Volfrac values match the relation value (if dense)
      if(volfrac_sparsity == SparsityLayout::DENSE && !m_dynamic_mode)
      {
        const StaticVariableRelationType& relPtr = m_staticRelations[(int)volfrac_layout];

        for(int i = 0; i < volfrac_map.firstSetSize(); ++i)
        {
          auto rel_iter = relPtr.begin(i);
          for(int j = 0; j < volfrac_map.secondSetSize(); ++j)
          {
            double volfrac = volfrac_map[i * volfrac_map.secondSetSize() + j];
            bool zero_volfrac = volfrac == 0.0;  //exact comp?
            if(rel_iter != relPtr.end(i) && *rel_iter == j)
            {  //cellmat rel is present
              if(zero_volfrac)
              {
                errStr << "\n\t*Volume fraction is zero for a material "
                       << " that exists in a cell";
                bValid = false;
              }
              ++rel_iter;
            }
            else
            {
              if(!zero_volfrac)
              {
                errStr << "\n\t*Volume fraction is non-zero for a material "
                       << "not presented in a cell.";
                bValid = false;
              }
            }
          }
        }
      }

      //Check Volfrac sums to up 1.0
      std::vector<double> volfrac_sum(getCellSet().size(), 0.0);
      for(int i = 0; i < volfrac_map.firstSetSize(); ++i)
      {
        const auto& submap = volfrac_map(i);
        for(int j = 0; j < submap.size(); ++j)
        {
          if(volfrac_layout == DataLayout::CELL_DOM)
          {
            volfrac_sum[i] += submap.value(j);
          }
          else
          {
            volfrac_sum[submap.index(j)] += submap.value(j);
          }
        }
      }
      for(unsigned int i = 0; i < volfrac_sum.size(); ++i)
      {
        if(std::abs(volfrac_sum[i] - 1.0) > 10e-9)
        {
          errStr << "\n\t*Volfrac does not sum to 1.0 in cell " << i;
          bValid = false;
        }
      }
    }
  }

  if(verboseOutput)
  {
    if(bValid)
    {
      errStr << "\n\t*MultiMat data was valid";
    }

    std::cout << errStr.str() << std::endl;
  }

  return bValid;
}

const MultiMat::BivariateSetType* MultiMat::get_mapped_biSet(DataLayout layout,
                                                             SparsityLayout sparsity) const
{
  const BivariateSetType* set_ptr = nullptr;

  if(sparsity == SparsityLayout::SPARSE)
  {
    set_ptr = &relSparseSet(layout);
  }
  else if(sparsity == SparsityLayout::DENSE)
  {
    set_ptr = &relDenseSet(layout);
  }

  return set_ptr;
}

const MultiMat::BivariateSetType* MultiMat::get_mapped_biSet(int field_idx) const
{
  SLIC_ASSERT(m_fieldMappingVec[field_idx] == FieldMapping::PER_CELL_MAT);
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_fieldNameVec.size()));

  DataLayout layout = m_fieldDataLayoutVec[field_idx];
  SparsityLayout sparsity = m_fieldSparsityLayoutVec[field_idx];

  return get_mapped_biSet(layout, sparsity);
}

MultiMat::StaticVariableRelationType* MultiMat::getRel(int field_idx)
{
  SLIC_ASSERT(m_fieldDataLayoutVec[field_idx] == DataLayout::CELL_DOM ||
              m_fieldDataLayoutVec[field_idx] == DataLayout::MAT_DOM);

  //get the sparse relation vector
  DataLayout layout = m_fieldDataLayoutVec[field_idx];
  StaticVariableRelationType* rel_ptr = &relStatic(layout);

  return rel_ptr;
}

namespace axom
{
namespace multimat
{
std::ostream& operator<<(std::ostream& os, FieldMapping mapping)
{
  switch(mapping)
  {
  case FieldMapping::PER_CELL:
    os << "PER_CELL";
    break;
  case FieldMapping::PER_MAT:
    os << "PER_MAT";
    break;
  case FieldMapping::PER_CELL_MAT:
    os << "PER_CELL_MAT";
    break;
  default:
    os << "Unknown";
    break;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, DataLayout layout)
{
  switch(layout)
  {
  case DataLayout::CELL_DOM:
    os << "CELL_DOM";
    break;
  case DataLayout::MAT_DOM:
    os << "MAT_DOM";
    break;
  default:
    os << "Unknown";
    break;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, SparsityLayout sparsity)
{
  switch(sparsity)
  {
  case SparsityLayout::DENSE:
    os << "DENSE";
    break;
  case SparsityLayout::SPARSE:
    os << "SPARSE";
    break;
  default:
    os << "Unknown";
    break;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, DataTypeSupported type)
{
  switch(type)
  {
  case DataTypeSupported::TypeUnknown:
    os << "TypeUnknown";
    break;
  case DataTypeSupported::TypeInt:
    os << "TypeInt";
    break;
  case DataTypeSupported::TypeDouble:
    os << "TypeDouble";
    break;
  case DataTypeSupported::TypeFloat:
    os << "TypeFloat";
    break;
  case DataTypeSupported::TypeUnsignChar:
    os << "TypeUnsignChar";
    break;
  default:
    os << "Unknown";
    break;
  }
  return os;
}
}  // end namespace multimat
}  // end namespace axom
