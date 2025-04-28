// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MMFIELD_H_
#define MMFIELD_H_

#include "axom/multimat/multimat.hpp"
#include "axom/multimat/mmsubfield.hpp"

namespace axom
{
namespace multimat
{
/**
 * Class for MultiMat 2D Fields
 * Ideally this class would have template specialization for each type of 
 * layout (dense/sparse, mat/cell dom, and maybe more). 
 */
template <typename DataType, typename BiSet = MultiMat::BivariateSetType>
class MMField2D : public MultiMat::BivariateMapType<DataType, BiSet>
{
public:
  using BiVarSetType = BiSet;
  using BiVarMapType = MultiMat::BivariateMapType<DataType, BiVarSetType>;
  using ProductSetType = MultiMat::ProductSetType;
  using RelationSetType = MultiMat::RelationSetType;

  using SetPosition = typename BiVarMapType::SetPosition;

  using Field2DType = MMField2D<DataType, BiVarSetType>;
  using SubFieldType = MMSubField2D<Field2DType>;
  using ConstSubFieldType = const MMSubField2D<const Field2DType>;

  //slam typedef
  using SubMapType = typename BiVarMapType::SubMapType;

  /**
  * \brief Constructor
  */
  MMField2D() = delete;
  // Constructor with specified layouts
  //MMField2D(MultiMat& mm, DataLayout data_layout, SparsityLayout sparsity_layout,
  //  const std::string& arr_name = "unnamed",
  //  const DataType* data_arr = nullptr,
  //  int stride = 1);
  // Constructor given the bivariate set
  MMField2D(const MultiMat& mm,
            const BiSet*,
            const int fieldIdx,
            axom::ArrayView<DataType> data_arr = {},
            int stride = 1);

  template <typename BiSetType = BiSet,
            typename Enable = std::enable_if_t<!std::is_abstract<BiSetType>::value>>
  MMField2D(const MultiMat& mm,
            const BiSetType&,
            const int fieldIdx,
            axom::ArrayView<DataType> data_arr = {},
            int stride = 1);

  bool operator==(const MMField2D& other) const
  {
    return ((m_mm == other.m_mm) && (this->set() == other.set()) &&
            (this->getMap()->data() == other.getMap()->data()));
  }

  using BiVarMapType::operator();  //why is this needed?

  //subfield (instead of SubMap)
  SubFieldType getSubfield(SetPosition firstIdx) { return operator()(firstIdx); }
  AXOM_HOST_DEVICE SubFieldType operator()(SetPosition firstIdx)
  {
    const bool hasInd = this->submapIndicesHaveIndirection();
    return SubFieldType(this, firstIdx, hasInd);
  }
  AXOM_HOST_DEVICE const ConstSubFieldType operator()(SetPosition firstIdx) const
  {
    const bool hasInd = this->submapIndicesHaveIndirection();
    return ConstSubFieldType(this, firstIdx, hasInd);
  }

  //Mimic BivariateMap operator(i) and return slam submap
  SubMapType getSlamSubMap(SetPosition firstIdx) { return BiVarMapType::operator()(firstIdx); }

  std::string getName() { return m_mm->getFieldName(m_fieldIdx); };

  MultiMat::IndexSet getSubfieldIndexingSet(int idx)
  {
    return m_mm->getSubfieldIndexingSet(idx, m_data_layout, m_sparsity_layout);
  }

  bool isDense() const { return m_sparsity_layout == SparsityLayout::DENSE; }
  bool isSparse() const { return m_sparsity_layout == SparsityLayout::SPARSE; }
  bool isCellDom() const { return m_data_layout == DataLayout::CELL_DOM; }
  bool isMatDom() const { return m_data_layout == DataLayout::MAT_DOM; }

private:
  const MultiMat* m_mm;

  DataLayout m_data_layout;
  SparsityLayout m_sparsity_layout;

  int m_fieldIdx;
};

//////////////////////// Implementation of MMField2D ////////////////////////////

// Constructor given layouts
//template<typename DataType, typename BiSet>
//inline MMField2D<DataType, BiSet>::MMField2D(MultiMat& mm,
//    DataLayout data_layout, SparsityLayout sparsity_layout,
//    const std::string& arr_name, const DataType* data_arr, int stride) :
//  MMField2D(mm, mm.get_mapped_biSet(data_layout, sparsity_layout),
//     arr_name, data_arr, stride)
//{ }

// Constructor given the biset
template <typename DataType, typename BiSet>
inline MMField2D<DataType, BiSet>::MMField2D(const MultiMat& mm,
                                             const BiSet* biset,
                                             const int fieldIdx,
                                             axom::ArrayView<DataType> data_arr,
                                             int stride)
  :  //call Bivariate map constructor
  BiVarMapType(biset, data_arr, stride)
  , m_mm(&mm)
  , m_fieldIdx(fieldIdx)
{
  SLIC_ASSERT(stride > 0);

  m_data_layout = mm.getFieldDataLayout(fieldIdx);
  m_sparsity_layout = mm.getFieldSparsityLayout(fieldIdx);
}

template <typename DataType, typename BiSet>
template <typename BiSetType, typename Enable>
MMField2D<DataType, BiSet>::MMField2D(const MultiMat& mm,
                                      const BiSetType& bisetValue,
                                      const int fieldIdx,
                                      axom::ArrayView<DataType> data_arr,
                                      int stride)
  : BiVarMapType(bisetValue, data_arr, stride)
  , m_mm(&mm)
  , m_fieldIdx(fieldIdx)
{
  SLIC_ASSERT(stride > 0);

  m_data_layout = mm.getFieldDataLayout(fieldIdx);
  m_sparsity_layout = mm.getFieldSparsityLayout(fieldIdx);
}

//////////////////////// MMField2D Templated ////////////////////////////

// Child class of MMField2D, typed with layout (cell/mat dom) and sparsity
template <typename DataType, DataLayout DataLayoutT, typename BiSet = MultiMat::BivariateSetType>
class MMField2DTemplated : public MMField2D<DataType, BiSet>
{
  using Field2DType = MMField2D<DataType, BiSet>;

public:
  MMField2DTemplated(MultiMat& mm, int fieldIdx, axom::ArrayView<DataType> data_arr = {}, int stride = 1)
    : Field2DType(mm, mm.getCompatibleBivarSet<BiSet>(fieldIdx), fieldIdx, data_arr, stride)
  { }
};

}  //end namespace multimat
}  //end namespace axom

#endif
