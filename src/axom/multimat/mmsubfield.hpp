// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MMSUBFIELD_H_
#define MMSUBFIELD_H_

#include "axom/multimat/multimat.hpp"

namespace axom
{
namespace multimat
{
/**
 * Class for MultiMat 2D SubFields
 */

template <typename Field2DType>
class MMSubField2D
  : public slam::SubMap<typename Field2DType::BiVarMapType, MultiMat::RangeSetType, slam::policies::ConcreteInterface>
{
public:
  using SubSetType = MultiMat::RangeSetType;
  using SubMapType =
    slam::SubMap<typename Field2DType::BiVarMapType, SubSetType, slam::policies::ConcreteInterface>;
  using SuperMapType = typename Field2DType::BiVarMapType;
  using BiVarSetType = typename Field2DType::BiVarSetType;

  // Default Constructor
  MMSubField2D() : SubMapType(), m_superfield(nullptr), firstSetIndex(-1) { };

  // Constructor
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE MMSubField2D(Field2DType* superfield,
                                int firstIndex,
                                bool indicesHaveIndirection = true)
    //why is without pointer type casting cause a compilation error?
    : SubMapType((SuperMapType*)superfield,
                 superfield->set()->elementRangeSet(firstIndex),
                 indicesHaveIndirection)
    , m_superfield(superfield)
    , firstSetIndex(firstIndex)
  { }

  int getOuterIndex() { return firstSetIndex; }

private:
  Field2DType* m_superfield;
  int firstSetIndex;

};  //end class MMSubField2D

template <typename Field2DType, DataLayout DataLayoutT>
class MMSubField2DWrap : public MMSubField2D<Field2DType>
{ };

// specialization for Cell Dom
template <typename Field2DType>
class MMSubField2DWrap<Field2DType, DataLayout::CELL_DOM> : public MMSubField2D<Field2DType>
{
public:
  using SFB = MMSubField2D<Field2DType>;
  MMSubField2DWrap(Field2DType* superfield, int firstIndex, bool indirection = true)
    : SFB(superfield, firstIndex, indirection)
  { }

  DataLayout getDataLayout() { return DataLayout::CELL_DOM; }
  int cellId() { return this->getOuterIndex(); }
  int matId(int i) { return this->index(i); }
};

// specialization for Mat Dom
template <typename Field2DType>
class MMSubField2DWrap<Field2DType, DataLayout::MAT_DOM> : public MMSubField2D<Field2DType>
{
public:
  using SFB = MMSubField2D<Field2DType>;
  MMSubField2DWrap(Field2DType* superfield, int firstIndex, bool indirection = true)
    : SFB(superfield, firstIndex, indirection)
  { }

  DataLayout getDataLayout() { return DataLayout::MAT_DOM; }
  int matId() { return this->getOuterIndex(); }
  int cellId(int i) { return this->index(i); }
};

}  //end namespace multimat
}  //end namespace axom

///////////////////////////////////////////////////////////////////////////

#endif
