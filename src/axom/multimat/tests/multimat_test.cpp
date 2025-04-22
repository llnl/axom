// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file multimat_test.cpp
 *
 * \brief Unit tests for MultiMat class
 */

#include "gtest/gtest.h"
#include <map>

#include "axom/multimat/multimat.hpp"

#include "axom/slic.hpp"

using namespace axom::multimat;

TEST(multimat, construct_empty_multimat_obj)
{
  MultiMat mm;

  EXPECT_TRUE(mm.isValid(true));
}

/* A structure to create test data for the MultiMat class. */
template <typename DataType>
struct MM_test_data
{
  std::vector<bool> fillBool_cellcen;
  std::vector<bool> fillBool_matcen;

  std::vector<int> matcount;

  axom::Array<double> volfrac_cellcen_dense;
  axom::Array<double> volfrac_matcen_dense;
  axom::Array<double> volfrac_cellcen_sparse;
  axom::Array<double> volfrac_matcen_sparse;

  axom::Array<DataType> cellmat_dense_arr;
  axom::Array<DataType> cellmat_sparse_arr;
  axom::Array<DataType> matcell_dense_arr;
  axom::Array<DataType> matcell_sparse_arr;

  std::vector<int> cellmat_beginvecs;
  std::vector<int> matcell_beginvecs;

  int num_cells, num_mats, stride;
  int nfilled = 0;

  DataType get_val(int c, int m, int s)
  {
    return static_cast<DataType>(c * 1000.0 + m * 1.0 + s * 0.01);
  }

  //constructor
  MM_test_data(int n_c = 20, int n_m = 10, int str = 1)
  {
    num_cells = n_c;
    num_mats = n_m;
    stride = str;

    nfilled = 0;
    fillBool_cellcen.resize(num_mats * num_cells, false);
    fillBool_matcen.resize(num_mats * num_cells, false);

    matcount.resize(num_cells, 0);
    for(int c = 0; c < num_cells; ++c)
    {
      for(int m = 0; m < num_mats; ++m)
      {
        int dense_cellcen_idx = c * num_mats + m;
        if(dense_cellcen_idx % 3 == 1)
        {
          int dense_matcen_idx = m * num_cells + c;
          fillBool_cellcen[dense_cellcen_idx] = true;
          fillBool_matcen[dense_matcen_idx] = true;
          nfilled++;
          matcount[c]++;
        }
      }
    }

    //create volfrac array
    int sparse_idx = 0;
    volfrac_cellcen_dense.resize(num_cells * num_mats, 0.0);
    volfrac_matcen_dense.resize(num_cells * num_mats, 0.0);
    volfrac_cellcen_sparse.resize(nfilled, 0.0);
    volfrac_matcen_sparse.resize(nfilled, 0.0);
    cellmat_beginvecs.resize(num_cells);
    matcell_beginvecs.resize(num_mats);
    for(auto i = 0; i < num_cells; ++i)
    {
      cellmat_beginvecs[i] = sparse_idx;
      for(auto m = 0; m < num_mats; ++m)
      {
        if(fillBool_cellcen[i * num_mats + m])
        {
          volfrac_cellcen_dense[i * num_mats + m] = 1.0 / (double)matcount[i];
          volfrac_matcen_dense[m * num_cells + i] = 1.0 / (double)matcount[i];
          volfrac_cellcen_sparse[sparse_idx++] = 1.0 / (double)matcount[i];
        }
      }
    }
    sparse_idx = 0;
    for(int m = 0; m < num_mats; ++m)
    {
      matcell_beginvecs[m] = sparse_idx;
      for(int c = 0; c < num_cells; ++c)
      {
        int dense_cellcen_idx = c * num_mats + m;
        if(fillBool_cellcen[dense_cellcen_idx])
        {
          volfrac_matcen_sparse[sparse_idx++] = 1.0 / (double)matcount[c];
        }
      }
    }

    //fill in the data
    cellmat_dense_arr.resize(num_mats * num_cells * stride);
    cellmat_sparse_arr.resize(nfilled * stride);
    matcell_dense_arr.resize(num_mats * num_cells * stride);
    matcell_sparse_arr.resize(nfilled * stride);
    sparse_idx = 0;
    for(int c = 0; c < num_cells; ++c)
    {
      for(int m = 0; m < num_mats; ++m)
      {
        int dense_cellcen_idx = c * num_mats + m;
        int dense_matcen_idx = m * num_cells + c;

        if(fillBool_cellcen[dense_cellcen_idx])
        {
          for(int s = 0; s < stride; ++s)
          {
            DataType val = get_val(c, m, s);
            cellmat_sparse_arr[sparse_idx * stride + s] = val;
            cellmat_dense_arr[dense_cellcen_idx * stride + s] = val;
          }
          ++sparse_idx;
        }
        if(fillBool_cellcen[dense_cellcen_idx])
        {
          for(int s = 0; s < stride; ++s)
          {
            DataType val = get_val(c, m, s);
            matcell_dense_arr[dense_matcen_idx * stride + s] = val;
          }
        }
      }
    }
    sparse_idx = 0;
    for(int m = 0; m < num_mats; ++m)
    {
      for(int c = 0; c < num_cells; ++c)
      {
        int dense_cellcen_idx = c * num_mats + m;
        if(fillBool_cellcen[dense_cellcen_idx])
        {
          for(int s = 0; s < stride; ++s)
          {
            DataType val = cellmat_dense_arr[dense_cellcen_idx * stride + s];
            matcell_sparse_arr[sparse_idx * stride + s] = val;
          }
          ++sparse_idx;
        }
      }
    }

  }  //end constructor

  //if val == 0, remove the entry. otherwise set entry using this volfrac
  //and the get_val function
  void setVal(int ci, int mi, double volfrac_val)
  {
    int cellcen_idx = ci * num_mats + mi;
    int matcen_idx = mi * num_cells + ci;

    bool already_filled = fillBool_cellcen[cellcen_idx];
    bool to_be_filled = volfrac_val != 0.0;

    if(!already_filled && !to_be_filled)
    {
      return;
    }

    volfrac_cellcen_dense[cellcen_idx] = volfrac_val;
    volfrac_matcen_dense[matcen_idx] = volfrac_val;

    if(already_filled && !to_be_filled)
    {
      matcount[ci] -= 1;
    }
    else if(!already_filled && to_be_filled)
    {
      matcount[ci] += 1;
    }

    for(int si = 0; si < stride; si++)
    {
      if(to_be_filled)
      {
        cellmat_dense_arr[cellcen_idx * stride + si] = get_val(ci, mi, si);
        matcell_dense_arr[matcen_idx * stride + si] = get_val(ci, mi, si);
      }
      else
      {
        cellmat_dense_arr[cellcen_idx * stride + si] = 0;
        matcell_dense_arr[matcen_idx * stride + si] = 0;
      }
    }

    int sparseidx = 0;
    for(int i = 0; i < mi; i++)
    {
      if(fillBool_cellcen[ci * num_mats + i])
      {
        sparseidx += 1;
      }
    }
    if(already_filled && to_be_filled)
    {  //just update
      volfrac_cellcen_sparse[cellmat_beginvecs[ci] + sparseidx] = volfrac_val;
      for(int si = 0; si < stride; si++)
      {
        int idx = (cellmat_beginvecs[ci] + sparseidx) * stride + si;
        cellmat_sparse_arr[idx] = get_val(ci, mi, si);
      }
    }
    else if(already_filled && !to_be_filled)
    {  //remove the entry
      volfrac_cellcen_sparse.erase(volfrac_cellcen_sparse.begin() + cellmat_beginvecs[ci] + sparseidx);
      for(int si = 0; si < stride; si++)
      {
        cellmat_sparse_arr.erase(cellmat_sparse_arr.begin() +
                                 (cellmat_beginvecs[ci] + sparseidx) * stride);
      }
      for(int i = ci + 1; i < num_cells; i++)
      {
        cellmat_beginvecs[i] -= 1;
      }
    }
    else if(!already_filled && to_be_filled)
    {  //adding an entry
      volfrac_cellcen_sparse.insert(
        volfrac_cellcen_sparse.begin() + (cellmat_beginvecs[ci] + sparseidx),
        volfrac_val);
      for(int si = 0; si < stride; si++)
      {
        cellmat_sparse_arr.insert(
          cellmat_sparse_arr.begin() + (cellmat_beginvecs[ci] + sparseidx) * stride + si,
          get_val(ci, mi, si));
      }
      for(int i = ci + 1; i < num_cells; i++)
      {
        cellmat_beginvecs[i] += 1;
      }
    }
    sparseidx = 0;
    for(int i = 0; i < ci; i++)
    {
      if(fillBool_matcen[mi * num_cells + i])
      {
        sparseidx += 1;
      }
    }
    if(already_filled && to_be_filled)
    {  //just update
      volfrac_matcen_sparse[matcell_beginvecs[mi] + sparseidx] = volfrac_val;
      for(int si = 0; si < stride; si++)
      {
        matcell_sparse_arr[(matcell_beginvecs[mi] + sparseidx) * stride + si] = get_val(ci, mi, si);
      }
    }
    else if(already_filled && !to_be_filled)
    {  //remove the entry
      volfrac_matcen_sparse.erase(volfrac_matcen_sparse.begin() + matcell_beginvecs[mi] + sparseidx);
      for(int si = 0; si < stride; si++)
      {
        matcell_sparse_arr.erase(matcell_sparse_arr.begin() +
                                 (matcell_beginvecs[mi] + sparseidx) * stride);
      }
      for(int i = mi + 1; i < num_mats; i++)
      {
        matcell_beginvecs[i] -= 1;
      }
    }
    else if(!already_filled && to_be_filled)
    {  //adding an entry
      volfrac_matcen_sparse.insert(volfrac_matcen_sparse.begin() + matcell_beginvecs[mi] + sparseidx,
                                   volfrac_val);
      for(int si = 0; si < stride; si++)
      {
        matcell_sparse_arr.insert(
          matcell_sparse_arr.begin() + (matcell_beginvecs[mi] + sparseidx) * stride + si,
          get_val(ci, mi, si));
      }
      for(int i = mi + 1; i < num_mats; i++)
      {
        matcell_beginvecs[i] += 1;
      }
    }

    fillBool_cellcen[cellcen_idx] = to_be_filled;
    fillBool_matcen[matcen_idx] = to_be_filled;
  }
};

/* Check the value using MultiMat findValue function */
template <typename DataType>
void check_values(MultiMat& mm, std::string arr_name, MM_test_data<DataType>& data)
{
  auto map_i = mm.getFieldIdx(arr_name);
  MultiMat::Field2D<DataType> map = mm.get2dField<DataType>(arr_name);
  EXPECT_TRUE(map.isValid());
  EXPECT_EQ(data.stride, map.stride());

  //check via findValue(...)
  for(int ci = 0; ci < data.num_cells; ++ci)
  {
    for(int mi = 0; mi < data.num_mats; ++mi)
    {
      for(int s = 0; s < data.stride; ++s)
      {
        DataType* d;
        if(mm.getFieldDataLayout(map_i) == DataLayout::CELL_DOM)
        {
          d = map.findValue(ci, mi, s);
        }
        else
        {
          d = map.findValue(mi, ci, s);
        }
        int dense_idx = ci * data.num_mats + mi;

        if(mm.getFieldSparsityLayout(map_i) == SparsityLayout::DENSE)
        {
          if(*d != data.cellmat_dense_arr[dense_idx * data.stride + s])
          {
            EXPECT_EQ(*d, data.cellmat_dense_arr[dense_idx * data.stride + s]);
          }
        }
        else
        {
          if(data.fillBool_cellcen[dense_idx])
          {
            EXPECT_NE(d, nullptr);
            EXPECT_EQ(*d, data.cellmat_dense_arr[dense_idx * data.stride + s]);
          }
          else
          {
            EXPECT_EQ(d, nullptr);
          }
        }
      }
    }
  }
}

/* Create a new multimat object of type T for testing. */
template <typename T>
MultiMat* newMM(MM_test_data<T>& data,
                DataLayout layout_used,
                SparsityLayout sparsity_used,
                std::string& array_name)
{
  MultiMat* mm_ptr = new MultiMat(layout_used, sparsity_used);
  MultiMat& mm = *mm_ptr;
  mm.setNumberOfCells(data.num_cells);
  mm.setNumberOfMaterials(data.num_mats);

  if(layout_used == DataLayout::CELL_DOM)
  {
    mm.setCellMatRel(data.fillBool_cellcen, DataLayout::CELL_DOM);
  }
  else
  {
    mm.setCellMatRel(data.fillBool_matcen,
                     DataLayout::MAT_DOM);  //todo change to MatCellRel
  }

  if(layout_used == DataLayout::CELL_DOM)
  {
    if(sparsity_used == SparsityLayout::DENSE)
    {
      mm.setVolfracField(data.volfrac_cellcen_dense, DataLayout::CELL_DOM, SparsityLayout::DENSE);
    }
    else
    {
      mm.setVolfracField(data.volfrac_cellcen_sparse, DataLayout::CELL_DOM, SparsityLayout::SPARSE);
    }
  }
  else
  {
    if(sparsity_used == SparsityLayout::DENSE)
    {
      mm.setVolfracField(data.volfrac_matcen_dense, DataLayout::MAT_DOM, SparsityLayout::DENSE);
    }
    else
    {
      mm.setVolfracField(data.volfrac_matcen_sparse, DataLayout::MAT_DOM, SparsityLayout::SPARSE);
    }
  }

  EXPECT_TRUE(mm.isValid());

  if(layout_used == DataLayout::CELL_DOM)
  {
    if(sparsity_used == SparsityLayout::DENSE)
    {
      mm.addField(array_name,
                  FieldMapping::PER_CELL_MAT,
                  layout_used,
                  sparsity_used,
                  data.cellmat_dense_arr.view(),
                  data.stride);
    }
    else
    {
      mm.addField(array_name,
                  FieldMapping::PER_CELL_MAT,
                  layout_used,
                  sparsity_used,
                  data.cellmat_sparse_arr.view(),
                  data.stride);
    }
  }
  else
  {
    if(sparsity_used == SparsityLayout::DENSE)
    {
      mm.addField(array_name,
                  FieldMapping::PER_CELL_MAT,
                  layout_used,
                  sparsity_used,
                  data.matcell_dense_arr.view(),
                  data.stride);
    }
    else
    {
      mm.addField(array_name,
                  FieldMapping::PER_CELL_MAT,
                  layout_used,
                  sparsity_used,
                  data.matcell_sparse_arr.view(),
                  data.stride);
    }
  }

  return mm_ptr;
}

//This doesn't test the MultiMat class, it tests the MM_test_data class
//to make sure the setValue function works correctly, so it can be used
TEST(multimat, test_data_setval)
{
  const int num_cells = 20;
  const int num_mats = 10;
  const int stride_val = 1;
  MM_test_data<double> data(num_cells, num_mats, stride_val);

  {  //remove everything, add everything back, and check
    auto data_cpy = data;

    for(int ci = 0; ci < data_cpy.num_cells; ci++)
    {
      for(int mi = 0; mi < data_cpy.num_mats; mi++)
      {
        data_cpy.setVal(ci, mi, 0);
      }
    }

    for(int ci = 0; ci < data.num_cells; ci++)
    {
      for(int mi = 0; mi < data.num_mats; mi++)
      {
        int dense_idx = ci * num_mats + mi;
        data_cpy.setVal(ci, mi, data.volfrac_cellcen_dense[dense_idx]);
      }
    }

    SLIC_ASSERT(data_cpy.fillBool_cellcen == data.fillBool_cellcen);
    SLIC_ASSERT(data_cpy.fillBool_matcen == data.fillBool_matcen);
    SLIC_ASSERT(data_cpy.matcount == data.matcount);
    SLIC_ASSERT(data_cpy.volfrac_cellcen_dense == data.volfrac_cellcen_dense);
    SLIC_ASSERT(data_cpy.volfrac_matcen_dense == data.volfrac_matcen_dense);
    SLIC_ASSERT(data_cpy.volfrac_cellcen_sparse == data.volfrac_cellcen_sparse);
    SLIC_ASSERT(data_cpy.volfrac_matcen_sparse == data.volfrac_matcen_sparse);
    SLIC_ASSERT(data_cpy.cellmat_dense_arr == data.cellmat_dense_arr);
    SLIC_ASSERT(data_cpy.cellmat_sparse_arr == data.cellmat_sparse_arr);
    SLIC_ASSERT(data_cpy.matcell_dense_arr == data.matcell_dense_arr);
    SLIC_ASSERT(data_cpy.matcell_sparse_arr == data.matcell_sparse_arr);
    SLIC_ASSERT(data_cpy.cellmat_beginvecs == data.cellmat_beginvecs);
    SLIC_ASSERT(data_cpy.matcell_beginvecs == data.matcell_beginvecs);
  }
}

template <typename DataType, int Stride>
void test_multimat_add_remove(std::pair<DataLayout, SparsityLayout> layout)
{
  const int num_cells = 20;
  const int num_mats = 10;
  const int stride_val = Stride;
  SLIC_INFO("--------------------------------------------------\n"
            << " Testing Multimat construction \n"
            << axom::fmt::format("Cells: {} Mats: {} Data type: {} Stride: {}",
                                 num_cells,
                                 num_mats,
                                 typeid(DataType).name(),
                                 Stride)
            << "--------------------------------------------------");
  MM_test_data<DataType> data(num_cells, num_mats, stride_val);

  const std::map<DataLayout, std::vector<bool>> relationMap {
    {DataLayout::CELL_DOM, data.fillBool_cellcen},
    {DataLayout::MAT_DOM, data.fillBool_matcen}};

  const std::map<std::pair<DataLayout, SparsityLayout>, axom::Array<double>> volFracMap {
    {{DataLayout::CELL_DOM, SparsityLayout::DENSE}, data.volfrac_cellcen_dense},
    {{DataLayout::CELL_DOM, SparsityLayout::SPARSE}, data.volfrac_cellcen_sparse},
    {{DataLayout::MAT_DOM, SparsityLayout::DENSE}, data.volfrac_matcen_dense},
    {{DataLayout::MAT_DOM, SparsityLayout::SPARSE}, data.volfrac_matcen_sparse},
  };

  const std::map<std::pair<DataLayout, SparsityLayout>, axom::Array<DataType>> fieldDataMap {
    {{DataLayout::CELL_DOM, SparsityLayout::DENSE}, data.cellmat_dense_arr},
    {{DataLayout::CELL_DOM, SparsityLayout::SPARSE}, data.cellmat_sparse_arr},
    {{DataLayout::MAT_DOM, SparsityLayout::DENSE}, data.matcell_dense_arr},
    {{DataLayout::MAT_DOM, SparsityLayout::SPARSE}, data.matcell_sparse_arr},
  };

  DataLayout layout_used = layout.first;
  SparsityLayout sparsity_used = layout.second;

  // Construct a Multimat object with a given layout
  MultiMat mm(layout.first, layout.second);
  mm.setNumberOfCells(data.num_cells);
  mm.setNumberOfMaterials(data.num_mats);
  mm.setCellMatRel(relationMap.find(layout_used)->second, layout_used);
  mm.setVolfracField(volFracMap.find(layout)->second.view(), layout_used, sparsity_used);

  SLIC_INFO(axom::fmt::format("Constructing Multimat object with layout {}/{}",
                              mm.getFieldDataLayoutAsString(0),
                              mm.getFieldSparsityLayoutAsString(0)));

  EXPECT_TRUE(mm.isValid(true));
  EXPECT_EQ(mm.getNumberOfCells(), data.num_cells);
  EXPECT_EQ(mm.getNumberOfMaterials(), data.num_mats);
  EXPECT_EQ(mm.getNumberOfFields(), 1);

  EXPECT_EQ(mm.getFieldIdx("OwnedField"), -1);

  // Add an owned field to the Multimat instance.
  axom::Array<DataType> owned_array = fieldDataMap.find(layout)->second;
  int field_idx = mm.addField("OwnedField",
                              FieldMapping::PER_CELL_MAT,
                              layout_used,
                              sparsity_used,
                              owned_array.view(),
                              data.stride);

  EXPECT_EQ(mm.getFieldDataLayout(field_idx), layout_used);
  EXPECT_EQ(mm.getFieldSparsityLayout(field_idx), sparsity_used);
  EXPECT_EQ(mm.getFieldName(field_idx), "OwnedField");
  EXPECT_EQ(mm.getFieldIdx("OwnedField"), field_idx);
  EXPECT_EQ(mm.getNumberOfFields(), 2);

  // Check field returned from the multimat instance.
  MMField2D<DataType> field = mm.get2dField<DataType>("OwnedField");
  EXPECT_EQ(field.isDense(), sparsity_used == SparsityLayout::DENSE);
  EXPECT_EQ(field.isCellDom(), layout_used == DataLayout::CELL_DOM);
  EXPECT_EQ(field.stride(), stride_val);

  // Check dense/sparse-typed fields as well.
  if(sparsity_used == SparsityLayout::DENSE)
  {
    MultiMat::DenseField2D<DataType> field_dense = mm.getDense2dField<DataType>("OwnedField");
    EXPECT_EQ(field_dense.isDense(), sparsity_used == SparsityLayout::DENSE);
    EXPECT_EQ(field_dense.isCellDom(), layout_used == DataLayout::CELL_DOM);
    EXPECT_EQ(field_dense.stride(), stride_val);
  }
  else
  {
    MultiMat::SparseField2D<DataType> field_sparse = mm.getSparse2dField<DataType>("OwnedField");
    EXPECT_EQ(field_sparse.isDense(), sparsity_used == SparsityLayout::DENSE);
    EXPECT_EQ(field_sparse.isCellDom(), layout_used == DataLayout::CELL_DOM);
    EXPECT_EQ(field_sparse.stride(), stride_val);
  }

  check_values<DataType>(mm, "OwnedField", data);

  // Remove the field.
  mm.removeField("OwnedField");

  EXPECT_EQ(mm.getFieldName(field_idx), "");
  EXPECT_EQ(mm.getFieldIdx("OwnedField"), -1);
  EXPECT_EQ(mm.getNumberOfFields(), 1);
}

TEST(multimat, construct_multimat_1_array)
{
  std::vector<std::pair<DataLayout, SparsityLayout>> layouts {
    {DataLayout::CELL_DOM, SparsityLayout::DENSE},
    {DataLayout::CELL_DOM, SparsityLayout::SPARSE},
    {DataLayout::MAT_DOM, SparsityLayout::DENSE},
    {DataLayout::MAT_DOM, SparsityLayout::SPARSE},
  };

  for(auto layout : layouts)
  {
    test_multimat_add_remove<int, 1>(layout);
    test_multimat_add_remove<float, 1>(layout);
    test_multimat_add_remove<double, 1>(layout);
  }
}

TEST(multimat, construct_multimat_3_array)
{
  std::vector<std::pair<DataLayout, SparsityLayout>> layouts {
    {DataLayout::CELL_DOM, SparsityLayout::DENSE},
    {DataLayout::CELL_DOM, SparsityLayout::SPARSE},
    {DataLayout::MAT_DOM, SparsityLayout::DENSE},
    {DataLayout::MAT_DOM, SparsityLayout::SPARSE},
  };

  for(auto layout : layouts)
  {
    test_multimat_add_remove<int, 3>(layout);
    test_multimat_add_remove<float, 3>(layout);
    test_multimat_add_remove<double, 3>(layout);
  }
}

template <typename DataType, int Stride>
void test_multimat_conversion(std::pair<DataLayout, SparsityLayout> from,
                              std::pair<DataLayout, SparsityLayout> to,
                              int allocatorID = axom::getDefaultAllocatorID())
{
  const int num_cells = 20;
  const int num_mats = 10;
  const int stride_val = Stride;
  SLIC_INFO("\n--------------------------------------------------\n"
            << " Testing Multimat construction and conversion \n"
            << axom::fmt::format(" Cells: {} Mats: {} Data type: {} Stride: {}\n",
                                 num_cells,
                                 num_mats,
                                 std::string(typeid(DataType).name()),
                                 Stride)
            << "--------------------------------------------------\n");
  MM_test_data<DataType> data(num_cells, num_mats, stride_val);

  const std::map<DataLayout, std::vector<bool>> relationMap {
    {DataLayout::CELL_DOM, data.fillBool_cellcen},
    {DataLayout::MAT_DOM, data.fillBool_matcen}};

  const std::map<std::pair<DataLayout, SparsityLayout>, axom::Array<double>> volFracMap {
    {{DataLayout::CELL_DOM, SparsityLayout::DENSE}, data.volfrac_cellcen_dense},
    {{DataLayout::CELL_DOM, SparsityLayout::SPARSE}, data.volfrac_cellcen_sparse},
    {{DataLayout::MAT_DOM, SparsityLayout::DENSE}, data.volfrac_matcen_dense},
    {{DataLayout::MAT_DOM, SparsityLayout::SPARSE}, data.volfrac_matcen_sparse},
  };

  const std::map<std::pair<DataLayout, SparsityLayout>, axom::Array<DataType>> fieldDataMap {
    {{DataLayout::CELL_DOM, SparsityLayout::DENSE}, data.cellmat_dense_arr},
    {{DataLayout::CELL_DOM, SparsityLayout::SPARSE}, data.cellmat_sparse_arr},
    {{DataLayout::MAT_DOM, SparsityLayout::DENSE}, data.matcell_dense_arr},
    {{DataLayout::MAT_DOM, SparsityLayout::SPARSE}, data.matcell_sparse_arr},
  };

  DataLayout layout_used = from.first;
  SparsityLayout sparsity_used = from.second;

  // Construct a Multimat object with a given layout
  MultiMat mm(from.first, from.second);
  mm.setNumberOfCells(data.num_cells);
  mm.setNumberOfMaterials(data.num_mats);
  mm.setCellMatRel(relationMap.find(layout_used)->second, layout_used);
  mm.setVolfracField(volFracMap.find(from)->second.view(), layout_used, sparsity_used);
  mm.setAllocatorID(allocatorID);

  SLIC_INFO(axom::fmt::format("Constructing Multimat object with layout {}/{}",
                              mm.getFieldDataLayoutAsString(0),
                              mm.getFieldSparsityLayoutAsString(0)));

  EXPECT_TRUE(mm.isValid(true));
  EXPECT_EQ(mm.getNumberOfCells(), data.num_cells);
  EXPECT_EQ(mm.getNumberOfMaterials(), data.num_mats);

  // Add an owned field and an unowned field to the Multimat instance.
  axom::Array<DataType> owned_array = fieldDataMap.find(from)->second;
  int int_field_idx = mm.addField("OwnedField",
                                  FieldMapping::PER_CELL_MAT,
                                  layout_used,
                                  sparsity_used,
                                  owned_array.view(),
                                  data.stride);

  int ext_field_idx = mm.addExternalField("UnownedField",
                                          FieldMapping::PER_CELL_MAT,
                                          layout_used,
                                          sparsity_used,
                                          owned_array.view(),
                                          data.stride);

  check_values<DataType>(mm, "OwnedField", data);
  check_values<DataType>(mm, "UnownedField", data);

  DataLayout to_layout = to.first;
  SparsityLayout to_sparsity = to.second;

  mm.convertLayout(to_layout, to_sparsity);

  SLIC_INFO(axom::fmt::format("Converted multimat instance to layout {}/{}",
                              mm.getFieldDataLayoutAsString(0),
                              mm.getFieldSparsityLayoutAsString(0)));

  // Check owned fields for conversion.
  EXPECT_EQ(mm.getFieldDataLayout(0), to_layout);
  EXPECT_EQ(mm.getFieldSparsityLayout(0), to_sparsity);
  EXPECT_EQ(mm.getFieldDataLayout(int_field_idx), to_layout);
  EXPECT_EQ(mm.getFieldSparsityLayout(int_field_idx), to_sparsity);

  // Special case: externally-owned fields can only be transposed between
  // material- and cell-dominant layouts; they cannot be converted between
  // sparse and dense layouts.
  EXPECT_EQ(mm.getFieldDataLayout(ext_field_idx), to_layout);
  EXPECT_EQ(mm.getFieldSparsityLayout(ext_field_idx), sparsity_used);

  check_values<DataType>(mm, "OwnedField", data);
  check_values<DataType>(mm, "UnownedField", data);
}

const std::vector<std::pair<DataLayout, SparsityLayout>> g_test_layouts {
  {DataLayout::CELL_DOM, SparsityLayout::DENSE},
  {DataLayout::CELL_DOM, SparsityLayout::SPARSE},
  {DataLayout::MAT_DOM, SparsityLayout::DENSE},
  {DataLayout::MAT_DOM, SparsityLayout::SPARSE},
};

TEST(multimat, convert_multimat_1_array)
{
  for(auto layout_from : g_test_layouts)
  {
    for(auto layout_to : g_test_layouts)
    {
      test_multimat_conversion<int, 1>(layout_from, layout_to);
      test_multimat_conversion<float, 1>(layout_from, layout_to);
      test_multimat_conversion<double, 1>(layout_from, layout_to);
    }
  }
}

TEST(multimat, convert_multimat_3_array)
{
  for(auto layout_from : g_test_layouts)
  {
    for(auto layout_to : g_test_layouts)
    {
      test_multimat_conversion<int, 3>(layout_from, layout_to);
      test_multimat_conversion<float, 3>(layout_from, layout_to);
      test_multimat_conversion<double, 3>(layout_from, layout_to);
    }
  }
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_CUDA) || defined(AXOM_USE_HIP)
TEST(multimat, convert_multimat_1_array_gpu)
{
  int unifiedAllocID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);

  for(auto layout_from : g_test_layouts)
  {
    for(auto layout_to : g_test_layouts)
    {
      test_multimat_conversion<int, 1>(layout_from, layout_to, unifiedAllocID);
      test_multimat_conversion<float, 1>(layout_from, layout_to, unifiedAllocID);
      test_multimat_conversion<double, 1>(layout_from, layout_to, unifiedAllocID);
    }
  }
}

TEST(multimat, convert_multimat_3_array_gpu)
{
  int unifiedAllocID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);

  for(auto layout_from : g_test_layouts)
  {
    for(auto layout_to : g_test_layouts)
    {
      test_multimat_conversion<int, 3>(layout_from, layout_to, unifiedAllocID);
      test_multimat_conversion<float, 3>(layout_from, layout_to, unifiedAllocID);
      test_multimat_conversion<double, 3>(layout_from, layout_to, unifiedAllocID);
    }
  }
}
  #endif  // defined(AXOM_USE_CUDA) || defined(AXOM_USE_HIP)
#endif    // defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

/* Test dynamic mode */
TEST(multimat, test_dynamic_multimat_1_array)
{
  const int num_cells = 20;
  const int num_mats = 10;
  const int stride_val = 1;

  std::vector<DataLayout> data_layouts = {DataLayout::CELL_DOM, DataLayout::MAT_DOM};
  std::vector<SparsityLayout> sparsity_layouts = {SparsityLayout::DENSE, SparsityLayout::SPARSE};

  std::string array_name = "Array 1";

  for(auto layout_used : data_layouts)
  {
    for(auto sparsity_used : sparsity_layouts)
    {
      MM_test_data<double> data(num_cells, num_mats, stride_val);

      SLIC_INFO("--------------------\nConstructing MultiMat object...");
      MultiMat* mm_ptr = newMM(data, layout_used, sparsity_used, array_name);
      MultiMat& mm = *mm_ptr;
      SLIC_INFO("Layout: " << mm.getFieldDataLayoutAsString(0) << " and "
                           << mm.getFieldSparsityLayoutAsString(0));

      EXPECT_TRUE(mm.isValid(true));
      EXPECT_EQ(mm.getNumberOfCells(), data.num_cells);
      EXPECT_EQ(mm.getNumberOfMaterials(), data.num_mats);

      check_values<double>(mm, array_name, data);

      SLIC_INFO("Convert to dynamic...");
      mm.convertToDynamic();

      auto volfrac_field = mm.getVolfracField();
      auto arr = mm.get2dField<double>(array_name);

      SLIC_INFO("Removing every other entries...");
      for(int ci = 0; ci < data.num_cells; ci++)
      {
        for(int mi = 0; mi < data.num_mats; mi++)
        {
          int cell_centric_dense_idx = ci * data.num_mats + mi;
          if(data.fillBool_cellcen[cell_centric_dense_idx])
          {
            data.setVal(ci, mi, 0);
            int idx1, idx2;
            if(layout_used == DataLayout::CELL_DOM)
            {
              idx1 = ci;
              idx2 = mi;
            }
            else
            {
              idx1 = mi;
              idx2 = ci;
            }

            EXPECT_TRUE(mm.removeEntry(ci, mi));
            volfrac_field(idx1, idx2) = 0.0;
            for(int s = 0; s < stride_val; s++)
            {
              arr(idx1, idx2, s) = 0.0;
            }
          }
        }
      }

      SLIC_INFO("Add material 0 into all cells");
      for(int ci = 0; ci < data.num_cells; ci++)
      {
        data.setVal(ci, 0, 1.0);
        int idx1, idx2;
        if(layout_used == DataLayout::CELL_DOM)
        {
          idx1 = ci;
          idx2 = 0;
        }
        else
        {
          idx1 = 0;
          idx2 = ci;
        }

        EXPECT_TRUE(mm.addEntry(ci, 0));
        volfrac_field(idx1, idx2) = 1.0;
        for(int s = 0; s < stride_val; s++)
        {
          arr(idx1, idx2, s) = data.cellmat_dense_arr[ci * data.num_mats * stride_val + s];
        }
      }

      EXPECT_TRUE(mm.isValid(true));
      EXPECT_EQ(mm.getNumberOfCells(), data.num_cells);
      EXPECT_EQ(mm.getNumberOfMaterials(), data.num_mats);
      for(int i = 0; i < mm.getNumberOfFields(); ++i)
      {
        //For now, during dynamic mode, layout is in dense
        EXPECT_EQ(mm.getFieldDataLayout(i), layout_used);
        EXPECT_EQ(mm.getFieldSparsityLayout(i), SparsityLayout::DENSE);
      }

      check_values<double>(mm, array_name, data);

      mm.convertToStatic();

      EXPECT_TRUE(mm.isValid(true));
      EXPECT_EQ(mm.getNumberOfCells(), data.num_cells);
      EXPECT_EQ(mm.getNumberOfMaterials(), data.num_mats);
      for(int i = 0; i < mm.getNumberOfFields(); ++i)
      {
        EXPECT_EQ(mm.getFieldDataLayout(i), layout_used);
        EXPECT_EQ(mm.getFieldSparsityLayout(i), sparsity_used);
      }
      check_values<double>(mm, array_name, data);

      delete mm_ptr;
    }
  }
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
