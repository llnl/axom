// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file traversal.cpp
 *
 * \brief Example of traversing data using various methods in MultiMat
 */

#include "axom/multimat/multimat.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <cstdlib>  // for std::rand(), RAND_MAX

using namespace axom::multimat;

double getRandomDouble(double low, double high)
{
  const double delta = high - low;
  const double c = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
  return (delta * c + low);
}

/**
 * \brief Show different ways to traverse the values inside MultiMat.
 *
 * The two fields used are "Cell Array" for 1d array
 * where there is an entry for each cell, and "CellMat Array" where
 * there is an entry for each cell and each material.
 *
 * Each traversal loops access the value and, for the sparse layouts, the index
 * of the material (assuming a cell-dominant layout).
 *
 * Note that some index access are commented out to avoid compiler warnings,
 * but are valid index accessing code.
 *
 */
void various_traversal_methods(int nmats, int ncells, int ncomp, bool use_sparse, double fill_percentage)
{
  axom::utilities::Timer timer;

  auto layout = DataLayout::CELL_DOM;
  auto sparsity = (use_sparse ? SparsityLayout::SPARSE : SparsityLayout::DENSE);

  MultiMat mm(layout, sparsity);

  int nfilled = 0;
  std::vector<bool> cellMatRel(nmats * ncells, false);
  for(int i = 0; i < nmats * ncells; i++)
  {
    if(getRandomDouble(0, 1) < fill_percentage)
    {
      cellMatRel[i] = true;
      nfilled++;
    }
  }

  mm.setNumberOfMaterials(nmats);
  mm.setNumberOfCells(ncells);
  mm.setCellMatRel(cellMatRel, layout);

  //create the std::vector data for the field arrays
  axom::Array<double> cell_arr(ncells * ncomp);
  double c_sum = 0;
  for(int i = 0; i < ncells; ++i)
  {
    for(int comp = 0; comp < ncomp; ++comp)
    {
      cell_arr[i * ncomp + comp] = (double)i * 2.0 + comp * 0.01;
      c_sum += cell_arr[i * ncomp + comp];
    }
  }

  axom::Array<double> cellmat_arr;
  cellmat_arr.resize((use_sparse ? nfilled : nmats * ncells) * ncomp);
  double x_sum = 0;
  for(axom::IndexType i = 0; i < cellmat_arr.size() / ncomp; i++)
  {
    if(use_sparse || cellMatRel[i])
    {
      for(int comp = 0; comp < ncomp; ++comp)
      {
        cellmat_arr[i * ncomp + comp] = (double)i * 1.1 + comp * 0.001;
        x_sum += cellmat_arr[i * ncomp + comp];
      }
    }
  }

  //create volfrac array
  axom::Array<double> volfrac_arr(ncells * nmats, 0);
  for(auto i = 0; i < ncells; ++i)
  {
    int matcount = 0;
    for(auto m = 0; m < nmats; ++m)
    {
      if(cellMatRel[i * nmats + m])
      {
        matcount += 1;
      }
    }
    for(auto m = 0; m < nmats; ++m)
    {
      if(cellMatRel[i * nmats + m])
      {
        volfrac_arr[i * nmats + m] = 1.0 / (double)matcount;
      }
    }
  }

  mm.setVolfracField(volfrac_arr, layout, SparsityLayout::DENSE);
  if(sparsity == SparsityLayout::SPARSE)
  {
    mm.convertFieldToSparse(0);
  }
  mm.addField("Cell Array", FieldMapping::PER_CELL, layout, sparsity, cell_arr.view(), ncomp);
  mm.addField("CellMat Array", FieldMapping::PER_CELL_MAT, layout, sparsity, cellmat_arr.view(), ncomp);

  double sum = 0;

  //Different accessing methods ...

  //get the volfrac field
  auto volfrac_map = mm.getVolfracField();
  auto volfrac_map2 = mm.get2dField<double>("Volfrac");
  SLIC_ASSERT(volfrac_map == volfrac_map2);
  //volfrac field is access the same way as a regular Field2d
  AXOM_UNUSED_VAR(volfrac_map);
  AXOM_UNUSED_VAR(volfrac_map2);

  // --------- returning SLAM map and submap -----------
  SLIC_INFO("\n -- Access from SLAM map (and submap) -- ");
  sum = 0;
  timer.reset();
  timer.start();
  {
    MultiMat::Field1D<double> map = mm.get1dField<double>("Cell Array");
    SLIC_ASSERT(ncomp == map.stride());
    SLIC_ASSERT(ncomp == map.numComp());
    for(int i = 0; i < mm.getNumberOfCells(); i++)
    {
      for(int comp = 0; comp < map.numComp(); comp++)
      {
        double val = map(i, comp);                          //<----
        SLIC_ASSERT(val == map[i * map.numComp() + comp]);  //1d bracket access
        sum += val;
      }
    }
  }
  timer.stop();
  SLIC_INFO("  Field1D: " << timer.elapsed() << " sec");
  SLIC_ASSERT(c_sum == sum);
  AXOM_UNUSED_VAR(c_sum);
  AXOM_UNUSED_VAR(sum);

  sum = 0;
  timer.reset();
  timer.start();
  {
    auto map2d = mm.get2dField<double>("CellMat Array");
    SLIC_ASSERT(ncomp == map2d.stride());
    SLIC_ASSERT(ncomp == map2d.numComp());
    for(int i = 0; i < map2d.firstSetSize(); i++)
    {
      const MultiMat::IdSet& rel_set = static_cast<const MultiMat::IdSet&>(map2d.indexSet(i));
      auto submap = map2d(i);
      SLIC_ASSERT(rel_set.size() == submap.size());
      for(int k = 0; k < submap.size(); k++)
      {
        int idx = submap.index(k);  //mat id
        int idx2 = rel_set[k];      //another way to get mat id
        SLIC_ASSERT(idx == idx2);
        AXOM_UNUSED_VAR(idx);
        AXOM_UNUSED_VAR(idx2);

        for(int c = 0; c < submap.numComp(); ++c)
        {
          double val = submap.value(k, c);                       //<----------
          SLIC_ASSERT(val == submap(k, c));                      //operator () access
          SLIC_ASSERT(val == submap[k * submap.numComp() + c]);  //bracket access
          sum += val;
        }
      }
    }
  }
  timer.stop();
  SLIC_INFO("  Field2D: " << timer.elapsed() << " sec");
  SLIC_ASSERT(x_sum == sum);
  AXOM_UNUSED_VAR(x_sum);

  // ------- Dense Access ----------
  SLIC_INFO("\n -- Dense Access via map-- ");
  sum = 0;
  timer.reset();
  timer.start();
  {
    // _multimat_using_fields_dense_start
    auto map = mm.get2dField<double>("CellMat Array");

    for(int i = 0; i < mm.getNumberOfCells(); i++)
    {
      for(int m = 0; m < mm.getNumberOfMaterials(); m++)
      {
        for(int c = 0; c < map.numComp(); ++c)
        {
          double* valptr = map.findValue(i, m, c);
          // ^ contains a hidden for-loop for sparse layouts, O(row_size) time
          if(valptr)
          {
            sum += *valptr;
          }
        }
      }
    }
    // _multimat_using_fields_dense_end
  }
  timer.stop();
  SLIC_INFO("  Field2D: " << timer.elapsed() << " sec");
  SLIC_ASSERT(x_sum == sum);

  if(sparsity == SparsityLayout::SPARSE)
  {
    // ------------ return index set --------------
    SLIC_INFO("\n -- Access by Map with indexing set-- ");
    sum = 0;
    timer.reset();
    timer.start();
    {
      //Note you can achieve the same thing by using
      // a Map pointer to point to a bivariateMap object
      //MultiMat::Field1D<double>& map = mm.get1dField<double>("CellMat Array");

      //_multimat_using_fields_indexset_start
      auto map = mm.get2dField<double>("CellMat Array");

      for(int i = 0; i < mm.getNumberOfCells(); i++)
      {
        //the materials (by id) in this cell
        MultiMat::IdSet setOfMaterialsInThisCell = mm.getMatInCell(i);
        //the indices into the maps
        MultiMat::IndexSet indexSet = mm.getIndexingSetOfCell(i, sparsity);

        SLIC_ASSERT(setOfMaterialsInThisCell.size() == indexSet.size());
        for(int j = 0; j < indexSet.size(); j++)
        {
          //int idx = setOfMaterialsInThisCell.at(j); //mat_id
          for(int comp = 0; comp < map.numComp(); ++comp)
          {
            double val = map[indexSet[j] * map.numComp() + comp];  //<-----

            //if 1dMap is used, this is also valid
            //SLIC_ASSERT(val == map(indexSet[j], comp));

            sum += val;
          }
        }
      }
      //_multimat_using_fields_indexset_end
    }
    timer.stop();
    SLIC_INFO("  Field2D: " << timer.elapsed() << " sec");
    SLIC_ASSERT(x_sum == sum);
  }

  // ---------- using iterator with Map and Submap -------------
  SLIC_INFO("\n -- With Map (and Submap) iterators -- ");
  sum = 0;
  timer.reset();
  timer.start();
  {
    auto map = mm.get1dField<double>("Cell Array");
    for(auto iter = map.set_begin(); iter != map.set_end(); iter++)
    {
      for(int comp = 0; comp < iter.numComp(); ++comp)
      {
        sum += iter(comp);  //<----
      }
    }
  }
  timer.stop();
  SLIC_INFO("  Field1D: " << timer.elapsed() << " sec");
  SLIC_ASSERT(c_sum == sum);

  sum = 0;
  timer.reset();
  timer.start();
  {
    auto map2d = mm.get2dField<double>("CellMat Array");
    for(int i = 0; i < map2d.firstSetSize(); i++)
    {
      auto submap = map2d(i);
      for(auto iter = submap.set_begin(); iter != submap.set_end(); iter++)
      {
        //int idx = iter.index();  //get the index
        for(int comp = 0; comp < map2d.numComp(); ++comp)
        {
          sum += iter(comp);                            //<----
          SLIC_ASSERT(iter(comp) == iter.value(comp));  //value()
        }
      }
    }
  }
  timer.stop();
  SLIC_INFO("  Field2D: " << timer.elapsed() << " sec");

  SLIC_ASSERT(x_sum == sum);

  // ---------- iterator for BivariateMap ------------
  SLIC_INFO("\n -- With Map iterators - begin(i) and end(i) -- ");
  sum = 0;
  timer.reset();
  timer.start();
  {
    //_multimat_using_fields_iterator_start
    auto map2d = mm.get2dField<double>("CellMat Array");
    for(int i = 0; i < mm.getNumberOfCells() /*map2d.firstSetSize()*/; i++)
    {
      for(auto iter = map2d.set_begin(i); iter != map2d.set_end(i); iter++)
      {
        // int idx = iter.index(); get the index
        for(int comp = 0; comp < map2d.numComp(); ++comp)
        {
          sum += iter(comp);                            //<----
          SLIC_ASSERT(iter(comp) == iter.value(comp));  //value()
        }
        SLIC_ASSERT(iter(0) == (*iter)[0]);  //2 ways to get the first component
        SLIC_ASSERT(iter(0) == iter.value(0));
      }
    }
    //_multimat_using_fields_iterator_end
  }
  timer.stop();
  SLIC_INFO("  Field2D: " << timer.elapsed() << " sec");
  SLIC_ASSERT(x_sum == sum);

  // -------- range-based for-loop, only works if there is 1 component ---------
  if(ncomp == 1)
  {
    SLIC_INFO("\n -- With range-based for-loop of Map -- ");
    sum = 0;
    timer.reset();
    timer.start();

    auto map2d = mm.get2dField<double>("CellMat Array");
    for(int i = 0; i < mm.getNumberOfCells(); i++)
    {
      for(double val : map2d(i))
      {
        sum += val;  //<----
      }
    }

    timer.stop();
    SLIC_INFO("  Field2D: " << timer.elapsed() << " sec");
    SLIC_ASSERT(x_sum == sum);
  }

  // ---------- flat iterator ------------
  SLIC_INFO("\n -- With BivariateMap flat iterators -- ");
  sum = 0;
  timer.reset();
  timer.start();
  {
    // _multimat_using_fields_flatiter_start
    auto map2d = mm.get2dField<double>("CellMat Array");
    for(auto iter = map2d.set_begin(); iter != map2d.set_end(); ++iter)
    {
      //get the indices
      //int cell_id = iter.firstIndex();
      //int mat_id = iter.secondIndex();
      for(int comp = 0; comp < map2d.numComp(); ++comp)
      {
        double val = iter(comp);               //<----
        SLIC_ASSERT(val == iter.value(comp));  //another way to get the value
        sum += val;
      }
    }
    // _multimat_using_fields_flatiter_end
  }
  timer.stop();
  SLIC_INFO("  Field2D: " << timer.elapsed() << " sec");
  SLIC_ASSERT(x_sum == sum);

  SLIC_INFO("\n");
}

void usage()
{
  SLIC_WARNING("Usage: ./multimat_traversal_ex " << "<num_cells> <num_mats> <num_comp> "
                                                 << "<bool use_sparse> <fill percentage>");
  SLIC_WARNING("  example: ./multimat_traversal_ex 10000 50 3 1 0.2");
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;

  //default options
  int ncells = 10000;
  int nmats = 50;
  int ncomp = 1;
  bool use_sparse = true;
  double fill_percentage = .2;

  if(argc != 1 && argc != 6)
  {
    usage();
    return 1;
  }

  if(argc == 6)
  {
    ncells = std::stoi(argv[1]);
    nmats = std::stoi(argv[2]);
    ncomp = std::stoi(argv[3]);
    int sparse = std::stoi(argv[4]);
    if(sparse == 0)
    {
      use_sparse = false;
    }
    else if(sparse == 1)
    {
      use_sparse = true;
    }
    else
    {
      usage();
      return 1;
    }
    fill_percentage = std::stod(argv[5]);
  }

  SLIC_INFO("Using arguments:");
  SLIC_INFO(" number of cells:      " << ncells);
  SLIC_INFO(" number of materials:  " << nmats);
  SLIC_INFO(" number of components: " << ncomp);
  SLIC_INFO(" use sparse:           " << (use_sparse ? "yes" : "no"));
  SLIC_INFO(" fill percentage:      " << fill_percentage);

  various_traversal_methods(nmats, ncells, ncomp, use_sparse, fill_percentage);

  return 0;
}
