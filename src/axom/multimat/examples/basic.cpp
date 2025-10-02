// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/CLI11.hpp"
#include "axom/core.hpp"
#include "axom/fmt.hpp"
#include "axom/multimat.hpp"
#include "axom/slic.hpp"

#ifdef AXOM_USE_CONDUIT
  #include "conduit.hpp"
  #include "conduit_relay.hpp"
#endif

#include <iostream>
#include <string>
#include <vector>

void addfields(axom::multimat::MultiMat &mm)
{
  // clang-format off

  //_multimat_using_fields_addfields_begin
  constexpr int ncells = 9;
  constexpr int nmats = 3;
  constexpr int nComponents = 1;
  // Add PER_CELL field.
  double perCellData[] = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
  axom::ArrayView<double> perCellAV(perCellData, ncells);
  mm.addField("perCell",
              axom::multimat::FieldMapping::PER_CELL,
              axom::multimat::DataLayout::CELL_DOM,
              axom::multimat::SparsityLayout::DENSE,
              perCellAV,
              nComponents);

  // Add PER_MAT field.
  double perMatData[] = {1., 2., 3.};
  axom::ArrayView<double> perMatAV(perMatData, nmats);
  mm.addField("perMat",
              axom::multimat::FieldMapping::PER_MAT,
              axom::multimat::DataLayout::MAT_DOM,
              axom::multimat::SparsityLayout::DENSE,
              perMatAV,
              nComponents);

  // Add PER_CELL_MAT DENSE field. 0's where there is no material.
  double perCellMatDense[ncells][nmats] = {
    {0.55,  0.45,  0.},   // cell 0
    {0.425, 0.425, 0.15}, // cell 1
    {0.3,   0.4,   0.3},  // cell 2
    {0.425, 0.425, 0.15}, // cell 3
    {0.,    0.2,   0.8},  // cell 4
    {0.,    0.,    1.},   // cell 5
    {0.3,   0.4,   0.3},  // cell 6
    {0.,    0.,    1.},   // cell 7
    {0.,    0.,    1.}    // cell 8
  };
  axom::ArrayView<double> perCellMatDenseAV(&perCellMatDense[0][0], ncells * nmats);
  mm.addField("perCellMatDense",
              axom::multimat::FieldMapping::PER_CELL_MAT,
              axom::multimat::DataLayout::CELL_DOM,
              axom::multimat::SparsityLayout::DENSE,
              perCellMatDenseAV,
              nComponents);

  // Add PER_CELL_MAT SPARSE field. 0's do not need to appear.
  double perCellMatSparse[] = {
    0.55,  0.45,        // cell 0
    0.425, 0.425, 0.15, // cell 1
    0.3,   0.4,   0.3,  // cell 2
    0.425, 0.425, 0.15, // cell 3
           0.2,   0.8,  // cell 4
                  1.,   // cell 5
    0.3,   0.4,   0.3,  // cell 6
                  1.,   // cell 7
                  1.    // cell 8
  };
  axom::ArrayView<double> perCellMatSparseAV(perCellMatSparse,
                                             sizeof(perCellMatSparse) / sizeof(double));
  mm.addField("perCellMatSparse",
              axom::multimat::FieldMapping::PER_CELL_MAT,
              axom::multimat::DataLayout::CELL_DOM,
              axom::multimat::SparsityLayout::SPARSE,
              perCellMatSparseAV,
              nComponents);
  //_multimat_using_fields_addfields_end

  // clang-format on
}

void multicomponent(axom::multimat::MultiMat &mm)
{
  // clang-format off

  //_multimat_using_fields_multicomponent_begin
  constexpr int nComponents = 2;
  double data[] = {0., 0.,  // cell 0 x,y components
                   1., 1.,  // cell 1 x,y components
                   2., 4.,  // cell 2 x,y components
                   3., 9.,  // cell 3 x,y components
                   4., 16., // cell 4 x,y components
                   5., 25., // cell 5 x,y components
                   6., 36., // cell 6 x,y components
                   7., 49., // cell 7 x,y components
                   8., 64.  // cell 8 x,y components
                  };
  axom::ArrayView<double> dataAV(data, sizeof(data) / sizeof(double));
  mm.addField("perCellMC",
              axom::multimat::FieldMapping::PER_CELL,
              axom::multimat::DataLayout::CELL_DOM,
              axom::multimat::SparsityLayout::DENSE,
              dataAV,
              nComponents);
  //_multimat_using_fields_multicomponent_end

  // clang-format on
}

void introspection(axom::multimat::MultiMat &mm)
{
  //_multimat_using_fields_introspection_begin

  // Print the field names in the MultiMat object mm.
  for(int i = 0; i < mm.getNumberOfFields(); i++)
  {
    // Get field properties
    auto name = mm.getFieldName(i);
    auto mapping = mm.getFieldMapping(i);
    auto layout = mm.getFieldDataLayout(i);
    auto sparsity = mm.getFieldSparsityLayout(i);
    auto dataType = mm.getFieldDataType(i);

    std::cout << name << ":"
              << "\n\tmapping: " << mapping << "\n\tlayout: " << layout
              << "\n\tsparsity: " << sparsity << "\n\tdataType: " << dataType << "\n";
  }
  std::cout << "Volfrac index: " << mm.getFieldIdx("Volfrac") << std::endl;

  //_multimat_using_fields_introspection_end
}

void using_fields_index_sets(axom::multimat::MultiMat &mm)
{
  //_multimat_using_fields_index_sets_begin
  // CELL_DOM data (iterate over cells then materials)
  const std::string fieldName("perCellMatSparse");
  auto f = mm.getSparse2dField<double>(fieldName);
  std::cout << "Field: " << fieldName << std::endl;
  for(int i = 0; i < mm.getNumberOfCells(); i++)
  {
    std::cout << "\tcell " << i << " values: ";
    for(const auto &idx : mm.getIndexingSetOfCell(i, axom::multimat::SparsityLayout::SPARSE))
    {
      std::cout << f[idx] << ", ";
    }
    std::cout << "\n";
  }
  //_multimat_using_fields_index_sets_end
}

void using_fields_1d(axom::multimat::MultiMat &mm)
{
  // _multimat_using_fields_1d_start
  // Sum all values in the field.
  double sum = 0.;
  auto f = mm.get1dField<double>("perCell");
  for(int i = 0; i < mm.getNumberOfCells(); i++)
  {
    sum += f[i];
  }
  // _multimat_using_fields_1d_end

  SLIC_INFO(axom::fmt::format("sum={}", sum));
}

void using_fields_multi_component(axom::multimat::MultiMat &mm)
{
  // _multimat_using_fields_1dmc_start
  double sum = 0.;
  auto f = mm.get1dField<double>("perCellMC");
  for(int i = 0; i < mm.getNumberOfCells(); i++)
  {
    for(int comp = 0; comp < f.numComp(); comp++)
    {
      sum += f(i, comp);
    }
  }
  // _multimat_using_fields_1dmc_end

  SLIC_INFO(axom::fmt::format("sum={}", sum));
}

void dynamic_mode(axom::multimat::MultiMat &mm)
{
  //_multimat_dynamic_mode_begin
  // mm is a MultiMat object.

  // Switch to dynamic mode
  mm.convertToDynamic();

  // Add material 1 in zone 7 that was not there before.
  mm.addEntry(7, 1);

  // Remove material 1 in zone 0
  mm.removeEntry(0, 1);

  // Volume fraction updates omitted (iterate Volfrac field, set new values)

  //_multimat_dynamic_mode_end
}

#ifdef AXOM_USE_CONDUIT
/**
 * \brief Turn a MultiMat into a Blueprint matset and add its fields
 *        to the Blueprint mesh too.
 *
 * \param mm The MultiMat object that contains the materials and fields.
 * \param mesh The node that contains the Blueprint mesh.
 */
void multimat_to_blueprint(axom::multimat::MultiMat &mm, conduit::Node &mesh)
{
  // Multimat to matset.
  const auto VF = mm.get2dField<double>("Volfrac");
  std::vector<int> material_ids, indices, sizes, offsets;
  std::vector<double> volume_fractions;
  int offset = 0, idx = 0;
  for(int c = 0; c < mm.getNumberOfCells(); c++)
  {
    auto vf_for_cell = VF(c);
    int size = 0;
    for(int i = 0; i < vf_for_cell.size(); i++)
    {
      double value = vf_for_cell(i);
      if(value > 0.)
      {
        material_ids.push_back(vf_for_cell.index(i));
        volume_fractions.push_back(vf_for_cell(i));
        indices.push_back(idx++);
        size++;
      }
    }
    sizes.push_back(size);
    offsets.push_back(offset);
    offset += size;
  }

  mesh["matsets/matset/topology"] = "main";
  for(int m = 0; m < mm.getNumberOfMaterials(); m++)
  {
    mesh[axom::fmt::format("matsets/matset/material_map/mat{}", m)] = m;
  }
  mesh["matsets/matset/material_ids"].set(material_ids);
  mesh["matsets/matset/volume_fractions"].set(volume_fractions);
  mesh["matsets/matset/indices"].set(indices);
  mesh["matsets/matset/sizes"].set(sizes);
  mesh["matsets/matset/offsets"].set(offsets);

  // MultiMat to fields.
  const std::string compNames[] = {"x", "y", "z"};
  for(int i = 0; i < mm.getNumberOfFields(); i++)
  {
    // Get field properties
    auto name = mm.getFieldName(i);
    auto mapping = mm.getFieldMapping(i);
    SLIC_ASSERT(mm.getFieldDataType(i) == axom::multimat::DataTypeSupported::TypeDouble);

    conduit::Node &n_f = mesh["fields/" + name];
    n_f["association"] = "element";
    n_f["topology"] = "main";

    if(mapping == axom::multimat::FieldMapping::PER_CELL)
    {
      auto f = mm.get1dField<double>(name);
      double *dptr = &f[0];

      if(f.numComp() == 1)
      {
        n_f["values"].set_external(dptr, mm.getNumberOfCells());
      }
      else
      {
        SLIC_ASSERT(f.numComp() <= 3);
        for(int c = 0; c < f.numComp(); c++)
        {
          n_f["values/" + compNames[c]].set_external(dptr,
                                                     mm.getNumberOfCells(),
                                                     0,
                                                     f.numComp() * sizeof(double));
          dptr++;
        }
      }
    }
    else if(mapping == axom::multimat::FieldMapping::PER_MAT)
    {
      auto f = mm.get1dField<double>(name);

      std::vector<double> values;
      for(int c = 0; c < mm.getNumberOfCells(); c++)
      {
        double value = 0.;
        const auto matsInCell = mm.getMatInCell(c);
        // Take the first material - could take the one with largest VF.
        if(matsInCell.size() > 0)
        {
          int mat = matsInCell[0];
          value = f[mat];
        }
        values.push_back(value);
      }

      n_f["values"].set(values);
    }
    else if(mapping == axom::multimat::FieldMapping::PER_CELL_MAT)
    {
      auto f = mm.get2dField<double>(name);

      if(f.numComp() == 1)
      {
        std::vector<double> values, matset_values;
        for(int c = 0; c < mm.getNumberOfCells(); c++)
        {
          const auto matsInCell = mm.getMatInCell(c);
          double avg = 0.;
          for(auto &m : matsInCell)
          {
            double *valptr = f.findValue(c, m);
            if(valptr != nullptr)
            {
              matset_values.push_back(*valptr);
              avg += *valptr;
            }
          }
          values.push_back(avg / std::max(double(matsInCell.size()), 1.));
        }

        n_f["matset"] = "matset";
        n_f["values"].set(values);
        n_f["matset_values"].set(matset_values);
      }
      else
      {
        std::vector<double> values[3], matset_values[3];
        for(int c = 0; c < mm.getNumberOfCells(); c++)
        {
          const auto matsInCell = mm.getMatInCell(c);
          double avg[3] = {0., 0., 0.};
          for(auto &m : matsInCell)
          {
            for(int comp = 0; comp < f.numComp(); comp++)
            {
              double *valptr = f.findValue(c, m, comp);
              if(valptr != nullptr)
              {
                matset_values[comp].push_back(*valptr);
                avg[comp] += *valptr;
              }
            }
          }
          for(int comp = 0; comp < f.numComp(); comp++)
          {
            values[comp].push_back(avg[comp] / std::max(double(matsInCell.size()), 1.));
          }
        }

        n_f["matset"] = "matset";
        for(int comp = 0; comp < f.numComp(); comp++)
        {
          n_f["values/" + compNames[comp]].set(values[comp]);
          n_f["matset_values/" + compNames[comp]].set(matset_values[comp]);
        }
      }
    }
  }
}

/**
 * \brief Save the MultiMat to a Blueprint output file.
 *
 * \param mm The MultiMat object to save. It is compatible with the mesh in this routine.
 */
void save_blueprint(axom::multimat::MultiMat &mm)
{
  const char *yaml = R"(
coordsets:
  coords:
    type: explicit
    values:
     x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
     y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 3., 3., 3., 3.]
topologies:
  main:
    type: unstructured
    coordset: coords
    elements:
      shape: quad
      connectivity: [0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10, 8,9,13,12, 9,10,14,13, 10,11,15,14]
      sizes: [4, 4, 4, 4, 4, 4, 4, 4, 4]
      offsets: [0, 4, 8, 12, 16, 20, 24, 28, 32]
  )";

  conduit::Node mesh;
  mesh.parse(yaml);

  multimat_to_blueprint(mm, mesh);

  SLIC_INFO("Saving MultiMat to Blueprint.");
  conduit::relay::io::blueprint::save_mesh(mesh, "basic", "hdf5");
  conduit::relay::io::save(mesh, "basic.yaml", "yaml");
}
#endif

int main(int argc, char *argv[])
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);
  axom::CLI::App app;
#ifdef AXOM_USE_CONDUIT
  bool output = false;
  app.add_flag("--output", output)
    ->description("Whether to write a Blueprint mesh of the MultiMat data.");
#endif
  // Parse command line options.
  app.parse(argc, argv);

  // NOTE: Construct the MultiMat object here so we can see the
  //       declaration in the docs.

  // clang-format off

  //_multimat_materials_cmr_begin
  constexpr int nmats = 3;
  constexpr int ncells = 9;

  // Create the MultiMat object mm
  axom::multimat::MultiMat mm;
  mm.setNumberOfMaterials(nmats);
  mm.setNumberOfCells(ncells);

  // Cell-Dominant data layout
  int rel[ncells][nmats] = {
    {1,1,0},
    {1,1,1},
    {1,1,1},
    {1,1,1},
    {0,1,1},
    {0,0,1},
    {1,1,1},
    {0,0,1},
    {0,0,1}
  };
  std::vector<bool> relation(ncells * nmats, false);
  for(int c = 0; c < ncells; c++)
  {
    for(int m = 0; m < nmats; m++)
    {
      relation[c * nmats + m] = rel[c][m] > 0;
    }
  }
  mm.setCellMatRel(relation, axom::multimat::DataLayout::CELL_DOM);
  //_multimat_materials_cmr_end

  //_multimat_materials_volfracs_begin
  // Cell-Dominant, DENSE data layout
  double volfracs[ncells][nmats] = {
    {0.55,  0.45,  0.},
    {0.425, 0.425, 0.15},
    {0.3,   0.4,   0.3},
    {0.425, 0.425, 0.15},
    {0.,    0.2,   0.8},
    {0.,    0.,    1.},
    {0.3,   0.4,   0.3},
    {0.,    0.,    1.},
    {0.,    0.,    1.}
  };
  axom::ArrayView<double> vfView(&volfracs[0][0], ncells * nmats);
  mm.setVolfracField(vfView,
                     axom::multimat::DataLayout::CELL_DOM,
                     axom::multimat::SparsityLayout::DENSE);
  //_multimat_materials_volfracs_end

  // clang-format on

  // Demonstrate other features.
  addfields(mm);
  multicomponent(mm);
  introspection(mm);
  using_fields_index_sets(mm);
  using_fields_1d(mm);
  using_fields_multi_component(mm);

#ifdef AXOM_USE_CONDUIT
  if(output)
  {
    save_blueprint(mm);
  }
#endif

  dynamic_mode(mm);
  introspection(mm);
  return 0;
}
