// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom utils
#include "axom/core/utilities/Utilities.hpp" /* for utilities::max */

// Mint includes
#include "axom/mint/config.hpp"                /* for IndexType, int64 */
#include "axom/mint/mesh/CellTypes.hpp"        /* for cell::vtk_types */
#include "axom/mint/mesh/CurvilinearMesh.hpp"  /* for CurvilinearMesh */
#include "axom/mint/mesh/Field.hpp"            /* for Field */
#include "axom/mint/mesh/FieldData.hpp"        /* for FieldData */
#include "axom/mint/mesh/FieldTypes.hpp"       /* for *_FIELD_TYPE */
#include "axom/mint/mesh/FieldVariable.hpp"    /* for FieldVariable */
#include "axom/mint/mesh/Mesh.hpp"             /* for Mesh */
#include "axom/mint/mesh/ParticleMesh.hpp"     /* for ParticleMesh */
#include "axom/mint/mesh/RectilinearMesh.hpp"  /* for RectilinearMesh */
#include "axom/mint/mesh/UniformMesh.hpp"      /* for UniformMesh */
#include "axom/mint/mesh/UnstructuredMesh.hpp" /* for UnstructuredMesh */
#include "axom/mint/utils/vtk_utils.hpp"       /* for write_vtk */
#include "mint_test_utilities.hpp"             /* for create_mesh */

// Slic includes
#include "axom/slic.hpp"

// C/C++ includes
#include <cmath>   /* for std::exp */
#include <fstream> /* for std::ifstream */
#include <iomanip> /* for std::setfill, std::setw */
#include <string>  /* for std::string */
#include <sstream> /* for std::stringstream */
#include <set>     /* for std::set */

// gtest includes
#include "gtest/gtest.h"

#ifndef DELETE_VTK_FILES
  #define DELETE_VTK_FILES 1
#endif

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
//  INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace internal
{
/*!
 * \brief Creates artificial 1 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != nullptr
 */
void create_scalar_data(Mesh* mesh)
{
  const IndexType mesh_num_nodes = mesh->getNumberOfNodes();
  const IndexType mesh_num_cells = mesh->getNumberOfCells();

  double* double_ptr = nullptr;
  int* int_ptr = nullptr;

  double_ptr = mesh->createField<double>("node_scalars_double", mint::NODE_CENTERED);
  int_ptr = mesh->createField<int>("node_scalars_int", mint::NODE_CENTERED);

  for(int idx = 0; idx < mesh_num_nodes; ++idx)
  {
    double coords[] = {0, 0, 0};
    mesh->getNode(idx, coords);

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    double r2 = (x * x) + (y * y) + (z * z);
    r2 = std::exp(-r2 / 100.0);
    double temp = (x * x) - 4 * y + (2 * y * z) - (5 * x * y * z + 1) * r2;
    double_ptr[idx] = temp;
    int_ptr[idx] = static_cast<int>(temp);
  }

  /* Particle meshes can only have node centered fields. */
  if(mesh->getMeshType() == PARTICLE_MESH)
  {
    return;
  }

  double_ptr = mesh->createField<double>("cell_scalars_double", mint::CELL_CENTERED);

  int_ptr = mesh->createField<int>("cell_scalars_int", mint::CELL_CENTERED);

  for(int idx = 0; idx < mesh_num_cells; ++idx)
  {
    double_ptr[idx] = idx;
    int_ptr[idx] = idx;
  }
}

/*!
 * \brief Creates artificial 3 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != nullptr
 */
void create_vector_data(Mesh* mesh)
{
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const IndexType num_cells = mesh->getNumberOfCells();

  double* double_ptr3 = nullptr;
  double* double_ptr2 = nullptr;
  int* int_ptr3 = nullptr;
  int* int_ptr2 = nullptr;

  double_ptr3 = mesh->createField<double>("node_vectors_3double", mint::NODE_CENTERED, 3);
  int_ptr3 = mesh->createField<int>("node_vectors_3int", mint::NODE_CENTERED, 3);
  double_ptr2 = mesh->createField<double>("node_vectors_2double", mint::NODE_CENTERED, 2);
  int_ptr2 = mesh->createField<int>("node_vectors_2int", mint::NODE_CENTERED, 2);

  for(int idx = 0; idx < num_nodes; ++idx)
  {
    double coords[] = {0, 0, 0};
    mesh->getNode(idx, coords);

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    double r2 = (x * x) + (y * y) + (z * z);
    r2 = std::exp(-r2 / 100.0);
    double v1 = (x * x) - 4 * y + (2 * y * z) - (5 * x * y * z + 1) * r2;
    double v2 = (y * y) - 4 * z + (2 * x * z) - (5 * x * y * z + 1) * r2;
    double v3 = (z * z) - 4 * x + (2 * x * y) - (5 * x * y * z + 1) * r2;

    double_ptr3[3 * idx] = v1;
    double_ptr3[3 * idx + 1] = v2;
    double_ptr3[3 * idx + 2] = v3;

    int_ptr3[3 * idx] = static_cast<int>(v1);
    int_ptr3[3 * idx + 1] = static_cast<int>(v2);
    int_ptr3[3 * idx + 2] = static_cast<int>(v3);

    double_ptr2[2 * idx] = v1;
    double_ptr2[2 * idx + 1] = v2;

    int_ptr2[2 * idx] = static_cast<int>(v1);
    int_ptr2[2 * idx + 1] = static_cast<int>(v2);
  }

  /* Particle meshes can only have node centered fields. */
  if(mesh->getMeshType() == PARTICLE_MESH)
  {
    return;
  }

  double_ptr3 = mesh->createField<double>("cell_vectors_3double", mint::CELL_CENTERED, 3);
  int_ptr3 = mesh->createField<int>("cell_vectors_3int", mint::CELL_CENTERED, 3);
  double_ptr2 = mesh->createField<double>("cell_vectors_2double", mint::CELL_CENTERED, 2);
  int_ptr2 = mesh->createField<int>("cell_vectors_2int", mint::CELL_CENTERED, 2);

  for(int idx = 0; idx < num_cells; ++idx)
  {
    double_ptr3[3 * idx] = idx;
    double_ptr3[3 * idx + 1] = idx + 1;
    double_ptr3[3 * idx + 2] = idx + 2;

    int_ptr3[3 * idx] = idx;
    int_ptr3[3 * idx + 1] = idx + 1;
    int_ptr3[3 * idx + 2] = idx + 2;

    double_ptr2[2 * idx] = idx;
    double_ptr2[2 * idx + 1] = idx + 1;

    int_ptr2[2 * idx] = idx;
    int_ptr2[2 * idx + 1] = idx + 1;
  }
}

/*!
 * \brief Creates artificial 4 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != nullptr
 */
void create_multidim_data(Mesh* mesh)
{
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const IndexType num_cells = mesh->getNumberOfCells();

  double* double_ptr = nullptr;
  int* int_ptr = nullptr;

  double_ptr = mesh->createField<double>("node_multidim_double", mint::NODE_CENTERED, 4);
  int_ptr = mesh->createField<int>("node_multidim_int", mint::NODE_CENTERED, 4);

  for(int idx = 0; idx < num_nodes; ++idx)
  {
    double coords[] = {0, 0, 0};
    mesh->getNode(idx, coords);

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    double r2 = (x * x) + (y * y) + (z * z);
    r2 = std::exp(-r2 / 100.0);
    double v1 = (x * x) - 4 * y + (2 * y * z) - (5 * x * y * z + 1) * r2;
    double v2 = (y * y) - 4 * z + (2 * x * z) - (5 * x * y * z + 1) * r2;
    double v3 = (z * z) - 4 * x + (2 * x * y) - (5 * x * y * z + 1) * r2;
    double v4 = (x * x + y * y + z * z) * r2;

    double_ptr[4 * idx + 0] = v1;
    double_ptr[4 * idx + 1] = v2;
    double_ptr[4 * idx + 2] = v3;
    double_ptr[4 * idx + 3] = v4;

    int_ptr[4 * idx + 0] = static_cast<int>(v1);
    int_ptr[4 * idx + 1] = static_cast<int>(v2);
    int_ptr[4 * idx + 2] = static_cast<int>(v3);
    int_ptr[4 * idx + 3] = static_cast<int>(v4);
  }

  /* Particle meshes can only have node centered fields. */
  if(mesh->getMeshType() == PARTICLE_MESH)
  {
    return;
  }

  double_ptr = mesh->createField<double>("cell_multidim_double", mint::CELL_CENTERED, 4);
  int_ptr = mesh->createField<int>("cell_multidim_int", mint::CELL_CENTERED, 4);

  for(int idx = 0; idx < num_cells; ++idx)
  {
    double_ptr[4 * idx + 0] = idx;
    double_ptr[4 * idx + 1] = idx + 1;
    double_ptr[4 * idx + 2] = idx + 2;
    double_ptr[4 * idx + 3] = idx + 3;

    int_ptr[4 * idx + 0] = idx;
    int_ptr[4 * idx + 1] = idx + 1;
    int_ptr[4 * idx + 2] = idx + 2;
    int_ptr[4 * idx + 3] = idx + 3;
  }
}

/*!
 * \brief Creates artificial mesh data and then writes the mesh out to disk.
 * \param [in] mesh the mesh to write out.
 * \param [in] path the path of the file to be written.
 * \pre mesh != nullptr
 */
void populate_and_write(Mesh* mesh, const std::string& path)
{
  create_scalar_data(mesh);
  create_vector_data(mesh);
  create_multidim_data(mesh);
  write_vtk(mesh, path);
}

/*!
 * \brief Checks that the VTK header was written correctly.
 * \param [in] file the file to parse.
 */
void check_header(std::ifstream& file)
{
  std::string buffer;
  std::getline(file, buffer);
  EXPECT_EQ(buffer, "# vtk DataFile Version 3.0");
  std::getline(file, buffer);
  EXPECT_EQ(buffer, "Mesh generated by axom::mint::write_vtk");
  std::getline(file, buffer);
  EXPECT_EQ(buffer, "ASCII");
}

/*!
 * \brief Checks that a scalar field was written correctly.
 * \param [in] field the field to check against.
 * \param [in] file the file to parse.
 * \param [in] offset the offset into the field to start at.
 * \pre field != nullptr
 */
void check_scalar(const Field* field, std::ifstream& file, unsigned int offset = 0)
{
  const int num_components = field->getNumComponents();
  const IndexType num_values = field->getNumTuples();

  std::string buffer;
  file >> buffer;
  EXPECT_EQ(buffer, "LOOKUP_TABLE");
  file >> buffer;
  EXPECT_EQ(buffer, "default");

  if(field->getType() == DOUBLE_FIELD_TYPE)
  {
    const double* field_data = Field::getDataPtr<double>(field);
    for(IndexType idx = 0; idx < num_values; ++idx)
    {
      double temp;
      file >> temp;
      EXPECT_DOUBLE_EQ(temp, field_data[num_components * idx + offset]);
    }
  }
  else if(field->getType() == INT32_FIELD_TYPE)
  {
    const int* field_data = Field::getDataPtr<int>(field);
    for(IndexType idx = 0; idx < num_values; ++idx)
    {
      int temp;
      file >> temp;
      EXPECT_EQ(field_data[num_components * idx + offset], temp);
    }
  }
}

/*!
 * \brief Checks that a vector field was written correctly.
 * \param [in] field the vector field to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 * \pre field->getNumComponents() == 2 || field->getNumComponents() == 3
 */
void check_vector_data(const Field* field, std::ifstream& file)
{
  const int num_components = field->getNumComponents();
  const IndexType num_values = field->getNumTuples();

  if(field->getType() == DOUBLE_FIELD_TYPE)
  {
    const double* field_data = Field::getDataPtr<double>(field);
    double temp;
    for(IndexType idx = 0; idx < num_values; ++idx)
    {
      for(int dim = 0; dim < num_components; ++dim)
      {
        file >> temp;
        EXPECT_DOUBLE_EQ(temp, field_data[num_components * idx + dim]);
      }

      if(num_components == 2)
      {
        file >> temp;
        EXPECT_EQ(temp, 0.0);
      }
    }
  }
  else if(field->getType() == INT32_FIELD_TYPE)
  {
    const int* field_data = Field::getDataPtr<int>(field);
    int temp;
    for(IndexType idx = 0; idx < num_values; ++idx)
    {
      for(int dim = 0; dim < num_components; ++dim)
      {
        file >> temp;
        EXPECT_EQ(temp, field_data[num_components * idx + dim]);
      }

      if(num_components == 2)
      {
        file >> temp;
        EXPECT_EQ(temp, 0);
      }
    }
  }
}

/*!
 * \brief Checks that a multidimensional field was written correctly.
 * \param [in] field the multidimensional field to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 * \pre field->getNumComponents() > 3
 */
void check_multidim_data(const Field* field, std::ifstream& file)
{
  const int num_components = field->getNumComponents();
  const int field_type = field->getType();
  check_scalar(field, file, 0);

  std::string type, name, d_type;
  for(int comp = 1; comp < num_components; ++comp)
  {
    file >> type >> name >> d_type;
    EXPECT_EQ(type, "SCALARS");

    std::stringstream temp;
    temp << field->getName() << "_";
    temp << std::setfill('0') << std::setw(3) << comp;
    EXPECT_EQ(name, temp.str());

    if(d_type == "double")
    {
      EXPECT_EQ(field_type, DOUBLE_FIELD_TYPE);
    }
    else if(d_type == "int")
    {
      EXPECT_EQ(field_type, INT32_FIELD_TYPE);
    }
    else
    {
      EXPECT_TRUE(false) << "Unknown field data type: " << d_type;
    }

    check_scalar(field, file, comp);
  }
}

/*!
 * \brief Checks that field data was written correctly.
 * \param [in] field_data the field data to check against.
 * \param [in] file the file to parse.
 * \pre field_data != nullptr
 */
void check_fieldData(const FieldData* field_data, std::ifstream& file)
{
  std::set<std::string> fields_read;
  std::string type, name, d_type;
  size_t cur_pos;
  while(file.good())
  {
    file >> type >> name >> d_type;

    if(field_data->hasField(name))
    {
      fields_read.insert(name);
      const Field* field = field_data->getField(name);

      if(d_type == "double")
      {
        EXPECT_EQ(field->getType(), DOUBLE_FIELD_TYPE);
      }
      else if(d_type == "int")
      {
        EXPECT_EQ(field->getType(), INT32_FIELD_TYPE);
      }
      else
      {
        EXPECT_TRUE(false) << "Unknown field data type: " << d_type;
      }

      if(type == "SCALARS")
      {
        EXPECT_TRUE(field->getNumComponents() == 1);
        check_scalar(field, file);
      }
      else if(type == "VECTORS")
      {
        EXPECT_TRUE(field->getNumComponents() == 2 || field->getNumComponents() == 3);
        check_vector_data(field, file);
      }
    }
    else
    {
      size_t underscore_pos = name.size() - 4;
      std::string true_name = name.substr(0, underscore_pos);
      EXPECT_EQ(name.substr(underscore_pos), "_000");
      ASSERT_TRUE(field_data->hasField(true_name)) << true_name;
      fields_read.insert(true_name);

      const Field* field = field_data->getField(true_name);

      if(d_type == "double")
      {
        EXPECT_EQ(field->getType(), DOUBLE_FIELD_TYPE);
      }
      else if(d_type == "int")
      {
        EXPECT_EQ(field->getType(), INT32_FIELD_TYPE);
      }
      else
      {
        EXPECT_TRUE(false) << "Unknown field data type: " << d_type;
      }

      EXPECT_TRUE(field->getNumComponents() > 3);
      check_multidim_data(field, file);
    }

    cur_pos = file.tellg();
    file >> type;
    file.seekg(cur_pos);
    if(type == "CELL_DATA" || type == "POINT_DATA")
    {
      break;
    }
  }

  EXPECT_EQ(static_cast<int>(fields_read.size()), field_data->getNumFields());
  for(std::set<std::string>::iterator it = fields_read.begin(); it != fields_read.end(); ++it)
  {
    EXPECT_TRUE(field_data->hasField(*it));
  }
}

/*!
 * \brief Checks that the mesh cell and node data were written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
void check_data(const Mesh* mesh, std::ifstream& file)
{
  std::string data_type;
  IndexType data_size;
  while(file.good())
  {
    file >> data_type >> data_size;
    if(data_type == "POINT_DATA")
    {
      EXPECT_EQ(data_size, mesh->getNumberOfNodes());
      const FieldData* node_data = mesh->getFieldData(mint::NODE_CENTERED);
      check_fieldData(node_data, file);
    }
    else if(data_type == "CELL_DATA")
    {
      EXPECT_EQ(data_size, mesh->getNumberOfCells());
      const FieldData* cell_data = mesh->getFieldData(mint::CELL_CENTERED);
      check_fieldData(cell_data, file);
    }
  }
}

/*!
 * \brief Checks that the StructuredMesh nodal dimensions were written 
 *  correctly.
 * \param [in] s_mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
void check_dimensions(const StructuredMesh* s_mesh, std::ifstream& file)
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(buffer, "DIMENSIONS");

  for(int i = 0; i < 3; ++i)
  {
    IndexType temp;
    file >> temp;
    if(i < s_mesh->getDimension())
    {
      EXPECT_EQ(temp, s_mesh->getNodeResolution(i));
    }
    else
    {
      EXPECT_EQ(temp, 1);
    }
  }
}

/*!
 * \brief Checks that the mesh node coordinates were written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
void check_points(const Mesh* mesh, std::ifstream& file)
{
  const IndexType num_nodes = mesh->getNumberOfNodes();

  std::string extracted_name, extracted_type;
  IndexType extracted_size;
  file >> extracted_name >> extracted_size >> extracted_type;
  EXPECT_EQ(extracted_name, "POINTS");
  EXPECT_EQ(extracted_size, num_nodes);
  EXPECT_EQ(extracted_type, "double");

  double extracted_coord = 0.0;
  for(IndexType idx = 0; idx < num_nodes; ++idx)
  {
    double coords[] = {0, 0, 0};
    mesh->getNode(idx, coords);

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    file >> extracted_coord;
    EXPECT_EQ(extracted_coord, x);

    file >> extracted_coord;
    EXPECT_EQ(extracted_coord, y);

    file >> extracted_coord;
    EXPECT_EQ(extracted_coord, z);
  }
}

/*!
 * \brief Checks that the mesh cells were written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
void check_cells(const Mesh* mesh, std::ifstream& file)
{
  const IndexType num_cells = mesh->getNumberOfCells();

  /* First need to get total size of the connectivity array. */
  /* If the mesh only has one cell type we can calculate this directly. */
  int max_cell_nodes = mesh->getNumberOfCellNodes(0);
  IndexType total_size = (max_cell_nodes + 1) * num_cells;

  /* If the mesh has mixed cells then we need to loop over the elements. */
  if(mesh->hasMixedCellTypes())
  {
    total_size = num_cells;
    for(IndexType cellIdx = 0; cellIdx < num_cells; ++cellIdx)
    {
      const int num_cell_nodes = mesh->getNumberOfCellNodes(cellIdx);
      max_cell_nodes = utilities::max(num_cell_nodes, max_cell_nodes);
      total_size += num_cell_nodes;
    }
  }

  std::string type;
  IndexType extracted_cells, extracted_size;
  file >> type >> extracted_cells >> extracted_size;
  EXPECT_EQ(type, "CELLS");
  EXPECT_EQ(extracted_cells, num_cells);
  EXPECT_EQ(extracted_size, total_size);

  /* Write out the mesh cell connectivity. */
  IndexType temp;
  IndexType* cell_nodes = new IndexType[max_cell_nodes];
  for(IndexType cellIdx = 0; cellIdx < num_cells; ++cellIdx)
  {
    const int num_cell_nodes = mesh->getNumberOfCellNodes(cellIdx);
    mesh->getCellNodeIDs(cellIdx, cell_nodes);

    file >> temp;
    EXPECT_EQ(temp, num_cell_nodes);
    for(int i = 0; i < num_cell_nodes; ++i)
    {
      file >> temp;
      EXPECT_EQ(temp, cell_nodes[i]);
    }
  }

  /* Write out the mesh cell types. */
  file >> type >> extracted_cells;
  EXPECT_EQ(type, "CELL_TYPES");
  EXPECT_EQ(extracted_cells, num_cells);
  for(IndexType cellIdx = 0; cellIdx < num_cells; ++cellIdx)
  {
    CellType cell_type = mesh->getCellType(cellIdx);
    file >> temp;
    EXPECT_EQ(temp, getCellInfo(cell_type).vtk_type);
  }

  delete[] cell_nodes;
  cell_nodes = nullptr;
}

/*!
 * \brief Checks that the unstructured or particle mesh header was written
 *  correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
template <class MeshType>
void check_mesh(const MeshType* mesh, std::ifstream& file)
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(buffer, "DATASET");
  file >> buffer;
  EXPECT_EQ(buffer, "UNSTRUCTURED_GRID");

  check_points(mesh, file);
  check_cells(mesh, file);
}

/*!
 * \brief Checks that the uniform mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
template <>
void check_mesh(const UniformMesh* mesh, std::ifstream& file)
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(buffer, "DATASET");
  file >> buffer;
  EXPECT_EQ(buffer, "STRUCTURED_POINTS");

  check_dimensions(mesh, file);

  const double* origin = mesh->getOrigin();
  file >> buffer;
  EXPECT_EQ(buffer, "ORIGIN");
  for(int i = 0; i < 3; ++i)
  {
    double temp;
    file >> temp;
    EXPECT_DOUBLE_EQ(temp, origin[i]);
  }

  const double* spacing = mesh->getSpacing();
  file >> buffer;
  EXPECT_EQ(buffer, "SPACING");
  for(int i = 0; i < 3; ++i)
  {
    double temp;
    file >> temp;
    EXPECT_DOUBLE_EQ(temp, spacing[i]);
  }
}

/*!
 * \brief Checks that the rectilinear mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
template <>
void check_mesh(const RectilinearMesh* r_mesh, std::ifstream& file)
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(buffer, "DATASET");
  file >> buffer;
  EXPECT_EQ(buffer, "RECTILINEAR_GRID");

  check_dimensions(r_mesh, file);

  std::string coord_names[3] = {"X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES"};
  std::string extracted_name, extracted_type;
  IndexType extracted_size;
  double extracted_coord;
  for(int dim = 0; dim < r_mesh->getDimension(); ++dim)
  {
    file >> extracted_name >> extracted_size >> extracted_type;
    EXPECT_EQ(extracted_name, coord_names[dim]);
    EXPECT_EQ(extracted_size, r_mesh->getNodeResolution(dim));
    EXPECT_EQ(extracted_type, "double");

    const double* coord_array = r_mesh->getCoordinateArray(dim);
    for(IndexType i = 0; i < r_mesh->getNodeResolution(dim); ++i)
    {
      file >> extracted_coord;
      EXPECT_EQ(extracted_coord, coord_array[i]);
    }
  }
  for(int dim = r_mesh->getDimension(); dim < 3; ++dim)
  {
    file >> extracted_name >> extracted_size >> extracted_type;
    file >> extracted_coord;

    EXPECT_EQ(extracted_name, coord_names[dim]);
    EXPECT_EQ(extracted_size, 1);
    EXPECT_EQ(extracted_type, "double");
    EXPECT_EQ(extracted_coord, 0.0);
  }
}

/*!
 * \brief Checks that the curvilinear mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != nullptr
 */
template <>
void check_mesh(const CurvilinearMesh* mesh, std::ifstream& file)
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(buffer, "DATASET");
  file >> buffer;
  EXPECT_EQ(buffer, "STRUCTURED_GRID");

  check_dimensions(mesh, file);

  check_points(mesh, file);
}

template <int MeshType, int Topology = SINGLE_SHAPE>
Mesh* build_mesh(int dimension)
{
  const IndexType Ni = 7;
  const IndexType Nj = (dimension >= 2) ? Ni + 1 : -1;
  const IndexType Nk = (dimension == 3) ? Ni + 2 : -1;

  const double lo[] = {-10, -9, -8};
  const double hi[] = {10, 9, 8};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  Mesh* test_mesh = create_mesh<MeshType, Topology>(uniform_mesh);
  EXPECT_TRUE(test_mesh != nullptr);

  return test_mesh;
}

template <class MeshType>
void test_mesh(MeshType* mesh, const std::string& path)
{
  populate_and_write(mesh, path);
  std::ifstream file(path.c_str());
  ASSERT_TRUE(file);
  check_header(file);
  check_mesh(mesh, file);
  check_data(mesh, file);

  file.close();
  delete mesh;
#if DELETE_VTK_FILES
  EXPECT_EQ(axom::utilities::filesystem::removeFile(path), 0);
#endif
}

} /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

/*!
 * \brief Creates a UniformMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST(mint_util_write_vtk, UniformMesh)
{
  for(int dim = 1; dim <= 3; ++dim)
  {
    const std::string path = "uniformMesh" + std::to_string(dim) + "D.vtk";
    UniformMesh* mesh = static_cast<UniformMesh*>(internal::build_mesh<STRUCTURED_UNIFORM_MESH>(dim));

    internal::test_mesh(mesh, path);
  }
}

/*!
 * \brief Creates a RectilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST(mint_util_write_vtk, RectilinearMesh)
{
  for(int dim = 1; dim <= 3; ++dim)
  {
    const std::string path = "rectilinearMesh" + std::to_string(dim) + "D.vtk";
    RectilinearMesh* mesh =
      static_cast<RectilinearMesh*>(internal::build_mesh<STRUCTURED_RECTILINEAR_MESH>(dim));

    internal::test_mesh(mesh, path);
  }
}

/*!
 * \brief Creates a CurvilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST(mint_util_write_vtk, CurvilinearMesh)
{
  for(int dim = 1; dim <= 3; ++dim)
  {
    const std::string path = "curvilinearMesh" + std::to_string(dim) + "D.vtk";
    CurvilinearMesh* mesh =
      static_cast<CurvilinearMesh*>(internal::build_mesh<STRUCTURED_CURVILINEAR_MESH>(dim));

    internal::test_mesh(mesh, path);
  }
}

/*!
 * \brief Creates a UnstructuredMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST(mint_util_write_vtk, UnstructuredMesh)
{
  for(int dim = 1; dim <= 3; ++dim)
  {
    const std::string path = "unstructuredMesh" + std::to_string(dim) + "D.vtk";
    UnstructuredMesh<SINGLE_SHAPE>* mesh = static_cast<UnstructuredMesh<SINGLE_SHAPE>*>(
      internal::build_mesh<UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim));

    internal::test_mesh(mesh, path);
  }
}

/*!
 * \brief Creates a UnstructuredMesh with mixed elements and writes it out to
 *  disk using mint::write_vtk and then reads the file back in to check for
 *  correctness. The UnstructuredMesh consists of one hexahedron with pyramids
 *  on the faces.
 */
TEST(mint_util_write_vtk, UnstructuredMixedMesh)
{
  for(int dim = 1; dim <= 3; ++dim)
  {
    const std::string path = "unstructuredMixedMesh" + std::to_string(dim) + "D.vtk";
    UnstructuredMesh<MIXED_SHAPE>* mesh = static_cast<UnstructuredMesh<MIXED_SHAPE>*>(
      internal::build_mesh<UNSTRUCTURED_MESH, MIXED_SHAPE>(dim));

    internal::test_mesh(mesh, path);
  }
}

/*!
 * \brief Creates a ParticleMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST(mint_util_write_vtk, ParticleMesh)
{
  for(int dim = 1; dim <= 3; ++dim)
  {
    const std::string path = "particleMesh" + std::to_string(dim) + "D.vtk";
    ParticleMesh* mesh = static_cast<ParticleMesh*>(internal::build_mesh<PARTICLE_MESH>(dim));

    internal::test_mesh(mesh, path);
  }
}

} /* end namespace mint */
} /* end namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();
  return result;
}
