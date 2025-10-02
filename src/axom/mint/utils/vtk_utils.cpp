// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Macros.hpp"

#include "axom/mint/utils/vtk_utils.hpp"  // file header

#include "axom/core/utilities/Utilities.hpp"  // for utilities::max
#include "axom/core/numerics/Matrix.hpp"      // for numerics::Matrix

#include "axom/mint/mesh/CellTypes.hpp"        // for cell::vtk_types
#include "axom/mint/mesh/Field.hpp"            // for Field
#include "axom/mint/mesh/FieldData.hpp"        // for FieldData
#include "axom/mint/mesh/FieldTypes.hpp"       // for *_FIELD_TYPE
#include "axom/mint/fem/FiniteElement.hpp"     // for mint::FiniteElement
#include "axom/mint/mesh/Mesh.hpp"             // for Mesh
#include "axom/mint/mesh/MeshTypes.hpp"        // for MINT_*_*_MESH
#include "axom/mint/mesh/ParticleMesh.hpp"     // for ParticleMesh
#include "axom/mint/mesh/RectilinearMesh.hpp"  // for RectilinearMesh
#include "axom/mint/mesh/StructuredMesh.hpp"   // for StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"      // for UniformMesh

#include "axom/slic/interface/slic.hpp"  // for slic macros

#include "axom/fmt.hpp"

// C/C++ includes
#include <fstream>  // for std::ofstream
#include <string>   // for std::string

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
 * \brief Computes the maximum number of nodes of a cell on the given mesh.
 *
 * \param [in] mesh the mesh object
 * \return max_cell_nodes the max number of nodes of a cell
 *
 * \pre mesh != nullptr
 * \post max_cell_nodes >= 1
 */
IndexType get_max_cell_nodes(const Mesh* mesh, IndexType& total_cell_nodes)
{
  SLIC_ASSERT(mesh != nullptr);

  int max_cell_nodes = 0;

  if(!mesh->hasMixedCellTypes())
  {
    // short-circuit
    max_cell_nodes = mesh->getNumberOfCellNodes();
    total_cell_nodes = max_cell_nodes * mesh->getNumberOfCells();
    return max_cell_nodes;
  }

  total_cell_nodes = 0;
  const axom::IndexType numCells = mesh->getNumberOfCells();
  for(axom::IndexType icell = 0; icell < numCells; ++icell)
  {
    CellType cell_type = mesh->getCellType(icell);
    const int num_nodes = mint::getCellInfo(cell_type).num_nodes;
    total_cell_nodes += num_nodes;
    if(num_nodes > max_cell_nodes)
    {
      max_cell_nodes = num_nodes;
    }
  }

  return max_cell_nodes;
}

/*!
 * \brief Writes mesh node locations to a VTK file in legacy ASCII format.
 * \param [in] mesh the mesh whose nodes will be written.
 * \param [in] file the stream to write to.
 * \pre mesh != nullptr
 */
void write_points(const Mesh* mesh, std::ofstream& file)
{
  SLIC_ASSERT(mesh != nullptr);
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const int mesh_dim = mesh->getDimension();

  const double* x = mesh->getCoordinateArray(X_COORDINATE);
  SLIC_ASSERT(x != nullptr);

  const double* y = (mesh_dim > 1) ? mesh->getCoordinateArray(Y_COORDINATE) : nullptr;
  const double* z = (mesh_dim > 2) ? mesh->getCoordinateArray(Z_COORDINATE) : nullptr;

  fmt::print(file, "POINTS {} double\n", num_nodes);
  for(IndexType nodeIdx = 0; nodeIdx < num_nodes; ++nodeIdx)
  {
    double xx = x[nodeIdx];
    double yy = (y != nullptr) ? y[nodeIdx] : 0.0;
    double zz = (z != nullptr) ? z[nodeIdx] : 0.0;
    fmt::print(file, "{} {} {}\n", xx, yy, zz);
  }
}

/*!
 * \brief Writes mesh cell connectivity and type to a VTK file
 *  using the legacy ASCII format.
 * \param [in] mesh the mesh whose cells will be written.
 * \param [in] file the stream to write to.
 * \pre mesh != nullptr
 */
void write_cells(const Mesh* mesh, std::ofstream& file)
{
  SLIC_ASSERT(mesh != nullptr);
  const IndexType num_cells = mesh->getNumberOfCells();

  /* First need to get total size of the connectivity array. */
  IndexType total_size;
  int max_cell_nodes = get_max_cell_nodes(mesh, total_size);
  total_size += num_cells;

  fmt::print(file, "CELLS {} {}\n", num_cells, total_size);

  /* Write out the mesh cell connectivity. */
  IndexType* cell_nodes = new IndexType[max_cell_nodes];
  for(IndexType cellIdx = 0; cellIdx < num_cells; ++cellIdx)
  {
    const int num_cell_nodes = mesh->getNumberOfCellNodes(cellIdx);
    mesh->getCellNodeIDs(cellIdx, cell_nodes);

    fmt::print(file,
               "{} {}\n",
               num_cell_nodes,
               fmt::join(cell_nodes, cell_nodes + num_cell_nodes, " "));
  }

  delete[] cell_nodes;

  /* Write out the mesh cell types. */
  fmt::print(file, "CELL_TYPES {}\n", num_cells);
  for(IndexType cellIdx = 0; cellIdx < num_cells; ++cellIdx)
  {
    CellType cell_type = mesh->getCellType(cellIdx);
    fmt::print(file, "{}\n", getCellInfo(cell_type).vtk_type);
  }
}

/*!
 * \brief Writes the dimensions of a structured mesh to a VTK file using the
 *  legacy ASCII format.
 * \param [in] mesh the structured mesh to write out.
 * \param [in] file the stream to write to.
 * \pre mesh != nullptr
 */
void write_dimensions(const StructuredMesh* mesh, std::ofstream& file)
{
  SLIC_ASSERT(mesh != nullptr);

  const int ndims = mesh->getDimension();
  SLIC_ASSERT(1 <= ndims && ndims <= 3);

  fmt::print(file, "DIMENSIONS ");
  if(ndims == 1)
  {
    fmt::print(file, "{} 1 1\n", mesh->getNodeResolution(0));
  }
  else if(ndims == 2)
  {
    fmt::print(file, "{} {} 1\n", mesh->getNodeResolution(0), mesh->getNodeResolution(1));
  }
  else
  {
    fmt::print(file,
               "{} {} {}\n",
               mesh->getNodeResolution(0),
               mesh->getNodeResolution(1),
               mesh->getNodeResolution(2));
  }
}

/*!
 * \brief Writes a rectilinear mesh to a VTK file using the legacy
 *  ASCII format.
 * \param [in] mesh the rectilinear mesh to write out.
 * \param [in] file the stream to write to.
 * \pre mesh != nullptr
 */
void write_rectilinear_mesh(const RectilinearMesh* mesh, std::ofstream& file)
{
  SLIC_ASSERT(mesh != nullptr);

  write_dimensions(mesh, file);

  std::string coord_names[3] = {"X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES"};

  for(int dim = 0; dim < mesh->getDimension(); ++dim)
  {
    fmt::print(file, "{} {} double\n", coord_names[dim], mesh->getNodeResolution(dim));
    const double* coords = mesh->getCoordinateArray(dim);
    fmt::print(file, "{}\n", fmt::join(coords, coords + mesh->getNodeResolution(dim), " "));
  }
  for(int dim = mesh->getDimension(); dim < 3; ++dim)
  {
    fmt::print(file, "{} 1 double\n0.0\n", coord_names[dim]);
  }
}

/*!
 * \brief Writes a uniform mesh to a VTK file using the legacy ASCII format.
 *
 * \param [in] mesh the uniform mesh to write out.
 * \param [in] file the stream to write to.
 * \pre mesh != nullptr
 */
void write_uniform_mesh(const UniformMesh* mesh, std::ofstream& file)
{
  SLIC_ASSERT(mesh != nullptr);

  write_dimensions(mesh, file);

  const double* temp = mesh->getOrigin();
  fmt::print(file, "ORIGIN {} {} {}\n", temp[0], temp[1], temp[2]);

  temp = mesh->getSpacing();
  fmt::print(file, "SPACING {} {} {}\n", temp[0], temp[1], temp[2]);
}

/*!
 * Templated helper function to write out scalar field of type \a T
 * See write_scalar_data for information about parameters
 */
template <typename T>
void write_scalar_helper(const std::string& type, const Field* field, std::ofstream& file)
{
  const T* data_ptr = Field::getDataPtr<T>(field);
  SLIC_ASSERT(data_ptr != nullptr);

  fmt::print(file, "SCALARS {} ", field->getName());
  fmt::print(file, "{}\n", type);
  fmt::print(file, "LOOKUP_TABLE default\n");
  const IndexType num_values = field->getNumTuples();
  fmt::print(file, "{}\n", fmt::join(data_ptr, data_ptr + num_values, "\n"));
}

/*!
 * \brief Writes a scalar field to a VTK file using the legacy ASCII format
 *
 * \param [in] field the scalar field to write out.
 * \param [in] file the stream to write to.
 * \pre field != nullptr
 * \pre field->getNumComponents() == 1
 */
void write_scalar_data(const Field* field, std::ofstream& file)
{
  SLIC_ASSERT(field != nullptr);
  SLIC_ASSERT(field->getNumComponents() == 1);

  switch(field->getType())
  {
  case FLOAT_FIELD_TYPE:
    write_scalar_helper<float>("float", field, file);
    break;
  case DOUBLE_FIELD_TYPE:
    write_scalar_helper<double>("double", field, file);
    break;
  case INT32_FIELD_TYPE:
    write_scalar_helper<std::int32_t>("int", field, file);
    break;
  case INT64_FIELD_TYPE:
    write_scalar_helper<std::int64_t>("long", field, file);
    break;
  default:
    SLIC_WARNING(fmt::format("Unsupported scalar field type ({}) for field '{}'",
                             field->getType(),
                             field->getName()));
    break;
  }
}

/*!
 * Templated helper function to write out vector field of type \a T
 * See write_vector_data for information about parameters
 */
template <typename T>
void write_vector_helper(const std::string& type, const Field* field, std::ofstream& file)
{
  const T* data_ptr = Field::getDataPtr<T>(field);
  SLIC_ASSERT(data_ptr != nullptr);

  fmt::print(file, "VECTORS {} ", field->getName());
  fmt::print(file, "{}\n", type);

  const int num_components = field->getNumComponents();
  const IndexType num_values = field->getNumTuples();
  for(IndexType i = 0; i < num_values; ++i)
  {
    fmt::print(file,
               "{} {} {}\n",
               data_ptr[num_components * i + 0],
               data_ptr[num_components * i + 1],
               num_components == 2 ? 0. : data_ptr[num_components * i + 2]);
  }
}

/*!
 * \brief Writes a vector field to a VTK file using the legacy ASCII format
 *
 * \param [in] field the vector field to write out.
 * \param [in] file the stream to write to.
 * \pre field != nullptr
 * \pre field->getNumComponents() == 2 || field->getNumComponents() == 3
 */
void write_vector_data(const Field* field, std::ofstream& file)
{
  SLIC_ASSERT(field != nullptr);
  const int num_components = field->getNumComponents();
  AXOM_UNUSED_VAR(num_components);  // silence warning in release configs
  SLIC_ASSERT(num_components == 2 || num_components == 3);

  switch(field->getType())
  {
  case FLOAT_FIELD_TYPE:
    write_vector_helper<float>("float", field, file);
    break;
  case DOUBLE_FIELD_TYPE:
    write_vector_helper<double>("double", field, file);
    break;
  case INT32_FIELD_TYPE:
    write_vector_helper<std::int32_t>("int", field, file);
    break;
  case INT64_FIELD_TYPE:
    write_vector_helper<std::int64_t>("long", field, file);
    break;
  default:
    SLIC_WARNING(fmt::format("Unsupported vector field type ({}) for field '{}'",
                             field->getType(),
                             field->getName()));
    break;
  }
}

/*!
 * Templated helper function to write out multidim field of type \a T
 * See write_multidim_data for information about parameters
 */
template <typename T>
void write_multidim_helper(const std::string& type, const Field* field, std::ofstream& file)
{
  const T* data_ptr = Field::getDataPtr<T>(field);
  SLIC_ASSERT(data_ptr != nullptr);

  const int num_components = field->getNumComponents();
  const IndexType num_values = field->getNumTuples();
  SLIC_ASSERT(num_components > 3);

  for(int cur_comp = 0; cur_comp < num_components; ++cur_comp)
  {
    fmt::print(file, "SCALARS {}_{:0>3} {}\n", field->getName(), cur_comp, type);
    fmt::print(file, "LOOKUP_TABLE default\n");

    for(IndexType i = 0; i < num_values; ++i)
    {
      fmt::print(file, "{}\n", data_ptr[num_components * i + cur_comp]);
    }
  }
}

/*!
 * \brief Writes a multidimensional field to a VTK file using the legacy
 *  ASCII format.
 * \param [in] field the multidimensional field to write out.
 * \param [in] file the stream to write to.
 * \pre field != nullptr
 * \pre field->getNumComponents > 3
 */
void write_multidim_data(const Field* field, std::ofstream& file)
{
  SLIC_ASSERT(field != nullptr);

  switch(field->getType())
  {
  case FLOAT_FIELD_TYPE:
    write_multidim_helper<float>("float", field, file);
    break;
  case DOUBLE_FIELD_TYPE:
    write_multidim_helper<double>("double", field, file);
    break;
  case INT32_FIELD_TYPE:
    write_multidim_helper<std::int32_t>("int", field, file);
    break;
  case INT64_FIELD_TYPE:
    write_multidim_helper<std::int64_t>("long", field, file);
    break;
  default:
    SLIC_WARNING(fmt::format("Unsupported multidim field type ({}) for field '{}'",
                             field->getType(),
                             field->getName()));
    break;
  }
}

/*!
 * \brief Writes mesh FieldData to a VTK file using the legacy ASCII format.
 * \param [in] field_data the data to write out.
 * \param [in] num_values the number of tuples each field is expected to have.
 * \param [in] file the stream to write to.
 * \pre field_data != nullptr
 */
void write_data(const FieldData* field_data, IndexType AXOM_DEBUG_PARAM(num_values), std::ofstream& file)
{
  const int numFields = field_data->getNumFields();
  for(int i = 0; i < numFields; ++i)
  {
    const Field* field = field_data->getField(i);
    SLIC_ASSERT(field != nullptr);
    const int num_components = field->getNumComponents();
    SLIC_ASSERT(field->getNumTuples() == num_values);

    const bool invalidType = (field->getType() >= NUMBER_OF_FIELD_TYPES);
    SLIC_ERROR_IF(invalidType, "Field [" << field->getName() << "] has invalid type");

    if(num_components == 1)
    {
      write_scalar_data(field, file);
    }
    else if(num_components == 2 || num_components == 3)
    {
      write_vector_data(field, file);
    }
    else if(num_components > 3)
    {
      write_multidim_data(field, file);
    }
    else
    {
      SLIC_WARNING("Field has an improper number of components.");
    }
  }
}

} /* namespace internal */

//------------------------------------------------------------------------------
int write_vtk(const Mesh* mesh, const std::string& file_path)
{
  SLIC_ASSERT(mesh != nullptr);
  int mesh_type = mesh->getMeshType();

  std::ofstream file(file_path.c_str());
  if(!file.good())
  {
    SLIC_WARNING("Could not open file at path " << file_path);
    return -1;
  }

  /* Write the VTK header */
  file << "# vtk DataFile Version 3.0\n";
  file << "Mesh generated by axom::mint::write_vtk\n";
  file << "ASCII\n";

  /* Write out the mesh node and cell coordinates. */
  if(mesh_type == mint::UNSTRUCTURED_MESH || mesh_type == mint::PARTICLE_MESH)
  {
    file << "DATASET UNSTRUCTURED_GRID\n";
    internal::write_points(mesh, file);
    internal::write_cells(mesh, file);
  }
  else if(mesh_type == mint::STRUCTURED_CURVILINEAR_MESH)
  {
    file << "DATASET STRUCTURED_GRID\n";
    const StructuredMesh* struc_mesh = dynamic_cast<const StructuredMesh*>(mesh);
    internal::write_dimensions(struc_mesh, file);
    internal::write_points(struc_mesh, file);
  }
  else if(mesh_type == mint::STRUCTURED_RECTILINEAR_MESH)
  {
    file << "DATASET RECTILINEAR_GRID\n";
    const RectilinearMesh* rect_mesh = dynamic_cast<const RectilinearMesh*>(mesh);
    internal::write_rectilinear_mesh(rect_mesh, file);
  }
  else if(mesh_type == mint::STRUCTURED_UNIFORM_MESH)
  {
    file << "DATASET STRUCTURED_POINTS\n";
    const UniformMesh* uniform_mesh = dynamic_cast<const UniformMesh*>(mesh);
    internal::write_uniform_mesh(uniform_mesh, file);
  }
  else
  {
    SLIC_WARNING("Mesh does not have a proper type (" << mesh_type << ") "
                                                      << "write aborted.");
    file.close();
    remove(file_path.c_str());
    return -1;
  }

  /* Write out the node data if any. */
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const FieldData* node_data = mesh->getFieldData(mint::NODE_CENTERED);
  if(node_data->getNumFields() > 0)
  {
    fmt::print(file, "POINT_DATA {}\n", num_nodes);
    internal::write_data(node_data, num_nodes, file);
  }

  /* Write out the cell data if any. */
  if(mesh->getMeshType() != mint::PARTICLE_MESH)
  {
    const IndexType num_cells = mesh->getNumberOfCells();
    const FieldData* cell_data = mesh->getFieldData(mint::CELL_CENTERED);
    if(cell_data->getNumFields() > 0)
    {
      fmt::print(file, "CELL_DATA {}\n", num_cells);
      internal::write_data(cell_data, num_cells, file);
    }
  }

  file.close();
  return 0;
}

//------------------------------------------------------------------------------
int write_vtk(mint::FiniteElement& fe, const std::string& file_path)
{
  if(file_path.empty())
  {
    return -1;
  }

  std::ofstream ofs(file_path.c_str());
  if(!ofs.is_open())
  {
    SLIC_WARNING("Could not open file at path " << file_path);
    return -1;
  }

  const bool zero_copy = true;
  const CellType cell_type = fe.getCellType();
  const int ndims = fe.getPhysicalDimension();
  const int nnodes = fe.getNumNodes();
  numerics::Matrix<double> nodes(ndims, nnodes, fe.getPhysicalNodes(), zero_copy);

  // write VTK header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " FiniteElement\n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";
  ofs << "POINTS " << nnodes << " double\n";

  // write the cell coordinates
  for(int i = 0; i < nnodes; ++i)
  {
    const double* pt = nodes.getColumn(i);
    const double x = pt[0];
    const double y = (ndims > 1) ? pt[1] : 0.0;
    const double z = (ndims > 2) ? pt[2] : 0.0;

    fmt::print(ofs, "{} {} {}\n", x, y, z);

  }  // END for all nodes

  // write cell connectivity
  ofs << "CELLS 1 " << nnodes + 1 << std::endl;
  ofs << nnodes << " ";
  for(int i = 0; i < nnodes; ++i)
  {
    ofs << i << " ";
  }
  ofs << std::endl;

  // write cell type information
  ofs << "CELL_TYPES 1\n";
  ofs << mint::getCellInfo(cell_type).vtk_type << std::endl;

  // close the file
  ofs.close();

  // return success
  return 0;
}

} /* namespace mint */
} /* namespace axom */
