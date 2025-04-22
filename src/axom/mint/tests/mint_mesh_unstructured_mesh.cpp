// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/config.hpp"                /* for mint defintions */
#include "axom/mint/mesh/UnstructuredMesh.hpp" /* for mint::Array */

#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic.hpp"

#include "mint_test_utilities.hpp" /* for create_mesh() */

#include "gtest/gtest.h"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
#endif

// C/C++ includes
#include <algorithm>
#include <string>
#include <sstream>
#include <unordered_map>

namespace axom
{
namespace mint
{
constexpr double PI = 3.14159265358979323846;
constexpr double E = 2.71828182845904523536;
const char IGNORE_OUTPUT[] = ".*";

#ifdef AXOM_MINT_USE_SIDRE
sidre::DataStore* ds = nullptr;
#endif

enum class CheckOption
{
  NONE,
  APPEND,
  SET
};

namespace internal
{
/*!
 * \brief Return the name of the field with the given association.
 *
 * \param [in] mesh the mesh in question.
 * \param [in] association the association.
 *
 * \pre mesh != nullptr
 * \pre association == NODE_CENTERED or association == CELL_CENTERED
 */
std::string getFieldName(const Mesh* mesh, int association)
{
  std::string field_name = "f1";
  if(association == CELL_CENTERED)
  {
    field_name = "f2";
  }

  if(mesh->hasSidreGroup())
  {
#ifdef AXOM_MINT_USE_SIDRE
    field_name = mesh->getTopologyName() + "_" + field_name;
#endif
  }

  return field_name;
}

/*!
 * \brief Create a field on the mesh with the given association.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] association the association.
 *
 * \pre mesh != nullptr
 * \pre association == NODE_CENTERED or association == CELL_CENTERED
 */
void createField(Mesh* mesh, int association)
{
  ASSERT_FALSE(mesh->isExternal());
  const std::string field_name = getFieldName(mesh, association);
  if(association == NODE_CENTERED)
  {
    mesh->createField<double>(field_name, NODE_CENTERED, 2);
    return;
  }
  mesh->createField<double>(field_name, CELL_CENTERED, 2);
}

/*!
 * \brief Create an external field on the mesh with the given association.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] association the association.
 * \param [in] capacity the capacity of the field.
 *
 * \pre mesh != nullptr
 * \pre association == NODE_CENTERED or association == CELL_CENTERED
 */
void createExternalField(Mesh* mesh, int association, IndexType capacity)
{
  ASSERT_TRUE(mesh->isExternal());
  double* data = new double[2 * capacity];
  const std::string field_name = getFieldName(mesh, association);
  if(association == NODE_CENTERED)
  {
    mesh->createField<double>(field_name, NODE_CENTERED, data, 2, capacity);
    return;
  }
  mesh->createField<double>(field_name, CELL_CENTERED, data, 2, capacity);
}

/*!
 * \brief Return the field variable of the field with the given association.
 *
 * \param [in] mesh the mesh in question.
 * \param [in] association the association.
 *
 * \pre mesh != nullptr
 * \pre association == NODE_CENTERED or association == CELL_CENTERED
 */
FieldVariable<double>* getFieldVar(Mesh* mesh, int association)
{
  const std::string field_name = getFieldName(mesh, association);
  Field* f = const_cast<Field*>(mesh->getFieldData(association)->getField(field_name));
  return dynamic_cast<FieldVariable<double>*>(f);
}

/*!
 * \brief Return the value of the field for the given tuple and component.
 *
 * \param [in] i the tuple index.
 * \param [in] j the component index.
 */
double getFieldValue(IndexType i, IndexType j) { return PI * i + E * j; }

/*!
 * \brief Return the new value of the field for the given tuple and component.
 *
 * \param [in] i the tuple index.
 * \param [in] j the component index.
 */
double getNewFieldValue(IndexType i, IndexType j) { return E * i + PI * j; }

/*!
 * \brief Set the given tuple of the given FieldVariable.
 *
 * \param [in/out] fv the FieldVariable in question.
 * \param [in] i the tuple index.
 * \param [in] final_i the index of the tuple in question once check_fields is
 *  called. Only needs to be specified when doing insertions.
 */
void setFieldTuple(FieldVariable<double>* fv,
                   IndexType n_tuples,
                   IndexType i,
                   IndexType final_i = USE_DEFAULT)
{
  if(final_i == USE_DEFAULT)
  {
    final_i = i;
  }

  EXPECT_EQ(n_tuples, fv->getNumTuples());

  double* data = fv->getFieldVariablePtr();
  const IndexType num_components = fv->getNumComponents();
  for(IndexType j = 0; j < num_components; ++j)
  {
    data[i * num_components + j] = getFieldValue(final_i, j);
  }
}

/*!
 * \brief Set the given tuple of the given FieldVariable to a new value.
 *
 * \param [in/out] fv the FieldVariable in question.
 * \param [in] i the tuple index.
 */
void setNewFieldTuple(FieldVariable<double>* fv, IndexType n_tuples, IndexType i)
{
  EXPECT_EQ(n_tuples, fv->getNumTuples());

  double* data = fv->getFieldVariablePtr();
  const IndexType num_components = fv->getNumComponents();
  for(IndexType j = 0; j < num_components; ++j)
  {
    data[i * num_components + j] = getNewFieldValue(i, j);
  }
}

/*!
 * \brief Check that the values of the given field are as expected.
 *
 * \param [in] fv the Mesh in question.
 * \param [in] assoc the association to check.
 * \param [in] newValues whether the field is expected to have new values.
 */
void check_fields(const Mesh* mesh, int assoc, bool newValues = false)
{
  double (*fieldValue)(IndexType, IndexType) = &getFieldValue;
  if(newValues)
  {
    fieldValue = &getNewFieldValue;
  }

  IndexType expected_n_tuples = mesh->getNumberOfNodes();
  IndexType expected_capacity = mesh->getNodeCapacity();
  if(assoc == CELL_CENTERED)
  {
    expected_n_tuples = mesh->getNumberOfCells();
    expected_capacity = mesh->getCellCapacity();
  }

  const FieldData* fd = mesh->getFieldData(assoc);
  ASSERT_TRUE(fd->checkConsistency(expected_n_tuples, expected_capacity));

  const std::string field_name = getFieldName(mesh, assoc);
  ASSERT_TRUE(fd->hasField(field_name));

  IndexType num_tuples = -1;
  IndexType num_components = -1;
  const double* data = fd->getFieldPtr<double>(field_name, num_tuples, num_components);

  ASSERT_EQ(num_tuples, expected_n_tuples);
  for(IndexType i = 0; i < num_tuples; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_DOUBLE_EQ(data[i * num_components + j], fieldValue(i, j));
    }
  }
}

/*!
 * \brief Return a pointer to a new UnstructuredMesh of single topology.
 *
 * \param [in] ndims the dimension of the mesh to create.
 * \param [in] cell_type the cell type of the mesh.
 * \param [in] cell_capacity space to allocate for cells.
 * \param [in] node_capacity space to allocate for nodes.
 */
UnstructuredMesh<SINGLE_SHAPE>* createExternalSingle(int ndims,
                                                     CellType cell_type,
                                                     IndexType cell_capacity,
                                                     IndexType node_capacity)
{
  double* x = new double[node_capacity];
  double* y = nullptr;
  double* z = nullptr;
  if(ndims > 1)
  {
    y = new double[node_capacity];
  }
  if(ndims > 2)
  {
    z = new double[node_capacity];
  }

  IndexType connec_capacity = getCellInfo(cell_type).num_nodes * cell_capacity;
  IndexType* connectivity = new IndexType[connec_capacity];

  return new UnstructuredMesh<SINGLE_SHAPE>(cell_type,
                                            0,
                                            cell_capacity,
                                            connectivity,
                                            0,
                                            node_capacity,
                                            x,
                                            y,
                                            z);
}

/*!
 * \brief Return a pointer to a new UnstructuredMesh of mixed topology.
 *
 * \param [in] ndims the dimension of the mesh to create.
 * \param [in] cell_capacity space to allocate for cells.
 * \param [in] node_capacity space to allocate for nodes.
 * \param [in] connec_capacity space to allocate for the connectivity array.
 */
UnstructuredMesh<MIXED_SHAPE>* createExternalMixed(int ndims,
                                                   IndexType cell_capacity,
                                                   IndexType node_capacity,
                                                   IndexType connec_capacity)
{
  double* x = new double[node_capacity];
  double* y = nullptr;
  double* z = nullptr;
  if(ndims > 1)
  {
    y = new double[node_capacity];
  }
  if(ndims > 2)
  {
    z = new double[node_capacity];
  }

  IndexType* connectivity = new IndexType[connec_capacity];
  IndexType* offsets = new IndexType[cell_capacity + 1];
  CellType* types = new CellType[cell_capacity];

  return new UnstructuredMesh<MIXED_SHAPE>(0,
                                           cell_capacity,
                                           connec_capacity,
                                           connectivity,
                                           offsets,
                                           types,
                                           0,
                                           node_capacity,
                                           x,
                                           y,
                                           z);
}

/*!
 * \brief Delete the given external mesh and its associated arrays.
 *
 * \param [in/out] mesh the mesh to delete.
 */
template <Topology TOPO>
void deleteExternalMesh(UnstructuredMesh<TOPO>*& mesh)
{
  ASSERT_TRUE(mesh->isExternal());

  const int ndims = mesh->getDimension();
  double* x = mesh->getCoordinateArray(X_COORDINATE);
  double* y = nullptr;
  double* z = nullptr;
  if(ndims > 1)
  {
    y = mesh->getCoordinateArray(Y_COORDINATE);
  }
  if(ndims > 2)
  {
    z = mesh->getCoordinateArray(Z_COORDINATE);
  }

  const IndexType* connectivity = mesh->getCellNodesArray();
  const IndexType* offsets = mesh->getCellNodesOffsetsArray();
  const CellType* types = mesh->getCellTypesArray();

  double* f1 = nullptr;
  double* f2 = nullptr;
  const std::string node_field_name = getFieldName(mesh, NODE_CENTERED);
  const std::string cell_field_name = getFieldName(mesh, CELL_CENTERED);
  if(mesh->hasField(node_field_name, NODE_CENTERED))
  {
    f1 = mesh->template getFieldPtr<double>(node_field_name, NODE_CENTERED);
  }
  if(mesh->hasField(cell_field_name, CELL_CENTERED))
  {
    f2 = mesh->template getFieldPtr<double>(cell_field_name, CELL_CENTERED);
  }

  delete mesh;
  mesh = nullptr;
  delete[] x;

  if(ndims > 1)
  {
    delete[] y;
  }
  if(ndims > 2)
  {
    delete[] z;
  }

  delete[] connectivity;
  delete[] offsets;
  delete[] types;

  if(f1 != nullptr)
  {
    delete[] f1;
  }
  if(f2 != nullptr)
  {
    delete[] f2;
  }
}

/*!
 * \brief Delete the given external mesh and create a copy from it's external
 *  arrays.
 *
 * \param [in/out] mesh the mesh to delete, it is replaced with the copy.
 */
/// @{

void deleteAndDuplicateExternalMesh(UnstructuredMesh<SINGLE_SHAPE>*& mesh)
{
  ASSERT_TRUE(mesh->isExternal());

  const int ndims = mesh->getDimension();
  const CellType cell_type = mesh->getCellType();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  IndexType* connectivity = mesh->getCellNodesArray();
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();
  double* x = mesh->getCoordinateArray(X_COORDINATE);
  double* y = nullptr;
  double* z = nullptr;
  if(ndims > 1)
  {
    y = mesh->getCoordinateArray(Y_COORDINATE);
  }
  if(ndims > 2)
  {
    z = mesh->getCoordinateArray(Z_COORDINATE);
  }

  double* f1 = nullptr;
  double* f2 = nullptr;
  const std::string node_field_name = getFieldName(mesh, NODE_CENTERED);
  const std::string cell_field_name = getFieldName(mesh, CELL_CENTERED);
  if(mesh->hasField(node_field_name, NODE_CENTERED))
  {
    f1 = mesh->getFieldPtr<double>(node_field_name, NODE_CENTERED);
  }
  if(mesh->hasField(cell_field_name, CELL_CENTERED))
  {
    f2 = mesh->getFieldPtr<double>(cell_field_name, CELL_CENTERED);
  }

  delete mesh;
  mesh = new UnstructuredMesh<SINGLE_SHAPE>(cell_type,
                                            n_cells,
                                            cell_capacity,
                                            connectivity,
                                            n_nodes,
                                            node_capacity,
                                            x,
                                            y,
                                            z);

  if(f1 != nullptr)
  {
    mesh->createField<double>(node_field_name, NODE_CENTERED, f1, 2, node_capacity);
  }
  if(f2 != nullptr)
  {
    mesh->createField<double>(cell_field_name, CELL_CENTERED, f2, 2, cell_capacity);
  }

  EXPECT_EQ(ndims, mesh->getDimension());
  EXPECT_EQ(cell_type, mesh->getCellType());
  EXPECT_EQ(n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cell_capacity, mesh->getCellCapacity());
  EXPECT_EQ(connectivity, mesh->getCellNodesArray());
  EXPECT_EQ(n_nodes, mesh->getNumberOfNodes());
  EXPECT_EQ(node_capacity, mesh->getNodeCapacity());
  EXPECT_EQ(x, mesh->getCoordinateArray(X_COORDINATE));
  if(ndims > 1)
  {
    EXPECT_EQ(y, mesh->getCoordinateArray(Y_COORDINATE));
  }
  if(ndims > 2)
  {
    EXPECT_EQ(z, mesh->getCoordinateArray(Z_COORDINATE));
  }
}

void deleteAndDuplicateExternalMesh(UnstructuredMesh<MIXED_SHAPE>*& mesh)
{
  ASSERT_TRUE(mesh->isExternal());

  const int ndims = mesh->getDimension();
  const CellType cell_type = mesh->getCellType();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  const IndexType connec_capacity = mesh->getCellNodesCapacity();
  IndexType* connectivity = mesh->getCellNodesArray();
  IndexType* offsets = const_cast<IndexType*>(mesh->getCellNodesOffsetsArray());
  CellType* types = const_cast<CellType*>(mesh->getCellTypesArray());
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();
  double* x = mesh->getCoordinateArray(X_COORDINATE);
  double* y = nullptr;
  double* z = nullptr;
  if(ndims > 1)
  {
    y = mesh->getCoordinateArray(Y_COORDINATE);
  }
  if(ndims > 2)
  {
    z = mesh->getCoordinateArray(Z_COORDINATE);
  }

  double* f1 = nullptr;
  double* f2 = nullptr;
  const std::string node_field_name = getFieldName(mesh, NODE_CENTERED);
  const std::string cell_field_name = getFieldName(mesh, CELL_CENTERED);
  if(mesh->hasField(node_field_name, NODE_CENTERED))
  {
    f1 = mesh->getFieldPtr<double>(node_field_name, NODE_CENTERED);
  }
  if(mesh->hasField(cell_field_name, CELL_CENTERED))
  {
    f2 = mesh->getFieldPtr<double>(cell_field_name, CELL_CENTERED);
  }

  delete mesh;
  mesh = new UnstructuredMesh<MIXED_SHAPE>(n_cells,
                                           cell_capacity,
                                           connec_capacity,
                                           connectivity,
                                           offsets,
                                           types,
                                           n_nodes,
                                           node_capacity,
                                           x,
                                           y,
                                           z);

  if(f1 != nullptr)
  {
    mesh->createField<double>(node_field_name, NODE_CENTERED, f1, 2, node_capacity);
  }
  if(f2 != nullptr)
  {
    mesh->createField<double>(cell_field_name, CELL_CENTERED, f2, 2, cell_capacity);
  }

  EXPECT_EQ(ndims, mesh->getDimension());
  EXPECT_EQ(cell_type, mesh->getCellType());
  EXPECT_EQ(n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cell_capacity, mesh->getCellCapacity());
  EXPECT_EQ(connectivity, mesh->getCellNodesArray());
  EXPECT_EQ(offsets, mesh->getCellNodesOffsetsArray());
  EXPECT_EQ(types, mesh->getCellTypesArray());
  EXPECT_EQ(n_nodes, mesh->getNumberOfNodes());
  EXPECT_EQ(node_capacity, mesh->getNodeCapacity());
  EXPECT_EQ(x, mesh->getCoordinateArray(X_COORDINATE));
  if(ndims > 1)
  {
    EXPECT_EQ(y, mesh->getCoordinateArray(Y_COORDINATE));
  }
  if(ndims > 2)
  {
    EXPECT_EQ(z, mesh->getCoordinateArray(Z_COORDINATE));
  }
}
/// @}

#ifdef AXOM_MINT_USE_SIDRE

/*!
 * \brief Delete the given sidre mesh and create a copy from it's group.
 *
 * \param [in/out] mesh the mesh to delete, it is replaced with the copy.
 */
template <Topology TOPO>
void deleteAndDuplicateSidreMesh(UnstructuredMesh<TOPO>*& mesh)
{
  ASSERT_TRUE(mesh->isInSidre());

  const int ndims = mesh->getDimension();
  const CellType cell_type = mesh->getCellType();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();

  sidre::Group* group = mesh->getSidreGroup();
  const std::string topo = mesh->getTopologyName();

  delete mesh;
  mesh = new UnstructuredMesh<TOPO>(group, topo);

  EXPECT_EQ(ndims, mesh->getDimension());
  EXPECT_EQ(cell_type, mesh->getCellType());
  EXPECT_EQ(n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cell_capacity, mesh->getCellCapacity());
  EXPECT_EQ(n_nodes, mesh->getNumberOfNodes());
  EXPECT_EQ(node_capacity, mesh->getNodeCapacity());
}

#endif

/*!
 * \brief Fill the given array with single topology meshes.
 *
 * \param [in/out] meshes the array to fill.
 * \param [in] n_nodes the nodal capacity for the external meshes.
 * \param [in] n_cells the cell capacity for the external meshes.
 * \param [in] connec_size the connectivity capacity for the external meshes.
 * \param [in] node_field whether to create a node centered field.
 * \param [in] cell_field whether to create a cell centered field.
 *
 * \note if using sidre 9 meshes are created otherwise only 6.
 */
/// @{

void createMeshes(UnstructuredMesh<SINGLE_SHAPE>** meshes,
                  IndexType n_nodes,
                  IndexType n_cells,
                  bool node_field,
                  bool cell_field)
{
#ifdef AXOM_MINT_USE_SIDRE
  SLIC_ERROR_IF(ds != nullptr, "Did not free DataStore.");
  ds = new sidre::DataStore();
  sidre::Group* root = ds->getRoot();
#endif

  int cur_mesh = 0;
  for(int dim = 1; dim <= 3; ++dim)
  {
    meshes[cur_mesh] = new UnstructuredMesh<SINGLE_SHAPE>(dim, QUAD);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;

    meshes[cur_mesh] = internal::createExternalSingle(dim, QUAD, n_cells, n_nodes);
    if(node_field)
    {
      internal::createExternalField(meshes[cur_mesh], NODE_CENTERED, n_nodes);
    }
    if(cell_field)
    {
      internal::createExternalField(meshes[cur_mesh], CELL_CENTERED, n_cells);
    }
    cur_mesh++;

#ifdef AXOM_MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string(dim);
    const std::string coordset = "c" + std::to_string(dim);
    meshes[cur_mesh] = new UnstructuredMesh<SINGLE_SHAPE>(dim, QUAD, root, topo, coordset);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;
#endif
  }
}

void createMeshes(UnstructuredMesh<MIXED_SHAPE>** meshes,
                  IndexType n_nodes,
                  IndexType n_cells,
                  IndexType connec_size,
                  bool node_field,
                  bool cell_field)
{
#ifdef AXOM_MINT_USE_SIDRE
  SLIC_ERROR_IF(ds != nullptr, "Did not free DataStore.");
  ds = new sidre::DataStore();
  sidre::Group* root = ds->getRoot();
#endif

  int cur_mesh = 0;
  for(int dim = 1; dim <= 3; ++dim)
  {
    meshes[cur_mesh] = new UnstructuredMesh<MIXED_SHAPE>(dim);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;

    meshes[cur_mesh] = internal::createExternalMixed(dim, n_cells, n_nodes, connec_size);
    if(node_field)
    {
      internal::createExternalField(meshes[cur_mesh], NODE_CENTERED, n_nodes);
    }
    if(cell_field)
    {
      internal::createExternalField(meshes[cur_mesh], CELL_CENTERED, n_cells);
    }
    cur_mesh++;

#ifdef AXOM_MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string(dim);
    const std::string coordset = "c" + std::to_string(dim);
    meshes[cur_mesh] = new UnstructuredMesh<MIXED_SHAPE>(dim, root, topo, coordset);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;
#endif
  }
}

/// @}

/*!
 * \brief Fill the given array with single topology meshes created to be used
 *  with the resize tests.
 *
 * \param [in/out] meshes the array to fill.
 * \param [in] n_nodes the nodal capacity.
 * \param [in] n_cells the cell capacity.
 * \param [in] connec_size the connectivity capacity.
 * \param [in] node_field whether to create a node centered field.
 * \param [in] cell_field whether to create a cell centered field.
 *
 * \note if using sidre 6 meshes are created otherwise only 4.
 */
/// @{

void createMeshesForResize(UnstructuredMesh<SINGLE_SHAPE>** meshes,
                           IndexType n_nodes,
                           IndexType n_cells,
                           bool node_field,
                           bool cell_field)
{
#ifdef AXOM_MINT_USE_SIDRE
  SLIC_ERROR_IF(ds != nullptr, "Did not free DataStore.");
  ds = new sidre::DataStore();
  sidre::Group* root = ds->getRoot();
#endif

  int cur_mesh = 0;
  for(int dim = 1; dim <= 3; ++dim)
  {
    meshes[cur_mesh] = new UnstructuredMesh<SINGLE_SHAPE>(dim, QUAD, n_nodes, n_cells);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;

#ifdef AXOM_MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string(dim);
    const std::string coordset = "c" + std::to_string(dim);
    meshes[cur_mesh] =
      new UnstructuredMesh<SINGLE_SHAPE>(dim, QUAD, root, topo, coordset, n_nodes, n_cells);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;
#endif
  }
}

void createMeshesForResize(UnstructuredMesh<MIXED_SHAPE>** meshes,
                           IndexType n_nodes,
                           IndexType n_cells,
                           bool node_field,
                           bool cell_field)
{
#ifdef AXOM_MINT_USE_SIDRE
  SLIC_ERROR_IF(ds != nullptr, "Did not free DataStore.");
  ds = new sidre::DataStore();
  sidre::Group* root = ds->getRoot();
#endif

  int cur_mesh = 0;
  for(int dim = 1; dim <= 3; ++dim)
  {
    meshes[cur_mesh] = new UnstructuredMesh<MIXED_SHAPE>(dim, n_nodes, n_cells);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;

#ifdef AXOM_MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string(dim);
    const std::string coordset = "c" + std::to_string(dim);
    meshes[cur_mesh] = new UnstructuredMesh<MIXED_SHAPE>(dim, root, topo, coordset, n_nodes, n_cells);
    if(node_field)
    {
      internal::createField(meshes[cur_mesh], NODE_CENTERED);
    }
    if(cell_field)
    {
      internal::createField(meshes[cur_mesh], CELL_CENTERED);
    }
    cur_mesh++;
#endif
  }
}

/// @}

/*!
 * \brief Create a UniformMesh and its equivalent
 * UnstructuredMesh<SINGLE_SHAPE> to be used with the face relation tests.
 *
 * \param [out] derived_mesh the UnstructuredMesh<SINGLE_SHAPE>
 * \param [out] source_mesh the UniformMesh
 * \param [in] resolution the resolution in the first dimension.
 *             Resolution in other dimensions is proportional to this
 *             parameter.
 * \param [in] dimension the dimensionality of the meshes (2 or 3).
 *
 * This routine uses create_mesh< UNSTRUCTURED_MESH >() from
 * mint_test_utilities.hpp to create UnstructuredMesh objects with
 * correct node and cell connectivity.  Thus, none of the objects
 * created by this routine use Sidre or external storage.
 */
void createMeshesForFace(UnstructuredMesh<SINGLE_SHAPE>*& test_mesh,
                         UniformMesh*& source_mesh,
                         IndexType resolution,
                         int dimension)
{
  const double lo[] = {0.0, 0.0, 0.0};
  const double hi[] = {2.0, 2.0, 2.0};

  constexpr double Y_RES_FACTOR = 0.75;
  constexpr double Z_RES_FACTOR = 1.4;

  switch(dimension)
  {
  case 2:
    source_mesh = new UniformMesh(lo, hi, resolution, (IndexType)(Y_RES_FACTOR * resolution));
    break;
  default:
    EXPECT_TRUE(dimension == 3);
    source_mesh = new UniformMesh(lo,
                                  hi,
                                  resolution,
                                  (IndexType)(Y_RES_FACTOR * resolution),
                                  (IndexType)(Z_RES_FACTOR * resolution));
  }

  Mesh* the_mesh = create_mesh<UNSTRUCTURED_MESH, SINGLE_SHAPE>(*source_mesh);
  test_mesh = static_cast<UnstructuredMesh<SINGLE_SHAPE>*>(the_mesh);
}

/*!
 * \brief Free the given array of meshes.
 *
 * \param [in/out] meshes the array of meshes.
 * \param [in] n_meshes the nuber of meshes in the array.
 */
template <Topology TOPO>
void deleteMeshes(UnstructuredMesh<TOPO>** meshes, int n_meshes)
{
  for(int i = 0; i < n_meshes; ++i)
  {
    if(meshes[i]->isExternal())
    {
      internal::deleteExternalMesh(meshes[i]);
    }
    else
    {
      delete meshes[i];
      meshes[i] = nullptr;
    }
  }

#ifdef AXOM_MINT_USE_SIDRE
  delete ds;
  ds = nullptr;
#endif
}

/*!
 * \brief Return the value of the nodal coordinate of the given node
 *  and dimension.
 *
 * \param [in] ndims the number of dimensions.
 * \param [in] i the nodal index.
 * \param [in] dim the dimension index.
 */
double getCoordValue(IndexType ndims, IndexType i, int dim) { return (ndims * i + dim) * PI; }

/*!
 * \brief Return the value of a new nodal coordinate of the given node
 *  and dimension.
 *
 * \param [in] ndims the number of dimensions.
 * \param [in] i the nodal index.
 * \param [in] dim the dimension index.
 */
double getNewCoordValue(IndexType ndims, IndexType i, int dim) { return (ndims * i + dim) * E; }

/*!
 * \brief Check that the nodal coordinates are correct.
 *
 * \param [in] mesh the mesh to check.
 * \param [in] newValues iff true the expected values are those corresponding to
 *  getNewCoordValue.
 */
template <Topology TOPO>
void check_append_nodes(const UnstructuredMesh<TOPO>* mesh, bool newValues = false)
{
  IndexType n_nodes = mesh->getNumberOfNodes();

  const int ndims = mesh->getDimension();
  const double* x = mesh->getCoordinateArray(X_COORDINATE);

  double (*coordValue)(IndexType, IndexType, int) = &getCoordValue;
  if(newValues)
  {
    coordValue = &getNewCoordValue;
  }

  /* Check using pointers */
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    EXPECT_EQ(x[i], coordValue(ndims, i, 0));
  }
  if(ndims > 1)
  {
    const double* y = mesh->getCoordinateArray(Y_COORDINATE);
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      EXPECT_EQ(y[i], coordValue(ndims, i, 1));
    }
  }
  if(ndims > 2)
  {
    const double* z = mesh->getCoordinateArray(Z_COORDINATE);
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      EXPECT_EQ(z[i], coordValue(ndims, i, 2));
    }
  }

  /* Check using getNode */
  double coords[3];
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    mesh->getNode(i, coords);
    for(int dim = 0; dim < ndims; ++dim)
    {
      EXPECT_EQ(coords[dim], coordValue(ndims, i, dim));
    }
  }

  /* Check using getNodeCoordinate */
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    for(int dim = 0; dim < ndims; ++dim)
    {
      EXPECT_EQ(mesh->getNodeCoordinate(i, dim), coordValue(ndims, i, dim));
    }
  }
}

/*!
 * \brief Append a single node to the mesh multiple times in a row.
 *
 * \param [in/out] mesh the mesh to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
template <Topology TOPO>
void append_node_single(UnstructuredMesh<TOPO>* mesh, IndexType n_nodes)
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, NODE_CENTERED);

  if(ndims == 1)
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      mesh->appendNode(getCoordValue(ndims, cur_n_nodes, 0));
      EXPECT_EQ(++cur_n_nodes, mesh->getNumberOfNodes());
      setFieldTuple(fv, cur_n_nodes, cur_n_nodes - 1);
    }
  }
  else if(ndims == 2)
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      mesh->appendNode(getCoordValue(ndims, cur_n_nodes, 0), getCoordValue(ndims, cur_n_nodes, 1));
      EXPECT_EQ(++cur_n_nodes, mesh->getNumberOfNodes());
      setFieldTuple(fv, cur_n_nodes, cur_n_nodes - 1);
    }
  }
  else
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      mesh->appendNode(getCoordValue(ndims, cur_n_nodes, 0),
                       getCoordValue(ndims, cur_n_nodes, 1),
                       getCoordValue(ndims, cur_n_nodes, 2));
      EXPECT_EQ(++cur_n_nodes, mesh->getNumberOfNodes());
      setFieldTuple(fv, cur_n_nodes, cur_n_nodes - 1);
    }
  }
}

/*!
 * \brief Append multiple nodes at once using the array of structs layout.
 *
 * \param [in/out] mesh the mesh to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
template <Topology TOPO>
void append_node_structs(UnstructuredMesh<TOPO>* mesh, IndexType n_nodes)
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, NODE_CENTERED);
  double* coords = new double[ndims * n_nodes];

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    for(IndexType dim = 0; dim < ndims; ++dim)
    {
      coords[ndims * i + dim] = getCoordValue(ndims, i + cur_n_nodes, dim);
    }
  }

  mesh->appendNodes(coords, n_nodes);
  cur_n_nodes += n_nodes;
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    setFieldTuple(fv, cur_n_nodes, cur_n_nodes - 1 - i);
  }

  delete[] coords;
}

/*!
 * \brief Append multiple nodes at once using the struct of arrays layout.
 *
 * \param [in/out] mesh the mesh to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
template <Topology TOPO>
void append_node_arrays(UnstructuredMesh<TOPO>* mesh, IndexType n_nodes)
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, NODE_CENTERED);
  double* coords = new double[ndims * n_nodes];

  for(IndexType dim = 0; dim < ndims; ++dim)
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      coords[dim * n_nodes + i] = getCoordValue(ndims, i + cur_n_nodes, dim);
    }
  }

  if(ndims == 1)
  {
    mesh->appendNodes(coords, n_nodes);
    cur_n_nodes += n_nodes;
    EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());
  }
  else if(ndims == 2)
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    mesh->appendNodes(x, y, n_nodes);
    cur_n_nodes += n_nodes;
    EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());
  }
  else
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    const double* z = coords + 2 * n_nodes;
    mesh->appendNodes(x, y, z, n_nodes);
    cur_n_nodes += n_nodes;
    EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());
  }

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    setFieldTuple(fv, cur_n_nodes, cur_n_nodes - 1 - i);
  }

  delete[] coords;
}

/*!
 * \brief Append nodes to the mesh using all the various methods.
 *
 * \param [in/out] mesh the mesh to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
template <Topology TOPO>
void append_nodes(UnstructuredMesh<TOPO>* mesh, IndexType n_nodes)
{
  ASSERT_EQ(n_nodes % 3, 0);

  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  ASSERT_EQ(cur_n_nodes, 0);
  EXPECT_TRUE(mesh->empty());

  /* Append one node at a time */
  append_node_single(mesh, n_nodes / 3);
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  /* Append multiple nodes at once using the array of structs layout. */
  append_node_structs(mesh, n_nodes / 3);
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ(cur_n_nodes, 2 * n_nodes / 3);
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  /* Append multiple nodes at once using the struct of arrays layout. */
  append_node_arrays(mesh, n_nodes / 3);
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ(cur_n_nodes, n_nodes);
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  check_append_nodes(mesh);
  check_fields(mesh, NODE_CENTERED);
}

/*!
 * \brief Return the node index of the given cell and vertex.
 *
 * \param [in] cur_cell cell in question.
 * \param [in] vertex the vertex in question.
 */
IndexType getCellConnecValue(IndexType cur_cell, IndexType vertex)
{
  return cur_cell * MAX_CELL_NODES + vertex;
}

/*!
 * \brief Return a new node index of the given cell and vertex.
 *
 * \param [in] cur_cell cell in question.
 * \param [in] vertex the vertex in question.
 */
IndexType getNewCellConnecValue(IndexType cur_n_cells, IndexType vertex)
{
  return (cur_n_cells * MAX_CELL_NODES + vertex) * vertex;
}

/*!
 * \brief Return the type of the given cell, only valid for mixed topology.
 *
 * \param [in] cur_cell cell in question.
 */
CellType getCellType(IndexType cur_cell)
{
  CellType type = QUAD;
  if(cur_cell % 2 != 0)
  {
    type = TRIANGLE;
  }
  return type;
}

/*!
 * \brief Return the total size of the connectivity array given the total number
 *  of cells, only valid for mixed topology.
 *
 * \param [in] n_cells the total number of cells.
 */
constexpr IndexType getConnectivitySize(IndexType n_cells) { return (n_cells * (4 + 3)) / 2; }

/*!
 * \brief Check that the cell connectivity is as expected.
 *
 * \param [in] mesh the mesh to check.
 * \param [in] newValues if true the expected values are those given by
 *  getNewCellConnecValue.
 */
/// @{

void check_append_cells(const UnstructuredMesh<SINGLE_SHAPE>* mesh, bool newValues = false)
{
  IndexType (*connecValue)(IndexType, IndexType) = &getCellConnecValue;
  if(newValues)
  {
    connecValue = &getNewCellConnecValue;
  }

  IndexType n_cells = mesh->getNumberOfCells();

  const CellType type = mesh->getCellType();
  const int nodes_per_cell = getCellInfo(type).num_nodes;

  /* Check using pointers */
  const IndexType* connectivity = mesh->getCellNodesArray();
  for(IndexType i = 0; i < n_cells; ++i)
  {
    EXPECT_EQ(type, mesh->getCellType(i));
    EXPECT_EQ(nodes_per_cell, mesh->getNumberOfCellNodes(i));

    for(IndexType j = 0; j < nodes_per_cell; ++j)
    {
      EXPECT_EQ(connectivity[i * nodes_per_cell + j], connecValue(i, j));
    }
  }

  /* Check using getCell */
  IndexType cell[mint::MAX_CELL_NODES];
  for(IndexType i = 0; i < n_cells; ++i)
  {
    mesh->getCellNodeIDs(i, cell);
    for(IndexType j = 0; j < nodes_per_cell; ++j)
    {
      EXPECT_EQ(cell[j], connecValue(i, j));
    }
  }

  /* Check using the other getCell */
  for(IndexType i = 0; i < n_cells; ++i)
  {
    const IndexType* cellPtr = mesh->getCellNodeIDs(i);
    for(IndexType j = 0; j < nodes_per_cell; ++j)
    {
      EXPECT_EQ(cellPtr[j], connecValue(i, j));
    }
  }
}

void check_append_cells(const UnstructuredMesh<MIXED_SHAPE>* mesh, bool newValues = false)
{
  IndexType (*connecValue)(IndexType, IndexType) = &getCellConnecValue;
  if(newValues)
  {
    connecValue = &getNewCellConnecValue;
  }

  IndexType n_cells = mesh->getNumberOfCells();

  /* Check using pointers */
  const IndexType* connectivity = mesh->getCellNodesArray();
  const IndexType* offsets = mesh->getCellNodesOffsetsArray();
  const CellType* types = mesh->getCellTypesArray();
  for(IndexType i = 0; i < n_cells; ++i)
  {
    CellType expected_type = getCellType(i);
    IndexType expected_n_nodes = getCellInfo(expected_type).num_nodes;
    ASSERT_EQ(expected_type, types[i]);
    ASSERT_EQ(expected_type, mesh->getCellType(i));
    ASSERT_EQ(expected_n_nodes, offsets[i + 1] - offsets[i]);
    ASSERT_EQ(expected_n_nodes, mesh->getNumberOfCellNodes(i));

    IndexType offset = offsets[i];
    for(IndexType j = 0; j < expected_n_nodes; ++j)
    {
      EXPECT_EQ(connectivity[offset + j], connecValue(i, j));
    }
  }

  /* Check using getCell */
  IndexType cell[MAX_CELL_NODES];
  for(IndexType i = 0; i < n_cells; ++i)
  {
    mesh->getCellNodeIDs(i, cell);
    IndexType num_nodes = mesh->getNumberOfCellNodes(i);
    for(IndexType j = 0; j < num_nodes; ++j)
    {
      EXPECT_EQ(cell[j], connecValue(i, j));
    }
  }

  /* Check using the other getCell */
  for(IndexType i = 0; i < n_cells; ++i)
  {
    const IndexType* cellNodeIds = mesh->getCellNodeIDs(i);
    IndexType num_nodes = mesh->getNumberOfCellNodes(i);
    for(IndexType j = 0; j < num_nodes; ++j)
    {
      EXPECT_EQ(cellNodeIds[j], connecValue(i, j));
    }
  }
}

/// @}

/*!
 * \brief Copy the connectivity of the given cell into the provided buffer.
 *
 * \param [in] cur_cell the cell in question.
 * \param [in] nodes_per_cell the number of nodes in this cell.
 * \param [out] connec the buffer to copy into.
 */
void getCellConnec(IndexType cur_cell, IndexType nodes_per_cell, IndexType* connec)
{
  for(IndexType i = 0; i < nodes_per_cell; ++i)
  {
    connec[i] = getCellConnecValue(cur_cell, i);
  }
}

/*!
 * \brief Copy the connectivity of the given cell into the provided buffer.
 *
 * \param [in] cur_cell the cell in question.
 * \param [out] connec the buffer to copy into.
 *
 * \note Only to be used with mixed topology.
 */
CellType getCellConnec(IndexType cur_cell, IndexType* connec)
{
  CellType type = getCellType(cur_cell);
  const IndexType nodes_per_cell = getCellInfo(type).num_nodes;
  for(IndexType i = 0; i < nodes_per_cell; ++i)
  {
    connec[i] = getCellConnecValue(cur_cell, i);
  }

  return type;
}

/*!
 * \brief Append a single cell to the mesh multiple times.
 *
 * \param [in/out] mesh the mesh to append to.
 * \param [in] n_cells the number of cells to append.
 */
/// @{

void append_cell_single(UnstructuredMesh<SINGLE_SHAPE>* mesh, IndexType n_cells)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  const IndexType nodes_per_cell = getCellInfo(mesh->getCellType()).num_nodes;
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  IndexType connec[mint::MAX_CELL_NODES];
  for(IndexType i = 0; i < n_cells; ++i)
  {
    getCellConnec(cur_n_cells, nodes_per_cell, connec);
    mesh->appendCell(connec);
    EXPECT_EQ(++cur_n_cells, mesh->getNumberOfCells());
    EXPECT_EQ(cur_n_cells * nodes_per_cell, mesh->getCellNodesSize());
    setFieldTuple(fv, cur_n_cells, cur_n_cells - 1);
  }
}

void append_cell_single(UnstructuredMesh<MIXED_SHAPE>* mesh, IndexType n_cells)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellNodesSize();
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  IndexType connec[mint::MAX_CELL_NODES];
  CellType type;
  for(IndexType i = 0; i < n_cells; ++i)
  {
    type = getCellConnec(cur_n_cells, connec);
    cur_connec_size += getCellInfo(type).num_nodes;
    mesh->appendCell(connec, type);
    EXPECT_EQ(++cur_n_cells, mesh->getNumberOfCells());
    EXPECT_EQ(cur_connec_size, mesh->getCellNodesSize());
    setFieldTuple(fv, cur_n_cells, cur_n_cells - 1);
  }
}

/// @}

/*!
 * \brief Append a multiple cells at once to the mesh.
 *
 * \param [in/out] mesh the mesh to append to.
 * \param [in] n_cells the number of cells to append.
 */
/// @{

void append_cell_multiple(UnstructuredMesh<SINGLE_SHAPE>* mesh, IndexType n_cells)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellNodesSize();
  const IndexType nodes_per_cell = mesh->getNumberOfCellNodes();
  IndexType* connectivity = new IndexType[nodes_per_cell * n_cells];
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  for(IndexType i = 0; i < n_cells; ++i)
  {
    IndexType* cur_connec = connectivity + i * nodes_per_cell;
    getCellConnec(cur_n_cells + i, nodes_per_cell, cur_connec);
  }

  mesh->appendCells(connectivity, n_cells);
  cur_n_cells += n_cells;
  cur_connec_size += n_cells * nodes_per_cell;
  EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cur_connec_size, mesh->getCellNodesSize());

  for(IndexType i = 0; i < n_cells; ++i)
  {
    setFieldTuple(fv, cur_n_cells, cur_n_cells - 1 - i);
  }

  delete[] connectivity;
}

void append_cell_multiple(UnstructuredMesh<MIXED_SHAPE>* mesh, IndexType n_cells)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellNodesSize();
  IndexType* connectivity = new IndexType[getConnectivitySize(n_cells)];
  IndexType* offsets = new IndexType[n_cells + 1];
  CellType* types = new CellType[n_cells];
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  offsets[0] = 0;
  for(IndexType i = 0; i < n_cells; ++i)
  {
    IndexType* cur_connec = connectivity + offsets[i];
    types[i] = getCellConnec(cur_n_cells + i, cur_connec);
    offsets[i + 1] = offsets[i] + getCellInfo(types[i]).num_nodes;
  }

  mesh->appendCells(connectivity, n_cells, offsets, types);
  cur_n_cells += n_cells;
  cur_connec_size += getConnectivitySize(n_cells);
  EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cur_connec_size, mesh->getCellNodesSize());

  for(IndexType i = 0; i < n_cells; ++i)
  {
    setFieldTuple(fv, cur_n_cells, cur_n_cells - 1 - i);
  }

  delete[] connectivity;
  delete[] offsets;
  delete[] types;
}

/// @}

/*!
 * \brief Append cells to the mesh.
 *
 * \param [in/out] mesh the mesh to append to.
 * \param [in] n_cells the number of cells to append.
 */
template <Topology TOPO>
void append_cells(UnstructuredMesh<TOPO>* mesh, IndexType n_cells)
{
  ASSERT_EQ(n_cells % 2, 0);

  IndexType cur_n_cells = mesh->getNumberOfCells();
  ASSERT_EQ(cur_n_cells, 0);

  /* Append cells one at a time */
  append_cell_single(mesh, n_cells / 2);
  cur_n_cells += n_cells / 2;
  EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());

  /* Append cells all at once */
  append_cell_multiple(mesh, n_cells / 2);
  cur_n_cells += n_cells / 2;
  EXPECT_EQ(cur_n_cells, n_cells);
  EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());

  check_append_cells(mesh);
  check_fields(mesh, CELL_CENTERED);
}

/*!
 * \brief Set the nodal coordinates to new values.
 *
 * \param [in/out] mesh the mesh in question.
 */
template <Topology TOPO>
void set_nodes(UnstructuredMesh<TOPO>* mesh)
{
  const int ndims = mesh->getDimension();
  const IndexType n_nodes = mesh->getNumberOfNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, NODE_CENTERED);

  double* x = mesh->getCoordinateArray(X_COORDINATE);
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    x[i] = getNewCoordValue(ndims, i, 0);
    setNewFieldTuple(fv, n_nodes, i);
  }

  if(ndims > 1)
  {
    double* y = mesh->getCoordinateArray(Y_COORDINATE);
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      y[i] = getNewCoordValue(ndims, i, 1);
      setNewFieldTuple(fv, n_nodes, i);
    }
  }
  if(ndims > 2)
  {
    double* z = mesh->getCoordinateArray(Z_COORDINATE);
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      z[i] = getNewCoordValue(ndims, i, 2);
      setNewFieldTuple(fv, n_nodes, i);
    }
  }

  EXPECT_EQ(n_nodes, mesh->getNumberOfNodes());
  check_append_nodes(mesh, true);
  check_fields(mesh, NODE_CENTERED, true);
}

/*!
 * \brief Set the cell connectivity to new values.
 *
 * \param [in/out] mesh the mesh in question.
 */
template <Topology TOPO>
void set_cells(UnstructuredMesh<TOPO>* mesh)
{
  const IndexType n_cells = mesh->getNumberOfCells();
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  for(IndexType i = 0; i < n_cells; ++i)
  {
    IndexType* cur_cell = mesh->getCellNodeIDs(i);
    const IndexType n_nodes = mesh->getNumberOfCellNodes(i);
    for(IndexType j = 0; j < n_nodes; ++j)
    {
      cur_cell[j] = getNewCellConnecValue(i, j);
    }
    setNewFieldTuple(fv, n_cells, i);
  }

  EXPECT_EQ(n_cells, mesh->getNumberOfCells());
  check_append_cells(mesh, true);
  check_fields(mesh, CELL_CENTERED, true);
}

/*!
 * \brief Insert a single node into the mesh.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this node after all
 *  insertions have taken place.
 * \param [in] update iff true will update the cell connectivity.
 */
template <Topology TOPO>
void insert_node_single(UnstructuredMesh<TOPO>* mesh,
                        IndexType pos,
                        IndexType final_pos,
                        bool update = true)
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, NODE_CENTERED);

  if(ndims == 1)
  {
    mesh->insertNode(pos, getCoordValue(ndims, final_pos, 0), update);
  }
  else if(ndims == 2)
  {
    mesh->insertNode(pos,
                     getCoordValue(ndims, final_pos, 0),
                     getCoordValue(ndims, final_pos, 1),
                     update);
  }
  else
  {
    mesh->insertNode(pos,
                     getCoordValue(ndims, final_pos, 0),
                     getCoordValue(ndims, final_pos, 1),
                     getCoordValue(ndims, final_pos, 2),
                     update);
  }

  EXPECT_EQ(++cur_n_nodes, mesh->getNumberOfNodes());
  setFieldTuple(fv, cur_n_nodes, pos, final_pos);
}

/*!
 * \brief Insert multiple nodes into the mesh using the array of structs layout.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] n_nodes the number of nodes to insert.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this node after all
 *  insertions have taken place.
 * \param [in] update iff true will update the cell connectivity.
 */
template <Topology TOPO>
void insert_node_structs(UnstructuredMesh<TOPO>* mesh,
                         IndexType n_nodes,
                         IndexType pos,
                         IndexType final_pos,
                         bool update = true)
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, NODE_CENTERED);
  double* coords = new double[ndims * n_nodes];

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    for(IndexType dim = 0; dim < ndims; ++dim)
    {
      coords[ndims * i + dim] = getCoordValue(ndims, final_pos + i, dim);
    }
  }

  mesh->insertNodes(pos, coords, n_nodes, update);
  cur_n_nodes += n_nodes;
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    setFieldTuple(fv, cur_n_nodes, pos + i, final_pos + i);
  }

  delete[] coords;
}

/*!
 * \brief Insert multiple nodes into the mesh using the struct of arrays layout.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] n_nodes the number of nodes to insert.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this node after all
 *  insertions have taken place.
 * \param [in] update iff true will update the cell connectivity.
 */
template <Topology TOPO>
void insert_node_arrays(UnstructuredMesh<TOPO>* mesh,
                        IndexType n_nodes,
                        IndexType pos,
                        IndexType final_pos,
                        bool update = true)
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, NODE_CENTERED);
  double* coords = new double[ndims * n_nodes];

  for(IndexType dim = 0; dim < ndims; ++dim)
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      coords[dim * n_nodes + i] = getCoordValue(ndims, final_pos + i, dim);
    }
  }

  if(ndims == 1)
  {
    mesh->insertNodes(pos, coords, n_nodes, update);
  }
  else if(ndims == 2)
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    mesh->insertNodes(pos, x, y, n_nodes, update);
  }
  else
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    const double* z = coords + 2 * n_nodes;
    mesh->insertNodes(pos, x, y, z, n_nodes, update);
  }

  cur_n_nodes += n_nodes;
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    setFieldTuple(fv, cur_n_nodes, pos + i, final_pos + i);
  }

  delete[] coords;
}

/*!
 * \brief Insert nodes into the mesh.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] n_nodes the number of nodes to insert.
 */
template <Topology TOPO>
void insert_nodes(UnstructuredMesh<TOPO>* mesh, IndexType n_nodes)
{
  ASSERT_EQ(n_nodes % 9, 0);

  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  ASSERT_EQ(cur_n_nodes, 0);
  EXPECT_TRUE(mesh->empty());

  const IndexType ninth_n_nodes = n_nodes / 9;
  const IndexType third_n_nodes = n_nodes / 3;

  /* The final positions of the nodes to insert. */
  IndexType final_pos[3] = {third_n_nodes - 1, third_n_nodes, 2 * third_n_nodes};

  /* The inital position of the nodes to insert. */
  IndexType insert_positions[3] = {0, 1, 2};

  /* Insert one node at a time */
  for(IndexType round = 0; round < ninth_n_nodes; ++round)
  {
    for(IndexType i = 0; i < 3; ++i)
    {
      insert_node_single(mesh, insert_positions[i], final_pos[i]);
      cur_n_nodes++;
      EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());
    }

    /* The middle insertion increases by two, the back by three. */
    insert_positions[1] += 2;
    insert_positions[2] += 3;

    /* The final position of the front decreases by one while the final position
     *  of the middle and back both increase by one. */
    final_pos[0]--;
    final_pos[1]++;
    final_pos[2]++;
  }
  EXPECT_EQ(cur_n_nodes, third_n_nodes);

  /* Insert multiple nodes at once using the array of structs layout. */
  insert_node_structs(mesh, ninth_n_nodes, 0, ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  insert_node_structs(mesh, ninth_n_nodes, third_n_nodes, 4 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  insert_node_structs(mesh, ninth_n_nodes, 5 * ninth_n_nodes, 7 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  EXPECT_EQ(cur_n_nodes, 2 * third_n_nodes);
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  /* Insert multiple nodes at once using the struct of arrays layout. */
  insert_node_arrays(mesh, ninth_n_nodes, 0, 0);
  cur_n_nodes += ninth_n_nodes;
  insert_node_arrays(mesh, ninth_n_nodes, 5 * ninth_n_nodes, 5 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  insert_node_arrays(mesh, ninth_n_nodes, 8 * ninth_n_nodes, 8 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  EXPECT_EQ(cur_n_nodes, n_nodes);
  EXPECT_EQ(cur_n_nodes, mesh->getNumberOfNodes());

  check_append_nodes(mesh);
  check_fields(mesh, NODE_CENTERED);
}

/*!
 * \brief Insert a single cell into the mesh
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this cell after all
 *  insertions have taken place.
 */
/// @{

void insert_cell_single(UnstructuredMesh<SINGLE_SHAPE>* mesh, IndexType pos, IndexType final_pos)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  const IndexType nodes_per_cell = mesh->getNumberOfCellNodes();
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  IndexType connec[mint::MAX_CELL_NODES];
  getCellConnec(final_pos, nodes_per_cell, connec);

  mesh->insertCell(connec, pos);
  EXPECT_EQ(++cur_n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cur_n_cells * nodes_per_cell, mesh->getCellNodesSize());
  setFieldTuple(fv, cur_n_cells, pos, final_pos);
}

void insert_cell_single(UnstructuredMesh<MIXED_SHAPE>* mesh, IndexType pos, IndexType final_pos)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellNodesSize();
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  IndexType connec[mint::MAX_CELL_NODES];
  CellType type = getCellConnec(final_pos, connec);
  cur_connec_size += getCellInfo(type).num_nodes;

  mesh->insertCell(connec, pos, type);
  EXPECT_EQ(++cur_n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cur_connec_size, mesh->getCellNodesSize());
  setFieldTuple(fv, cur_n_cells, pos, final_pos);
}

/// @}

/*!
 * \brief Insert multiple cells into the mesh
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] n_cells the number of cells to insert.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this cell after all
 *  insertions have taken place.
 */
/// @{

void insert_cell_multiple(UnstructuredMesh<SINGLE_SHAPE>* mesh,
                          IndexType n_cells,
                          IndexType pos,
                          IndexType final_pos)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellNodesSize();
  const IndexType nodes_per_cell = mesh->getNumberOfCellNodes();
  IndexType* connectivity = new IndexType[nodes_per_cell * n_cells];
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  for(IndexType i = 0; i < n_cells; ++i)
  {
    IndexType* cur_connec = connectivity + i * nodes_per_cell;
    getCellConnec(final_pos + i, nodes_per_cell, cur_connec);
  }

  mesh->insertCells(connectivity, pos, n_cells);
  cur_n_cells += n_cells;
  cur_connec_size += n_cells * nodes_per_cell;
  EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cur_connec_size, mesh->getCellNodesSize());

  for(IndexType i = 0; i < n_cells; ++i)
  {
    setFieldTuple(fv, cur_n_cells, pos + i, final_pos + i);
  }

  delete[] connectivity;
}

void insert_cell_multiple(UnstructuredMesh<MIXED_SHAPE>* mesh,
                          IndexType n_cells,
                          IndexType pos,
                          IndexType final_pos)
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellNodesSize();
  IndexType* connectivity = new IndexType[getConnectivitySize(n_cells)];
  IndexType* offsets = new IndexType[n_cells + 1];
  CellType* types = new CellType[n_cells];
  FieldVariable<double>* fv = getFieldVar(mesh, CELL_CENTERED);

  offsets[0] = 0;
  for(IndexType i = 0; i < n_cells; ++i)
  {
    IndexType* cur_connec = connectivity + offsets[i];
    types[i] = getCellConnec(final_pos + i, cur_connec);
    offsets[i + 1] = offsets[i] + getCellInfo(types[i]).num_nodes;
  }

  mesh->insertCells(connectivity, pos, n_cells, offsets, types);
  cur_n_cells += n_cells;
  cur_connec_size += getConnectivitySize(n_cells);
  EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());
  EXPECT_EQ(cur_connec_size, mesh->getCellNodesSize());

  for(IndexType i = 0; i < n_cells; ++i)
  {
    setFieldTuple(fv, cur_n_cells, pos + i, final_pos + i);
  }

  delete[] connectivity;
  delete[] offsets;
  delete[] types;
}

/// @}

/*!
 * \brief Insert multiple cells into the mesh
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] n_cells the number of cells to insert.
 */
template <Topology TOPO>
void insert_cells(UnstructuredMesh<TOPO>* mesh, IndexType n_cells)
{
  ASSERT_EQ(n_cells % 6, 0);

  IndexType cur_n_cells = mesh->getNumberOfCells();
  ASSERT_EQ(cur_n_cells, 0);

  const IndexType sixth_n_cells = n_cells / 6;
  const IndexType third_n_cells = n_cells / 3;
  const IndexType half_n_cells = n_cells / 2;

  /* The final positions of the cells to insert. */
  IndexType final_pos[3] = {third_n_cells - 1, third_n_cells, 2 * third_n_cells};

  /* The inital position of the cells to insert. */
  IndexType insert_positions[3] = {0, 1, 2};

  /* Insert one cell at a time */
  for(IndexType round = 0; round < sixth_n_cells; ++round)
  {
    for(IndexType i = 0; i < 3; ++i)
    {
      insert_cell_single(mesh, insert_positions[i], final_pos[i]);
      cur_n_cells++;
      EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());
    }

    /* The middle insertion increases by two, the back by three. */
    insert_positions[1] += 2;
    insert_positions[2] += 3;

    /* The final position of the front decreases by one while the final position
     *  of the middle and back both increase by one. */
    final_pos[0]--;
    final_pos[1]++;
    final_pos[2]++;
  }
  EXPECT_EQ(cur_n_cells, half_n_cells);

  /* Insert multiple cells at once using the struct of arrays layout. */
  insert_cell_multiple(mesh, sixth_n_cells, 0, 0);
  cur_n_cells += sixth_n_cells;
  insert_cell_multiple(mesh, sixth_n_cells, half_n_cells, half_n_cells);
  cur_n_cells += sixth_n_cells;
  insert_cell_multiple(mesh, sixth_n_cells, 5 * sixth_n_cells, 5 * sixth_n_cells);
  cur_n_cells += sixth_n_cells;
  EXPECT_EQ(cur_n_cells, n_cells);
  EXPECT_EQ(cur_n_cells, mesh->getNumberOfCells());

  check_append_cells(mesh);
  check_fields(mesh, CELL_CENTERED);
}

/*!
 * \brief Insert nodes into the mesh without updating the cell connectivity.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] n_nodes the number of nodes to insert.
 */
template <Topology TOPO>
void insert_nodes_no_update(UnstructuredMesh<TOPO>* mesh, IndexType n_nodes)
{
  const IndexType initial_n_nodes = mesh->getNumberOfNodes();
  const IndexType half_way = initial_n_nodes / 2;
  const IndexType n_cells = mesh->getNumberOfCells();
  ASSERT_GT(initial_n_nodes, 0);
  ASSERT_GT(n_cells, 0);

  /* Insert a single node. */
  insert_node_single(mesh, 0, 0, false);

  /* Check that the connectivity remains uneffected. */
  check_append_cells(mesh);

  /* Insert multiple nodes. */
  insert_node_arrays(mesh, n_nodes - 1, half_way, half_way, false);

  /* check that the connectivity remains uneffected. */
  check_append_cells(mesh);

  EXPECT_EQ(initial_n_nodes + n_nodes, mesh->getNumberOfNodes());
  EXPECT_EQ(n_cells, mesh->getNumberOfCells());
}

/*!
 * \brief Insert nodes into the mesh while updating the cell connectivity.
 *
 * \param [in/out] mesh the mesh in question.
 * \param [in] n_nodes the number of nodes to insert.
 */
template <Topology TOPO>
void insert_nodes_update(UnstructuredMesh<TOPO>* mesh, IndexType n_nodes)
{
  const IndexType initial_n_nodes = mesh->getNumberOfNodes();
  const IndexType half_way = initial_n_nodes / 2;
  const IndexType n_cells = mesh->getNumberOfCells();
  ASSERT_GT(initial_n_nodes, 0);
  ASSERT_GT(n_cells, 0);

  /* Insert a single node. */
  insert_node_single(mesh, 0, 0);

  /* Check that the connectivity of each cell is incremented by one. */
  for(IndexType i = 0; i < n_cells; ++i)
  {
    const IndexType* cell = mesh->getCellNodeIDs(i);
    const IndexType num_nodes = mesh->getNumberOfCellNodes(i);
    for(IndexType j = 0; j < num_nodes; ++j)
    {
      const IndexType normal_value = getCellConnecValue(i, j);
      EXPECT_EQ(cell[j], normal_value + 1);
    }
  }

  /* Insert the rest of the nodes at the half way point. */
  insert_node_arrays(mesh, n_nodes - 1, half_way, half_way);

  /* check that the connectivity of each cell is changed appropriately. */
  for(IndexType i = 0; i < n_cells; ++i)
  {
    const IndexType* cell = mesh->getCellNodeIDs(i);
    const IndexType num_nodes = mesh->getNumberOfCellNodes(i);
    for(IndexType j = 0; j < num_nodes; ++j)
    {
      const IndexType normal_value = getCellConnecValue(i, j);
      if(normal_value < half_way - 1)
      {
        EXPECT_EQ(cell[j], normal_value + 1);
      }
      else
      {
        EXPECT_EQ(cell[j], normal_value + n_nodes);
      }
    }
  }

  EXPECT_EQ(initial_n_nodes + n_nodes, mesh->getNumberOfNodes());
  EXPECT_EQ(n_cells, mesh->getNumberOfCells());
}

/*!
 * \brief Append nodes to the mesh testing the resizing/capacity functionality.
 *
 * \param [in/out] mesh the mesh in question.
 */
template <Topology TOPO>
void resize_nodes(UnstructuredMesh<TOPO>* mesh)
{
  ASSERT_FALSE(mesh->isExternal());
  ASSERT_TRUE(mesh->empty());

  const int ndims = mesh->getDimension();
  IndexType n_nodes = mesh->getNumberOfNodes();
  IndexType node_capacity = mesh->getNodeCapacity();
  double resize_ratio = mesh->getNodeResizeRatio();
  const double* x_coords = mesh->getCoordinateArray(X_COORDINATE);
  const double* y_coords = nullptr;
  const double* z_coords = nullptr;
  if(ndims > 1)
  {
    y_coords = mesh->getCoordinateArray(Y_COORDINATE);
  }
  if(ndims > 2)
  {
    z_coords = mesh->getCoordinateArray(Z_COORDINATE);
  }

  /* Fill up the array */
  append_node_arrays(mesh, node_capacity - n_nodes);
  n_nodes = node_capacity;
  ASSERT_EQ(n_nodes, mesh->getNumberOfNodes());
  ASSERT_EQ(node_capacity, mesh->getNodeCapacity());
  EXPECT_EQ(x_coords, mesh->getCoordinateArray(X_COORDINATE));
  if(ndims > 1)
  {
    EXPECT_EQ(y_coords, mesh->getCoordinateArray(Y_COORDINATE));
  }
  if(ndims > 2)
  {
    EXPECT_EQ(z_coords, mesh->getCoordinateArray(Z_COORDINATE));
  }

  /* Append one more, should trigger a resize. */
  append_node_single(mesh, 1);
  n_nodes++;
  node_capacity = static_cast<IndexType>(n_nodes * resize_ratio + 0.5);
  ASSERT_EQ(n_nodes, mesh->getNumberOfNodes());
  ASSERT_EQ(node_capacity, mesh->getNodeCapacity());

  /* Shrink. */
  mesh->shrinkNodes();
  node_capacity = n_nodes;
  ASSERT_EQ(n_nodes, mesh->getNumberOfNodes());
  ASSERT_EQ(node_capacity, mesh->getNodeCapacity());

  /* Set the resize ratio and append 100 more nodes, should trigger a resize. */
  resize_ratio = 1.5;
  mesh->setNodeResizeRatio(resize_ratio);
  append_node_structs(mesh, 100);
  n_nodes += 100;
  node_capacity = static_cast<IndexType>(resize_ratio * n_nodes + 0.5);
  ASSERT_EQ(n_nodes, mesh->getNumberOfNodes());
  ASSERT_EQ(node_capacity, mesh->getNodeCapacity());

  /* Reserve 100 more than the current capacity and fill up the array. */
  node_capacity += 100;
  mesh->reserveNodes(node_capacity);
  x_coords = mesh->getCoordinateArray(X_COORDINATE);
  if(ndims > 1)
  {
    y_coords = mesh->getCoordinateArray(Y_COORDINATE);
  }
  if(ndims > 2)
  {
    z_coords = mesh->getCoordinateArray(Z_COORDINATE);
  }
  append_node_arrays(mesh, node_capacity - n_nodes);
  n_nodes = node_capacity;
  ASSERT_EQ(n_nodes, mesh->getNumberOfNodes());
  ASSERT_EQ(node_capacity, mesh->getNodeCapacity());
  EXPECT_EQ(x_coords, mesh->getCoordinateArray(X_COORDINATE));
  if(ndims > 1)
  {
    EXPECT_EQ(y_coords, mesh->getCoordinateArray(Y_COORDINATE));
  }
  if(ndims > 2)
  {
    EXPECT_EQ(z_coords, mesh->getCoordinateArray(Z_COORDINATE));
  }

  /* Set the resize ratio to 1 and append ten values, should trigger a resize.
   */
  resize_ratio = 1.0;
  mesh->setNodeResizeRatio(resize_ratio);
  append_node_structs(mesh, 10);
  n_nodes += 10;
  node_capacity = n_nodes;
  ASSERT_EQ(n_nodes, mesh->getNumberOfNodes());
  ASSERT_EQ(node_capacity, mesh->getNodeCapacity());

  check_append_nodes(mesh);
}

/*!
 * \brief Append cells to the mesh testing the resizing/capacity functionality.
 *
 * \param [in/out] mesh the mesh in question.
 */
template <Topology TOPO>
void resize_cells(UnstructuredMesh<TOPO>* mesh)
{
  ASSERT_FALSE(mesh->isExternal());
  ASSERT_TRUE(mesh->empty());

  IndexType n_cells = mesh->getNumberOfCells();
  IndexType cell_capacity = mesh->getCellCapacity();
  IndexType connec_capacity = mesh->getCellNodesCapacity();
  double resize_ratio = mesh->getCellResizeRatio();
  const IndexType* connectivity = mesh->getCellNodesArray();
  const IndexType* offsets = mesh->getCellNodesOffsetsArray();
  const CellType* types = mesh->getCellTypesArray();

  /* Fill up the array */
  append_cell_multiple(mesh, cell_capacity - n_cells);
  n_cells = cell_capacity;
  ASSERT_EQ(n_cells, mesh->getNumberOfCells());
  ASSERT_EQ(cell_capacity, mesh->getCellCapacity());
  EXPECT_EQ(offsets, mesh->getCellNodesOffsetsArray());
  EXPECT_EQ(types, mesh->getCellTypesArray());
  if(connec_capacity == mesh->getCellNodesCapacity())
  {
    EXPECT_EQ(connectivity, mesh->getCellNodesArray());
  }
  else
  {
    connec_capacity = mesh->getCellNodesCapacity();
  }

  /* Append one more, should trigger a resize. */
  append_cell_single(mesh, 1);
  n_cells++;
  cell_capacity = axom::utilities::max<IndexType>(cell_capacity * resize_ratio + 0.5, n_cells);
  connec_capacity = mesh->getCellNodesCapacity();
  ASSERT_EQ(n_cells, mesh->getNumberOfCells());
  ASSERT_EQ(cell_capacity, mesh->getCellCapacity());

  /* Shrink. */
  mesh->shrinkCells();
  cell_capacity = n_cells;
  connec_capacity = mesh->getCellNodesCapacity();
  ASSERT_EQ(n_cells, mesh->getNumberOfCells());
  ASSERT_EQ(cell_capacity, mesh->getCellCapacity());

  /* Set the resize ratio and append 100 more cells, should trigger a resize. */
  resize_ratio = 1.5;
  mesh->setCellResizeRatio(resize_ratio);
  append_cell_multiple(mesh, 100);
  n_cells += 100;
  cell_capacity = axom::utilities::max<IndexType>(cell_capacity * resize_ratio + 0.5, n_cells);
  ASSERT_EQ(n_cells, mesh->getNumberOfCells());
  ASSERT_EQ(cell_capacity, mesh->getCellCapacity());

  /* Reserve 100 more than the current capacity and fill up the array. */
  cell_capacity += 100;
  mesh->reserveCells(cell_capacity);
  connec_capacity = mesh->getCellNodesCapacity();
  connectivity = mesh->getCellNodesArray();
  offsets = mesh->getCellNodesOffsetsArray();
  types = mesh->getCellTypesArray();
  append_cell_multiple(mesh, cell_capacity - n_cells);
  n_cells = cell_capacity;
  ASSERT_EQ(n_cells, mesh->getNumberOfCells());
  ASSERT_EQ(cell_capacity, mesh->getCellCapacity());
  EXPECT_EQ(offsets, mesh->getCellNodesOffsetsArray());
  EXPECT_EQ(types, mesh->getCellTypesArray());
  if(connec_capacity == mesh->getCellNodesCapacity())
  {
    EXPECT_EQ(connectivity, mesh->getCellNodesArray());
  }
  else
  {
    connec_capacity = mesh->getCellNodesCapacity();
  }

  /* Set the resize ratio to 1 and append ten cells, should trigger a resize. */
  resize_ratio = 1.0;
  mesh->setCellResizeRatio(resize_ratio);
  append_cell_multiple(mesh, 10);
  n_cells += 10;
  cell_capacity = n_cells;
  ASSERT_EQ(n_cells, mesh->getNumberOfCells());
  ASSERT_EQ(cell_capacity, mesh->getCellCapacity());

  check_append_cells(mesh);
}

/*!
 * \brief Check that meshes restored from sidre or external buffers are
 *  equivalent to the original mesh.
 *
 * \param [in/out] meshes the array of meshes.
 * \param [in] n_meshes the number of meshes.
 * \param [in] nodeOpt the check to apply on the nodes.
 * \param [in] cellOpt the check to apply on the cells.
 */
template <Topology TOPO>
void check_restoration(UnstructuredMesh<TOPO>** meshes,
                       int n_meshes,
                       CheckOption nodeOpt,
                       CheckOption cellOpt)
{
  for(int i = 0; i < n_meshes; ++i)
  {
    if(meshes[i]->isExternal())
    {
      deleteAndDuplicateExternalMesh(meshes[i]);
    }
#ifdef AXOM_MINT_USE_SIDRE
    else if(meshes[i]->isInSidre())
    {
      deleteAndDuplicateSidreMesh(meshes[i]);
    }
#endif
    else
    {
      continue;
    }

    if(nodeOpt == CheckOption::APPEND)
    {
      check_append_nodes(meshes[i]);
      check_fields(meshes[i], NODE_CENTERED);
    }
    else if(nodeOpt == CheckOption::SET)
    {
      check_append_nodes(meshes[i], true);
      check_fields(meshes[i], NODE_CENTERED, true);
    }

    if(cellOpt == CheckOption::APPEND)
    {
      check_append_cells(meshes[i]);
      check_fields(meshes[i], CELL_CENTERED);
    }
    else if(cellOpt == CheckOption::SET)
    {
      check_append_cells(meshes[i], true);
      check_fields(meshes[i], CELL_CENTERED, true);
    }
  }
}

struct FaceInfo
{
  FaceInfo() : key(), nodeCount(0), faceType(UNDEFINED_CELL), opposingCellID(-1) { }
  FaceInfo(std::string theKey, IndexType theNodeCount, CellType theFaceType, IndexType oppCellID)
    : key(theKey)
    , nodeCount(theNodeCount)
    , faceType(theFaceType)
    , opposingCellID(oppCellID)
  { }

  std::string key;
  IndexType nodeCount;
  CellType faceType;
  IndexType opposingCellID;
};

bool operator==(const FaceInfo& a, const FaceInfo& b)
{
  return a.key == b.key && a.nodeCount == b.nodeCount && a.faceType == b.faceType &&
    a.opposingCellID == b.opposingCellID;
}

/*!
 * \brief Check that a cell's faces (nodes and types) match between
 *        a UniformMesh and its derived UnstructuredMesh.
 *
 * \param [in] s_mesh the standard mesh, a UniformMesh
 * \param [in] d_mesh the derived mesh, an UnstructuredMesh<SINGLE_SHAPE>
 * \param [in] s_fcount the number of faces in s_mesh
 * \param [in] d_fcount the number of faces in d_mesh
 * \param [in] c the cellID to check faces
 * \param [out] errMesg an error message indicating the nature of any mismatch
 *
 * \returns match true if cell c's faces match between s_mesh and d_mesh,
 *                false otherwise
 *
 * This routine compares the faces of a cell existing in s_mesh and d_mesh.
 * These meshes are assumed to be equivalent, as produced by
 * create_mesh< UNSTRUCTURED_MESH >().  Specifically, cell c refers to the
 * same cell in both meshes.  However, face ID, face ordering, and face
 * node ordering may vary between the two meshes.
 *
 * For each face of cell c of s_mesh, store in a struct FaceInfo
 * - node count
 * - face type (exercising the cell-to-face relation)
 * - opposing cell ID (exercising the face-to-cell and cell-to-face relation)
 * - key, composed of concatenated, sorted face node IDs (exercising the
 *   face-to-node relation)
 * - Put the FaceInfo into a map hashed by its key.
 *
 * For each face of cell c of d_mesh, create the same FaceInfo.  If its
 * key exists for the s_mesh FaceInfo map, and the FaceInfo structs are
 * equal, it matches: remove the entry from the s_mesh map.  Otherwise,
 * d_mesh has an extra face: store the unmatched d_mesh FaceInfo in its own
 * map.
 *
 * At the end of the function, if there are no extra d_mesh faces and no
 * unmatched s_mesh faces, the faces of cell c are equivalent between s_mesh
 * and d_mesh.
 */
bool check_cellFaceNodeType(UniformMesh* s_mesh,
                            UnstructuredMesh<SINGLE_SHAPE>* d_mesh,
                            IndexType s_fcount,
                            IndexType d_fcount,
                            IndexType c,
                            std::string& errMesg)
{
  using FaceInfoMap = std::unordered_map<std::string, FaceInfo>;

  FaceInfoMap s_info, d_info;
  IndexType faceIDs[MAX_CELL_FACES];
  IndexType faceNodeIDs[MAX_FACE_NODES];

  // Create the map for the standard mesh
  s_mesh->getCellFaceIDs(c, faceIDs);
  for(IndexType f = 0; f < s_fcount; ++f)
  {
    IndexType nodeCount = s_mesh->getFaceNodeIDs(faceIDs[f], faceNodeIDs);
    CellType faceType = s_mesh->getFaceType(faceIDs[f]);
    IndexType faceCell1, faceCell2;
    s_mesh->getFaceCellIDs(faceIDs[f], faceCell1, faceCell2);
    IndexType oppCell = (c == faceCell2 ? faceCell1 : faceCell2);
    std::string key = make_face_key(nodeCount, faceNodeIDs, '.');
    s_info[key] = FaceInfo(key, nodeCount, faceType, oppCell);
  }

  // Create the map for the derived mesh
  d_mesh->getCellFaceIDs(c, faceIDs);
  for(IndexType f = 0; f < d_fcount; ++f)
  {
    IndexType nodeCount = d_mesh->getFaceNodeIDs(faceIDs[f], faceNodeIDs);
    CellType faceType = d_mesh->getFaceType(faceIDs[f]);
    IndexType faceCell1, faceCell2;
    d_mesh->getFaceCellIDs(faceIDs[f], faceCell1, faceCell2);
    IndexType oppCell = (c == faceCell2 ? faceCell1 : faceCell2);

    std::string key = make_face_key(nodeCount, faceNodeIDs, '.');
    FaceInfo info(key, nodeCount, faceType, oppCell);

    // If this face exists in s_mesh, and matches,
    bool faceMatches = (s_info.count(key) > 0 && info == s_info[key]);

    // erase it from the s_mesh map.
    if(faceMatches)
    {
      s_info.erase(key);
    }
    // Otherwise, there is a mismatch!!  Store it (as an extra) in d_mesh map.
    else
    {
      d_info[key] = info;
    }
  }

  std::stringstream msg;

  // Compose an error message if necessary, adding extra d_mesh faces.
  if(d_info.size() > 0)
  {
    msg << d_info.size() << " extra UnstructuredMesh faces:\n";
  }
  for(auto d_face : d_info)
  {
    FaceInfo& fi = d_face.second;
    msg << "  type " << getCellInfo(fi.faceType).name << " node count " << fi.nodeCount
        << " nodes (" << fi.key << ") opposite cell ID " << fi.opposingCellID << std::endl;
  }

  // Add any s_mesh faces to the error message.
  if(s_info.size() > 0)
  {
    msg << d_info.size() << " missing UniformMesh faces:\n";
  }
  for(auto s_face : s_info)
  {
    FaceInfo& fi = s_face.second;
    msg << "  type " << getCellInfo(fi.faceType).name << " node count " << fi.nodeCount
        << " nodes (" << fi.key << ") opposite cell ID " << fi.opposingCellID << std::endl;
  }
  errMesg = msg.str();

  // Success is defined as no un-matched and no extra faces for cell c.
  return d_info.size() < 1 && s_info.size() < 1;
}

/*!
 * \brief Check for consistency between faces of an UnstructuredMesh
 * and faces of the UniformMesh it is based on.
 *
 * \param [in] d_mesh the UnstructuredMesh
 * \param [in] s_mesh the UniformMesh
 */
void check_faces(UnstructuredMesh<SINGLE_SHAPE>* d_mesh, UniformMesh* s_mesh)
{
  // check total number of faces
  EXPECT_EQ(d_mesh->getNumberOfFaces(), s_mesh->getNumberOfFaces());

  IndexType cellcount = d_mesh->getNumberOfCells();

  for(IndexType c = 0; c < cellcount; ++c)
  {
    // c is the cell ID in d_mesh.  Assume the same cell ID
    // in s_mesh due to create_mesh()'s implementation.

    // check face count
    IndexType d_facecount = d_mesh->getNumberOfCellFaces(c);
    IndexType s_facecount = s_mesh->getNumberOfCellFaces(c);
    EXPECT_EQ(d_facecount, s_facecount);

    // check face nodes and type, and neighbor across face
    std::string errMesg;
    bool facesMatch = check_cellFaceNodeType(s_mesh, d_mesh, s_facecount, d_facecount, c, errMesg);
    EXPECT_TRUE(facesMatch) << errMesg;
  }
}

} /* end namespace internal */

/*******************************************************************************
 *                          Append tests                                       *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, appendNodesSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, appendNodesMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;
  constexpr IndexType CONNEC_SIZE = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, appendCellsSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::append_cells(meshes[i], N_CELLS);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::APPEND);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, appendCellsMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize(N_CELLS);

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::append_cells(meshes[i], N_CELLS);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::APPEND);

  internal::deleteMeshes(meshes, N_MESHES);
}

/*******************************************************************************
 *                             Set tests                                       *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, setNodesSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::set_nodes(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::SET, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, setNodesMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;
  constexpr IndexType CONNEC_SIZE = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::set_nodes(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::SET, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, setCellsSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::append_cells(meshes[i], N_CELLS);
    internal::set_cells(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::SET);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, setCellsMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize(N_CELLS);

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::append_cells(meshes[i], N_CELLS);
    internal::set_cells(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::SET);

  internal::deleteMeshes(meshes, N_MESHES);
}

/*******************************************************************************
 *                             Insert tests                                    *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertNodesSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 900;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::insert_nodes(meshes[i], N_NODES);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertNodesMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 900;
  constexpr IndexType N_CELLS = 0;
  constexpr IndexType CONNEC_SIZE = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::insert_nodes(meshes[i], N_NODES);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertCellsSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::insert_cells(meshes[i], N_CELLS);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::APPEND);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertCellsMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize(N_CELLS);

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES);
    internal::insert_cells(meshes[i], N_CELLS);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::APPEND);

  internal::deleteMeshes(meshes, N_MESHES);
}

/*******************************************************************************
 *                             Insert and update tests                         *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertNodesNoUpdateSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_NODES_INSERT = 300;
  constexpr IndexType N_CELLS = 900;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES - N_NODES_INSERT);
    internal::append_cells(meshes[i], N_CELLS);
    internal::insert_nodes_no_update(meshes[i], N_NODES_INSERT);
  }

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertNodesNoUpdateMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_NODES_INSERT = 300;
  constexpr IndexType N_CELLS = 900;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize(N_CELLS);

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES - N_NODES_INSERT);
    internal::append_cells(meshes[i], N_CELLS);
    internal::insert_nodes_no_update(meshes[i], N_NODES_INSERT);
  }

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertNodesUpdateSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_NODES_INSERT = 300;
  constexpr IndexType N_CELLS = 900;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES - N_NODES_INSERT);
    internal::append_cells(meshes[i], N_CELLS);
    internal::insert_nodes_update(meshes[i], N_NODES_INSERT);
  }

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, insertNodesUpdateMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_NODES_INSERT = 300;
  constexpr IndexType N_CELLS = 900;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize(N_CELLS);

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshes(meshes, N_NODES, N_CELLS, CONNEC_SIZE, true, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::append_nodes(meshes[i], N_NODES - N_NODES_INSERT);
    internal::append_cells(meshes[i], N_CELLS);
    internal::insert_nodes_update(meshes[i], N_NODES_INSERT);
  }

  internal::deleteMeshes(meshes, N_MESHES);
}

/*******************************************************************************
 *                             Resize tests                                    *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, resizeNodesSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 2;
#else
  constexpr int STRIDE = 1;
#endif

  constexpr IndexType N_NODES = 100;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshesForResize(meshes, N_NODES, N_CELLS, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::resize_nodes(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, resizeNodesMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 2;
#else
  constexpr int STRIDE = 1;
#endif

  constexpr IndexType N_NODES = 100;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshesForResize(meshes, N_NODES, N_CELLS, true, false);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::resize_nodes(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::APPEND, CheckOption::NONE);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, resizeCellsSingle)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 2;
#else
  constexpr int STRIDE = 1;
#endif

  constexpr IndexType N_NODES = 100;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  internal::createMeshesForResize(meshes, N_NODES, N_CELLS, false, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::resize_cells(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::NONE, CheckOption::APPEND);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, resizeCellsMixed)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 2;
#else
  constexpr int STRIDE = 1;
#endif

  constexpr IndexType N_NODES = 0;
  constexpr IndexType N_CELLS = 100;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  internal::createMeshesForResize(meshes, N_NODES, N_CELLS, false, true);

  for(int i = 0; i < N_MESHES; ++i)
  {
    internal::resize_cells(meshes[i]);
  }

  internal::check_restoration(meshes, N_MESHES, CheckOption::NONE, CheckOption::APPEND);

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, resize_mesh_single_topology)
{
  constexpr int NDIMS = 2;
  constexpr IndexType N_NODES = 6;
  constexpr IndexType N_CELLS = 2;
  constexpr double MAGIC_DOUBLE = 42.0;
  constexpr IndexType MAGIC_INT = 42;

  using UnstructuredMeshType = mint::UnstructuredMesh<SINGLE_SHAPE>;

  UnstructuredMeshType mesh(NDIMS, mint::QUAD);
  EXPECT_EQ(mesh.getNumberOfNodes(), 0);
  EXPECT_EQ(mesh.getNumberOfCells(), 0);
  EXPECT_FALSE(mesh.hasMixedCellTypes());
  EXPECT_TRUE(mesh.getCellNodesOffsetsArray() == nullptr);
  EXPECT_TRUE(mesh.getCellTypesArray() == nullptr);

  mesh.resize(N_NODES, N_CELLS);
  EXPECT_EQ(mesh.getNumberOfNodes(), N_NODES);
  EXPECT_EQ(mesh.getNumberOfCells(), N_CELLS);

  double* x = mesh.getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh.getCoordinateArray(mint::Y_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);

  // touch memory associated with nodes to ensure it is allocated properly
  for(IndexType inode = 0; inode < N_NODES; ++inode)
  {
    x[inode] = MAGIC_DOUBLE;
    y[inode] = MAGIC_DOUBLE;
  }

  // touch memory associated with the cell connectivity
  IndexType stride = mesh.getNumberOfCellNodes();
  EXPECT_EQ(stride, 4);
  IndexType* conn = mesh.getCellNodesArray();
  for(IndexType icell = 0; icell < N_CELLS; ++icell)
  {
    conn[icell * 4] = MAGIC_INT;
    conn[icell * 4 + 1] = MAGIC_INT;
    conn[icell * 4 + 2] = MAGIC_INT;
    conn[icell * 4 + 3] = MAGIC_INT;
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, resize_mesh_mixed_topology)
{
  constexpr int NDIMS = 2;
  constexpr IndexType N_NODES = 6;
  constexpr IndexType N_CELLS = 3;
  constexpr double MAGIC_DOUBLE = 42.0;

  using UnstructuredMeshType = mint::UnstructuredMesh<MIXED_SHAPE>;

  UnstructuredMeshType mesh(NDIMS);
  EXPECT_EQ(mesh.getNumberOfNodes(), 0);
  EXPECT_EQ(mesh.getNumberOfCells(), 0);
  EXPECT_TRUE(mesh.hasMixedCellTypes());

  mesh.resize(N_NODES, N_CELLS);
  EXPECT_EQ(mesh.getNumberOfNodes(), N_NODES);
  EXPECT_EQ(mesh.getNumberOfCells(), N_CELLS);
  EXPECT_TRUE(mesh.getCellNodesOffsetsArray() != nullptr);
  EXPECT_TRUE(mesh.getCellTypesArray() != nullptr);

  double* x = mesh.getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh.getCoordinateArray(mint::Y_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);

  // touch memory associated with nodes to ensure it is allocated properly
  for(IndexType inode = 0; inode < N_NODES; ++inode)
  {
    x[inode] = MAGIC_DOUBLE;
    y[inode] = MAGIC_DOUBLE;
  }

  // touch memory associated with the cell connectivity
  IndexType* conn = mesh.getCellNodesArray();
  IndexType* offsets = mesh.getCellNodesOffsetsArray();
  CellType* types = mesh.getCellTypesArray();

  for(IndexType icell = 0; icell < N_CELLS; ++icell)
  {
    if(icell < 2)
    {
      types[icell] = mint::TRIANGLE;
      offsets[icell] = (icell == 0) ? 0 : offsets[icell];
      offsets[icell + 1] = offsets[icell] + 3;

      IndexType* triangle = &conn[offsets[icell]];
      triangle[0] = icell;
      triangle[1] = icell + 1;
      triangle[2] = icell + 2;
    }
    else
    {
      types[icell] = mint::QUAD;
      offsets[icell + 1] = offsets[icell] + 4;

      IndexType* quad = &conn[offsets[icell]];
      quad[0] = icell;
      quad[1] = icell + 1;
      quad[2] = icell + 2;
      quad[3] = icell + 3;
    }
  }
}

/*******************************************************************************
 *                             External death tests                            *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh_DeathTest, invalidExternalOpsSingle)
{
  constexpr IndexType N_NODES = 100;
  constexpr IndexType N_CELLS = 100;
  constexpr int N_MESHES = 3;
  UnstructuredMesh<SINGLE_SHAPE>* meshes[N_MESHES];

  int cur_mesh = 0;
  for(int dim = 1; dim <= 3; ++dim)
  {
    meshes[cur_mesh] = internal::createExternalSingle(dim, QUAD, N_NODES, N_CELLS);
    internal::createExternalField(meshes[cur_mesh], NODE_CENTERED, N_NODES);
    internal::createExternalField(meshes[cur_mesh], CELL_CENTERED, N_CELLS);

    internal::append_node_arrays(meshes[cur_mesh], N_NODES);
    internal::append_cell_multiple(meshes[cur_mesh], N_CELLS);
    cur_mesh++;
  }

  for(int i = 0; i < N_MESHES; ++i)
  {
    /* These are all valid operations. */
    meshes[i]->shrinkNodes();
    meshes[i]->shrinkCells();
    meshes[i]->reserveNodes(N_NODES);
    meshes[i]->reserveCells(N_CELLS);

    /* These are invalid */
    EXPECT_DEATH_IF_SUPPORTED(internal::append_node_single(meshes[i], 1), IGNORE_OUTPUT);
    EXPECT_DEATH_IF_SUPPORTED(internal::append_cell_single(meshes[i], 1), IGNORE_OUTPUT);
    EXPECT_DEATH_IF_SUPPORTED(meshes[i]->reserveNodes(N_NODES * 2), IGNORE_OUTPUT);
    EXPECT_DEATH_IF_SUPPORTED(meshes[i]->reserveCells(N_NODES * 2), IGNORE_OUTPUT);
  }

  internal::deleteMeshes(meshes, N_MESHES);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh_DeathTest, invalidExternalOpsMixed)
{
  constexpr IndexType N_NODES = 100;
  constexpr IndexType N_CELLS = 100;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize(N_CELLS);
  constexpr int N_MESHES = 3;
  UnstructuredMesh<MIXED_SHAPE>* meshes[N_MESHES];

  int cur_mesh = 0;
  for(int dim = 1; dim <= 3; ++dim)
  {
    meshes[cur_mesh] = internal::createExternalMixed(dim, N_NODES, N_CELLS, CONNEC_SIZE);
    internal::createExternalField(meshes[cur_mesh], NODE_CENTERED, N_NODES);
    internal::createExternalField(meshes[cur_mesh], CELL_CENTERED, N_CELLS);

    internal::append_node_arrays(meshes[cur_mesh], N_NODES);
    internal::append_cell_multiple(meshes[cur_mesh], N_CELLS);
    cur_mesh++;
  }

  for(int i = 0; i < N_MESHES; ++i)
  {
    /* These are all valid operations. */
    meshes[i]->shrinkNodes();
    meshes[i]->shrinkCells();
    meshes[i]->reserveNodes(N_NODES);
    meshes[i]->reserveCells(N_CELLS);

    /* These are invalid */
    EXPECT_DEATH_IF_SUPPORTED(internal::append_node_single(meshes[i], 1), IGNORE_OUTPUT);
    EXPECT_DEATH_IF_SUPPORTED(internal::append_cell_single(meshes[i], 1), IGNORE_OUTPUT);
    EXPECT_DEATH_IF_SUPPORTED(meshes[i]->reserveNodes(N_NODES * 2), IGNORE_OUTPUT);
    EXPECT_DEATH_IF_SUPPORTED(meshes[i]->reserveCells(N_CELLS * 2), IGNORE_OUTPUT);
  }

  internal::deleteMeshes(meshes, N_MESHES);
}

/*******************************************************************************
 *                             Face relation tests                             *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_mesh_unstructured_mesh, check_face_connectivity)
{
  constexpr IndexType N_NODES = 5;

  constexpr int SIZE_STEPS = 2;
  constexpr IndexType SIZE_FACTOR[] = {1, 3};

  constexpr int DIM_STEPS = 2;
  constexpr int DIMENSIONS[] = {2, 3};

  for(int isize = 0; isize < SIZE_STEPS; ++isize)
  {
    const IndexType resolution = SIZE_FACTOR[isize] * N_NODES;
    for(int idim = 0; idim < DIM_STEPS; ++idim)
    {
      const int dimension = DIMENSIONS[idim];

      UniformMesh* source_mesh = nullptr;
      UnstructuredMesh<SINGLE_SHAPE>* test_mesh = nullptr;

      internal::createMeshesForFace(test_mesh, source_mesh, resolution, dimension);

      internal::check_faces(test_mesh, source_mesh);

      SLIC_INFO("Tested " << dimension << "D mesh with base resolution " << resolution);

      delete source_mesh;
      source_mesh = nullptr;
      delete test_mesh;
      test_mesh = nullptr;
    }
  }
}

} /* end namespace mint */
} /* end namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
