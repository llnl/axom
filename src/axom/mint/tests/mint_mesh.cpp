// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/config.hpp"  // for compile-time type definitions

// Mint includes
#include "axom/mint/mesh/CurvilinearMesh.hpp"   // for CurvilinearMesh
#include "axom/mint/mesh/Mesh.hpp"              // for Mesh base class
#include "axom/mint/mesh/ParticleMesh.hpp"      // for ParticleMesh definition
#include "axom/mint/mesh/RectilinearMesh.hpp"   // for RectilinearMesh
#include "axom/mint/mesh/UniformMesh.hpp"       // for UniformMesh
#include "axom/mint/mesh/UnstructuredMesh.hpp"  // for UnstructuredMesh

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "axom/slic.hpp"
#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>    // for cosh, etc.
#include <cstring>  // for memcpy

namespace axom
{
namespace mint
{
constexpr double E = 2.71828182845904523536;
const char* IGNORE_OUTPUT = ".*";

namespace internal
{
/*!
 * \brief Append a structured grid of nodes to the mesh.
 *
 * \param [in/out] mesh the mesh to append the nodes to.
 * \param [in] x_extent the number of nodes in the x direction.
 * \param [in] y_extent the number of nodes in the y direction.
 * \param [in] spacing the spacing between the nodes.
 */
template <Topology TOPO>
void append_nodes(UnstructuredMesh<TOPO>* mesh, IndexType x_extent, IndexType y_extent, double spacing)
{
  for(IndexType j = 0; j < y_extent; ++j)
  {
    for(IndexType i = 0; i < x_extent; ++i)
    {
      mesh->appendNode(i * spacing, j * spacing);
    }
  }
}

/*!
 * \brief Check that the nodes were appended correctly.
 *
 * \param [in] mesh the mesh to check.
 * \param [in] x_extent the number of nodes in the x direction.
 * \param [in] y_extent the number of nodes in the y direction.
 * \param [in] spacing the spacing between the nodes.
 */
void check_append_nodes(const Mesh* mesh, IndexType x_extent, IndexType y_extent, double spacing)
{
  ASSERT_EQ(x_extent * y_extent, mesh->getNumberOfNodes());
  ASSERT_EQ(mesh->getDimension(), 2);

  IndexType node_ID = 0;
  double node[2];
  for(IndexType j = 0; j < y_extent; ++j)
  {
    for(IndexType i = 0; i < x_extent; ++i)
    {
      mesh->getNode(node_ID++, node);
      EXPECT_EQ(node[0], i * spacing);
      EXPECT_EQ(node[1], j * spacing);
    }
  }
}

/*!
 * \brief Append a structured grid of cells to the mesh.
 *
 * \param [in/out] mesh the mesh to append the cells to.
 * \param [in] x_extent the number of nodes in the x direction.
 * \param [in] y_extent the number of nodes in the y direction.
 */
void append_cells(UnstructuredMesh<SINGLE_SHAPE>* mesh, IndexType x_extent, IndexType y_extent)
{
  for(IndexType j = 0; j < y_extent - 1; ++j)
  {
    for(IndexType i = 0; i < x_extent - 1; ++i)
    {
      const IndexType bottom_left = j * x_extent + i;
      const IndexType bottom_right = bottom_left + 1;
      const IndexType top_right = bottom_right + x_extent;
      const IndexType top_left = bottom_left + x_extent;
      const IndexType cell[4] = {bottom_left, bottom_right, top_right, top_left};
      mesh->appendCell(cell);
    }
  }
}

/*!
 * \brief Check that the cells were appended correctly.
 *
 * \param [in] mesh the mesh to check.
 * \param [in] x_extent the number of nodes in the x direction.
 * \param [in] y_extent the number of nodes in the y direction.
 */
void check_append_cells_single(const Mesh* mesh, IndexType x_extent, IndexType y_extent)
{
  ASSERT_FALSE(mesh->hasMixedCellTypes());
  ASSERT_EQ(mesh->getNumberOfCellNodes(), 4);
  ASSERT_EQ(mesh->getCellType(), QUAD);

  IndexType cell_ID = 0;
  IndexType cell[4];
  for(IndexType j = 0; j < y_extent - 1; ++j)
  {
    for(IndexType i = 0; i < x_extent - 1; ++i)
    {
      mesh->getCellNodeIDs(cell_ID++, cell);

      const IndexType bottom_left = j * x_extent + i;
      const IndexType bottom_right = bottom_left + 1;
      const IndexType top_right = bottom_right + x_extent;
      const IndexType top_left = bottom_left + x_extent;

      EXPECT_EQ(cell[0], bottom_left);
      EXPECT_EQ(cell[1], bottom_right);
      EXPECT_EQ(cell[2], top_right);
      EXPECT_EQ(cell[3], top_left);
    }
  }
}

/*!
 * \brief Return true iff the cell (i, j) should be a quad.
 */
inline bool appendQuad(axom::IndexType i, axom::IndexType j)
{
  return (i % 2 == 0 && j % 2 == 0) || (i % 2 == 1 && j % 2 == 1);
}

/*!
 * \brief Append a structured grid of cells to the mesh. the cells alternate
 *  between quads and triangles.
 *
 * \param [in/out] mesh the mesh to append the cells to.
 * \param [in] x_extent the number of nodes in the x direction.
 * \param [in] y_extent the number of nodes in the y direction.
 */
void append_cells(UnstructuredMesh<MIXED_SHAPE>* mesh, IndexType x_extent, IndexType y_extent)
{
  for(IndexType j = 0; j < y_extent - 1; ++j)
  {
    for(IndexType i = 0; i < x_extent - 1; ++i)
    {
      const IndexType bottom_left = j * x_extent + i;
      const IndexType bottom_right = bottom_left + 1;
      const IndexType top_right = bottom_right + x_extent;
      const IndexType top_left = bottom_left + x_extent;

      if(appendQuad(i, j))
      {
        /* Append a quad. */
        const IndexType quad[4] = {bottom_left, bottom_right, top_right, top_left};
        mesh->appendCell(quad, QUAD);
      }
      else
      {
        /* Append two triangles. */
        const IndexType tri0[3] = {bottom_left, bottom_right, top_right};
        mesh->appendCell(tri0, TRIANGLE);

        const IndexType tri1[3] = {top_right, top_left, bottom_left};
        mesh->appendCell(tri1, TRIANGLE);
      }
    }
  }
}

/*!
 * \brief Check that the cells were appended correctly.
 *
 * \param [in] mesh the mesh to check.
 * \param [in] x_extent the number of nodes in the x direction.
 * \param [in] y_extent the number of nodes in the y direction.
 */
void check_append_cells_mixed(const Mesh* mesh, IndexType x_extent, IndexType y_extent)
{
  ASSERT_TRUE(mesh->hasMixedCellTypes());

  IndexType cell_ID = 0;
  IndexType cell[MAX_CELL_NODES];
  for(IndexType j = 0; j < y_extent - 1; ++j)
  {
    for(IndexType i = 0; i < x_extent - 1; ++i)
    {
      CellType type = mesh->getCellType(cell_ID);
      IndexType n_nodes = mesh->getCellNodeIDs(cell_ID++, cell);

      const IndexType bottom_left = j * x_extent + i;
      const IndexType bottom_right = bottom_left + 1;
      const IndexType top_right = bottom_right + x_extent;
      const IndexType top_left = bottom_left + x_extent;

      if(appendQuad(i, j))
      {
        EXPECT_EQ(n_nodes, 4);
        EXPECT_EQ(type, QUAD);
        EXPECT_EQ(cell[0], bottom_left);
        EXPECT_EQ(cell[1], bottom_right);
        EXPECT_EQ(cell[2], top_right);
        EXPECT_EQ(cell[3], top_left);
      }
      else
      {
        EXPECT_EQ(n_nodes, 3);
        EXPECT_EQ(type, TRIANGLE);
        EXPECT_EQ(cell[0], bottom_left);
        EXPECT_EQ(cell[1], bottom_right);
        EXPECT_EQ(cell[2], top_right);

        type = mesh->getCellType(cell_ID);
        n_nodes = mesh->getCellNodeIDs(cell_ID++, cell);

        EXPECT_EQ(n_nodes, 3);
        EXPECT_EQ(type, TRIANGLE);
        EXPECT_EQ(cell[0], top_right);
        EXPECT_EQ(cell[1], top_left);
        EXPECT_EQ(cell[2], bottom_left);
      }
    }
  }
}

/*!
 * \brief Set the node fields.
 *
 * \param [in] n_nodes the number of nodes.
 * \param [out] vx the first field.
 * \param [out] vy the second field.
 */
void set_node_fields(IndexType n_nodes, double* vx, double* vy)
{
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    vx[i] = std::cos(E * i);
    vy[i] = std::sin(E * E * i);
  }
}

/*!
 * \brief Check that the node fields were set correctly.
 *
 * \param [in] n_nodes the number of nodes.
 * \param [in] vx the first field.
 * \param [in] vy the second field.
 */
void check_node_fields(IndexType n_nodes, const double* vx, const double* vy)
{
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    EXPECT_DOUBLE_EQ(vx[i], std::cos(E * i));
    EXPECT_DOUBLE_EQ(vy[i], std::sin(E * E * i));
  }
}

/*!
 * \brief Set the cell fields.
 *
 * \param [in] n_nodes the number of nodes.
 * \param [out] p the first field.
 */
void set_cell_fields(IndexType n_cells, double* p)
{
  for(IndexType i = 0; i < n_cells; ++i)
  {
    p[i] = std::cosh(E * i);
  }
}

/*!
 * \brief Check that the cell fields were set correctly.
 *
 * \param [in] n_nodes the number of nodes.
 * \param [out] p the first field.
 */
void check_cell_fields(IndexType n_cells, const double* p)
{
  for(IndexType i = 0; i < n_cells; ++i)
  {
    EXPECT_DOUBLE_EQ(p[i], std::cosh(E * i));
  }
}

} /* namespace internal */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_DeathTest, enforce_unique_field_names)
{
  mint::CurvilinearMesh m(5, 5);
  m.createField<double>("foo", mint::NODE_CENTERED);

  EXPECT_DEATH_IF_SUPPORTED(m.createField<double>("foo", mint::CELL_CENTERED), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE

TEST(mint_mesh_DeathTest, get_mesh_null_group)
{
  EXPECT_DEATH_IF_SUPPORTED(mint::getMesh(nullptr), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(mint_mesh, get_curvilinear_mesh_from_sidre)
{
  constexpr int DIMENSION = 3;
  constexpr int BLOCKID = 9;
  constexpr int PARTID = 10;

  /* STEP 0: get empty Sidre group where to store the mesh. */
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* STEP 1: populate group with a CurvilinearMesh */
  CurvilinearMesh* cm = new CurvilinearMesh(root, 20, 21, 22);
  cm->setBlockId(BLOCKID);
  cm->setPartitionId(PARTID);
  double* foo = cm->createField<double>("foo", NODE_CENTERED);
  double* x = cm->getCoordinateArray(X_COORDINATE);
  double* y = cm->getCoordinateArray(Y_COORDINATE);
  double* z = cm->getCoordinateArray(Z_COORDINATE);
  const IndexType numNodes = cm->getNumberOfNodes();
  const IndexType numCells = cm->getNumberOfCells();

  delete cm;

  /* STEP 2: get a mesh object from the group */
  Mesh* m = mint::getMesh(root);

  /* STEP 3: check mesh */
  EXPECT_EQ(m->getMeshType(), STRUCTURED_CURVILINEAR_MESH);
  EXPECT_EQ(m->getBlockId(), BLOCKID);
  EXPECT_EQ(m->getPartitionId(), PARTID);
  EXPECT_EQ(m->getDimension(), DIMENSION);
  EXPECT_EQ(m->getNumberOfNodes(), numNodes);
  EXPECT_EQ(m->getNumberOfCells(), numCells);
  EXPECT_TRUE(m->hasField("foo", NODE_CENTERED));
  EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
  EXPECT_EQ(m->getCoordinateArray(Y_COORDINATE), y);
  EXPECT_EQ(m->getCoordinateArray(Z_COORDINATE), z);

  double* foo_test = m->getFieldPtr<double>("foo", NODE_CENTERED);
  EXPECT_EQ(foo, foo_test);

  const CurvilinearMesh* M = dynamic_cast<const CurvilinearMesh*>(m);
  EXPECT_TRUE(M != nullptr);

  delete m;
}

//------------------------------------------------------------------------------
TEST(mint_mesh, get_rectilinear_mesh_from_sidre)
{
  constexpr int DIMENSION = 3;
  constexpr int BLOCKID = 9;
  constexpr int PARTID = 10;

  /* STEP 0: get empty Sidre group where to store the mesh. */
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* STEP 1: populate group with a RectilinearMesh */
  RectilinearMesh* rm = new RectilinearMesh(root, 20, 21, 22);
  rm->setBlockId(BLOCKID);
  rm->setPartitionId(PARTID);
  double* foo = rm->createField<double>("foo", NODE_CENTERED);
  double* x = rm->getCoordinateArray(X_COORDINATE);
  double* y = rm->getCoordinateArray(Y_COORDINATE);
  double* z = rm->getCoordinateArray(Z_COORDINATE);
  const IndexType numNodes = rm->getNumberOfNodes();
  const IndexType numCells = rm->getNumberOfCells();

  delete rm;

  /* STEP 2: get a mesh object from the group */
  Mesh* m = mint::getMesh(root);

  /* STEP 3: check mesh */
  EXPECT_EQ(m->getMeshType(), STRUCTURED_RECTILINEAR_MESH);
  EXPECT_EQ(m->getBlockId(), BLOCKID);
  EXPECT_EQ(m->getPartitionId(), PARTID);
  EXPECT_EQ(m->getDimension(), DIMENSION);
  EXPECT_EQ(m->getNumberOfNodes(), numNodes);
  EXPECT_EQ(m->getNumberOfCells(), numCells);
  EXPECT_TRUE(m->hasField("foo", NODE_CENTERED));
  EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
  EXPECT_EQ(m->getCoordinateArray(Y_COORDINATE), y);
  EXPECT_EQ(m->getCoordinateArray(Z_COORDINATE), z);

  double* foo_test = m->getFieldPtr<double>("foo", NODE_CENTERED);
  EXPECT_EQ(foo, foo_test);

  const RectilinearMesh* M = dynamic_cast<const RectilinearMesh*>(m);
  EXPECT_TRUE(M != nullptr);

  delete m;
}

//------------------------------------------------------------------------------
TEST(mint_mesh, get_uniform_mesh_from_sidre)
{
  const double LO[] = {-2.0, -2.0, -2.0};
  const double HI[] = {2.0, 2.0, 2.0};
  constexpr int DIMENSION = 3;
  constexpr int BLOCKID = 9;
  constexpr int PARTID = 10;

  /* STEP 0: get empty Sidre group where to store the mesh. */
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* STEP 1: populate group with a UniformMesh */
  UniformMesh* um = new UniformMesh(root, LO, HI, 20, 21, 22);
  um->setBlockId(BLOCKID);
  um->setPartitionId(PARTID);
  double* foo = um->createField<double>("foo", NODE_CENTERED);

  double SPACING[3];
  memcpy(SPACING, um->getSpacing(), 3 * sizeof(double));

  const IndexType numNodes = um->getNumberOfNodes();
  const IndexType numCells = um->getNumberOfCells();

  delete um;

  /* STEP 2: get a mesh object from the group */
  Mesh* m = mint::getMesh(root);

  /* STEP 3: check mesh */
  EXPECT_EQ(m->getMeshType(), STRUCTURED_UNIFORM_MESH);
  EXPECT_EQ(m->getBlockId(), BLOCKID);
  EXPECT_EQ(m->getPartitionId(), PARTID);
  EXPECT_EQ(m->getDimension(), DIMENSION);
  EXPECT_EQ(m->getNumberOfNodes(), numNodes);
  EXPECT_EQ(m->getNumberOfCells(), numCells);
  EXPECT_TRUE(m->hasField("foo", NODE_CENTERED));

  double* foo_test = m->getFieldPtr<double>("foo", NODE_CENTERED);
  EXPECT_EQ(foo, foo_test);

  /* STEP 4: down-cast and test the UniformMesh object. */
  const UniformMesh* M = dynamic_cast<const UniformMesh*>(m);
  EXPECT_TRUE(M != nullptr);

  const double* h = M->getSpacing();
  const double* lo = M->getOrigin();
  for(int idim = 0; idim < DIMENSION; ++idim)
  {
    EXPECT_DOUBLE_EQ(SPACING[idim], h[idim]);
    EXPECT_DOUBLE_EQ(LO[idim], lo[idim]);
  }

  delete m;
}

//------------------------------------------------------------------------------
TEST(mint_mesh, get_particle_mesh_from_sidre)
{
  constexpr int DIMENSION = 3;
  constexpr IndexType NUM_PARTICLES = 10;
  constexpr double MAGIC_NUMBER = 42.0;
  constexpr int BLOCKID = 9;
  constexpr int PARTID = 10;

  /* STEP 0: get empty Sidre group where to store the mesh. */
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* STEP 1: populate the group with a particle mesh (positions + fields). */
  ParticleMesh* particles = new ParticleMesh(DIMENSION, NUM_PARTICLES, root);
  particles->setBlockId(BLOCKID);
  particles->setPartitionId(PARTID);
  double* x = particles->getCoordinateArray(X_COORDINATE);
  double* y = particles->getCoordinateArray(Y_COORDINATE);
  double* z = particles->getCoordinateArray(Z_COORDINATE);

  double* phi = particles->createField<double>("phi", NODE_CENTERED);
  IndexType* id = particles->createField<IndexType>("id", NODE_CENTERED);

  for(IndexType ipart = 0; ipart < NUM_PARTICLES; ++ipart)
  {
    const double val = ipart + 1;
    x[ipart] = y[ipart] = z[ipart] = val * val;
    phi[ipart] = MAGIC_NUMBER;
    id[ipart] = ipart;
  }

  delete particles;

  /* STEP 2: get a mesh object from the group. */
  Mesh* m = mint::getMesh(root);

  /* STEP 3: test the object. */
  EXPECT_EQ(m->getMeshType(), PARTICLE_MESH);
  EXPECT_EQ(m->getBlockId(), BLOCKID);
  EXPECT_EQ(m->getPartitionId(), PARTID);
  EXPECT_EQ(m->getDimension(), DIMENSION);
  EXPECT_EQ(m->getNumberOfNodes(), NUM_PARTICLES);
  EXPECT_TRUE(m->hasField("phi", NODE_CENTERED));
  EXPECT_TRUE(m->hasField("id", NODE_CENTERED));
  EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
  EXPECT_EQ(m->getCoordinateArray(Y_COORDINATE), y);
  EXPECT_EQ(m->getCoordinateArray(Z_COORDINATE), z);

  double* phi_test = m->getFieldPtr<double>("phi", NODE_CENTERED);
  EXPECT_EQ(phi, phi_test);

  IndexType* id_test = m->getFieldPtr<IndexType>("id", NODE_CENTERED);
  EXPECT_EQ(id, id_test);

  /* STEP 4: de-allocate. */
  delete m;
}

//------------------------------------------------------------------------------
TEST(mint_mesh, get_single_topology_unstructured_from_sidre)
{
  constexpr int DIMENSION = 2;
  constexpr IndexType X_EXTENT = 11;
  constexpr IndexType Y_EXTENT = 11;
  constexpr double SPACING = 1.0;
  constexpr CellType CELL_TYPE = QUAD;
  constexpr int BLOCKID = 9;
  constexpr int PARTID = 10;

  // STEP 0: get empty Sidre group where to store the mesh
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* STEP 1: create the UnstructuredMesh */
  UnstructuredMesh<SINGLE_SHAPE>* mesh =
    new UnstructuredMesh<SINGLE_SHAPE>(DIMENSION, CELL_TYPE, root);
  mesh->setBlockId(BLOCKID);
  mesh->setPartitionId(PARTID);

  /* STEP 2: Add the nodes and cells. */
  internal::append_nodes(mesh, X_EXTENT, Y_EXTENT, SPACING);
  internal::append_cells(mesh, X_EXTENT, Y_EXTENT);

  /* STEP 3: Add and fill the fields. */
  double* vx = mesh->createField<double>("vx", NODE_CENTERED);
  double* vy = mesh->createField<double>("vy", NODE_CENTERED);
  double* p = mesh->createField<double>("pressure", CELL_CENTERED);
  internal::set_node_fields(mesh->getNumberOfNodes(), vx, vy);
  internal::set_cell_fields(mesh->getNumberOfCells(), p);

  /* STEP 4: Get info and delete. */
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  const double* x = mesh->getCoordinateArray(X_COORDINATE);
  const double* y = mesh->getCoordinateArray(Y_COORDINATE);
  const IndexType* connec = mesh->getCellNodesArray();

  delete mesh;
  mesh = nullptr;

  /* STEP 5: get a mesh object from the group. */
  const Mesh* m = mint::getMesh(root);

  /* STEP 6: test the object. */
  EXPECT_TRUE(m->hasSidreGroup());
  EXPECT_EQ(m->getMeshType(), UNSTRUCTURED_MESH);
  EXPECT_TRUE(m->hasExplicitCoordinates());
  EXPECT_TRUE(m->hasExplicitConnectivity());
  EXPECT_FALSE(m->hasMixedCellTypes());
  EXPECT_EQ(m->getBlockId(), BLOCKID);
  EXPECT_EQ(m->getPartitionId(), PARTID);
  EXPECT_EQ(m->getDimension(), DIMENSION);
  EXPECT_EQ(m->getNumberOfNodes(), n_nodes);
  EXPECT_EQ(m->getNodeCapacity(), node_capacity);
  EXPECT_EQ(m->getNumberOfCells(), n_cells);
  EXPECT_EQ(m->getCellCapacity(), cell_capacity);
  EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
  EXPECT_EQ(m->getCoordinateArray(Y_COORDINATE), y);

  internal::check_append_nodes(m, X_EXTENT, Y_EXTENT, SPACING);
  internal::check_append_cells_single(m, X_EXTENT, Y_EXTENT);

  const double* vx_cpy = m->getFieldPtr<double>("vx", NODE_CENTERED);
  const double* vy_cpy = m->getFieldPtr<double>("vy", NODE_CENTERED);
  const double* p_cpy = m->getFieldPtr<double>("pressure", CELL_CENTERED);

  EXPECT_EQ(vx_cpy, vx);
  EXPECT_EQ(vy_cpy, vy);
  EXPECT_EQ(p_cpy, p);

  internal::check_node_fields(n_nodes, vx_cpy, vy_cpy);
  internal::check_cell_fields(n_cells, p_cpy);

  /* STEP 7: down-cast and test the UnstructuredMesh object. */
  const UnstructuredMesh<SINGLE_SHAPE>* M = dynamic_cast<const UnstructuredMesh<SINGLE_SHAPE>*>(m);
  EXPECT_EQ(M->getCellNodesArray(), connec);

  /* STEP 8: de-allocate. */
  delete m;
}

//------------------------------------------------------------------------------
TEST(mint_mesh, get_mixed_topology_unstructured_from_sidre)
{
  constexpr int DIMENSION = 2;
  constexpr IndexType X_EXTENT = 11;
  constexpr IndexType Y_EXTENT = 11;
  constexpr double SPACING = 1.0;
  constexpr int BLOCKID = 9;
  constexpr int PARTID = 10;

  // STEP 0: get empty Sidre group where to store the mesh
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* STEP 1: create the UnstructuredMesh */
  UnstructuredMesh<MIXED_SHAPE>* mesh = new UnstructuredMesh<MIXED_SHAPE>(DIMENSION, root);
  mesh->setBlockId(BLOCKID);
  mesh->setPartitionId(PARTID);

  /* STEP 2: Add the nodes and cells. */
  internal::append_nodes(mesh, X_EXTENT, Y_EXTENT, SPACING);
  internal::append_cells(mesh, X_EXTENT, Y_EXTENT);

  /* STEP 3: Add and fill the fields. */
  double* vx = mesh->createField<double>("vx", NODE_CENTERED);
  double* vy = mesh->createField<double>("vy", NODE_CENTERED);
  double* p = mesh->createField<double>("pressure", CELL_CENTERED);
  internal::set_node_fields(mesh->getNumberOfNodes(), vx, vy);
  internal::set_cell_fields(mesh->getNumberOfCells(), p);

  /* STEP 4: Get info and delete. */
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  const double* x = mesh->getCoordinateArray(X_COORDINATE);
  const double* y = mesh->getCoordinateArray(Y_COORDINATE);
  const IndexType connec_size = mesh->getCellNodesSize();
  const IndexType connec_capacity = mesh->getCellNodesCapacity();
  const IndexType* connec = mesh->getCellNodesArray();
  const IndexType* offsets = mesh->getCellNodesOffsetsArray();
  const CellType* types = mesh->getCellTypesArray();

  delete mesh;
  mesh = nullptr;

  /* STEP 5: get a mesh object from the group. */
  const Mesh* m = mint::getMesh(root);

  /* STEP 6: test the object. */
  EXPECT_TRUE(m->hasSidreGroup());
  EXPECT_EQ(m->getMeshType(), UNSTRUCTURED_MESH);
  EXPECT_TRUE(m->hasExplicitCoordinates());
  EXPECT_TRUE(m->hasExplicitConnectivity());
  EXPECT_TRUE(m->hasMixedCellTypes());
  EXPECT_EQ(m->getBlockId(), BLOCKID);
  EXPECT_EQ(m->getPartitionId(), PARTID);
  EXPECT_EQ(m->getDimension(), DIMENSION);
  EXPECT_EQ(m->getNumberOfNodes(), n_nodes);
  EXPECT_EQ(m->getNodeCapacity(), node_capacity);
  EXPECT_EQ(m->getNumberOfCells(), n_cells);
  EXPECT_EQ(m->getCellCapacity(), cell_capacity);
  EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
  EXPECT_EQ(m->getCoordinateArray(Y_COORDINATE), y);

  internal::check_append_nodes(m, X_EXTENT, Y_EXTENT, SPACING);
  internal::check_append_cells_mixed(m, X_EXTENT, Y_EXTENT);

  const double* vx_cpy = m->getFieldPtr<double>("vx", NODE_CENTERED);
  const double* vy_cpy = m->getFieldPtr<double>("vy", NODE_CENTERED);
  const double* p_cpy = m->getFieldPtr<double>("pressure", CELL_CENTERED);

  EXPECT_EQ(vx_cpy, vx);
  EXPECT_EQ(vy_cpy, vy);
  EXPECT_EQ(p_cpy, p);

  internal::check_node_fields(n_nodes, vx_cpy, vy_cpy);
  internal::check_cell_fields(n_cells, p_cpy);

  /* STEP 7: down-cast and test the UnstructuredMesh object. */
  const UnstructuredMesh<MIXED_SHAPE>* M = dynamic_cast<const UnstructuredMesh<MIXED_SHAPE>*>(m);
  EXPECT_EQ(M->getCellNodesSize(), connec_size);
  EXPECT_EQ(M->getCellNodesCapacity(), connec_capacity);
  EXPECT_EQ(M->getCellNodesArray(), connec);
  EXPECT_EQ(M->getCellNodesOffsetsArray(), offsets);
  EXPECT_EQ(M->getCellTypesArray(), types);

  /* STEP 8: de-allocate. */
  delete m;
}

//------------------------------------------------------------------------------
TEST(mint_mesh, immediate_ug_reserve)
{
  axom::sidre::DataStore objectDS;
  axom::sidre::Group* meshGroup = objectDS.getRoot()->createGroup("myGroup");
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> contourMesh(2,
                                                                     axom::mint::CellType::SEGMENT,
                                                                     meshGroup);
  contourMesh.reserveCells(10);  // This may crash.
}

#endif /* AXOM_MINT_USE_SIDRE */

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
