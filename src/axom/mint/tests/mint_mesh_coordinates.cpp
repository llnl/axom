// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic.hpp"

#include "axom/mint/config.hpp"
#include "axom/mint/mesh/MeshCoordinates.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
#endif

#include "gtest/gtest.h"

namespace utilities = axom::utilities;

namespace axom
{
namespace mint
{
#ifdef AXOM_MINT_USE_SIDRE
sidre::DataStore* ds = nullptr;
#endif

// constants used in tests
constexpr double PI = 3.14159265358979323846;
constexpr IndexType DEFAULT_CAPACITY = 100;
constexpr IndexType ZERO_NUM_NODES = 0;
constexpr IndexType SMALL_NUM_NODES = 4;
constexpr IndexType LARGE_NUM_NODES = 256;
constexpr IndexType IGNORE_CAPACITY = -1;
constexpr IndexType SMALL_NODE_CAPACITY = 5;
constexpr IndexType LARGE_NODE_CAPACITY = 256;

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace internal
{
/*!
 * \brief Checks that the values of two arrays are the same.
 * \param [in] actual pointer to the buffer of values to check
 * \param [in] expected pointer to the buffer consisting the expected values
 * \param [in] N the length of the buffer.
 *
 * \pre actual != nullptr
 * \pre expected != nullptr
 * \pre N > 0
 *
 */
void check_array_values(const double* actual, const double* expected, IndexType N)
{
  SLIC_ASSERT(actual != nullptr);
  SLIC_ASSERT(expected != nullptr);
  SLIC_ASSERT(N > 0);

  for(IndexType i = 0; i < N; ++i)
  {
    EXPECT_DOUBLE_EQ(actual[i], expected[i]);
  }
}

#ifdef AXOM_MINT_USE_SIDRE
/*!
 * \brief Create a valid sidre datastore to use in tests
 *
 * \param [in,out] ds pointer to the datastore
 * \param [in] dimension the requested dimension
 *
 * \pre 1 <= dimension <= 3
 * \note Caller must delete the datastore instance and group.
 */
void create_sidre_data(sidre::DataStore& ds, int dimension)
{
  SLIC_ASSERT((dimension >= 1) && (dimension <= 3));

  double x[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double y[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double z[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double* ptrs[3] = {x, y, z};

  sidre::Group* gp = ds.getRoot();
  SLIC_ASSERT(gp != nullptr);

  gp->createView("type")->setString("explicit");
  sidre::Group* values = gp->createGroup("values");

  const char* coord_names[3] = {"x", "y", "z"};

  for(int idim = 0; idim < dimension; ++idim)
  {
    const char* name = coord_names[idim];
    sidre::View* coord_view = values->createView(std::string(name));

    // NOTE: even though the array goes out-of-scope here, the data
    // remains persistent in sidre
    sidre::deprecated::MCArray<double> coord_array(coord_view, SMALL_NUM_NODES, 1, SMALL_NUM_NODES);

    coord_array.set(ptrs[idim], SMALL_NUM_NODES, 0);
  }  // END for all dimensions
}
#endif

/*!
 * \brief Tests that the coordinate arrays are not nullptr
 * \param [in] coords const pointer to the MeshCoordinates object.
 */
void check_coordinate_arrays(const MeshCoordinates* coords)
{
  SLIC_ASSERT(coords != nullptr);

  const int ndims = coords->dimension();
  for(int i = 0; i < ndims; ++i)
  {
    const double* coordsptr = coords->getCoordinateArray(i);
    EXPECT_TRUE(coordsptr != nullptr);
  }
}

/*!
 * \brief Tests construction of a MeshCoordinates object of specified dimension.
 * \param [in] dimension the mesh dimension
 */
void check_constructor(int dimension)
{
  EXPECT_TRUE(dimension >= 1 && dimension <= 3);

  // STEP 0: construct MeshCoordinates object
  MeshCoordinates coords(dimension);

  // STEP 1: check post-conditions
  EXPECT_EQ(dimension, coords.dimension());
  EXPECT_EQ(0, coords.numNodes());
  EXPECT_TRUE(coords.capacity() > 0);
  EXPECT_EQ(coords.capacity(), DEFAULT_CAPACITY);
  EXPECT_TRUE(coords.empty());

  // STEP 2: check coordinate arrays
  check_coordinate_arrays(&coords);
}

/*!
 * \brief Tests construction of a MeshCoordinates instance of specified
 *  dimension, number of nodes and optionally initial max capacity.
 *
 * \param [in] dimension the mesh dimension
 * \param [in] numNodes the number of nodes in the mesh
 * \param [in] capacity max initial capacity (optional)
 */
void check_constructor(int dimension, IndexType numNodes, IndexType capacity = IGNORE_CAPACITY)
{
  EXPECT_TRUE(dimension >= 1 && dimension <= 3);

  // STEP 0: construct the MeshCoordinates object
  MeshCoordinates* coords = nullptr;
  if(capacity == IGNORE_CAPACITY)
  {
    coords = new MeshCoordinates(dimension, numNodes);
  }
  else
  {
    coords = new MeshCoordinates(dimension, numNodes, capacity);
  }
  EXPECT_TRUE(coords != nullptr);

  // STEP 1: check post-conditions
  EXPECT_EQ(dimension, coords->dimension());
  EXPECT_EQ(numNodes, coords->numNodes());
  EXPECT_TRUE(coords->numNodes() <= coords->capacity());

  if(numNodes == ZERO_NUM_NODES)
  {
    EXPECT_TRUE(coords->empty());
  }

  // STEP 2: check actual capacity
  const IndexType actual_capacity = coords->capacity();

  const double ratio = deprecated::MCArray<double>::DEFAULT_RESIZE_RATIO;
  const IndexType expected_computed_capacity =
    utilities::max(DEFAULT_CAPACITY, static_cast<IndexType>(numNodes * ratio + 0.5));

  if(capacity == IGNORE_CAPACITY)
  {
    EXPECT_EQ(actual_capacity, expected_computed_capacity);
  }
  else
  {
    EXPECT_EQ(actual_capacity, capacity);
  }

  // STEP 3: check coordinate arrays
  check_coordinate_arrays(coords);

  // STEP 4: delete MeshCoordinates object
  delete coords;
  coords = nullptr;
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
 * \brief Check that the nodal coordinates are correct.
 *
 * \param [in] mesh_coords the MeshCoordinates to check.
 */
void check_append(const MeshCoordinates* mesh_coords)
{
  IndexType n_nodes = mesh_coords->numNodes();

  const int ndims = mesh_coords->dimension();
  const double* x = mesh_coords->getCoordinateArray(X_COORDINATE);

  /* Check using pointers */
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    EXPECT_EQ(x[i], getCoordValue(ndims, i, 0));
  }
  if(ndims > 1)
  {
    const double* y = mesh_coords->getCoordinateArray(Y_COORDINATE);
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      EXPECT_EQ(y[i], getCoordValue(ndims, i, 1));
    }
  }
  if(ndims > 2)
  {
    const double* z = mesh_coords->getCoordinateArray(Z_COORDINATE);
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      EXPECT_EQ(z[i], getCoordValue(ndims, i, 2));
    }
  }

  /* Check using getCoordinates */
  double coords[3];
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    mesh_coords->getCoordinates(i, coords);
    for(int dim = 0; dim < ndims; ++dim)
    {
      EXPECT_EQ(coords[dim], getCoordValue(ndims, i, dim));
    }
  }

  /* Check using getCoordinate */
  for(IndexType i = 0; i < n_nodes; ++i)
  {
    for(int dim = 0; dim < ndims; ++dim)
    {
      EXPECT_EQ(mesh_coords->getCoordinate(i, dim), getCoordValue(ndims, i, dim));
    }
  }
}

/*!
 * \brief Append a single node multiple times in a row.
 *
 * \param [in/out] mesh_coords the MeshCoordinates to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
void append_node_single(MeshCoordinates* mesh_coords, IndexType n_nodes)
{
  const int ndims = mesh_coords->dimension();
  IndexType cur_n_nodes = mesh_coords->numNodes();

  if(ndims == 1)
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      mesh_coords->append(getCoordValue(ndims, cur_n_nodes, 0));
      EXPECT_EQ(++cur_n_nodes, mesh_coords->numNodes());
    }
  }
  else if(ndims == 2)
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      mesh_coords->append(getCoordValue(ndims, cur_n_nodes, 0), getCoordValue(ndims, cur_n_nodes, 1));
      EXPECT_EQ(++cur_n_nodes, mesh_coords->numNodes());
    }
  }
  else
  {
    for(IndexType i = 0; i < n_nodes; ++i)
    {
      mesh_coords->append(getCoordValue(ndims, cur_n_nodes, 0),
                          getCoordValue(ndims, cur_n_nodes, 1),
                          getCoordValue(ndims, cur_n_nodes, 2));
      EXPECT_EQ(++cur_n_nodes, mesh_coords->numNodes());
    }
  }
}

/*!
 * \brief Append multiple nodes at once using the array of structs layout.
 *
 * \param [in/out] mesh_coords the MeshCoordinates to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
void append_node_structs(MeshCoordinates* mesh_coords, IndexType n_nodes)
{
  const int ndims = mesh_coords->dimension();
  IndexType cur_n_nodes = mesh_coords->numNodes();
  double* coords = new double[ndims * n_nodes];

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    for(IndexType dim = 0; dim < ndims; ++dim)
    {
      coords[ndims * i + dim] = getCoordValue(ndims, i + cur_n_nodes, dim);
    }
  }

  mesh_coords->append(coords, n_nodes);
  cur_n_nodes += n_nodes;
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  delete[] coords;
}

/*!
 * \brief Append multiple nodes at once using the struct of arrays layout.
 *
 * \param [in/out] mesh_coords the MeshCoordinates to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
void append_node_arrays(MeshCoordinates* mesh_coords, IndexType n_nodes)
{
  const int ndims = mesh_coords->dimension();
  IndexType cur_n_nodes = mesh_coords->numNodes();
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
    mesh_coords->append(coords, n_nodes);
    cur_n_nodes += n_nodes;
    EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());
  }
  else if(ndims == 2)
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    mesh_coords->append(x, y, n_nodes);
    cur_n_nodes += n_nodes;
    EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());
  }
  else
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    const double* z = coords + 2 * n_nodes;
    mesh_coords->append(x, y, z, n_nodes);
    cur_n_nodes += n_nodes;
    EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());
  }

  delete[] coords;
}

/*!
 * \brief Append nodes using all the various methods.
 *
 * \param [in/out] mesh_coords the MeshCoordinates to append to.
 * \param [in] n_nodes the number of nodes to append.
 */
void append_nodes(MeshCoordinates* mesh_coords, IndexType n_nodes)
{
  ASSERT_EQ(n_nodes % 3, 0);

  IndexType cur_n_nodes = mesh_coords->numNodes();
  ASSERT_EQ(cur_n_nodes, 0);
  EXPECT_TRUE(mesh_coords->empty());

  /* Append one node at a time */
  append_node_single(mesh_coords, n_nodes / 3);
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  /* Append multiple nodes at once using the array of structs layout. */
  append_node_structs(mesh_coords, n_nodes / 3);
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ(cur_n_nodes, 2 * n_nodes / 3);
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  /* Append multiple nodes at once using the struct of arrays layout. */
  append_node_arrays(mesh_coords, n_nodes / 3);
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ(cur_n_nodes, n_nodes);
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  check_append(mesh_coords);
}

/*!
 * \brief Insert a single node.
 *
 * \param [in/out] mesh_coords the MeshCoordinates in question.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this node after all
 *  insertions have taken place.
 */
void insert_single(MeshCoordinates* mesh_coords, IndexType pos, IndexType final_pos)
{
  const int ndims = mesh_coords->dimension();
  IndexType cur_n_nodes = mesh_coords->numNodes();

  if(ndims == 1)
  {
    mesh_coords->insert(pos, getCoordValue(ndims, final_pos, 0));
  }
  else if(ndims == 2)
  {
    mesh_coords->insert(pos, getCoordValue(ndims, final_pos, 0), getCoordValue(ndims, final_pos, 1));
  }
  else
  {
    mesh_coords->insert(pos,
                        getCoordValue(ndims, final_pos, 0),
                        getCoordValue(ndims, final_pos, 1),
                        getCoordValue(ndims, final_pos, 2));
  }

  EXPECT_EQ(++cur_n_nodes, mesh_coords->numNodes());
}

/*!
 * \brief Insert multiple nodes using the array of structs layout.
 *
 * \param [in/out] mesh_coords the MeshCoordinates in question.
 * \param [in] n_nodes the number of nodes to insert.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this node after all
 *  insertions have taken place.
 */
void insert_structs(MeshCoordinates* mesh_coords, IndexType n_nodes, IndexType pos, IndexType final_pos)
{
  const int ndims = mesh_coords->dimension();
  IndexType cur_n_nodes = mesh_coords->numNodes();
  double* coords = new double[ndims * n_nodes];

  for(IndexType i = 0; i < n_nodes; ++i)
  {
    for(IndexType dim = 0; dim < ndims; ++dim)
    {
      coords[ndims * i + dim] = getCoordValue(ndims, final_pos + i, dim);
    }
  }

  mesh_coords->insert(pos, coords, n_nodes);
  cur_n_nodes += n_nodes;
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  delete[] coords;
}

/*!
 * \brief Insert multiple nodes using the struct of arrays layout.
 *
 * \param [in/out] mesh_coords the MeshCoordinates in question.
 * \param [in] n_nodes the number of nodes to insert.
 * \param [in] pos the insertion position.
 * \param [in] final_pos the expected final position of this node after all
 *  insertions have taken place.
 */
void insert_arrays(MeshCoordinates* mesh_coords, IndexType n_nodes, IndexType pos, IndexType final_pos)
{
  const int ndims = mesh_coords->dimension();
  IndexType cur_n_nodes = mesh_coords->numNodes();
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
    mesh_coords->insert(pos, coords, n_nodes);
  }
  else if(ndims == 2)
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    mesh_coords->insert(pos, x, y, n_nodes);
  }
  else
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    const double* z = coords + 2 * n_nodes;
    mesh_coords->insert(pos, x, y, z, n_nodes);
  }

  cur_n_nodes += n_nodes;
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  delete[] coords;
}

/*!
 * \brief Insert nodes into the MeshCoordinates.
 *
 * \param [in/out] mesh_coords the MeshCoordinates in question.
 * \param [in] n_nodes the number of nodes to insert.
 */
void insert_nodes(MeshCoordinates* mesh_coords, IndexType n_nodes)
{
  ASSERT_EQ(n_nodes % 9, 0);

  IndexType cur_n_nodes = mesh_coords->numNodes();
  ASSERT_EQ(cur_n_nodes, 0);
  EXPECT_TRUE(mesh_coords->empty());

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
      insert_single(mesh_coords, insert_positions[i], final_pos[i]);
      cur_n_nodes++;
      EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());
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
  insert_structs(mesh_coords, ninth_n_nodes, 0, ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  insert_structs(mesh_coords, ninth_n_nodes, third_n_nodes, 4 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  insert_structs(mesh_coords, ninth_n_nodes, 5 * ninth_n_nodes, 7 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  EXPECT_EQ(cur_n_nodes, 2 * third_n_nodes);
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  /* Insert multiple nodes at once using the struct of arrays layout. */
  insert_arrays(mesh_coords, ninth_n_nodes, 0, 0);
  cur_n_nodes += ninth_n_nodes;
  insert_arrays(mesh_coords, ninth_n_nodes, 5 * ninth_n_nodes, 5 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  insert_arrays(mesh_coords, ninth_n_nodes, 8 * ninth_n_nodes, 8 * ninth_n_nodes);
  cur_n_nodes += ninth_n_nodes;
  EXPECT_EQ(cur_n_nodes, n_nodes);
  EXPECT_EQ(cur_n_nodes, mesh_coords->numNodes());

  check_append(mesh_coords);
}

/*!
 * \brief Test set/get a node from the given MeshCoordinates object.
 * \param [in] mc the mesh coordinates object to test
 * \pre mc != nullptr
 */
void check_set_and_get(MeshCoordinates* mc)
{
  EXPECT_TRUE(mc != nullptr);

  constexpr double TEST_VALUE = 7;
  constexpr IndexType targetIdx = 3;

  const int nnodes = mc->numNodes();
  const int ndims = mc->dimension();

  ASSERT_TRUE(targetIdx < nnodes);

  constexpr double set_value = TEST_VALUE + 0.5;
  if(ndims == 1)
  {
    mc->set(targetIdx, set_value);
  }
  else if(ndims == 2)
  {
    mc->set(targetIdx, set_value, set_value);
  }
  else
  {
    mc->set(targetIdx, set_value, set_value, set_value);
  }

  for(int i = 0; i < ndims; ++i)
  {
    EXPECT_DOUBLE_EQ(mc->getCoordinate(targetIdx, i), set_value);
  }
}

/*!
 * \brief Return a pointer to a new MeshCoordinates.
 *
 * \param [in] ndims the dimension of the mesh to create.
 * \param [in] n_nodes space to allocate for nodes.
 */
MeshCoordinates* createExternal(int dim, IndexType n_nodes)
{
  double* x = new double[n_nodes];
  double* y = nullptr;
  double* z = nullptr;

  if(dim > 1)
  {
    y = new double[n_nodes];
  }
  if(dim > 2)
  {
    z = new double[n_nodes];
  }

  return new MeshCoordinates(0, n_nodes, x, y, z);
}

/*!
 * \brief Fill the given array with MeshCoordinates.
 *
 * \param [in/out] mesh_coords the array to fill.
 * \param [in] n_nodes the nodal capacity for the external meshes.
 *
 * \note if using sidre 9 meshes are created otherwise only 6.
 */
void createCoords(MeshCoordinates** mesh_coords, IndexType n_nodes)
{
#ifdef AXOM_MINT_USE_SIDRE
  SLIC_ERROR_IF(ds != nullptr, "Did not free DataStore.");
  ds = new sidre::DataStore();
  sidre::Group* root = ds->getRoot();
#endif

  int cur_coords = 0;
  for(int dim = 1; dim <= 3; ++dim)
  {
    mesh_coords[cur_coords++] = new MeshCoordinates(dim);

    mesh_coords[cur_coords++] = createExternal(dim, n_nodes);

#ifdef AXOM_MINT_USE_SIDRE
    const std::string coordset_name = "c" + std::to_string(dim);
    sidre::Group* group = root->createGroup(coordset_name);
    mesh_coords[cur_coords++] = new MeshCoordinates(group, dim);
#endif
  }
}

/*!
 * \brief Delete the given external MeshCoordinates and its associated arrays.
 *
 * \param [in/out] mesh_coords the MeshCoordinates to delete.
 */
void deleteExternal(MeshCoordinates*& mesh_coords)
{
  const int ndims = mesh_coords->dimension();
  const double* x = mesh_coords->getCoordinateArray(X_COORDINATE);
  const double* y = nullptr;
  const double* z = nullptr;

  if(ndims > 1)
  {
    y = mesh_coords->getCoordinateArray(Y_COORDINATE);
  }
  if(ndims > 2)
  {
    z = mesh_coords->getCoordinateArray(Z_COORDINATE);
  }

  delete mesh_coords;
  mesh_coords = nullptr;

  delete[] x;
  delete[] y;
  delete[] z;
}

/*!
 * \brief Free the given array of MeshCoordinates.
 *
 * \param [in/out] mesh_coords the array of MeshCoordinates.
 * \param [in] n_coords the nuber of MeshCoordinates in the array.
 */
void deleteCoords(MeshCoordinates** mesh_coords, int n_coords)
{
  for(int i = 0; i < n_coords; ++i)
  {
    if(mesh_coords[i]->isExternal())
    {
      internal::deleteExternal(mesh_coords[i]);
    }
    else
    {
      delete mesh_coords[i];
      mesh_coords[i] = nullptr;
    }
  }

#ifdef AXOM_MINT_USE_SIDRE
  delete ds;
  ds = nullptr;
#endif
}

/*!
 * \brief Delete the given external MeshCoordinates and create a copy from it's
 *  external arrays.
 *
 * \param [in/out] mesh_coords the MeshCoordinates to delete, it is replaced
 *  with the copy.
 */
void deleteAndDuplicateExternal(MeshCoordinates*& mesh_coords)
{
  ASSERT_TRUE(mesh_coords->isExternal());
  const int ndims = mesh_coords->dimension();
  const int n_nodes = mesh_coords->numNodes();
  const int capacity = mesh_coords->capacity();
  double* x = mesh_coords->getCoordinateArray(X_COORDINATE);
  double* y = nullptr;
  double* z = nullptr;

  if(ndims > 1)
  {
    y = mesh_coords->getCoordinateArray(Y_COORDINATE);
  }
  if(ndims > 2)
  {
    z = mesh_coords->getCoordinateArray(Z_COORDINATE);
  }

  delete mesh_coords;

  mesh_coords = new MeshCoordinates(n_nodes, capacity, x, y, z);
  EXPECT_EQ(ndims, mesh_coords->dimension());
  EXPECT_EQ(n_nodes, mesh_coords->numNodes());
  EXPECT_EQ(capacity, mesh_coords->capacity());
  EXPECT_EQ(x, mesh_coords->getCoordinateArray(X_COORDINATE));
  if(ndims > 1)
  {
    EXPECT_EQ(y, mesh_coords->getCoordinateArray(Y_COORDINATE));
  }
  if(ndims > 2)
  {
    EXPECT_EQ(z, mesh_coords->getCoordinateArray(Z_COORDINATE));
  }
}

#ifdef AXOM_MINT_USE_SIDRE

/*!
 * \brief Delete the given sidre MeshCoordinates and create a copy from
 *  it's group.
 *
 * \param [in/out] mesh_coords the MeshCoordinates to delete, it is replaced
 *  with the copy.
 */
void deleteAndDuplicateSidre(MeshCoordinates*& mesh_coords)
{
  ASSERT_TRUE(mesh_coords->isInSidre());

  const int ndims = mesh_coords->dimension();
  const int n_nodes = mesh_coords->numNodes();
  const int capacity = mesh_coords->capacity();

  sidre::Group* group = mesh_coords->getSidreGroup();

  delete mesh_coords;
  mesh_coords = new MeshCoordinates(group);

  EXPECT_EQ(ndims, mesh_coords->dimension());
  EXPECT_EQ(n_nodes, mesh_coords->numNodes());
  EXPECT_EQ(capacity, mesh_coords->capacity());
}

#endif

/*!
 * \brief Check that meshes restored from sidre or external buffers are
 *  equivalent to the original mesh.
 *
 * \param [in/out] meshes the array of meshes.
 * \param [in] n_coords the number of meshes.
 */
void check_restoration(MeshCoordinates** mesh_coords, int n_coords)
{
  for(int i = 0; i < n_coords; ++i)
  {
    if(mesh_coords[i]->isExternal())
    {
      deleteAndDuplicateExternal(mesh_coords[i]);
    }
#ifdef AXOM_MINT_USE_SIDRE
    else if(mesh_coords[i]->isInSidre())
    {
      deleteAndDuplicateSidre(mesh_coords[i]);
    }
#endif
    else
    {
      continue;
    }

    check_append(mesh_coords[i]);
  }
}

}  // namespace internal

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, append)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;

  constexpr int N_COORDS = STRIDE * 3;
  MeshCoordinates* mesh_coords[N_COORDS];

  internal::createCoords(mesh_coords, N_NODES);

  for(int i = 0; i < N_COORDS; ++i)
  {
    internal::append_nodes(mesh_coords[i], N_NODES);
  }

  internal::check_restoration(mesh_coords, N_COORDS);

  internal::deleteCoords(mesh_coords, N_COORDS);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, insert)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 900;

  constexpr int N_COORDS = STRIDE * 3;
  MeshCoordinates* mesh_coords[N_COORDS];

  internal::createCoords(mesh_coords, N_NODES);

  for(int i = 0; i < N_COORDS; ++i)
  {
    internal::insert_nodes(mesh_coords[i], N_NODES);
  }

  internal::check_restoration(mesh_coords, N_COORDS);

  internal::deleteCoords(mesh_coords, N_COORDS);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, set_and_get)
{
  constexpr int NDIMS = 3;

  double x[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double y[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double z[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};

  for(int dim = 1; dim <= NDIMS; ++dim)
  {
    // MeshCoordinates object with native store
    MeshCoordinates mc1(dim, SMALL_NUM_NODES);
    internal::check_set_and_get(&mc1);

#ifdef AXOM_MINT_USE_SIDRE
    // MeshCoordinates object from sidre store
    sidre::DataStore ds;
    MeshCoordinates mc2(ds.getRoot(), dim, SMALL_NUM_NODES, SMALL_NUM_NODES);
    internal::check_set_and_get(&mc2);
#endif

    // MeshCoordinates object tied to external buffers
    if(dim == 1)
    {
      MeshCoordinates mc3(SMALL_NUM_NODES, SMALL_NUM_NODES, x);
      internal::check_set_and_get(&mc3);
    }
    else if(dim == 2)
    {
      MeshCoordinates mc3(SMALL_NUM_NODES, SMALL_NUM_NODES, x, y);
      internal::check_set_and_get(&mc3);
    }
    else
    {
      EXPECT_EQ(dim, 3);
      MeshCoordinates mc3(SMALL_NUM_NODES, SMALL_NUM_NODES, x, y, z);
      internal::check_set_and_get(&mc3);
    }
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, reserve)
{
  constexpr int NDIMS = 3;

  for(int dim = 1; dim <= NDIMS; ++dim)
  {
    // Construct an empty MeshCoordinates object
    MeshCoordinates mc1(dim);
    IndexType numNodes1 = mc1.numNodes();

    mc1.reserve(LARGE_NODE_CAPACITY);

    EXPECT_EQ(mc1.numNodes(), numNodes1);
    EXPECT_EQ(mc1.capacity(), LARGE_NODE_CAPACITY);

#ifdef AXOM_MINT_USE_SIDRE
    // Construct a MeshCoordinates object from sidre
    sidre::DataStore ds;
    internal::create_sidre_data(ds, dim);
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE(coords_group != nullptr);

    MeshCoordinates mc2(coords_group);
    IndexType numNodes2 = mc2.numNodes();

    mc2.reserve(LARGE_NODE_CAPACITY);

    EXPECT_EQ(mc2.numNodes(), numNodes2);
    EXPECT_EQ(mc2.capacity(), LARGE_NODE_CAPACITY);
#endif
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coorindates, resize)
{
  constexpr int NDIMS = 3;

  for(int dim = 1; dim <= NDIMS; ++dim)
  {
    // Construct an empty MeshCoordinates object
    MeshCoordinates mc1(dim);
    mc1.resize(LARGE_NUM_NODES);
    EXPECT_EQ(mc1.numNodes(), LARGE_NUM_NODES);
    EXPECT_TRUE(mc1.numNodes() <= mc1.capacity());

#ifdef AXOM_MINT_USE_SIDRE
    // Construct a MeshCoordinates object from sidre
    sidre::DataStore ds;
    internal::create_sidre_data(ds, dim);
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE(coords_group != nullptr);

    MeshCoordinates mc2(coords_group);
    mc2.resize(LARGE_NUM_NODES);
    EXPECT_EQ(mc2.numNodes(), LARGE_NUM_NODES);
    EXPECT_TRUE(mc2.numNodes() <= mc2.capacity());
#endif
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, shrink)
{
  constexpr int NDIMS = 3;

  for(int dim = 1; dim <= NDIMS; ++dim)
  {
    // test shrink() when constructed using native store
    MeshCoordinates mc1(dim, SMALL_NUM_NODES);
    EXPECT_EQ(mc1.numNodes(), SMALL_NUM_NODES);
    EXPECT_TRUE(mc1.numNodes() < mc1.capacity());

    mc1.shrink();

    EXPECT_EQ(mc1.numNodes(), mc1.capacity());

#ifdef AXOM_MINT_USE_SIDRE
    // test shrink() when constructed using sidre
    sidre::DataStore ds;
    internal::create_sidre_data(ds, dim);
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE(coords_group != nullptr);

    MeshCoordinates mc2(coords_group);

    mc2.shrink();

    EXPECT_EQ(mc2.numNodes(), mc2.capacity());
#endif
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, change_resize_ratio)
{
  constexpr int NDIMS = 3;
  constexpr double DEFAULT_RESIZE_RATIO = deprecated::MCArray<double>::DEFAULT_RESIZE_RATIO;
  constexpr double NEW_RESIZE_RATIO = 2.5;

  MeshCoordinates mc(NDIMS);
  EXPECT_DOUBLE_EQ(mc.getResizeRatio(), DEFAULT_RESIZE_RATIO);

  mc.setResizeRatio(NEW_RESIZE_RATIO);
  EXPECT_DOUBLE_EQ(mc.getResizeRatio(), NEW_RESIZE_RATIO);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, native_constructors)
{
  for(int dim = 1; dim <= 3; ++dim)
  {
    internal::check_constructor(dim);

    internal::check_constructor(dim, ZERO_NUM_NODES);
    internal::check_constructor(dim, ZERO_NUM_NODES, SMALL_NODE_CAPACITY);
    internal::check_constructor(dim, ZERO_NUM_NODES, LARGE_NODE_CAPACITY);

    internal::check_constructor(dim, SMALL_NUM_NODES);
    internal::check_constructor(dim, SMALL_NUM_NODES, SMALL_NODE_CAPACITY);
    internal::check_constructor(dim, SMALL_NUM_NODES, LARGE_NODE_CAPACITY);

    internal::check_constructor(dim, LARGE_NUM_NODES);
    internal::check_constructor(dim, LARGE_NUM_NODES, LARGE_NODE_CAPACITY);
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, external_constructor)
{
  double x[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double y[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double z[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};

  // Test 1-D
  {
    constexpr int EXPECTED_DIMENSION = 1;
    MeshCoordinates coords(SMALL_NUM_NODES, SMALL_NUM_NODES, x);
    internal::check_coordinate_arrays(&coords);
    EXPECT_EQ(coords.dimension(), EXPECTED_DIMENSION);
    EXPECT_EQ(coords.numNodes(), SMALL_NUM_NODES);
    EXPECT_EQ(coords.capacity(), SMALL_NUM_NODES);
    EXPECT_EQ(coords.getCoordinateArray(X_COORDINATE), x);
  }
  for(int i = 0; i < SMALL_NUM_NODES; ++i)
  {
    EXPECT_DOUBLE_EQ(x[i], i + 1);
    EXPECT_DOUBLE_EQ(x[i], i + 1);
    EXPECT_DOUBLE_EQ(x[i], i + 1);
  }

  // Test 2-D
  {
    constexpr int EXPECTED_DIMENSION = 2;
    MeshCoordinates coords(SMALL_NUM_NODES, SMALL_NUM_NODES, x, y);
    internal::check_coordinate_arrays(&coords);
    EXPECT_EQ(coords.dimension(), EXPECTED_DIMENSION);
    EXPECT_EQ(coords.numNodes(), SMALL_NUM_NODES);
    EXPECT_EQ(coords.capacity(), SMALL_NUM_NODES);
    EXPECT_EQ(coords.getCoordinateArray(X_COORDINATE), x);
    EXPECT_EQ(coords.getCoordinateArray(Y_COORDINATE), y);
  }
  for(int i = 0; i < SMALL_NUM_NODES; ++i)
  {
    EXPECT_DOUBLE_EQ(x[i], i + 1);
    EXPECT_DOUBLE_EQ(x[i], i + 1);
    EXPECT_DOUBLE_EQ(x[i], i + 1);
  }

  // Test 3-D
  {
    constexpr int EXPECTED_DIMENSION = 3;
    MeshCoordinates coords(SMALL_NUM_NODES, SMALL_NUM_NODES, x, y, z);
    internal::check_coordinate_arrays(&coords);
    EXPECT_EQ(coords.dimension(), EXPECTED_DIMENSION);
    EXPECT_EQ(coords.numNodes(), SMALL_NUM_NODES);
    EXPECT_EQ(coords.capacity(), SMALL_NUM_NODES);
    EXPECT_EQ(coords.getCoordinateArray(X_COORDINATE), x);
    EXPECT_EQ(coords.getCoordinateArray(Y_COORDINATE), y);
    EXPECT_EQ(coords.getCoordinateArray(Z_COORDINATE), z);
  }
  for(int i = 0; i < SMALL_NUM_NODES; ++i)
  {
    EXPECT_DOUBLE_EQ(x[i], i + 1);
    EXPECT_DOUBLE_EQ(x[i], i + 1);
    EXPECT_DOUBLE_EQ(x[i], i + 1);
  }
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE

TEST(mint_mesh_coordinates, sidre_pull_constructor)
{
  double x[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double y[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double z[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double* expected_data[3] = {x, y, z};
  const char* coord_names[3] = {"x", "y", "z"};

  for(int dim = 1; dim <= 3; ++dim)
  {
    // STEP 0: create a test datastore with a coordinates groupp
    sidre::DataStore ds;
    internal::create_sidre_data(ds, dim);
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE(coords_group != nullptr);

    // STEP 1: create a MeshCoordinates instance within a scope from
    // the sidre group.
    // BEGIN SCOPE
    {
      MeshCoordinates coords(coords_group);
      EXPECT_EQ(coords.dimension(), dim);
      internal::check_coordinate_arrays(&coords);

      EXPECT_FALSE(coords.empty());
      EXPECT_EQ(coords.numNodes(), SMALL_NUM_NODES);
      EXPECT_TRUE(coords.numNodes() <= coords.capacity());

      internal::check_coordinate_arrays(&coords);

      for(int j = 0; j < dim; ++j)
      {
        internal::check_array_values(coords.getCoordinateArray(j), expected_data[j], SMALL_NUM_NODES);
      }
    }
    // END SCOPE

    // STEP 2: ensure data is persistent in the data-store after the
    // MeshCoordinates object goes out-of-scope
    EXPECT_TRUE(coords_group != nullptr);

    // Ensure that the data remains persistent in the data-store
    sidre::Group* values_group = coords_group->getGroup("values");
    EXPECT_TRUE(values_group != nullptr);
    for(int j = 0; j < dim; ++j)
    {
      sidre::View* coords_view = values_group->getView(std::string(coord_names[j]));

      double* coord_data = static_cast<double*>(coords_view->getVoidPtr());

      internal::check_array_values(coord_data, expected_data[j], SMALL_NUM_NODES);
    }

  }  // END for all dimensions
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates, sidre_push_constructor)
{
  double x[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double y[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double z[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double* data[3] = {x, y, z};
  const char* coord_names[3] = {"x", "y", "z"};

  for(int dim = 1; dim <= 3; ++dim)
  {
    sidre::DataStore ds;
    sidre::Group* coords_group = ds.getRoot();
    SLIC_ASSERT(coords_group != nullptr);

    EXPECT_TRUE(coords_group->getNumGroups() == 0);
    EXPECT_TRUE(coords_group->getNumViews() == 0);

    // BEGIN SCOPE
    {
      MeshCoordinates mesh_coords(coords_group, dim, SMALL_NUM_NODES);
      internal::check_coordinate_arrays(&mesh_coords);

      // ensure the sidre::Group is populated accordingly
      EXPECT_TRUE(coords_group->getNumGroups() == 1);
      EXPECT_TRUE(coords_group->getNumViews() == 1);
      EXPECT_TRUE(coords_group->hasView("type"));
      EXPECT_TRUE(coords_group->hasChildGroup("values"));

      sidre::Group* values_group = coords_group->getGroup("values");
      EXPECT_NE(values_group, nullptr);
      EXPECT_EQ(values_group->getNumViews(), dim);

      for(int j = 0; j < dim; ++j)
      {
        EXPECT_TRUE(values_group->hasChildView(coord_names[j]));

        sidre::View* coord_view = values_group->getView(coord_names[j]);
        EXPECT_TRUE(coord_view != nullptr);

        EXPECT_EQ(mesh_coords.getCoordinateArray(j), coord_view->getVoidPtr());
      }

      EXPECT_EQ(mesh_coords.dimension(), dim);
      EXPECT_EQ(mesh_coords.numNodes(), SMALL_NUM_NODES);
      EXPECT_TRUE(mesh_coords.numNodes() <= mesh_coords.capacity());

      IndexType capacity = SMALL_NUM_NODES * mesh_coords.getResizeRatio() + 0.5;
      if(capacity < MCArray<IndexType>::MIN_DEFAULT_CAPACITY)
      {
        capacity = MCArray<IndexType>::MIN_DEFAULT_CAPACITY;
      }
      EXPECT_EQ(mesh_coords.capacity(), capacity);

      // populate the coordinates, writes to the corresponding sidre views
      axom::deprecated::MCArray<double> xx(dim, 1, dim);
      for(int inode = 0; inode < SMALL_NUM_NODES; ++inode)
      {
        for(int j = 0; j < dim; ++j)
        {
          xx(j) = data[j][inode];
        }

        if(dim == 1)
        {
          mesh_coords.set(inode, xx[0]);
        }
        else if(dim == 2)
        {
          mesh_coords.set(inode, xx[0], xx[1]);
        }
        else
        {
          mesh_coords.set(inode, xx[0], xx[1], xx[2]);
        }
      }
    }
    // END SCOPE

    // Ensure the data is persistent in sidre
    EXPECT_TRUE(coords_group->getNumGroups() == 1);
    EXPECT_TRUE(coords_group->getNumViews() == 1);
    EXPECT_TRUE(coords_group->hasView("type"));
    EXPECT_TRUE(coords_group->hasChildGroup("values"));

    sidre::Group* values_group = coords_group->getGroup("values");
    EXPECT_NE(values_group, nullptr);
    EXPECT_EQ(values_group->getNumViews(), dim);

    for(int j = 0; j < dim; ++j)
    {
      EXPECT_TRUE(values_group->hasChildView(coord_names[j]));

      sidre::View* coord_view = values_group->getView(coord_names[j]);
      EXPECT_TRUE(coord_view != nullptr);

      double* actual_data = static_cast<double*>(coord_view->getVoidPtr());
      internal::check_array_values(actual_data, data[j], SMALL_NUM_NODES);
    }

  }  // END for all dimensions
}

#endif /* AXOM_MINT_USE_SIDRE */

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates_DeathTest, invalid_operations)
{
  // NOTE: this test ensures that when constructing a MeshCoordinates
  // object with

  const char* IGNORE_OUTPUT = ".*";
  constexpr IndexType N = 4;
  double x[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double y[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  double z[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};

  MeshCoordinates coords(N, N, x, y, z);

  EXPECT_DEATH_IF_SUPPORTED(coords.append(1, 1, 1), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(coords.resize(10), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(coords.reserve(10), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_coordinates_DeathTest, invalid_construction)
{
  const char* IGNORE_OUTPUT = ".*";

  // STEP 0: test construction with invalid dimension
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(0), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(4), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(0, SMALL_NUM_NODES, SMALL_NODE_CAPACITY), IGNORE_OUTPUT);

  // STEP 1: test construction with invalid numNodes,max_capacity settings
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(2, LARGE_NUM_NODES, SMALL_NODE_CAPACITY), IGNORE_OUTPUT);

  // STEP 2: test invalid construction with null external buffers
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(10, 10, nullptr), IGNORE_OUTPUT);

  // STEP 3: test invalid construction from external buffers with 0 nodes
  double x[SMALL_NUM_NODES] = {1.0, 2.0, 3.0, 4.0};
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(ZERO_NUM_NODES, ZERO_NUM_NODES, x), IGNORE_OUTPUT);

#ifdef AXOM_MINT_USE_SIDRE

  // STEP 4: test construction with a null Sidre group
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(nullptr), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(nullptr, 2, 10, 10), IGNORE_OUTPUT);

  // STEP 5: test sidre pull-constructor that does not conform to blueprint
  sidre::DataStore ds;
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(ds.getRoot()), IGNORE_OUTPUT);

  // STEP 6: test sidre push-constructor with with a non-empty group
  internal::create_sidre_data(ds, 2);
  EXPECT_DEATH_IF_SUPPORTED(MeshCoordinates(ds.getRoot(), 2, 5, 5), IGNORE_OUTPUT);

#endif /* AXOM_MINT_USE_SIDRE */
}

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
