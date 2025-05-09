// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/config.hpp"  // for compile-time type definitions

// Mint includes
#include "axom/mint/mesh/blueprint.hpp"             // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"             // for CellTypes enum definition
#include "axom/mint/mesh/CurvilinearMesh.hpp"       // for CurivilinearMesh
#include "axom/mint/mesh/ParticleMesh.hpp"          // for ParticleMesh
#include "axom/mint/mesh/internal/MeshHelpers.hpp"  // for axom::mint::internal::dim
#include "StructuredMesh_helpers.hpp"               // for StructuredMesh test helpers

// Slic includes
#include "axom/slic.hpp"

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"

using namespace axom::mint;
using IndexType = axom::IndexType;

// globals
const char* IGNORE_OUTPUT = ".*";

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
//------------------------------------------------------------------------------
void set_coordinates(IndexType N, double* x, double* y = nullptr, double* z = nullptr)
{
  const IndexType ndims = axom::mint::internal::dim(x, y, z);
  double* const coords[3] = {x, y, z};

  double factor = axom::mint::internal::PI;
  for(int dim = 0; dim < ndims; ++dim)
  {
    SLIC_ASSERT(coords[dim] != nullptr);
    for(IndexType i = 0; i < N; ++i)
    {
      coords[dim][i] = factor * i;
    }

    factor *= axom::mint::internal::PI;
  }
}

//------------------------------------------------------------------------------
void set_coordinates(CurvilinearMesh* m)
{
  const int ndims = m->getDimension();
  const IndexType numNodes = m->getNumberOfNodes();
  double* x = m->getCoordinateArray(X_COORDINATE);

  if(ndims == 1)
  {
    set_coordinates(numNodes, x);
  }
  else if(ndims == 2)
  {
    double* y = m->getCoordinateArray(Y_COORDINATE);
    set_coordinates(numNodes, x, y);
  }
  else
  {
    double* y = m->getCoordinateArray(Y_COORDINATE);
    double* z = m->getCoordinateArray(Z_COORDINATE);
    set_coordinates(numNodes, x, y, z);
  }
}

//------------------------------------------------------------------------------
void check_coordinates(IndexType N,
                       int ndims,
                       const double* x,
                       const double* y = nullptr,
                       const double* z = nullptr)
{
  const double* const coords[3] = {x, y, z};

  double factor = axom::mint::internal::PI;
  for(int dim = 0; dim < ndims; ++dim)
  {
    ASSERT_NE(coords[dim], nullptr);
    for(IndexType i = 0; i < N; ++i)
    {
      EXPECT_DOUBLE_EQ(coords[dim][i], factor * i);
    }

    factor *= axom::mint::internal::PI;
  }
}

//------------------------------------------------------------------------------
void check_coordinates(const CurvilinearMesh* m)
{
  const int ndims = m->getDimension();
  const IndexType numNodes = m->getNumberOfNodes();
  const double* x = m->getCoordinateArray(X_COORDINATE);

  if(ndims == 1)
  {
    check_coordinates(numNodes, ndims, x);
  }
  else if(ndims == 2)
  {
    const double* y = m->getCoordinateArray(Y_COORDINATE);
    check_coordinates(numNodes, ndims, x, y);
  }
  else
  {
    const double* y = m->getCoordinateArray(Y_COORDINATE);
    const double* z = m->getCoordinateArray(Z_COORDINATE);
    check_coordinates(numNodes, ndims, x, y, z);
  }
}

}  // END namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_curvilinear_mesh_DeathTest, invalid_construction)
{
  const IndexType N[] = {5, 5, 5};
  double x[] = {0, 1, 2, 3, 4, 5};

  // check 2nd native constructor
  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(-1, N[1], N[2]), IGNORE_OUTPUT);

  // check external constructor
  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(5, nullptr), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(-1, x), IGNORE_OUTPUT);

#ifdef AXOM_MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  sidre::Group* valid_group = root->createGroup("mesh");
  sidre::Group* particle_mesh = root->createGroup("particle_mesh");
  ParticleMesh(3, 10, particle_mesh);

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(nullptr, ""), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(root, ""), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(particle_mesh, ""), IGNORE_OUTPUT);

  // check 2nd push constructor
  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(nullptr, N[0]), IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(CurvilinearMesh(valid_group, -1), IGNORE_OUTPUT);

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_curvilinear_mesh, native_constructor)
{
  constexpr int NDIMS = 3;
  const IndexType N[] = {5, 6, 7};
  const int64 extent[] = {0, 4, 10, 15, 7, 13};

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    CurvilinearMesh* m;
    switch(idim)
    {
    case 1:
      m = new CurvilinearMesh(N[0]);
      break;
    case 2:
      m = new CurvilinearMesh(N[0], N[1]);
      break;
    default:
      EXPECT_EQ(idim, 3);
      m = new CurvilinearMesh(N[0], N[1], N[2]);
    }  // END switch

    axom::mint::internal::check_constructor(m, STRUCTURED_CURVILINEAR_MESH, idim, N);
    EXPECT_FALSE(m->isExternal());
    EXPECT_FALSE(m->hasSidreGroup());
    m->setExtent(idim, extent);
    axom::mint::internal::check_node_extent(m, extent);
    set_coordinates(m);
    axom::mint::internal::check_create_fields(m);
    delete m;
  }  // END for all dimensions
}

//------------------------------------------------------------------------------
TEST(mint_mesh_curvilinear_mesh, external_constructor)
{
  constexpr int NDIMS = 3;
  const IndexType N[] = {5, 6, 7};
  const int64 extent[] = {0, 4, 10, 15, 7, 13};
  const IndexType maxNumNodes = N[0] * N[1] * N[2];

  double* x = new double[maxNumNodes];
  double* y = new double[maxNumNodes];
  double* z = new double[maxNumNodes];
  set_coordinates(maxNumNodes, x, y, z);

  IndexType curNumNodes = 1;
  IndexType curNumCells = 1;
  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    curNumNodes *= N[idim - 1];
    curNumCells *= N[idim - 1] - 1;
    CurvilinearMesh* m = nullptr;

    switch(idim)
    {
    case 1:
    {
      m = new CurvilinearMesh(N[0], x);
      EXPECT_EQ(x, m->getCoordinateArray(X_COORDINATE));
    }  // END 1D
    break;
    case 2:
    {
      m = new CurvilinearMesh(N[0], x, N[1], y);
      EXPECT_EQ(x, m->getCoordinateArray(X_COORDINATE));
      EXPECT_EQ(y, m->getCoordinateArray(Y_COORDINATE));
    }  // END 2D
    break;
    default:
    {
      m = new CurvilinearMesh(N[0], x, N[1], y, N[2], z);
      EXPECT_EQ(x, m->getCoordinateArray(X_COORDINATE));
      EXPECT_EQ(y, m->getCoordinateArray(Y_COORDINATE));
      EXPECT_EQ(z, m->getCoordinateArray(Z_COORDINATE));
    }  // END 3D
    }  // END switch

    check_coordinates(m);
    axom::mint::internal::check_constructor(m, STRUCTURED_CURVILINEAR_MESH, idim, N);
    m->setExtent(idim, extent);
    axom::mint::internal::check_node_extent(m, extent);

    EXPECT_FALSE(m->hasSidreGroup());
    EXPECT_TRUE(m->isExternal());
    EXPECT_EQ(m->getNumberOfNodes(), curNumNodes);
    EXPECT_EQ(m->getNumberOfCells(), curNumCells);

    delete m;
    m = nullptr;

    // ensure array buffers are persistent
    check_coordinates(curNumNodes, idim, x, y, z);
  }  // END for all dimensions

  delete[] x;
  delete[] y;
  delete[] z;
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE

TEST(mint_mesh_curvilinear_mesh, sidre_constructor)
{
  constexpr int NDIMS = 3;
  const IndexType N[] = {5, 6, 7};
  const int64 extent[] = {0, 4, 10, 15, 7, 13};

  IndexType numNodes = 1;
  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    numNodes *= N[idim - 1];

    // STEP 0: create a data-store with an groups
    sidre::DataStore ds;
    sidre::Group* root = ds.getRoot();
    sidre::Group* meshGroup = root->createGroup("mesh");

    // STEP 1: populate the mesh in sidre.
    CurvilinearMesh* m;
    switch(idim)
    {
    case 1:
      m = new CurvilinearMesh(meshGroup, N[I_DIRECTION]);
      break;
    case 2:
      m = new CurvilinearMesh(meshGroup, N[I_DIRECTION], N[J_DIRECTION]);
      break;
    default:
      EXPECT_EQ(idim, 3);
      m = new CurvilinearMesh(meshGroup, N[I_DIRECTION], N[J_DIRECTION], N[K_DIRECTION]);
    }  // END switch

    EXPECT_TRUE(m->hasSidreGroup());
    EXPECT_FALSE(m->isExternal());
    axom::mint::internal::check_constructor(m, STRUCTURED_CURVILINEAR_MESH, idim, N);
    m->setExtent(idim, extent);
    axom::mint::internal::check_node_extent(m, extent);
    set_coordinates(m);
    axom::mint::internal::check_create_fields(m);

    delete m;
    m = nullptr;

    // STEP 2: pull the mesh from sidre in to new instances.
    m = new CurvilinearMesh(meshGroup);
    EXPECT_TRUE(m->hasSidreGroup());
    EXPECT_FALSE(m->isExternal());
    axom::mint::internal::check_constructor(m, STRUCTURED_CURVILINEAR_MESH, idim, N);
    axom::mint::internal::check_fields(m, true);
    EXPECT_EQ(idim, m->getDimension());
    EXPECT_EQ(numNodes, m->getNumberOfNodes());
    check_coordinates(m);

    delete m;

    EXPECT_TRUE(blueprint::isValidRootGroup(meshGroup));

  }  // END for all dimensions
}

#endif /* AXOM_MINT_USE_SIDRE */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
