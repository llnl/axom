// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/config.hpp"  // for compile-time definitions

// Mint includes
#include "axom/mint/mesh/blueprint.hpp"        // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"        // for CellType enum
#include "axom/mint/mesh/RectilinearMesh.hpp"  // for RectilinearMesh
#include "axom/mint/mesh/ParticleMesh.hpp"     // for ParticleMesh
#include "StructuredMesh_helpers.hpp"          // for StructuredMesh test helpers

// Slic includes
#include "axom/slic.hpp"

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>  // for exp
using namespace axom::mint;
using IndexType = axom::IndexType;

// globals
const char* IGNORE_OUTPUT = ".*";

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
void exponential_distribution(double origin, IndexType N, double* x)
{
  EXPECT_TRUE(x != nullptr);

  constexpr double beta = 0.05;
  const double expbeta = exp(beta);
  const double invf = 1 / (expbeta - 1.0);

  x[0] = origin;
  for(int i = 1; i < N; ++i)
  {
    const double prev = x[i - 1];
    const double dx = (exp(i * beta) - 1.0) * invf;
    x[i] = prev + dx;
  }
}

//------------------------------------------------------------------------------
void check_coordinate(const double* x, const double* expected, IndexType N)
{
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(expected != nullptr);
  EXPECT_TRUE(N >= 0);

  for(IndexType i = 0; i < N; ++i)
  {
    EXPECT_DOUBLE_EQ(x[i], expected[i]);
  }
}

//------------------------------------------------------------------------------
void check_fill_coords(RectilinearMesh* m)
{
  EXPECT_TRUE(m != nullptr);

  const int ndims = m->getDimension();
  for(int idim = 0; idim < ndims; ++idim)
  {
    const IndexType N = m->getNodeResolution(idim);
    double* x = m->getCoordinateArray(idim);
    exponential_distribution(42.0, N, x);
  }
}

//------------------------------------------------------------------------------

}  // END anonymous namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_rectilinear_mesh_DeathTest, invalid_construction)
{
  const IndexType N[] = {5, 5, 5};
  double x[5];

  // check 1st native constructor
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(2, nullptr), IGNORE_OUTPUT);

  // check 2nd native constructor
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(-1, N[1], N[2]), IGNORE_OUTPUT);

  // check external constructor
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(5, nullptr), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(-1, x), IGNORE_OUTPUT);

#ifdef AXOM_MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  sidre::Group* valid_group = root->createGroup("mesh");
  sidre::Group* particle_mesh = root->createGroup("particle_mesh");
  ParticleMesh(3, 10, particle_mesh);

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(nullptr, ""), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(root, ""), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(particle_mesh, ""), IGNORE_OUTPUT);

  // check 2nd push constructor
  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(nullptr, N[0]), IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(RectilinearMesh(valid_group, -1), IGNORE_OUTPUT);

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_rectilinear_mesh, native_constructor)
{
  constexpr int NDIMS = 3;
  const IndexType N[] = {5, 6, 7};
  const int64 extent[] = {0, 4, 10, 15, 7, 13};

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    RectilinearMesh* m;
    switch(idim)
    {
    case 1:
      m = new RectilinearMesh(N[0]);
      break;
    case 2:
      m = new RectilinearMesh(N[0], N[1]);
      break;
    default:
      EXPECT_EQ(idim, 3);
      m = new RectilinearMesh(N[0], N[1], N[2]);
    }

    axom::mint::internal::check_constructor(m, STRUCTURED_RECTILINEAR_MESH, idim, N);
    EXPECT_FALSE(m->isExternal());
    EXPECT_FALSE(m->hasSidreGroup());
    m->setExtent(idim, extent);
    axom::mint::internal::check_node_extent(m, extent);
    check_fill_coords(m);
    axom::mint::internal::check_create_fields(m);

    delete m;
  }  // END for all dimensions
}

//------------------------------------------------------------------------------
TEST(mint_mesh_rectilinear_mesh, external_costructor)
{
  constexpr int NDIMS = 3;
  const IndexType N[] = {5, 6, 7};
  const int64 extent[] = {0, 4, 10, 15, 7, 13};
  constexpr double MAGIC_VAL = 42.0;

  double* X = new double[N[0]];
  double* Y = new double[N[1]];
  double* Z = new double[N[2]];
  exponential_distribution(MAGIC_VAL, N[0], X);
  exponential_distribution(MAGIC_VAL, N[1], Y);
  exponential_distribution(MAGIC_VAL, N[2], Z);

  double* x = new double[N[0]];
  double* y = new double[N[1]];
  double* z = new double[N[2]];
  exponential_distribution(MAGIC_VAL, N[0], x);
  exponential_distribution(MAGIC_VAL, N[1], y);
  exponential_distribution(MAGIC_VAL, N[2], z);

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    RectilinearMesh* m = nullptr;

    switch(idim)
    {
    case 1:
      m = new RectilinearMesh(N[0], x);
      EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
      check_coordinate(m->getCoordinateArray(X_COORDINATE), x, N[0]);
      break;
    case 2:
      m = new RectilinearMesh(N[0], x, N[1], y);
      EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
      EXPECT_EQ(m->getCoordinateArray(Y_COORDINATE), y);
      check_coordinate(m->getCoordinateArray(X_COORDINATE), x, N[0]);
      check_coordinate(m->getCoordinateArray(Y_COORDINATE), y, N[1]);
      break;
    default:
      m = new RectilinearMesh(N[0], x, N[1], y, N[2], z);
      EXPECT_EQ(m->getCoordinateArray(X_COORDINATE), x);
      EXPECT_EQ(m->getCoordinateArray(Y_COORDINATE), y);
      EXPECT_EQ(m->getCoordinateArray(Z_COORDINATE), z);
      check_coordinate(m->getCoordinateArray(X_COORDINATE), x, N[0]);
      check_coordinate(m->getCoordinateArray(Y_COORDINATE), y, N[1]);
      check_coordinate(m->getCoordinateArray(Z_COORDINATE), z, N[2]);

    }  // END switch

    EXPECT_TRUE(m != nullptr);
    axom::mint::internal::check_constructor(m, STRUCTURED_RECTILINEAR_MESH, idim, N);
    EXPECT_FALSE(m->hasSidreGroup());
    EXPECT_TRUE(m->isExternal());
    m->setExtent(idim, extent);
    axom::mint::internal::check_node_extent(m, extent);

    // deallocate
    delete m;
    m = nullptr;

    // ensure coordinates are not changed after mesh gets deleted
    EXPECT_TRUE(x != nullptr);
    EXPECT_TRUE(y != nullptr);
    EXPECT_TRUE(z != nullptr);
    check_coordinate(x, X, N[0]);
    check_coordinate(y, Y, N[1]);
    check_coordinate(z, Z, N[2]);

  }  // END for all dimensions

  delete[] x;
  delete[] X;
  delete[] y;
  delete[] Y;
  delete[] z;
  delete[] Z;
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE

TEST(mint_mesh_rectilinear_mesh, sidre_constructor)
{
  constexpr int NDIMS = 3;
  const IndexType N[] = {5, 6, 7};
  const IndexType maxDim = 7;
  const int64 extent[] = {0, 4, 10, 15, 7, 13};
  constexpr double MAGIC_VAL = 42.0;

  double* expected_coords = new double[maxDim];
  exponential_distribution(MAGIC_VAL, maxDim, expected_coords);

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    // STEP 0: create a data-store with two groups
    sidre::DataStore ds;
    sidre::Group* root = ds.getRoot();
    sidre::Group* meshGroup = root->createGroup("mesh");

    // STEP 1: populate meshes in Sidre using the 2 sidre constructors
    RectilinearMesh* m;
    switch(idim)
    {
    case 1:
      m = new RectilinearMesh(meshGroup, N[I_DIRECTION]);
      break;
    case 2:
      m = new RectilinearMesh(meshGroup, N[I_DIRECTION], N[J_DIRECTION]);
      break;
    default:
      EXPECT_EQ(idim, 3);
      m = new RectilinearMesh(meshGroup, N[I_DIRECTION], N[J_DIRECTION], N[K_DIRECTION]);
    }  // END switch

    axom::mint::internal::check_constructor(m, STRUCTURED_RECTILINEAR_MESH, idim, N);
    EXPECT_TRUE(m->hasSidreGroup());
    EXPECT_FALSE(m->isExternal());
    axom::mint::internal::check_constructor(m, STRUCTURED_RECTILINEAR_MESH, idim, N);
    m->setExtent(idim, extent);
    axom::mint::internal::check_node_extent(m, extent);
    check_fill_coords(m);
    axom::mint::internal::check_create_fields(m);

    delete m;
    m = nullptr;

    // STEP 2: pull the mesh from sidre and check correctness
    m = new RectilinearMesh(meshGroup);
    axom::mint::internal::check_constructor(m, STRUCTURED_RECTILINEAR_MESH, idim, N);
    EXPECT_TRUE(m->hasSidreGroup());
    EXPECT_FALSE(m->isExternal());

    for(int ii = 0; ii < idim; ++ii)
    {
      check_coordinate(m->getCoordinateArray(ii), expected_coords, m->getNodeResolution(ii));
    }

    axom::mint::internal::check_fields(m, true);
    delete m;

    // STEP 3: ensure data is persistent in Sidre
    EXPECT_TRUE(blueprint::isValidRootGroup(meshGroup));

  }  // END for all dimensions

  delete[] expected_coords;
}

#endif /* ENDIF AXOM_MINT_USE_SIDRE */

//------------------------------------------------------------------------------
TEST(mint_mesh_rectilinear_mesh, get_node)
{
  // initialize coordinates
  const IndexType N[] = {5, 5, 5};
  double x[5];
  double y[5];
  double z[5];
  exponential_distribution(42.0, 5, x);
  exponential_distribution(0.0, 5, y);
  exponential_distribution(0.0, 5, z);

  RectilinearMesh m(N[0], x, N[1], y, N[2], z);
  axom::mint::internal::check_constructor(&m, STRUCTURED_RECTILINEAR_MESH, 3, N);

  const IndexType kp = m.nodeKp();
  const IndexType jp = m.nodeJp();
  const double* xx = m.getCoordinateArray(X_COORDINATE);
  const double* yy = m.getCoordinateArray(Y_COORDINATE);
  const double* zz = m.getCoordinateArray(Z_COORDINATE);
  EXPECT_EQ(x, xx);
  EXPECT_EQ(y, yy);
  EXPECT_EQ(z, zz);

  for(IndexType k = 0; k < N[K_DIRECTION]; ++k)
  {
    const IndexType k_offset = k * kp;
    for(IndexType j = 0; j < N[J_DIRECTION]; ++j)
    {
      const IndexType j_offset = j * jp;
      for(IndexType i = 0; i < N[I_DIRECTION]; ++i)
      {
        const IndexType nodeID = i + j_offset + k_offset;

        double node[3];
        m.getNode(nodeID, node);

        EXPECT_DOUBLE_EQ(xx[i], node[X_COORDINATE]);
        EXPECT_DOUBLE_EQ(yy[j], node[Y_COORDINATE]);
        EXPECT_DOUBLE_EQ(zz[k], node[Z_COORDINATE]);

      }  // END for all i
    }  // END for all j
  }  // END for all k
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();
  return result;
}
