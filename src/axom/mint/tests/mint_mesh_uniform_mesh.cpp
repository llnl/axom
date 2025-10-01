// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/config.hpp"             // for compile-time type definitions
#include "axom/mint/mesh/blueprint.hpp"     // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"     // for CellTypes enum definition
#include "axom/mint/mesh/ParticleMesh.hpp"  // for ParticleMesh
#include "axom/mint/mesh/UniformMesh.hpp"   // for UniformMesh
#include "StructuredMesh_helpers.hpp"       // for StructuredMesh test helpers

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
#ifdef AXOM_MINT_USE_SIDRE

void check_sidre_group(sidre::Group* root_group,
                       int expected_dimension,
                       const double* expected_origin,
                       const double* expected_spacing)
{
  EXPECT_TRUE(blueprint::isValidRootGroup(root_group));

  const sidre::Group* topology = blueprint::getTopologyGroup(root_group);
  EXPECT_TRUE(blueprint::isValidTopologyGroup(topology));
  EXPECT_EQ(topology->getParent()->getNumGroups(), 1);
  EXPECT_EQ(topology->getParent()->getNumViews(), 0);

  const sidre::Group* coordset = blueprint::getCoordsetGroup(root_group, topology);
  EXPECT_TRUE(blueprint::isValidCoordsetGroup(coordset));
  EXPECT_EQ(coordset->getParent()->getNumGroups(), 1);
  EXPECT_EQ(coordset->getParent()->getNumViews(), 0);

  int mesh_type = UNDEFINED_MESH;
  int dimension = -1;
  blueprint::getMeshTypeAndDimension(mesh_type, dimension, root_group);
  EXPECT_EQ(mesh_type, STRUCTURED_UNIFORM_MESH);
  EXPECT_EQ(dimension, expected_dimension);

  double mesh_origin[3];
  double mesh_spacing[3];
  blueprint::getUniformMeshProperties(dimension, mesh_origin, mesh_spacing, coordset);

  for(int i = 0; i < dimension; ++i)
  {
    EXPECT_DOUBLE_EQ(expected_origin[i], mesh_origin[i]);
    EXPECT_DOUBLE_EQ(expected_spacing[i], mesh_spacing[i]);
  }
}

#endif

}  // namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_uniform_mesh_DeathTest, invalid_construction)
{
  const double lo[] = {0.0, 0.0, 0.0};
  const double hi[] = {2.0, 2.0, 2.0};
  const IndexType Ni = 5;
  const IndexType Nj = 5;
  const IndexType Nk = 5;

  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(nullptr, hi, Ni, Nj, Nk), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(lo, nullptr, Ni, Nj, Nk), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(lo, hi, -1), IGNORE_OUTPUT);

#ifdef AXOM_MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  sidre::Group* valid_group = root->createGroup("mesh");
  sidre::Group* particle_mesh = root->createGroup("particle_mesh");
  ParticleMesh(3, 10, particle_mesh);

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(root, ""), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(particle_mesh, ""), IGNORE_OUTPUT);

  // check push constructor
  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(particle_mesh, lo, hi, Ni, Nj, Nk), IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(nullptr, lo, hi, Ni, Nj, Nk), IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(valid_group, lo, nullptr, Ni, Nj, Nk), IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(valid_group, nullptr, hi, Ni, Nj, Nk), IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(UniformMesh(valid_group, lo, hi, -1), IGNORE_OUTPUT);

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_uniform_mesh_DeathTest, invalid_operations)
{
  const double low[] = {0.0, 0.0, 0.0};
  const double high[] = {0.5, 0.5, 0.5};
  const IndexType Ni = 4;
  const IndexType Nj = 4;
  const IndexType Nk = 4;

  UniformMesh m(low, high, Ni, Nj, Nk);
  EXPECT_DEATH_IF_SUPPORTED(m.getCoordinateArray(0), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(m.getCoordinateArray(1), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(m.getCoordinateArray(2), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_uniform_mesh, native_constructor)
{
  constexpr int NDIMS = 3;

  const double lo[] = {0.0, 0.0, 0.0};
  const double hi[] = {2.0, 2.0, 2.0};
  const double h[] = {0.5, 0.5, 0.5};
  const IndexType N[] = {5, 5, 5};
  const int64 node_ext[] = {-55, -51, 10, 14, 0, 4};

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    UniformMesh* m;
    switch(idim)
    {
    case 1:
      m = new UniformMesh(lo, hi, N[0]);
      break;
    case 2:
      m = new UniformMesh(lo, hi, N[0], N[1]);
      break;
    default:
      EXPECT_TRUE(idim == 3);
      m = new UniformMesh(lo, hi, N[0], N[1], N[2]);
    }

    axom::mint::internal::check_constructor(m, idim, lo, h, N);
    EXPECT_FALSE(m->hasSidreGroup());

    axom::mint::internal::check_create_fields(m);
    axom::mint::internal::check_fields(m, false);
    m->setExtent(idim, node_ext);
    axom::mint::internal::check_node_extent(m, node_ext);
    delete m;
  }  // END for all dimensions
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
TEST(mint_mesh_uniform_mesh, sidre_constructor)
{
  constexpr int NDIMS = 3;
  const double lo[] = {0.0, 0.0, 0.0};
  const double hi[] = {2.0, 2.0, 2.0};
  const double h[] = {0.5, 0.5, 0.5};
  const IndexType N[] = {5, 5, 5};
  const int64 node_ext[] = {-55, -51, 10, 14, 0, 4};

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    // STEP 0: create a datastore to store the mesh.
    sidre::DataStore ds;
    sidre::Group* root = ds.getRoot();
    sidre::Group* meshGroup = root->createGroup("mesh");

    // STEP 1: create a uniform mesh in Sidre.
    UniformMesh* m;
    switch(idim)
    {
    case 1:
      m = new UniformMesh(meshGroup, lo, hi, N[0]);
      break;
    case 2:
      m = new UniformMesh(meshGroup, lo, hi, N[0], N[1]);
      break;
    default:
      EXPECT_TRUE(idim == 3);
      m = new UniformMesh(meshGroup, lo, hi, N[0], N[1], N[2]);
    }  // END switch

    EXPECT_TRUE(m->hasSidreGroup());
    axom::mint::internal::check_constructor(m, idim, lo, h, N);
    axom::mint::internal::check_create_fields(m);
    m->setExtent(idim, node_ext);
    axom::mint::internal::check_node_extent(m, node_ext);

    delete m;
    m = nullptr;

    // STEP 2: pull the mesh from the sidre groups, and check with expected
    m = new UniformMesh(meshGroup);
    axom::mint::internal::check_constructor(m, idim, lo, h, N);
    axom::mint::internal::check_fields(m, true);
    axom::mint::internal::check_node_extent(m, node_ext);
    delete m;

    // STEP 3: ensure the data is persistent in sidre
    check_sidre_group(meshGroup, idim, lo, h);
  }
}

#endif

//------------------------------------------------------------------------------
TEST(mint_mesh_uniform_mesh, check_evaluate_coordinate)
{
  constexpr int NDIMS = 3;
  const double lo[] = {0.0, 0.0, 0.0};
  const double hi[] = {2.0, 2.0, 2.0};
  const double h[] = {0.5, 0.5, 0.5};
  const IndexType N[] = {5, 5, 5};

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    switch(idim)
    {
    case 1:
    {
      UniformMesh m(lo, hi, N[0]);
      axom::mint::internal::check_constructor(&m, idim, lo, h, N);
      EXPECT_FALSE(m.hasSidreGroup());

      const IndexType Ni = m.getNodeResolution(I_DIRECTION);
      EXPECT_EQ(Ni, N[I_DIRECTION]);

      for(IndexType i = 0; i < Ni; ++i)
      {
        const double expected_x = lo[X_COORDINATE] + i * h[X_COORDINATE];
        const double x = m.evaluateCoordinate(i, I_DIRECTION);
        EXPECT_DOUBLE_EQ(expected_x, x);
      }  // END for all i

    }  // END 1D
    break;
    case 2:
    {
      UniformMesh m(lo, hi, N[0], N[1]);
      axom::mint::internal::check_constructor(&m, idim, lo, h, N);
      EXPECT_FALSE(m.hasSidreGroup());

      const IndexType Ni = m.getNodeResolution(I_DIRECTION);
      const IndexType Nj = m.getNodeResolution(J_DIRECTION);
      EXPECT_EQ(Ni, N[I_DIRECTION]);
      EXPECT_EQ(Nj, N[J_DIRECTION]);

      for(IndexType j = 0; j < Nj; ++j)
      {
        for(IndexType i = 0; i < Ni; ++i)
        {
          const double expected_x = lo[X_COORDINATE] + i * h[X_COORDINATE];
          const double expected_y = lo[Y_COORDINATE] + j * h[Y_COORDINATE];

          const double x = m.evaluateCoordinate(i, I_DIRECTION);
          const double y = m.evaluateCoordinate(j, J_DIRECTION);

          EXPECT_DOUBLE_EQ(expected_x, x);
          EXPECT_DOUBLE_EQ(expected_y, y);
        }  // END for all i
      }  // END for all j

    }  // END 2D
    break;
    default:
      EXPECT_EQ(idim, 3);
      {
        UniformMesh m(lo, hi, N[0], N[1], N[2]);
        axom::mint::internal::check_constructor(&m, idim, lo, h, N);
        EXPECT_FALSE(m.hasSidreGroup());

        const IndexType Ni = m.getNodeResolution(I_DIRECTION);
        const IndexType Nj = m.getNodeResolution(J_DIRECTION);
        const IndexType Nk = m.getNodeResolution(K_DIRECTION);
        EXPECT_EQ(Ni, N[I_DIRECTION]);
        EXPECT_EQ(Nj, N[J_DIRECTION]);
        EXPECT_EQ(Nk, N[K_DIRECTION]);

        for(IndexType k = 0; k < Nk; ++k)
        {
          for(IndexType j = 0; j < Nj; ++j)
          {
            for(IndexType i = 0; i < Ni; ++i)
            {
              const double expected_x = lo[X_COORDINATE] + i * h[X_COORDINATE];
              const double expected_y = lo[Y_COORDINATE] + j * h[Y_COORDINATE];
              const double expected_z = lo[Z_COORDINATE] + k * h[Z_COORDINATE];

              const double x = m.evaluateCoordinate(i, I_DIRECTION);
              const double y = m.evaluateCoordinate(j, J_DIRECTION);
              const double z = m.evaluateCoordinate(k, K_DIRECTION);

              EXPECT_DOUBLE_EQ(expected_x, x);
              EXPECT_DOUBLE_EQ(expected_y, y);
              EXPECT_DOUBLE_EQ(expected_z, z);
            }  // END for all i
          }  // END for all j
        }  // END for all k

      }  // END 3D

    }  // END switch
  }
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
