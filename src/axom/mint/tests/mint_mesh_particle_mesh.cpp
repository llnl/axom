// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/config.hpp"  // for compile-time type definitions

#include "axom/mint/mesh/blueprint.hpp"     // for mint::blueprint() functions
#include "axom/mint/mesh/ParticleMesh.hpp"  // for ParticleMesh definition

#include "axom/slic.hpp"

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"

// namespace aliases
namespace mint = axom::mint;

namespace
{
// globals
const char* IGNORE_OUTPUT = ".*";

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
void check_constructor(mint::ParticleMesh* particles,
                       int expected_dimension,
                       axom::IndexType expected_num_particles)
{
  EXPECT_TRUE(particles != nullptr);
  EXPECT_EQ(particles->getDimension(), expected_dimension);
  EXPECT_EQ(particles->getNumberOfNodes(), expected_num_particles);
  EXPECT_EQ(particles->getNumberOfCells(), expected_num_particles);
  EXPECT_TRUE(particles->getNumberOfNodes() <= particles->getNodeCapacity());
  EXPECT_EQ(particles->getMeshType(), mint::PARTICLE_MESH);
  EXPECT_TRUE(particles->hasExplicitCoordinates());
  EXPECT_FALSE(particles->hasExplicitConnectivity());
  EXPECT_FALSE(particles->hasMixedCellTypes());
  EXPECT_EQ(particles->getCellType(), mint::VERTEX);

  axom::IndexType ncells = particles->getNumberOfCells();
  axom::IndexType cell[1];
  for(axom::IndexType icell = 0; icell < ncells; ++icell)
  {
    EXPECT_EQ(particles->getCellType(icell), mint::VERTEX);

    particles->getCellNodeIDs(icell, cell);
    EXPECT_EQ(cell[0], icell);
  }

  const int ndims = particles->getDimension();
  for(int idim = 0; idim < ndims; ++idim)
  {
    const axom::IndexType numParticles = particles->getNumberOfNodes();
    const axom::IndexType lastParticle = numParticles - 1;
    double* pos = particles->getCoordinateArray(idim);
    EXPECT_TRUE(pos != nullptr);

    pos[0] = pos[lastParticle] = 42.0;
    EXPECT_DOUBLE_EQ(pos[0], 42.0);
    EXPECT_DOUBLE_EQ(pos[lastParticle], 42.0);
  }
}

//------------------------------------------------------------------------------
void check_resize(mint::ParticleMesh* particles)
{
  EXPECT_TRUE(particles != nullptr);
  constexpr axom::IndexType NEW_SIZE = 512;

  particles->resize(NEW_SIZE);
  EXPECT_EQ(particles->getNumberOfNodes(), NEW_SIZE);
  EXPECT_TRUE(particles->getNumberOfNodes() <= particles->getNodeCapacity());

  // ensure all fiels are resized as well
  const mint::FieldData* fd = particles->getFieldData(mint::NODE_CENTERED);
  EXPECT_TRUE(fd->getNumFields() > 0);
  for(int ifield = 0; ifield < fd->getNumFields(); ++ifield)
  {
    const mint::Field* f = fd->getField(ifield);
    EXPECT_EQ(f->getNumTuples(), particles->getNumberOfNodes());
  }
}

//------------------------------------------------------------------------------
void check_reserve(mint::ParticleMesh* particles)
{
  EXPECT_TRUE(particles != nullptr);
  constexpr axom::IndexType NEW_CAPACITY = 512;

  particles->reserve(NEW_CAPACITY);
  EXPECT_EQ(particles->getNodeCapacity(), NEW_CAPACITY);
  EXPECT_TRUE(particles->getNumberOfNodes() <= particles->getNodeCapacity());

  // ensure all fiels are resized as well
  const mint::FieldData* fd = particles->getFieldData(mint::NODE_CENTERED);
  EXPECT_TRUE(fd->getNumFields() > 0);
  for(int ifield = 0; ifield < fd->getNumFields(); ++ifield)
  {
    const mint::Field* f = fd->getField(ifield);
    EXPECT_EQ(f->getCapacity(), particles->getNodeCapacity());
  }
}

//------------------------------------------------------------------------------
void check_shrink(mint::ParticleMesh* particles, axom::IndexType NUM_PARTICLES, axom::IndexType CAPACITY)
{
  EXPECT_TRUE(particles != nullptr);
  EXPECT_TRUE(NUM_PARTICLES > 0);
  EXPECT_EQ(particles->getNumberOfNodes(), NUM_PARTICLES);
  EXPECT_EQ(particles->getNodeCapacity(), CAPACITY);

  // ensure all fields have the specified capacity
  const mint::FieldData* fd = particles->getFieldData(mint::NODE_CENTERED);
  EXPECT_TRUE(fd->getNumFields() > 0);
  for(int ifield = 0; ifield < fd->getNumFields(); ++ifield)
  {
    const mint::Field* f = fd->getField(ifield);
    EXPECT_EQ(f->getCapacity(), particles->getNodeCapacity());
    EXPECT_EQ(f->getNumTuples(), particles->getNumberOfNodes());
  }

  // call shrink
  particles->shrink();

  // check particles
  EXPECT_EQ(particles->getNodeCapacity(), particles->getNumberOfNodes());

  for(int ifield = 0; ifield < fd->getNumFields(); ++ifield)
  {
    const mint::Field* f = fd->getField(ifield);
    EXPECT_EQ(f->getCapacity(), f->getNumTuples());
    EXPECT_EQ(f->getNumTuples(), particles->getNumberOfNodes());
  }
}

//------------------------------------------------------------------------------
void check_append(mint::ParticleMesh* particles)
{
  EXPECT_TRUE(particles != nullptr);

  const mint::FieldData* fd = particles->getFieldData(mint::NODE_CENTERED);
  EXPECT_TRUE(fd->getNumFields() > 0);

  constexpr int NUM_APPENDS = 3;
  constexpr double MAGIC_NUMBER = 42;

  const double* x = nullptr;
  const double* y = nullptr;
  const double* z = nullptr;
  axom::IndexType lidx = -1;

  const int dimension = particles->getDimension();

  for(int iter = 0; iter < NUM_APPENDS; ++iter)
  {
    axom::IndexType current_num_particles = particles->getNumberOfNodes();

    switch(dimension)
    {
    case 1:
      particles->append(MAGIC_NUMBER);
      x = particles->getCoordinateArray(mint::X_COORDINATE);
      lidx = particles->getNumberOfNodes() - 1;
      EXPECT_EQ(x[lidx], MAGIC_NUMBER);
      break;
    case 2:
      particles->append(MAGIC_NUMBER, MAGIC_NUMBER);
      lidx = particles->getNumberOfNodes() - 1;
      x = particles->getCoordinateArray(mint::X_COORDINATE);
      y = particles->getCoordinateArray(mint::Y_COORDINATE);
      EXPECT_EQ(x[lidx], MAGIC_NUMBER);
      EXPECT_EQ(y[lidx], MAGIC_NUMBER);
      break;
    default:
      EXPECT_TRUE(dimension == 3);
      particles->append(MAGIC_NUMBER, MAGIC_NUMBER, MAGIC_NUMBER);
      lidx = particles->getNumberOfNodes() - 1;
      x = particles->getCoordinateArray(mint::X_COORDINATE);
      y = particles->getCoordinateArray(mint::Y_COORDINATE);
      z = particles->getCoordinateArray(mint::Z_COORDINATE);
      EXPECT_EQ(x[lidx], MAGIC_NUMBER);
      EXPECT_EQ(y[lidx], MAGIC_NUMBER);
      EXPECT_EQ(z[lidx], MAGIC_NUMBER);
    }  // END switch

    EXPECT_EQ(particles->getNumberOfNodes(), current_num_particles + 1);

    // ensure the fields are also resized accordingly when a new particle
    // is appended
    for(int ifield = 0; ifield < fd->getNumFields(); ++ifield)
    {
      const mint::Field* f = fd->getField(ifield);
      EXPECT_TRUE(f != nullptr);
      EXPECT_EQ(f->getNumTuples(), particles->getNumberOfNodes());
    }

    // check invariant
    EXPECT_TRUE(particles->getNumberOfNodes() <= particles->getNodeCapacity());

    // shrink so that next append triggers a realloc
    particles->shrink();
    EXPECT_EQ(particles->getNumberOfNodes(), particles->getNodeCapacity());
  }
}
//------------------------------------------------------------------------------
void check_create_field(mint::ParticleMesh* particles, const std::string& name, int numComponents)
{
  EXPECT_TRUE(particles != nullptr);

  const mint::Field* f = nullptr;
  const int assoc = mint::NODE_CENTERED;

  double* vel = particles->createField<double>(name, assoc, numComponents);
  EXPECT_TRUE(vel != nullptr);
  EXPECT_TRUE(particles->hasField(name, assoc));

  f = particles->getFieldData(mint::NODE_CENTERED)->getField(name);
  EXPECT_EQ(particles->getNumberOfNodes(), f->getNumTuples());
  EXPECT_EQ(f->getNumComponents(), numComponents);
  EXPECT_EQ(vel, mint::Field::getDataPtr<double>(f));
}

}  // namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh_DeathTest, invalid_construction)
{
  const int INVALID_DIMENSION = -1;
  const int INVALID_NUM_PARTICLES = -10;
  EXPECT_DEATH_IF_SUPPORTED(mint::ParticleMesh(INVALID_DIMENSION, 10), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(mint::ParticleMesh(3, INVALID_NUM_PARTICLES), IGNORE_OUTPUT);

  // capacity cannot be smaller than specified num particles
  EXPECT_DEATH_IF_SUPPORTED(mint::ParticleMesh(3, 10, 5), IGNORE_OUTPUT);

#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  // construct from an empty group should fail!
  EXPECT_DEATH_IF_SUPPORTED(mint::ParticleMesh(root, ""), IGNORE_OUTPUT);

  // construct from a group that is not blueprint-conforming also fails!
  root->createGroup("foo")->createView("bar");
  EXPECT_DEATH_IF_SUPPORTED(mint::ParticleMesh(root, ""), IGNORE_OUTPUT);

#endif
}

//-----------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh_DeathTest, invalid_operations)
{
  const axom::IndexType numParticles = 10;
  mint::ParticleMesh particles(2, numParticles);
  EXPECT_DEATH_IF_SUPPORTED(particles.getCoordinateArray(mint::Z_COORDINATE), IGNORE_OUTPUT);

  // creating/accessing anything other than node-centered fields on a particle
  // mesh does not make sense and should fail
  EXPECT_DEATH_IF_SUPPORTED(particles.getFieldData(mint::CELL_CENTERED), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(particles.createField<double>("foobar", mint::CELL_CENTERED),
                            IGNORE_OUTPUT);

  // append a particle to a 2-D ParticleMesh using a 1-D append()
  EXPECT_DEATH_IF_SUPPORTED(particles.append(42.0), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(particles.append(1, 2, 3), IGNORE_OUTPUT);

  double x[3] = {1.0, 2.0, 3.0};
  mint::ParticleMesh particles_external(3, x);
  EXPECT_DEATH_IF_SUPPORTED(particles_external.append(2), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(particles_external.resize(10), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(particles_external.reserve(20), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh, native_constructor)
{
  const axom::IndexType numParticles = 10;
  for(int dim = 1; dim <= 3; ++dim)
  {
    mint::ParticleMesh particles(dim, numParticles);
    check_constructor(&particles, dim, numParticles);
    check_create_field(&particles, "foo", 3);
    check_create_field(&particles, "bar", 1);
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh, external_constructor)
{
  const axom::IndexType numParticles = 4;
  double x[] = {1.0, 2.0, 3.0, 4.0};
  double y[] = {1.0, 2.0, 3.0, 4.0};
  double z[] = {1.0, 2.0, 3.0, 4.0};

  // BEGIN SCOPE
  {
    mint::ParticleMesh particles1d(numParticles, x);
    check_constructor(&particles1d, 1, numParticles);
    check_create_field(&particles1d, "foobar", 4);
    EXPECT_EQ(particles1d.getCoordinateArray(mint::X_COORDINATE), x);

    mint::ParticleMesh particles2d(numParticles, x, y);
    check_constructor(&particles2d, 2, numParticles);
    check_create_field(&particles2d, "foobar", 5);
    EXPECT_EQ(particles2d.getCoordinateArray(mint::X_COORDINATE), x);
    EXPECT_EQ(particles2d.getCoordinateArray(mint::Y_COORDINATE), y);

    mint::ParticleMesh particles3d(numParticles, x, y, z);
    check_constructor(&particles3d, 3, numParticles);
    check_create_field(&particles3d, "foobar", 1);
    EXPECT_EQ(particles3d.getCoordinateArray(mint::X_COORDINATE), x);
    EXPECT_EQ(particles3d.getCoordinateArray(mint::Y_COORDINATE), y);
    EXPECT_EQ(particles3d.getCoordinateArray(mint::Z_COORDINATE), z);
  }
  // END SCOPE

  for(int i = 0; i < numParticles; ++i)
  {
    const double expected_val =
      ((i == 0) || (i == numParticles - 1)) ? 42.0 : static_cast<double>(i + 1);

    EXPECT_DOUBLE_EQ(x[i], expected_val);
    EXPECT_DOUBLE_EQ(y[i], expected_val);
    EXPECT_DOUBLE_EQ(z[i], expected_val);
  }
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh, sidre_constructor)
{
  const axom::IndexType numParticles = 4;
  double x[] = {1.0, 2.0, 3.0, 4.0};
  double y[] = {1.0, 2.0, 3.0, 4.0};
  double z[] = {1.0, 2.0, 3.0, 4.0};
  double* data[3] = {x, y, z};
  const double MAGIC_NUMBER = 42.0;
  const int BLOCK_ID = 9;
  const int PART_ID = 10;

  for(int dim = 1; dim <= 3; ++dim)
  {
    // create a datastore
    sidre::DataStore ds;
    sidre::Group* root = ds.getRoot();

    // BEGIN SCOPE
    {
      // create a particle mesh on sidre
      mint::ParticleMesh particles(dim, numParticles, root);
      particles.setBlockId(BLOCK_ID);
      particles.setPartitionId(PART_ID);
      check_constructor(&particles, dim, numParticles);
      check_create_field(&particles, "foo", 3);
      EXPECT_TRUE(particles.hasSidreGroup());
      EXPECT_EQ(particles.getBlockId(), BLOCK_ID);
      EXPECT_EQ(particles.getPartitionId(), PART_ID);

      for(int idim = 0; idim < dim; ++idim)
      {
        double* pos = particles.getCoordinateArray(idim);
        memcpy(pos, data[idim], numParticles * sizeof(double));
      }

      double* foo = particles.getFieldPtr<double>("foo", mint::NODE_CENTERED);
      EXPECT_TRUE(foo != nullptr);

      for(axom::IndexType ipart = 0; ipart < numParticles; ++ipart)
      {
        foo[ipart * 3] = MAGIC_NUMBER;
        foo[ipart * 3 + 1] = MAGIC_NUMBER;
        foo[ipart * 3 + 2] = MAGIC_NUMBER;
      }
    }
    // END SCOPE

    // BEGIN SCOPE
    {
      // create a particle mesh from the sidre group (ensure data is persistent)
      mint::ParticleMesh particles(root);
      check_constructor(&particles, dim, numParticles);
      EXPECT_TRUE(particles.hasSidreGroup());
      EXPECT_TRUE(particles.hasField("foo", mint::NODE_CENTERED));
      EXPECT_EQ(particles.getBlockId(), BLOCK_ID);
      EXPECT_EQ(particles.getPartitionId(), PART_ID);

      axom::IndexType numComp = -1;
      const double* foo = particles.getFieldPtr<double>("foo", mint::NODE_CENTERED, numComp);
      EXPECT_TRUE(foo != nullptr);
      EXPECT_EQ(numComp, 3);

      for(axom::IndexType i = 0; i < numParticles; ++i)
      {
        EXPECT_DOUBLE_EQ(foo[i * 3], MAGIC_NUMBER);
        EXPECT_DOUBLE_EQ(foo[i * 3 + 1], MAGIC_NUMBER);
        EXPECT_DOUBLE_EQ(foo[i * 3 + 2], MAGIC_NUMBER);
      }

      for(int idim = 0; idim < dim; ++idim)
      {
        double* pos = particles.getCoordinateArray(idim);
        EXPECT_TRUE(pos != nullptr);

        for(axom::IndexType i = 0; i < numParticles; ++i)
        {
          const double expected_val =
            ((i == 0) || (i == numParticles - 1)) ? 42.0 : static_cast<double>(i + 1);
          EXPECT_DOUBLE_EQ(pos[i], expected_val);
        }
      }
    }
    // END SCOPE

    // ensure data is persistent in sidre
    EXPECT_TRUE(mint::blueprint::isValidRootGroup(root));

    sidre::Group* coordsets = root->getGroup("coordsets");
    sidre::Group* topologies = root->getGroup("topologies");
    sidre::Group* fields = root->getGroup("fields");
    EXPECT_TRUE(coordsets != nullptr);
    EXPECT_TRUE(topologies != nullptr);
    EXPECT_TRUE(fields != nullptr);
    EXPECT_EQ(coordsets->getNumGroups(), 1);
    EXPECT_EQ(topologies->getNumGroups(), 1);
    EXPECT_EQ(fields->getNumGroups(), 1);
    EXPECT_TRUE(fields->hasChildGroup("foo"));
  }  // END for all dimensions
}

#endif

//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh, append)
{
  constexpr int NDIMS = 3;
  constexpr axom::IndexType NUM_PARTICLES = 10;

  for(int idim = 1; idim <= NDIMS; ++idim)
  {
    mint::ParticleMesh p1(idim, NUM_PARTICLES);
    p1.createField<double>("vel", mint::NODE_CENTERED, 3);
    p1.createField<int>("id", mint::NODE_CENTERED);
    check_append(&p1);

#ifdef AXOM_MINT_USE_SIDRE

    // create a datastore
    sidre::DataStore ds;
    sidre::Group* root = ds.getRoot();

    mint::ParticleMesh p2(idim, 0, root);
    p2.createField<double>("vel", mint::NODE_CENTERED, 3);
    p2.createField<int>("id", mint::NODE_CENTERED);
    check_append(&p2);

#endif
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh, resize)
{
  constexpr int NDIMS = 3;
  mint::ParticleMesh p1(NDIMS, 10);
  p1.createField<double>("vel", mint::NODE_CENTERED, 3);
  p1.createField<int>("id", mint::NODE_CENTERED);
  check_resize(&p1);

#ifdef AXOM_MINT_USE_SIDRE

  // create a datastore
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  mint::ParticleMesh p2(NDIMS, 0, root);
  p2.createField<double>("vel", mint::NODE_CENTERED, 3);
  p2.createField<int>("id", mint::NODE_CENTERED);
  check_resize(&p2);

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh, reserve)
{
  constexpr int NDIMS = 3;
  mint::ParticleMesh p1(NDIMS, 10);
  p1.createField<double>("vel", mint::NODE_CENTERED, 3);
  p1.createField<int>("id", mint::NODE_CENTERED);
  check_reserve(&p1);

#ifdef AXOM_MINT_USE_SIDRE

  // create a datastore
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  mint::ParticleMesh p2(NDIMS, 0, root);
  p2.createField<double>("vel", mint::NODE_CENTERED, 3);
  p2.createField<int>("id", mint::NODE_CENTERED);
  check_reserve(&p2);

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_particle_mesh, shrink)
{
  constexpr int NDIMS = 3;
  constexpr axom::IndexType NUM_PARTICLES = 10;
  constexpr axom::IndexType CAPACITY = 512;
  mint::ParticleMesh p1(NDIMS, NUM_PARTICLES, CAPACITY);
  p1.createField<double>("vel", mint::NODE_CENTERED, 3, true);
  p1.createField<int>("id", mint::NODE_CENTERED, 1, true);
  check_shrink(&p1, NUM_PARTICLES, CAPACITY);

#ifdef AXOM_MINT_USE_SIDRE

  // create a datastore
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  mint::ParticleMesh p2(NDIMS, NUM_PARTICLES, root, CAPACITY);
  p2.createField<double>("vel", mint::NODE_CENTERED, 3, true);
  p2.createField<int>("id", mint::NODE_CENTERED, 1, true);
  check_shrink(&p2, NUM_PARTICLES, CAPACITY);

#endif
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
