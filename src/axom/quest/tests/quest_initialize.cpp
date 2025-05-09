// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/mint.hpp"
#include "quest_test_utilities.hpp"
#if defined AXOM_USE_SIDRE
  #include "axom/sidre.hpp"
#endif

#include "axom/quest/interface/inout.hpp"
#include "axom/quest/interface/signed_distance.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

// Google test include
#include "gtest/gtest.h"

namespace mint = axom::mint;
namespace quest = axom::quest;

// Test initializing quest inout from a preloaded mesh
TEST(quest_initialize, inout_pointer_initialize)
{
  int rc = quest::QUEST_INOUT_SUCCESS;

  std::shared_ptr<mint::Mesh> input_mesh {axom::quest::utilities::make_tetrahedron_mesh()};

  // Note: the following call updates the input_mesh pointer
#ifdef AXOM_USE_MPI
  rc = quest::inout_init(input_mesh, MPI_COMM_WORLD);
#else
  rc = quest::inout_init(input_mesh);
#endif
  EXPECT_EQ(quest::QUEST_INOUT_SUCCESS, rc);

  EXPECT_TRUE(quest::inout_initialized());

  EXPECT_TRUE(quest::inout_evaluate(3, 2, 0));
  EXPECT_TRUE(quest::inout_evaluate(-1, 2, -1));
  EXPECT_FALSE(quest::inout_evaluate(4, 4, -7));

  rc = quest::inout_finalize();
  EXPECT_EQ(quest::QUEST_INOUT_SUCCESS, rc);
}

// Test initializing quest signed_distance from a preloaded mesh
TEST(quest_initialize, signed_distance_pointer_initialize)
{
  int rc = 0;

  mint::Mesh* input_mesh = axom::quest::utilities::make_tetrahedron_mesh();

#ifdef AXOM_USE_MPI
  rc = quest::signed_distance_init(input_mesh, MPI_COMM_WORLD);
#else
  rc = quest::signed_distance_init(input_mesh);
#endif
  EXPECT_EQ(0, rc);

  EXPECT_TRUE(quest::signed_distance_initialized());

  EXPECT_GT(0., quest::signed_distance_evaluate(3, 2, 0));
  EXPECT_GT(0., quest::signed_distance_evaluate(-1, 2, -1));
  EXPECT_LT(0., quest::signed_distance_evaluate(4, 4, -7));

  quest::signed_distance_finalize();

  delete input_mesh;
}

/*
  The *_reallocations tests and the immediate_ug_reserve are
  reproducers for the case of multiple Umpire instances that happens
  when Umpire is built with static lib and Axom with shared libs.
  The failures only appear when there are several shared libs that
  include the static Umpire, e.g. axom::core, axom::sidre, axom::mint
  and axom::quest.  This issue is specific to the linker commands in
  quest.
*/
// Test allocation/reallocation using axom::allocate and axom::reallocate bytes
TEST(quest_initialize, byte_reallocations)
{
  std::int8_t* b1 = axom::allocate<std::int8_t>(40);
  EXPECT_NE(b1, nullptr);

  b1 = axom::reallocate<std::int8_t>(b1, 80);
  EXPECT_NE(b1, nullptr);

  axom::deallocate(b1);
}

// Test allocation/reallocation using axom::allocate and axom::reallocate ints
TEST(quest_initialize, int_reallocations)
{
  int* b3 = axom::allocate<int>(10);
  EXPECT_NE(b3, nullptr);

  b3 = axom::reallocate<int>(b3, 20);
  EXPECT_NE(b3, nullptr);

  axom::deallocate(b3);
}

#if defined(AXOM_USE_SIDRE)
// Test allocation/reallocation using sidre::Buffer.
TEST(quest_initialize, buffer_reallocations)
{
  axom::sidre::DataStore dataStore;

  axom::sidre::Buffer* b1 = dataStore.createBuffer(axom::sidre::DataTypeId::INT32_ID, 10);
  b1->allocate();
  EXPECT_NE(b1->getVoidPtr(), nullptr);

  b1->reallocate(20);
  EXPECT_NE(b1->getVoidPtr(), nullptr);

  b1->deallocate();
}

// Test allocation/reallocation using sidre::View.
TEST(quest_initialize, view_reallocations)
{
  axom::sidre::DataStore dataStore;
  axom::sidre::Group* group = dataStore.getRoot();

  axom::sidre::View* v3 = group->createView("v3", conduit::DataType::int32(10));
  v3->allocate();
  EXPECT_NE(v3->getVoidPtr(), nullptr);

  v3->reallocate(20);
  EXPECT_NE(v3->getVoidPtr(), nullptr);

  v3->deallocate();
}

  #ifdef AXOM_MINT_USE_SIDRE
// Test immediately reserving space in UnstructuredMesh.
TEST(quest_initialize, immediate_ug_reserve)
{
  axom::sidre::DataStore dataStore;
  axom::sidre::Group* meshGroup = dataStore.getRoot()->createGroup("myGroup");
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> contourMesh(2,
                                                                     axom::mint::CellType::SEGMENT,
                                                                     meshGroup);
  contourMesh.reserveCells(10);  // This may unexpectedly crash.
}
  #endif
#endif

int main(int argc, char** argv)
{
#ifdef AXOM_USE_MPI
  // Initialize MPI
  MPI_Init(&argc, &argv);
#endif

  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return result;
}
