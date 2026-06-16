// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \brief Blueprint-focused unit tests for quest's SamplingShaper class.
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/klee.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/quest/SamplingShaper.hpp"
#include "axom/quest/util/mesh_helpers.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"

#if !defined(AXOM_USE_CONDUIT) || !defined(AXOM_USE_BUMP) || !defined(AXOM_USE_SIDRE)
  #error "Quest's Blueprint SamplingShaper tests require Conduit, Bump, and Sidre."
#endif

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#include <map>
#include <string>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;
namespace slic = axom::slic;

TEST(SamplingShaperBlueprintTest, sidre_blueprint_quadrature_persists)
{
  sidre::DataStore dataStore;
  auto* meshGroup = dataStore.getRoot()->createGroup("mesh");

  const primal::BoundingBox<double, 2> bbox {{0., 0.}, {1., 1.}};
  const axom::NumericArray<int, 2> res {{2, 2}};
  quest::util::make_unstructured_blueprint_box_mesh_2d(meshGroup, bbox, res, "mesh", "coords");

  constexpr axom::IndexType cellCount = 4;
  const std::string backgroundVolFracName = quest::shaping::volumeFractionFieldName("background");
  const std::string backgroundMatInOutName = quest::shaping::materialInOutFieldName("background");
  auto* fieldGroup = meshGroup->createGroup(axom::fmt::format("fields/{}", backgroundVolFracName));
  fieldGroup->createViewString("association", "element");
  fieldGroup->createViewString("topology", "mesh");
  auto* valuesView =
    fieldGroup->createViewAndAllocate("values", axom::sidre::DataTypeId::FLOAT64_ID, cellCount);
  auto* values = static_cast<double*>(valuesView->getVoidPtr());
  for(axom::IndexType i = 0; i < cellCount; ++i)
  {
    values[i] = 1.;
  }

  klee::ShapeSet shapeSet;
  quest::SamplingShaper shaper(axom::runtime_policy::Policy::seq,
                               axom::policyToDefaultAllocatorID(axom::runtime_policy::Policy::seq),
                               shapeSet,
                               meshGroup,
                               "mesh");
  shaper.setSamplingResolution(2);

  auto* bpMeshNode = shaper.getBlueprintMeshNode();
  ASSERT_NE(bpMeshNode, nullptr);

  std::map<std::string, conduit::Node*> initialVolumeFractions;
  initialVolumeFractions["background"] =
    &bpMeshNode->fetch_existing(axom::fmt::format("fields/{}", backgroundVolFracName));
  shaper.importInitialVolumeFractions(initialVolumeFractions);

  EXPECT_TRUE(meshGroup->hasGroup("coordsets/quadrature_points"));
  EXPECT_TRUE(meshGroup->hasGroup("topologies/quadrature_points"));
  EXPECT_TRUE(meshGroup->hasGroup("fields/originalElements"));
  EXPECT_TRUE(meshGroup->hasGroup("fields/quadratureWeights"));
  EXPECT_TRUE(meshGroup->hasGroup(axom::fmt::format("fields/{}", backgroundMatInOutName)));

  conduit::Node refreshedMesh;
  meshGroup->createNativeLayout(refreshedMesh);
  EXPECT_TRUE(refreshedMesh.has_path("coordsets/quadrature_points"));
  EXPECT_TRUE(refreshedMesh.has_path("topologies/quadrature_points"));
  EXPECT_TRUE(refreshedMesh.has_path("fields/originalElements/values"));
  EXPECT_TRUE(refreshedMesh.has_path("fields/quadratureWeights/values"));
  EXPECT_TRUE(
    refreshedMesh.has_path(axom::fmt::format("fields/{}/values", backgroundMatInOutName)));
}

int main(int argc, char* argv[])
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  slic::SimpleLogger logger(slic::message::Info);

  const int result = RUN_ALL_TESTS();
  return result;
}
