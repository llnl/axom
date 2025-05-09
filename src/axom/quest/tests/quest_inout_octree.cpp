// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"

// _quest_inout_cpp_include_start
#include "axom/quest/InOutOctree.hpp"
// _quest_inout_cpp_include_end

#include "quest_test_utilities.hpp"

namespace
{
const int NUM_PT_TESTS = 50000;
const int DIM = 3;
}  // namespace

// _quest_inout_cpp_typedef_start
using Octree3D = axom::quest::InOutOctree<DIM>;

using GeometricBoundingBox = Octree3D::GeometricBoundingBox;
using SpacePt = Octree3D::SpacePt;
// _quest_inout_cpp_typedef_end

using SpaceVector = Octree3D::SpaceVector;
using GridPt = Octree3D::GridPt;
using BlockIndex = Octree3D::BlockIndex;

#include <cstdlib>
#include <limits>

// Uncomment the line below for true randomized points
#ifndef INOUT_OCTREE_TESTER_SHOULD_SEED
//  #define INOUT_OCTREE_TESTER_SHOULD_SEED
#endif

#ifdef INOUT_OCTREE_TESTER_SHOULD_SEED
  #include <ctime>  // for time() used by srand()
#endif

/// Returns a SpacePt corresponding to the given vertex id \a vIdx  in \a mesh
SpacePt getVertex(axom::mint::Mesh& mesh, int vIdx)
{
  SpacePt pt;
  mesh.getNode(vIdx, pt.data());

  return pt;
}

GeometricBoundingBox computeBoundingBox(axom::mint::Mesh& mesh)
{
  GeometricBoundingBox bbox;
  for(int i = 0; i < mesh.getNumberOfNodes(); ++i)
  {
    bbox.addPoint(getVertex(mesh, i));
  }

  return bbox;
}

/// Runs randomized inout queries on an octahedron mesh
void queryOctahedronMesh(std::shared_ptr<axom::mint::Mesh>& mesh, const GeometricBoundingBox& bbox)
{
  const double bbMin = bbox.getMin()[0];
  const double bbMax = bbox.getMax()[0];

  // _quest_inout_cpp_init_start
  Octree3D octree(bbox, mesh);
  octree.generateIndex();
  // _quest_inout_cpp_init_end

  SLIC_INFO("Testing point containment on an octahedron surface mesh.");
  SLIC_INFO("Note: Points on the surface might issue a warning, "
            << "but should be considered outside the surface.");
  SLIC_INFO("--[==[");

  // Query the mesh containment
  axom::utilities::Timer timer(true);
  for(int i = 0; i < NUM_PT_TESTS; ++i)
  {
    SpacePt pt;

    // test a few special points and a lot of random points
    switch(i)
    {
    case 0:
    case 1:
    case 2:  // Test point at mesh vertices
    case 3:
    case 4:
    case 5:
      pt = getVertex(*mesh, i);
      break;
    case 6:
    case 7:
    case 8:
    case 9:  // Test point at triangle centers
    case 10:
    case 11:
    case 12:
    case 13:
    {
      axom::IndexType tIdx = (i - 6);
      GridPt vertInds;
      mesh->getCellNodeIDs(tIdx, vertInds.data());

      pt = axom::quest::utilities::getCentroid(getVertex(*mesh, vertInds[0]),
                                               getVertex(*mesh, vertInds[1]),
                                               getVertex(*mesh, vertInds[2]));
    }
    break;

    case 14:
    case 15:
    case 16:
    case 17:  // Test point at edge centers
    case 18:
    case 19:
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
    {
      // Define explicit vertex indices for edges of octahedron
      const int v1[] = {0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5};
      const int v2[] = {1, 2, 3, 4, 2, 3, 4, 1, 1, 2, 3, 4};

      int eIdx = (i - 14);
      pt =
        axom::quest::utilities::getCentroid(getVertex(*mesh, v1[eIdx]), getVertex(*mesh, v2[eIdx]));
    }
    break;

    case 26:  // origin
      pt = SpacePt();
      break;
    case 27:  // outside bounding box
      pt = SpacePt(2 * bbMin);
      break;
    case 28:  // outside bounding box
      pt = SpacePt(2 * bbMax);
      break;
    default:  // random points in bounding box
      pt = axom::quest::utilities::randomSpacePt<DIM>(bbMin, bbMax);
      break;
    }

    // Unit octahedron is a unit sphere under the L1 metric. Points are
    // inside when the sum of coordinate magnitudes is less than one
    double absCoordSum = std::abs(pt[0]) + std::abs(pt[1]) + std::abs(pt[2]);

    // For the time being, we allow the within() test to fail when the
    // query point is sufficiently close to the surface
    bool expectInside = absCoordSum < 1.;
    EXPECT_TRUE(octree.within(pt) == expectInside || axom::utilities::isNearlyEqual(absCoordSum, 1.))
      << "Point " << pt << " was not " << (expectInside ? "inside" : "outside")
      << " surface of octahedron as expected."
      << " Sum of absolute values of coords was: " << absCoordSum
      << " and point is inside when this is less than 1.";
  }
  timer.stop();
  SLIC_INFO("--]==]");

  SLIC_INFO("-- querying octahedron with " << NUM_PT_TESTS << " points took " << timer.elapsed()
                                           << " seconds.");
  SLIC_INFO("***");
}

TEST(quest_inout_octree, octahedron_mesh)
{
  SLIC_INFO("*** This test creates a simple mesh of an octahedron "
            << " and tests point containment.\n");

  // Generate the InOutOctree
  std::shared_ptr<axom::mint::Mesh> mesh(axom::quest::utilities::make_octahedron_mesh());
  // axom::mint::write_vtk(mesh, "octahedron.vtk");

  ///
  SpacePt ptNeg1(-1.);
  SpacePt ptPos1(1.);
  GeometricBoundingBox bbox1(ptNeg1, ptPos1);
  SLIC_INFO("Testing InOutOctree on octahedron mesh with bounding box " << bbox1);
  queryOctahedronMesh(mesh, bbox1);

  ///
  SpacePt ptNeg2(-2.);
  SpacePt ptPos2(2.);
  GeometricBoundingBox bbox2(ptNeg2, ptPos2);
  SLIC_INFO("Testing InOutOctree on octahedron mesh with bounding box " << bbox2);
  queryOctahedronMesh(mesh, bbox2);

  bbox2.shift(SpaceVector(0.01));
  SLIC_INFO("Testing InOutOctree on octahedron mesh with shifted bounding box " << bbox2);
  queryOctahedronMesh(mesh, bbox2);
}

TEST(quest_inout_octree, tetrahedron_mesh)
{
  SLIC_INFO("*** Exercises InOutOctree queries for several thresholds.\n");

  namespace mint = axom::mint;
  namespace quest = axom::quest;

  std::vector<double> thresholds = {1E-9, 1E-5, 1E-1, 0.};

  for(auto thresh : thresholds)
  {
    std::shared_ptr<mint::Mesh> mesh {quest::utilities::make_tetrahedron_mesh()};
    GeometricBoundingBox bbox = computeBoundingBox(*mesh);

    Octree3D octree(bbox, mesh);
    octree.setVertexWeldThreshold(thresh);

    octree.generateIndex();

    SpacePt queryInside = quest::utilities::getCentroid(getVertex(*mesh, 0),
                                                        getVertex(*mesh, 1),
                                                        getVertex(*mesh, 2),
                                                        getVertex(*mesh, 3));
    SpacePt queryOutside = SpacePt(2. * bbox.getMax().array());

    EXPECT_TRUE(octree.within(queryInside));
    EXPECT_FALSE(octree.within(queryOutside));
  }
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

#ifdef INOUT_OCTREE_TESTER_SHOULD_SEED
  std::srand(std::time(0));
#else
  std::srand(105);
#endif

  int result = RUN_ALL_TESTS();
  return result;
}
