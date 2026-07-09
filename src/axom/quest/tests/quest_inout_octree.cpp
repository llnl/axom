// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
TEST(quest_inout_octree, on_surface_points)
{
  // Regression test for https://github.com/LLNL/axom/issues/611 (3D)
  //
  // Query point lying exactly on the surface should be marked as inside the surface.
  // This test builds a unit cube as a closed triangle surface (12 triangles) and
  // checks that points sampled exactly on the surface are reported as 'within'.
  // It compares against generalized winding number summed over the 12 triangles.

  namespace mint = axom::mint;
  namespace quest = axom::quest;
  namespace primal = axom::primal;

  using Point3D = primal::Point<double, 3>;
  using Triangle3D = primal::Triangle<double, 3>;

  constexpr double x_lo = 0.;
  constexpr double y_lo = 0.;
  constexpr double z_lo = 0.;
  constexpr double x_mid = 0.5;
  constexpr double y_mid = 0.5;
  constexpr double z_mid = 0.5;
  constexpr double x_hi = 1.;
  constexpr double y_hi = 1.;
  constexpr double z_hi = 1.;

  // 8 corners of the unit cube.
  axom::Array<Point3D> V {Point3D {x_lo, y_lo, z_lo},
                          Point3D {x_hi, y_lo, z_lo},
                          Point3D {x_hi, y_hi, z_lo},
                          Point3D {x_lo, y_hi, z_lo},
                          Point3D {x_lo, y_lo, z_hi},
                          Point3D {x_hi, y_lo, z_hi},
                          Point3D {x_hi, y_hi, z_hi},
                          Point3D {x_lo, y_hi, z_hi}};

  // 12 triangles (2 per face), wound counter-clockwise as seen from outside (outward normals).
  axom::Array<std::tuple<int, int, int>> TRI {{0, 3, 2},  // z=0 bottom
                                              {0, 2, 1},
                                              {4, 5, 6},  // z=1 top
                                              {4, 6, 7},
                                              {0, 1, 5},  // y=0 front
                                              {0, 5, 4},
                                              {3, 7, 6},  // y=1 back
                                              {3, 6, 2},
                                              {0, 4, 7},  // x=0 left
                                              {0, 7, 3},
                                              {1, 2, 6},  // x=1 right
                                              {1, 6, 5}};

  // Build a mesh over the triangles and an octree over the mesh
  std::shared_ptr<mint::Mesh> mesh = [&V, &TRI]() {
    auto m = std::make_shared<mint::UnstructuredMesh<mint::SINGLE_SHAPE>>(DIM, mint::TRIANGLE);
    for(const auto& v : V)
    {
      m->appendNode(v[0], v[1], v[2]);
    }
    for(const auto& [t0, t1, t2] : TRI)
    {
      axom::IndexType cell[3] = {t0, t1, t2};
      m->appendCell(cell);
    }
    return m;
  }();

  // Use an explicit vertex-weld threshold so the test's on-surface tolerance
  // and the octree's on-surface tolerance are the same quantity.
  const double weldThresh = 1e-6;
  GeometricBoundingBox bbox = computeBoundingBox(*mesh);
  Octree3D octree(bbox, mesh);
  octree.setVertexWeldThreshold(weldThresh);
  octree.generateIndex();

  // The octree treats points within its weld threshold of the surface as 'inside'
  // The oracle below uses the same tolerance for consistency.
  const double edgeTol = octree.getVertexWeldThreshold();

  // Matching triangle list for the winding-number oracle.
  axom::Array<Triangle3D> tris;
  for(const auto& [t0, t1, t2] : TRI)
  {
    tris.emplace_back(V[t0], V[t1], V[t2]);
  }

  // Oracle: on-surface (by distance) => within; else sign of rounded GWN sum.
  // The default tolerance matches the octree's vertex-weld threshold.
  auto expectedWithin = [&tris, edgeTol](const Point3D& q, double edge_tol = -1.0) -> bool {
    if(edge_tol < 0.0)
    {
      edge_tol = edgeTol;
    }
    const double edge_tol_2 = edge_tol * edge_tol;
    double wn = 0.0;
    for(const auto& tri : tris)
    {
      bool onThis {};
      wn += primal::winding_number(q, tri, onThis, edge_tol, edge_tol);
      if(onThis || primal::squared_distance(q, tri) <= edge_tol_2)
      {
        return true;  // on the surface
      }
    }
    return std::lround(wn) != 0;
  };

  // Adds some hand-picked on-surface points: face centers, edge midpoints, corners,
  // and points on the shared diagonal edge between the two triangles of a face.
  axom::Array<SpacePt> onSurface {
    // face centers
    SpacePt {x_mid, y_mid, z_lo},
    SpacePt {x_mid, y_mid, z_hi},
    SpacePt {x_mid, y_lo, z_mid},
    SpacePt {x_mid, y_hi, z_mid},
    SpacePt {x_lo, y_mid, z_mid},
    SpacePt {x_hi, y_mid, z_mid},
    // edge midpoints
    SpacePt {x_mid, y_lo, z_lo},
    SpacePt {x_lo, y_mid, z_lo},
    SpacePt {x_hi, y_hi, z_mid},
    // corners
    SpacePt {x_lo, y_lo, z_lo},
    SpacePt {x_hi, y_hi, z_hi},
    SpacePt {x_hi, y_lo, z_hi},
    // off-center on faces
    SpacePt {0.3, 0.7, z_hi},
    SpacePt {x_hi, 0.25, 0.6},
    SpacePt {0.4, 0.4, z_lo}  // on a face diagonal edge
  };

  // Also add dense set of samples on each triangle
  const int bres = 8;
  for(const auto& tri : tris)
  {
    for(int a = 0; a <= bres; ++a)
    {
      const double u = static_cast<double>(a) / bres;
      for(int b = 0; a + b <= bres; ++b)
      {
        const double v = static_cast<double>(b) / bres;
        onSurface.push_back(tri.baryToPhysical(SpacePt {u, v, 1. - u - v}));
      }
    }
  }

  // Run the on-surface comparisons
  for(const auto& q : onSurface)
  {
    EXPECT_TRUE(expectedWithin(q)) << "Oracle: point " << q << " should be on/within the surface";
    EXPECT_TRUE(octree.within(q)) << "On-surface point " << q << " should be within the surface";
  }

  // Sanity check for several interior query points
  for(const auto& q_interior : {SpacePt {0.5, 0.5, 0.5}, SpacePt {0.25, 0.75, 0.5}})
  {
    EXPECT_TRUE(expectedWithin(q_interior));
    EXPECT_TRUE(octree.within(q_interior));
  }

  // Sanity check for several exterior query points
  for(const auto& q_exterior :
      {SpacePt {0.5, 0.5, 1.5}, SpacePt {1.5, 0.5, 0.5}, SpacePt {2.0, 2.0, 2.0}})
  {
    EXPECT_FALSE(expectedWithin(q_exterior));
    EXPECT_FALSE(octree.within(q_exterior));
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
