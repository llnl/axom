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

#include "axom/quest/InOutOctree.hpp"

#include "quest_test_utilities.hpp"

#include <cstdlib>
#include <limits>

// Uncomment the define below for true randomized points
#ifndef INOUT_OCTREE_TESTER_SHOULD_SEED
//  #define INOUT_OCTREE_TESTER_SHOULD_SEED
#endif

#ifdef INOUT_OCTREE_TESTER_SHOULD_SEED
  #include <ctime>  // for time() used by srand()
#endif

namespace
{
const int NUM_PT_TESTS = 10000;
const int DIM = 2;
using Octree2D = axom::quest::InOutOctree<DIM>;
using GeometricBoundingBox = Octree2D::GeometricBoundingBox;
using SpacePt = Octree2D::SpacePt;

using SpaceVector = Octree2D::SpaceVector;
using GridPt = Octree2D::GridPt;
using BlockIndex = Octree2D::BlockIndex;
}  // namespace

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

TEST(quest_inout_quadtree, triangle_boundary_mesh)
{
  SLIC_INFO("*** Exercises InOutOctree queries for several thresholds.\n");

  namespace mint = axom::mint;
  namespace quest = axom::quest;

  std::vector<double> thresholds = {1E-9, 1E-5, 1E-1, 0.};

  for(auto thresh : thresholds)
  {
    // create a simple mesh on the boundary of an equilateral triangle
    std::shared_ptr<mint::Mesh> mesh = [=]() {
      auto mesh = std::make_shared<mint::UnstructuredMesh<mint::SINGLE_SHAPE>>(DIM, mint::SEGMENT);
      mesh->appendNode(1, 1);
      mesh->appendNode(4, 1);
      mesh->appendNode(2.5, 3 * sqrt(3) / 2.);

      axom::IndexType cell[4] = {0, 1, 2, 0};
      mesh->appendCell(&cell[0]);
      mesh->appendCell(&cell[1]);
      mesh->appendCell(&cell[2]);

      return mesh;
    }();

    GeometricBoundingBox bbox = computeBoundingBox(*mesh);

    Octree2D octree(bbox, mesh);
    octree.setVertexWeldThreshold(thresh);

    octree.generateIndex();

    SpacePt queryInside =
      quest::utilities::getCentroid(getVertex(*mesh, 0), getVertex(*mesh, 1), getVertex(*mesh, 2));
    SpacePt queryOutside = SpacePt(2. * bbox.getMax().array());

    EXPECT_TRUE(octree.within(queryInside));
    EXPECT_FALSE(octree.within(queryOutside));
  }
}

TEST(quest_inout_quadtree, circle_mesh)
{
  SLIC_INFO("*** Exercises InOutOctree over the boundary of a circle.\n");

  namespace mint = axom::mint;
  namespace quest = axom::quest;

  for(int num_segments : {3, 100, 1000, 10000})
  {
    for(double radius : {1. / 3., 1., sqrt(2.), 1234.5678})
    {
      ASSERT_TRUE(num_segments >= 3);
      std::shared_ptr<mint::Mesh> mesh {quest::utilities::make_circle_mesh_2d(radius, num_segments)};
      //mint::write_vtk(mesh,axom::fmt::format("circle_mesh_r{:.3f}_s{:06}.vtk", radius, num_segments));

      GeometricBoundingBox bbox = computeBoundingBox(*mesh).scale(1.2);

      Octree2D octree(bbox, mesh);
      octree.generateIndex();

      SpacePt queryInside = SpacePt {0, 0};
      SpacePt queryOutside = SpacePt {2. * radius, 2. * radius};

      EXPECT_TRUE(octree.within(queryInside));
      EXPECT_FALSE(octree.within(queryOutside));

      // Regression: Test a point that was incorrectly categorized when num_segments was 3
      if(radius <= 1.)
      {
        SpacePt additionalQuery {1.1494147076163739, 0.51644760397625789};
        EXPECT_FALSE(octree.within(additionalQuery));
      }

      // Determine a confidence interval for status of query point against discretized circle
      // Uses the midpoint of a segment for the inner radius and a vertex for the outer one
      SpacePt a, b;
      mesh->getNode(0, a.data());
      mesh->getNode(1, b.data());
      auto segmentCentroid = quest::utilities::getCentroid(a, b);
      double innerConfidence = 0.95 * SpaceVector(segmentCentroid).norm();
      double outerConfidence = 1.05 * radius;

      int insideCount = 0;
      int outsideCount = 0;
      int uncertainCount = 0;

      for(int i = 0; i < NUM_PT_TESTS; ++i)
      {
        SpacePt queryPt = quest::utilities::randomSpacePt<2>(0, 1.25 * radius);
        const double mag = SpaceVector(queryPt).norm();
        const bool expectInside = mag < innerConfidence;
        const bool expectOutside = mag > outerConfidence;

        if(expectInside)
        {
          EXPECT_TRUE(octree.within(queryPt)) << "Query point: " << queryPt << "; norm: " << mag;
          ++insideCount;
        }
        else if(expectOutside)
        {
          EXPECT_FALSE(octree.within(queryPt)) << "Query point: " << queryPt << "; norm: " << mag;
          ++outsideCount;
        }
        else
        {
          // Not sure if point should be inside or outside
          ++uncertainCount;
        }
      }

      // Output some stats about the query
      SLIC_INFO(
        axom::fmt::format("Queried quadtree over circle mesh of radius {}"
                          " defined by {} segments using {} query points. \n "
                          "Of which: "
                          " {:.2f}% were known to be inside; {:.2f}% were known to be outside; "
                          " and {:.2f}% were too close to the boundary for our simple model.",
                          radius,
                          num_segments,
                          NUM_PT_TESTS,
                          100. * (static_cast<double>(insideCount) / NUM_PT_TESTS),
                          100. * (static_cast<double>(outsideCount) / NUM_PT_TESTS),
                          100. * (static_cast<double>(uncertainCount) / NUM_PT_TESTS)));
    }
  }
}

TEST(quest_inout_quadtree, on_surface_points)
{
  // Regression test for https://github.com/LLNL/axom/issues/611 (2D)
  // Query points on the surface should be marked as inside the surface.
  // This test also checks for consistency with the winding number results.

  namespace mint = axom::mint;
  namespace quest = axom::quest;
  namespace primal = axom::primal;

  using Polygon2D = primal::Polygon<double, 2>;

  // lambda to linearly interpolate a point on an edge of a polygon at parameter 0 <= t <= 1
  auto lerp_edge = [](const Polygon2D& poly, int edge, double t) {
    SLIC_ASSERT(edge >= 0 && edge <= poly.numVertices());
    SLIC_ASSERT(t >= 0. && t <= 1.);
    return SpacePt::lerp(poly[edge], poly[(edge + 1) == poly.numVertices() ? 0 : edge + 1], t);
  };

  const Polygon2D unitSquare({SpacePt {0., 0.}, SpacePt {1., 0.}, SpacePt {1., 1.}, SpacePt {0., 1.}});

  // Create mesh of unit square, with several edge refinements
  for(int segsPerSide : {1, 2, 3, 5})
  {
    // Build the perimeter vertices (CCW; corners not duplicated).
    axom::Array<SpacePt> verts;
    for(int edge = 0; edge < unitSquare.numVertices(); ++edge)
    {
      for(int s = 0; s < segsPerSide; ++s)
      {
        const double t = static_cast<double>(s) / segsPerSide;
        verts.push_back(lerp_edge(unitSquare, edge, t));
      }
    }
    const int nverts = static_cast<int>(verts.size());

    // Segment mesh for the quadtree.
    std::shared_ptr<mint::Mesh> mesh = [&]() {
      auto m = std::make_shared<mint::UnstructuredMesh<mint::SINGLE_SHAPE>>(DIM, mint::SEGMENT);
      for(const auto& v : verts)
      {
        m->appendNode(v[0], v[1]);
      }
      for(int i = 0; i < nverts; ++i)
      {
        axom::IndexType cell[2] = {i, (i + 1) % nverts};
        m->appendCell(cell);
      }
      return m;
    }();

    // Build quadtree over the linearized unit square
    GeometricBoundingBox bbox = computeBoundingBox(*mesh);
    Octree2D octree(bbox, mesh);
    octree.generateIndex();

    // sanity check on some interior and exterior points
    EXPECT_TRUE(octree.within(SpacePt {0.5, 0.5}));
    EXPECT_TRUE(octree.within(SpacePt {0.25, 0.75}));
    EXPECT_FALSE(octree.within(SpacePt {0.5, 1.5}));
    EXPECT_FALSE(octree.within(SpacePt {1.5, 0.5}));
    EXPECT_FALSE(octree.within(SpacePt {2.0, 2.0}));

    // We will use the winding number over the polygon as an oracle
    Polygon2D poly;
    for(const auto& v : verts)
    {
      poly.addVertex(v);
    }

    // Create a set of query points on the boundary, starting from unit square vertices and edge midpoints
    axom::Array<SpacePt> queryPoints {SpacePt {0.5, 1.0},
                                      SpacePt {0.5, 0.0},
                                      SpacePt {0.0, 0.5},
                                      SpacePt {1.0, 0.5},
                                      SpacePt {0.0, 0.0},
                                      SpacePt {1.0, 1.0}};

    // Add regression case from the issue
    queryPoints.push_back(SpacePt {0.370667, 1.0});

    // Add a dense set of points along the edges
    const int edgeSamples = 500 / poly.numVertices();
    for(int edge = 0; edge < poly.numVertices(); ++edge)
    {
      for(int sample = 0; sample < edgeSamples; ++sample)
      {
        queryPoints.push_back(lerp_edge(poly, edge, static_cast<double>(sample) / edgeSamples));
      }
    }

    // Run the tests
    for(const auto& q : queryPoints)
    {
      bool isOnEdge {};
      const int wn = primal::winding_number(q, poly, isOnEdge, /*includeBoundary=*/true);

      // Sanity check on the oracle: these points really are on the boundary.
      EXPECT_TRUE(isOnEdge) << axom::fmt::format("Oracle: point {} should be on the boundary", q);
      EXPECT_NE(0, wn) << axom::fmt::format("Oracle: point {} should have nonzero winding number", q);

      // The actual regression assertion: on-surface points are 'within'.
      EXPECT_TRUE(octree.within(q))
        << axom::fmt::format("Boundary point {} should be within the surface (segsPerSide={})",
                             q,
                             segsPerSide);
    }
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
