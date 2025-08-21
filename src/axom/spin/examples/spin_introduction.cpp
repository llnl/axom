// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file spin_introduction.cpp
 *  \brief This example code is a demonstration of Axom's spin component.
 *
 *  This file shows how to use spin with primal to instantiate, populate,
 *  and use a spatial index.  Running the executable produced from
 *  this file will produce a collection of Asymptote source files.  When
 *  compiled, the Asymptote files produce the figures that accompany spin's
 *  Sphinx documentation.
 */

/* This example code contains snippets used in the spin Sphinx documentation.
 * They begin and end with comments such as
 *
 * prims_header_start
 * prims_header_end
 * clip_header_start
 * clip_header_end
 * closest_point_header_start
 * closest_point_header_end
 *
 * each prepended with an underscore.
 */

#include "axom/config.hpp"

// Axom primitives
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"

// Axom operations
#include "axom/primal/operators/compute_bounding_box.hpp"
#include "axom/primal/operators/intersect.hpp"

// C++ headers
#include <cmath>  // do we need this?
#include <iostream>
#include <fstream>
#include <vector>

#include "axom/fmt.hpp"

// almost all our examples are in 3D
constexpr int in3D = 3;
constexpr int in2D = 2;

// primitives represented by doubles in 3D
using BoundingBoxType = axom::primal::BoundingBox<double, in3D>;
using PointType = axom::primal::Point<double, in3D>;
using TriangleType = axom::primal::Triangle<double, in3D>;

// Axom spatial index
// _ugrid_triintersect_header_start
#include "axom/spin/UniformGrid.hpp"
// the UniformGrid will store ints ("thing" indexes) in 3D
using UniformGridType = axom::spin::UniformGrid<int, in3D>;
// _ugrid_triintersect_header_end

// _bvh_header_start
#include "axom/spin/BVH.hpp"
// the BVH is in 2D, storing an index to 2D triangles
using BVH2DType = axom::spin::BVH<in2D>;
// supporting classes
using BoundingBox2DType = axom::primal::BoundingBox<double, in2D>;
using Point2DType = axom::primal::Point<double, in2D>;
using Triangle2DType = axom::primal::Triangle<double, in2D>;
// _bvh_header_end

// _igrid_header_start
#include "axom/spin/ImplicitGrid.hpp"

// the ImplicitGrid will be in 2D
using IGridT = axom::spin::ImplicitGrid<in2D>;

// useful derived types
using IGridCell = typename IGridT::GridCell;
using ISpacePt = typename IGridT::SpacePoint;
using IBBox = typename IGridT::SpatialBoundingBox;
using IBitsetType = typename IGridT::BitsetType;
using IBitsetIndexType = typename IBitsetType::Index;

// some functions we'll use
bool expensiveTest(ISpacePt& query, Triangle2DType& tri);
void makeTriangles(std::vector<Triangle2DType>& tris);
// _igrid_header_end

// _rectlattice_header_start
#include "axom/spin/RectangularLattice.hpp"
// We'll be using a 2D lattice with space coordinates of type double
// and cell coordinates of type int.
using RectLatticeType = axom::spin::RectangularLattice<in2D, double, int>;
// Get the types of coordinates and bounding box from the RectangularLattice
using RLGridCell = RectLatticeType::GridCell;
using RLSpacePt = RectLatticeType::SpacePoint;
using RLSpaceVec = RectLatticeType::SpaceVector;
using RLBBox = RectLatticeType::SpatialBoundingBox;
// _rectlattice_header_end

// _morton_header_start
#include "axom/spin/MortonIndex.hpp"
#include <unordered_map>
// The PointHash will allow us to use integral N-D coordinates as a hash key.
// This example will use RectangularLattice grid cell coordinates as keys to
// a std::unordered_map.
using PointHashType = axom::spin::PointHash<RLGridCell::CoordType>;
// Here's a class defined elsewhere that will do some work on a point.
class DataContainer;
using MapType = std::unordered_map<RLGridCell, DataContainer, PointHashType>;
// _morton_header_end

// _octree_header_start
#include "axom/spin/SpatialOctree.hpp"

using LeafNodeType = axom::spin::BlockData;

using OctreeType = axom::spin::SpatialOctree<in3D, LeafNodeType>;
using OctBlockIndex = OctreeType::BlockIndex;
using OctSpacePt = OctreeType::SpacePt;
using OctBBox = OctreeType::GeometricBoundingBox;
// _octree_header_end

class DataContainer
{
public:
  DataContainer() : count(0) { }
  void registerPoint(RLSpacePt p)
  {
    (void)p;
    count += 1;
  }

  int count;
};

std::vector<RLSpacePt> generatePoints()
{
  return std::vector<RLSpacePt> {RLSpacePt {2.8, 3.1},
                                 RLSpacePt {2.9, 1.9},
                                 RLSpacePt {1.0, 3.6},
                                 RLSpacePt {.28, 310},
                                 RLSpacePt {-460, -0.9}};
}

void demoMorton()
{
  std::cout << "----- demoMorton -----" << std::endl;

  // _morton_use_start
  // Make a RectangularLattice to bin query points.
  double origin[] = {-0.6, -0.2};
  double spacing[] = {1.2, 0.8};
  RectLatticeType lat(origin, spacing);

  // Make the map from grid point to DataContainer
  MapType map;

  // For several query points, create a DataContainer if necessary and register the point.
  std::vector<RLSpacePt> pts = generatePoints();
  for(RLSpacePt p : pts)
  {
    RLGridCell g = lat.gridCell(p);
    DataContainer dat;
    if(map.count(g) > 0)
    {
      dat = map[g];
    }
    dat.registerPoint(p);
    map[g] = dat;
  }

  // Report on what was registered.
  for(auto iter : map)
  {
    RLGridCell g = iter.first;
    DataContainer dat = iter.second;
    std::cout << "Grid cell " << g << " holds " << dat.count << " points." << std::endl;
  }
  // _morton_use_end

  std::cout << std::endl;
}

void demoRectangularLattice()
{
  std::cout << "----- demoRectangularLattice -----" << std::endl;

  // _rectlattice_use_start
  // Origin and spacing
  double origin[] = {-0.6, -0.2};
  double spacing[] = {1.2, 0.8};

  // Instantiate a RectangularLattice.
  // Other constructors allow the use of Point and Vector objects.
  RectLatticeType lat(origin, spacing);

  // Query point (2.0, 1.2) should be in grid cell (2, 1)
  RLSpacePt pA {2.0, 1.2};
  RLGridCell cA = lat.gridCell(pA);
  std::cout << "Point " << pA << " is in grid cell " << cA << " (should be (2, 1))" << std::endl;

  // Query point (2.3, 0.8) should also be in grid cell (2, 1)
  RLSpacePt pB {2.3, 0.8};
  RLGridCell cB = lat.gridCell(pB);
  std::cout << "Point " << pB << " is in grid cell " << cB << " (should be (2, 1))" << std::endl;

  // What is the lowest corner and bounding box of the shared cell?
  RLSpacePt cellcorner = lat.spacePoint(cB);
  RLBBox cellbbox = lat.cellBounds(cB);
  std::cout << "The lower corner of the grid cell is " << cellcorner << " and its bounding box is "
            << cellbbox << std::endl;
  // _rectlattice_use_end

  std::cout << std::endl;
}

// _naive_triintersect_start
void findTriIntersectionsNaively(std::vector<TriangleType>& tris,
                                 std::vector<std::pair<int, int>>& clashes)
{
  int tcount = static_cast<int>(tris.size());

  for(int i = 0; i < tcount; ++i)
  {
    TriangleType& t1 = tris[i];
    for(int j = i + 1; j < tcount; ++j)
    {
      TriangleType& t2 = tris[j];
      if(intersect(t1, t2))
      {
        clashes.push_back(std::make_pair(i, j));
      }
    }
  }
}
// _naive_triintersect_end

BoundingBoxType findBbox(std::vector<TriangleType>& tris)
{
  BoundingBoxType bbox;

  for(size_t i = 0; i < tris.size(); ++i)
  {
    bbox.addPoint(tris[i][0]);
    bbox.addPoint(tris[i][1]);
    bbox.addPoint(tris[i][2]);
  }

  return bbox;
}

BoundingBoxType findBbox(TriangleType& tri)
{
  BoundingBoxType bbox;

  bbox.addPoint(tri[0]);
  bbox.addPoint(tri[1]);
  bbox.addPoint(tri[2]);

  return bbox;
}

BoundingBox2DType findBbox(Triangle2DType& tri)
{
  BoundingBox2DType bbox;

  bbox.addPoint(tri[0]);
  bbox.addPoint(tri[1]);
  bbox.addPoint(tri[2]);

  return bbox;
}

// _ugrid_build_start
BoundingBoxType findBbox(std::vector<TriangleType>& tris);
BoundingBoxType findBbox(TriangleType& tri);

UniformGridType* buildUniformGrid(std::vector<TriangleType>& tris)
{
  // Prepare to construct the UniformGrid.
  BoundingBoxType allbbox = findBbox(tris);
  const PointType& minBBPt = allbbox.getMin();
  const PointType& maxBBPt = allbbox.getMax();

  int tcount = static_cast<int>(tris.size());

  // The number of bins along one side of the UniformGrid.
  // This is a heuristic.
  int res = (int)(1 + std::pow(tcount, 1 / 3.));
  int ress[3] = {res, res, res};

  // Construct the UniformGrid with minimum point, maximum point,
  // and number of bins along each side.  Then insert the triangles.
  UniformGridType* ugrid = new UniformGridType(minBBPt.data(), maxBBPt.data(), ress);
  for(int i = 0; i < tcount; ++i)
  {
    TriangleType& t1 = tris[i];
    BoundingBoxType bbox = findBbox(t1);
    ugrid->insert(bbox, i);
  }

  return ugrid;
}
// _ugrid_build_end

// _ugrid_candidate_start
void findNeighborCandidates(TriangleType& t1, int i, UniformGridType* ugrid, std::vector<int>& neighbors)
{
  BoundingBoxType bbox = findBbox(t1);

  // Get all the bins t1 occupies
  const std::vector<int> bToCheck = ugrid->getBinsForBbox(bbox);
  size_t checkcount = bToCheck.size();

  // Load all the triangles in these bins whose indices are greater than i into a vector.
  for(size_t curb = 0; curb < checkcount; ++curb)
  {
    axom::ArrayView<int> ntlist = ugrid->getBinContents(bToCheck[curb]);
    for(axom::IndexType j = 0; j < ntlist.size(); ++j)
    {
      if(ntlist[j] > i)
      {
        neighbors.push_back(ntlist[j]);
      }
    }
  }

  // Sort the neighboring triangles, and throw out duplicates.
  // This is not strictly necessary but saves some calls to intersect().
  std::sort(neighbors.begin(), neighbors.end());
  std::vector<int>::iterator jend = std::unique(neighbors.begin(), neighbors.end());
  neighbors.erase(jend, neighbors.end());
}
// _ugrid_candidate_end

// _ugrid_triintersect_start
void findTriIntersectionsAccel(std::vector<TriangleType>& tris,
                               UniformGridType* ugrid,
                               std::vector<std::pair<int, int>>& clashes)
{
  int tcount = static_cast<int>(tris.size());

  // For each triangle t1,
  for(int i = 0; i < tcount; ++i)
  {
    TriangleType& t1 = tris[i];
    std::vector<int> neighbors;
    findNeighborCandidates(t1, i, ugrid, neighbors);

    // Test for intersection between t1 and each of its neighbors.
    int ncount = static_cast<int>(neighbors.size());
    for(int n = 0; n < ncount; ++n)
    {
      int j = neighbors[n];
      TriangleType& t2 = tris[j];
      if(axom::primal::intersect(t1, t2))
      {
        clashes.push_back(std::make_pair(i, j));
      }
    }
  }
}
// _ugrid_triintersect_end

bool expensiveTest(ISpacePt& query, Triangle2DType& tri)
{
  // This is supposed to be really expensive.  Thankfully it isn't quite.
  return tri.contains(query);
}

void showImplicitGrid()
{
  // _igrid_build_start
  // here are the triangles.
  std::vector<Triangle2DType> tris;
  makeTriangles(tris);

  // Set up the ImplicitGrid: ten bins on an axis
  IGridCell res(10);
  // establish the domain of the ImplicitGrid.
  IBBox bbox(ISpacePt::zero(), ISpacePt::ones());
  // room for one hundred elements in the index
  const int numElts = static_cast<int>(tris.size());
  IGridT grid(bbox, &res, numElts);

  // load the bounding box of each triangle, along with its index, into the ImplicitGrid.
  for(int i = 0; i < numElts; ++i)
  {
    grid.insert(findBbox(tris[i]), i);
  }
  // _igrid_build_end

  // _igrid_query_start
  // Here is our query point
  ISpacePt qpt {0.63, 0.42};

  // Which triangles might it intersect?
  IBitsetType candidates = grid.getCandidates(qpt);
  int totalTrue = 0;

  // Iterate over the bitset and test the candidates expensively.
  IBitsetIndexType index = candidates.find_first();
  while(index != IBitsetType::npos)
  {
    if(expensiveTest(qpt, tris[index]))
    {
      totalTrue += 1;
    }
    index = candidates.find_next(index);
  }
  // _igrid_query_end

  // Report on intersection tests
  std::cout << "----- showImplicitGrid -----" << std::endl;
  std::cout << numElts << " total triangles, " << candidates.count()
            << " triangles expensively tested against query point, " << totalTrue << " found true."
            << std::endl;
}

void makeTriangles(std::vector<TriangleType>& tris)
{
  PointType p[8] = {PointType {0.3, 0.93, 0.03},    // 0
                    PointType {0.1, 0.85, 0.01},    // 1
                    PointType {0.3, 0.78, 0.03},    // 2
                    PointType {0.18, 0.36, 0.018},  // 3
                    PointType {0.8, 0.58, 0.08},    // 4
                    PointType {0.6, 0.5, 0.06},     // 5
                    PointType {0.55, 0.42, 0.055},  // 6
                    PointType {0.61, 0.1, 0.061}};  // 7

  tris.emplace_back(p[0], p[1], p[2]);
  tris.emplace_back(p[2], p[1], p[3]);
  tris.emplace_back(p[2], p[3], p[6]);
  tris.emplace_back(p[6], p[3], p[7]);
  tris.emplace_back(p[4], p[2], p[6]);
  tris.emplace_back(p[4], p[5], p[7]);
}

void makeTriangles(std::vector<Triangle2DType>& tris)
{
  Point2DType p[8] = {Point2DType {0.3, 0.93},   // 0
                      Point2DType {0.1, 0.85},   // 1
                      Point2DType {0.3, 0.78},   // 2
                      Point2DType {0.18, 0.36},  // 3
                      Point2DType {0.8, 0.58},   // 4
                      Point2DType {0.6, 0.5},    // 5
                      Point2DType {0.55, 0.42},  // 6
                      Point2DType {0.61, 0.1}};  // 7

  tris.emplace_back(p[0], p[1], p[2]);
  tris.emplace_back(p[2], p[1], p[3]);
  tris.emplace_back(p[2], p[3], p[6]);
  tris.emplace_back(p[6], p[3], p[7]);
  tris.emplace_back(p[4], p[2], p[6]);
  tris.emplace_back(p[4], p[5], p[7]);
}

void printPairs(std::string title, std::vector<std::pair<int, int>>& clashes)
{
  int ccount = static_cast<int>(clashes.size());
  std::cout << ccount << title << std::endl;
  for(int i = 0; i < ccount; ++i)
  {
    std::cout << clashes[i].first << "   " << clashes[i].second << std::endl;
  }
}

void driveUniformGrid()
{
  std::vector<TriangleType> tris;
  makeTriangles(tris);

  std::vector<std::pair<int, int>> naiveclashes;
  findTriIntersectionsNaively(tris, naiveclashes);

  UniformGridType* ugrid = buildUniformGrid(tris);
  std::vector<std::pair<int, int>> accelclashes;
  findTriIntersectionsAccel(tris, ugrid, accelclashes);
  delete ugrid;

  std::cout << "----- driveUniformGrid -----" << std::endl;
  printPairs(" clashes found by naive algorithm:", naiveclashes);
  printPairs(" clashes found by accelerated algorithm:", accelclashes);
  std::cout << std::endl;
}

// _bvh_build_start
BoundingBox2DType findBbox(Triangle2DType& tri);

BVH2DType* buildBVHTree(std::vector<Triangle2DType>& tris)
{
  // Create a BVH
  BVH2DType* tree = new BVH2DType;

  std::vector<BoundingBox2DType> bboxes(tris.size());

  // Get bounding boxes of each element
  for(size_t i = 0; i < tris.size(); ++i)
  {
    bboxes[i] = findBbox(tris[i]);
  }

  // Build bounding volume hierarchy from bounding boxes
  tree->initialize(bboxes.data(), bboxes.size());

  return tree;
}
// _bvh_build_end

// _bvh_candidate_start
void findCandidateBVHTreeBins(BVH2DType* tree, Point2DType ppoint, std::vector<int>& candidates)
{
  axom::IndexType offsets;
  axom::IndexType counts;
  // Get the candidates for a given probe point:
  // BVH::findPoints takes an array of points, and allocates and fills an array
  // for all the candidate intersections with the points in a packed manner.
  axom::Array<axom::IndexType> candidatesArray;
  tree->findPoints(axom::ArrayView<axom::IndexType>(&offsets, 1),
                   axom::ArrayView<axom::IndexType>(&counts, 1),
                   candidatesArray,
                   1,
                   &ppoint);

  // Since we are only querying one point, offsets == 0 and
  // len(candidatesPtr) == counts
  candidates = std::vector<int>(candidatesArray.data(), candidatesArray.data() + counts);
}
// _bvh_candidate_end

// _bvh_cand_int_start
void findIntersectionsWithCandidates(std::vector<Triangle2DType>& tris,
                                     std::vector<int>& candidates,
                                     Point2DType ppoint,
                                     std::vector<int>& intersections)
{
  // Test if ppoint lands in any of its neighbor triangles.
  int csize = static_cast<int>(candidates.size());
  for(int i = 0; i < csize; ++i)
  {
    Triangle2DType& t = tris[candidates[i]];
    if(t.contains(ppoint))
    {
      intersections.push_back(candidates[i]);
    }
  }
}
// _bvh_cand_int_end

void makeTreeTriangles(std::vector<Triangle2DType>& tris)
{
  Point2DType p[19] = {Point2DType {.13, .88},   //  0
                       Point2DType {.26, .87},   //  1
                       Point2DType {.11, .77},   //  2
                       Point2DType {.18, .78},   //  3
                       Point2DType {.13, .74},   //  4
                       Point2DType {.37, .75},   //  5
                       Point2DType {.12, .61},   //  6
                       Point2DType {.25, .51},   //  7
                       Point2DType {.11, .44},   //  8
                       Point2DType {.26, .40},   //  9
                       Point2DType {.12, .25},   // 10
                       Point2DType {.85, .38},   // 11
                       Point2DType {.94, .37},   // 12
                       Point2DType {.84, .26},   // 13
                       Point2DType {.92, .27},   // 14
                       Point2DType {.96, .28},   // 15
                       Point2DType {.84, .16},   // 16
                       Point2DType {.92, .16},   // 17
                       Point2DType {.93, .09}};  // 18

  tris.emplace_back(p[1], p[0], p[3]);
  tris.emplace_back(p[0], p[2], p[3]);
  tris.emplace_back(p[1], p[3], p[5]);
  tris.emplace_back(p[3], p[2], p[4]);
  tris.emplace_back(p[3], p[4], p[5]);
  tris.emplace_back(p[2], p[6], p[4]);
  tris.emplace_back(p[5], p[4], p[6]);
  tris.emplace_back(p[5], p[6], p[7]);
  tris.emplace_back(p[7], p[6], p[8]);
  tris.emplace_back(p[7], p[8], p[9]);
  tris.emplace_back(p[9], p[8], p[10]);
  tris.emplace_back(p[11], p[13], p[14]);
  tris.emplace_back(p[12], p[11], p[14]);
  tris.emplace_back(p[12], p[14], p[15]);
  tris.emplace_back(p[14], p[13], p[16]);
  tris.emplace_back(p[14], p[16], p[17]);
  tris.emplace_back(p[16], p[18], p[17]);
}

void driveBVHTree()
{
  std::vector<Triangle2DType> tris;
  makeTreeTriangles(tris);
  Point2DType ppoint {0.45, 0.25};
  std::vector<int> intersections, candidates;

  BVH2DType* tree = buildBVHTree(tris);
  findCandidateBVHTreeBins(tree, ppoint, candidates);
  findIntersectionsWithCandidates(tris, candidates, ppoint, intersections);
  tree->writeVtkFile("BVHTree.out.vtk");

  std::cout << "----- driveBVHTree -----" << std::endl;
  std::cout << "Point " << ppoint << " hit the following " << intersections.size()
            << " triangles (0 expected):" << std::endl;
  for(size_t i = 0; i < intersections.size(); ++i)
  {
    std::cout << intersections[i] << std::endl;
  }

  delete tree;
}

void driveOctree()
{
  // _octree_start
  OctBBox bb(OctSpacePt(10), OctSpacePt(20));

  // Generate a point within the bounding box
  double alpha = 2. / 3.;
  OctSpacePt queryPt = OctSpacePt::lerp(bb.getMin(), bb.getMax(), alpha);

  // Instantiate the Octree
  OctreeType octree(bb);

  // Find the block containing the query point
  OctBlockIndex leafBlock = octree.findLeafBlock(queryPt);
  // and the bounding box of the block.
  OctBBox leafBB = octree.blockBoundingBox(leafBlock);

  for(int i = 0; i < octree.maxInternalLevel(); ++i)
  {
    // SpatialOctree allows a code to refine (subdivide) a block
    octree.refineLeaf(leafBlock);
    // and locate the (new) child block containing the query point.
    leafBlock = octree.findLeafBlock(queryPt);
  }
  // _octree_end

  AXOM_UNUSED_VAR(leafBB);
}

int main(int argc, char** argv)
{
  // Deal with unused variables
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  demoRectangularLattice();
  demoMorton();
  driveUniformGrid();
  showImplicitGrid();
  driveBVHTree();
  driveOctree();

  return 0;
}
