// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/quest/MeshTester.hpp"

#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
inline detail::SpatialBoundingBox compute_bounds(detail::UMesh* mesh)
{
  SLIC_ASSERT(mesh != nullptr);

  detail::SpatialBoundingBox meshBB;
  detail::Point3 pt;

  const double* x = mesh->getCoordinateArray(mint::X_COORDINATE);
  const double* y = mesh->getCoordinateArray(mint::Y_COORDINATE);
  const double* z = mesh->getCoordinateArray(mint::Z_COORDINATE);

  const axom::IndexType numNodes = mesh->getNumberOfNodes();
  for(axom::IndexType i = 0; i < numNodes; ++i)
  {
    pt[0] = x[i];
    pt[1] = y[i];
    pt[2] = z[i];

    meshBB.addPoint(pt);
  }

  if(numNodes > 0)
  {
    SLIC_ASSERT(meshBB.isValid());
  }

  return meshBB;
}

inline bool areTriangleIndicesDistinct(axom::IndexType* indices)
{
  SLIC_ASSERT(indices != nullptr);

  return indices[0] != indices[1] && indices[1] != indices[2] && indices[2] != indices[0];
}

/// Find and report self-intersections and degenerate triangles in a triangle surface mesh using a Uniform Grid.
void findTriMeshIntersections(detail::UMesh* surface_mesh,
                              std::vector<std::pair<int, int>>& intersections,
                              std::vector<int>& degenerateIndices,
                              int spatialIndexResolution,
                              double intersectionThreshold)
{
  AXOM_ANNOTATE_SCOPE("quest::findTriMeshIntersections");

  detail::Triangle3 t1 {};
  detail::Triangle3 t2 {};
  SLIC_INFO("Running mesh_tester with UniformGrid index");

  const int ncells = surface_mesh->getNumberOfCells();

  // find the specified resolution.  If we're passed a number less than one, use the cube root of the number of triangles.
  if(spatialIndexResolution < 1)
  {
    spatialIndexResolution = (int)(1 + std::pow(ncells, 1 / 3.));
  }
  NumericArray<int, 3> resolutions(spatialIndexResolution);

  SLIC_INFO("Building UniformGrid index...");
  std::vector<int> nondegenerateIndices;
  nondegenerateIndices.reserve(ncells);

  axom::Array<int> triIdxs(ncells);
  axom::Array<detail::SpatialBoundingBox> triBboxes(ncells);
  for(int i = 0; i < ncells; i++)
  {
    t1 = detail::getMeshTriangle(i, surface_mesh);

    if(t1.degenerate())
    {
      degenerateIndices.push_back(i);
    }
    else
    {
      nondegenerateIndices.push_back(i);

      triIdxs[i] = i;
      triBboxes[i] = compute_bounding_box(t1);
    }
  }
  detail::UniformGrid3 ugrid(resolutions, triBboxes.view(), triIdxs.view());

  // Iterate through triangle indices *idx.
  // Check against each other triangle with index greater than the index *idx that also shares a UniformGrid bin.
  SLIC_INFO(
    axom::fmt::format(axom::utilities::locale(), "Checking mesh with a total of {:L} cells", ncells));

  std::vector<int>::iterator idx = nondegenerateIndices.begin(), ndgend = nondegenerateIndices.end();
  for(; idx != ndgend; ++idx)
  {
    // Retrieve the triangle at *idx and construct a bounding box around it
    t1 = detail::getMeshTriangle(*idx, surface_mesh);
    detail::SpatialBoundingBox triBB2 = compute_bounding_box(t1);

    // Get a list of all triangles in bins this triangle will touch, whose indices are greater than this triangle's index
    std::vector<int> neighborTriangles;
    const std::vector<int> binsToCheck = ugrid.getBinsForBbox(triBB2);
    size_t checkcount = binsToCheck.size();
    for(size_t curbin = 0; curbin < checkcount; ++curbin)
    {
      axom::ArrayView<int> ntlist = ugrid.getBinContents(binsToCheck[curbin]);
      for(const int nbr : ntlist)
      {
        if(nbr > *idx)
        {
          neighborTriangles.push_back(nbr);
        }
      }
    }

    std::sort(neighborTriangles.begin(), neighborTriangles.end());
    std::vector<int>::iterator nend = std::unique(neighborTriangles.begin(), neighborTriangles.end());
    std::vector<int>::iterator nit = neighborTriangles.begin();

    // test any remaining neighbor tris for intersection
    while(nit != nend)
    {
      t2 = detail::getMeshTriangle(*nit, surface_mesh);
      if(primal::intersect(t1, t2, false, intersectionThreshold))
      {
        intersections.push_back(std::make_pair(*idx, *nit));
      }
      ++nit;
    }
  }
}

/// Check a surface mesh for holes using its face relation.
WatertightStatus isSurfaceMeshWatertight(detail::UMesh* surface_mesh)
{
  AXOM_ANNOTATE_SCOPE("quest::isSurfaceMeshWatertight");

  // Make sure the mesh is reasonable
  SLIC_ASSERT_MSG(surface_mesh != nullptr,
                  "surface_mesh must be a valid pointer to a triangle mesh");

  // Calculate the face relations---this can take awhile
  bool success = surface_mesh->initializeFaceConnectivity();

  if(!success)
  {
    return WatertightStatus::CHECK_FAILED;
  }

  WatertightStatus retval = WatertightStatus::WATERTIGHT;

  constexpr int ON_BOUNDARY = 1;
  constexpr int INTERNAL = 0;

  // Create fields to store boundary
  int* bndry_face = surface_mesh->createField<int>("bndry_face", mint::FACE_CENTERED);
  int* boundary = surface_mesh->createField<int>("boundary", mint::CELL_CENTERED);

  // Mark boundary faces
  const IndexType numFaces = surface_mesh->getNumberOfFaces();
  for(IndexType iface = 0; iface < numFaces; ++iface)
  {
    IndexType c1, c2;
    surface_mesh->getFaceCellIDs(iface, c1, c2);
    SLIC_ASSERT(c1 != static_cast<IndexType>(mint::UNDEFINED_CELL));

    if(c2 == static_cast<IndexType>(mint::UNDEFINED_CELL))
    {
      bndry_face[iface] = ON_BOUNDARY;
      retval = WatertightStatus::NOT_WATERTIGHT;
    }
    else
    {
      bndry_face[iface] = INTERNAL;
    }
  }

  if(retval == WatertightStatus::WATERTIGHT)
  {
    /* short-circuit */
    const IndexType numCells = surface_mesh->getNumberOfCells();
    std::memset(boundary, INTERNAL, sizeof(int) * numCells);
    return retval;
  }

  // Mark boundary cells
  const IndexType numCells = surface_mesh->getNumberOfCells();
  for(IndexType icell = 0; icell < numCells; ++icell)
  {
    // NOTE: this check currently assumes triangles
    SLIC_ASSERT(surface_mesh->getNumberOfCellFaces(icell) == 3);
    const IndexType* faceids = surface_mesh->getCellFaceIDs(icell);

    if((bndry_face[faceids[0]] == ON_BOUNDARY) || (bndry_face[faceids[1]] == ON_BOUNDARY) ||
       (bndry_face[faceids[2]] == ON_BOUNDARY))
    {
      boundary[icell] = ON_BOUNDARY;
    }
    else
    {
      boundary[icell] = INTERNAL;
    }
  }

  return retval;
}

/// Weld vertices of a triangle mesh that are closer than \a eps
void weldTriMeshVertices(detail::UMesh** surface_mesh, double eps)
{
  AXOM_ANNOTATE_SCOPE("quest::weldTriMeshVertices");

  // Note: Use 64-bit index to accomodate small values of epsilon
  using IdxType = std::int64_t;
  using Lattice3 = spin::RectangularLattice<3, double, IdxType>;
  using GridCell = Lattice3::GridCell;

  // Define a lambda for hashing points
  // implementation of hash combiner is from boost's hash_combine()
  auto point_hash = [](const GridCell& pt) {
    auto seed = std::hash<IdxType> {}(pt[0]);
    for(int i = 1; i < GridCell::DIMENSION; ++i)
    {
      seed ^= std::hash<IdxType> {}(pt[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  };

  using GridCellToIndexMap = std::unordered_map<GridCell, IdxType, decltype(point_hash)>;

  /// Implementation notes:
  ///
  /// This function welds vertices in the triangle mesh by
  /// quantizing them to an integer lattice with spacing \a eps.
  ///
  /// To avoid searching neighbor elements on the grid, we run the
  /// algorithm twice. Once on the original lattice, and once on a
  /// lattice that is shifted by half the grid spacing.
  ///
  /// Due to running this algorithm twice, it is possible in extreme cases
  /// for vertices that are at distances up to (1.5 * eps)
  /// (under the max norm) to be identified in this process.
  /// This can correspond to distances of at most 1.5 * eps * \sqrt(2)
  /// in the Euclidean norm.

  SLIC_ASSERT_MSG(eps > 0., "Epsilon must be greater than 0. Passed in value was " << eps);
  SLIC_ASSERT_MSG(surface_mesh != nullptr && *surface_mesh != nullptr,
                  "surface_mesh must be a valid pointer to a pointer to a triangle mesh");

  int const DIM = 3;
  detail::UMesh* oldMesh = *surface_mesh;

  detail::SpatialBoundingBox meshBB = compute_bounds(oldMesh).expand(eps);

  // Run the algorithm twice -- on the original grid and a translated grid
  for(const auto& offset : {detail::Point3(0.), detail::Point3(eps / 2.)})
  {
    // We will build up a new triangle mesh with the welded indices
    detail::UMesh* newMesh = new detail::UMesh(DIM, mint::TRIANGLE);

    // Set up the lattice for quantizing points to an integer lattice
    detail::Point3 origin(meshBB.getMin().array() - offset.array());
    Lattice3 lattice(origin, Lattice3::SpaceVector(detail::Point3(eps)));

    // First, find unique indices for the welded vertices
    const int numVerts = oldMesh->getNumberOfNodes();
    int uniqueVertCount = 0;

    // A map from GridCells to the new vertex indices
    GridCellToIndexMap vertexIndexMap(numVerts, point_hash);

    std::vector<int> vertex_remap;  // stores the new vertex indices
    vertex_remap.resize(numVerts);  // for each old vertex

    detail::Point3 vert;
    const double* x = oldMesh->getCoordinateArray(mint::X_COORDINATE);
    const double* y = oldMesh->getCoordinateArray(mint::Y_COORDINATE);
    const double* z = oldMesh->getCoordinateArray(mint::Z_COORDINATE);

    for(int i = 0; i < numVerts; ++i)
    {
      // get the vertex from the mesh
      vert[0] = x[i];
      vert[1] = y[i];
      vert[2] = z[i];

      // find the new vertex index; if not present, insert vertex into new mesh
      auto res = vertexIndexMap.insert(std::make_pair(lattice.gridCell(vert), uniqueVertCount));
      if(res.second == true)
      {
        uniqueVertCount++;
        newMesh->appendNodes(vert.data());
      }
      vertex_remap[i] = static_cast<int>(res.first->second);
    }

    // Next, add triangles into the new mesh using the unique vertex indices
    const int NUM_TRI_VERTS = 3;
    axom::IndexType triInds[NUM_TRI_VERTS];
    const axom::IndexType numTris = oldMesh->getNumberOfCells();
    for(axom::IndexType i = 0; i < numTris; ++i)
    {
      memcpy(triInds, oldMesh->getCellNodeIDs(i), NUM_TRI_VERTS * sizeof(axom::IndexType));

      for(int d = 0; d < NUM_TRI_VERTS; ++d)
      {
        triInds[d] = vertex_remap[triInds[d]];
      }

      // Degeneracy check -- vertices need to be distinct
      if(areTriangleIndicesDistinct(triInds))
      {
        newMesh->appendCell(triInds, mint::TRIANGLE);
      }
    }

    // Finally, delete old mesh and swap pointers
    delete oldMesh;
    oldMesh = newMesh;
  }

  // Update the original mesh pointer to the new mesh
  *surface_mesh = oldMesh;
}

}  // end namespace quest
}  // end namespace axom
