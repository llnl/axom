// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file DelaunayImpl.hpp
 *
 * \brief Defines the main incremental insertion and mesh-update routines for `quest::Delaunay`.
 */

#ifndef AXOM_QUEST_DETAIL_DELAUNAY_IMPL_HPP_
#define AXOM_QUEST_DETAIL_DELAUNAY_IMPL_HPP_

namespace axom
{
namespace quest
{

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE bool Delaunay<DIM>::isSearchableElement(IndexType element_idx) const
{
  return m_mesh.isValidElement(element_idx);
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE typename Delaunay<DIM>::CircumsphereEval
Delaunay<DIM>::evaluateCircumsphereOnMesh(const IAMeshType& mesh, IndexType element_idx)
{
  using axom::numerics::determinant;

  const auto verts = mesh.boundaryVertices(element_idx);

  if constexpr(DIM == 2)
  {
    const PointType& p0 = mesh.getVertexPosition(verts[0]);
    const PointType& p1 = mesh.getVertexPosition(verts[1]);
    const PointType& p2 = mesh.getVertexPosition(verts[2]);

    const double vx0 = p1[0] - p0[0];
    const double vx1 = p2[0] - p0[0];
    const double vy0 = p1[1] - p0[1];
    const double vy1 = p2[1] - p0[1];

    const double sq0 = vx0 * vx0 + vy0 * vy0;
    const double sq1 = vx1 * vx1 + vy1 * vy1;

    const double a = determinant(vx0, vx1, vy0, vy1);
    const double eps = (a >= 0.) ? primal::PRIMAL_TINY : -primal::PRIMAL_TINY;
    const double ood = 1. / (2. * a + eps);

    const double center_offset_x = determinant(sq0, sq1, vy0, vy1) * ood;
    const double center_offset_y = -determinant(sq0, sq1, vx0, vx1) * ood;

    return CircumsphereEval(p0, VectorType {center_offset_x, center_offset_y});
  }
  else
  {
    const PointType& p0 = mesh.getVertexPosition(verts[0]);
    const PointType& p1 = mesh.getVertexPosition(verts[1]);
    const PointType& p2 = mesh.getVertexPosition(verts[2]);
    const PointType& p3 = mesh.getVertexPosition(verts[3]);

    const double vx0 = p1[0] - p0[0];
    const double vx1 = p2[0] - p0[0];
    const double vx2 = p3[0] - p0[0];
    const double vy0 = p1[1] - p0[1];
    const double vy1 = p2[1] - p0[1];
    const double vy2 = p3[1] - p0[1];
    const double vz0 = p1[2] - p0[2];
    const double vz1 = p2[2] - p0[2];
    const double vz2 = p3[2] - p0[2];

    const double sq0 = vx0 * vx0 + vy0 * vy0 + vz0 * vz0;
    const double sq1 = vx1 * vx1 + vy1 * vy1 + vz1 * vz1;
    const double sq2 = vx2 * vx2 + vy2 * vy2 + vz2 * vz2;

    const double a = determinant(vx0, vx1, vx2, vy0, vy1, vy2, vz0, vz1, vz2);
    const double eps = (a >= 0.) ? primal::PRIMAL_TINY : -primal::PRIMAL_TINY;
    const double ood = 1. / (2. * a + eps);

    const double center_offset_x = determinant(sq0, sq1, sq2, vy0, vy1, vy2, vz0, vz1, vz2) * ood;
    const double center_offset_y = determinant(sq0, sq1, sq2, vz0, vz1, vz2, vx0, vx1, vx2) * ood;
    const double center_offset_z = determinant(sq0, sq1, sq2, vx0, vx1, vx2, vy0, vy1, vy2) * ood;

    return CircumsphereEval(p0, VectorType {center_offset_x, center_offset_y, center_offset_z});
  }
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE double Delaunay<DIM>::sphereSquaredDistanceTolerance(
  const CircumsphereEval& sphere,
  const PointType& x,
  double distance_sq)
{
  double scale = 1.;

  for(int dim = 0; dim < DIM; ++dim)
  {
    scale = axom::utilities::max(scale, axom::utilities::abs(sphere.center[dim]));
    scale = axom::utilities::max(scale, axom::utilities::abs(x[dim]));
  }

  const double local_span =
    axom::utilities::max(1., 2. * std::sqrt(axom::utilities::max(distance_sq, sphere.radius_sq)));

  return 256. * std::numeric_limits<double>::epsilon() * scale * local_span;
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE int Delaunay<DIM>::inSphereOrientationOnMesh(const IAMeshType& mesh,
                                                                              const PointType& q,
                                                                              IndexType element_idx)
{
  const auto sphere = evaluateCircumsphereOnMesh(mesh, element_idx);
  const double distance_sq = primal::squared_distance(sphere.center, q);

  const double delta_sq = distance_sq - sphere.radius_sq;
  const double tol = sphereSquaredDistanceTolerance(sphere, q, distance_sq);
  if(axom::utilities::abs(delta_sq) <= tol)
  {
    return primal::ON_BOUNDARY;
  }

  return delta_sq < 0. ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE bool Delaunay<DIM>::isPointInSphereOnMesh(const IAMeshType& mesh,
                                                                           const PointType& q,
                                                                           IndexType element_idx,
                                                                           bool includeBoundary)
{
  const int res = inSphereOrientationOnMesh(mesh, q, element_idx);
  return includeBoundary ? (res != primal::ON_POSITIVE_SIDE) : (res == primal::ON_NEGATIVE_SIDE);
}

template <int DIM>
inline bool Delaunay<DIM>::shouldCompactMesh() const
{
  constexpr int MIN_REMOVED_ELEMENTS = DIM == 2 ? 512 : 2048;
  constexpr double REMOVED_ELEMENT_FRACTION = DIM == 2 ? 0.25 : 0.35;
  return static_cast<int>(m_deleted_elements.size()) > MIN_REMOVED_ELEMENTS &&
    (static_cast<double>(m_deleted_elements.size()) >
     REMOVED_ELEMENT_FRACTION * static_cast<double>(m_mesh.elements().size()));
}

template <int DIM>
inline void Delaunay<DIM>::compactMesh()
{
  m_mesh.compact();
  m_deleted_elements.clear();
  m_element_finder.recomputeGrid(m_mesh, m_bounding_box);
  if(m_next_regrid_vertex_count > 0)
  {
    while(m_next_regrid_vertex_count <= m_mesh.vertices().size())
    {
      m_next_regrid_vertex_count *= 2;
    }
  }
}

template <int DIM>
inline void Delaunay<DIM>::writeToVTKFile(const std::string& filename)
{
  const auto CELL_TYPE = DIM == 2 ? mint::TRIANGLE : mint::TET;
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mint_mesh(DIM, CELL_TYPE);

  this->compactMesh();

  for(auto v : m_mesh.vertices().positions())
  {
    mint_mesh.appendNodes(m_mesh.getVertexPosition(v).data(), 1);
  }

  for(auto e : m_mesh.elements().positions())
  {
    mint_mesh.appendCell(&(m_mesh.boundaryVertices(e)[0]), CELL_TYPE);
  }

  mint::write_vtk(&mint_mesh, filename);
}

template <int DIM>
inline void Delaunay<DIM>::removeBoundary()
{
  if(m_has_boundary)
  {
    const int num_boundary_pts = 1 << DIM;

    for(auto e : m_mesh.elements().positions())
    {
      if(m_mesh.isValidElement(e))
      {
        const auto verts = m_mesh.boundaryVertices(e);
        bool touches_boundary = false;
        for(int i = 0; i < VERT_PER_ELEMENT; ++i)
        {
          touches_boundary |= (verts[i] >= 0 && verts[i] < num_boundary_pts);
        }

        if(touches_boundary)
        {
          m_mesh.removeElement(e);
        }
      }
    }

    for(int v = 0; v < num_boundary_pts; ++v)
    {
      m_mesh.removeVertex(v);
    }

    for(auto e : m_mesh.elements().positions())
    {
      if(!m_mesh.isValidElement(e))
      {
        continue;
      }

      const auto verts = m_mesh.boundaryVertices(e);
      bool has_invalid_vertex = false;
      for(int i = 0; i < VERT_PER_ELEMENT; ++i)
      {
        has_invalid_vertex |= !m_mesh.isValidVertex(verts[i]);
      }

      if(has_invalid_vertex)
      {
        m_mesh.removeElement(e);
      }
    }

    this->compactMesh();
    m_has_boundary = false;
  }
}

template <int DIM>
inline void Delaunay<DIM>::initializeBoundary(const BoundingBox& bb)
{
  std::vector<DataType> points;
  IndexArray elem;

  generateInitialMesh(points, elem, bb);

  m_mesh = IAMeshType(points, elem);
  m_element_finder.recomputeGrid(m_mesh, bb);
  m_next_regrid_vertex_count = 1024;
  m_walk_visited = slam::BitSet(static_cast<int>(m_mesh.elements().size()));
  m_total_removed_elements = 0;
  m_max_removed_elements = 0;
  m_num_insertions = 0;
  m_num_walk_calls = 0;
  m_num_walk_found = 0;
  m_num_walk_outside = 0;
  m_num_walk_failed = 0;
  m_total_walk_steps = 0;
  m_max_walk_steps = 0;

  m_candidate_elements_scratch.clear();
  m_candidate_elements_scratch.reserve(16);
  m_walked_elements_scratch.clear();
  m_walked_elements_scratch.reserve(256);
  m_walk_local_elements_scratch.clear();
  m_walk_local_elements_scratch.reserve(256);
  m_initial_vertices_scratch.clear();
  m_initial_vertices_scratch.reserve(8);
  m_deleted_elements.clear();
  m_deleted_elements.reserve(DIM == 2 ? 128 : 512);

  if(!m_insertion_helper)
  {
    m_insertion_helper = std::make_unique<InsertionHelper>(m_mesh);
  }

  m_bounding_box = bb;
  m_has_boundary = true;
}

template <int DIM>
inline void Delaunay<DIM>::insertPoint(const PointType& new_pt)
{
  SLIC_ASSERT_MSG(m_has_boundary, "Error: Need a predefined boundary box prior to adding points.");
  SLIC_ASSERT_MSG(m_insertion_helper != nullptr,
                  "Error: Insertion helper was not initialized. "
                  "Delaunay::initializeBoundary() needs to be called first.");
  SLIC_ASSERT_MSG(m_bounding_box.contains(new_pt),
                  "Error: new point is outside of the boundary box.");

  IndexType element_i = findContainingElement(new_pt);

  if(element_i == INVALID_INDEX)
  {
    SLIC_WARNING(
      fmt::format("Could not insert point {} into Delaunay triangulation: "
                  "Element containing that point was not found",
                  new_pt));
    return;
  }

  const BaryCoordType bary_coord = getBaryCoords(element_i, new_pt);
  const IndexArray seed_elements = getSeedElements(element_i, bary_coord);
  validateInsertionSeed(element_i, new_pt, bary_coord, seed_elements);

  auto& insertionHelper = *m_insertion_helper;
  insertionHelper.reset();
  insertionHelper.containing_element = element_i;
  insertionHelper.containing_bary = bary_coord;
  insertionHelper.seed_elements_debug = seed_elements;
  insertionHelper.findCavityElements(
    new_pt,
    seed_elements,
    [&](const PointType& query_pt, IndexType element_idx) {
      return isPointInSphereOnMesh(m_mesh, query_pt, element_idx, /*includeBoundary=*/true);
    });
  ++m_num_insertions;
  m_total_removed_elements += static_cast<std::uint64_t>(insertionHelper.numRemovedElements());
  m_max_removed_elements =
    axom::utilities::max(m_max_removed_elements,
                         static_cast<std::uint64_t>(insertionHelper.numRemovedElements()));
  validateCavityBoundary(insertionHelper);
  insertionHelper.createCavity(m_deleted_elements);
  IndexType new_pt_i = m_mesh.addVertex(new_pt);
  insertionHelper.delaunayBall(new_pt_i, m_deleted_elements);
  validateInsertedBall(new_pt_i, insertionHelper);
  validateInsertionResult();

  m_element_finder.updateBin(new_pt, new_pt_i);
  if(m_next_regrid_vertex_count > 0 && m_mesh.vertices().size() >= m_next_regrid_vertex_count)
  {
    m_element_finder.recomputeGrid(m_mesh, m_bounding_box);
    while(m_next_regrid_vertex_count <= m_mesh.vertices().size())
    {
      m_next_regrid_vertex_count *= 2;
    }
  }

  if(shouldCompactMesh())
  {
    this->compactMesh();
  }
}

template <int DIM>
inline void Delaunay<DIM>::generateInitialMesh(std::vector<DataType>& points,
                                               std::vector<IndexType>& elem,
                                               const BoundingBox& bb)
{
  const PointType& mins = bb.getMin();
  const PointType& maxs = bb.getMax();

  if constexpr(DIM == 2)
  {
    std::vector<DataType> pt {mins[0], mins[1], mins[0], maxs[1], maxs[0], mins[1], maxs[0], maxs[1]};
    std::vector<IndexType> el {0, 2, 1, 3, 1, 2};

    points.swap(pt);
    elem.swap(el);
  }
  else
  {
    std::vector<DataType> pt {mins[0], mins[1], mins[2], mins[0], mins[1], maxs[2],
                              mins[0], maxs[1], mins[2], mins[0], maxs[1], maxs[2],
                              maxs[0], mins[1], mins[2], maxs[0], mins[1], maxs[2],
                              maxs[0], maxs[1], mins[2], maxs[0], maxs[1], maxs[2]};

    std::vector<IndexType> el {3, 2, 4, 0, 3, 4, 1, 0, 3, 2, 6, 4,
                               3, 6, 7, 4, 3, 5, 1, 4, 3, 7, 5, 4};

    points.swap(pt);
    elem.swap(el);
  }
}

}  // namespace quest
}  // namespace axom

#endif
