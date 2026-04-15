// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_DETAIL_DELAUNAY_POINT_LOCATION_HPP_
#define AXOM_QUEST_DETAIL_DELAUNAY_POINT_LOCATION_HPP_

namespace axom
{
namespace quest
{

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE typename Delaunay<DIM>::BaryCoordType Delaunay<DIM>::getBaryCoords(
  IndexType element_idx,
  const PointType& query_pt) const
{
  const auto verts = m_mesh.boundaryVertices(element_idx);

  if constexpr(DIM == 2)
  {
    const ElementType tri(m_mesh.getVertexPosition(verts[0]),
                          m_mesh.getVertexPosition(verts[1]),
                          m_mesh.getVertexPosition(verts[2]));

    return tri.physToBarycentric(query_pt);
  }
  else
  {
    const ElementType tet(m_mesh.getVertexPosition(verts[0]),
                          m_mesh.getVertexPosition(verts[1]),
                          m_mesh.getVertexPosition(verts[2]),
                          m_mesh.getVertexPosition(verts[3]));

    return tet.physToBarycentric(query_pt);
  }
}

template <int DIM>
inline typename Delaunay<DIM>::BaryCoordType Delaunay<DIM>::getRawBarycentricDeterminants(
  IndexType element_idx,
  const PointType& query_pt) const
{
  const auto verts = m_mesh.boundaryVertices(element_idx);

  if constexpr(DIM == 2)
  {
    const ElementType tri(m_mesh.getVertexPosition(verts[0]),
                          m_mesh.getVertexPosition(verts[1]),
                          m_mesh.getVertexPosition(verts[2]));
    return tri.physToBarycentric(query_pt, /*skipNormalization=*/true);
  }
  else
  {
    const ElementType tet(m_mesh.getVertexPosition(verts[0]),
                          m_mesh.getVertexPosition(verts[1]),
                          m_mesh.getVertexPosition(verts[2]),
                          m_mesh.getVertexPosition(verts[3]));
    return tet.physToBarycentric(query_pt, /*skipNormalization=*/true);
  }
}

template <int DIM>
inline double Delaunay<DIM>::rawBarycentricDeterminantTolerance(IndexType element_idx,
                                                                const PointType& query_pt) const
{
  const auto verts = m_mesh.boundaryVertices(element_idx);

  double scale = 1.;
  for(int i = 0; i < VERT_PER_ELEMENT; ++i)
  {
    const auto diff = m_mesh.getVertexPosition(verts[i]) - query_pt;
    scale = axom::utilities::max(scale, diff.norm());
  }

  const double k = 64.;
  if constexpr(DIM == 2)
  {
    return k * std::numeric_limits<double>::epsilon() * scale * scale;
  }
  else
  {
    return k * std::numeric_limits<double>::epsilon() * scale * scale * scale;
  }
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE bool Delaunay<DIM>::isPointInsideForLocation(
  IndexType element_idx,
  const PointType& query_pt,
  const BaryCoordType& bary_coord,
  ModularFaceIndex* exit_face) const
{
  ModularFaceIndex min_face(bary_coord.array().argMin());

  for(int i = 0; i < VERT_PER_ELEMENT; ++i)
  {
    if(bary_coord[i] < -BARY_EPS)
    {
      if(exit_face != nullptr)
      {
        *exit_face = min_face;
      }
      return false;
    }
  }

  int first_determinant_negative = -1;
  if constexpr(DIM == 3)
  {
    bool has_near_zero = false;
    for(int i = 0; i < VERT_PER_ELEMENT; ++i)
    {
      has_near_zero |= axom::utilities::abs(bary_coord[i]) <= BARY_EPS;
    }

    if(has_near_zero)
    {
      const BaryCoordType raw = getRawBarycentricDeterminants(element_idx, query_pt);
      const double tol = rawBarycentricDeterminantTolerance(element_idx, query_pt);
      for(int i = 0; i < VERT_PER_ELEMENT; ++i)
      {
        if(axom::utilities::abs(bary_coord[i]) <= BARY_EPS && signWithTolerance(raw[i], tol) < 0)
        {
          first_determinant_negative = i;
          break;
        }
      }
    }
  }

  if(first_determinant_negative >= 0)
  {
    if(exit_face != nullptr)
    {
      *exit_face = ModularFaceIndex(first_determinant_negative);
    }
    return false;
  }

  return true;
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE typename Delaunay<DIM>::PointLocationResult
Delaunay<DIM>::walkToContainingElement(const PointType& query_pt,
                                       IndexType start_element,
                                       std::vector<IndexType>* visited_elements_out) const
{
  constexpr IndexType invalid_element = IAMeshType::INVALID_ELEMENT_INDEX;

  if(!m_mesh.isValidElement(start_element))
  {
    return {};
  }

  static constexpr int MAX_WALK_STEPS = 256;
  std::vector<IndexType>& visited_elements =
    visited_elements_out != nullptr ? *visited_elements_out : m_walk_local_elements_scratch;
  visited_elements.clear();
  if(static_cast<int>(visited_elements.capacity()) < MAX_WALK_STEPS)
  {
    visited_elements.reserve(MAX_WALK_STEPS);
  }
  IndexType element_i = start_element;

  if(m_walk_visited.size() < static_cast<int>(m_mesh.elements().size()))
  {
    m_walk_visited = slam::BitSet(static_cast<int>(m_mesh.elements().size()));
  }

  auto clearVisitedBits = [&]() {
    for(const IndexType visited : visited_elements)
    {
      m_walk_visited.clear(static_cast<int>(visited));
    }
  };

  int step_count = 0;

  auto recordWalk = [&](PointLocationStatus status) {
    if(!m_collect_location_stats)
    {
      return;
    }

    ++m_num_walk_calls;
    m_total_walk_steps += static_cast<std::uint64_t>(step_count);
    m_max_walk_steps = axom::utilities::max(m_max_walk_steps, static_cast<std::uint64_t>(step_count));
    switch(status)
    {
    case PointLocationStatus::Found:
      ++m_num_walk_found;
      break;
    case PointLocationStatus::Outside:
      ++m_num_walk_outside;
      break;
    default:
      ++m_num_walk_failed;
      break;
    }
  };

  while(1)
  {
    ++step_count;
    if(m_walk_visited.test(static_cast<int>(element_i)))
    {
      recordWalk(PointLocationStatus::Failed);
      clearVisitedBits();
      return {};
    }
    m_walk_visited.set(static_cast<int>(element_i));
    visited_elements.push_back(element_i);

    const BaryCoordType bary_coord = getBaryCoords(element_i, query_pt);
    ModularFaceIndex modular_idx(0);
    if(isPointInsideForLocation(element_i, query_pt, bary_coord, &modular_idx))
    {
      recordWalk(PointLocationStatus::Found);
      clearVisitedBits();
      return {element_i, PointLocationStatus::Found};
    }

    if(static_cast<int>(visited_elements.size()) >= MAX_WALK_STEPS)
    {
      recordWalk(PointLocationStatus::Failed);
      clearVisitedBits();
      return {};
    }

    const IndexType next_element = m_mesh.adjacentElements(element_i)[modular_idx + 1];
    if(next_element == invalid_element)
    {
      recordWalk(PointLocationStatus::Outside);
      clearVisitedBits();
      return {INVALID_INDEX, PointLocationStatus::Outside};
    }

    SLIC_ASSERT(m_mesh.isValidElement(next_element));
    element_i = next_element;
  }
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE void Delaunay<DIM>::appendCandidateElement(
  std::vector<IndexType>& candidate_elements,
  IndexType vertex_i) const
{
  const IndexType element_i = m_mesh.coboundaryElement(vertex_i);
  if(isSearchableElement(element_i) &&
     std::find(candidate_elements.begin(), candidate_elements.end(), element_i) ==
       candidate_elements.end())
  {
    candidate_elements.push_back(element_i);
  }
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE void Delaunay<DIM>::appendCandidateElementsFromVertices(
  std::vector<IndexType>& candidate_elements,
  const std::vector<IndexType>& candidate_vertices) const
{
  for(const auto vertex_i : candidate_vertices)
  {
    appendCandidateElement(candidate_elements, vertex_i);
  }
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE void Delaunay<DIM>::getInitialCandidateElements(
  const PointType& query_pt,
  std::vector<IndexType>& candidate_elements) const
{
  candidate_elements.clear();
  candidate_elements.reserve(16);

  m_initial_vertices_scratch.clear();
  m_element_finder.getNearbyVertices(m_mesh,
                                     query_pt,
                                     m_initial_vertices_scratch,
                                     /*search_radius=*/1,
                                     /*max_candidates=*/8);
  appendCandidateElementsFromVertices(candidate_elements, m_initial_vertices_scratch);

  if(candidate_elements.empty())
  {
    if(m_collect_location_stats)
    {
      ++m_num_empty_seed_fallbacks;
    }
    for(auto elem : m_mesh.elements().positions())
    {
      if(isSearchableElement(elem))
      {
        candidate_elements.push_back(elem);
        break;
      }
    }
  }
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE typename Delaunay<DIM>::PointLocationResult
Delaunay<DIM>::walkCandidateElements(const PointType& query_pt,
                                     const std::vector<IndexType>& candidate_elements,
                                     std::size_t start_idx,
                                     std::vector<IndexType>* walked_elements) const
{
  for(std::size_t idx = start_idx; idx < candidate_elements.size(); ++idx)
  {
    std::vector<IndexType>* visited_elements =
      (walked_elements != nullptr && idx == start_idx) ? walked_elements : nullptr;
    PointLocationResult walk_result =
      walkToContainingElement(query_pt, candidate_elements[idx], visited_elements);
    if(walk_result.status != PointLocationStatus::Failed)
    {
      return walk_result;
    }
  }

  return {};
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE typename Delaunay<DIM>::PointLocationResult
Delaunay<DIM>::findContainingElementWithQueryFallbacks(const PointType& query_pt,
                                                       std::vector<IndexType>& candidate_elements,
                                                       const std::vector<IndexType>& walked_elements) const
{
  const IndexType walk_region_elem = findContainingElementFromNeighbors(query_pt, walked_elements);
  if(walk_region_elem != INVALID_INDEX)
  {
    return {walk_region_elem, PointLocationStatus::Found};
  }

  m_fallback_vertices_scratch.clear();
  m_element_finder.getNearbyVertices(m_mesh,
                                     query_pt,
                                     m_fallback_vertices_scratch,
                                     QUERY_SEARCH_RADIUS,
                                     QUERY_CANDIDATE_LIMIT);
  const std::size_t initial_candidate_count = candidate_elements.size();
  candidate_elements.reserve(candidate_elements.size() + m_fallback_vertices_scratch.size());
  appendCandidateElementsFromVertices(candidate_elements, m_fallback_vertices_scratch);

  PointLocationResult walk_result =
    walkCandidateElements(query_pt, candidate_elements, initial_candidate_count);
  if(walk_result.status != PointLocationStatus::Failed)
  {
    return walk_result;
  }

  const IndexType nearby_elem = findContainingElementNearby(query_pt, m_fallback_vertices_scratch);
  if(nearby_elem != INVALID_INDEX)
  {
    return {nearby_elem, PointLocationStatus::Found};
  }

  return {};
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE typename Delaunay<DIM>::IndexType
Delaunay<DIM>::findContainingElementFromNeighbors(const PointType& query_pt,
                                                  const std::vector<IndexType>& seed_elements) const
{
  if(seed_elements.empty())
  {
    return INVALID_INDEX;
  }

  std::vector<IndexType> nearby_elements;
  nearby_elements.reserve(seed_elements.size() * (1 + WALK_NEIGHBORHOOD_LAYERS * VERT_PER_ELEMENT));

  auto appendUniqueElement = [&](IndexType element_idx, std::vector<IndexType>& frontier) {
    if(isSearchableElement(element_idx) &&
       std::find(nearby_elements.begin(), nearby_elements.end(), element_idx) == nearby_elements.end())
    {
      nearby_elements.push_back(element_idx);
      frontier.push_back(element_idx);
    }
  };

  std::vector<IndexType> frontier;
  frontier.reserve(seed_elements.size());
  for(const IndexType element_idx : seed_elements)
  {
    appendUniqueElement(element_idx, frontier);
  }

  for(int layer = 0; layer < WALK_NEIGHBORHOOD_LAYERS && !frontier.empty(); ++layer)
  {
    std::vector<IndexType> next_frontier;
    next_frontier.reserve(frontier.size() * VERT_PER_ELEMENT);

    for(const IndexType element_idx : frontier)
    {
      const auto neighbors = m_mesh.adjacentElements(element_idx);
      for(int i = 0; i < VERT_PER_ELEMENT; ++i)
      {
        appendUniqueElement(neighbors[ModularFaceIndex(i) + 1], next_frontier);
      }
    }

    frontier.swap(next_frontier);
  }

  for(const IndexType element_idx : nearby_elements)
  {
    const BaryCoordType bary_coord = getBaryCoords(element_idx, query_pt);
    const DataType min_bary = bary_coord[bary_coord.array().argMin()];
    if(min_bary >= -BARY_EPS)
    {
      return element_idx;
    }
  }

  return INVALID_INDEX;
}

template <int DIM>
inline typename Delaunay<DIM>::IndexType Delaunay<DIM>::findContainingElementLinear(
  const PointType& query_pt,
  bool warnOnInvalid) const
{
  IndexType best_element = INVALID_INDEX;
  DataType best_min_bary = -std::numeric_limits<DataType>::max();

  for(auto element_idx : m_mesh.elements().positions())
  {
    if(!isSearchableElement(element_idx))
    {
      continue;
    }

    const BaryCoordType bary_coord = getBaryCoords(element_idx, query_pt);
    const DataType min_bary = bary_coord[bary_coord.array().argMin()];
    if(min_bary > best_min_bary)
    {
      best_min_bary = min_bary;
      best_element = element_idx;
    }

    if(min_bary >= -BARY_EPS)
    {
      return element_idx;
    }
  }

  SLIC_WARNING_IF(warnOnInvalid,
                  fmt::format("Unable to locate containing element for point {} after exhaustive "
                              "neighbor search; returning closest candidate with min barycentric "
                              "coordinate {:.17g}",
                              query_pt,
                              best_min_bary));

  return best_element;
}

template <int DIM>
AXOM_QUEST_DELAUNAY_FORCE_INLINE typename Delaunay<DIM>::IndexType
Delaunay<DIM>::findContainingElementNearby(const PointType& query_pt,
                                           const std::vector<IndexType>& nearby_vertices) const
{
  std::vector<IndexType> nearby_elements;
  for(const IndexType vertex_idx : nearby_vertices)
  {
    if(!m_mesh.isValidVertex(vertex_idx))
    {
      continue;
    }

    const auto star = m_mesh.vertexStar(vertex_idx);
    for(const IndexType elem : star)
    {
      if(isSearchableElement(elem))
      {
        nearby_elements.push_back(elem);
      }
    }
  }

  std::sort(nearby_elements.begin(), nearby_elements.end());
  nearby_elements.erase(std::unique(nearby_elements.begin(), nearby_elements.end()),
                        nearby_elements.end());

  for(const IndexType element_idx : nearby_elements)
  {
    const BaryCoordType bary_coord = getBaryCoords(element_idx, query_pt);
    const DataType min_bary = bary_coord[bary_coord.array().argMin()];
    if(min_bary >= -BARY_EPS)
    {
      return element_idx;
    }
  }

  return INVALID_INDEX;
}

template <int DIM>
inline typename Delaunay<DIM>::IndexType Delaunay<DIM>::findContainingElement(const PointType& query_pt,
                                                                              bool warnOnInvalid) const
{
  if(m_mesh.isEmpty())
  {
    SLIC_ERROR_IF(warnOnInvalid,
                  "Attempting to insert point into empty Delaunay triangulation."
                  "Delaunay::initializeBoundary() needs to be called first");
    return INVALID_INDEX;
  }
  if(!m_bounding_box.contains(query_pt))
  {
    SLIC_WARNING_IF(warnOnInvalid, "Attempting to locate element at location outside valid domain");
    return INVALID_INDEX;
  }

  m_candidate_elements_scratch.clear();
  getInitialCandidateElements(query_pt, m_candidate_elements_scratch);
  m_walked_elements_scratch.clear();
  PointLocationResult walk_result =
    walkCandidateElements(query_pt, m_candidate_elements_scratch, 0, &m_walked_elements_scratch);

  if(walk_result.status == PointLocationStatus::Found)
  {
    return walk_result.element_idx;
  }

  if(walk_result.status == PointLocationStatus::Outside)
  {
    return INVALID_INDEX;
  }

  if(!m_candidate_elements_scratch.empty())
  {
    const PointLocationResult fallback_result =
      findContainingElementWithQueryFallbacks(query_pt,
                                              m_candidate_elements_scratch,
                                              m_walked_elements_scratch);
    if(fallback_result.status == PointLocationStatus::Found)
    {
      return fallback_result.element_idx;
    }
  }

  if(m_collect_location_stats)
  {
    ++m_num_linear_fallbacks;
  }
  return findContainingElementLinear(query_pt, warnOnInvalid);
}

}  // namespace quest
}  // namespace axom

#endif
