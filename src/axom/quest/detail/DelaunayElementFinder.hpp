// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_DETAIL_DELAUNAY_ELEMENT_FINDER_HPP_
#define AXOM_QUEST_DETAIL_DELAUNAY_ELEMENT_FINDER_HPP_

#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include <utility>
#include <vector>

namespace axom
{
namespace quest
{
namespace detail
{

template <int DIM, typename PointType, typename IAMeshType, typename BoundingBox, typename IndexType>
class DelaunayElementFinder
{
public:
  using NumericArrayType = NumericArray<IndexType, DIM>;
  using LatticeType = spin::RectangularLattice<DIM, double, IndexType>;

  explicit DelaunayElementFinder() = default;

  void recomputeGrid(const IAMeshType& mesh, const BoundingBox& bb)
  {
    const auto& verts = mesh.vertices();

    // Choose the grid resolution so that each bin contains ~O(1) points per
    // dimension on average. This keeps the "nearby bin" seed used by point
    // insertion and query walks close (in terms of simplex adjacency hops)
    // even for very large point sets.
    //
    // Target occupancy is ~ 4^DIM points per bin (16 in 2D, 64 in 3D).
    constexpr double BIN_SIDE_SPACING = 4.0;
    const double res_root = std::pow(static_cast<double>(verts.size()), 1.0 / DIM);
    const IndexType res =
      axom::utilities::max(IndexType {2},
                           static_cast<IndexType>(std::ceil(res_root / BIN_SIDE_SPACING)));

    auto expandedBB = BoundingBox(bb).scale(1.05);

    m_lattice = spin::rectangular_lattice_from_bounding_box(expandedBB, NumericArrayType(res));

    resizeArray<DIM>(res);
    m_bins.fill(INVALID_INDEX);

    for(auto idx : verts.positions())
    {
      if(!mesh.isValidVertex(idx))
      {
        continue;
      }

      const IndexType coboundary = mesh.coboundaryElement(idx);
      if(!mesh.isValidElement(coboundary))
      {
        continue;
      }

      const auto& pos = mesh.getVertexPosition(idx);
      const auto cell = m_lattice.gridCell(pos);
      IndexType& slot = flatIndex(cell);
      if(!mesh.isValidVertex(slot) || !mesh.isValidElement(mesh.coboundaryElement(slot)))
      {
        slot = idx;
      }
    }
  }

  inline void getNearbyVertices(const IAMeshType& mesh,
                                const PointType& pt,
                                std::vector<IndexType>& nearby_vertices,
                                int search_radius = 1,
                                int max_candidates = 1) const
  {
    const auto cell = m_lattice.gridCell(pt);
    m_candidate_scratch.clear();
    const int span = 2 * search_radius + 1;
    const int max_bins = (DIM == 2) ? (span * span) : (span * span * span);
    m_candidate_scratch.reserve(static_cast<std::size_t>(max_bins));

    auto tryCandidate = [&](const typename LatticeType::GridCell& candidate_cell) {
      const IndexType vertex_idx = flatIndex(candidate_cell);
      if(mesh.isValidVertex(vertex_idx) && mesh.isValidElement(mesh.coboundaryElement(vertex_idx)))
      {
        const double sq_dist = primal::squared_distance(mesh.getVertexPosition(vertex_idx), pt);
        m_candidate_scratch.emplace_back(sq_dist, vertex_idx);
      }
    };

    if constexpr(DIM == 2)
    {
      for(int dj = -search_radius; dj <= search_radius; ++dj)
      {
        const IndexType j = cell[1] + dj;
        if(j < 0 || j >= m_bins.shape()[1])
        {
          continue;
        }

        for(int di = -search_radius; di <= search_radius; ++di)
        {
          const IndexType i = cell[0] + di;
          if(i < 0 || i >= m_bins.shape()[0])
          {
            continue;
          }

          tryCandidate(typename LatticeType::GridCell {{i, j}});
        }
      }
    }
    else
    {
      for(int dk = -search_radius; dk <= search_radius; ++dk)
      {
        const IndexType k = cell[2] + dk;
        if(k < 0 || k >= m_bins.shape()[2])
        {
          continue;
        }

        for(int dj = -search_radius; dj <= search_radius; ++dj)
        {
          const IndexType j = cell[1] + dj;
          if(j < 0 || j >= m_bins.shape()[1])
          {
            continue;
          }

          for(int di = -search_radius; di <= search_radius; ++di)
          {
            const IndexType i = cell[0] + di;
            if(i < 0 || i >= m_bins.shape()[0])
            {
              continue;
            }

            tryCandidate(typename LatticeType::GridCell {{i, j, k}});
          }
        }
      }
    }

    if(static_cast<int>(m_candidate_scratch.size()) > max_candidates)
    {
      auto kth = m_candidate_scratch.begin() + max_candidates;
      std::nth_element(m_candidate_scratch.begin(),
                       kth,
                       m_candidate_scratch.end(),
                       [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
      m_candidate_scratch.resize(static_cast<std::size_t>(max_candidates));
    }

    std::sort(m_candidate_scratch.begin(),
              m_candidate_scratch.end(),
              [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    nearby_vertices.clear();
    nearby_vertices.reserve(
      axom::utilities::min(max_candidates, static_cast<int>(m_candidate_scratch.size())));
    for(const auto& candidate : m_candidate_scratch)
    {
      nearby_vertices.push_back(candidate.second);
    }
  }

  inline IndexType getNearbyVertex(const PointType& pt) const
  {
    const auto cell = m_lattice.gridCell(pt);
    return flatIndex(cell);
  }

  inline void updateBin(const PointType& pt, IndexType vertex_id)
  {
    const auto cell = m_lattice.gridCell(pt);
    flatIndex(cell) = vertex_id;
  }

private:
  static constexpr IndexType INVALID_INDEX = IndexType {-1};

  inline IndexType& flatIndex(const typename LatticeType::GridCell& cell)
  {
    const IndexType idx = numerics::dot_product(cell.data(), m_bins.strides().begin(), DIM);
    return m_bins.flatIndex(idx);
  }

  inline const IndexType& flatIndex(const typename LatticeType::GridCell& cell) const
  {
    const IndexType idx = numerics::dot_product(cell.data(), m_bins.strides().begin(), DIM);
    return m_bins.flatIndex(idx);
  }

  template <int TDIM>
  typename std::enable_if<TDIM == 2, void>::type resizeArray(IndexType res)
  {
    m_bins.resize(res, res);
  }

  template <int TDIM>
  typename std::enable_if<TDIM == 3, void>::type resizeArray(IndexType res)
  {
    m_bins.resize(res, res, res);
  }

private:
  axom::Array<IndexType, DIM> m_bins;
  LatticeType m_lattice;
  mutable std::vector<std::pair<double, IndexType>> m_candidate_scratch;
};

}  // namespace detail
}  // namespace quest
}  // namespace axom

#endif
