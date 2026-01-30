// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_QUEST_FAST_WINDING_NUMBER_HPP
#define AXOM_QUEST_FAST_WINDING_NUMBER_HPP
#include "axom/primal.hpp"

#include <type_traits>

namespace axom
{
namespace quest
{
namespace detail
{
constexpr int get_num_moment_entries(int ndim, int ord)
{
  int n_entries = 0;
  for(int i = 0; i < ord + 1; ++i)
  {
    int power = 1;
    for(int j = 0; j < i + 1; ++j)
    {
      power *= ndim;
    }
    n_entries += power;
  }
  return n_entries;
}
}  // namespace detail

template <typename T, typename NDIMS, typename ORD>
class GWNMomentData
{
  static_assert((NDIMS == 2 || NDIMS == 3), "Must be defined in 2 or 3 dimensions");
  static_assert(0 <= ORD && ORD <= 2, "Only supported for orders 0, 1, or 2");

  static constexpr int NumberOfEntries = get_num_moment_entries(NDIMS, ORD);

  /// Addition overload to find the sum of two sets of raw moments
  friend GWNMomentData operator+(GWNMomentData& b1, GWNMomentData& b2)
  {
    GWNMomentData<T, NDIMS, ORD> b_out;

    b_out.a = b1.a + b2.a;
    b_out.ax = b1.ax + b2.ax;
    b_out.ay = b1.ay + b2.ay;
    b_out.az = b1.az + b2.az;

    for(int i = 0; i < NumberOfEntries; ++i) b_out.rm[i] = b1.rm[i] + b2.rm[i];

    b_out.compute_coefficients();

    return b_out;
  }

public:
  GWNMomentData() = default;

  SurfaceMomentData(const axom::primal::Triangle<T, 3>& a_tri)
  {
    static_assert(NDIMS == 3, "GWN Moments for triangles are defined only for 3D");

    // Track the centroid across the tree, and return the rest of the data
    auto centroid = a_tri.centroid();
    a = a_tri.area();
    ax = a * centroid[0];
    ay = a * centroid[1];
    az = a * centroid[2];

    auto normal = 0.5 * a_tri.normal();
    rm[0] = normal[0];
    rm[1] = normal[1];
    rm[2] = normal[2];

    if constexpr(ORD >= 1)
    {
      int m = 3;
      for(int i = 0; i < 9; ++i, ++m)
      {
        // In tensor product notation, equal to
        // centroid \otimes normal
        rm[m] = normal[i % 3] * centroid[i / 3];
      }
    }

    if constexpr(ORD >= 2)
    {
      constexpr auto twlv = 0.083333333333333333333333333333;
      const auto ab = axom::primal::Vector<T, 3> {a_tri[0].array() + a_tri[1].array()};
      const auto bc = axom::primal::Vector<T, 3> {a_tri[1].array() + a_tri[2].array()};
      const auto ac = axom::primal::Vector<T, 3> {a_tri[0].array() + a_tri[2].array()};

      for(int i = 0; i < 27; ++i, ++m)
      {
        // 1/12 * (ab \otimes ab + bc \otimes bc + ac \otimes ac) \otimes normal
        rm[m] = twlv * normal[i % 3] *
          (ab[i / 9] * ab[(i / 3) % 3] + bc[i / 9] * bc[(i / 3) % 3] + ac[i / 9] * ac[(i / 3) % 3]);
      }
    }

    compute_coefficients();
  }

  double winding_number(axom::primal::Point<T, NDIMS> query) const
  {
    if(axom::utilities::isNearlyEqual(std::abs(a), 0.0)) return 0.0;

    T terms[3] = {0.0, 0.0, 0.0};
    axom::primal::Vector<T, NDIMS> pq = ap / a - query;
    const double norm = pq.norm();

    if constexpr(NDIMS == 2)
    {
    }

    if constexpr(NDIMS == 3)
    {
      const axom::primal::Vector<T, 3> F1 {ec[0], ec[1], ec[2]};
      axom::primal::Vector<T, 3> G1 {pq[0], pq[1], pq[2]};

      // zeroth order
      const auto norm_to_3 = norm * norm * norm;
      terms[0] = F1.dot(G1) / norm_to_3;

      if constexpr(ORD >= 1)
      {
        const axom::primal::Vector<T, 9>
          F2 {ec[3], ec[4], ec[5], ec[6], ec[7], ec[8], ec[9], ec[10], ec[11]};
        // clang-format off
				axom::primal::Vector<T, 9> G2_1{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
				axom::primal::Vector<T, 9> G2_2{ pq[0] * pq[0], pq[1] * pq[0], pq[2] * pq[0],
																				 pq[0] * pq[1], pq[1] * pq[1], pq[2] * pq[1],
																				 pq[0] * pq[2], pq[1] * pq[2], pq[2] * pq[2] };
        // clang-format on

        const auto norm_to_5 = norm_to_3 * norm * norm;
        terms[1] = F2.dot(G2_1 * (1. / norm_to_3) - G2_2 * (3. / norm_to_5));

        if constexpr(ORD >= 2)
        {
          const axom::primal::Vector<T, 27> F3 {
            ec[12], ec[13], ec[14], ec[15], ec[16], ec[17], ec[18], ec[19], ec[20],
            ec[21], ec[22], ec[23], ec[24], ec[25], ec[26], ec[27], ec[28], ec[29],
            ec[30], ec[31], ec[32], ec[33], ec[34], ec[35], ec[36], ec[37], ec[38]};

          // clang-format off
				  axom::primal::Vector<T, 27> G3_1{ 3 * pq[0], pq[1], pq[2], pq[1], pq[0], 0, pq[2], 0, pq[0], pq[1], pq[0], 0, pq[0],3 * pq[1], pq[2], 0, pq[2], pq[1], pq[2], 0, pq[0], 0, pq[2], pq[1], pq[0], pq[1], 3 * pq[2] };
				  axom::primal::Vector<T, 27> G3_2{ pq[0] * pq[0] * pq[0], pq[0] * pq[0] * pq[1], pq[0] * pq[0] * pq[2], pq[0] * pq[1] * pq[0], pq[0] * pq[1] * pq[1], pq[0] * pq[1] * pq[2], pq[0] * pq[2] * pq[0], pq[0] * pq[2] * pq[1], pq[0] * pq[2] * pq[2],
																					  pq[1] * pq[0] * pq[0], pq[1] * pq[0] * pq[1], pq[1] * pq[0] * pq[2], pq[1] * pq[1] * pq[0], pq[1] * pq[1] * pq[1], pq[1] * pq[1] * pq[2], pq[1] * pq[2] * pq[0], pq[1] * pq[2] * pq[1], pq[1] * pq[2] * pq[2],
																					  pq[2] * pq[0] * pq[0], pq[2] * pq[0] * pq[1], pq[2] * pq[0] * pq[2], pq[2] * pq[1] * pq[0], pq[2] * pq[1] * pq[1], pq[2] * pq[1] * pq[2], pq[2] * pq[2] * pq[0], pq[2] * pq[2] * pq[1], pq[2] * pq[2] * pq[2] };
          // clang-format on

          const auto norm_to_7 = norm_to_5 * norm * norm;
          terms[2] = F3.dot(G3_1 * (-3. / norm_to_5) + G3_2 * (15. / norm_to_7));
        }
      }

      return 0.25 * M_1_PI * (terms[0] + terms[1] + terms[2]);
    }
  }

  axom::primal::Point<T, NDIMS> getSource() const
  {
    return axom::primal::Point<T, NDIMS>((ap / a).array());
  }

private:
  void compute_coefficients()
  {
    auto px = getSource();
    
    for (int i = 0; i < NDIMS; ++i)
    {
      ec[i] = rm[i];
    }

    if constexpr(ORD >= 1)
    {
      int m = NDIMS;
      for (int i = 0; i < NDIMS * NDIMS; ++i, ++m)
      {
        ec[m] = rm[m] - px[i / NDIMS] * rm[i % NDIMS];
      }

      if constexpr(ORD == 2)
      {
        //int m = NDIMS * NDIMS;
        //for(int i = 0; i < NDIMS * NDIMS * NDIMS; ++i, ++m)
        //{
        //  ec[m] = 0.5 * (rm[m] - );
        //}

        // Tensor notation forthcoming...
        ec[12] = 0.5 * (rm[12] - px[0] * rm[3] - px[0] * rm[3] + px[0] * px[0] * rm[0]);
        ec[13] = 0.5 * (rm[13] - px[0] * rm[4] - px[0] * rm[4] + px[0] * px[0] * rm[1]);
        ec[14] = 0.5 * (rm[14] - px[0] * rm[5] - px[0] * rm[5] + px[0] * px[0] * rm[2]);
        ec[15] = 0.5 * (rm[15] - px[0] * rm[6] - px[1] * rm[3] + px[0] * px[1] * rm[0]);
        ec[16] = 0.5 * (rm[16] - px[0] * rm[7] - px[1] * rm[4] + px[0] * px[1] * rm[1]);
        ec[17] = 0.5 * (rm[17] - px[0] * rm[8] - px[1] * rm[5] + px[0] * px[1] * rm[2]);
        ec[18] = 0.5 * (rm[18] - px[0] * rm[9] - px[2] * rm[3] + px[0] * px[2] * rm[0]);
        ec[19] = 0.5 * (rm[19] - px[0] * rm[10] - px[2] * rm[4] + px[0] * px[2] * rm[1]);
        ec[20] = 0.5 * (rm[20] - px[0] * rm[11] - px[2] * rm[5] + px[0] * px[2] * rm[2]);
        ec[21] = 0.5 * (rm[21] - px[1] * rm[3] - px[0] * rm[6] + px[1] * px[0] * rm[0]);
        ec[22] = 0.5 * (rm[22] - px[1] * rm[4] - px[0] * rm[7] + px[1] * px[0] * rm[1]);
        ec[23] = 0.5 * (rm[23] - px[1] * rm[5] - px[0] * rm[8] + px[1] * px[0] * rm[2]);
        ec[24] = 0.5 * (rm[24] - px[1] * rm[6] - px[1] * rm[6] + px[1] * px[1] * rm[0]);
        ec[25] = 0.5 * (rm[25] - px[1] * rm[7] - px[1] * rm[7] + px[1] * px[1] * rm[1]);
        ec[26] = 0.5 * (rm[26] - px[1] * rm[8] - px[1] * rm[8] + px[1] * px[1] * rm[2]);
        ec[27] = 0.5 * (rm[27] - px[1] * rm[9] - px[2] * rm[6] + px[1] * px[2] * rm[0]);
        ec[28] = 0.5 * (rm[28] - px[1] * rm[10] - px[2] * rm[7] + px[1] * px[2] * rm[1]);
        ec[29] = 0.5 * (rm[29] - px[1] * rm[11] - px[2] * rm[8] + px[1] * px[2] * rm[2]);
        ec[30] = 0.5 * (rm[30] - px[2] * rm[3] - px[0] * rm[9] + px[2] * px[0] * rm[0]);
        ec[31] = 0.5 * (rm[31] - px[2] * rm[4] - px[0] * rm[10] + px[2] * px[0] * rm[1]);
        ec[32] = 0.5 * (rm[32] - px[2] * rm[5] - px[0] * rm[11] + px[2] * px[0] * rm[2]);
        ec[33] = 0.5 * (rm[33] - px[2] * rm[6] - px[1] * rm[9] + px[2] * px[1] * rm[0]);
        ec[34] = 0.5 * (rm[34] - px[2] * rm[7] - px[1] * rm[10] + px[2] * px[1] * rm[1]);
        ec[35] = 0.5 * (rm[35] - px[2] * rm[8] - px[1] * rm[11] + px[2] * px[1] * rm[2]);
        ec[36] = 0.5 * (rm[36] - px[2] * rm[9] - px[2] * rm[9] + px[2] * px[2] * rm[0]);
        ec[37] = 0.5 * (rm[37] - px[2] * rm[10] - px[2] * rm[10] + px[2] * px[2] * rm[1]);
        ec[38] = 0.5 * (rm[38] - px[2] * rm[11] - px[2] * rm[11] + px[2] * px[2] * rm[2]);
      }
    }
  }

public:
  // Store accumulated values across the node's children
  axom::primal::Vector<T, NDIMS> ap;  // Scaled centroid
  double a {};

  // Raw moments {
  axom::StackArray<T, NumberOfEntries> rm {};

  // Expansion coefficients
  axom::StackArray<T, NumberOfEntries> ec {};

  // Store the order of the approximation for more efficient processing later
  int m_order;

#ifdef JACOBS_DEBUG
  axom::Array<int> components;
#endif
};

template <typename T>
axom::Array<primal::NURBSPatch<T, 3>> process_surfaces(
  const axom::Array<primal::NURBSPatch<T, 3>>& input_surfaces,
  double threshold,
  int npasses = 10)
{
  using BoxType = primal::BoundingBox<T, 3>;
  using NURBSType = primal::NURBSPatch<T, 3>;

  // Create initial array of processed patches,
  //  beginning by clipping each patch parameter space
  //  to a bounding box of its trimming curves
  // Compute a bounding box of all the surfaces
  axom::Array<NURBSType> candidates;
  candidates.reserve(input_surfaces.size());
  BoxType total_bbox;
  for(auto& surf : input_surfaces)
  {
    // This is where we would do Bezier extraction, if it worked better :(
    //auto beziers = surf.extractTrimmedBezier();
    //for(auto& bez : beziers)
    {
      auto the_patch = surf;

      if(the_patch.getNumTrimmingCurves() == 0) continue;

      the_patch.normalize();
      the_patch.clipToCurves();

      // Re-check if the patch is empty after clipping to curve
      if(the_patch.getNumTrimmingCurves() == 0) continue;

      candidates.push_back(the_patch);
      total_bbox.addBox(the_patch.boundingBox());
    }
  }

  // Iterate over all the surfaces until no patch has a bounding box
  //  bigger than threshold * (total_bbox's size)
  for(int i = 0; i < npasses; ++i)
  {
    axom::Array<NURBSType> subdivisions;
    subdivisions.reserve(candidates.size() * 3 / 2);

    BoxType new_bbox;

    // If any patch is bigger than the threshold, subdivide it, clip it,
    //  and add it to the next level. Repeat as needed.
    const double max_range_norm = threshold * total_bbox.range().norm();
    for(const auto& candidate : candidates)
    {
      if(candidate.boundingBox().range().norm() < max_range_norm)
      {
        new_bbox.addBox(candidate.boundingBox());
        subdivisions.push_back(candidate);
        continue;
      }

      NURBSType subpatches[4];
      candidate.nearBisect(subpatches[0], subpatches[1], subpatches[2], subpatches[3]);
      for(int si = 0; si < 4; si++)
      {
        if(subpatches[si].getNumTrimmingCurves() == 0) continue;

        subdivisions.emplace_back(std::move(subpatches[si]));
        subdivisions.back().clipToCurves();
        new_bbox.addBox(subdivisions.back().boundingBox());
      }
    }

    // Break if no additional subdivisions are made
    if(candidates.size() == subdivisions.size()) break;

    candidates.swap(subdivisions);
    total_bbox = new_bbox;
  }

  return candidates;
}

// I know this method is inefficient, but it's only for figure visualization
template <typename T>
axom::Array<primal::Triangle<T, 3>> process_triangles(
  const axom::Array<primal::Triangle<T, 3>>& input_triangles,
  int npasses = 3)
{
  axom::Array<primal::Triangle<T, 3>> candidates = input_triangles;
  for(int i = 0; i < npasses; ++i)
  {
    axom::Array<primal::Triangle<T, 3>> subdivisions;
    subdivisions.reserve(candidates.size() * 3 / 2);

    double min_area = candidates[0].area();
    double max_area = candidates[0].area();

    for(auto& tri : candidates)
    {
      min_area = axom::utilities::min(min_area, tri.area());
      max_area = axom::utilities::max(max_area, tri.area());
    }

    // Any triangle with area bigger than 0.8 * (max_area - min_area), subdivide into 4
    double area_threshold = 0.8 * (max_area - min_area);
    for(auto& tri : candidates)
    {
      if(tri.area() > area_threshold)
      {
        const auto& A = tri[0];
        const auto& B = tri[1];
        const auto& C = tri[2];
        const primal::Point<T, 3> P1(0.5 * (A.array() + B.array()));
        const primal::Point<T, 3> P2(0.5 * (B.array() + C.array()));
        const primal::Point<T, 3> P3(0.5 * (C.array() + A.array()));

        subdivisions.push_back(primal::Triangle<T, 3>(A, P1, P3));
        subdivisions.push_back(primal::Triangle<T, 3>(P1, B, P2));
        subdivisions.push_back(primal::Triangle<T, 3>(P1, P2, P3));
        subdivisions.push_back(primal::Triangle<T, 3>(P2, C, P3));
      }
      else
        subdivisions.push_back(tri);
    }

    candidates = subdivisions;
  }

  return candidates;
}

template <typename T, typename LeafGeometry, typename TraverserType>
double winding_number(const primal::Point<T, 3>& query,
                      const TraverserType& traverser,
                      const ArrayView<LeafGeometry>& leaf_caches,
                      const ArrayView<SurfaceMomentData<T>>& interior_moments,
                      bool useDirect = false,
                      const double edge_tol = 1e-8,
                      const double ls_tol = 1e-8,
                      const double quad_tol = 1e-8,
                      const double disk_size = 0.01,
                      const double EPS = 1e-8)
{
  double gwn = 0.0;

  const auto leaf_caches_view = leaf_caches;
  const auto interior_moments_view = interior_moments;
  using LeafGeom = std::decay_t<LeafGeometry>;

  // NOTE: "query" capture was removed here because the parameter list contains a "query"
  //       variable. They seem like they would be the same thing but not 100% sure.
  auto bbContain = [&gwn, &interior_moments_view, &useDirect](const primal::Point<T, 3>& query,
                                                              const primal::BoundingBox<T, 3>& bvhBbox,
                                                              std::int32_t node_index) -> bool {
    constexpr double beta = 2.0;
    const bool near_tree =
      axom::primal::squared_distance(query, interior_moments_view[node_index].getSource()) <
      beta * beta * bvhBbox.range().squared_norm() / 4;

    // If we're inside an internal node bbox, need to keep recurring
    if(useDirect || near_tree)
    {
      return true;
    }
    // If we're outside an internal node bbox, we can add its contribution to the GWN
    else
    {
      gwn += interior_moments_view[node_index].winding_number(query);
      return false;
    }
  };

  if constexpr(std::is_same_v<LeafGeom, axom::primal::Triangle<T, 3>>)
  {
    auto leaf_gwn = [&query, &gwn, leaf_caches_view, edge_tol, EPS](
                      std::int32_t currentNode,
                      const std::int32_t* leafNodes) -> void {
      const auto idx = leafNodes[currentNode];
      gwn += axom::primal::winding_number(query, leaf_caches_view[idx], edge_tol, EPS);
    };

    traverser.traverse_tree(query, leaf_gwn, bbContain);
  }
  else if constexpr(std::is_same_v<LeafGeom, axom::primal::Dipole<T>>)
  {
    auto leaf_gwn = [&query, &gwn, leaf_caches_view, edge_tol, EPS](
                      std::int32_t currentNode,
                      const std::int32_t* leafNodes) -> void {
      const auto idx = leafNodes[currentNode];
      gwn += axom::primal::winding_number(query, leaf_caches_view[idx], edge_tol, EPS);
    };

    traverser.traverse_tree(query, leaf_gwn, bbContain);
  }
  else
  {
    auto leaf_gwn = [&query, &gwn, leaf_caches_view, edge_tol, ls_tol, quad_tol, disk_size, EPS](
                      std::int32_t currentNode,
                      const std::int32_t* leafNodes) -> void {
      const auto idx = leafNodes[currentNode];
      gwn += axom::primal::winding_number(query,
                                          leaf_caches_view[idx],
                                          edge_tol,
                                          ls_tol,
                                          quad_tol,
                                          disk_size,
                                          EPS);
    };

    traverser.traverse_tree(query, leaf_gwn, bbContain);
  }

  return gwn;
}

template <typename T>
void normalize_volume_by_bbox(ArrayView<axom::primal::NURBSPatch<T, 3>> surfaces,
                              axom::primal::BoundingBox<double, 3> bbox)
{
  auto centroid = bbox.getCentroid();
  auto longest_dim = bbox.getLongestDimension();
  auto scale = bbox.getMax()[longest_dim] - bbox.getMin()[longest_dim];

  for(auto& surf : surfaces)
  {
    for(int i = 0; i < surf.getNumControlPoints_u(); ++i)
      for(int j = 0; j < surf.getNumControlPoints_v(); ++j)
        surf(i, j).array() = (surf(i, j).array() - centroid.array()) / scale;
  }
}

template <typename T>
void normalize_volume_by_bbox(ArrayView<axom::primal::Triangle<T, 3>> triangles,
                              axom::primal::BoundingBox<double, 3> bbox)
{
  auto centroid = bbox.getCentroid();
  auto longest_dim = bbox.getLongestDimension();
  auto scale = bbox.getMax()[longest_dim] - bbox.getMin()[longest_dim];

  for(auto& tri : triangles)
  {
    tri[0].array() = (tri[0].array() - centroid.array()) / scale;
    tri[1].array() = (tri[1].array() - centroid.array()) / scale;
    tri[2].array() = (tri[2].array() - centroid.array()) / scale;
  }
}

}  // end namespace quest
}  // end namespace axom

#endif
