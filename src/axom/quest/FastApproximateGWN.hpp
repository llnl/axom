// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_QUEST_FAST_APPROXIMATE_GWN_HPP
#define AXOM_QUEST_FAST_APPROXIMATE_GWN_HPP
#include "axom/primal.hpp"

#include <type_traits>

namespace axom
{
namespace quest
{
namespace detail
{

/// Computes the total number of moments for a given degree and order, i.e.
/// ndim + ndim^2 + ... + ndim^(ord + 1)
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

/*
 * \class GWN Moment Data
 *
 * \brief A class to compute and store the geometric moment data which parameterizes the 
 *  Taylor expansion of a GWN approximation for a cluster of geometric primitives. 
 * \tparam T The numeric type
 * \tparam NDIMS The number of spatial dimensions of the geometric primitives (2 or 3)
 * \tparam ORD The order of the Taylor expansion (0, 1, or 2)
 * 
 * Stores arrays `rm` of raw moments and `ec` of Taylor expansion coefficients, as only the 
 * latter is used to compute the approximated GWN. Transforming raw moments to expansion coefficients
 * requires knowing the center of the expansion, which we take as the centroid of the collection.
 * 
 * Stores double `a` and vector `ap` which are the total unsigned area and weighted centroid of
 * the collection of geometric objects. 
 */
template <typename T, int NDIMS, int ORD>
class GWNMomentData
{
  static_assert((NDIMS == 2 || NDIMS == 3), "Must be defined in 2 or 3 dimensions");
  static_assert(0 <= ORD && ORD <= 2, "Only supported for orders 0, 1, or 2");

  static constexpr int NumberOfEntries = detail::get_num_moment_entries(NDIMS, ORD);

  /// Addition overload to find the sum of two sets of raw moments.
  /// TODO: Technically, the raw moments for b1 and b2 can be deallocated after this
  ///  function is called, which would decrease the memory footprint
  friend GWNMomentData operator+(GWNMomentData& b1, GWNMomentData& b2)
  {
    GWNMomentData<T, NDIMS, ORD> b_out;

    b_out.a = b1.a + b2.a;
    b_out.ap = b1.ap + b2.ap;

    for(int i = 0; i < NumberOfEntries; ++i)
    {
      b_out.rm[i] = b1.rm[i] + b2.rm[i];
    }

    b_out.compute_coefficients();

    return b_out;
  }

public:
  GWNMomentData() = default;

  /// Construct moments from a 3D triangle
  explicit GWNMomentData(const axom::primal::Triangle<T, 3>& a_tri)
  {
    static_assert(NDIMS == 3, "GWN Moments for triangles are defined only for 3D");

    // Track the centroid across the tree, and return the rest of the data
    auto centroid = a_tri.centroid();
    a = a_tri.area();
    ap[0] = a * centroid[0];
    ap[1] = a * centroid[1];
    ap[2] = a * centroid[2];

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

      if constexpr(ORD >= 2)
      {
        constexpr auto twlv = 1.0 / 12.0;
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
    }

    compute_coefficients();
  }

  /// Construct moments from a trimmed NURBS surface
  explicit GWNMomentData(const axom::primal::NURBSPatch<T, 3>& a_patch)
  {
    const auto patch_data = a_patch.calculateSurfaceMoments<ORD>();

    a = patch_data[0];
    ap[0] = patch_data[1];
    ap[1] = patch_data[2];
    ap[2] = patch_data[3];

    for(int i = 0; i < NumberOfEntries; ++i) rm[i] = patch_data[i + 4];

    compute_coefficients();
  }

  /// Construct moments from the endpoints of a 2D segment
  explicit GWNMomentData(const axom::primal::NURBSCurve<T, 2>& c)
    : GWNMomentData(c.getInitPoint(), c.getEndPoint())
  { }

  /// Construct moments from a 2D Segment
  explicit GWNMomentData(const axom::primal::Segment<T, 2>& s)
    : GWNMomentData(s.source(), s.target())
  { }

  /// Construct moments from the endpoints of a 2D segment
  explicit GWNMomentData(const axom::primal::Point<T, 2>& p0, const axom::primal::Point<T, 2>& p1)
  {
    static_assert(NDIMS == 2, "GWN Moments for segments are defined only for 2D");

    const auto x0 = p0[0];
    const auto y0 = p0[1];
    const auto x1 = p1[0];
    const auto y1 = p1[1];

    // Needed to track the centroid across the tree
    const auto dx = x1 - x0;
    const auto dy = y0 - y1;  // actually -dy since it was more useful.
    a = sqrt(dx * dx + dy * dy);
    ap[0] = a * 0.5 * (x0 + x1);
    ap[1] = a * 0.5 * (y0 + y1);

    const auto x0x0 = x0 * x0;
    const auto x0x1 = x0 * x1;
    const auto x1x1 = x1 * x1;

    const auto y0y0 = y0 * y0;
    const auto y0y1 = y0 * y1;
    const auto y1y1 = y1 * y1;

    const auto x0y0 = x0 * y0;
    const auto x0y1 = x0 * y1;
    const auto x1y0 = x1 * y0;
    const auto x1y1 = x1 * y1;

    rm[0] = dy;
    rm[1] = dx;

    if constexpr(ORD >= 1)
    {
      rm[2] = 0.5 * dy * (x0 + x1);
      rm[3] = 0.5 * (x1x1 - x0x0);
      rm[4] = 0.5 * (y0y0 - y1y1);
      rm[5] = 0.5 * (dx) * (y0 + y1);

      if constexpr(ORD == 2)
      {
        const auto A = (x0x0 + x0x1 + x1x1) / 3.0;
        rm[6] = dy * A;
        const auto B = (x0y0 + 0.5 * (x0y1 + x1y0) + x1y1) / 3.0;
        rm[7] = dy * B;
        rm[8] = rm[7];
        const auto C = (y0y0 + y0y1 + y1y1) / 3.0;
        rm[9] = dy * C;
        rm[10] = dx * A;
        rm[11] = dx * B;
        rm[12] = rm[11];
        rm[13] = dx * C;
      }
    }

    compute_coefficients();
  }

  /// Computes the approximated GWN field at the given query.
  ///  Formulae are taken from "Fast Winding Numbers for Soups and Clouds" by
  ///  Barill et al. (2018)
  double approx_winding_number(axom::primal::Point<T, NDIMS> query) const
  {
    if(axom::utilities::isNearlyEqual(std::abs(a), 0.0)) return 0.0;

    T terms[3] = {0.0, 0.0, 0.0};
    axom::primal::Vector<T, NDIMS> pq(query, getCenter());
    const double norm = pq.norm();

    if constexpr(NDIMS == 2)
    {
      const axom::primal::Vector<double, 2> F1 {ec[0], ec[1]};
      axom::primal::Vector<double, 2> G1 {pq[0], pq[1]};

      const auto norm_to_2 = norm * norm;
      terms[0] = F1.dot(G1) / norm_to_2;

      if constexpr(ORD >= 1)
      {
        const axom::primal::Vector<double, 4> F2 {ec[2], ec[3], ec[4], ec[5]};
        axom::primal::Vector<double, 4> G2 {pq[0] * pq[0] - pq[1] * pq[1], 2 * pq[0] * pq[1], 0.0, 0.0};
        G2[2] = G2[1];
        G2[3] = -G2[0];

        const auto norm_to_4 = norm_to_2 * norm * norm;
        terms[1] = -F2.dot(G2) / norm_to_4;

        if constexpr(ORD == 2)
        {
          const axom::primal::Vector<double, 8>
            F3 {ec[6], ec[7], ec[8], ec[9], ec[10], ec[11], ec[12], ec[13]};
          const double g0 = 2 * pq[0] * (pq[0] * pq[0] - 3 * pq[1] * pq[1]);
          const double g1 = 2 * pq[1] * (3 * pq[0] * pq[0] - pq[1] * pq[1]);
          axom::primal::Vector<double, 8> G3 {g0, g1, g1, -g0, g1, -g0, -g0, -g1};

          const auto norm_to_6 = norm_to_4 * norm * norm;
          terms[2] = F3.dot(G3) / norm_to_6;
        }
      }

      return -0.5 * M_1_PI * (terms[0] + terms[1] + terms[2]);
    }

    if constexpr(NDIMS == 3)
    {
      const axom::primal::Vector<T, 3> F1 {ec[0], ec[1], ec[2]};
      axom::primal::Vector<T, 3> G1 {pq[0], pq[1], pq[2]};

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

          // The formula in [Barill 2018] incorrectly lists (-1. / norm_to_5)
          const auto norm_to_7 = norm_to_5 * norm * norm;
          terms[2] = F3.dot(G3_1 * (-3. / norm_to_5) + G3_2 * (15. / norm_to_7));
        }
      }

      return 0.25 * M_1_PI * (terms[0] + terms[1] + terms[2]);
    }
  }

  /// Return the center of the Taylor expansion
  axom::primal::Point<T, NDIMS> getCenter() const
  {
    if(a == 0)
    {
      return axom::primal::Point<T, NDIMS> {};
    }

    return axom::primal::Point<T, NDIMS>((ap / a).array());
  }

  /// Return the normal if computed from 3D data
  axom::primal::Vector<T, 3> getNormal() const
  {
    static_assert(NDIMS == 3, "GWN Moments for triangles are defined only for 3D");
    return axom::primal::Vector<T, 3> {ec[0], ec[1], ec[2]};
  }

private:
  /// Transform raw moments into expansion coefficients
  void compute_coefficients()
  {
    auto p = getCenter();

    for(int i = 0; i < NDIMS; ++i)
    {
      ec[i] = rm[i];
    }

    if constexpr(ORD >= 1)
    {
      int m = NDIMS;
      for(int i = 0; i < NDIMS * NDIMS; ++i, ++m)
      {
        ec[m] = rm[m] - p[i / NDIMS] * rm[i % NDIMS];
      }

      if constexpr(ORD == 2)
      {
        for(int i = 0; i < NDIMS * NDIMS * NDIMS; ++i, ++m)
        {
          // Example values for NDIMS = 2
          const int A = (i / NDIMS / NDIMS) % NDIMS;  // 0, 0, 0, 0, 1, 1, 1, 1
          const int B = (i / NDIMS) % NDIMS;          // 0, 0, 1, 1, 0, 0, 1, 1
          const int C = (i / 1) % NDIMS;              // 0, 1, 0, 1, 0, 1, 0, 1
          const int D = i % (NDIMS * NDIMS);          // 0, 1, 2, 3, 0, 1, 2, 3

          ec[m] = 0.5 *
            (rm[m] - p[A] * rm[NDIMS + D] - p[B] * rm[(A + 1) * NDIMS + C] + p[A] * p[B] * rm[C]);
        }
      }
    }
  }

public:
  // Store accumulated values across the node's children
  axom::primal::Vector<T, NDIMS> ap;  // Weighted centroid
  double a {};

  // Raw moments
  axom::StackArray<T, NumberOfEntries> rm {};

  // Expansion coefficients
  axom::StackArray<T, NumberOfEntries> ec {};
};

/*!
 * \brief Evaluate a hierarchical approximation of the GWN for the shape in 
 *   `leaf_objects`, encompassed by `traverser`'s BVH
 *
 * Traverse the BVH at the input `query`. 
 *  For any node which is considered "far-away", approximate the GWN for that node
 *  For any node which is a leaf (i.e. close to the query), compute the GWN directly
 * 
 * \tparam T Numeric type of geometric primitives
 * \tparam NDIMS The spatial dimensions of the shape (2 or 3)
 * \tparam ORD The number of terms in the Taylor expansion
 * \tparam LeafGeometry The type of the objects which could require direct evaluation
 * \tparam TraverserType The derived type of the traverser for the BVH tree
 * \param [in] query The query point at which to evaluate the GWN
 * \param [in] traverser A traverser for the BVH tree 
 * \param [in] leaf_objects A view of all individual geometric objects in the shape,
 *                           indexed by position in the BVH tree
 * \param [in] internal_moments A view of GWNMomentData objects for each internal node 
 *                               of the BVH tree, each representing a cluster of leaf objects
 * \param [in] wt A structure of possible tolerances for GWN evaluation to permit
 *                 flexible evaluation for many different type of leaf objects
 * \param [in] beta An "accuracy parameter" which scales at what distance from an AABB a query
 *                   point is considered to be "far-away".
 * 
 * The default parameter beta = 2.0 is suggested by the work 
 * "Fast Winding Numbers for Soups and Clouds" by Barill et al. (2018)
 * 
 * \return The approximated GWN at the query point
 */
template <typename T, int NDIMS, int ORD, typename LeafGeometry, typename TraverserType>
double fast_approximate_winding_number(const primal::Point<T, NDIMS>& query,
                                       const TraverserType& traverser,
                                       const ArrayView<LeafGeometry>& leaf_objects,
                                       const ArrayView<GWNMomentData<T, NDIMS, ORD>>& internal_moments,
                                       const primal::WindingTolerances& wt,
                                       double beta = 2.0)
{
  double gwn = 0.0;

  const auto leaf_objects_view = leaf_objects;
  const auto internal_moments_view = internal_moments;
  using LeafGeom = std::decay_t<LeafGeometry>;

  auto bbContain = [&gwn, &internal_moments_view, &beta](const primal::Point<T, NDIMS>& query,
                                                         const primal::BoundingBox<T, NDIMS>& bvhBbox,
                                                         std::int32_t node_index) -> bool {
    const bool near_tree =
      axom::primal::squared_distance(query, internal_moments_view[node_index].getCenter()) <
      beta * beta * bvhBbox.range().squared_norm() / 4;

    // If we're inside an internal node bbox, need to keep recurring
    if(near_tree)
    {
      return true;
    }
    // If we're outside an internal node bbox, we can add its contribution to the GWN
    else
    {
      gwn += internal_moments_view[node_index].approx_winding_number(query);
      return false;
    }
  };

  if constexpr(std::is_same_v<LeafGeom, axom::primal::Triangle<T, 3>> ||
               std::is_same_v<LeafGeom, axom::primal::NURBSCurve<T, 2>> ||
               std::is_same_v<LeafGeom, axom::primal::detail::NURBSCurveGWNCache<T>>)
  {
    auto leaf_gwn = [&query, &gwn, leaf_objects_view, &wt](std::int32_t currentNode,
                                                           const std::int32_t* leafNodes) -> void {
      const auto idx = leafNodes[currentNode];
      gwn += axom::primal::winding_number(query, leaf_objects_view[idx], wt.edge_tol, wt.EPS);
    };

    traverser.traverse_tree(query, leaf_gwn, bbContain);
  }

  if constexpr(std::is_same_v<LeafGeom, axom::primal::Segment<T, 2>>)
  {
    auto leaf_gwn = [&query, &gwn, leaf_objects_view, &wt](std::int32_t currentNode,
                                                           const std::int32_t* leafNodes) -> void {
      const auto idx = leafNodes[currentNode];
      gwn += axom::primal::winding_number(query, leaf_objects_view[idx], wt.edge_tol);
    };

    traverser.traverse_tree(query, leaf_gwn, bbContain);
  }

  if constexpr(std::is_same_v<LeafGeom, axom::primal::NURBSPatch<T, 3>> ||
               std::is_same_v<LeafGeom, axom::primal::detail::NURBSPatchGWNCache<T>>)
  {
    auto leaf_gwn = [&query, &gwn, leaf_objects_view, &wt](std::int32_t currentNode,
                                                           const std::int32_t* leafNodes) -> void {
      const auto idx = leafNodes[currentNode];
      gwn += axom::primal::winding_number(query,
                                          leaf_objects_view[idx],
                                          wt.edge_tol,
                                          wt.ls_tol,
                                          wt.quad_tol,
                                          wt.disk_size,
                                          wt.EPS);
    };

    traverser.traverse_tree(query, leaf_gwn, bbContain);
  }
  // Support for other leaf types forthcoming...

  return gwn;
}

template <typename T>
axom::Array<primal::NURBSCurve<T, 2>> subdivide_curves(
  const axom::ArrayView<const primal::NURBSCurve<T, 2>>& input_curves_view,
  double bbox_threshold,
  int npasses = 10)
{
  using BoxType = primal::BoundingBox<T, 2>;
  using NURBSType = primal::NURBSCurve<T, 2>;
  using BezierType = primal::BezierCurve<T, 2>;

  // Compute a bounding box of all the curves
  axom::Array<BezierType> candidates;
  BoxType total_bbox;

  // For NURBSCurves, first do a pass of Bezier extraction
  for(auto& curv : input_curves_view)
  {
    for(auto& bez : curv.extractBezier())
    {
      candidates.push_back(bez);
      total_bbox.addBox(bez.boundingBox());
    }
  }

  // Iterate over all the curves until none have a bounding box
  //  bigger than threshold * (total_bbox's size)
  for(int i = 0; i < npasses; ++i)
  {
    axom::Array<BezierType> subdivisions;
    subdivisions.reserve(candidates.size() * 3 / 2);

    BoxType new_bbox;

    // If any patch is bigger than the threshold, subdivide it,
    //  and add it to the next level. Repeat as needed.
    const double max_range_norm = bbox_threshold * total_bbox.range().norm();
    for(const auto& candidate : candidates)
    {
      if(candidate.boundingBox().range().norm() < max_range_norm)
      {
        new_bbox.addBox(candidate.boundingBox());
        subdivisions.push_back(candidate);
        continue;
      }

      BezierType subcurves[2];
      candidate.split(0.5, subcurves[0], subcurves[1]);
      for(int si = 0; si < 2; si++)
      {
        subdivisions.emplace_back(std::move(subcurves[si]));
        new_bbox.addBox(subdivisions.back().boundingBox());
      }
    }

    // Break if no additional subdivisions are made
    if(candidates.size() == subdivisions.size()) break;

    candidates.swap(subdivisions);
    total_bbox = new_bbox;
  }

  // Do one final pass to turn the array of candidates into NURBS
  axom::Array<NURBSType> candidates_nurbs(0, candidates.size());
  for(auto& c : candidates)
  {
    candidates_nurbs.emplace_back(NURBSType(c));
  }

  return candidates_nurbs;
}

template <typename T>
axom::Array<primal::NURBSPatch<T, 3>> subdivide_patches(
  const axom::ArrayView<const primal::NURBSPatch<T, 3>>& input_patches_view,
  double bbox_threshold,
  int npasses = 10)
{
  using BoxType = primal::BoundingBox<T, 3>;
  using NURBSType = primal::NURBSPatch<T, 3>;

  axom::Array<NURBSType> candidates;
  candidates.reserve(input_patches_view.size() * 3 / 2);
  BoxType total_bbox;

  // Create initial array of processed patches,
  //  beginning by clipping each patch parameter space
  //  to a bounding box of its trimming curves
  // Then compute a bounding box of all the surfaces
  for(auto& surf : input_patches_view)
  {
    // This is where we would do Bezier extraction, if the curve-curve intersection
    //  routine were more robust :(
    //for(auto& bez : surf.extractTrimmedBezier())
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
    const double max_range_norm = bbox_threshold * total_bbox.range().norm();
    for(const auto& candidate : candidates)
    {
      if(candidate.boundingBox().range().norm() < max_range_norm)
      {
        new_bbox.addBox(candidate.boundingBox());
        subdivisions.push_back(candidate);
        continue;
      }

      NURBSType subpatches[2];
      candidate.nearBisectOnLongestAxis(subpatches[0], subpatches[1]);
      for(int si = 0; si < 2; si++)
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

}  // end namespace quest
}  // end namespace axom

#endif
