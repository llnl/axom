// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Delaunay.hpp
 *
 * \brief Defines an incremental 2D/3D Delaunay triangulation.
 */

#ifndef QUEST_DELAUNAY_H_
#define QUEST_DELAUNAY_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "detail/DelaunayElementFinder.hpp"
#include "detail/DelaunayInsertionHelper.hpp"

#include "axom/fmt.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace axom
{
namespace quest
{

/**
 * \brief A class for incremental generation of a 2D or 3D Delaunay triangulation
 *
 * Construct a Delaunay triangulation incrementally by inserting points one by one.
 * A bounding box of the points needs to be defined first via \a initializeBoundary(...)
 */
template <int DIM = 2>
class Delaunay
{
public:
  AXOM_STATIC_ASSERT_MSG(DIM == 2 || DIM == 3, "The template parameter DIM can only be 2 or 3. ");
  static constexpr double BARY_EPS = 1e-12;

  using DataType = double;

  using PointType = primal::Point<DataType, DIM>;
  using VectorType = primal::Vector<DataType, DIM>;
  using ElementType =
    typename std::conditional<DIM == 2, primal::Triangle<DataType, 2>, primal::Tetrahedron<DataType, 3>>::type;
  using BaryCoordType = primal::Point<DataType, DIM + 1>;
  using BoundingBox = primal::BoundingBox<DataType, DIM>;

  using IAMeshType = slam::IAMesh<DIM, DIM, PointType>;
  using IndexType = typename IAMeshType::IndexType;
  using IndexArray = typename IAMeshType::IndexArray;

  using IndexPairType = std::pair<IndexType, IndexType>;

  static constexpr int VERT_PER_ELEMENT = DIM + 1;
  static constexpr int VERTS_PER_FACET = VERT_PER_ELEMENT - 1;
  static constexpr IndexType INVALID_INDEX = -1;

  enum class InsertionValidationMode
  {
    /// No additional insertion-time validation beyond existing asserts.
    None,
    /// Checks only insertion-local seed/cavity/ball invariants.
    Local,
    /// Checks local cavity/ball invariants and that the resulting IA mesh is conforming.
    /// Intended for debugging and unit tests.
    ConformingMesh,
    /// Additionally runs a global empty-circumsphere check after each insertion (very expensive).
    Full
  };

  enum class PointLocationStatus
  {
    Found,
    Outside,
    Failed
  };

  struct PointLocationResult
  {
    IndexType element_idx {INVALID_INDEX};
    PointLocationStatus status {PointLocationStatus::Failed};
  };

  struct InsertionStats
  {
    std::uint64_t insertions {0};
    std::uint64_t total_removed {0};
    std::uint64_t max_removed {0};

    double mean_removed() const
    {
      return insertions > 0 ? static_cast<double>(total_removed) / static_cast<double>(insertions)
                            : 0.0;
    }
  };

  struct PointLocationStats
  {
    std::uint64_t walk_calls {0};
    std::uint64_t walk_found {0};
    std::uint64_t walk_outside {0};
    std::uint64_t walk_failed {0};
    std::uint64_t total_walk_steps {0};
    std::uint64_t max_walk_steps {0};
    std::uint64_t linear_fallbacks {0};
    std::uint64_t empty_seed_fallbacks {0};

    double mean_walk_steps() const
    {
      return walk_calls > 0 ? static_cast<double>(total_walk_steps) / static_cast<double>(walk_calls)
                            : 0.0;
    }
  };

  struct CircumsphereEval
  {
    PointType center {};
    double radius_sq {0.};

    CircumsphereEval() = default;

    CircumsphereEval(const PointType& origin, const VectorType& center_offset)
      : center(origin + center_offset)
      , radius_sq(center_offset.squared_norm())
    { }
  };

  struct OrientationEval
  {
    double det {0.};
    double tol {0.};
    int orientation {primal::ON_BOUNDARY};
  };

private:
  using ModularFaceIndex =
    slam::ModularInt<slam::policies::CompileTimeSize<IndexType, VERT_PER_ELEMENT>>;

  using FacetKey = typename IAMeshType::FacetKey;
  using FacetRecord = typename IAMeshType::FacetRecord;

  using ElementFinder =
    detail::DelaunayElementFinder<DIM, PointType, IAMeshType, BoundingBox, IndexType>;
  using InsertionHelper =
    detail::DelaunayInsertionHelper<DIM, PointType, BaryCoordType, IndexType, IndexArray, IAMeshType>;

  // These broader fallbacks are only used for 3D query point location after the
  // initial directed walk fails. Insertions stay on the cheaper local path.
  static constexpr int QUERY_SEARCH_RADIUS = 6;
  static constexpr int QUERY_CANDIDATE_LIMIT = 128;
  static constexpr int WALK_NEIGHBORHOOD_LAYERS = 2;

  /**
   * \brief Lightweight LIFO pool of invalid simplex slots that can be reused
   * during cavity retriangulation before growing the mesh element array.
   */
  struct RecycledElementPool
  {
    void reserve(IndexType count) { m_slots.reserve(count); }

    void clear() { m_slots.clear(); }

    void release(IndexType element_idx) { m_slots.push_back(element_idx); }

    bool empty() const { return m_slots.empty(); }

    IndexType size() const { return m_slots.size(); }

    IndexType capacity() const { return m_slots.capacity(); }

    IndexType acquire()
    {
      SLIC_ASSERT(!m_slots.empty());
      const IndexType element_idx = m_slots.back();
      m_slots.resize(m_slots.size() - 1);
      return element_idx;
    }

  private:
    axom::Array<IndexType> m_slots;
  };

  IAMeshType m_mesh;
  BoundingBox m_bounding_box;
  bool m_has_boundary;
  InsertionValidationMode m_insertion_validation_mode;
  RecycledElementPool m_deleted_elements;

  ElementFinder m_element_finder;
  IndexType m_next_regrid_vertex_count {0};
  mutable slam::BitSet m_walk_visited;
  std::uint64_t m_total_removed_elements {0};
  std::uint64_t m_max_removed_elements {0};
  std::uint64_t m_num_insertions {0};
  bool m_collect_location_stats {false};
  mutable std::uint64_t m_num_walk_calls {0};
  mutable std::uint64_t m_num_walk_found {0};
  mutable std::uint64_t m_num_walk_outside {0};
  mutable std::uint64_t m_num_walk_failed {0};
  mutable std::uint64_t m_total_walk_steps {0};
  mutable std::uint64_t m_max_walk_steps {0};
  mutable std::uint64_t m_num_linear_fallbacks {0};
  mutable std::uint64_t m_num_empty_seed_fallbacks {0};

  // Scratch buffers used by point location to avoid per-call heap allocations.
  // Delaunay is not thread-safe, so these are safe to reuse between calls.
  mutable std::vector<IndexType> m_candidate_elements_scratch;
  mutable std::vector<IndexType> m_walked_elements_scratch;
  mutable std::vector<IndexType> m_walk_local_elements_scratch;
  mutable std::vector<IndexType> m_initial_vertices_scratch;
  mutable std::vector<IndexType> m_fallback_vertices_scratch;
  std::unique_ptr<InsertionHelper> m_insertion_helper;

public:
  /**
   * \brief Default constructor
   * \note User must call initializeBoundary(BoundingBox) before adding points.
   */
  Delaunay() : m_has_boundary(false), m_insertion_validation_mode(InsertionValidationMode::None) { }

  InsertionStats getInsertionStats() const
  {
    return {m_num_insertions, m_total_removed_elements, m_max_removed_elements};
  }

  void setCollectPointLocationStats(bool enabled) { m_collect_location_stats = enabled; }

  PointLocationStats getPointLocationStats() const
  {
    return {m_num_walk_calls,
            m_num_walk_found,
            m_num_walk_outside,
            m_num_walk_failed,
            m_total_walk_steps,
            m_max_walk_steps,
            m_num_linear_fallbacks,
            m_num_empty_seed_fallbacks};
  }

  /// \brief Controls the amount of validation performed around each point insertion
  ///
  /// \note This is intended for debugging. `InsertionValidationMode::Full` is a diagnostic mode
  /// and should not be enabled in performance-sensitive runs.
  void setInsertionValidationMode(InsertionValidationMode mode)
  {
    m_insertion_validation_mode = mode;
  }

  /// \brief Returns the current insertion validation mode
  InsertionValidationMode getInsertionValidationMode() const { return m_insertion_validation_mode; }

  /**
   * \brief Defines the boundary of the triangulation.
   * \details subsequent points added to the triangulation must not be outside of this boundary.
   */
  void initializeBoundary(const BoundingBox& bb);

  /// \brief Reserve storage for an expected number of inserted points.
  ///
  /// Uses a dimension-specific heuristic for the total simplex count so large
  /// bulk-builds can avoid repeated mesh-container reallocations.
  void reserveForPointCount(IndexType num_points)
  {
    if(!m_has_boundary || num_points <= 0)
    {
      return;
    }

    const IndexType expected_vertices = m_mesh.vertices().size() + num_points;
    constexpr double ELEMENTS_PER_POINT = DIM == 2 ? 3.0 : 9.0;
    const IndexType expected_elements = m_mesh.elements().size() +
      static_cast<IndexType>(std::ceil(ELEMENTS_PER_POINT * static_cast<double>(num_points)));

    m_mesh.reserveVertices(expected_vertices);
    m_mesh.reserveElements(expected_elements);
    if(m_walk_visited.size() < static_cast<int>(expected_elements))
    {
      m_walk_visited = slam::BitSet(static_cast<int>(expected_elements));
    }
    if(static_cast<IndexType>(m_deleted_elements.capacity()) < expected_elements / 8)
    {
      m_deleted_elements.reserve(expected_elements / 8);
    }
  }

  /**
   * \brief Adds a new point and locally re-triangulates the mesh to ensure that it stays Delaunay
   *
   * This function will traverse the mesh to find the element that contains
   * this point, creates the Delaunay cavity, which takes out all the elements
   * that contains the point in its sphere, and fill it with a Delaunay ball.
   *
   * \pre The current mesh must already be Delaunay.
   */
  void insertPoint(const PointType& new_pt);

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2, ElementType>::type getElement(int element_index) const
  {
    const auto verts = m_mesh.boundaryVertices(element_index);
    const PointType& p0 = m_mesh.getVertexPosition(verts[0]);
    const PointType& p1 = m_mesh.getVertexPosition(verts[1]);
    const PointType& p2 = m_mesh.getVertexPosition(verts[2]);

    return ElementType(p0, p1, p2);
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3, ElementType>::type getElement(int element_index) const
  {
    const auto verts = m_mesh.boundaryVertices(element_index);
    const PointType& p0 = m_mesh.getVertexPosition(verts[0]);
    const PointType& p1 = m_mesh.getVertexPosition(verts[1]);
    const PointType& p2 = m_mesh.getVertexPosition(verts[2]);
    const PointType& p3 = m_mesh.getVertexPosition(verts[3]);

    return ElementType(p0, p1, p2, p3);
  }

  /**
   * \brief Prints out mesh details, for debugging purpose.
   */
  void printMesh() { m_mesh.print_all(); }

  /**
   * \brief Write the m_mesh to a legacy VTK file
   *
   * \param filename The name of the file to write to,
   * \note The suffix ".vtk" will be appended to the provided filename
   * \details This function uses mint to write the m_mesh to VTK format.
   */
  void writeToVTKFile(const std::string& filename);

  /**
   * \brief Removes the vertices that defines the boundary of the mesh,
   * and the elements attached to them.
   *
   * \details After this function is called, no more points can be added to the m_mesh.
   */
  void removeBoundary();

  /// \brief Get the IA mesh data pointer
  const IAMeshType* getMeshData() const { return &m_mesh; }

  /**
   * \brief Checks that the underlying mesh is a valid Delaunay triangulation of the point set
   *
   * A Delaunay triangulation is valid when none of the vertices are inside the circumspheres
   * of any of the elements of the mesh
   */
  bool isValid(bool verboseOutput = false) const;

  /**
   * \brief Checks that the underlying mesh is conforming and consistently oriented
   *
   * Topological conformity (manifold facets and reciprocal adjacencies) is
   * delegated to `slam::IAMesh::isConforming()`. This routine additionally
   * verifies that the simplices remain positively oriented, and (when the
   * initial bounding-box boundary is still present) that boundary facets lie on
   * that bounding box.
   */
  bool isConforming(bool verboseOutput = false) const;

  /// \brief Returns true when an element slot is active and can participate in point-location
  ///
  /// Point-location only needs to reject tombstones here. During incremental
  /// insertion all active simplices retain valid vertices, and `removeBoundary()`
  /// compacts before post-build query use.
  bool isSearchableElement(IndexType element_idx) const;

  /// \brief Find the index of the element that contains the query point, or the element closest to the point.
  IndexType findContainingElement(const PointType& query_pt, bool warnOnInvalid = true) const;

  /**
   * \brief helper function to retrieve the barycentric coordinate of the query point in the element
   */
  BaryCoordType getBaryCoords(IndexType element_idx, const PointType& q_pt) const;

  /// \brief Returns cavity seed elements based on the simplex feature containing the query point
  IndexArray getSeedElements(IndexType element_idx, const BaryCoordType& bary_coord) const
  {
    IndexArray seed_elements;
    seed_elements.push_back(element_idx);
    constexpr IndexType invalid_element = IAMeshType::INVALID_ELEMENT_INDEX;

    for(int i = 0; i < VERT_PER_ELEMENT; ++i)
    {
      if(axom::utilities::abs(bary_coord[i]) <= BARY_EPS)
      {
        const IndexType nbr = m_mesh.adjacentElements(element_idx)[ModularFaceIndex(i) + 1];
        if(nbr != invalid_element)
        {
          SLIC_ASSERT(m_mesh.isValidElement(nbr));
          seed_elements.push_back(nbr);
        }
      }
    }

    return seed_elements;
  }

  static int signWithTolerance(double value, double tolerance)
  {
    return value > tolerance ? 1 : (value < -tolerance ? -1 : 0);
  }

  template <std::size_t N>
  static double getPointMagnitudeScale(const std::array<PointType, N>& pts)
  {
    double max_abs_coord = 1.;
    for(const auto& pt : pts)
    {
      for(int dim = 0; dim < DIM; ++dim)
      {
        max_abs_coord = axom::utilities::max(max_abs_coord, axom::utilities::abs(pt[dim]));
      }
    }

    return max_abs_coord;
  }

  static double orientationTolerance(const std::array<PointType, 4>& pts);

  static double determinant3(const PointType& p0, const PointType& p1, const PointType& p2);

  static double orientationDeterminant(const std::array<PointType, 4>& pts);

  static int symbolicOrientationSign(const std::array<PointType, 4>& pts,
                                     const std::array<IndexType, 4>& ranks);

  static CircumsphereEval evaluateCircumsphereOnMesh(const IAMeshType& mesh, IndexType element_idx);

  static double sphereSquaredDistanceTolerance(const CircumsphereEval& sphere,
                                               const PointType& x,
                                               double distance_sq);

  template <typename SphereType>
  static double sphereSignedDistanceTolerance(const SphereType& sphere, const PointType& x);

  static int inSphereOrientationOnMesh(const IAMeshType& mesh,
                                       const PointType& q,
                                       IndexType element_idx);

  static double inSphereDeterminantOnMesh(const IAMeshType& mesh,
                                          const PointType& q,
                                          IndexType element_idx);

  static const char* orientationResultName(int result);

  static bool isPointInSphereOnMesh(const IAMeshType& mesh,
                                    const PointType& q,
                                    IndexType element_idx,
                                    bool includeBoundary);

  static int classifyOrientationDeterminant(double det, double tol);

  OrientationEval evaluateElementOrientationDeterminant(IndexType element_idx) const;

private:
  template <typename FacetSubsetType>
  static FacetKey makeSortedFaceKey(const FacetSubsetType& facet);

  static std::string facetKeyString(const FacetKey& facet_key);

  double getBoundaryCoordinateTolerance() const;

  bool isFacetOnBoundingBox(const FacetKey& facet_key) const;

  double getElementSignedMeasure(IndexType element_idx) const;

  double getElementMeasureTolerance() const;

  /**
   * \brief Validates the point-location result and cavity seed selection before building the cavity
   *
   * This is a light-weight check meant to catch point-location failures early:
   * - the reported containing element must be searchable and contain the point
   * - the seed set must include that element (and may include face-adjacent neighbors)
   */
  void validateInsertionSeed(IndexType element_idx,
                             const PointType& query_pt,
                             const BaryCoordType& bary_coord,
                             const IndexArray& seed_elements) const;

  /**
   * \brief Validates that the cavity boundary facets match the faces between cavity and non-cavity elements
   *
   * The `InsertionHelper` collects a boundary facet for every cavity face that borders either:
   * - an element outside the cavity, or
   * - the temporary bounding-box boundary (invalid neighbor).
   */
  void validateCavityBoundary(const InsertionHelper& insertion_helper) const;

  /**
   * \brief Validates that the inserted ball covers the cavity boundary and is internally stitched
   *
   * Each inserted simplex must contain `new_pt_i`. Boundary facets must appear
   * exactly once among the inserted elements and point to the recorded outside
   * neighbor. Internal facets must be paired with reciprocal adjacency.
   */
  void validateInsertedBall(IndexType new_pt_i, const InsertionHelper& insertion_helper) const;

  /**
   * \brief Validates that the global mesh invariants still hold after insertion
   *
   * `ConformingMesh` checks IA validity and topological conformity, plus
   * Delaunay's geometry-specific checks (positive simplex orientation and
   * bounding-box boundary consistency while the fake boundary exists).
   * `Full` additionally runs the global Delaunay empty-circumsphere validation.
   */
  void validateInsertionResult() const;

  /// \brief Predicate for when to compact internal mesh data structures after removing elements
  bool shouldCompactMesh() const;

  /// \brief Compacts the underlying mesh
  void compactMesh();

  BaryCoordType getRawBarycentricDeterminants(IndexType element_idx, const PointType& query_pt) const;

  double rawBarycentricDeterminantTolerance(IndexType element_idx, const PointType& query_pt) const;

  bool isPointInsideForLocation(IndexType element_idx,
                                const PointType& query_pt,
                                const BaryCoordType& bary_coord,
                                ModularFaceIndex* exit_face = nullptr) const;

  PointLocationResult walkToContainingElement(
    const PointType& query_pt,
    IndexType start_element,
    std::vector<IndexType>* visited_elements_out = nullptr) const;

  void appendCandidateElement(std::vector<IndexType>& candidate_elements, IndexType vertex_i) const;

  void appendCandidateElementsFromVertices(std::vector<IndexType>& candidate_elements,
                                           const std::vector<IndexType>& candidate_vertices) const;

  void getInitialCandidateElements(const PointType& query_pt,
                                   std::vector<IndexType>& candidate_elements) const;

  PointLocationResult walkCandidateElements(const PointType& query_pt,
                                            const std::vector<IndexType>& candidate_elements,
                                            std::size_t start_idx = 0,
                                            std::vector<IndexType>* walked_elements = nullptr) const;

  PointLocationResult findContainingElementWithQueryFallbacks(
    const PointType& query_pt,
    std::vector<IndexType>& candidate_elements,
    const std::vector<IndexType>& walked_elements) const;

  /// \brief Scan a small adjacency region around a failed directed walk before falling back to a full scan
  IndexType findContainingElementFromNeighbors(const PointType& query_pt,
                                               const std::vector<IndexType>& seed_elements) const;

  /// \brief Last-resort exhaustive scan used when the cheaper local search path cannot classify the point
  IndexType findContainingElementLinear(const PointType& query_pt, bool warnOnInvalid) const;

  /// \brief Scan the stars of nearby seed vertices as a cheaper local fallback for query point location
  IndexType findContainingElementNearby(const PointType& query_pt,
                                        const std::vector<IndexType>& nearby_vertices) const;

  /**
   * \brief Helper function to fill the array with the initial mesh.
   * \details create a rectangle for 2D, cube for 3D, and fill the array with the mesh data.
   */
  void generateInitialMesh(std::vector<DataType>& points,
                           std::vector<IndexType>& elem,
                           const BoundingBox& bb);
};

template <int DIM>
constexpr typename Delaunay<DIM>::IndexType Delaunay<DIM>::INVALID_INDEX;

}  // namespace quest
}  // namespace axom

#include "detail/DelaunayPointLocation.hpp"
#include "detail/DelaunayImpl.hpp"

#endif  // QUEST_DELAUNAY_H_
