// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Delaunay.hpp
 *
 * \brief Declares the public `quest::Delaunay` incremental 2D/3D triangulation API.
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

#if defined(_MSC_VER)
  #define AXOM_QUEST_DELAUNAY_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
  #define AXOM_QUEST_DELAUNAY_FORCE_INLINE inline __attribute__((always_inline))
#else
  #define AXOM_QUEST_DELAUNAY_FORCE_INLINE inline
#endif

namespace axom
{
namespace quest
{

/**
   * \brief A class for incremental generation of a 2D or 3D Delaunay triangulation
   *
   * Construct a Delaunay triangulation incrementally by inserting points one by one.
   * The public API lives in this header while the larger insertion, point-location,
   * and validation routines are split into companion `detail/` headers.
   *
   * A bounding box of the points needs to be defined first via \a initializeBoundary(...).
   * The algorithm uses the Bowyer-Watson incremental insertion approach with robust
   * geometric predicates for handling degenerate cases including regular grids and
   * co-spherical point configurations.
   *
   * \note This class is not thread-safe. Multiple concurrent insertions are not supported.
   *
   * \tparam DIM The spatial dimension (2 for triangulation, 3 for tetrahedralization)
   */
template <int DIM = 2>
class Delaunay
{
public:
  AXOM_STATIC_ASSERT_MSG(DIM == 2 || DIM == 3, "The template parameter DIM can only be 2 or 3. ");

  /// Tolerance for barycentric coordinate comparisons (distinguishes interior from boundary)
  /// Value chosen to handle typical floating-point error accumulation in coordinate computation
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

  /// Controls the level of validation performed during point insertion
  enum class InsertionValidationMode
  {
    /// No additional insertion-time validation beyond existing asserts (production mode).
    None,
    /// Checks only insertion-local seed/cavity/ball invariants.
    Local,
    /// Checks local cavity/ball invariants and that the resulting IA mesh is conforming.
    /// Intended for debugging and unit tests.
    ConformingMesh,
    /// Additionally runs a global empty-circumsphere check after each insertion (very expensive).
    /// Use only for deep debugging; unsuitable for large meshes.
    Full
  };

  /// Result status from point location queries
  enum class PointLocationStatus
  {
    Found,    ///< Query point was successfully located inside an element
    Outside,  ///< Query point lies outside the triangulation boundary
    Failed    ///< Location failed (internal error or degenerate case)
  };

  /// Combined result from point location including element index and status
  struct PointLocationResult
  {
    IndexType element_idx {INVALID_INDEX};  ///< Index of containing element (or INVALID_INDEX)
    PointLocationStatus status {PointLocationStatus::Failed};  ///< Status of the location query
  };

  /// Statistics about point insertion operations (cavity size during Bowyer-Watson insertion)
  struct InsertionStats
  {
    std::uint64_t insertions {0};     ///< Total number of points inserted
    std::uint64_t total_removed {0};  ///< Cumulative elements removed across all insertions
    std::uint64_t max_removed {0};    ///< Maximum elements removed in a single insertion

    /// Average number of elements removed per insertion (cavity size)
    double mean_removed() const
    {
      return insertions > 0 ? static_cast<double>(total_removed) / static_cast<double>(insertions)
                            : 0.0;
    }
  };

  /// Statistics about point location query performance
  struct PointLocationStats
  {
    std::uint64_t walk_calls {0};    ///< Number of walk attempts
    std::uint64_t walk_found {0};    ///< Walks that successfully found containing element
    std::uint64_t walk_outside {0};  ///< Walks that terminated outside boundary
    std::uint64_t walk_failed {0};  ///< Walks that failed (should be rare with robust implementation)
    std::uint64_t total_walk_steps {0};  ///< Cumulative simplex-to-simplex steps across all walks
    std::uint64_t max_walk_steps {0};    ///< Maximum steps in a single walk

    /// Average number of simplex traversals per walk
    double mean_walk_steps() const
    {
      return walk_calls > 0 ? static_cast<double>(total_walk_steps) / static_cast<double>(walk_calls)
                            : 0.0;
    }
  };

  /// Precomputed circumsphere center and squared radius for in-sphere tests
  struct CircumsphereEval
  {
    PointType center {};    ///< Circumsphere center point
    double radius_sq {0.};  ///< Squared radius (avoids sqrt in distance comparisons)

    CircumsphereEval() = default;

    /// Construct from origin point and offset vector to center
    CircumsphereEval(const PointType& origin, const VectorType& center_offset)
      : center(origin + center_offset)
      , radius_sq(center_offset.squared_norm())
    { }
  };

  /// Orientation determinant evaluation with tolerance for robust geometric tests
  struct OrientationEval
  {
    double det {0.};  ///< Raw determinant value
    double tol {0.};  ///< Context-aware tolerance for this determinant
    int orientation {
      primal::ON_BOUNDARY};  ///< Classified orientation (ON_NEGATIVE_SIDE/ON_BOUNDARY/ON_POSITIVE_SIDE)
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

  /// Number of adjacency layers to search around visited simplices when a directed walk fails.
  /// If a walk cycles or reaches its step budget, we probe this many layers of neighbors
  /// before giving up. Value of 2 balances coverage vs. cost for typical meshes.
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

  // Scratch buffers used by point location to avoid per-call heap allocations.
  // Delaunay is not thread-safe, so these are safe to reuse between calls.
  mutable std::vector<IndexType> m_candidate_elements_scratch;
  mutable std::vector<IndexType> m_walked_elements_scratch;
  mutable std::vector<IndexType> m_walk_local_elements_scratch;
  mutable std::vector<IndexType> m_initial_vertices_scratch;
  std::unique_ptr<InsertionHelper> m_insertion_helper;

public:
  /**
   * \brief Default constructor
   * \note User must call initializeBoundary(BoundingBox) before adding points.
   */
  Delaunay() : m_has_boundary(false), m_insertion_validation_mode(InsertionValidationMode::None) { }

  /// \brief Returns statistics about insertion operations
  /// \return InsertionStats containing cavity size metrics
  InsertionStats getInsertionStats() const
  {
    return {m_num_insertions, m_total_removed_elements, m_max_removed_elements};
  }

  /// \brief Enable or disable collection of point location statistics
  /// \param enabled If true, track walk performance metrics (adds minimal overhead)
  /// \note Stats collection is disabled by default
  void setCollectPointLocationStats(bool enabled) { m_collect_location_stats = enabled; }

  /// \brief Returns statistics about point location query performance
  /// \return PointLocationStats containing walk performance metrics
  /// \note Returns zeros if setCollectPointLocationStats(true) was not called
  PointLocationStats getPointLocationStats() const
  {
    return {m_num_walk_calls,
            m_num_walk_found,
            m_num_walk_outside,
            m_num_walk_failed,
            m_total_walk_steps,
            m_max_walk_steps};
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
   * \details Subsequent points added to the triangulation must lie within this boundary.
   *          Creates an initial bounding box triangulation (2 triangles in 2D, 6 tetrahedra in 3D).
   * \param bb The bounding box that will contain all inserted points
   * \pre Must be called before any points are inserted
   */
  void initializeBoundary(const BoundingBox& bb);

  /// \brief Reserve storage for an expected number of inserted points.
  ///
  /// Uses a dimension-specific heuristic for the total simplex count so large
  /// bulk-builds can avoid repeated mesh-container reallocations.
  /// \param num_points Expected number of points to be inserted
  /// \note Based on Euler characteristic: 2D expects ~2n triangles for n points,
  ///       3D expects ~6n tetrahedra. Heuristics account for bounding box and over-tessellation.
  void reserveForPointCount(IndexType num_points)
  {
    if(!m_has_boundary || num_points <= 0)
    {
      return;
    }

    const IndexType expected_vertices = m_mesh.vertices().size() + num_points;

    // Heuristic element count: 2D triangle count ≈ 2n (Euler), 3D tet count ≈ 6n (Euler).
    // Use conservative estimates (3.0 and 9.0) to account for boundary and over-tessellation.
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
   * \brief Adds a new point and locally re-triangulates the mesh to ensure it remains Delaunay
   *
   * Uses the Bowyer-Watson incremental insertion algorithm:
   * 1. Locate the element containing the new point via directed walk
   * 2. Expand a cavity of all elements whose circumspheres contain the new point
   * 3. Retriangulate the cavity by connecting the new point to the cavity boundary
   *
   * \param new_pt The point to insert
   * \pre initializeBoundary() must have been called
   * \pre new_pt must lie within the bounding box
   * \pre The current mesh must be a valid Delaunay triangulation
   * \post The mesh remains a valid Delaunay triangulation
   * \note If insertion fails (rare), a warning is logged and the mesh is unchanged
   */
  void insertPoint(const PointType& new_pt);

  /// \brief Retrieves the geometric element (triangle or tetrahedron) at the given index
  /// \param element_index The index of the element to retrieve
  /// \return Triangle (2D) or Tetrahedron (3D) with vertex coordinates
  ElementType getElement(int element_index) const
  {
    const auto verts = m_mesh.boundaryVertices(element_index);
    if constexpr(DIM == 2)
    {
      return ElementType(m_mesh.getVertexPosition(verts[0]),
                         m_mesh.getVertexPosition(verts[1]),
                         m_mesh.getVertexPosition(verts[2]));
    }
    else
    {
      return ElementType(m_mesh.getVertexPosition(verts[0]),
                         m_mesh.getVertexPosition(verts[1]),
                         m_mesh.getVertexPosition(verts[2]),
                         m_mesh.getVertexPosition(verts[3]));
    }
  }

  /**
   * \brief Prints out mesh details, for debugging purpose.
   */
  void printMesh() { m_mesh.print_all(); }

  /**
   * \brief Write the mesh to a legacy VTK file for visualization
   *
   * \param filename The name of the file to write to
   * \note The suffix ".vtk" will be appended to the provided filename
   * \note This method compacts the mesh before export to eliminate deleted elements
   */
  void writeToVTKFile(const std::string& filename);

  /**
   * \brief Removes the bounding box vertices and all attached elements
   *
   * The initial bounding box (created by initializeBoundary) consists of 2^DIM vertices
   * that define a rectangular boundary. This method removes those vertices and all
   * simplices incident to them, leaving only the Delaunay triangulation of the inserted points.
   *
   * \post No more points can be inserted after calling this method
   * \post The mesh is compacted (deleted element slots are removed)
   */
  void removeBoundary();

  /// \brief Get the underlying IA mesh data structure
  /// \return Pointer to the IAMesh (for advanced users or ScatteredInterpolation)
  const IAMeshType* getMeshData() const { return &m_mesh; }

  /**
   * \brief Checks that the mesh satisfies the Delaunay empty-circumsphere property
   *
   * A Delaunay triangulation is valid when no vertex lies strictly inside the
   * circumsphere of any simplex. This method checks all vertex-simplex pairs.
   *
   * \param verboseOutput If true, prints detailed diagnostic information
   * \return true if the Delaunay property holds, false otherwise
   * \note This is an O(n·m) check where n = vertices, m = elements. Use for validation.
   */
  bool isValid(bool verboseOutput = false) const;

  /**
   * \brief Checks mesh conformity, orientation, and boundary consistency
   *
   * Verifies three properties:
   * 1. Topological conformity (manifold facets, reciprocal adjacencies) via IAMesh::isConforming()
   * 2. All simplices are positively oriented (positive signed volume/area)
   * 3. If boundary exists, all boundary facets lie on the bounding box faces
   *
   * \param verboseOutput If true, prints detailed diagnostic information
   * \return true if all checks pass, false otherwise
   */
  bool isConforming(bool verboseOutput = false) const;

  /// \brief Returns true if an element is active and can participate in point-location queries
  ///
  /// \param element_idx The element index to check
  /// \return true if element is valid (not a deleted/recycled slot)
  /// \note During insertion, all active elements have valid vertices. After removeBoundary(),
  ///       the mesh is compacted so all element indices are valid.
  bool isSearchableElement(IndexType element_idx) const;

  /// \brief Locate the element containing a query point using directed walk
  ///
  /// \param query_pt The point to locate
  /// \param warnOnInvalid If true, log a warning if location fails
  /// \return Element index containing the point, or INVALID_INDEX if not found
  /// \note Uses grid-seeded directed walk for O(n^(1/DIM)) expected performance
  IndexType findContainingElement(const PointType& query_pt, bool warnOnInvalid = true) const;

  /**
   * \brief Compute barycentric coordinates of a query point within an element
   *
   * \param element_idx The element containing (or near) the query point
   * \param q_pt The query point
   * \return Barycentric coordinates (DIM+1 values that sum to 1)
   * \note If the point is inside, all coordinates are non-negative (within tolerance)
   */
  BaryCoordType getBaryCoords(IndexType element_idx, const PointType& q_pt) const;

  /// \brief Returns initial seed elements for cavity expansion based on point location
  ///
  /// If the query point lies on a face/edge of the containing element (indicated by a
  /// barycentric coordinate near zero), the neighbor across that face is also included
  /// as a seed. This ensures the cavity expansion starts from all elements that might
  /// contain the point in their circumsphere.
  ///
  /// \param element_idx The element containing (or nearest to) the query point
  /// \param bary_coord Barycentric coordinates of the query point in that element
  /// \return Array of seed element indices (at least 1, possibly more if point is on boundary feature)
  IndexArray getSeedElements(IndexType element_idx, const BaryCoordType& bary_coord) const
  {
    IndexArray seed_elements;
    seed_elements.push_back(element_idx);
    constexpr IndexType invalid_element = IAMeshType::INVALID_ELEMENT_INDEX;

    // If query point lies on a facet (bary coord approximately 0), include the neighbor across that facet
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

  /// \brief Classify a value as positive, negative, or zero within tolerance
  /// \param value The value to classify
  /// \param tolerance The tolerance for considering the value as zero
  /// \return +1 if value > tolerance, -1 if value < -tolerance, 0 otherwise
  static int signWithTolerance(double value, double tolerance)
  {
    return value > tolerance ? 1 : (value < -tolerance ? -1 : 0);
  }

  static CircumsphereEval evaluateCircumsphereOnMesh(const IAMeshType& mesh, IndexType element_idx);

  static double sphereSquaredDistanceTolerance(const CircumsphereEval& sphere,
                                               const PointType& x,
                                               double distance_sq);

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

  /// \brief Scan a small adjacency region around a failed directed walk before falling back to a full scan
  IndexType findContainingElementFromNeighbors(const PointType& query_pt,
                                               const std::vector<IndexType>& seed_elements) const;

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
#include "detail/DelaunayValidation.hpp"
#include "detail/DelaunayImpl.hpp"

#undef AXOM_QUEST_DELAUNAY_FORCE_INLINE

#endif  // QUEST_DELAUNAY_H_
