// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DELAUNAY_H_
#define QUEST_DELAUNAY_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "axom/fmt.hpp"

#include <list>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <array>
#include <memory>
#include <cstdlib>
#include <cmath>
#include <limits>

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

private:
  using ModularFaceIndex =
    slam::ModularInt<slam::policies::CompileTimeSize<IndexType, VERT_PER_ELEMENT>>;

  using FacetKey = typename IAMeshType::FacetKey;
  using FacetRecord = typename IAMeshType::FacetRecord;

private:
  struct ElementFinder;
  struct InsertionHelper;

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

  InsertionStats getInsertionStats() const
  {
    return {m_num_insertions, m_total_removed_elements, m_max_removed_elements};
  }

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

public:
  /**
   * \brief Default constructor
   * \note User must call initializeBoundary(BoundingBox) before adding points.
   */
  Delaunay() : m_has_boundary(false), m_insertion_validation_mode(InsertionValidationMode::None) { }

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
  void writeToVTKFile(const std::string& filename)
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

  /**
   * \brief Removes the vertices that defines the boundary of the mesh,
   * and the elements attached to them.
   *
   * \details After this function is called, no more points can be added to the m_mesh.
   */
  void removeBoundary()
  {
    if(m_has_boundary)
    {
      //remove the boundary box, which will be the first 4 points for triangles, first 8 for tetrahedron
      const int num_boundary_pts = 1 << DIM;

      // Remove all elements incident to boundary vertices. Avoid relying on
      // `vertexStar()` here since it may be incomplete when the mesh is not
      // manifold around the temporary boundary.
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

      // Defensive cleanup: ensure no valid element references a removed vertex.
      // This can happen if the boundary-vertex star is non-manifold and
      // `removeVertex()` cannot discover all incident elements via adjacency.
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

  /// \brief Get the IA mesh data pointer
  const IAMeshType* getMeshData() const { return &m_mesh; }

  /**
   * \brief Checks that the underlying mesh is a valid Delaunay triangulation of the point set
   *
   * A Delaunay triangulation is valid when none of the vertices are inside the circumspheres 
   * of any of the elements of the mesh
   */
  bool isValid(bool verboseOutput = false) const
  {
    // Implementation note: We use an UniformGrid spatial index to find the candidate elements
    // whose circumsphere might contain the vertices of the mesh
    // To build this faster, we bootstrap the UniformGrid with an ImplicitGrid

    using ImplicitGridType = spin::ImplicitGrid<DIM, axom::SEQ_EXEC, IndexType>;
    using UniformGridType = spin::UniformGrid<IndexType, DIM>;
    using NumericArrayType = NumericArray<DataType, DIM>;
    using axom::numerics::dot_product;

    bool valid = true;

    std::vector<std::pair<IndexType, IndexType>> invalidEntries;
    std::vector<IndexType> invalidElements;

    const IndexType totalVertices = m_mesh.vertices().size();
    const IndexType totalElements = m_mesh.elements().size();
    const IndexType res = axom::utilities::ceil(0.33 * std::pow(totalVertices, 1. / DIM));
    UniformGridType grid(m_bounding_box, NumericArray<int, DIM>(res).data());

    // An array to cache the circumspheres associated with each element
    axom::Array<typename ElementType::SphereType> circumspheres(totalElements);

    auto vertexInsideCircumsphere = [&](const PointType& vertex, IndexType element_idx) {
      // Use the same mesh-side geometric circumsphere classifier as insertion so
      // the empty-circumsphere validation follows the exact same inside/boundary
      // decisions. Boundary cases are treated as "not inside" for global checks.
      return isPointInSphereOnMesh(m_mesh, vertex, element_idx, /*includeBoundary=*/false);
    };

    auto circumsphereSignedDistanceTol = [](const typename ElementType::SphereType& sphere,
                                            const PointType& x) {
      const auto& center = sphere.getCenter();
      double scale = axom::utilities::max(1., sphere.getRadius());
      for(int dim = 0; dim < DIM; ++dim)
      {
        scale = axom::utilities::max(scale, axom::utilities::abs(center[dim]));
        scale = axom::utilities::max(scale, axom::utilities::abs(x[dim]));
      }

      // Signed distance to a sphere involves a subtraction after a norm. Use a
      // tolerance proportional to the coordinate/radius scale to give a
      // meaningful diagnostic in the verbose output.
      return 256. * std::numeric_limits<double>::epsilon() * scale;
    };

    // bootstrap the uniform grid using an implicit grid
    {
      using GridCell = typename ImplicitGridType::GridCell;

      // Add (bounding boxes of) element circumspheres to temporary implicit grid
      const auto resCell = GridCell(res);
      ImplicitGridType implicitGrid(m_bounding_box, &resCell, totalElements);
      for(auto element_idx : m_mesh.elements().positions())
      {
        if(m_mesh.isValidElement(element_idx))
        {
          const auto verts = m_mesh.boundaryVertices(element_idx);
          bool has_all_vertices = true;
          for(int i = 0; i < VERT_PER_ELEMENT; ++i)
          {
            has_all_vertices &= m_mesh.isValidVertex(verts[i]);
          }

          if(!has_all_vertices)
          {
            valid = false;
            if(verboseOutput)
            {
              invalidElements.push_back(element_idx);
            }
            continue;
          }

          circumspheres[element_idx] = this->getElement(element_idx).circumsphere();
          const auto& sphere = circumspheres[element_idx];
          const auto& center = sphere.getCenter().array();
          const auto offset = NumericArrayType(sphere.getRadius());

          BoundingBox bb;
          bb.addPoint(PointType(center - offset));
          bb.addPoint(PointType(center + offset));

          implicitGrid.insert(bb, element_idx);  // insert valid entries into grid
        }
      }

      // copy candidates from implicit grid directly into uniform grid
      const int kUpper = (DIM == 2) ? 0 : res;
      const IndexType stride[3] = {1, res, (DIM == 2) ? 0 : res * res};
      for(IndexType k = 0; k < kUpper; ++k)
      {
        for(IndexType j = 0; j < res; ++j)
        {
          for(IndexType i = 0; i < res; ++i)
          {
            const IndexType vals[3] = {i, j, k};
            const GridCell cell(vals);
            const auto idx = dot_product(cell.data(), stride, DIM);
            const auto binValues = implicitGrid.getCandidatesAsArray(cell);
            grid.getBinContents(idx).insert(0, binValues.size(), binValues.data());
          }
        }
      }
    }

    // for each vertex -- check in_sphere condition for candidate element
    for(auto vertex_idx : m_mesh.vertices().positions())
    {
      // skip if vertex at this index is not valid
      if(!m_mesh.isValidVertex(vertex_idx))
      {
        continue;
      }

      const auto& vertex = m_mesh.getVertexPosition(vertex_idx);
      for(const auto element_idx : grid.getBinContents(grid.getBinIndex(vertex)))
      {
        // no need to check for invalid elements -- only valid elements were added to grid

        // skip if this is a vertex of the element
        if(slam::is_subset(vertex_idx, m_mesh.boundaryVertices(element_idx)))
        {
          continue;
        }

        // check insphere condition
        if(vertexInsideCircumsphere(vertex, element_idx))
        {
          valid = false;

          if(verboseOutput)
          {
            invalidEntries.push_back(std::make_pair(vertex_idx, element_idx));
          }
        }
      }
    }

    if(verboseOutput)
    {
      if(valid)
      {
        SLIC_INFO("Delaunay complex was valid");
      }
      else
      {
        fmt::memory_buffer out;
        if(!invalidElements.empty())
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\t{} valid elements referenced invalid vertices: {}",
                         invalidElements.size(),
                         fmt::join(invalidElements, ", "));
        }
        for(const auto& pr : invalidEntries)
        {
          const auto vertex_idx = pr.first;
          const auto element_idx = pr.second;
          const auto& pos = m_mesh.getVertexPosition(vertex_idx);
          const auto element = this->getElement(element_idx);
          const auto circumsphere = element.circumsphere();
          const double tol = circumsphereSignedDistanceTol(circumsphere, pos);
          fmt::format_to(std::back_inserter(out),
                         "\n\tVertex {} @ {}"
                         "\n\tElement {}: {} w/ circumsphere: {}"
                         "\n\tDistance to circumcenter: {} (tol={})",
                         vertex_idx,
                         pos,
                         element_idx,
                         element,
                         circumsphere,
                         circumsphere.computeSignedDistance(pos),
                         tol);
        }

        SLIC_INFO(
          fmt::format("Delaunay complex was NOT valid. There were {} "
                      "vertices in the circumsphere of an element. {}",
                      invalidEntries.size(),
                      fmt::to_string(out)));
      }
    }

    return valid;
  }

  /**
   * \brief Checks that the underlying mesh is conforming and consistently oriented
   *
   * Topological conformity (manifold facets and reciprocal adjacencies) is
   * delegated to `slam::IAMesh::isConforming()`. This routine additionally
   * verifies that the simplices remain positively oriented, and (when the
   * initial bounding-box boundary is still present) that boundary facets lie on
   * that bounding box.
   */
  bool isConforming(bool verboseOutput = false) const
  {
    fmt::memory_buffer out;

    bool valid = m_mesh.isConforming(verboseOutput);

    // Geometry-specific checks that do not belong in slam::IAMesh.
    for(auto element_idx : m_mesh.elements().positions())
    {
      if(!m_mesh.isValidElement(element_idx))
      {
        continue;
      }

      // check orientations of top simplices
      const auto orient = evaluateElementOrientationDeterminant(element_idx);
      if(orient.orientation != primal::ON_POSITIVE_SIDE)
      {
        if(verboseOutput)
        {
          fmt::format_to(
            std::back_inserter(out),
            "\n\tElement {} has non-positive orientation determinant {:.17g} (tol={:.3g})",
            element_idx,
            orient.det,
            orient.tol);
        }
        valid = false;
      }

      // check that boundary elements are not on the Delaunay bounding box
      if(m_has_boundary)
      {
        const auto neighbors = m_mesh.adjacentElements(element_idx);
        for(int facet_idx = 0; facet_idx < VERT_PER_ELEMENT; ++facet_idx)
        {
          if(m_mesh.isValidElement(neighbors[facet_idx]))
          {
            continue;
          }

          const FacetKey facet_key = m_mesh.getSortedFacetKey(element_idx, facet_idx);
          if(!isFacetOnBoundingBox(facet_key))
          {
            if(verboseOutput)
            {
              fmt::format_to(std::back_inserter(out),
                             "\n\tBoundary facet {} is not on the initial bounding box",
                             facetKeyString(facet_key));
            }
            valid = false;
          }
        }
      }
    }

    if(verboseOutput)
    {
      if(valid)
      {
        SLIC_INFO("Delaunay mesh was conforming");
      }
      else
      {
        SLIC_INFO("Delaunay mesh was NOT conforming. Summary: " << fmt::to_string(out));
      }
    }

    return valid;
  }

private:
  template <typename FacetSubsetType>
  static FacetKey makeSortedFaceKey(const FacetSubsetType& facet)
  {
    // Canonical key for a simplex facet: store the vertex ids in sorted order
    // so the same topological facet is mapped identically regardless of
    // orientation or local indexing.
    FacetKey key {};
    for(int i = 0; i < VERTS_PER_FACET; ++i)
    {
      key[i] = facet[i];
    }
    std::sort(key.begin(), key.end());
    return key;
  }

  static std::string facetKeyString(const FacetKey& facet_key)
  {
    return fmt::format("[{}]", fmt::join(facet_key, ", "));
  }

  double getBoundaryCoordinateTolerance() const
  {
    const auto min_pt = m_bounding_box.getMin();
    const auto max_pt = m_bounding_box.getMax();

    double max_extent = 1.;
    for(int dim = 0; dim < DIM; ++dim)
    {
      max_extent = axom::utilities::max(max_extent, max_pt[dim] - min_pt[dim]);
    }

    return 64. * std::numeric_limits<double>::epsilon() * max_extent;
  }

  bool isFacetOnBoundingBox(const FacetKey& facet_key) const
  {
    // During incremental construction, the triangulation includes an initial
    // bounding box (or cube) to guarantee the domain is closed. When that
    // temporary boundary is present, any facet with an invalid neighbor must
    // lie on one of the bounding box planes within a small coordinate tolerance.
    const auto& min_pt = m_bounding_box.getMin();
    const auto& max_pt = m_bounding_box.getMax();
    const double tol = getBoundaryCoordinateTolerance();

    for(int dim = 0; dim < DIM; ++dim)
    {
      bool on_min_face = true;
      bool on_max_face = true;
      for(const IndexType vertex_idx : facet_key)
      {
        const auto& vertex = m_mesh.getVertexPosition(vertex_idx);
        on_min_face &= axom::utilities::abs(vertex[dim] - min_pt[dim]) <= tol;
        on_max_face &= axom::utilities::abs(vertex[dim] - max_pt[dim]) <= tol;
      }

      if(on_min_face || on_max_face)
      {
        return true;
      }
    }

    return false;
  }

  double getElementSignedMeasure(IndexType element_idx) const
  {
    if constexpr(DIM == 2)
    {
      return getElement(element_idx).signedArea();
    }
    else
    {
      return getElement(element_idx).signedVolume();
    }
  }

  double getElementMeasureTolerance() const
  {
    const double scale = axom::utilities::max(1., m_bounding_box.range().norm());
    return 256. * std::numeric_limits<double>::epsilon() * std::pow(scale, static_cast<double>(DIM));
  }

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
                             const IndexArray& seed_elements) const
  {
    if(m_insertion_validation_mode == InsertionValidationMode::None)
    {
      return;
    }

    fmt::memory_buffer out;
    bool valid = true;

    if(!isSearchableElement(element_idx))
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tContaining element {} is not searchable",
                     element_idx);
      valid = false;
    }
    else
    {
      const auto verts = m_mesh.boundaryVertices(element_idx);
      for(auto idx : verts.positions())
      {
        if(!m_mesh.isValidVertex(verts[idx]))
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\tContaining element {} references invalid vertex {}",
                         element_idx,
                         verts[idx]);
          valid = false;
          break;
        }
      }
    }

    if(!isPointInsideForLocation(element_idx, query_pt, bary_coord))
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tPoint {} is not inside containing element {}",
                     query_pt,
                     element_idx);
      valid = false;
    }

    if(seed_elements.empty())
    {
      fmt::format_to(std::back_inserter(out), "\n\tNo cavity seed elements were generated");
      valid = false;
    }

    bool found_containing_element = false;
    for(const IndexType seed_element : seed_elements)
    {
      found_containing_element |= (seed_element == element_idx);
      if(!m_mesh.isValidElement(seed_element))
      {
        fmt::format_to(std::back_inserter(out), "\n\tSeed element {} is invalid", seed_element);
        valid = false;
      }
    }

    if(!found_containing_element)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tContaining element {} is missing from the seed set",
                     element_idx);
      valid = false;
    }

    SLIC_ERROR_IF(!valid, "Delaunay insertion seed validation failed:" << fmt::to_string(out));
  }

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
   *  Delaunay's geometry-specific checks (positive simplex orientation and
   *  bounding-box boundary consistency while the fake boundary exists).
   * `Full` additionally runs the global Delaunay empty-circumsphere validation.
   */
  void validateInsertionResult() const
  {
    if(m_insertion_validation_mode == InsertionValidationMode::None)
    {
      return;
    }

    if(m_insertion_validation_mode == InsertionValidationMode::Local)
    {
      return;
    }

    // Note: These checks are intentionally global and can be expensive. They
    // are only enabled when the caller opts into insertion validation.
    if(!m_mesh.isValid(false))
    {
      m_mesh.isValid(true);
      SLIC_ERROR("Delaunay insertion produced an invalid IAMesh");
    }

    if(!isConforming(false))
    {
      isConforming(true);
      SLIC_ERROR("Delaunay insertion produced a non-conforming triangulation");
    }

    if(m_insertion_validation_mode == InsertionValidationMode::Full)
    {
      if(!isValid(false))
      {
        isValid(true);
        SLIC_ERROR(
          "Delaunay insertion produced a triangulation that violates the empty-circumsphere "
          "condition");
      }
    }
  }

public:
  /// \brief Returns true when an element slot is active and can participate in point-location
  ///
  /// Point-location only needs to reject tombstones here. During incremental
  /// insertion all active simplices retain valid vertices, and `removeBoundary()`
  /// compacts before post-build query use.
  bool isSearchableElement(IndexType element_idx) const
  {
    return m_mesh.isValidElement(element_idx);
  }

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

  // These broader fallbacks are only used for 3D query point location after the
  // initial directed walk fails. Insertions stay on the cheaper local path.
  static constexpr int QUERY_SEARCH_RADIUS = 6;
  static constexpr int QUERY_CANDIDATE_LIMIT = 128;
  static constexpr int WALK_NEIGHBORHOOD_LAYERS = 2;

  /// \brief Find the index of the element that contains the query point, or the element closest to the point.
  IndexType findContainingElement(const PointType& query_pt, bool warnOnInvalid = true) const
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
      SLIC_WARNING_IF(warnOnInvalid,
                      "Attempting to locate element at location outside valid domain");
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

    // Local recovery before falling back to a global scan. These steps are
    // intentionally conservative: they only succeed when they find a simplex
    // whose barycentric coordinates are non-negative (up to tolerance).
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

  static double orientationTolerance(const std::array<PointType, 4>& pts)
  {
    const double scale = getPointMagnitudeScale(pts);
    if constexpr(DIM == 2)
    {
      return 64. * std::numeric_limits<double>::epsilon() * scale * scale;
    }
    else
    {
      return 64. * std::numeric_limits<double>::epsilon() * scale * scale * scale;
    }
  }

  static double determinant3(const PointType& p0, const PointType& p1, const PointType& p2)
  {
    return axom::numerics::determinant(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
  }

  static double orientationDeterminant(const std::array<PointType, 4>& pts)
  {
    return axom::numerics::determinant(pts[0][0],
                                       pts[0][1],
                                       pts[0][2],
                                       1.,
                                       pts[1][0],
                                       pts[1][1],
                                       pts[1][2],
                                       1.,
                                       pts[2][0],
                                       pts[2][1],
                                       pts[2][2],
                                       1.,
                                       pts[3][0],
                                       pts[3][1],
                                       pts[3][2],
                                       1.);
  }

  static int symbolicOrientationSign(const std::array<PointType, 4>& pts,
                                     const std::array<IndexType, 4>& ranks)
  {
    const double det = orientationDeterminant(pts);
    const int det_sign = signWithTolerance(det, orientationTolerance(pts));
    if(det_sign != 0)
    {
      return det_sign;
    }

    const std::array<double, 4> cofactors {-determinant3(pts[1], pts[2], pts[3]),
                                           determinant3(pts[0], pts[2], pts[3]),
                                           -determinant3(pts[0], pts[1], pts[3]),
                                           determinant3(pts[0], pts[1], pts[2])};

    std::array<int, 4> order {{0, 1, 2, 3}};
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) { return ranks[lhs] < ranks[rhs]; });

    const double cofactor_tol = 64. * std::numeric_limits<double>::epsilon() *
      axom::utilities::max(1., orientationTolerance(pts));
    for(const int row : order)
    {
      const int sign = signWithTolerance(cofactors[row], cofactor_tol);
      if(sign != 0)
      {
        return sign;
      }
    }

    return 0;
  }

  //-----------------------------------------------------------------------------
  // In-sphere predicate helpers
  //
  // Delaunay intentionally classifies points against element circumspheres
  // using a shared geometric circumsphere classifier rather than primal's
  // determinant-based in-sphere classifier. This keeps cavity growth,
  // insertion validation, and the global empty-circumsphere check consistent
  // on ill-conditioned slivers while still letting the hot path avoid the
  // explicit sphere signed-distance sqrt. Raw determinant values are still
  // exposed separately for diagnostics.
  //-----------------------------------------------------------------------------
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

  static CircumsphereEval evaluateCircumsphereOnMesh(const IAMeshType& mesh, IndexType element_idx)
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

  static double sphereSquaredDistanceTolerance(const CircumsphereEval& sphere,
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

    // Since d^2 - r^2 = (d - r)(d + r), map the old signed-distance tolerance to
    // squared distance with a local bound on (d + r) rather than another global
    // coordinate-scale factor. This preserves the large-coordinate sliver cases
    // that motivated the geometric classifier in the first place.
    return 256. * std::numeric_limits<double>::epsilon() * scale * local_span;
  }

  template <typename SphereType>
  static double sphereSignedDistanceTolerance(const SphereType& sphere, const PointType& x)
  {
    const auto& center = sphere.getCenter();
    double scale = axom::utilities::max(1., sphere.getRadius());
    for(int dim = 0; dim < DIM; ++dim)
    {
      scale = axom::utilities::max(scale, axom::utilities::abs(center[dim]));
      scale = axom::utilities::max(scale, axom::utilities::abs(x[dim]));
    }

    return 256. * std::numeric_limits<double>::epsilon() * scale;
  }

  static int inSphereOrientationOnMesh(const IAMeshType& mesh, const PointType& q, IndexType element_idx)
  {
    const auto sphere = evaluateCircumsphereOnMesh(mesh, element_idx);
    const double distance_sq = VectorType(sphere.center, q).squared_norm();

    const double delta_sq = distance_sq - sphere.radius_sq;
    const double tol = sphereSquaredDistanceTolerance(sphere, q, distance_sq);
    if(axom::utilities::abs(delta_sq) <= tol)
    {
      return primal::ON_BOUNDARY;
    }

    return delta_sq < 0. ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
  }

  static double inSphereDeterminantOnMesh(const IAMeshType& mesh,
                                          const PointType& q,
                                          IndexType element_idx)
  {
    const auto verts = mesh.boundaryVertices(element_idx);

    if constexpr(DIM == 2)
    {
      return primal::in_sphere_determinant(q,
                                           mesh.getVertexPosition(verts[0]),
                                           mesh.getVertexPosition(verts[1]),
                                           mesh.getVertexPosition(verts[2]));
    }
    else
    {
      return primal::in_sphere_determinant(q,
                                           mesh.getVertexPosition(verts[0]),
                                           mesh.getVertexPosition(verts[1]),
                                           mesh.getVertexPosition(verts[2]),
                                           mesh.getVertexPosition(verts[3]));
    }
  }

  static const char* orientationResultName(int result)
  {
    switch(result)
    {
    case primal::ON_NEGATIVE_SIDE:
      return "inside";
    case primal::ON_BOUNDARY:
      return "boundary";
    case primal::ON_POSITIVE_SIDE:
      return "outside";
    default:
      return "unknown";
    }
  }

  static bool isPointInSphereOnMesh(const IAMeshType& mesh,
                                    const PointType& q,
                                    IndexType element_idx,
                                    bool includeBoundary)
  {
    const int res = inSphereOrientationOnMesh(mesh, q, element_idx);
    return includeBoundary ? (res != primal::ON_POSITIVE_SIDE) : (res == primal::ON_NEGATIVE_SIDE);
  }

  //-----------------------------------------------------------------------------
  // Simplex orientation helpers
  //
  // Validate that simplices are positively oriented using the same
  // determinant/tolerance pattern (determinants are 2*area in 2D and 6*volume in 3D).
  //-----------------------------------------------------------------------------
  struct OrientationEval
  {
    double det {0.};
    double tol {0.};
    int orientation {primal::ON_BOUNDARY};
  };

  static int classifyOrientationDeterminant(double det, double tol)
  {
    const int sign = signWithTolerance(det, tol);
    return sign < 0 ? primal::ON_NEGATIVE_SIDE
                    : (sign > 0 ? primal::ON_POSITIVE_SIDE : primal::ON_BOUNDARY);
  }

  OrientationEval evaluateElementOrientationDeterminant(IndexType element_idx) const
  {
    const double scale = (DIM == 2) ? 2. : 6.;
    const double det = scale * getElementSignedMeasure(element_idx);
    const double tol = scale * getElementMeasureTolerance();
    return {det, tol, classifyOrientationDeterminant(det, tol)};
  }

  BaryCoordType getRawBarycentricDeterminants(IndexType element_idx, const PointType& query_pt) const
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

  double rawBarycentricDeterminantTolerance(IndexType element_idx, const PointType& query_pt) const
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

  bool isPointInsideForLocation(IndexType element_idx,
                                const PointType& query_pt,
                                const BaryCoordType& bary_coord,
                                ModularFaceIndex* exit_face = nullptr) const
  {
    ModularFaceIndex min_face(bary_coord.array().argMin());

    // Fast path: if any barycentric coordinate is clearly negative, the point
    // lies outside the simplex across the most-negative face.
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

    // Ambiguous path: for near-zero barycentric coordinates, fall back to the
    // underlying (unnormalized) determinants to decide the sign consistently.
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

  /// \brief Walk from a starting element until the containing element is found or the walk cycles
  PointLocationResult walkToContainingElement(const PointType& query_pt,
                                              IndexType start_element,
                                              std::vector<IndexType>* visited_elements_out = nullptr) const
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
      m_max_walk_steps =
        axom::utilities::max(m_max_walk_steps, static_cast<std::uint64_t>(step_count));
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

  void appendCandidateElement(std::vector<IndexType>& candidate_elements, IndexType vertex_i) const
  {
    const IndexType element_i = m_mesh.coboundaryElement(vertex_i);
    if(isSearchableElement(element_i) &&
       std::find(candidate_elements.begin(), candidate_elements.end(), element_i) ==
         candidate_elements.end())
    {
      candidate_elements.push_back(element_i);
    }
  }

  void appendCandidateElementsFromVertices(std::vector<IndexType>& candidate_elements,
                                           const std::vector<IndexType>& candidate_vertices) const
  {
    for(const auto vertex_i : candidate_vertices)
    {
      appendCandidateElement(candidate_elements, vertex_i);
    }
  }

  void getInitialCandidateElements(const PointType& query_pt,
                                   std::vector<IndexType>& candidate_elements) const
  {
    candidate_elements.clear();
    candidate_elements.reserve(16);

    // Prefer a small set of vertices from nearby bins rather than a single bin
    // representative. For large meshes, a single cached vertex can be far (in
    // terms of simplex-to-simplex walks) from the query point even if it lies
    // in the same bin; providing a few local candidates keeps directed walks
    // short and avoids expensive fallbacks.
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

  PointLocationResult walkCandidateElements(const PointType& query_pt,
                                            const std::vector<IndexType>& candidate_elements,
                                            std::size_t start_idx = 0,
                                            std::vector<IndexType>* walked_elements = nullptr) const
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

  PointLocationResult findContainingElementWithQueryFallbacks(
    const PointType& query_pt,
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

  /// \brief Scan a small adjacency region around a failed directed walk before falling back to a full scan
  IndexType findContainingElementFromNeighbors(const PointType& query_pt,
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
         std::find(nearby_elements.begin(), nearby_elements.end(), element_idx) ==
           nearby_elements.end())
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
      // Keep this recovery step strict: it only succeeds on a true containing
      // element. Broader "closest element" behavior stays in the final linear
      // fallback so outside-hull queries can still exit cleanly.
      if(min_bary >= -BARY_EPS)
      {
        return element_idx;
      }
    }

    return INVALID_INDEX;
  }

  /// \brief Last-resort exhaustive scan used when the cheaper local search path cannot classify the point
  IndexType findContainingElementLinear(const PointType& query_pt, bool warnOnInvalid) const
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

  /// \brief Scan the stars of nearby seed vertices as a cheaper local fallback for query point location
  IndexType findContainingElementNearby(const PointType& query_pt,
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

private:
  /// \brief Predicate for when to compact internal mesh data structures after removing elements
  bool shouldCompactMesh() const
  {
    // Note: This auto-compacting feature is hard coded.
    // It may be good to let user have control of this option in the future.
    constexpr int MIN_REMOVED_ELEMENTS = DIM == 2 ? 512 : 2048;
    constexpr double REMOVED_ELEMENT_FRACTION = DIM == 2 ? 0.25 : 0.35;
    return static_cast<int>(m_deleted_elements.size()) > MIN_REMOVED_ELEMENTS &&
      (static_cast<double>(m_deleted_elements.size()) >
       REMOVED_ELEMENT_FRACTION * static_cast<double>(m_mesh.elements().size()));
  }

  /// \brief Compacts the underlying mesh
  void compactMesh()
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

  /**
   * \brief Helper function to fill the array with the initial mesh.
   * \details create a rectangle for 2D, cube for 3D, and fill the array with the mesh data.
   */
  void generateInitialMesh(std::vector<DataType>& points,
                           std::vector<IndexType>& elem,
                           const BoundingBox& bb);

private:
  /// Helper struct to find the first element near a point to be inserted
  struct ElementFinder
  {
    using NumericArrayType = NumericArray<IndexType, DIM>;
    using LatticeType = spin::RectangularLattice<DIM, double, IndexType>;

    explicit ElementFinder() = default;

    /**
     * \brief Resizes the grid and reinserts vertices
     *
     * Resizes using a heuristic based on the number of vertices in the mesh.
     */
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

      // regenerate lattice
      m_lattice = spin::rectangular_lattice_from_bounding_box(expandedBB, NumericArrayType(res));

      // resize m_bins
      resizeArray<DIM>(res);
      m_bins.fill(INVALID_INDEX);

      // insert vertices into lattice
      for(auto idx : verts.positions())
      {
        if(!mesh.isValidVertex(idx))
        {
          continue;
        }

        // Skip vertices that are no longer incident to any valid element (can
        // occur for interior vertices of removed cavities). Using them as
        // point-location seeds forces expensive global fallbacks.
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

    /**
     * \brief Returns the index of the vertex in the bin containing point \a pt
     *
     * \param pt The position in space of the vertex that we're checking
     * \note Some bins might not point to a vertex, so users should check
     * that the returned index is a valid vertex, e.g. using \a mesh.isValidVertex(vertex_id)
     */
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

    /// \brief Returns the index of the vertex in the bin containing point \a pt
    inline IndexType getNearbyVertex(const PointType& pt) const
    {
      const auto cell = m_lattice.gridCell(pt);
      return flatIndex(cell);
    }

    /// \brief Updates the cached value of the bin containing point \a pt to \a vertex_id
    inline void updateBin(const PointType& pt, IndexType vertex_id)
    {
      const auto cell = m_lattice.gridCell(pt);
      flatIndex(cell) = vertex_id;
    }

  private:
    /// Returns a reference to the index in the array for the ND point with grid index \a cell
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

    /// Dimension-specific helper for resizing the ND array in 2D
    template <int TDIM>
    typename std::enable_if<TDIM == 2, void>::type resizeArray(IndexType res)
    {
      m_bins.resize(res, res);
    }

    /// Dimension-specific helper for resizing the ND array in 3D
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

  /// Helper struct to locally insert a new point into a Delaunay complex while keeping the mesh Delaunay
  struct InsertionHelper
  {
  public:
    struct BoundaryFacet
    {
      std::array<IndexType, VERTS_PER_FACET> vertices {};
      IndexType neighbor {IAMeshType::ElementAdjacencyRelation::INVALID_INDEX};
    };

    InsertionHelper(IAMeshType& mesh) : m_mesh(mesh) { }

    void reset()
    {
      for(const IndexType element_idx : cavity_elems)
      {
        if(element_idx >= 0 && static_cast<std::size_t>(element_idx) < m_cavity_membership.size())
        {
          m_cavity_membership[static_cast<std::size_t>(element_idx)] = 0;
        }
      }

      boundary_facets.clear();
      cavity_elems.clear();
      inserted_elems.clear();
      containing_element = INVALID_INDEX;
      seed_elements_debug.clear();
      m_stack.clear();
    }

    /**
   * \brief Find the Delaunay cavity: the elements whose circumspheres contain the query point
   *
   * \details This function starts from an element \a element_i and searches through
   * neighboring elements for a list of element indices whose circumspheres
   * contain or touch the query point.
   * It also finds the faces on the boundaries of the cavity to help with filling the cavity
   * in the \a delaunayBall function.  
   *
   * \param query_pt the query point
   * \param element_i the element to start the search at
   */
    void findCavityElements(const PointType& query_pt, const IndexArray& seed_elements)
    {
      constexpr int reserveSize = (DIM == 2) ? 16 : 64;
      ensureCavityMembershipCapacity();

      if(m_stack.capacity() < reserveSize)
      {
        m_stack.reserve(reserveSize);
      }
      if(cavity_elems.capacity() < reserveSize)
      {
        cavity_elems.reserve(reserveSize);
      }
      if(boundary_facets.capacity() < reserveSize)
      {
        boundary_facets.reserve(reserveSize);
      }
      if(inserted_elems.capacity() < reserveSize)
      {
        inserted_elems.reserve(reserveSize);
      }
      constexpr IndexType invalid_element = IAMeshType::INVALID_ELEMENT_INDEX;

      // Seed the cavity with the containing element, and with any face-adjacent
      // neighbors when the insertion point lies on the containing simplex
      // boundary. This avoids repeatedly retriangulating structured inputs that
      // insert directly onto existing edges/faces.
      for(const IndexType element_i : seed_elements)
      {
        SLIC_ASSERT(m_mesh.isValidElement(element_i));
        addCavityElement(element_i);
      }

      while(!m_stack.empty())
      {
        const IndexType element_idx = m_stack.back();
        m_stack.pop_back();

        // Invariant: this element is valid, was checked and is in the cavity
        // Each neighbor is either in the cavity or the shared face is on the cavity boundary
        const auto neighbors = m_mesh.adjacentElements(element_idx);
        for(auto n_idx : neighbors.positions())
        {
          const IndexType nbr = neighbors[n_idx];

          // invalid neighbor means face is on domain boundary, and thus on cavity boundary
          if(nbr != invalid_element)
          {
            SLIC_ASSERT(m_mesh.isValidElement(nbr));

            // If the neighbor is already in the cavity, the shared face is
            // internal and not part of the cavity boundary.
            if(containsCavityElement(nbr))
            {
              continue;  // both elem and neighbor along face are in cavity
            }

            // neighbor is valid but not in cavity; check circumsphere and add when appropriate
            if(isPointInCircumsphere(query_pt, nbr))
            {
              addCavityElement(nbr);
              continue;  // face is internal to cavity, nothing left to do for this face
            }
          }

          // if we got here, the face is on the boundary of the Delaunay cavity
          // add it to the boundary facet list
          {
            const auto bdry = m_mesh.boundaryVertices(element_idx);

            typename IAMeshType::ModularVertexIndex mod_idx(n_idx);
            BoundaryFacet facet;
            for(int i = 0; i < VERTS_PER_FACET; i++)
            {
              facet.vertices[i] = bdry[mod_idx++];
            }
            //For tetrahedron, if the element face is odd, reverse vertex order
            if(DIM == 3 && n_idx % 2 == 1)
            {
              axom::utilities::swap(facet.vertices[1], facet.vertices[2]);
            }

            facet.neighbor = nbr;
            boundary_facets.push_back(facet);
          }
        }
      }

      SLIC_ASSERT_MSG(!cavity_elems.empty(), "Error: New point is not contained in the mesh");
      SLIC_ASSERT(!boundary_facets.empty());
    }

    /**
    * \brief Remove the elements in the Delaunay cavity
    */
    void createCavity(RecycledElementPool& deleted_elements)
    {
      for(const auto elem : cavity_elems)
      {
        m_mesh.removeElement(elem);
        deleted_elements.release(elem);
      }
    }

    /// \brief Fill in the Delaunay cavity with new elements containing the insertion point
    void delaunayBall(IndexType new_pt_i, RecycledElementPool& deleted_elements)
    {
      const int numFaces = static_cast<int>(boundary_facets.size());
      const IndexType invalid_neighbor = IAMeshType::ElementAdjacencyRelation::INVALID_INDEX;
      const PointType& new_pt = m_mesh.getVertexPosition(new_pt_i);

      IndexType vlist[VERT_PER_ELEMENT] {};
      IndexType neighbors[VERT_PER_ELEMENT] {};
      for(int i = 0; i < numFaces; ++i)
      {
        // Create a new element from the face and the inserted point
        for(int d = 0; d < VERTS_PER_FACET; ++d)
        {
          vlist[d] = boundary_facets[static_cast<std::size_t>(i)].vertices[d];
        }
        vlist[VERTS_PER_FACET] = new_pt_i;

        // Preserve positive simplex orientation regardless of the boundary-face
        // ordering used to describe the cavity facet.
        if constexpr(DIM == 2)
        {
          const PointType& p0 = m_mesh.getVertexPosition(vlist[0]);
          const PointType& p1 = m_mesh.getVertexPosition(vlist[1]);
          if(primal::orientation_determinant(p0, p1, new_pt) < 0.)
          {
            axom::utilities::swap(vlist[0], vlist[1]);
          }
        }
        else
        {
          const PointType& p0 = m_mesh.getVertexPosition(vlist[0]);
          const PointType& p1 = m_mesh.getVertexPosition(vlist[1]);
          const PointType& p2 = m_mesh.getVertexPosition(vlist[2]);
          if(primal::orientation_determinant(p0, p1, p2, new_pt) < 0.)
          {
            axom::utilities::swap(vlist[1], vlist[2]);
          }
        }

        // Face 0 is the cavity boundary face opposite the inserted point.
        // The remaining faces stay invalid until fixVertexNeighborhood() stitches
        // the new ball together around the inserted vertex.
        const auto nID = boundary_facets[static_cast<std::size_t>(i)].neighbor;
        for(int d = 0; d < VERT_PER_ELEMENT; ++d)
        {
          neighbors[d] = invalid_neighbor;
        }
        neighbors[0] = nID;

        IndexType new_el = invalid_neighbor;
        if(!deleted_elements.empty())
        {
          new_el = deleted_elements.acquire();
          m_mesh.reuseElement(new_el, vlist, neighbors);
        }
        else
        {
          new_el = m_mesh.addElement(vlist, neighbors);
        }
        inserted_elems.push_back(new_el);
      }

      // Fix neighborhood around the new point
      m_mesh.fixVertexNeighborhood(new_pt_i, inserted_elems);
    }

    /// \brief Returns the number of elements removed during this insertion
    int numRemovedElements() const { return static_cast<int>(cavity_elems.size()); }

    bool containsCavityElement(IndexType element_idx) const
    {
      return element_idx >= 0 && static_cast<std::size_t>(element_idx) < m_cavity_membership.size() &&
        m_cavity_membership[static_cast<std::size_t>(element_idx)] != 0;
    }

    /// \brief Returns true when the query point is inside or on the element circumsphere
    bool isPointInCircumsphere(const PointType& query_pt, IndexType element_idx) const;

  public:
    IAMeshType& m_mesh;

    std::vector<BoundaryFacet> boundary_facets;
    std::vector<IndexType> cavity_elems;
    std::vector<IndexType> inserted_elems;
    IndexType containing_element {INVALID_INDEX};
    BaryCoordType containing_bary;
    IndexArray seed_elements_debug;

    IndexArray m_stack;

  private:
    void ensureCavityMembershipCapacity()
    {
      const std::size_t required_size = static_cast<std::size_t>(m_mesh.elements().size());
      if(m_cavity_membership.size() < required_size)
      {
        m_cavity_membership.resize(required_size, 0);
      }
    }

    void addCavityElement(IndexType element_idx)
    {
      ensureCavityMembershipCapacity();
      if(containsCavityElement(element_idx))
      {
        return;
      }

      m_cavity_membership[static_cast<std::size_t>(element_idx)] = 1;
      cavity_elems.push_back(element_idx);
      m_stack.push_back(element_idx);
    }

    std::vector<unsigned char> m_cavity_membership;
  };
};

template <int DIM>
constexpr typename Delaunay<DIM>::IndexType Delaunay<DIM>::INVALID_INDEX;

template <int DIM>
void Delaunay<DIM>::initializeBoundary(const BoundingBox& bb)
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
  m_num_linear_fallbacks = 0;
  m_num_empty_seed_fallbacks = 0;

  m_candidate_elements_scratch.clear();
  m_candidate_elements_scratch.reserve(QUERY_CANDIDATE_LIMIT);
  m_walked_elements_scratch.clear();
  m_walked_elements_scratch.reserve(256);
  m_walk_local_elements_scratch.clear();
  m_walk_local_elements_scratch.reserve(256);
  m_initial_vertices_scratch.clear();
  m_initial_vertices_scratch.reserve(8);
  m_fallback_vertices_scratch.clear();
  m_fallback_vertices_scratch.reserve(QUERY_CANDIDATE_LIMIT);
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
void Delaunay<DIM>::insertPoint(const PointType& new_pt)
{
  //Make sure initializeBoundary(...) is called first
  SLIC_ASSERT_MSG(m_has_boundary, "Error: Need a predefined boundary box prior to adding points.");
  SLIC_ASSERT_MSG(m_insertion_helper != nullptr,
                  "Error: Insertion helper was not initialized. "
                  "Delaunay::initializeBoundary() needs to be called first.");

  //Make sure the new point is inside the boundary box
  SLIC_ASSERT_MSG(m_bounding_box.contains(new_pt),
                  "Error: new point is outside of the boundary box.");

  // Find the mesh element containing the insertion point
  IndexType element_i = findContainingElement(new_pt);

  if(element_i == INVALID_INDEX)
  {
    SLIC_WARNING(
      fmt::format("Could not insert point {} into Delaunay triangulation: "
                  "Element containing that point was not found",
                  new_pt));
    return;
  }

  // Run the insertion operation by finding invalidated elements around the point (the "cavity")
  // and replacing them with new valid elements (the Delaunay "ball")
  const BaryCoordType bary_coord = getBaryCoords(element_i, new_pt);
  const IndexArray seed_elements = getSeedElements(element_i, bary_coord);
  validateInsertionSeed(element_i, new_pt, bary_coord, seed_elements);

  auto& insertionHelper = *m_insertion_helper;
  insertionHelper.reset();
  insertionHelper.containing_element = element_i;
  insertionHelper.containing_bary = bary_coord;
  insertionHelper.seed_elements_debug = seed_elements;
  insertionHelper.findCavityElements(new_pt, seed_elements);
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

  // Compact the mesh if there are too many removed elements
  if(shouldCompactMesh())
  {
    this->compactMesh();
  }
}

template <int DIM>
void Delaunay<DIM>::validateCavityBoundary(const InsertionHelper& insertion_helper) const
{
  if(m_insertion_validation_mode == InsertionValidationMode::None)
  {
    return;
  }

  // Cavity boundary invariant:
  // - Every cavity face that borders a non-cavity neighbor (or the temporary
  //   bounding-box boundary) must appear exactly once in `boundary_facets`.
  // - The facet's recorded neighbor must match the mesh adjacency.
  struct FacetInfo
  {
    IndexType neighbor_idx {INVALID_INDEX};
    bool matched {false};
  };

  fmt::memory_buffer out;
  bool valid = true;
  std::map<FacetKey, FacetInfo> cavity_boundary;

  for(const auto& facet : insertion_helper.boundary_facets)
  {
    const FacetKey facet_key = makeSortedFaceKey(facet.vertices);
    const auto insert_status = cavity_boundary.insert({facet_key, {facet.neighbor, false}});
    if(!insert_status.second)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tCavity boundary facet {} was recorded more than once",
                     facetKeyString(facet_key));
      valid = false;
    }
  }

  for(const IndexType cavity_element : insertion_helper.cavity_elems)
  {
    const auto neighbors = m_mesh.adjacentElements(cavity_element);
    for(int facet_idx = 0; facet_idx < VERT_PER_ELEMENT; ++facet_idx)
    {
      const IndexType neighbor_idx = neighbors[facet_idx];
      if(m_mesh.isValidElement(neighbor_idx) && insertion_helper.containsCavityElement(neighbor_idx))
      {
        continue;
      }

      const FacetKey facet_key = m_mesh.getSortedFacetKey(cavity_element, facet_idx);
      auto facet_it = cavity_boundary.find(facet_key);
      if(facet_it == cavity_boundary.end())
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCavity facet {} on element {} face {} is missing from the "
                       "boundary facet set",
                       facetKeyString(facet_key),
                       cavity_element,
                       facet_idx);
        valid = false;
        continue;
      }

      if(facet_it->second.neighbor_idx != neighbor_idx)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCavity face {} expects neighbor {} but facet relation stores {}",
                       facetKeyString(facet_key),
                       neighbor_idx,
                       facet_it->second.neighbor_idx);
        valid = false;
      }

      if(!m_mesh.isValidElement(neighbor_idx) && m_has_boundary && !isFacetOnBoundingBox(facet_key))
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCavity boundary face {} exits the mesh away from the bounding box",
                       facetKeyString(facet_key));
        valid = false;
      }

      facet_it->second.matched = true;
    }
  }

  for(const auto& facet_entry : cavity_boundary)
  {
    if(!facet_entry.second.matched)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tFacet {} does not correspond to a cavity boundary face",
                     facetKeyString(facet_entry.first));
      valid = false;
    }
  }

  SLIC_ERROR_IF(!valid,
                "Delaunay cavity validation failed after insertion " << m_num_insertions << ":"
                                                                     << fmt::to_string(out));
}

template <int DIM>
void Delaunay<DIM>::validateInsertedBall(IndexType new_pt_i,
                                         const InsertionHelper& insertion_helper) const
{
  if(m_insertion_validation_mode == InsertionValidationMode::None)
  {
    return;
  }

  // Ball invariant:
  // - Every inserted element contains `new_pt_i` and is positively oriented.
  // - Each cavity boundary facet is covered by exactly one inserted element and
  //   points to the recorded outside neighbor.
  // - All remaining facets are internal to the ball and are paired with reciprocal adjacencies.
  struct BoundaryFacetInfo
  {
    IndexType neighbor_idx {INVALID_INDEX};
    bool matched {false};
  };

  fmt::memory_buffer out;
  bool valid = true;
  std::map<FacetKey, BoundaryFacetInfo> cavity_boundary;
  std::map<FacetKey, std::vector<FacetRecord>> inserted_faces;
  auto findOppositeVertex = [&](IndexType element_idx, const FacetKey& facet_key) {
    const auto verts = m_mesh.boundaryVertices(element_idx);
    for(int i = 0; i < VERT_PER_ELEMENT; ++i)
    {
      bool on_face = false;
      for(int j = 0; j < VERTS_PER_FACET; ++j)
      {
        on_face |= (verts[i] == facet_key[j]);
      }
      if(!on_face)
      {
        return verts[i];
      }
    }

    return INVALID_INDEX;
  };

  if(insertion_helper.inserted_elems.size() != insertion_helper.boundary_facets.size())
  {
    fmt::format_to(
      std::back_inserter(out),
      "\n\tInserted ball element count {} does not match cavity boundary facet count {}",
      insertion_helper.inserted_elems.size(),
      insertion_helper.boundary_facets.size());
    valid = false;
  }

  if(insertion_helper.containing_element != INVALID_INDEX)
  {
    fmt::format_to(std::back_inserter(out),
                   "\n\tContaining element {} barycentric coordinates {}",
                   insertion_helper.containing_element,
                   insertion_helper.containing_bary);
    if(!insertion_helper.seed_elements_debug.empty())
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tInsertion seeds: [{}]",
                     fmt::join(insertion_helper.seed_elements_debug, ", "));
    }
  }

  for(const auto& facet : insertion_helper.boundary_facets)
  {
    cavity_boundary.insert({makeSortedFaceKey(facet.vertices), {facet.neighbor, false}});
  }

  for(const IndexType element_idx : insertion_helper.inserted_elems)
  {
    if(!m_mesh.isValidElement(element_idx))
    {
      fmt::format_to(std::back_inserter(out), "\n\tInserted element {} is invalid", element_idx);
      valid = false;
      continue;
    }

    const auto verts = m_mesh.boundaryVertices(element_idx);
    if(!slam::is_subset(new_pt_i, verts))
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tInserted element {} does not contain the new vertex {}",
                     element_idx,
                     new_pt_i);
      valid = false;
    }

    const auto orient = evaluateElementOrientationDeterminant(element_idx);
    if(orient.orientation != primal::ON_POSITIVE_SIDE)
    {
      fmt::format_to(
        std::back_inserter(out),
        "\n\tInserted element {} has non-positive orientation determinant {:.17g} (tol={:.3g})",
        element_idx,
        orient.det,
        orient.tol);
      valid = false;
    }

    const auto neighbors = m_mesh.adjacentElements(element_idx);
    for(int facet_idx = 0; facet_idx < VERT_PER_ELEMENT; ++facet_idx)
    {
      inserted_faces[m_mesh.getSortedFacetKey(element_idx, facet_idx)].push_back(
        {element_idx, facet_idx, neighbors[facet_idx]});
    }
  }

  for(auto& face_entry : inserted_faces)
  {
    auto cavity_it = cavity_boundary.find(face_entry.first);
    auto& records = face_entry.second;
    if(cavity_it != cavity_boundary.end())
    {
      if(records.size() != 1)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tBoundary face {} of the inserted ball is used by {} new elements",
                       facetKeyString(face_entry.first),
                       records.size());
        valid = false;
        continue;
      }

      if(records.front().neighbor_idx != cavity_it->second.neighbor_idx)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tInserted boundary face {} points to neighbor {} instead of {}",
                       facetKeyString(face_entry.first),
                       records.front().neighbor_idx,
                       cavity_it->second.neighbor_idx);
        valid = false;
      }

      if(m_mesh.isValidElement(cavity_it->second.neighbor_idx))
      {
        const IndexType inserted_opposite =
          findOppositeVertex(records.front().element_idx, face_entry.first);
        const IndexType neighbor_opposite =
          findOppositeVertex(cavity_it->second.neighbor_idx, face_entry.first);

        if(inserted_opposite == INVALID_INDEX || neighbor_opposite == INVALID_INDEX)
        {
          fmt::format_to(
            std::back_inserter(out),
            "\n\tCould not identify opposite vertices across inserted boundary face {}",
            facetKeyString(face_entry.first));
          valid = false;
        }
        else
        {
          const PointType& inserted_point = m_mesh.getVertexPosition(inserted_opposite);
          const PointType& neighbor_point = m_mesh.getVertexPosition(neighbor_opposite);

          if(isPointInSphereOnMesh(m_mesh,
                                   inserted_point,
                                   cavity_it->second.neighbor_idx,
                                   /*includeBoundary=*/false))
          {
            const int query_in_neighbor =
              inSphereOrientationOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx);
            const int inserted_in_neighbor =
              inSphereOrientationOnMesh(m_mesh, inserted_point, cavity_it->second.neighbor_idx);
            fmt::format_to(
              std::back_inserter(out),
              "\n\tInserted boundary face {} leaves new vertex {} inside neighbor {} "
              "circumsphere (query {}={}, det={:.17g}; opposite {}={}, det={:.17g})",
              facetKeyString(face_entry.first),
              inserted_opposite,
              cavity_it->second.neighbor_idx,
              new_pt_i,
              orientationResultName(query_in_neighbor),
              inSphereDeterminantOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx),
              inserted_opposite,
              orientationResultName(inserted_in_neighbor),
              inSphereDeterminantOnMesh(m_mesh, inserted_point, cavity_it->second.neighbor_idx));
            valid = false;
          }

          if(isPointInSphereOnMesh(m_mesh,
                                   neighbor_point,
                                   records.front().element_idx,
                                   /*includeBoundary=*/false))
          {
            const int query_in_neighbor =
              inSphereOrientationOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx);
            const int neighbor_in_inserted =
              inSphereOrientationOnMesh(m_mesh, neighbor_point, records.front().element_idx);
            fmt::format_to(
              std::back_inserter(out),
              "\n\tInserted boundary face {} leaves neighbor vertex {} inside new element {} "
              "circumsphere (query {} in neighbor {} is {}, det={:.17g}; opposite {} in new "
              "element is {}, det={:.17g})",
              facetKeyString(face_entry.first),
              neighbor_opposite,
              records.front().element_idx,
              new_pt_i,
              cavity_it->second.neighbor_idx,
              orientationResultName(query_in_neighbor),
              inSphereDeterminantOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx),
              neighbor_opposite,
              orientationResultName(neighbor_in_inserted),
              inSphereDeterminantOnMesh(m_mesh, neighbor_point, records.front().element_idx));
            valid = false;
          }
        }
      }

      cavity_it->second.matched = true;
    }
    else if(records.size() != 2 || records[0].neighbor_idx != records[1].element_idx ||
            records[1].neighbor_idx != records[0].element_idx)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tInternal face {} of the inserted ball has inconsistent adjacency",
                     facetKeyString(face_entry.first));
      valid = false;
    }
    else
    {
      const IndexType lhs_opposite = findOppositeVertex(records[0].element_idx, face_entry.first);
      const IndexType rhs_opposite = findOppositeVertex(records[1].element_idx, face_entry.first);

      if(lhs_opposite == INVALID_INDEX || rhs_opposite == INVALID_INDEX)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCould not identify opposite vertices across inserted internal face {}",
                       facetKeyString(face_entry.first));
        valid = false;
      }
      else
      {
        const PointType& lhs_point = m_mesh.getVertexPosition(lhs_opposite);
        const PointType& rhs_point = m_mesh.getVertexPosition(rhs_opposite);

        if(isPointInSphereOnMesh(m_mesh,
                                 lhs_point,
                                 records[1].element_idx,
                                 /*includeBoundary=*/false))
        {
          const int query_in_rhs = inSphereOrientationOnMesh(m_mesh,
                                                             m_mesh.getVertexPosition(new_pt_i),
                                                             records[1].element_idx);
          const int lhs_in_rhs = inSphereOrientationOnMesh(m_mesh, lhs_point, records[1].element_idx);
          fmt::format_to(
            std::back_inserter(out),
            "\n\tInserted internal face {} leaves vertex {} inside adjacent element {} "
            "circumsphere (query {}={}, det={:.17g}; opposite {}={}, det={:.17g})",
            facetKeyString(face_entry.first),
            lhs_opposite,
            records[1].element_idx,
            new_pt_i,
            orientationResultName(query_in_rhs),
            inSphereDeterminantOnMesh(m_mesh,
                                      m_mesh.getVertexPosition(new_pt_i),
                                      records[1].element_idx),
            lhs_opposite,
            orientationResultName(lhs_in_rhs),
            inSphereDeterminantOnMesh(m_mesh, lhs_point, records[1].element_idx));
          valid = false;
        }

        if(isPointInSphereOnMesh(m_mesh,
                                 rhs_point,
                                 records[0].element_idx,
                                 /*includeBoundary=*/false))
        {
          const int query_in_lhs = inSphereOrientationOnMesh(m_mesh,
                                                             m_mesh.getVertexPosition(new_pt_i),
                                                             records[0].element_idx);
          const int rhs_in_lhs = inSphereOrientationOnMesh(m_mesh, rhs_point, records[0].element_idx);
          fmt::format_to(
            std::back_inserter(out),
            "\n\tInserted internal face {} leaves vertex {} inside adjacent element {} "
            "circumsphere (query {}={}, det={:.17g}; opposite {}={}, det={:.17g})",
            facetKeyString(face_entry.first),
            rhs_opposite,
            records[0].element_idx,
            new_pt_i,
            orientationResultName(query_in_lhs),
            inSphereDeterminantOnMesh(m_mesh,
                                      m_mesh.getVertexPosition(new_pt_i),
                                      records[0].element_idx),
            rhs_opposite,
            orientationResultName(rhs_in_lhs),
            inSphereDeterminantOnMesh(m_mesh, rhs_point, records[0].element_idx));
          valid = false;
        }
      }
    }
  }

  for(const auto& facet_entry : cavity_boundary)
  {
    if(!facet_entry.second.matched)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tCavity boundary face {} was not covered by the inserted ball",
                     facetKeyString(facet_entry.first));
      valid = false;
    }
  }

  SLIC_ERROR_IF(!valid,
                "Delaunay ball validation failed after insertion " << m_num_insertions << ":"
                                                                   << fmt::to_string(out));
}

//--------------------------------------------------------------------------------
// Below are 2D and 3D specializations for methods in the Delaunay class
//--------------------------------------------------------------------------------

// 2D specialization for generateInitialMesh(...)
template <>
inline void Delaunay<2>::generateInitialMesh(std::vector<DataType>& points,
                                             std::vector<IndexType>& elem,
                                             const BoundingBox& bb)
{
  //Set up the initial IA mesh of 2 triangles forming a rectangle

  const PointType& mins = bb.getMin();
  const PointType& maxs = bb.getMax();

  // clang-format off
  std::vector<DataType> pt { mins[0], mins[1], 
                             mins[0], maxs[1], 
                             maxs[0], mins[1], 
                             maxs[0], maxs[1] };

  std::vector<IndexType> el { 0, 2, 1, 
                              3, 1, 2 };
  // clang-format on

  points.swap(pt);
  elem.swap(el);
}

// 3D specialization for generateInitialMesh(...)
template <>
inline void Delaunay<3>::generateInitialMesh(std::vector<DataType>& points,
                                             std::vector<IndexType>& elem,
                                             const BoundingBox& bb)
{
  //Set up the initial IA mesh of 6 tetrahedrons forming a cube
  const PointType& mins = bb.getMin();
  const PointType& maxs = bb.getMax();

  // clang-format off
  std::vector<DataType> pt { mins[0], mins[1], mins[2], 
                             mins[0], mins[1], maxs[2], 
                             mins[0], maxs[1], mins[2], 
                             mins[0], maxs[1], maxs[2],
                             maxs[0], mins[1], mins[2], 
                             maxs[0], mins[1], maxs[2], 
                             maxs[0], maxs[1], mins[2], 
                             maxs[0], maxs[1], maxs[2] };

  std::vector<IndexType> el { 3, 2, 4, 0, 
                              3, 4, 1, 0, 
                              3, 2, 6, 4,
                              3, 6, 7, 4, 
                              3, 5, 1, 4, 
                              3, 7, 5, 4 };
  // clang-format on

  points.swap(pt);
  elem.swap(el);
}

template <int DIM>
inline typename Delaunay<DIM>::BaryCoordType Delaunay<DIM>::getBaryCoords(IndexType element_idx,
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
inline bool Delaunay<DIM>::InsertionHelper::isPointInCircumsphere(const PointType& query_pt,
                                                                  IndexType element_idx) const
{
  // The cavity is defined by elements whose circumspheres contain or touch the
  // insertion point. Returning "true on boundary" ensures the cavity is
  // topologically closed for co-spherical inputs (e.g. regular grids).
  return Delaunay<DIM>::isPointInSphereOnMesh(m_mesh, query_pt, element_idx, /*includeBoundary=*/true);
}

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
