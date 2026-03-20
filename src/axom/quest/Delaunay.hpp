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

  IAMeshType m_mesh;
  BoundingBox m_bounding_box;
  bool m_has_boundary;
  int m_num_removed_elements_since_last_compact;
  InsertionValidationMode m_insertion_validation_mode;

  ElementFinder m_element_finder;

public:
  /**
   * \brief Default constructor
   * \note User must call initializeBoundary(BoundingBox) before adding points.
   */
  Delaunay()
    : m_has_boundary(false)
    , m_num_removed_elements_since_last_compact(0)
    , m_insertion_validation_mode(InsertionValidationMode::None)
  { }

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
  void initializeBoundary(const BoundingBox& bb)
  {
    std::vector<DataType> points;
    IndexArray elem;

    generateInitialMesh(points, elem, bb);

    m_mesh = IAMeshType(points, elem);
    m_element_finder.recomputeGrid(m_mesh, bb);

    m_bounding_box = bb;
    m_has_boundary = true;
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
  void insertPoint(const PointType& new_pt)
  {
    //Make sure initializeBoundary(...) is called first
    SLIC_ASSERT_MSG(m_has_boundary, "Error: Need a predefined boundary box prior to adding points.");

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

    InsertionHelper insertionHelper(m_mesh);
    insertionHelper.findCavityElements(new_pt, seed_elements);
    validateCavityBoundary(insertionHelper);
    insertionHelper.createCavity();
    IndexType new_pt_i = m_mesh.addVertex(new_pt);
    insertionHelper.delaunayBall(new_pt_i);
    validateInsertedBall(new_pt_i, insertionHelper);
    validateInsertionResult();

    m_element_finder.updateBin(new_pt, new_pt_i);
    m_num_removed_elements_since_last_compact += insertionHelper.numRemovedElements();

    // Compact the mesh if there are too many removed elements
    if(shouldCompactMesh())
    {
      this->compactMesh();
    }
  }

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
      // `Sphere::getOrientation()` depends on an explicitly constructed
      // circumsphere (center + radius). For sliver tetrahedra this construction
      // is ill-conditioned and can produce false positives when validating the
      // empty-circumsphere property. Use the determinant-based predicate used
      // during cavity construction and treat boundary cases as "not inside" for
      // global validation.
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
  /// \brief Returns true when an element and all of its vertices are valid for point-location predicates
  bool isSearchableElement(IndexType element_idx) const
  {
    if(!m_mesh.isValidElement(element_idx))
    {
      return false;
    }

    const auto verts = m_mesh.boundaryVertices(element_idx);
    for(auto idx : verts.positions())
    {
      if(!m_mesh.isValidVertex(verts[idx]))
      {
        return false;
      }
    }

    return true;
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

    // Query mode (`warnOnInvalid == false`) accepts points outside the convex hull,
    // so it uses broader 3D recovery steps before falling back to a full scan.
    const bool use_query_fallbacks = !warnOnInvalid && DIM == 3;
    std::vector<IndexType> candidate_elements = getInitialCandidateElements(query_pt);
    std::vector<IndexType> walked_elements;
    PointLocationResult walk_result =
      walkCandidateElements(query_pt,
                            candidate_elements,
                            0,
                            use_query_fallbacks ? &walked_elements : nullptr);

    if(walk_result.status == PointLocationStatus::Found)
    {
      return walk_result.element_idx;
    }

    if(walk_result.status == PointLocationStatus::Outside)
    {
      return INVALID_INDEX;
    }

    if(use_query_fallbacks)
    {
      walk_result =
        findContainingElementWithQueryFallbacks(query_pt, candidate_elements, walked_elements);
      if(walk_result.status == PointLocationStatus::Found)
      {
        return walk_result.element_idx;
      }
      if(walk_result.status == PointLocationStatus::Outside)
      {
        return INVALID_INDEX;
      }
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

    for(int i = 0; i < VERT_PER_ELEMENT; ++i)
    {
      if(axom::utilities::abs(bary_coord[i]) <= BARY_EPS)
      {
        const IndexType nbr = m_mesh.adjacentElements(element_idx)[ModularFaceIndex(i) + 1];
        if(m_mesh.isValidElement(nbr))
        {
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

  //-----------------------------------------------------------------------------
  // In-sphere predicate helpers
  //
  // Delaunay uses determinant-based in-sphere tests for both cavity growth and
  // global empty-circumsphere validation. We rely on primal::robust::in_sphere()
  // to classify {inside, outside, on boundary} and only compute a
  // scale-dependent determinant tolerance here.
  //-----------------------------------------------------------------------------
  static double inSphereDeterminantTolerance(double scale)
  {
    // Determinant magnitude scales like length^(DIM+2): L^4 in 2D, L^5 in 3D.
    const double k = 128.;
    if constexpr(DIM == 2)
    {
      return k * std::numeric_limits<double>::epsilon() * scale * scale * scale * scale;
    }
    else
    {
      return k * std::numeric_limits<double>::epsilon() * scale * scale * scale * scale * scale;
    }
  }

  static int inSphereOrientationOnMesh(const IAMeshType& mesh, const PointType& q, IndexType element_idx)
  {
    const auto verts = mesh.boundaryVertices(element_idx);

    if constexpr(DIM == 2)
    {
      const PointType& p0 = mesh.getVertexPosition(verts[0]);
      const PointType& p1 = mesh.getVertexPosition(verts[1]);
      const PointType& p2 = mesh.getVertexPosition(verts[2]);

      const auto ba = p1 - p0;
      const auto ca = p2 - p0;
      const auto qa = q - p0;
      const double scale = axom::utilities::max(
        1.,
        axom::utilities::max(ba.norm(), axom::utilities::max(ca.norm(), qa.norm())));
      const double eps = inSphereDeterminantTolerance(scale);

      return primal::robust::in_sphere(q, p0, p1, p2, eps);
    }
    else
    {
      const PointType& p0 = mesh.getVertexPosition(verts[0]);
      const PointType& p1 = mesh.getVertexPosition(verts[1]);
      const PointType& p2 = mesh.getVertexPosition(verts[2]);
      const PointType& p3 = mesh.getVertexPosition(verts[3]);

      const auto ba = p1 - p0;
      const auto ca = p2 - p0;
      const auto da = p3 - p0;
      const auto qa = q - p0;
      const double scale = axom::utilities::max(
        1.,
        axom::utilities::max(
          ba.norm(),
          axom::utilities::max(ca.norm(), axom::utilities::max(da.norm(), qa.norm()))));
      const double eps = inSphereDeterminantTolerance(scale);

      return primal::robust::in_sphere(q, p0, p1, p2, p3, eps);
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
    if(!isSearchableElement(start_element))
    {
      return {};
    }

    static constexpr int MAX_WALK_STEPS = 256;
    std::vector<IndexType> local_visited_elements;
    std::vector<IndexType>& visited_elements =
      visited_elements_out != nullptr ? *visited_elements_out : local_visited_elements;
    visited_elements.clear();
    visited_elements.reserve(MAX_WALK_STEPS);
    IndexType element_i = start_element;

    while(1)
    {
      if(std::find(visited_elements.begin(), visited_elements.end(), element_i) !=
         visited_elements.end())
      {
        return {};
      }
      visited_elements.push_back(element_i);

      const BaryCoordType bary_coord = getBaryCoords(element_i, query_pt);
      ModularFaceIndex modular_idx(0);
      if(isPointInsideForLocation(element_i, query_pt, bary_coord, &modular_idx))
      {
        return {element_i, PointLocationStatus::Found};
      }

      if(static_cast<int>(visited_elements.size()) >= MAX_WALK_STEPS)
      {
        return {};
      }

      const IndexType next_element = m_mesh.adjacentElements(element_i)[modular_idx + 1];
      if(!m_mesh.isValidElement(next_element))
      {
        return {INVALID_INDEX, PointLocationStatus::Outside};
      }

      element_i = next_element;
      if(!isSearchableElement(element_i))
      {
        return {};
      }
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

  std::vector<IndexType> getInitialCandidateElements(const PointType& query_pt) const
  {
    std::vector<IndexType> candidate_elements;
    candidate_elements.reserve(1);

    const IndexType initial_vertex = m_element_finder.getNearbyVertex(query_pt);
    if(m_mesh.isValidVertex(initial_vertex))
    {
      appendCandidateElement(candidate_elements, initial_vertex);
    }

    if(candidate_elements.empty())
    {
      for(auto elem : m_mesh.elements().positions())
      {
        if(isSearchableElement(elem))
        {
          candidate_elements.push_back(elem);
          break;
        }
      }
    }

    return candidate_elements;
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

    const auto fallback_vertices =
      m_element_finder.getNearbyVertices(m_mesh, query_pt, QUERY_SEARCH_RADIUS, QUERY_CANDIDATE_LIMIT);
    const std::size_t initial_candidate_count = candidate_elements.size();
    appendCandidateElementsFromVertices(candidate_elements, fallback_vertices);

    PointLocationResult walk_result =
      walkCandidateElements(query_pt, candidate_elements, initial_candidate_count);
    if(walk_result.status != PointLocationStatus::Failed)
    {
      return walk_result;
    }

    const IndexType nearby_elem = findContainingElementNearby(query_pt, fallback_vertices);
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
    return m_num_removed_elements_since_last_compact > 512 &&
      (m_num_removed_elements_since_last_compact > .2 * m_mesh.elements().size());
  }

  /// \brief Compacts the underlying mesh
  void compactMesh()
  {
    m_mesh.compact();
    m_num_removed_elements_since_last_compact = 0;
    m_element_finder.recomputeGrid(m_mesh, m_bounding_box);
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

      // Use heuristic for resolution in each dimension to minimize storage
      // Use 2*square root of nth root (n==DIM)
      // e.g. for 1,000,000 verts in 2D, sqrt root is 1000, leading to ~ 60^2 grid w/ ~250 verts per bin
      // e.g. for 1,000,000 verts in 3D, cube root is 100, leading to a 20^3 grid w/ ~125 verts per bin
      const double res_root = std::pow(verts.size(), 1.0 / DIM);
      const IndexType res = axom::utilities::max(2, 2 * static_cast<int>(std::sqrt(res_root)));

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

        const auto& pos = mesh.getVertexPosition(idx);
        const auto cell = m_lattice.gridCell(pos);
        flatIndex(cell) = idx;
      }
    }

    /**
     * \brief Returns the index of the vertex in the bin containing point \a pt
     *
     * \param pt The position in space of the vertex that we're checking
     * \note Some bins might not point to a vertex, so users should check
     * that the returned index is a valid vertex, e.g. using \a mesh.isValidVertex(vertex_id)
     */
    inline std::vector<IndexType> getNearbyVertices(const IAMeshType& mesh,
                                                    const PointType& pt,
                                                    int search_radius = 1,
                                                    int max_candidates = 1) const
    {
      const auto cell = m_lattice.gridCell(pt);
      std::vector<std::pair<double, IndexType>> candidates;

      auto tryCandidate = [&](const typename LatticeType::GridCell& candidate_cell) {
        const IndexType vertex_idx = flatIndex(candidate_cell);
        if(mesh.isValidVertex(vertex_idx))
        {
          const double sq_dist = primal::squared_distance(mesh.getVertexPosition(vertex_idx), pt);
          candidates.emplace_back(sq_dist, vertex_idx);
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

      std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
      });

      std::vector<IndexType> nearby_vertices;
      nearby_vertices.reserve(
        axom::utilities::min(max_candidates, static_cast<int>(candidates.size())));
      for(const auto& candidate : candidates)
      {
        nearby_vertices.push_back(candidate.second);
        if(static_cast<int>(nearby_vertices.size()) == max_candidates)
        {
          break;
        }
      }

      return nearby_vertices;
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
  };

  /// Helper struct to locally insert a new point into a Delaunay complex while keeping the mesh Delaunay
  struct InsertionHelper
  {
  public:
    InsertionHelper(IAMeshType& mesh)
      : m_mesh(mesh)
      , facet_set(0)
      , fv_rel(&facet_set, &m_mesh.vertices())
      , fc_rel(&facet_set, &m_mesh.elements())
      , cavity_elems(0)
      , inserted_elems(0)
    { }

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

      IndexArray stack;
      stack.reserve(reserveSize);

      // Seed the cavity with the containing element, and with any face-adjacent
      // neighbors when the insertion point lies on the containing simplex
      // boundary. This avoids repeatedly retriangulating structured inputs that
      // insert directly onto existing edges/faces.
      for(const IndexType element_i : seed_elements)
      {
        if(m_mesh.isValidElement(element_i) && m_checked_element_set.insert(element_i).second)
        {
          cavity_elems.insert(element_i);
          stack.push_back(element_i);
        }
      }

      while(!stack.empty())
      {
        const IndexType element_idx = stack.back();
        stack.pop_back();

        // Invariant: this element is valid, was checked and is in the cavity
        // Each neighbor is either in the cavity or the shared face is on the cavity boundary
        const auto neighbors = m_mesh.adjacentElements(element_idx);
        for(auto n_idx : neighbors.positions())
        {
          const IndexType nbr = neighbors[n_idx];

          // invalid neighbor means face is on domain boundary, and thus on cavity boundary
          if(m_mesh.isValidElement(nbr))
          {
            // neighbor is valid; check circumsphere (if necesary), and add to cavity as appropriate
            if(m_checked_element_set.insert(nbr).second)
            {
              if(isPointInCircumsphere(query_pt, nbr))
              {
                cavity_elems.insert(nbr);
                stack.push_back(nbr);
                continue;  // face is internal to cavity, nothing left to do for this face
              }
            }
            // check if neighbor is already in the cavity
            else if(cavity_elems.findIndex(nbr) != ElementSet::INVALID_ENTRY)
            {
              continue;  // both elem and neighbor along face are in cavity
            }
          }

          // if we got here, the face is on the boundary of the Delaunay cavity
          // add it to facet sets and associated relations
          {
            auto fIdx = facet_set.insert();
            fv_rel.updateSizes();
            fc_rel.updateSizes();

            const auto bdry = m_mesh.boundaryVertices(element_idx);

            auto faceVerts = fv_rel[fIdx];
            typename IAMeshType::ModularVertexIndex mod_idx(n_idx);
            for(int i = 0; i < VERTS_PER_FACET; i++)
            {
              faceVerts[i] = bdry[mod_idx++];
            }
            //For tetrahedron, if the element face is odd, reverse vertex order
            if(DIM == 3 && n_idx % 2 == 1)
            {
              axom::utilities::swap(faceVerts[1], faceVerts[2]);
            }

            fc_rel.insert(fIdx, nbr);
          }
        }
      }

      SLIC_ASSERT_MSG(!cavity_elems.empty(), "Error: New point is not contained in the mesh");
      SLIC_ASSERT(!facet_set.empty());
    }

    /**
    * \brief Remove the elements in the Delaunay cavity
    */
    void createCavity()
    {
      for(auto elem : cavity_elems)
      {
        m_mesh.removeElement(elem);
      }
    }

    /// \brief Fill in the Delaunay cavity with new elements containing the insertion point
    void delaunayBall(IndexType new_pt_i)
    {
      const int numFaces = facet_set.size();
      const IndexType invalid_neighbor = IAMeshType::ElementAdjacencyRelation::INVALID_INDEX;

      IndexType vlist[VERT_PER_ELEMENT] {};
      IndexType neighbors[VERT_PER_ELEMENT] {};
      for(int i = 0; i < numFaces; ++i)
      {
        // Create a new element from the face and the inserted point
        for(int d = 0; d < VERTS_PER_FACET; ++d)
        {
          vlist[d] = fv_rel[i][d];
        }
        vlist[VERTS_PER_FACET] = new_pt_i;

        // Face 0 is the cavity boundary face opposite the inserted point.
        // The remaining faces stay invalid until fixVertexNeighborhood() stitches
        // the new ball together around the inserted vertex.
        const auto nID = fc_rel[i][0];
        for(int d = 0; d < VERT_PER_ELEMENT; ++d)
        {
          neighbors[d] = invalid_neighbor;
        }
        neighbors[0] = nID;

        IndexType new_el = m_mesh.addElement(vlist, neighbors);
        inserted_elems.insert(new_el);
      }

      // Fix neighborhood around the new point
      m_mesh.fixVertexNeighborhood(new_pt_i, inserted_elems.data());
    }

    /// \brief Returns the number of elements removed during this insertion
    int numRemovedElements() const { return cavity_elems.size(); }

    /// \brief Returns true when the query point is inside or on the element circumsphere
    bool isPointInCircumsphere(const PointType& query_pt, IndexType element_idx) const;

  public:
    // we create a surface mesh
    // sets: vertex, facet
    using PositionType = typename IAMeshType::PositionType;
    using ElementType = typename IAMeshType::ElementType;

    using ElementSet = typename IAMeshType::ElementSet;
    using VertexSet = typename IAMeshType::VertexSet;
    using FacetSet = slam::DynamicSet<PositionType, ElementType>;

    // relations: facet->vertex, facet->cell
    static constexpr int VERTS_PER_FACET = IAMeshType::VERTS_PER_ELEM - 1;
    using FacetBoundaryRelation =
      typename IAMeshType::template IADynamicConstantRelation<VERTS_PER_FACET>;
    using FacetCoboundaryRelation = typename IAMeshType::template IADynamicConstantRelation<1>;

  public:
    IAMeshType& m_mesh;

    FacetSet facet_set;
    FacetBoundaryRelation fv_rel;
    FacetCoboundaryRelation fc_rel;

    ElementSet cavity_elems;
    ElementSet inserted_elems;

    std::set<IndexType> m_checked_element_set;
  };
};

template <int DIM>
constexpr typename Delaunay<DIM>::IndexType Delaunay<DIM>::INVALID_INDEX;

template <int DIM>
void Delaunay<DIM>::validateCavityBoundary(const InsertionHelper& insertion_helper) const
{
  if(m_insertion_validation_mode == InsertionValidationMode::None)
  {
    return;
  }

  // Cavity boundary invariant:
  // - Every cavity face that borders a non-cavity neighbor (or the temporary
  //   bounding-box boundary) must appear exactly once in `facet_set`.
  // - The facet's recorded neighbor (`fc_rel`) must match the mesh adjacency.
  struct FacetInfo
  {
    IndexType neighbor_idx {INVALID_INDEX};
    bool matched {false};
  };

  fmt::memory_buffer out;
  bool valid = true;
  std::map<FacetKey, FacetInfo> cavity_boundary;

  for(auto facet_idx : insertion_helper.facet_set.positions())
  {
    const FacetKey facet_key = makeSortedFaceKey(insertion_helper.fv_rel[facet_idx]);
    const auto insert_status =
      cavity_boundary.insert({facet_key, {insertion_helper.fc_rel[facet_idx][0], false}});
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
      if(m_mesh.isValidElement(neighbor_idx) &&
         insertion_helper.cavity_elems.findIndex(neighbor_idx) !=
           InsertionHelper::ElementSet::INVALID_ENTRY)
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

  SLIC_ERROR_IF(!valid, "Delaunay cavity validation failed:" << fmt::to_string(out));
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

  if(insertion_helper.inserted_elems.size() != insertion_helper.facet_set.size())
  {
    fmt::format_to(
      std::back_inserter(out),
      "\n\tInserted ball element count {} does not match cavity boundary facet count {}",
      insertion_helper.inserted_elems.size(),
      insertion_helper.facet_set.size());
    valid = false;
  }

  for(auto facet_idx : insertion_helper.facet_set.positions())
  {
    cavity_boundary.insert({makeSortedFaceKey(insertion_helper.fv_rel[facet_idx]),
                            {insertion_helper.fc_rel[facet_idx][0], false}});
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

  SLIC_ERROR_IF(!valid, "Delaunay ball validation failed:" << fmt::to_string(out));
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
