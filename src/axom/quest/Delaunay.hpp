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
  static constexpr IndexType INVALID_INDEX = -1;

private:
  using ModularFaceIndex =
    slam::ModularInt<slam::policies::CompileTimeSize<IndexType, VERT_PER_ELEMENT>>;

private:
  struct ElementFinder;

  IAMeshType m_mesh;
  BoundingBox m_bounding_box;
  bool m_has_boundary;
  int m_num_removed_elements_since_last_compact;

  ElementFinder m_element_finder;

public:
  /**
   * \brief Default constructor
   * \note User must call initializeBoundary(BoundingBox) before adding points.
   */
  Delaunay() : m_has_boundary(false), m_num_removed_elements_since_last_compact(0) { }

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

    InsertionHelper insertionHelper(m_mesh);
    insertionHelper.findCavityElements(new_pt, seed_elements);
    insertionHelper.createCavity();
    IndexType new_pt_i = m_mesh.addVertex(new_pt);
    insertionHelper.delaunayBall(new_pt_i);

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

      //Collect a list of elements to remove first, because
      //the list may be incomplete if generated during the removal.
      IndexArray elements_to_remove;
      for(int v = 0; v < num_boundary_pts; ++v)
      {
        IndexArray elems = m_mesh.vertexStar(v);
        elements_to_remove.insert(elements_to_remove.end(), elems.begin(), elems.end());
      }
      for(auto e : elements_to_remove)
      {
        if(m_mesh.isValidElement(e))
        {
          m_mesh.removeElement(e);
        }
      }

      for(int v = 0; v < num_boundary_pts; ++v)
      {
        m_mesh.removeVertex(v);
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

    const IndexType totalVertices = m_mesh.vertices().size();
    const IndexType totalElements = m_mesh.elements().size();
    const IndexType res = axom::utilities::ceil(0.33 * std::pow(totalVertices, 1. / DIM));
    UniformGridType grid(m_bounding_box, NumericArray<int, DIM>(res).data());

    // An array to cache the circumspheres associated with each element
    axom::Array<typename ElementType::SphereType> circumspheres(totalElements);

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
        if(circumspheres[element_idx].getOrientation(vertex) == primal::ON_NEGATIVE_SIDE)
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
        for(const auto& pr : invalidEntries)
        {
          const auto vertex_idx = pr.first;
          const auto element_idx = pr.second;
          const auto& pos = m_mesh.getVertexPosition(vertex_idx);
          const auto element = this->getElement(element_idx);
          const auto circumsphere = element.circumsphere();
          fmt::format_to(std::back_inserter(out),
                         "\n\tVertex {} @ {}"
                         "\n\tElement {}: {} w/ circumsphere: {}"
                         "\n\tDistance to circumcenter: {}",
                         vertex_idx,
                         pos,
                         element_idx,
                         element,
                         circumsphere,
                         circumsphere.computeSignedDistance(pos));
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

  static double getPointMagnitudeScale(const std::array<PointType, 4>& pts)
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
    return 64. * std::numeric_limits<double>::epsilon() * scale * scale * scale;
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
    // Simulation-of-simplicity style tie-break: if the 4x4 orientation
    // determinant is effectively zero, use the earliest nonzero cofactor in a
    // fixed symbolic rank order to choose one consistent sign.
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

  int getBarycentricSign(IndexType element_idx,
                         const PointType& query_pt,
                         const BaryCoordType& bary_coord,
                         int bary_idx) const
  {
    const double value = bary_coord[bary_idx];
    if(axom::utilities::abs(value) > BARY_EPS || DIM == 2)
    {
      return signWithTolerance(value, BARY_EPS);
    }

    // The only ambiguous case is a near-zero barycentric coordinate. Interpret
    // that face test symbolically by replacing the corresponding tetrahedron
    // vertex with the query point and evaluating the signed orientation.
    const auto verts = m_mesh.boundaryVertices(element_idx);
    const std::array<PointType, 4> pts {
      {bary_idx == 0 ? query_pt : m_mesh.getVertexPosition(verts[0]),
       bary_idx == 1 ? query_pt : m_mesh.getVertexPosition(verts[1]),
       bary_idx == 2 ? query_pt : m_mesh.getVertexPosition(verts[2]),
       bary_idx == 3 ? query_pt : m_mesh.getVertexPosition(verts[3])}};

    const IndexType query_rank = 1 +
      axom::utilities::max(axom::utilities::max(verts[0], verts[1]),
                           axom::utilities::max(verts[2], verts[3]));
    // Give the query point the highest symbolic rank so zero-case face tests
    // resolve deterministically without changing the stored simplex ordering.
    const std::array<IndexType, 4> ranks {{bary_idx == 0 ? query_rank : verts[0],
                                           bary_idx == 1 ? query_rank : verts[1],
                                           bary_idx == 2 ? query_rank : verts[2],
                                           bary_idx == 3 ? query_rank : verts[3]}};

    const std::array<PointType, 4> tet_pts {{m_mesh.getVertexPosition(verts[0]),
                                             m_mesh.getVertexPosition(verts[1]),
                                             m_mesh.getVertexPosition(verts[2]),
                                             m_mesh.getVertexPosition(verts[3])}};
    const std::array<IndexType, 4> tet_ranks {{verts[0], verts[1], verts[2], verts[3]}};

    const int numerator_sign = symbolicOrientationSign(pts, ranks);
    const int denominator_sign = symbolicOrientationSign(tet_pts, tet_ranks);
    return numerator_sign * denominator_sign;
  }

  bool isPointInsideForLocation(IndexType element_idx,
                                const PointType& query_pt,
                                const BaryCoordType& bary_coord,
                                ModularFaceIndex* exit_face = nullptr) const
  {
    int first_symbolic_negative = -1;
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

      if(axom::utilities::abs(bary_coord[i]) <= BARY_EPS &&
         getBarycentricSign(element_idx, query_pt, bary_coord, i) < 0)
      {
        if(first_symbolic_negative < 0)
        {
          first_symbolic_negative = i;
        }
      }
    }

    if(first_symbolic_negative >= 0)
    {
      if(exit_face != nullptr)
      {
        *exit_face = ModularFaceIndex(first_symbolic_negative);
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

        // set all neighbors to nID; they'll be fixed in the fixVertexNeighborhood function below
        const auto nID = fc_rel[i][0];
        for(int d = 0; d < VERTS_PER_FACET; ++d)
        {
          neighbors[d] = nID;
        }

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
  const auto verts = m_mesh.boundaryVertices(element_idx);

  if constexpr(DIM == 2)
  {
    const PointType& p0 = m_mesh.getVertexPosition(verts[0]);
    const PointType& p1 = m_mesh.getVertexPosition(verts[1]);
    const PointType& p2 = m_mesh.getVertexPosition(verts[2]);
    return primal::in_sphere(query_pt, p0, p1, p2, primal::PRIMAL_TINY, false);
  }
  else
  {
    const PointType& p0 = m_mesh.getVertexPosition(verts[0]);
    const PointType& p1 = m_mesh.getVertexPosition(verts[1]);
    const PointType& p2 = m_mesh.getVertexPosition(verts[2]);
    const PointType& p3 = m_mesh.getVertexPosition(verts[3]);

    const auto ba = p1 - p0;
    const auto ca = p2 - p0;
    const auto da = p3 - p0;
    const auto qa = query_pt - p0;

    const double det = axom::numerics::determinant(ba[0],
                                                   ba[1],
                                                   ba[2],
                                                   ba.squared_norm(),
                                                   ca[0],
                                                   ca[1],
                                                   ca[2],
                                                   ca.squared_norm(),
                                                   da[0],
                                                   da[1],
                                                   da[2],
                                                   da.squared_norm(),
                                                   qa[0],
                                                   qa[1],
                                                   qa[2],
                                                   qa.squared_norm());

    const double scale = axom::utilities::max(
      1.,
      axom::utilities::max(
        ba.norm(),
        axom::utilities::max(ca.norm(), axom::utilities::max(da.norm(), qa.norm()))));
    const double det_tol =
      128. * std::numeric_limits<double>::epsilon() * scale * scale * scale * scale;
    const int det_sign = Delaunay<DIM>::signWithTolerance(det, det_tol);
    if(det_sign != 0)
    {
      return det_sign < 0;
    }

    // Resolve exact co-spherical ties symbolically from the lifted determinant
    // cofactors, ordered by the fixed vertex/query ranks. The query point gets
    // the highest rank so exact ties choose one deterministic inclusive result
    // without perturbing coordinates.
    const IndexType query_rank = 1 +
      axom::utilities::max(axom::utilities::max(verts[0], verts[1]),
                           axom::utilities::max(verts[2], verts[3]));
    const std::array<PointType, 5> lifted_pts {{p0, p1, p2, p3, query_pt}};
    const std::array<IndexType, 5> ranks {{verts[0], verts[1], verts[2], verts[3], query_rank}};
    std::array<int, 5> order {{0, 1, 2, 3, 4}};
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) { return ranks[lhs] < ranks[rhs]; });

    const double cofactor_tol = 128. * std::numeric_limits<double>::epsilon() * scale * scale * scale;
    for(const int row : order)
    {
      std::array<PointType, 4> cofactor_pts;
      int next_idx = 0;
      for(int pt_idx = 0; pt_idx < 5; ++pt_idx)
      {
        if(pt_idx == row)
        {
          continue;
        }
        cofactor_pts[next_idx++] = lifted_pts[pt_idx];
      }

      const double cofactor =
        ((row + 3) % 2 == 0 ? 1. : -1.) * Delaunay<DIM>::orientationDeterminant(cofactor_pts);
      const int sign = Delaunay<DIM>::signWithTolerance(cofactor, cofactor_tol);
      if(sign != 0)
      {
        return sign < 0;
      }
    }

    return true;
  }
}

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
