// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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
#include <cstdlib>
#include <cmath>

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
    InsertionHelper insertionHelper(m_mesh);
    insertionHelper.findCavityElements(new_pt, element_i);
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

    // Find a starting element using ElementFinder helper class
    IndexType element_i = INVALID_INDEX;
    {
      const auto vertex_i = m_element_finder.getNearbyVertex(query_pt);
      if(m_mesh.isValidVertex(vertex_i))
      {
        element_i = m_mesh.coboundaryElement(vertex_i);
      }

      // Fallback -- start from last valid element that was inserted
      if(!m_mesh.isValidElement(element_i))
      {
        element_i = m_mesh.getValidElementIndex();
      }

      SLIC_ASSERT(m_mesh.isValidElement(element_i));
    }

    while(1)
    {
      const BaryCoordType bary_coord = getBaryCoords(element_i, query_pt);

      //Find the index of the most negative barycentric coord
      //Use modular index since it could wrap around to 0
      ModularFaceIndex modular_idx(bary_coord.array().argMin());

      if(bary_coord[modular_idx] >= 0)  // inside if smallest bary coord positive
      {
        return element_i;
      }

      // else, move to that neighbor
      element_i = m_mesh.adjacentElements(element_i)[modular_idx + 1];

      // Either there is a hole in the m_mesh, or the point is outside of the m_mesh.
      // Logically, this should never happen.
      if(!m_mesh.isValidElement(element_i))
      {
        SLIC_WARNING_IF(warnOnInvalid,
                        fmt::format("Entered invalid element in "
                                    "Delaunay::findContainingElement(). Underlying mesh {} valid",
                                    m_mesh.isValid() ? "is" : "is not"));
        return INVALID_INDEX;
      }
    }
  }

  /**
   * \brief helper function to retrieve the barycentric coordinate of the query point in the element
   */
  BaryCoordType getBaryCoords(IndexType element_idx, const PointType& q_pt) const;

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
   * neighboring elements for a list of element indices whose circumspheres contains the query point.
   * It also finds the faces on the boundaries of the cavity to help with filling the cavity
   * in the \a delaunayBall function.  
   *
   * \param query_pt the query point
   * \param element_i the element to start the search at
   */
    void findCavityElements(const PointType& query_pt, IndexType element_i)
    {
      constexpr int reserveSize = (DIM == 2) ? 16 : 64;

      IndexArray stack;
      stack.reserve(reserveSize);

      // add first element (if valid and point is in its circumsphere)
      if(m_mesh.isValidElement(element_i) && isPointInSphere(query_pt, element_i))
      {
        m_checked_element_set.insert(element_i);
        cavity_elems.insert(element_i);
        stack.push_back(element_i);
      }

      while(!stack.empty())
      {
        IndexType element_idx = stack.back();
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
              if(isPointInSphere(query_pt, nbr))
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

    /// \brief Helper function returns true if the query point is in the sphere formed by the element vertices
    bool isPointInSphere(const PointType& query_pt, IndexType element_idx) const;

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

// 2D specialization for getBaryCoords(...)
template <>
inline Delaunay<2>::BaryCoordType Delaunay<2>::getBaryCoords(IndexType element_idx,
                                                             const PointType& query_pt) const
{
  const auto verts = m_mesh.boundaryVertices(element_idx);
  const ElementType tri(m_mesh.getVertexPosition(verts[0]),
                        m_mesh.getVertexPosition(verts[1]),
                        m_mesh.getVertexPosition(verts[2]));

  return tri.physToBarycentric(query_pt);
}

// 3D specialization for getBaryCoords(...)
template <>
inline Delaunay<3>::BaryCoordType Delaunay<3>::getBaryCoords(IndexType element_idx,
                                                             const PointType& query_pt) const
{
  const auto verts = m_mesh.boundaryVertices(element_idx);
  const ElementType tet(m_mesh.getVertexPosition(verts[0]),
                        m_mesh.getVertexPosition(verts[1]),
                        m_mesh.getVertexPosition(verts[2]),
                        m_mesh.getVertexPosition(verts[3]));

  return tet.physToBarycentric(query_pt);
}

// 2D specialization for isPointInSphere(...)
template <>
inline bool Delaunay<2>::InsertionHelper::isPointInSphere(const PointType& query_pt,
                                                          IndexType element_idx) const
{
  const auto verts = m_mesh.boundaryVertices(element_idx);
  const PointType& p0 = m_mesh.getVertexPosition(verts[0]);
  const PointType& p1 = m_mesh.getVertexPosition(verts[1]);
  const PointType& p2 = m_mesh.getVertexPosition(verts[2]);
  return primal::in_sphere(query_pt, p0, p1, p2, 0.);
}

// 3D specialization for isPointInSphere(...)
template <>
inline bool Delaunay<3>::InsertionHelper::isPointInSphere(const PointType& query_pt,
                                                          IndexType element_idx) const
{
  const auto verts = m_mesh.boundaryVertices(element_idx);
  const PointType& p0 = m_mesh.getVertexPosition(verts[0]);
  const PointType& p1 = m_mesh.getVertexPosition(verts[1]);
  const PointType& p2 = m_mesh.getVertexPosition(verts[2]);
  const PointType& p3 = m_mesh.getVertexPosition(verts[3]);
  return primal::in_sphere(query_pt, p0, p1, p2, p3, 0.);
}

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
