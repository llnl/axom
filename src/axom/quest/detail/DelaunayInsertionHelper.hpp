// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file DelaunayInsertionHelper.hpp
 *
 * \brief Defines the cavity-construction and retriangulation helper used by
 * incremental Delaunay insertion.
 *
 * Implements the Bowyer-Watson algorithm for incremental point insertion:
 * 1. findCavityElements(): Expands cavity from seed elements using circumsphere test
 * 2. createCavity(): Removes cavity elements from the mesh
 * 3. delaunayBall(): Retriangulates by connecting new point to cavity boundary
 *
 * The helper is reused across insertions to avoid repeated allocations.
 */

#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/slam.hpp"

#include <array>
#include <vector>

namespace axom
{
namespace quest
{
namespace detail
{

/**
 * \brief Tracks the local cavity state for a single Bowyer-Watson insertion.
 *
 * The owning `quest::Delaunay` instance reuses one helper across insertions so
 * cavity membership, boundary facets, and inserted-element scratch storage can
 * be cleared without reallocation on every point.
 *
 * Public members track insertion state for validation:
 * - cavity_elems: Elements whose circumspheres contain the new point
 * - boundary_facets: Faces between cavity and non-cavity elements
 * - inserted_elems: New simplices created by connecting new point to boundary
 * - containing_element, containing_bary, seed_elements_debug: For diagnostics
 *
 * \tparam DIM The spatial dimension (2 for triangles, 3 for tetrahedra)
 * \tparam PointType The point type used for the new (query) point, e.g. primal::Point
 * \tparam BaryCoordType The barycentric-coordinate type used during point location
 * \tparam IndexType The integer type indexing vertices and elements of the mesh
 * \tparam IndexArray The container type used to pass a set of seed elements
 * \tparam IAMeshType The topological mesh data structure type (slam::IAMesh) being modified
 */
template <int DIM, typename PointType, typename BaryCoordType, typename IndexType, typename IndexArray, typename IAMeshType>
class DelaunayInsertionHelper
{
public:
  static constexpr int VERT_PER_ELEMENT = DIM + 1;
  static constexpr int VERTS_PER_FACET = VERT_PER_ELEMENT - 1;
  static constexpr IndexType INVALID_INDEX = IndexType {-1};

  /**
   * \brief A facet on the boundary of the cavity, paired with its adjacent non-cavity element
   *
   *  Used to re-stitch the triangulation after the cavity is removed.
   */
  struct BoundaryFacet
  {
    std::array<IndexType, VERTS_PER_FACET> vertices {};
    IndexType neighbor {INVALID_INDEX};
  };

  /**
   * \brief Construct the helper bound to the mesh it will modify in place.
   *
   * \param mesh The IA mesh whose cavity will be carved and re-triangulated
   */
  explicit DelaunayInsertionHelper(IAMeshType& mesh) : m_mesh(mesh) { }

  /**
   * \brief Clears all per-insertion scratch state so the helper can be reused
   *  for the next point without reallocating its internal buffers.
   */
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
   * \brief Computes the Bowyer-Watson cavity: the set of elements whose circumspheres
   *  contain the query point, found by a flood fill outward from the seed elements across shared facets.
   *
   * Populates \a cavity_elems (the elements to delete) and their \a boundary_facets
   * 
   * \tparam CircumspherePredicate Callable `(IndexType elem) -> bool` returning
   *   true when the query point lies inside element \a elem's circumsphere
   *
   * \param query_pt The point being inserted
   * \param seed_elements One or more elements known to contain \a query_pt in
   *   their circumsphere, used to seed the flood fill
   * \param isPointInCircumsphere The in-circumsphere predicate (see tparam)
   */
  template <typename CircumspherePredicate>
  void findCavityElements(const PointType& query_pt,
                          const IndexArray& seed_elements,
                          CircumspherePredicate&& isPointInCircumsphere)
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

    for(const IndexType element_i : seed_elements)
    {
      SLIC_ASSERT(m_mesh.isValidElement(element_i));
      addCavityElement(element_i);
    }

    while(!m_stack.empty())
    {
      const IndexType element_idx = m_stack.back();
      m_stack.pop_back();

      const auto neighbors = m_mesh.adjacentElements(element_idx);
      for(auto n_idx : neighbors.positions())
      {
        const IndexType nbr = neighbors[n_idx];

        if(nbr != invalid_element)
        {
          SLIC_ASSERT(m_mesh.isValidElement(nbr));

          if(containsCavityElement(nbr))
          {
            continue;
          }

          if(isPointInCircumsphere(query_pt, nbr))
          {
            addCavityElement(nbr);
            continue;
          }
        }

        const auto bdry = m_mesh.boundaryVertices(element_idx);

        typename IAMeshType::ModularVertexIndex mod_idx(static_cast<int>(n_idx));
        BoundaryFacet facet;
        for(int i = 0; i < VERTS_PER_FACET; i++)
        {
          facet.vertices[i] = bdry[mod_idx++];
        }
        if(DIM == 3 && n_idx % 2 == 1)
        {
          axom::utilities::swap(facet.vertices[1], facet.vertices[2]);
        }

        facet.neighbor = nbr;
        boundary_facets.push_back(facet);
      }
    }

    SLIC_ASSERT_MSG(!cavity_elems.empty(), "Error: New point is not contained in the mesh");
    SLIC_ASSERT(!boundary_facets.empty());
  }

  /**
   * \brief Removes the cavity elements from the mesh, releasing their slots to
   *  the deleted-element pool for reuse by the new simplices.
   *
   * \tparam DeletedElementPool A pool type exposing `release(IndexType)`
   * \param deleted_elements The pool that receives the freed element slots
   */
  template <typename DeletedElementPool>
  void createCavity(DeletedElementPool& deleted_elements)
  {
    for(const auto elem : cavity_elems)
    {
      m_mesh.removeElement(elem);
      deleted_elements.release(elem);
    }
  }

  /**
   * \brief Retriangulates the cavity by connecting the new point to each
   *  boundary facet (forming the "Delaunay ball"), reusing freed slots first.
   *
   * \tparam DeletedElementPool A pool type exposing reusable element slots
   * \param new_pt_i The vertex index of the just-inserted point
   * \param deleted_elements The pool of element slots freed by createCavity()
   */
  template <typename DeletedElementPool>
  void delaunayBall(IndexType new_pt_i, DeletedElementPool& deleted_elements)
  {
    const int numFaces = static_cast<int>(boundary_facets.size());
    const IndexType invalid_neighbor = IAMeshType::ElementAdjacencyRelation::INVALID_INDEX;
    const PointType& new_pt = m_mesh.getVertexPosition(new_pt_i);

    IndexType vlist[VERT_PER_ELEMENT] {};
    IndexType neighbors[VERT_PER_ELEMENT] {};
    for(int i = 0; i < numFaces; ++i)
    {
      for(int d = 0; d < VERTS_PER_FACET; ++d)
      {
        vlist[d] = boundary_facets[static_cast<std::size_t>(i)].vertices[d];
      }
      vlist[VERTS_PER_FACET] = new_pt_i;

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

    m_mesh.fixVertexNeighborhood(new_pt_i, inserted_elems);
  }

  int numRemovedElements() const { return static_cast<int>(cavity_elems.size()); }

  /**
   * \brief Returns true if \a element_idx is a member of the current cavity.
   *
   * \param element_idx The element to test for cavity membership
   */
  bool containsCavityElement(IndexType element_idx) const
  {
    return element_idx >= 0 && static_cast<std::size_t>(element_idx) < m_cavity_membership.size() &&
      m_cavity_membership[static_cast<std::size_t>(element_idx)] != 0;
  }

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

}  // namespace detail
}  // namespace quest
}  // namespace axom
