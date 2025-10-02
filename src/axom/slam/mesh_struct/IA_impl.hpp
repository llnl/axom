// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_IA_IMPL_H_
#define SLAM_IA_IMPL_H_

/*
 * \file IA_impl.hpp
 *
 * \brief Contains the implementation of the IAMesh class and helper functions
 */

#include "axom/core/Macros.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/ModularInt.hpp"

#include "axom/fmt.hpp"

#include <vector>
#include <map>

namespace axom
{
namespace slam
{
/**
 * Helper function on std::vectors
 */

namespace /*anonymous*/
{
// Checks if \a v is in the list \a s
template <typename T, typename IterableT>
bool is_subset(T v, const IterableT& iterable)
{
  for(auto item : iterable)
  {
    if(item == v)
    {
      return true;
    }
  }
  return false;
}

// Formatted output of a relation or map to an array of strings
// helper function for IAMesh::print_all()
template <typename RelOrMap, typename SetType>
std::vector<std::string> entries_as_vec(const RelOrMap& outer, const SetType& s)
{
  const int sz = outer.size();
  std::vector<std::string> strs(sz);
  for(auto pos : s.positions())
  {
    strs[pos] = s.isValidEntry(pos) ? fmt::format("{}: {}", axom::fmt::streamed(pos), outer[pos])
                                    : fmt::format("{}: {{}}", axom::fmt::streamed(pos));
  }

  return strs;
}

}  //end anonymous namespace

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::ElementAndFaceIdxType IAMesh<TDIM, SDIM, P>::ElemNbrFinder(
  V2EMapType& vertpair_to_elem_map,
  IndexType element_i,
  IndexType side_i)
{
  // NOTE: V2EMapType maps a sorted tuple of vertex IDs to a face on a given
  //       mesh element. It is used to find the element index of the opposite
  //       face within the mesh

  IndexArray vlist = getElementFace(element_i, side_i);
  std::sort(vlist.begin(), vlist.end());

  ElementAndFaceIdxType zs_pair(element_i, side_i);

  auto map_ret2 = vertpair_to_elem_map.insert(std::make_pair(vlist, zs_pair));
  if(!map_ret2.second)  //if this pair is in the map, we've found our match
  {
    auto orig_pair = map_ret2.first->second;
    vertpair_to_elem_map.erase(map_ret2.first);
    return orig_pair;
  }

  //No matching pair is found. Return an invalid pair
  return ElementAndFaceIdxType(-1, -1);
}

template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::print_all() const
{
  fmt::memory_buffer out;
  fmt::format_to(std::back_inserter(out),
                 "IA mesh: {} mesh in {}d with {} valid vertices (of {}) "
                 "and {} valid elements (of {})\n",
                 (TDIM == 2 ? "triangle" : "tetrahedral"),
                 SDIM,
                 vertex_set.numberOfValidEntries(),
                 vertex_set.size(),
                 element_set.numberOfValidEntries(),
                 element_set.size());

  //print out the element and vertex sets
  fmt::format_to(std::back_inserter(out),
                 "  element_set ({}/{}): [{}]\n",
                 element_set.numberOfValidEntries(),
                 element_set.size(),
                 fmt::join(element_set, ", "));

  fmt::format_to(std::back_inserter(out),
                 "  vertex_set ({}/{}): [{}]\n",
                 vertex_set.numberOfValidEntries(),
                 vertex_set.size(),
                 fmt::join(vertex_set, ", "));

  //print out the relations on the sets (ev, ve and ee)
  fmt::format_to(std::back_inserter(out),
                 "  ev_rel ({}/{}): [{}]\n",
                 ev_rel.numberOfValidEntries(),
                 ev_rel.size(),
                 fmt::join(entries_as_vec(ev_rel, element_set), "; "));

  fmt::format_to(std::back_inserter(out),
                 "  ve_rel ({}/{}): [{}]\n",
                 ve_rel.numberOfValidEntries(),
                 ve_rel.size(),
                 fmt::join(entries_as_vec(ve_rel, vertex_set), "; "));

  fmt::format_to(std::back_inserter(out),
                 "  ee_rel ({}/{}): [{}]\n",
                 ee_rel.numberOfValidEntries(),
                 ee_rel.size(),
                 fmt::join(entries_as_vec(ee_rel, element_set), "; "));

  //print out the coordinate map (i.e the positions)
  fmt::format_to(std::back_inserter(out),
                 "  vertex coord ({}/{}): [{}]\n",
                 vcoord_map.numberOfValidEntries(),
                 vcoord_map.size(),
                 fmt::join(entries_as_vec(vcoord_map, vertex_set), "; "));

  SLIC_INFO(fmt::to_string(out));
}

//-----------------------------------------------------------------------------

template <int TDIM, int SDIM, typename P>
IAMesh<TDIM, SDIM, P>::IAMesh()
  : vertex_set(0)
  , element_set(0)
  , ev_rel(&element_set, &vertex_set)
  , ve_rel(&vertex_set, &element_set)
  , ee_rel(&element_set, &element_set)
  , vcoord_map(&vertex_set)
{ }

template <int TDIM, int SDIM, typename P>
IAMesh<TDIM, SDIM, P>::IAMesh(const IAMesh& m)
{
  operator=(m);
}

template <int TDIM, int SDIM, typename P>
IAMesh<TDIM, SDIM, P>& IAMesh<TDIM, SDIM, P>::operator=(const IAMesh& m)
{
  if(&m != this)
  {
    vertex_set = m.vertex_set;
    element_set = m.element_set;
    ev_rel = ElementBoundaryRelation(&element_set, &vertex_set);
    ve_rel = VertexCoboundaryRelation(&vertex_set, &element_set);
    ee_rel = ElementAdjacencyRelation(&element_set, &element_set);
    vcoord_map = PositionMap(&vertex_set);

    ev_rel.data() = m.ev_rel.data();
    ve_rel.data() = m.ve_rel.data();
    ee_rel.data() = m.ee_rel.data();
    vcoord_map.data() = m.vcoord_map.data();
  }

  return *this;
}

template <int TDIM, int SDIM, typename P>
IAMesh<TDIM, SDIM, P>::IAMesh(std::vector<double>& points, std::vector<IndexType>& tri)
  : vertex_set(points.size() / COORDS_PER_VERT)
  , element_set(tri.size() / VERTS_PER_ELEM)
  , ev_rel(&element_set, &vertex_set)
  , ve_rel(&vertex_set, &element_set)
  , ee_rel(&element_set, &element_set)
  , vcoord_map(&vertex_set)
{
  // Relation, element to vertex boundary relation
  for(auto e : element_set)
  {
    const int offset = VERTS_PER_ELEM * e;
    for(int idx = 0; idx < VERTS_PER_ELEM; ++idx)
    {
      ev_rel.insert(e, tri[offset + idx]);
    }
  }
  SLIC_ASSERT_MSG(ev_rel.isValid(), "Error creating (dynamic) relation from elements to vertices!");

  // The map, vertex to coordinates
  for(auto v : vertex_set)
  {
    vcoord_map[v] = Point(&(points[v * COORDS_PER_VERT]));
  }
  SLIC_ASSERT_MSG(vcoord_map.isValid(true), "Error creating map from vertex to coords!");

  //Vertex element relation. 1->1 mapping only 1 element per vertex.
  for(auto e : element_set)
  {
    for(auto v : ev_rel[e])
    {
      ve_rel.modify(v, 0, e);
    }
  }
  SLIC_ASSERT_MSG(ve_rel.isValid(true),
                  "Error creating (dynamic) relation from vertices to elements!\n");

  //Before making element to element relation, construct the data.
  // For every cell, find the union of triangles for each pair of vertices
  IndexArray element_element_vec(element_set.size() * VERTS_PER_ELEM,
                                 ElementBoundaryRelation::INVALID_INDEX);

  V2EMapType vertpair_to_elem_map;

  for(auto element_i : element_set)
  {
    for(IndexType side_i = 0; side_i < VERTS_PER_ELEM; ++side_i)
    {
      ElementAndFaceIdxType nst = ElemNbrFinder(vertpair_to_elem_map, element_i, side_i);

      IndexType other_element_idx = nst.first;
      IndexType other_side_idx = nst.second;

      if(element_set.isValidEntry(other_element_idx))
      {
        int idx0 = element_i * VERTS_PER_ELEM + side_i;
        element_element_vec[idx0] = other_element_idx;

        int idx1 = other_element_idx * VERTS_PER_ELEM + other_side_idx;
        element_element_vec[idx1] = element_i;
      }
    }
  }

  //Element adjacency relation along facets
  for(auto i : element_set)
  {
    for(int j = 0; j < VERTS_PER_ELEM; j++)
    {
      ee_rel.modify(i, j, element_element_vec[i * VERTS_PER_ELEM + j]);
    }
  }
  SLIC_ASSERT_MSG(ee_rel.isValid(true),
                  "Error creating (dynamic) relation from elements to elements!");
}

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexArray IAMesh<TDIM, SDIM, P>::vertexStar(IndexType vertex_idx) const
{
  // reasonable expected size of vertex star in triangle and tet meshes
  constexpr int EXP_SZ = (TDIM == 2) ? 8 : 32;

  IndexArray ret;
  ret.reserve(EXP_SZ);

  if(!ve_rel.isValidEntry(vertex_idx))
  {
    //this vertex is not connected to any elements
    SLIC_WARNING_IF(!vertex_set.isValidEntry(vertex_idx),
                    "Attempting to retrieve data with an invalid vertex id: " << vertex_idx);
    return ret;
  }

  IndexType starting_element_idx = ve_rel[vertex_idx][0];

  ret.push_back(starting_element_idx);
  IndexArray element_traverse_queue;
  element_traverse_queue.push_back(starting_element_idx);

  while(!element_traverse_queue.empty())
  {
    IndexType element_i = element_traverse_queue.back();
    element_traverse_queue.pop_back();

    for(auto nbr : ee_rel[element_i])
    {
      // If nbr is valid, has not already been found and contains the vertex in question,
      // add it and enqueue to check neighbors.
      if(element_set.isValidEntry(nbr)  //
         && !is_subset(nbr, ret)        //
         && is_subset(vertex_idx, ev_rel[nbr]))
      {
        ret.push_back(nbr);
        element_traverse_queue.push_back(nbr);
      }
    }
  }

  return ret;
}

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexArray IAMesh<TDIM, SDIM, P>::getElementFace(IndexType element_idx,
                                                                                 IndexType face_idx) const
{
  constexpr int VERTS_PER_FACET = VERTS_PER_ELEM - 1;

  ModularVertexIndex mod_face(face_idx);

  IndexArray ret;
  ret.reserve(VERTS_PER_FACET);

  if(!element_set.isValidEntry(element_idx))
  {
    SLIC_WARNING("Attempting to retrieve data with an invalid element: " << element_idx);

    return ret;
  }

  SLIC_ASSERT_MSG(0 <= face_idx && face_idx < VERTS_PER_ELEM, "Face index is invalid.");

  auto ev = ev_rel[element_idx];
  for(int i = 0; i < VERTS_PER_FACET; ++i)
  {
    ret.push_back(ev[mod_face + i]);
  }

  return ret;
}

template <int TDIM, int SDIM, typename P>
const typename IAMesh<TDIM, SDIM, P>::Point& IAMesh<TDIM, SDIM, P>::getVertexPosition(
  IndexType vertex_idx) const
{
  SLIC_ASSERT(isValidVertex(vertex_idx));

  return vcoord_map[vertex_idx];
}

template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::removeVertex(IndexType vertex_idx)
{
  if(!vertex_set.isValidEntry(vertex_idx))
  {
    SLIC_WARNING("Attempting to remove an invalid vertex");
    return;
  }

  //check if any element uses this vertex. If so, remove them too.
  for(auto incident_element : vertexStar(vertex_idx))
  {
    removeElement(incident_element);
  }

  ve_rel.remove(vertex_idx);
  vertex_set.remove(vertex_idx);
  //Note: after set entry is removed, its corresponding map entry will be invalid
}

template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::removeElement(IndexType element_idx)
{
  if(!element_set.isValidEntry(element_idx))
  {
    SLIC_WARNING("Attempting to remove an invalid element");
    return;
  }

  // update vertex coboundary relation for vertices of the removed cell (when necessary)
  for(auto vertex_i : ev_rel[element_idx])
  {
    // update VE relation for vertex_i when it points to deleted element
    IndexType& cbdry = coboundaryElement(vertex_i);
    if(cbdry == element_idx)
    {
      IndexType new_elem = ElementSet::INVALID_ENTRY;
      for(auto nbr : ee_rel[element_idx])
      {
        // update to a valid neighbor that is incident in vertex_i
        if(element_set.isValidEntry(nbr) && is_subset(vertex_i, ev_rel[nbr]))
        {
          new_elem = nbr;
          break;
        }
      }
      cbdry = new_elem;
    }
  }

  //erase this element and it boundary relation
  ev_rel.remove(element_idx);

  //erase neighbor element's adjacency data pointing to deleted element
  for(auto nbr : ee_rel[element_idx])
  {
    if(isValidElement(nbr))
    {
      auto nbr_ee = ee_rel[nbr];
      for(auto idx : nbr_ee.positions())
      {
        if(nbr_ee[idx] == element_idx)
        {
          nbr_ee[idx] = ElementBoundaryRelation::INVALID_INDEX;
          break;
        }
      }
    }
  }
  ee_rel.remove(element_idx);

  element_set.remove(element_idx);
}

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::addVertex(const Point& p)
{
  IndexType vertex_idx = vertex_set.insert();
  vcoord_map.insert(vertex_idx, p);
  ve_rel.insert(vertex_idx, VertexCoboundaryRelation::INVALID_INDEX);

  return vertex_idx;
}

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::addElement(IndexType v0,
                                                                            IndexType v1,
                                                                            IndexType v2,
                                                                            IndexType v3)
{
  SLIC_ASSERT(VERTS_PER_ELEM <= 4);
  IndexType vlist[] = {v0, v1, v2, v3};
  return addElement(vlist);
}

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::addElement(const IndexType* vlist)
{
  // Implementation note:
  //   This function reconstructs the vertex-element relation
  //    on each vertex ID of the new element
  // Can we optimize this function?

  for(int i = 0; i < VERTS_PER_ELEM; i++)
  {
    SLIC_WARNING_IF(!vertex_set.isValidEntry(vlist[i]),
                    "Trying to add an element with invalid vertex index:" << vlist[i]);
  }

  IndexType element_idx = element_set.insert();
  ev_rel.updateSizes();
  ee_rel.updateSizes();

  auto bdry = ev_rel[element_idx];
  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    bdry[i] = vlist[i];
  }

  //make sure the space is allocated in ee_rel
  ee_rel.insert(element_idx, ElementAdjacencyRelation::INVALID_INDEX);

  V2EMapType vertpair_to_elem_map;

  //First add each face of this new element into the map
  for(int side_i = 0; side_i < VERTS_PER_ELEM; ++side_i)
  {
    ElementAndFaceIdxType zs_pair = ElemNbrFinder(vertpair_to_elem_map, element_idx, side_i);
    SLIC_ASSERT(zs_pair.first == -1);
    AXOM_UNUSED_VAR(zs_pair);
  }

  //Make a list of elements that shares at least 1 vertex of the new element
  std::set<IndexType> elem_list;
  for(int n = 0; n < VERTS_PER_ELEM; ++n)
  {
    for(auto elem : vertexStar(vlist[n]))
    {
      elem_list.insert(elem);
    }
  }

  //Check if any of the elements share a face with the new element.
  // If so, modify ee_rel to reflect that.
  for(auto otherElementIdx : elem_list)
  {
    if(otherElementIdx < 0 || otherElementIdx == element_idx)
    {
      continue;
    }
    for(IndexType s = 0; s < VERTS_PER_ELEM; s++)
    {
      IndexType otherSideIdx = s;
      // insert the pair
      ElementAndFaceIdxType zs_pair =
        ElemNbrFinder(vertpair_to_elem_map, otherElementIdx, otherSideIdx);

      //If zs_pair returned is the new element, save this nbr to set later
      IndexType foundElementIdx = zs_pair.first;
      IndexType foundSideIdx = zs_pair.second;

      if(foundElementIdx == element_idx)
      {
        // if there is already a neighbor on the save list, this mesh is not a
        // manifold.  Example: Having an edge with 3 faces...

        SLIC_ASSERT(ee_rel[otherElementIdx][otherSideIdx] < 0);

        ee_rel.modify(foundElementIdx, foundSideIdx, otherElementIdx);
        ee_rel.modify(otherElementIdx, otherSideIdx, foundElementIdx);

        //put new element pair back in queue to check if mesh is manifold
        ElemNbrFinder(vertpair_to_elem_map, foundElementIdx, foundSideIdx);
      }
    }
  }

  //update ve_rel
  for(int i = 0; i < VERTS_PER_ELEM; i++)
  {
    if(!ve_rel.isValidEntry(vlist[i]))
    {
      ve_rel.modify(vlist[i], 0, element_idx);
    }
  }

  return element_idx;
}

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::addElement(const IndexType* vlist,
                                                                            const IndexType* neighbors)
{
  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    SLIC_ASSERT_MSG(vertex_set.isValidEntry(vlist[i]),
                    "Trying to add an element with invalid vertex index:" << vlist[i]);
  }

  IndexType element_idx = element_set.insert();
  ev_rel.updateSizes();
  ee_rel.updateSizes();

  // set the vertices in this element's ev relation
  auto bdry = ev_rel[element_idx];
  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    bdry[i] = vlist[i];
  }

  // set the neighbor elements in this element's ee relation
  auto adj = ee_rel[element_idx];
  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    adj[i] = neighbors[i];
  }

  // update coboundary relation of this element's vertices, if necessary
  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    const IndexType v = vlist[i];
    IndexType& cbdry = coboundaryElement(v);

    if(!element_set.isValidEntry(cbdry))
    {
      cbdry = element_idx;
    }
  }

  return element_idx;
}

template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::fixVertexNeighborhood(IndexType vertex_idx,
                                                  const std::vector<IndexType>& new_elements)
{
  using FaceLinkVerts = axom::StackArray<IndexType, TDIM - 1>;

  // Struct used to associate a collection of vertex indices with a face of the mesh
  struct FaceLinkMapping
  {
    IndexType boundary_hash;  // unique identifier for collection of face link verts
    IndexType element_idx;    // the element containing this face
    IndexType face_idx;       // the index of the face w.r.t. the element
  };

  const IndexType totalVerts = elements().size();
  std::vector<FaceLinkMapping> mapping;
  mapping.reserve(TDIM * new_elements.size());

  // helper lambda for determining if a face on one element (given by boundary verts nbr_verts)
  // is shared with another element (given by boundary verts elem_verts)
  auto isSharedFace =
    [](const BoundarySubset& nbr_verts, IndexType face_idx, const BoundarySubset& elem_verts) {
      // v_skip is not a vertex in this face
      const auto v_skip = (face_idx == 0) ? VERTS_PER_ELEM - 1 : face_idx - 1;
      for(int v_idx = 0; v_idx < VERTS_PER_ELEM; ++v_idx)
      {
        if(v_idx != v_skip &&                         // v_idx is a vertex in this face
           !is_subset(nbr_verts[v_idx], elem_verts))  //and is not in other elem
        {
          return false;
        }
      }
      return true;
    };

  for(auto el : new_elements)
  {
    const auto bdry = ev_rel[el];
    for(int face_i = 0; face_i < VERTS_PER_ELEM; ++face_i)
    {
      // This face is either a boundary facet of the star...
      const auto vert_i = (face_i == 0) ? VERTS_PER_ELEM - 1 : face_i - 1;
      if(bdry[vert_i] == vertex_idx)
      {
        // figure out which face this is on the neighbor
        // and update neighbor's adjacency to point to current element
        const IndexType nbr = ee_rel[el][face_i];
        if(element_set.isValidEntry(nbr))
        {
          const auto nbr_ev = ev_rel[nbr];
          auto nbr_ee = ee_rel[nbr];

          for(auto face_j : nbr_ev.positions())
          {
            if(!element_set.isValidEntry(nbr_ee[face_j]) && isSharedFace(nbr_ev, face_j, bdry))
            {
              nbr_ee[face_j] = el;
              break;
            }
          }
        }
      }
      // ... or it is incident in the common vertex: vertex_idx
      else
      {
        // Add all element-face associations to an array;
        // we'll update the adjacencies outside this loop
        FaceLinkVerts face_link_verts {};
        for(int i = 0, idx = 0; i < TDIM; ++i)
        {
          if(i != vert_i &&          // i is a vertex in face_i
             bdry[i] != vertex_idx)  // and not the common vertex
          {
            face_link_verts[idx++] = bdry[i];
          }
        }

        FaceLinkMapping m;
        if(TDIM == 2)
        {
          m.boundary_hash = face_link_verts[0];
        }
        else  // TDIM == 3
        {
          // compute unique identifier based on sorted face link vertices
          m.boundary_hash = (face_link_verts[0] < face_link_verts[1])
            ? face_link_verts[0] * totalVerts + face_link_verts[1]
            : face_link_verts[1] * totalVerts + face_link_verts[0];
        }

        m.element_idx = el;
        m.face_idx = face_i;
        mapping.emplace_back(m);
      }
    }
  }

  SLIC_ASSERT(mapping.size() == new_elements.size() * TDIM);

  // Sort by face vertices
  std::sort(mapping.begin(), mapping.end(), [](const FaceLinkMapping& lhs, const FaceLinkMapping& rhs) {
    return lhs.boundary_hash < rhs.boundary_hash;
  });

  // Apply neighbor data from matching face pairs
  const int SZ = mapping.size();
  for(int idx = 1; idx < SZ; idx += 2)
  {
    const auto& left = mapping[idx - 1];
    const auto& right = mapping[idx];
    ee_rel.modify(left.element_idx, left.face_idx, right.element_idx);
    ee_rel.modify(right.element_idx, right.face_idx, left.element_idx);
  }
}

// Remove all the invalid entries in the IA structure
template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::compact()
{
  constexpr IndexType INVALID_VERTEX = VertexSet::INVALID_ENTRY;
  constexpr IndexType INVALID_ELEMENT = ElementSet::INVALID_ENTRY;

  //Construct an array that maps original set indices to new compacted indices
  IndexArray vertex_set_map(vertex_set.size(), INVALID_VERTEX);
  IndexArray element_set_map(element_set.size(), INVALID_ELEMENT);

  int v_count = 0;
  for(auto v : vertex_set.positions())
  {
    if(vertex_set.isValidEntry(v))
    {
      vertex_set_map[v] = v_count++;
    }
  }

  int e_count = 0;
  for(auto e : element_set.positions())
  {
    if(element_set.isValidEntry(e))
    {
      element_set_map[e] = e_count++;
    }
  }

  //update the EV boundary relation
  for(auto e : element_set.positions())
  {
    const auto new_e = element_set_map[e];
    if(new_e != INVALID_ELEMENT)
    {
      const auto ev_old = ev_rel[e];
      auto ev_new = ev_rel[new_e];
      for(auto i : ev_new.positions())
      {
        const auto old = ev_old[i];
        ev_new[i] = (old != INVALID_VERTEX) ? vertex_set_map[old] : INVALID_VERTEX;
      }
    }
  }

  //update the VE coboundary relation
  for(auto v : vertex_set.positions())
  {
    const auto new_v = vertex_set_map[v];
    if(new_v != INVALID_VERTEX)
    {
      // cardinality of VE relation is 1
      const auto old = ve_rel[v][0];
      ve_rel[new_v][0] = (old != INVALID_ELEMENT) ? element_set_map[old] : INVALID_ELEMENT;
    }
  }

  //update the EE adjacency relation
  for(auto e : element_set.positions())
  {
    int new_e = element_set_map[e];
    if(new_e != INVALID_ELEMENT)
    {
      const auto ee_old = ee_rel[e];
      auto ee_new = ee_rel[new_e];
      for(auto i : ee_new.positions())
      {
        const auto old = ee_old[i];
        ee_new[i] = (old != INVALID_ELEMENT) ? element_set_map[old] : INVALID_ELEMENT;
      }
    }
  }

  //Update the coordinate positions map
  for(auto v : vertex_set.positions())
  {
    int new_entry_index = vertex_set_map[v];
    if(new_entry_index != INVALID_VERTEX)
    {
      vcoord_map[new_entry_index] = vcoord_map[v];
    }
  }

  //update the sets
  vertex_set.reset(v_count);
  element_set.reset(e_count);

  ev_rel.updateSizes();
  ve_rel.updateSizes();
  ee_rel.updateSizes();
  vcoord_map.resize(v_count);
}

template <int TDIM, int SDIM, typename P>
bool IAMesh<TDIM, SDIM, P>::isEmpty() const
{
  return vertex_set.size() == 0;
}

template <int TDIM, int SDIM, typename P>
bool IAMesh<TDIM, SDIM, P>::isValid(bool verboseOutput) const
{
  fmt::memory_buffer out;

  bool bValid = true;

  //Check that sizes for vertices match
  if(vertex_set.size() != ve_rel.size() || vertex_set.size() != vcoord_map.size())
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\t vertex set and relation size(s) don't match.\n\t"
                     "vertex size: {}, ve_rel size: {}, vcoord size: {}",
                     vertex_set.size(),
                     ve_rel.size(),
                     vcoord_map.size());
    }
    bValid = false;
  }

  //Check that sizes for elements match
  if(element_set.size() != ev_rel.size() || element_set.size() != ee_rel.size())
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\t element set and relation size(s) don't match.\n\t"
                     "element_set size: {}, ev_rel size: {}, ee_rel size: {}",
                     element_set.size(),
                     ev_rel.size(),
                     ee_rel.size());
    }
    bValid = false;
  }

  // Check sets, relations and maps for validity
  if(!vertex_set.isValid(verboseOutput))
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out), "\n\t Vertex set invalid");
    }
    bValid = false;
  }
  if(!element_set.isValid(verboseOutput))
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out), "\n\t Element set invalid");
    }
    bValid = false;
  }
  if(!ev_rel.isValid(verboseOutput))
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out), "\n\t Boundary relation invalid");
    }
    bValid = false;
  }
  if(!ve_rel.isValid(verboseOutput))
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out), "\n\t Coboundary relation invalid");
    }
    bValid = false;
  }
  if(!ee_rel.isValid(verboseOutput))
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out), "\n\t Adjacency relation invalid");
    }
    bValid = false;
  }
  if(!vcoord_map.isValid(verboseOutput))
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out), "\n\t Coordinate map is invalid");
    }
    bValid = false;
  }

  if(verboseOutput)
  {
    if(bValid)
    {
      SLIC_INFO("IA mesh was valid");
    }
    else
    {
      SLIC_INFO("IA mesh was not valid.\n Summary: " << fmt::to_string(out));
    }
  }

  return bValid;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_IA_IMPL_H_
