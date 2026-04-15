// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include "axom/core/StaticArray.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/ModularInt.hpp"

#include "axom/fmt.hpp"

#include <vector>
#include <map>
#include <array>
#include <algorithm>
#include <memory>
#include <cstddef>

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
      if(nbr != INVALID_ELEMENT_INDEX  //
         && !is_subset(nbr, ret)       //
         && is_subset(vertex_idx, ev_rel[nbr]))
      {
        SLIC_ASSERT(element_set.isValidEntry(nbr));
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
        if(nbr != INVALID_ELEMENT_INDEX && is_subset(vertex_i, ev_rel[nbr]))
        {
          SLIC_ASSERT(element_set.isValidEntry(nbr));
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
    if(nbr != INVALID_ELEMENT_INDEX)
    {
      SLIC_ASSERT(isValidElement(nbr));
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
typename IAMesh<TDIM, SDIM, P>::IndexType IAMesh<TDIM, SDIM, P>::reuseElement(IndexType element_idx,
                                                                              const IndexType* vlist,
                                                                              const IndexType* neighbors)
{
  SLIC_ASSERT_MSG(element_idx >= 0 && element_idx < element_set.size(),
                  "Trying to reuse an out-of-range element index:" << element_idx);
  SLIC_ASSERT_MSG(!element_set.isValidEntry(element_idx),
                  "Trying to reuse an element slot that is already valid:" << element_idx);

  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    SLIC_ASSERT_MSG(vertex_set.isValidEntry(vlist[i]),
                    "Trying to reuse an element with invalid vertex index:" << vlist[i]);
  }

  element_set[element_idx] = element_idx;

  auto bdry = ev_rel[element_idx];
  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    bdry[i] = vlist[i];
  }

  auto adj = ee_rel[element_idx];
  for(int i = 0; i < VERTS_PER_ELEM; ++i)
  {
    adj[i] = neighbors[i];
  }

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
  constexpr IndexType EMPTY_SLOT = INVALID_ELEMENT_INDEX;
  constexpr IndexType TOMBSTONE_SLOT = INVALID_ELEMENT_INDEX - 1;

  struct PendingFace
  {
    IndexType key0 {INVALID_VERTEX_INDEX};
    IndexType key1 {INVALID_VERTEX_INDEX};
    IndexType element_idx {EMPTY_SLOT};
    IndexType face_idx {INVALID_ELEMENT_INDEX};
  };

  static thread_local std::vector<PendingFace> pending_faces;
  static thread_local std::vector<std::size_t> used_slots;
  for(const auto slot : used_slots)
  {
    pending_faces[slot].element_idx = EMPTY_SLOT;
  }
  used_slots.clear();

  std::size_t table_size = 8;
  const std::size_t target_slots = std::max<std::size_t>(8, 4 * new_elements.size());
  while(table_size < target_slots)
  {
    table_size <<= 1;
  }
  if(pending_faces.size() < table_size)
  {
    pending_faces.resize(table_size);
  }

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

  auto getPendingFaceKey = [vertex_idx](const BoundarySubset& bdry, IndexType vert_i) {
    PendingFace key;
    for(int i = 0, idx = 0; i < TDIM; ++i)
    {
      if(i != vert_i && bdry[i] != vertex_idx)
      {
        if(idx == 0)
        {
          key.key0 = bdry[i];
        }
        else
        {
          key.key1 = bdry[i];
        }
        ++idx;
      }
    }

    if constexpr(TDIM == 3)
    {
      if(key.key1 < key.key0)
      {
        axom::utilities::swap(key.key0, key.key1);
      }
    }

    return key;
  };

  auto pendingFaceHash = [](IndexType key0, IndexType key1) -> std::size_t {
    std::size_t seed = static_cast<std::size_t>(key0);
    seed ^= static_cast<std::size_t>(key1) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  };

  int num_pending_faces = 0;
  int num_incident_faces = 0;

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
        if(nbr != INVALID_ELEMENT_INDEX)
        {
          SLIC_ASSERT(element_set.isValidEntry(nbr));
          const auto nbr_ev = ev_rel[nbr];
          auto nbr_ee = ee_rel[nbr];

          for(auto face_j : nbr_ev.positions())
          {
            const IndexType nbr_face = nbr_ee[face_j];
            SLIC_ASSERT(nbr_face == INVALID_ELEMENT_INDEX || element_set.isValidEntry(nbr_face));
            if(nbr_face == INVALID_ELEMENT_INDEX && isSharedFace(nbr_ev, face_j, bdry))
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
        const PendingFace face = getPendingFaceKey(bdry, vert_i);
        const std::size_t mask = pending_faces.size() - 1;
        std::size_t slot = pendingFaceHash(face.key0, face.key1) & mask;
        std::size_t insert_slot = pending_faces.size();

        for(;; slot = (slot + 1) & mask)
        {
          auto& pending = pending_faces[slot];
          if(pending.element_idx == EMPTY_SLOT)
          {
            if(insert_slot == pending_faces.size())
            {
              insert_slot = slot;
            }

            auto& insert_entry = pending_faces[insert_slot];
            insert_entry = face;
            insert_entry.element_idx = el;
            insert_entry.face_idx = face_i;
            used_slots.push_back(insert_slot);
            ++num_pending_faces;
            break;
          }

          if(pending.element_idx == TOMBSTONE_SLOT)
          {
            if(insert_slot == pending_faces.size())
            {
              insert_slot = slot;
            }
            continue;
          }

          if(pending.key0 == face.key0 && pending.key1 == face.key1)
          {
            SLIC_ASSERT_MSG(pending.element_idx != el,
                            "Each face in the inserted star should be shared by two elements");
            ee_rel.modify(pending.element_idx, pending.face_idx, el);
            ee_rel.modify(el, face_i, pending.element_idx);
            pending.element_idx = TOMBSTONE_SLOT;
            --num_pending_faces;
            break;
          }
        }
        ++num_incident_faces;
      }
    }
  }

  AXOM_UNUSED_VAR(num_incident_faces);
  AXOM_UNUSED_VAR(num_pending_faces);
  SLIC_ASSERT(num_incident_faces == static_cast<int>(new_elements.size()) * TDIM);
  SLIC_ASSERT_MSG(num_pending_faces == 0,
                  "All faces in the inserted star should be paired exactly once");
}

// Remove all the invalid entries in the IA structure
template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::compact()
{
  constexpr IndexType INVALID_VERTEX = VertexSet::INVALID_ENTRY;
  constexpr IndexType INVALID_ELEMENT = ElementSet::INVALID_ENTRY;
  const IndexType vertex_size = vertex_set.size();
  const IndexType element_size = element_set.size();
  const auto& vertex_data = vertex_set.data();
  const auto& element_data = element_set.data();
  auto& ev_data = ev_rel.data();
  auto& ve_data = ve_rel.data();
  auto& ee_data = ee_rel.data();

  bool has_invalid_vertices = false;
  IndexType v_count = 0;
  for(IndexType v = 0; v < vertex_size; ++v)
  {
    has_invalid_vertices |= vertex_data[v] == INVALID_VERTEX;
    v_count += (vertex_data[v] != INVALID_VERTEX) ? 1 : 0;
  }

  bool has_invalid_elements = false;
  std::unique_ptr<IndexType[]> element_set_map(new IndexType[element_size]);
  IndexType e_count = 0;
  for(IndexType e = 0; e < element_size; ++e)
  {
    if(element_data[e] != INVALID_ELEMENT)
    {
      element_set_map[e] = e_count++;
    }
    else
    {
      has_invalid_elements = true;
    }
  }

  if(!has_invalid_vertices && !has_invalid_elements)
  {
    return;
  }

  auto remapElement = [&](IndexType old_element) {
    return (old_element >= 0 && old_element < element_size &&
            element_data[old_element] != INVALID_ELEMENT)
      ? element_set_map[old_element]
      : INVALID_ELEMENT;
  };

  if(!has_invalid_vertices)
  {
    for(IndexType e = 0; e < element_size; ++e)
    {
      if(element_data[e] == INVALID_ELEMENT)
      {
        continue;
      }

      const IndexType new_e = element_set_map[e];
      const IndexType old_base = e * VERTS_PER_ELEM;
      const IndexType new_base = new_e * VERTS_PER_ELEM;
      for(int i = 0; i < VERTS_PER_ELEM; ++i)
      {
        ev_data[new_base + i] = ev_data[old_base + i];
        ee_data[new_base + i] = remapElement(ee_data[old_base + i]);
      }
    }

    for(IndexType v = 0; v < vertex_size; ++v)
    {
      ve_data[v] = remapElement(ve_data[v]);
    }

    element_set.reset(e_count);
    ev_rel.updateSizes();
    ee_rel.updateSizes();
    return;
  }

  std::unique_ptr<IndexType[]> vertex_set_map(new IndexType[vertex_size]);
  v_count = 0;
  for(IndexType v = 0; v < vertex_size; ++v)
  {
    if(vertex_data[v] != INVALID_VERTEX)
    {
      vertex_set_map[v] = v_count++;
    }
  }

  auto remapVertex = [&](IndexType old_vertex) {
    return (old_vertex >= 0 && old_vertex < vertex_size && vertex_data[old_vertex] != INVALID_VERTEX)
      ? vertex_set_map[old_vertex]
      : INVALID_VERTEX;
  };

  for(IndexType e = 0; e < element_size; ++e)
  {
    if(element_data[e] == INVALID_ELEMENT)
    {
      continue;
    }

    const IndexType new_e = element_set_map[e];
    const IndexType old_base = e * VERTS_PER_ELEM;
    const IndexType new_base = new_e * VERTS_PER_ELEM;
    for(int i = 0; i < VERTS_PER_ELEM; ++i)
    {
      ev_data[new_base + i] = remapVertex(ev_data[old_base + i]);
      ee_data[new_base + i] = remapElement(ee_data[old_base + i]);
    }
  }

  auto& coord_data = vcoord_map.data();
  for(IndexType v = 0; v < vertex_size; ++v)
  {
    if(vertex_data[v] == INVALID_VERTEX)
    {
      continue;
    }

    const IndexType new_v = vertex_set_map[v];
    ve_data[new_v] = remapElement(ve_data[v]);
    coord_data[new_v] = coord_data[v];
  }

  vertex_set.reset(v_count);
  element_set.reset(e_count);

  ev_rel.updateSizes();
  ve_rel.updateSizes();
  ee_rel.updateSizes();
  vcoord_map.resize(v_count);
}

template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::reserveVertices(IndexType vertex_capacity)
{
  vertex_set.reserve(vertex_capacity);
  ve_rel.reserve(vertex_capacity);
  vcoord_map.reserve(vertex_capacity);
}

template <int TDIM, int SDIM, typename P>
void IAMesh<TDIM, SDIM, P>::reserveElements(IndexType element_capacity)
{
  element_set.reserve(element_capacity);
  ev_rel.reserve(element_capacity);
  ee_rel.reserve(element_capacity);
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

template <int TDIM, int SDIM, typename P>
typename IAMesh<TDIM, SDIM, P>::FacetKey IAMesh<TDIM, SDIM, P>::getSortedFacetKey(
  IndexType element_idx,
  IndexType facet_idx) const
{
  FacetKey key {};
  if(!element_set.isValidEntry(element_idx))
  {
    return key;
  }

  SLIC_ASSERT_MSG(0 <= facet_idx && facet_idx < VERTS_PER_ELEM, "Face index is invalid.");

  const auto verts = ev_rel[element_idx];
  ModularVertexIndex mod_face(facet_idx);
  for(int i = 0; i < VERTS_PER_ELEM - 1; ++i)
  {
    key[i] = verts[mod_face + i];
  }

  std::sort(key.begin(), key.end());
  return key;
}

template <int TDIM, int SDIM, typename P>
bool IAMesh<TDIM, SDIM, P>::isConforming(bool verboseOutput) const
{
  fmt::memory_buffer out;

  bool valid = isValid(verboseOutput);

  struct FacetBucket
  {
    axom::StaticArray<FacetRecord, 2> records;
    int incident_count {0};
  };

  std::map<FacetKey, FacetBucket> facet_records;

  auto facetKeyString = [](const FacetKey& facet_key) {
    return fmt::format("[{}]", fmt::join(facet_key, ", "));
  };

  for(auto element_idx : elements().positions())
  {
    if(!isValidElement(element_idx))
    {
      continue;
    }

    // check that element vertices are all valid and non-repeating
    const auto verts = boundaryVertices(element_idx);
    for(int i = 0; i < VERTS_PER_ELEM; ++i)
    {
      if(!isValidVertex(verts[i]))
      {
        if(verboseOutput)
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\tElement {} references invalid vertex {}",
                         element_idx,
                         verts[i]);
        }
        valid = false;
      }

      for(int j = i + 1; j < VERTS_PER_ELEM; ++j)
      {
        if(verts[i] == verts[j])
        {
          if(verboseOutput)
          {
            fmt::format_to(std::back_inserter(out),
                           "\n\tElement {} repeats vertex {}",
                           element_idx,
                           verts[i]);
          }
          valid = false;
        }
      }
    }

    // build facet-element co-boundary in facet_records
    const auto neighbors = adjacentElements(element_idx);
    for(int facet_idx = 0; facet_idx < VERTS_PER_ELEM; ++facet_idx)
    {
      const FacetKey facet_key = getSortedFacetKey(element_idx, facet_idx);
      FacetBucket& bucket = facet_records[facet_key];
      bucket.incident_count++;
      if(bucket.incident_count <= 2)
      {
        bucket.records.push_back({element_idx, facet_idx, neighbors[facet_idx]});
      }
      else
      {
        if(verboseOutput && bucket.incident_count == 3)
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\tFacet {} is non-manifold (>=3 incident elements); "
                         "first two are ({}:{}) and ({}:{})",
                         facetKeyString(facet_key),
                         bucket.records[0].element_idx,
                         bucket.records[0].facet_idx,
                         bucket.records[1].element_idx,
                         bucket.records[1].facet_idx);
        }
        valid = false;
      }
    }
  }

  // check for valid facet-coboundary relation
  // each facet is referenced once if it is on mesh boundary;
  // and twice otherwise, with consistent adjacencies
  for(const auto& [facet_key, bucket] : facet_records)
  {
    if(bucket.incident_count == 1)
    {
      const FacetRecord& record = bucket.records[0];
      if(isValidElement(record.neighbor_idx))
      {
        if(verboseOutput)
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\tBoundary face {} on element {} facet {} points to neighbor {}",
                         facetKeyString(facet_key),
                         record.element_idx,
                         record.facet_idx,
                         record.neighbor_idx);
        }
        valid = false;
      }
    }
    else if(bucket.incident_count == 2)
    {
      const FacetRecord& lhs = bucket.records[0];
      const FacetRecord& rhs = bucket.records[1];
      if(lhs.element_idx == rhs.element_idx || lhs.neighbor_idx != rhs.element_idx ||
         rhs.neighbor_idx != lhs.element_idx)
      {
        if(verboseOutput)
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\tInterior facet {} has inconsistent adjacency: "
                         "({}:{}) -> {}, ({}:{}) -> {}",
                         facetKeyString(facet_key),
                         lhs.element_idx,
                         lhs.facet_idx,
                         lhs.neighbor_idx,
                         rhs.element_idx,
                         rhs.facet_idx,
                         rhs.neighbor_idx);
        }
        valid = false;
      }
    }
    else
    {
      // bucket.incident_count > 2 was already reported, if requested
    }
  }

  if(verboseOutput)
  {
    if(valid)
    {
      SLIC_INFO("IA mesh was conforming");
    }
    else
    {
      SLIC_INFO("IA mesh was not conforming.\n Summary: " << fmt::to_string(out));
    }
  }

  return valid;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_IA_IMPL_H_
