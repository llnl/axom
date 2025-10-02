// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IA.hpp
 *
 * \brief Contains the header information of IA class
 */

#ifndef SLAM_IA_H_
#define SLAM_IA_H_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/StaticRelation.hpp"
#include "axom/slam/DynamicConstantRelation.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/DynamicMap.hpp"
#include "axom/slam/FieldRegistry.hpp"

namespace axom
{
namespace slam
{
/**
 * \class IA
 *
 * \brief Implements the Indexed mesh data structure with Adjacencies (IA)
 *
 * \tparam TDIM The topological dimension (2=triangle, 3=tetrahedra).
 * \tparam SDIM The dimension of the embedding space.
 * \tparam PointType A point type for the mesh position data
 *
 * \details The IAMesh class is an adjacency-based topological mesh data
 * structure for simplicial complexes.
 * It encodes:
 * - the boundary relation from elements to their vertices
 * - the partial coboundary relation from vertices to one incident element, and
 * - the adjacency relation between elements along their facets (faces of dimension TDIM-1).
 *
 * PointType is required to have the following interface:
 * - .ctor(T*)          -- constructor from an array of T
 * - operator[](index)  -- subscript operator to access the coordinates
 * - data() -> T*       -- direct access to the data
 * - operator==()       -- equality operator
 *
 * The basic form of this data structure was introduced in:
 *   - A. Paoluzzi, F. Bernardini, C. Cattani, and V. Ferrucci. Dimension-independent modeling with simplicial complexes.
 *     ACM Trans. Graph. 12, 1 (Jan. 1993), 56-102. DOI: https://doi.org/10.1145/169728.169719A. 
 *   - G. M. Nielson. Tools for Triangulations and Tetrahedrizations.
 *     Scientific Visualization: Overviews, Methodologies , Techniques, CS Press, pp.429-525, 1994.
 * We're using the extended formulation that includes a partial vertex coboundary relation as described in:
 *   - De Floriani L, Hui A. Data Structures for Simplicial Complexes: An Analysis And A Comparison. 
 *     Symposium on Geometry Processing 2005 Jul 5 (pp. 119-128). DOI: https://doi.org/10.2312/SGP/SGP05/119-128
 */
template <int TDIM = 2, int SDIM = 3, typename PointType = axom::slam::util::Point3<double>>
class IAMesh
{
public:
  static constexpr int COORDS_PER_VERT = SDIM;     //2 or 3 dimensional space
  static constexpr int VERTS_PER_ELEM = TDIM + 1;  //3 = triangle, 4 = tetrahedron

  using IndexType = axom::IndexType;
  using DataType = double;

  using IndexArray = std::vector<IndexType>;

  using Point = PointType;
  using PositionType = IndexType;
  using ElementType = IndexType;

  /// types for Sets
  using SetBase = slam::Set<PositionType, ElementType>;
  using VertexSet = slam::DynamicSet<PositionType, ElementType>;
  using ElementSet = slam::DynamicSet<PositionType, ElementType>;

  static constexpr int INVALID_VERTEX_INDEX = VertexSet::INVALID_ENTRY;
  static constexpr int INVALID_ELEMENT_INDEX = ElementSet::INVALID_ENTRY;

  /// types for Relations

  // A templated type alias to help with defining relations
  template <int SIZE>
  using IADynamicConstantRelation = slam::DynamicConstantRelation<
    PositionType,
    ElementType,
    slam::policies::ConstantCardinality<PositionType, slam::policies::CompileTimeStride<PositionType, SIZE>>>;

  using ElementBoundaryRelation = IADynamicConstantRelation<VERTS_PER_ELEM>;
  using VertexCoboundaryRelation = IADynamicConstantRelation<1>;
  using ElementAdjacencyRelation = IADynamicConstantRelation<VERTS_PER_ELEM>;

  using BoundarySubset = typename ElementBoundaryRelation::RelationSubset;
  using AdjacencySubset = typename ElementAdjacencyRelation::RelationSubset;

  /// types for Maps
  using PositionMap = DynamicMap<VertexSet, Point>;

  using IndexBuf = slam::FieldRegistry<SetBase, IndexType>;

  using ModularVertexIndex =
    slam::ModularInt<slam::policies::CompileTimeSize<IndexType, VERTS_PER_ELEM>>;
  using ModularFacetIndex =
    slam::ModularInt<slam::policies::CompileTimeSize<IndexType, VERTS_PER_ELEM>>;

public:
  /// \brief Default Constructor for an empty mesh
  IAMesh();

  /// \brief Constructs an IA mesh with the given point coordinates and vertex indices for each element
  IAMesh(std::vector<double>& points, std::vector<IndexType>& ev_vec);

  /// \brief Copy constructor
  IAMesh(const IAMesh&);

  IAMesh& operator=(const IAMesh&);

  /**
   * \brief check that the mesh data stored is valid.
   *
   * \note Valid meshes are not necessarily manifold.
   */
  bool isValid(bool verboseOutput = false) const;

  /// \name Accessors for encoded Sets, Relations and Maps
  /// @{

  /// \brief Returns the set of vertices in the mesh
  VertexSet& vertices() { return vertex_set; }

  /// \brief Returns the set of vertices in the mesh
  const VertexSet& vertices() const { return vertex_set; }

  /// \brief Returns the set of elements in the mesh
  ElementSet& elements() { return element_set; }

  /// \brief Returns the set of elements in the mesh
  const ElementSet& elements() const { return element_set; }

  /** 
   * \brief Returns the set of vertices in the boundary of element \a element_index
   * \pre Assumes \a element_index is a valid index in the set, i.e. 0 <= element_index < elements().size(),
   *      but does not require that the element is in the mesh.
   *      To check if the element is valid, call \a this->isValidElement(element_index)
   */
  BoundarySubset boundaryVertices(IndexType element_index) { return ev_rel[element_index]; }
  const BoundarySubset boundaryVertices(IndexType element_index) const
  {
    return ev_rel[element_index];
  }

  /** 
   * \brief Returns an element in the coboundary (star) of vertex \a vertex_index
   * \pre Assumes \a vertex_index is a valid index in the set, i.e. 0 <= vertex_index < vertices().size(),
   *      but does not require that the vertex is in the mesh.
   *      To check if the vertex is valid, call \a this->isValidVertex(vertex_index)
   * \note To return all elements in the coboundary, call \a this->vertexStar(vertex_index)
   */
  IndexType& coboundaryElement(IndexType vertex_index) { return ve_rel[vertex_index][0]; }
  const IndexType& coboundaryElement(IndexType vertex_index) const
  {
    return ve_rel[vertex_index][0];
  }

  /** 
   * \brief Returns the set of elements adjacent to element \a element_index
   * \pre Assumes \a element_index is a valid index in the set, i.e. 0 <= element_index < elements().size(),
   *      but does not require that the element is in the mesh.
   *      To check if the element is valid, call \a this->isValidElement(element_index)
   */
  AdjacencySubset adjacentElements(IndexType element_index) { return ee_rel[element_index]; }
  const AdjacencySubset adjacentElements(IndexType element_index) const
  {
    return ee_rel[element_index];
  }

  /**
   * \brief Given a vertex index, return its Point coordinate.
   */
  const Point& getVertexPosition(IndexType vertex_idx) const;

  /// @}

  /**
   * \brief Given a vertex index, return a list of incident elements, i.e. the star of the vertex
   *
   * \note If the index is invalid or out of bounds,
   * an empty list is returned. This function may return incorrect/partial
   * results if the mesh is not a manifold mesh.
   * \note This function performs a search of the local neighborhood
   * in the vicinity of the vertex using the element adjacency relation
   */
  IndexArray vertexStar(IndexType vertex_idx) const;

  /**
   * \brief Given an element index, and a face index i,
   * return a list of the vertices on the i^th face of that element
   *
   * \details If the element index is invalid or out of bounds,
   * an empty list is returned.
   */
  IndexArray getElementFace(IndexType element_idx, IndexType face_idx) const;

  /**
   * \brief Return true if the mesh has no vertex or element
   */
  bool isEmpty() const;

  /**
   * \brief Returns the number of elements in the mesh.
   * \note To get the total number of elements, including invalid and deleted elements, 
   * call \a this->elements().size()
   */
  IndexType getNumberOfValidElements() const { return element_set.numberOfValidEntries(); }

  /**
   * \brief Returns the number of vertices in the mesh
   * \note To get the total number of vertices, including invalid and deleted vertices, 
   * call \a this->vertices().size()
   */
  IndexType getNumberOfValidVertices() const { return vertex_set.numberOfValidEntries(); }

  /**
   * \brief Returns true if \a element_idx is a valid element index in the mesh
   *
   * An element index is valid when the element is not deleted.
   */
  inline bool isValidElement(IndexType element_idx) const
  {
    return element_set.isValidEntry(element_idx);
  }

  /**
   * \brief Returns true if \a vertex_idx is a valid vertex index in the mesh
   *
   * A vertex index is valid when the vertex is not deleted.
   */
  inline bool isValidVertex(IndexType vertex_idx) const
  {
    return vertex_set.isValidEntry(vertex_idx);
  }

  /**
   * \brief Fix the element adjacency relation in the neighborhood of a vertex
   *
   * \details Given a vertex index and a list of all the elements incident in
   * that vertex, fix the element->element relation data.
   * Sometimes when modifying the mesh, the mesh becomes non-manifold.
   * Adding elements may result in incorrect element->element data.
   */
  void fixVertexNeighborhood(IndexType vertex_idx, const std::vector<IndexType>& new_elements);

  /**
   * \brief Return a valid element index
   */
  IndexType getValidElementIndex() const
  {
    for(int i = element_set.size() - 1; i >= 0; --i)
    {
      if(isValidElement(i))
      {
        return i;
      }
    }
    return INVALID_ELEMENT_INDEX;
  }

  /**
   * \brief Add a vertex to the mesh,
   *
   * \param p Coordinates for the new vertex
   */
  IndexType addVertex(const Point& p);

  /**
   * \brief Add an element to the mesh.
   *
   * This is a convenience function that only uses the first VERTS_PER_ELEM
   * vertex identifiers.
   *
   * \param v0 - first vertex id of the element
   * \param v1 - second vertex id of the element
   * \param v2 - third vertex id of the element
   * \param v3 - fourth vertex id of the element
   */
  IndexType addElement(IndexType v0,
                       IndexType v1,
                       IndexType v2 = INVALID_VERTEX_INDEX,
                       IndexType v3 = INVALID_VERTEX_INDEX);

  /**
   * \brief Add an element to the mesh.
   *
   * \param vlist - A pointer to the vertex indices of the new element.
   * The array size should be at least VERTS_PER_ELEM
   */
  IndexType addElement(const IndexType* vlist);

  /**
   * \brief Add an element to the mesh. Overload to also set neighbors
   *
   * \param vlist - A pointer to the vertex indices of the new element.
   * The array size should be at least VERTS_PER_ELEM
   * \param neighbors - A pointer to the neighbor indices of the new element.
   * The array size should be at least VERTS_PER_ELEM
   */
  IndexType addElement(const IndexType* vlist, const IndexType* neighbors);

  /**
   * \brief Removes an element from the mesh
   *
   * \details If the index is invalid or out of bounds,
   * no changes are made to the mesh.
   *
   * \param element_idx The index of the element to remove.
   * \warning Removing an element could make one of its vertices non-manifold
   * (its link will have more than one connected components) since the
   * IA only retains 1 element reference per vertex.
   */
  void removeElement(IndexType element_idx);

  /**
   * \brief Removes a vertex and all incident elements from the mesh
   *
   * \details If the index is invalid or out of bounds,
   * no changes are made to the mesh.
   * \param vertex_idx The index of the vertex to remove.
   */
  void removeVertex(IndexType vertex_idx);

  /**
   * \brief Removes all the invalid entries in the mesh and reduce memory used
   * \details This function may invalidates all indices in user code.
   */
  void compact();

  /**
   * \brief Prints the IA mesh structure, for debug purpose.
   */
  void print_all() const;

private:
  //helper functions to find element to element relations
  using ElementIndexType = PositionType;
  using FaceIndexType = PositionType;
  using ElementAndFaceIdxType = std::pair<ElementIndexType, FaceIndexType>;
  using V2EMapType = std::map<IndexArray, ElementAndFaceIdxType>;

  /**
   * \brief Helper function to find adjacent elements when adding a new element
   *
   * \param vertpair_to_elem_map A map from facets encoded using a sorted tuple of
   * vertex IDs to the face encoded as a element index and a face index
   * \param element_i The index of the element whose neighbor we are seeking
   * \param side_i The local index of the desired face w.r.t. \a element_i
   */
  ElementAndFaceIdxType ElemNbrFinder(V2EMapType& vertpair_to_elem_map,
                                      IndexType element_i,
                                      IndexType side_i);

private:
  VertexSet vertex_set;             //Set of vertices
  ElementSet element_set;           //Set of elements
  ElementBoundaryRelation ev_rel;   //Element to vertex relation.
  VertexCoboundaryRelation ve_rel;  //Vertex to one element partial relation.
  ElementAdjacencyRelation ee_rel;  //Element to neighboring element relation
  PositionMap vcoord_map;           //map of coordinates per vertex.
};

template <int TDIM, int SDIM, typename P>
constexpr int IAMesh<TDIM, SDIM, P>::COORDS_PER_VERT;

template <int TDIM, int SDIM, typename P>
constexpr int IAMesh<TDIM, SDIM, P>::VERTS_PER_ELEM;

}  // end namespace slam
}  // end namespace axom

#include "axom/slam/mesh_struct/IA_impl.hpp"

#endif  //  SLAM_IA_H_
