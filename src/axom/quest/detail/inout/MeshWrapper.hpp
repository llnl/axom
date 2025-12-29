// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MeshWrapper.hpp
 *
 * \brief Defines a templated mesh wrapper class for the InOutOctree.
 */

#ifndef AXOM_QUEST_INOUT_OCTREE_MESH_WRAPPER__HPP_
#define AXOM_QUEST_INOUT_OCTREE_MESH_WRAPPER__HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
/**
   * \brief A utility class that wraps access to the mesh data of an InOutOctree
   *
   * This class helps separate the specifics of accessing the underlying mesh
   * for an InOutOctree. It is customized for unstructured Segment meshes in 2D
   * and Triangle meshes in 3D. 
   *
   * If we want to support other surface mesh types (e.g. quad meshes in 3D),
   * we'll have to customize it a bit more.
   *
   * \note Uses the CRTP pattern to allow the dimension-specific derived classes
   * to share the dimension-independent implementation
   */
template <int DIM>
class MeshWrapper;

template <int DIM, typename Derived>
class SimplexMeshWrapper
{
public:
  using VertexIndex = axom::IndexType;
  using CellIndex = axom::IndexType;

  using MeshVertexSet = slam::PositionSet<>;
  using MeshElementSet = slam::PositionSet<>;
  using SurfaceMesh = mint::Mesh;

  using SpacePt = primal::Point<double, DIM>;
  using SpaceVector = primal::Vector<double, DIM>;
  using GeometricBoundingBox = primal::BoundingBox<double, DIM>;

  using VertexIndexMap = slam::Map<VertexIndex>;
  using VertexPositionMap = slam::Map<SpacePt>;

  /// Always DIM verts since we're representing a d-dimensional simplicial mesh in dimension d
  static constexpr int NUM_CELL_VERTS = DIM;
  using STLIndirection = slam::policies::STLVectorIndirection<VertexIndex, VertexIndex>;
  using TVStride = slam::policies::CompileTimeStride<VertexIndex, NUM_CELL_VERTS>;
  using ConstantCardinality = slam::policies::ConstantCardinality<VertexIndex, TVStride>;
  using CellVertexRelation =
    slam::StaticRelation<VertexIndex, VertexIndex, ConstantCardinality, STLIndirection, MeshElementSet, MeshVertexSet>;
  using CellVertIndices = typename CellVertexRelation::RelationSubset;

  // \brief A vertex index to indicate that there is no associated vertex
  static constexpr VertexIndex NO_VERTEX = -1;

  /// \brief A vertex index to indicate that there is no associated vertex
  static constexpr CellIndex NO_CELL = -1;

protected:
  SimplexMeshWrapper(std::shared_ptr<SurfaceMesh>& meshPtr)
    : m_surfaceMesh(meshPtr)
    , m_vertexPositions(&m_vertexSet)
    , m_cellToVertexRelation()
  { }

private:
  /// Utility functions to get a pointer to the derived type (part of CRTP pattern)
  Derived* getDerived() { return static_cast<Derived*>(this); }
  /// Utility functions to get a const pointer to the derived type (part of CRTP pattern)
  const Derived* getDerived() const { return static_cast<const Derived*>(this); }

public:
  /** Predicate to determine if the wrapped surface mesh has been reindexed */
  bool meshWasReindexed() const { return m_meshWasReindexed; }

  /** Const accessor to the vertex set of the wrapped surface mesh */
  const MeshVertexSet& vertexSet() const { return m_vertexSet; }

  /** Accessor to the vertex set of the wrapped surface mesh */
  MeshVertexSet& vertexSet() { return m_vertexSet; }

  /** Const accessor to the element set of the wrapped surface mesh */
  const MeshElementSet& elementSet() const { return m_elementSet; }

  /** Accessor to the element set of the wrapped surface mesh */
  MeshElementSet& elementSet() { return m_elementSet; }

  /** Accessor for the number of vertices in the wrapped surface mesh */
  int numMeshVertices() const
  {
    if(m_meshWasReindexed)
    {
      return m_vertexSet.size();
    }
    else
    {
      return m_surfaceMesh->getNumberOfNodes();
    }
  }

  /** Accessor for the number of elements in the wrapped surface mesh */
  int numMeshCells() const
  {
    return m_meshWasReindexed ? static_cast<int>(m_elementSet.size())
                              : static_cast<int>(m_surfaceMesh->getNumberOfCells());
  }

  /**
   * \brief Returns position of vertex with index \a idx within the surface mesh
   *
   * \note Use this function instead of vertexPosition() if calling in a context
   * where the mesh might not yet have been reindexed
   * \sa vertexPosition()
   */
  SpacePt getMeshVertexPosition(VertexIndex idx) const
  {
    if(m_meshWasReindexed)
    {
      return m_vertexPositions[idx];
    }
    else
    {
      SLIC_ASSERT(m_surfaceMesh->getDimension() == SpacePt::dimension());

      SpacePt pt;
      m_surfaceMesh->getNode(idx, pt.data());

      return pt;
    }
  }

  /**
   * \brief Returns spatial position of vertex with index \a idx from wrapped surface mesh
   * 
   * \note Use after mesh has been reindexed
   */
  const SpacePt& vertexPosition(VertexIndex idx) const { return m_vertexPositions[idx]; }

  /**
     * \brief Returns the indices of the boundary vertices of the element
     * of the wrapped surface mesh with the given index
     *
     * \param idx The index of an element within the surface mesh
     */
  CellVertIndices cellVertexIndices(CellIndex idx) const { return m_cellToVertexRelation[idx]; }

  /**
   * \brief Finds the index of a vertex in cell \a c1 that is not in cell \a c0
   *
   * \param c0 The index of the first cell
   * \param c1 The index of the second cell
   * \pre \a c0 and \a c1 must be distinct cells
   * \return The index of a vertex in \a c1 that is not in \a c0; NO_VERTEX if one does not exist
   */
  VertexIndex distinctVertex(CellIndex c0, CellIndex c1) const
  {
    SLIC_ASSERT_MSG(c0 != c1,
                    "Expected two different cell indices in "
                    "quest::MeshWrapper::findDistinctVertex,"
                      << " got cell index " << c0 << " twice.");

    // Find a vertex from the local surface that is not incident in the first cell
    CellVertIndices cvRel0 = cellVertexIndices(c0);
    CellVertIndices cvRel1 = cellVertexIndices(c1);

    for(int i = 0; i < NUM_CELL_VERTS; ++i)
    {
      if(!incidentInVertex(cvRel0, cvRel1[i]))
      {
        return cvRel1[i];
      }
    }

    SLIC_ASSERT_MSG(
      false,
      fmt::format("There should be a vertex in cell {} that was not in cell {}", c1, c0));
    return NO_VERTEX;
  }

  /**
   * \brief Determine if the two given cells have a vertex in common
   *
   * \param c0 The index of the first cell
   * \param c1 The index of the second cell
   * \param [out] sharedVert The index of the shared vertex, if it exists
   * \return true if the two cells have a vertex
   * in common (returned in sharedVert), false otherwise
   */
  bool haveSharedVertex(CellIndex c0, CellIndex c1, VertexIndex& sharedVert) const
  {
    // There are two cells  -- check that they have at least one common vertex
    CellVertIndices cvRel0 = cellVertexIndices(c0);
    CellVertIndices cvRel1 = cellVertexIndices(c1);

    for(int i = 0; i < NUM_CELL_VERTS; ++i)
    {
      if(getDerived()->incidentInVertex(cvRel0, cvRel1[i]))
      {
        sharedVert = cvRel1[i];
        return true;
      }
    }
    return false;
  }

  /**
   * \brief Determine if the three given cells have a vertex in common
   *
   * \param c0 The index of the first cell
   * \param c1 The index of the second cell
   * \param c2 The index of the second cell
   * \param [out] sharedVert The index of the shared vertex, if it exists
   * \return true if the three cells have a vertex in common
   *  (returned in sharedVert), false otherwise
   */
  bool haveSharedVertex(CellIndex c0, CellIndex c1, CellIndex c2, VertexIndex& sharedVert) const
  {
    CellVertIndices c0Verts = cellVertexIndices(c0);
    CellVertIndices c1Verts = cellVertexIndices(c1);
    CellVertIndices c2Verts = cellVertexIndices(c2);

    for(int i = 0; i < NUM_CELL_VERTS; ++i)
    {
      // check if a vertex from the third cell is in the first and second
      if(getDerived()->incidentInVertex(c0Verts, c2Verts[i]) &&
         getDerived()->incidentInVertex(c1Verts, c2Verts[i]))
      {
        sharedVert = c2Verts[i];
        return true;
      }
    }

    return false;
  }

protected:
  std::shared_ptr<SurfaceMesh>& m_surfaceMesh;  // ref to pointer to allow changing the mesh

  MeshVertexSet m_vertexSet {0};
  MeshElementSet m_elementSet {0};

  VertexPositionMap m_vertexPositions;

  std::vector<VertexIndex> m_cv_data;
  CellVertexRelation m_cellToVertexRelation;

  bool m_meshWasReindexed {false};
};

template <>
class MeshWrapper<2> : public SimplexMeshWrapper<2, MeshWrapper<2>>
{
public:
  static constexpr int DIM = 2;
  using Base = SimplexMeshWrapper<DIM, MeshWrapper<DIM>>;
  using Base::CellIndex;
  using Base::GeometricBoundingBox;
  using Base::SpacePt;
  using Base::SpaceVector;
  using Base::SurfaceMesh;
  using Base::VertexIndex;
  using Base::VertexIndexMap;
  using Base::VertexPositionMap;

  /// \brief A constant for the number of boundary vertices in an edge */
  static constexpr int NUM_EDGE_VERTS = DIM;
  using SpaceCell = primal::Segment<double, NUM_EDGE_VERTS>;

public:
  /// \brief Constructor for a mesh wrapper */
  MeshWrapper(std::shared_ptr<SurfaceMesh>& meshPtr) : Base(meshPtr) { }

  /**
   * \brief Helper function to compute the bounding box of an edge
   *
   * \param idx The edges's index within the surface mesh
   */
  GeometricBoundingBox cellBoundingBox(CellIndex idx) const
  {
    // Get the ids of the verts bounding this edge
    CellVertIndices vertIds = cellVertexIndices(idx);
    return GeometricBoundingBox(vertexPosition(vertIds[0]), vertexPosition(vertIds[1]));
  }

  /**
   * \brief Utility function to retrieve the positions of the edge's vertices
   *
   * \return A SpaceCell (Segment) whose vertices are positioned in space
   */
  SpaceCell cellPositions(CellIndex idx) const
  {
    CellVertIndices verts = cellVertexIndices(idx);
    return SpaceCell(vertexPosition(verts[0]), vertexPosition(verts[1]));
  }

  /// \brief Checks whether the indexed segment contains a reference to the given vertex
  bool incidentInVertex(const CellVertIndices& ids, VertexIndex vIdx) const
  {
    return (ids[0] == vIdx) || (ids[1] == vIdx);
  }

  /**
   * \brief Returns the normal vector of the surface for the cell with index \a cidx
   * 
   * If we are at an endpoint of the segment (i.e. if the \a segmentParameter is close
   * to 0 or 1), we compute the average normal of its incident segments
   */
  template <typename CellIndexSet>
  SpaceVector surfaceNormal(CellIndex cidx, double segmentParameter, const CellIndexSet& otherCells) const
  {
    SpaceVector vec = this->cellPositions(cidx).normal();

    // Check if the point is at the first vertex of the segment
    if(axom::utilities::isNearlyEqual(segmentParameter, 0.))
    {
      vec = vec.unitVector();
      for(auto idx : otherCells)
      {
        auto vidx = cellVertexIndices(cidx)[0];
        if(idx != cidx && incidentInVertex(cellVertexIndices(idx), vidx))
        {
          vec += this->cellPositions(idx).normal().unitVector();
        }
      }
    }
    // Check if the point is at the second vertex of the segment
    else if(axom::utilities::isNearlyEqual(segmentParameter, 1.))
    {
      vec = vec.unitVector();
      for(auto idx : otherCells)
      {
        auto vidx = cellVertexIndices(cidx)[1];
        if(idx != cidx && incidentInVertex(cellVertexIndices(idx), vidx))
        {
          vec += this->cellPositions(idx).normal().unitVector();
        }
      }
    }

    return vec.unitVector();
  }

  /**
   * \brief Reindexes the mesh vertices and edge indices using the given map
   *
   * \param numVertices The number of vertices in the new mesh
   * \param vertexIndexMap A mapping from the old vertex indices to the new ones
   * \note This step clears out the original mesh,
   * which can be reconstructed using the regenerateSurfaceMesh() function
   */
  void reindexMesh(int numVertices, const VertexIndexMap& vertexIndexMap)
  {
    // Create a vertex set on the new vertices and grab coordinates from the old ones
    m_vertexSet = MeshVertexSet(numVertices);
    m_vertexPositions = VertexPositionMap(&m_vertexSet);

    int numOrigVertices = numMeshVertices();
    for(int i = 0; i < numOrigVertices; ++i)
    {
      const VertexIndex& vInd = vertexIndexMap[i];
      m_vertexPositions[vInd] = getMeshVertexPosition(i);
    }

    // Update the vertex IDs of the triangles to the new vertices
    // and create a SLAM relation on these
    int numOrigEdges = numMeshCells();

    m_cv_data.clear();
    m_cv_data.reserve(NUM_EDGE_VERTS * numOrigEdges);
    for(axom::IndexType i = 0; i < numOrigEdges; ++i)
    {
      // Grab relation from mesh
      using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
      axom::IndexType* vertIds = std::static_pointer_cast<UMesh>(m_surfaceMesh)->getCellNodeIDs(i);

      // Remap the vertex IDs
      for(int j = 0; j < NUM_EDGE_VERTS; ++j)
      {
        vertIds[j] = vertexIndexMap[vertIds[j]];
      }

      // Add to relation if not degenerate edge
      // (namely, we need 2 unique vertex IDs)
      if((vertIds[0] != vertIds[1]))
      {
        m_cv_data.push_back(vertIds[0]);
        m_cv_data.push_back(vertIds[1]);
      }
    }

    m_elementSet = MeshElementSet(static_cast<int>(m_cv_data.size()) / NUM_EDGE_VERTS);
    m_cellToVertexRelation = CellVertexRelation(&m_elementSet, &m_vertexSet);
    m_cellToVertexRelation.bindIndices(static_cast<int>(m_cv_data.size()), &m_cv_data);

    // Delete old mesh
    m_surfaceMesh.reset();

    m_meshWasReindexed = true;
  }

  /// \brief Add the vertex positions and edge boundary relations to the surface mesh
  void regenerateSurfaceMesh()
  {
    if(m_surfaceMesh != nullptr)
    {
      m_surfaceMesh.reset();
    }

    using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
    UMesh* edgeMesh = new UMesh(DIM, mint::SEGMENT, m_vertexSet.size(), m_elementSet.size());

    // Add vertices to the mesh (i.e. vertex positions)
    for(int i = 0; i < m_vertexSet.size(); ++i)
    {
      const SpacePt& pt = vertexPosition(i);
      edgeMesh->appendNode(pt[0], pt[1]);
    }

    // Add edges to the mesh (i.e. boundary vertices)
    for(int i = 0; i < m_elementSet.size(); ++i)
    {
      const CellIndex* tv = &cellVertexIndices(i)[0];
      edgeMesh->appendCell(tv);
    }

    m_surfaceMesh.reset(edgeMesh);
  }
};

template <>
class MeshWrapper<3> : public SimplexMeshWrapper<3, MeshWrapper<3>>
{
public:
  static constexpr int DIM = 3;
  using Base = SimplexMeshWrapper<DIM, MeshWrapper<DIM>>;
  using Base::CellIndex;
  using Base::GeometricBoundingBox;
  using Base::SpacePt;
  using Base::SurfaceMesh;
  using Base::VertexIndex;
  using Base::VertexIndexMap;
  using Base::VertexPositionMap;

  /// \brief A constant for the number of boundary vertices in a triangle */
  static constexpr int NUM_TRI_VERTS = 3;
  using SpaceCell = primal::Triangle<double, DIM>;

public:
  /// \brief Constructor for a mesh wrapper */
  MeshWrapper(std::shared_ptr<SurfaceMesh>& meshPtr) : Base(meshPtr) { }

  /**
   * \brief Helper function to compute the bounding box of a triangle
   *
   * \param idx The triangle's index within the surface mesh
   */
  GeometricBoundingBox cellBoundingBox(CellIndex idx) const
  {
    // Get the ids of the verts bounding this triangle
    CellVertIndices vertIds = cellVertexIndices(idx);

    GeometricBoundingBox bb(vertexPosition(vertIds[0]));
    bb.addPoint(vertexPosition(vertIds[1]));
    bb.addPoint(vertexPosition(vertIds[2]));

    return bb;
  }

  /**
   * \brief Utility function to retrieve the positions of the triangle's vertices
   *
   * \return A SpaceCell (Triangle) whose vertices are positioned in space
   */
  SpaceCell cellPositions(CellIndex idx) const
  {
    CellVertIndices verts = cellVertexIndices(idx);
    return SpaceCell(vertexPosition(verts[0]), vertexPosition(verts[1]), vertexPosition(verts[2]));
  }

  /// \brief Checks whether the indexed triangle contains a reference to the given vertex
  bool incidentInVertex(const CellVertIndices& ids, VertexIndex vIdx) const
  {
    return (ids[0] == vIdx) || (ids[1] == vIdx) || (ids[2] == vIdx);
  }

  /**
   * \brief Reindexes the mesh vertices and triangle indices using the given map
   *
   * \param numVertices The number of vertices in the new mesh
   * \param vertexIndexMap A mapping from the old vertex indices to the new ones
   * \note This step clears out the original mesh,
   * which can be reconstructed using the regenerateSurfaceMesh() function
   */
  void reindexMesh(int numVertices, const VertexIndexMap& vertexIndexMap)
  {
    // Create a vertex set on the new vertices and grab coordinates from the old ones
    m_vertexSet = MeshVertexSet(numVertices);
    m_vertexPositions = VertexPositionMap(&m_vertexSet);

    int numOrigVertices = numMeshVertices();
    for(int i = 0; i < numOrigVertices; ++i)
    {
      const VertexIndex& vInd = vertexIndexMap[i];
      m_vertexPositions[vInd] = getMeshVertexPosition(i);
    }

    // Update the vertex IDs of the triangles to the new vertices
    // and create a SLAM relation on these
    int numOrigTris = numMeshCells();

    m_cv_data.clear();
    m_cv_data.reserve(NUM_TRI_VERTS * numOrigTris);
    for(axom::IndexType i = 0; i < numOrigTris; ++i)
    {
      // Grab relation from mesh
      using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
      axom::IndexType* vertIds = std::static_pointer_cast<UMesh>(m_surfaceMesh)->getCellNodeIDs(i);

      // Remap the vertex IDs
      for(int j = 0; j < NUM_TRI_VERTS; ++j)
      {
        vertIds[j] = vertexIndexMap[vertIds[j]];
      }

      // Add to relation if not degenerate triangles
      // (namely, we need 3 unique vertex IDs)
      if((vertIds[0] != vertIds[1]) && (vertIds[1] != vertIds[2]) && (vertIds[2] != vertIds[0]))
      {
        m_cv_data.push_back(vertIds[0]);
        m_cv_data.push_back(vertIds[1]);
        m_cv_data.push_back(vertIds[2]);
      }
    }

    m_elementSet = MeshElementSet(static_cast<int>(m_cv_data.size()) / NUM_TRI_VERTS);
    m_cellToVertexRelation = CellVertexRelation(&m_elementSet, &m_vertexSet);
    m_cellToVertexRelation.bindIndices(static_cast<int>(m_cv_data.size()), &m_cv_data);

    // Delete old mesh
    m_surfaceMesh.reset();

    m_meshWasReindexed = true;
  }

  /// \brief Add the vertex positions and triangle boundary relations to the surface mesh
  void regenerateSurfaceMesh()
  {
    if(m_surfaceMesh != nullptr)
    {
      m_surfaceMesh.reset();
    }

    using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
    UMesh* triMesh = new UMesh(DIM, mint::TRIANGLE, m_vertexSet.size(), m_elementSet.size());

    // Add vertices to the mesh (i.e. vertex positions)
    for(int i = 0; i < m_vertexSet.size(); ++i)
    {
      const SpacePt& pt = vertexPosition(i);
      triMesh->appendNode(pt[0], pt[1], pt[2]);
    }

    // Add triangles to the mesh (i.e. boundary vertices)
    for(int i = 0; i < m_elementSet.size(); ++i)
    {
      const CellIndex* tv = &cellVertexIndices(i)[0];
      triMesh->appendCell(tv);
    }

    m_surfaceMesh.reset(triMesh);
  }
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_INOUT_OCTREE_MESH_WRAPPER__HPP_
