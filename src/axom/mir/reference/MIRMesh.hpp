// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MIRMesh.hpp
 * 
 * \brief Contains the specification for the MIRMesh class.
 * 
 */

#ifndef __MIR_MESH_H__
#define __MIR_MESH_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"  // unified header for slam classes and functions

#include "axom/mir/reference/MIRMeshTypes.hpp"
#include "axom/mir/reference/CellData.hpp"

// C/C++ includes
#include <cmath>  // for definition of M_PI, exp()
#include <vector>
#include <string>
#include <fstream>

namespace numerics = axom::numerics;
namespace slam = axom::slam;

//--------------------------------------------------------------------------------
namespace axom
{
namespace mir
{
const int NULL_MAT = -1;

//--------------------------------------------------------------------------------

/**
   * \class MIRMesh
   * 
   * \brief The MIRMesh class represents a finite element mesh containing element volume fractions.
   * 
   * \detail This class is meant to be used in conjunction with the InterfaceReconstructor class
   *         to process the mesh.
   * 
   */
class MIRMesh
{
public:
  /**
       * \brief Default constructor.
       */
  MIRMesh();

  /**
       * \brief Copy constructor.
       */
  MIRMesh(const MIRMesh& other);

  /**
       * \brief Default destructor.
       */
  ~MIRMesh() = default;

  /**
       * \brief Copy assignment.
       */
  MIRMesh& operator=(const MIRMesh& other);

  /**
       * \brief Initializes a mesh with the provided data and topology.
       * 
       * \param _verts  The set of vertices of the mesh.
       * \param _elems  The set of elements of the mesh.
       * \param _numMaterials  The number of materials present in the mesh.
       * \param _topology  The topology/connectivity of the mesh.
       * \param _mapData  The data used to initialized the maps associated with the vertex and element sets.
       * \param _elementVF  The volume fractions of each element. Note that this is an optional parameter.
       */
  void initializeMesh(const VertSet _verts,
                      const ElemSet _elems,
                      const int _numMaterials,
                      const CellTopologyData& _topology,
                      const CellMapData& _mapData,
                      const std::vector<std::vector<axom::float64>>& _elementVF = {});

  /**
       * \brief Constructs the mesh boundary and coboundary relations.
       * 
       * \note The boundary relation is from each element to its vertices.
       * \note The coboundary relation is from each vertex to its elements.
       */
  void constructMeshRelations();

  /**
       * \brief Constructs the element and vertex volume fraction maps.
       * 
       * \param elementVF  The volume fractions of each element.
       */
  void constructMeshVolumeFractionsMaps(const std::vector<std::vector<axom::float64>>& elementVF);

  /**
       * \brief Constructs the vertex volume fraction map.
       * 
       * \param vertexVF  The volume fractions of each vertex.
       */
  void constructMeshVolumeFractionsVertex(const std::vector<std::vector<axom::float64>>& vertexVF);

  /**
       * \brief Prints out the data contained within this mesh in a nice format.
       */
  void print();

  /**
       * \brief Reads in a mesh specified by the given file.
       * 
       * \param filename  The location where the mesh file will be read from.
       */
  void readMeshFromFile(const std::string filename);

  /**
       * \brief Writes out the mesh to the given file.
       * 
       * \param filename  The location where the mesh file will be written.
       * 
       * \note Currently reads in an ASCII, UNSTRUCTURED_GRID .vtk file.
       */
  void writeMeshToFile(const std::string& dirName,
                       const std::string& fileName,
                       const std::string& separator = "/");

  /**
       * \brief  Computes the volume fractions of the elements of the original mesh.
       */
  std::vector<std::vector<axom::float64>> computeOriginalElementVolumeFractions();

  /**
       * \brief Checks that this instance of MIRMesh is valid
       *
       * \return True if this instance is valid, false otherwise
       */
  bool isValid(bool verbose) const;

private:
  /// Utility predicate to check if the mesh's topology is valid
  bool isTopologyValid(bool verbose) const;
  /// Utility predicate to check if the mesh's volume fractions are valid
  bool areVolumeFractionsValid(bool verbose) const;
  /// Utility predicate to check if the mesh's maps are valid
  bool areMapsValid(bool verbose) const;

  /**
       * \brief Constructs the positions map on the vertices.
       * 
       * \param data  The vector of position data for each vertex.
       */
  void constructVertexPositionMap(const std::vector<Point2>& data);

  /**
       * \brief Constructs the map of elements to their original element parent.
       * 
       * \param cellParents  The vector of parent IDs for each element of the mesh.
       */
  void constructElementParentMap(const std::vector<int>& elementParents);

  /**
       * \brief Constructs the map of elements to their dominant materials.
       * 
       * \param dominantMaterials  A vector of material ids that are the dominant material of each element.
       */
  void constructElementDominantMaterialMap(const std::vector<int>& dominantMaterials);

  /**
       * \brief Constructs the map of elements to their shape types.
       * 
       * \param shapeTypes  A vector of shape enumerators that are the shape type of each element.
       */
  void constructElementShapeTypesMap(const std::vector<mir::Shape>& shapeTypes);

  /**
       * \brief Computes the area of the triangle defined by the given three vertex positions using Heron's formula.
       * 
       * \param p0  The position of the first vertex.
       * \param p1  The position of the second vertex.
       * \param p2  The position of the third vertex.
       */
  axom::float64 computeTriangleArea(Point2 p0, Point2 p1, Point2 p2);

  /**
       * \brief  Computes the area of the quad defined by the given four vertex positions.
       * 
       * \param p0  The position of the first vertex.
       * \param p1  The position of the second vertex.
       * \param p2  The position of the third vertex.
       * \param p3  The position of the fourth vertex.
       * 
       * \note It is assumed that the points are given in consecutive, counter-clockwise order.
       */
  axom::float64 computeQuadArea(Point2 p0, Point2 p1, Point2 p2, Point2 p3);

  /****************************************************************
     *                        VARIABLES
     ****************************************************************/
public:
  // Mesh Set Definitions
  VertSet m_verts;  // the set of vertices in the mesh
  ElemSet m_elems;  // the set of elements in the mesh

  // Mesh Relation Definitions
  ElemToVertRelation m_bdry;    // Boundary relation from elements to vertices
  VertToElemRelation m_cobdry;  // Coboundary relation from vertices to elements

  // Mesh Map Definitions
  PointMap m_vertexPositions;                               // vertex position for each vertex
  std::vector<ScalarMap> m_materialVolumeFractionsElement;  // the volume fractions of each material for each element
  std::vector<ScalarMap> m_materialVolumeFractionsVertex;  // the volume fractions of each material for each vertex
  IntMap m_elementParentIDs;  // the ID of the parent element from the original mesh
  IntMap m_elementDominantMaterials;  // the dominant material of the cell (a processed mesh should have only one material per cell)
  IntMap m_shapeTypes;  // the int enumerator of what type of shape each element is

  int m_numMaterials;               // the number of materials present in the mesh
  CellTopologyData m_meshTopology;  // the topology/connectivity of the mesh
};

//--------------------------------------------------------------------------------
}  // namespace mir
}  // namespace axom
#endif
