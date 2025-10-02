// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file CellData.hpp
 * 
 * \brief Contains the specifications for the CellData class 
 *        and CellTopologyData and CellMapData structs.
 */

#ifndef __CELL_DATA_H__
#define __CELL_DATA_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"

#include "axom/mir/reference/MIRMeshTypes.hpp"

namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{
/**
   * \struct CellTopologyData
   * 
   * \brief Struct for collecting data that specifies a mesh or cell's connectivity/topology.
   */
struct CellTopologyData
{
  std::vector<PosType> m_evInds;
  std::vector<PosType> m_evBegins;
  std::vector<PosType> m_veInds;
  std::vector<PosType> m_veBegins;
};

/**
   * \struct CellMapData
   * 
   * \brief Struct for collecting the data that will populate a mesh's map data structures.
   */
struct CellMapData
{
  std::vector<mir::Point2> m_vertexPositions;  // Data that goes into MIRMesh's vertexPositions PointMap
  std::vector<std::vector<axom::float64>>
    m_vertexVolumeFractions;  // Data that goes into MIRMesh's materialVolumeFractionsVertex ScalarMap
  std::vector<int> m_elementDominantMaterials;  // Data that goes into MIRMesh's elementDominantColors IntMap
  std::vector<int> m_elementParents;     // Data that goes into MIRMesh's elementParentIDs IntMap
  std::vector<mir::Shape> m_shapeTypes;  // Data that goes into MIRMesh's shapeType IntMap
};

/**
   * \class CellData
   * 
   * \brief The CellData class represents an arbitrary number of cells that are 
   *        within the same local coordinate system (i.e. share a set of vertices 
   *        and elements).
   * 
   * \detail This class is intended to be used as a helper class to hold intermediate 
   *         data while processing a mesh, and to be used as input to the MIRMesh class
   *         to fully initialize it.
   */
class CellData
{
public:
  /**
      * \brief Default constructor.
      */
  CellData();

  /**
       * \brief Default destructor.
       */
  ~CellData() = default;

  /**
       * \brief Merges the cell data from the given cell into this cell.
       * 
       * \param cellToMerge  The cell whose data will be taken and merged in.
       */
  void mergeCell(const CellData& cellToMerge);

public:
  int m_numVerts;
  int m_numElems;

  CellTopologyData m_topology;
  CellMapData m_mapData;
};
}  // namespace mir
}  // namespace axom

#endif
