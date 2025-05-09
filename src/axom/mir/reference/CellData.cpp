// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/reference/CellData.hpp"

namespace axom
{
namespace mir
{
//--------------------------------------------------------------------------------

CellData::CellData() : m_numVerts(0), m_numElems(0) { }

//--------------------------------------------------------------------------------

void CellData::mergeCell(const CellData& cellToMerge)
{
  // Initialize index offsets
  int evBeginsOffset = m_topology.m_evInds.size();
  int veBeginsOffset = m_topology.m_veInds.size();

  int vertexIndexOffset = m_numVerts;
  int elementIndexOffset = m_numElems;

  // Merge the cell topology information
  for(unsigned long i = 0; i < cellToMerge.m_topology.m_evInds.size(); ++i)
  {
    m_topology.m_evInds.push_back(cellToMerge.m_topology.m_evInds[i] + vertexIndexOffset);
  }

  for(unsigned long i = 1; i < cellToMerge.m_topology.m_evBegins.size(); ++i)
  {
    m_topology.m_evBegins.push_back(cellToMerge.m_topology.m_evBegins[i] + evBeginsOffset);
  }

  for(unsigned long i = 0; i < cellToMerge.m_topology.m_veInds.size(); ++i)
  {
    m_topology.m_veInds.push_back(cellToMerge.m_topology.m_veInds[i] + elementIndexOffset);
  }

  for(unsigned long i = 1; i < cellToMerge.m_topology.m_veBegins.size(); ++i)
  {
    m_topology.m_veBegins.push_back(cellToMerge.m_topology.m_veBegins[i] + veBeginsOffset);
  }

  // Merge the vertex positions
  for(unsigned long i = 0; i < cellToMerge.m_mapData.m_vertexPositions.size(); ++i)
  {
    m_mapData.m_vertexPositions.push_back(cellToMerge.m_mapData.m_vertexPositions[i]);
  }

  // Merge the vertex volume fractions
  for(unsigned long matID = 0; matID < m_mapData.m_vertexVolumeFractions.size(); ++matID)
  {
    for(unsigned long vID = 0; vID < cellToMerge.m_mapData.m_vertexVolumeFractions[matID].size();
        ++vID)
    {
      m_mapData.m_vertexVolumeFractions[matID].push_back(
        cellToMerge.m_mapData.m_vertexVolumeFractions[matID][vID]);
    }
  }

  // Merge the elements' dominant materials
  for(unsigned long i = 0; i < cellToMerge.m_mapData.m_elementDominantMaterials.size(); ++i)
  {
    m_mapData.m_elementDominantMaterials.push_back(
      cellToMerge.m_mapData.m_elementDominantMaterials[i]);
  }

  // Merge the elements' parent ids
  for(unsigned long i = 0; i < cellToMerge.m_mapData.m_elementParents.size(); ++i)
  {
    m_mapData.m_elementParents.push_back(cellToMerge.m_mapData.m_elementParents[i]);
  }

  // Merge the elements' shape types
  for(unsigned long i = 0; i < cellToMerge.m_mapData.m_shapeTypes.size(); ++i)
  {
    m_mapData.m_shapeTypes.push_back(cellToMerge.m_mapData.m_shapeTypes[i]);
  }

  // Merge the total number of verts and elems in the resulting cell
  m_numVerts += cellToMerge.m_numVerts;
  m_numElems += cellToMerge.m_numElems;
}

//--------------------------------------------------------------------------------

}  // namespace mir
}  // namespace axom
