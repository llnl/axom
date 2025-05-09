// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/reference/CellGenerator.hpp"

namespace axom
{
namespace mir
{
//--------------------------------------------------------------------------------

CellGenerator::CellGenerator() { }

//--------------------------------------------------------------------------------

CellGenerator::~CellGenerator() { }

//--------------------------------------------------------------------------------

void CellGenerator::generateTopologyData(const std::map<int, std::vector<int>>& newElements,
                                         const std::map<int, std::vector<int>>& newVertices,
                                         CellData& out_cellData)
{
  // Store the evInds and evBegins data in the output vectors
  int currentEVBeginIndex = 0;
  for(auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    // Push the start index of the next element
    out_cellData.m_topology.m_evBegins.push_back(currentEVBeginIndex);

    // Push the next element's vertices
    for(unsigned int vIndex = 0; vIndex < itr->second.size(); ++vIndex)
    {
      out_cellData.m_topology.m_evInds.push_back(itr->second[vIndex]);
      ++currentEVBeginIndex;
    }
  }

  // Push the index that occurs after the last vertex
  out_cellData.m_topology.m_evBegins.push_back(currentEVBeginIndex);

  // Store the veInds and veBegins data in the output vectors
  int currentVEBeginIndex = 0;
  for(auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    // Push the start index of the vertex's elements
    out_cellData.m_topology.m_veBegins.push_back(currentVEBeginIndex);

    // Push the next vertex's elements
    for(unsigned int eIndex = 0; eIndex < itr->second.size(); eIndex++)
    {
      out_cellData.m_topology.m_veInds.push_back(itr->second[eIndex]);
      ++currentVEBeginIndex;
    }
  }

  // Push the index that occurs after the last element
  out_cellData.m_topology.m_veBegins.push_back(currentVEBeginIndex);
}

//--------------------------------------------------------------------------------

void CellGenerator::generateVertexPositions(const mir::Shape shapeType,
                                            const std::map<int, std::vector<int>>& newVertices,
                                            const std::vector<mir::Point2>& vertexPositions,
                                            axom::float64* tValues,
                                            CellData& out_cellData)
{
  for(auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    int vID = itr->first;

    if(vID < mir::utilities::numVerts(shapeType))
    {
      // This vertex is one of the shape's original vertices
      out_cellData.m_mapData.m_vertexPositions.push_back(vertexPositions[vID]);
    }
    else if(mir::utilities::isCenterVertex(shapeType, vID))
    {
      // Average the vertex position values at the corners of the shape
      mir::Point2 centroid = mir::utilities::computeAveragePoint(vertexPositions);

      out_cellData.m_mapData.m_vertexPositions.push_back(centroid);
    }
    else
    {
      // This vertex is between two of the shape's original vertices
      int vIDFrom = mir::utilities::getEdgeEndpoint(shapeType, vID, true);
      int vIDTo = mir::utilities::getEdgeEndpoint(shapeType, vID, false);

      out_cellData.m_mapData.m_vertexPositions.push_back(
        mir::Point2::lerp(vertexPositions[vIDFrom], vertexPositions[vIDTo], tValues[vID]));
    }
  }
}

//--------------------------------------------------------------------------------

void CellGenerator::generateVertexVolumeFractions(const mir::Shape shapeType,
                                                  const std::map<int, std::vector<int>>& newVertices,
                                                  const std::vector<std::vector<axom::float64>>& vertexVF,
                                                  axom::float64* tValues,
                                                  CellData& out_cellData)
{
  out_cellData.m_mapData.m_vertexVolumeFractions.resize(
    vertexVF.size());  // vertexVF size is the number of materials

  for(auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    int vID = itr->first;

    for(unsigned long matID = 0; matID < vertexVF.size(); ++matID)
    {
      if(vID < mir::utilities::numVerts(shapeType))
      {
        // This vertex is one of the shape's original vertices
        out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(vertexVF[matID][vID]);
      }
      else if(mir::utilities::isCenterVertex(shapeType, vID))
      {
        // Average the vertex volume fractions values at the corners of the shape
        axom::float64 averageValue = mir::utilities::computeAverageFloat(vertexVF[matID]);

        out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(averageValue);
      }
      else
      {
        // This vertex is between two of the shape's original vertices
        int vIDFrom = mir::utilities::getEdgeEndpoint(shapeType, vID, true);
        int vIDTo = mir::utilities::getEdgeEndpoint(shapeType, vID, false);

        out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(
          axom::utilities::lerp(vertexVF[matID][vIDFrom], vertexVF[matID][vIDTo], tValues[vID]));
      }
    }
  }
}

//--------------------------------------------------------------------------------

int CellGenerator::determineCleanCellMaterial(const Shape elementShape,
                                              const std::vector<int>& vertexIDs,
                                              const int matOne,
                                              const int matTwo,
                                              const std::vector<std::vector<axom::float64>>& vertexVF)
{
  int dominantMaterial = matOne;

  axom::float64 matOneVF = -1.0;
  axom::float64 matTwoVF = -1.0;

  for(unsigned long it = 0; it < vertexIDs.size(); ++it)
  {
    int vID = vertexIDs[it];

    if(vID >= 0 && vID < mir::utilities::numVerts(elementShape))
    {
      if(matOne != NULL_MAT)
      {
        matOneVF = vertexVF[matOne][vID];
      }
      if(matTwo != NULL_MAT)
      {
        matTwoVF = vertexVF[matTwo][vID];
      }

      dominantMaterial = (matOneVF > matTwoVF) ? matOne : matTwo;
    }
  }

  return dominantMaterial;
}

//--------------------------------------------------------------------------------

}  // namespace mir
}  // namespace axom
