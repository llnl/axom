// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/reference/CellClipper.hpp"

namespace axom
{
namespace mir
{
//--------------------------------------------------------------------------------

CellClipper::CellClipper() { }

//--------------------------------------------------------------------------------

CellClipper::~CellClipper() { }

//--------------------------------------------------------------------------------

// Computes the t-values where each edge is clipped, as well as the topology of the new output cells after clipping the original cell
// Outputs the newElements, newVertices maps and the verticesClippingTValue array[]
void CellClipper::computeClippingPoints(const mir::Shape shapeType,
                                        const std::vector<std::vector<axom::float64>>& vertexVF,
                                        std::map<int, std::vector<int>>& newElements,
                                        std::map<int, std::vector<int>>& newVertices,
                                        axom::float64* tValues)
{
  // Determine the clipping case for the current element
  unsigned int caseIndex = determineClippingCase(shapeType, vertexVF[0], vertexVF[1]);

  std::vector<std::vector<int>> clipTable = getClipTable(shapeType);

  // Create the new polygons based on the clipping case
  int currentElementIndex = 0;  // the next available element index
  int i = 0;
  int numVertices = clipTable[caseIndex][i];

  // for each new element in the current clipping case
  while(numVertices != -1)
  {
    // for each vertex of the new element
    for(int j = 0; j < numVertices; ++j)
    {
      // Find the id of the next vertex of the new element
      int vID = clipTable[caseIndex][i + (j + 1)];

      // Associate the vertex and element together
      newElements[currentElementIndex].push_back(vID);
      newVertices[vID].push_back(currentElementIndex);

      if(vID >= mir::utilities::numVerts(shapeType) &&
         vID < (mir::utilities::maxPossibleNumVerts(shapeType) - 1))
      {
        int vertexOneID = mir::utilities::getEdgeEndpoint(shapeType, vID, true);
        int vertexTwoID = mir::utilities::getEdgeEndpoint(shapeType, vID, false);

        tValues[vID] = computeTValueOnEdge(vertexVF[0][vertexOneID],
                                           vertexVF[1][vertexOneID],
                                           vertexVF[0][vertexTwoID],
                                           vertexVF[1][vertexTwoID]);
      }
    }

    // Increment the element index counter, marking the current element as being finished processed
    currentElementIndex++;

    // Increase index into lookup table to the next element
    i += (numVertices + 1);
    numVertices = clipTable[caseIndex][i];
  }
}

//--------------------------------------------------------------------------------

unsigned int CellClipper::determineClippingCase(const mir::Shape shapeType,
                                                const std::vector<axom::float64>& matOneVF,
                                                const std::vector<axom::float64>& matTwoVF)
{
  unsigned int caseIndex = 0;

  int numVertices = mir::utilities::numVerts(shapeType);

  for(int vID = 0; vID < numVertices; ++vID)
  {
    if(matOneVF[vID] > matTwoVF[vID])
    {
      unsigned int bitIndex = (numVertices - 1) - vID;

      caseIndex |= (1 << bitIndex);
    }
  }

  return caseIndex;
}

//--------------------------------------------------------------------------------

axom::float64 CellClipper::computeTValueOnEdge(axom::float64 vfMatOneVertexOne,
                                               axom::float64 vfMatTwoVertexOne,
                                               axom::float64 vfMatOneVertexTwo,
                                               axom::float64 vfMatTwoVertexTwo)
{
  axom::float64 ret = 0.0;

  // TODO: Perhaps just handle NULL_MAT by return 0, since that is what will happen anyways?
  // Handle NULL_MAT, which has a vf of -1.0, but which needs to be 0.0 for the purposes of computing the clipping point
  if(vfMatOneVertexOne < 0.0)
  {
    vfMatOneVertexOne = 0.0;
  };
  if(vfMatTwoVertexOne < 0.0)
  {
    vfMatTwoVertexOne = 0.0;
  };
  if(vfMatOneVertexTwo < 0.0)
  {
    vfMatOneVertexTwo = 0.0;
  };
  if(vfMatTwoVertexTwo < 0.0)
  {
    vfMatTwoVertexTwo = 0.0;
  };

  axom::float64 numerator = vfMatTwoVertexOne - vfMatOneVertexOne;
  axom::float64 denominator =
    -vfMatOneVertexOne + vfMatOneVertexTwo + vfMatTwoVertexOne - vfMatTwoVertexTwo;

  if(denominator != 0.0)
  {
    ret = numerator / denominator;
  }

  if(ret > 1.0 || ret < 0.0)
  {
    // This shouldn't happen...
    printf("    OUT OF BOUNDS T VALUE: %f\n", ret);

    // Clamp the t value
    ret = fmin(1.0, ret);
    ret = fmax(0.0, ret);
  }

  return ret;
}

//--------------------------------------------------------------------------------

const std::vector<std::vector<int>>& CellClipper::getClipTable(const mir::Shape shapeType)
{
  switch(shapeType)
  {
  case mir::Shape::Triangle:
    return triangleClipTableVec;
  case mir::Shape::Quad:
    return quadClipTableVec;
  case mir::Shape::Tetrahedron:
    return tetrahedronClipTableVec;
  case mir::Shape::Pyramid:
    return pyramidClipTableVec;
  case mir::Shape::Triangular_Prism:
    return triangularPrismClipTableVec;
  case mir::Shape::Hexahedron:
    return hexahedronClipTableVec;
  default:
    printf("No clipping table for this shape type.\n");
    return triangleClipTableVec;
  }
}

//--------------------------------------------------------------------------------

}  // namespace mir
}  // namespace axom
