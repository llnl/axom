// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MIRUtilities.hpp
 * 
 * \brief Contains a set of functions that provide general utility
 *        within Axom's MIR component.
 * 
 */

#ifndef __MIR_UTILITIES_HPP__
#define __MIR_UTILITIES_HPP__

#include "axom/mir/reference/ZooClippingTables.hpp"

//--------------------------------------------------------------------------------

namespace axom
{
namespace mir
{
namespace utilities
{
//--------------------------------------------------------------------------------

/**
   * \brief Determines the number of vertices of the given shape in the finite element zoo.
   * 
   * \param shape  The shape type from the finite element zoo.
   * 
   * \return  THe number of vertices of the shape.
   */
inline int numVerts(mir::Shape shape)
{
  int numVertices = 0;
  switch(shape)
  {
  case mir::Shape::Triangle:
    numVertices = 3;
    break;
  case mir::Shape::Quad:
    numVertices = 4;
    break;
  case mir::Shape::Tetrahedron:
    numVertices = 4;
    break;
  case mir::Shape::Pyramid:
    numVertices = 5;
    break;
  case mir::Shape::Triangular_Prism:
    numVertices = 6;
    break;
  case mir::Shape::Hexahedron:
    numVertices = 8;
    break;
  default:
    printf("Invalid shape. Cannot determine numVerts().\n");
  }
  return numVertices;
}

//--------------------------------------------------------------------------------

/**
   * \brief Determines the maximum number of possible vertices of the given shape in the finite element zoo.
   *        This number includes the midpoint vertices between each of the original shape's vertices.
   * 
   * \param shape  The shape type from the finite element zoo.
   * 
   * \return  THe number of vertices of the shape.
   */
inline int maxPossibleNumVerts(mir::Shape shape)
{
  int numVertices = -1;
  switch(shape)
  {
  case mir::Shape::Triangle:
    numVertices = 6 + 1;  // add one for the central vertex (not used)
    break;
  case mir::Shape::Quad:
    numVertices = 8 + 1;  // add one for the central vertex (not used)
    break;
  case mir::Shape::Tetrahedron:
    numVertices = 10 + 1;  // add one for the central vertex (not used)
    break;
  case mir::Shape::Pyramid:
    numVertices = 13 + 1;  // add one for the central vertex
    break;
  case mir::Shape::Triangular_Prism:
    numVertices = 15 + 1;  // add one for the central vertex
    break;
  case mir::Shape::Hexahedron:
    numVertices = 20 + 1;  // add one for the central vertex
    break;
  default:
    printf("Invalid shape. Cannot determine maxPossibleNumVerts().\n");
  }
  return numVertices;
}

//--------------------------------------------------------------------------------

/**
   * \brief  Returns the local vertex ID of the from/to vertex that is one of 
   *         the two endpoints that the edge the given midpoint is on.
   * 
   * \param shapeType  The shape type from the finite element zoo.
   * \param midpointVertexID  The ID of the vertex between the two endpoints.
   * \param isFromVertex  A flag denoting which of the two edge endpoints to return.
   * 
   * \return The vertex ID of one of the endpoints.
   */
inline int getEdgeEndpoint(const mir::Shape shapeType,
                           const int midpointVertexID,
                           const bool isFromVertex)
{
  switch(shapeType)
  {
  case mir::Shape::Triangle:
    if(midpointVertexID == 3 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 3 && !isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 4 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 4 && !isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 5 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 5 && !isFromVertex)
    {
      return 0;
    }
    break;
  case mir::Shape::Quad:
    if(midpointVertexID == 4 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 4 && !isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 5 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 5 && !isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 6 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 6 && !isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 7 && isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 7 && !isFromVertex)
    {
      return 0;
    }
    break;
  case mir::Shape::Tetrahedron:
    if(midpointVertexID == 4 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 4 && !isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 5 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 5 && !isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 6 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 6 && !isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 7 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 7 && !isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 8 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 8 && !isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 9 && isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 9 && !isFromVertex)
    {
      return 1;
    }
    break;
  case mir::Shape::Pyramid:
    if(midpointVertexID == 5 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 5 && !isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 6 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 6 && !isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 7 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 7 && !isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 8 && isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 8 && !isFromVertex)
    {
      return 0;
    }

    if(midpointVertexID == 9 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 9 && !isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 10 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 10 && !isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 11 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 11 && !isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 12 && isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 12 && !isFromVertex)
    {
      return 4;
    }
    break;
  case mir::Shape::Triangular_Prism:
    if(midpointVertexID == 6 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 6 && !isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 7 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 7 && !isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 8 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 8 && !isFromVertex)
    {
      return 0;
    }

    if(midpointVertexID == 9 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 9 && !isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 10 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 10 && !isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 11 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 11 && !isFromVertex)
    {
      return 5;
    }

    if(midpointVertexID == 12 && isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 12 && !isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 13 && isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 13 && !isFromVertex)
    {
      return 5;
    }
    if(midpointVertexID == 14 && isFromVertex)
    {
      return 5;
    }
    if(midpointVertexID == 14 && !isFromVertex)
    {
      return 3;
    }
    break;
  case mir::Shape::Hexahedron:
    if(midpointVertexID == 8 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 8 && !isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 9 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 9 && !isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 10 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 10 && !isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 11 && isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 11 && !isFromVertex)
    {
      return 0;
    }

    if(midpointVertexID == 12 && isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 12 && !isFromVertex)
    {
      return 5;
    }
    if(midpointVertexID == 13 && isFromVertex)
    {
      return 5;
    }
    if(midpointVertexID == 13 && !isFromVertex)
    {
      return 6;
    }
    if(midpointVertexID == 14 && isFromVertex)
    {
      return 6;
    }
    if(midpointVertexID == 14 && !isFromVertex)
    {
      return 7;
    }
    if(midpointVertexID == 15 && isFromVertex)
    {
      return 7;
    }
    if(midpointVertexID == 15 && !isFromVertex)
    {
      return 4;
    }

    if(midpointVertexID == 16 && isFromVertex)
    {
      return 0;
    }
    if(midpointVertexID == 16 && !isFromVertex)
    {
      return 4;
    }
    if(midpointVertexID == 17 && isFromVertex)
    {
      return 1;
    }
    if(midpointVertexID == 17 && !isFromVertex)
    {
      return 5;
    }
    if(midpointVertexID == 18 && isFromVertex)
    {
      return 2;
    }
    if(midpointVertexID == 18 && !isFromVertex)
    {
      return 6;
    }
    if(midpointVertexID == 19 && isFromVertex)
    {
      return 3;
    }
    if(midpointVertexID == 19 && !isFromVertex)
    {
      return 7;
    }
    break;
  default:
    printf("Edge endpoint case not implemented.\n");
    return -1;
    break;
  }
  return -1;
}

//--------------------------------------------------------------------------------

/**
 * \brief Determines if the shape type is a three dimensional shape or not.
 * 
 * \param shapeType  The shape type from the finite element zoo.
 * 
 * \return True, if the shape is a tetrahedron, pyramid, triangular prism, or a hexahedron.
 */
inline bool isShapeThreeDimensional(const mir::Shape shapeType)
{
  return (shapeType == mir::Tetrahedron || shapeType == mir::Pyramid ||
          shapeType == mir::Triangular_Prism || shapeType == mir::Hexahedron);
}

//--------------------------------------------------------------------------------

/**
 * \brief  Get the local vertex ID of the central vertex used when decomposing a 3D shape.
 * 
 * \param shapeType  The shape type from the finite element zoo.
 * 
 * \return  The local vertex ID of the central vertex of a 3D shape.
 * 
 * \note 2D shapes do not have a central vertex because they do not decompose.
 */
inline int getCenterVertex(const mir::Shape shapeType)
{
  switch(shapeType)
  {
  case mir::Tetrahedron:
    return 10;
  case mir::Pyramid:
    return 13;
  case mir::Triangular_Prism:
    return 15;
  case mir::Hexahedron:
    return 20;
  default:
    return -1;
  }
}

//--------------------------------------------------------------------------------

/**
 * \brief Determines if the given vertex is the central vertex of a 3D shape.
 * 
 * \param shapeType  The shape type from the finite element zoo.
 * \param vID  The local vertex ID.
 * 
 * \return True, if vertex ID is that of the central vertex of a 3D shape.
 * 
 * \note 2D shapes do not have a central vertex because they do not decompose.
 */
inline bool isCenterVertex(const mir::Shape shapeType, const int vID)
{
  if(shapeType == mir::Tetrahedron && vID == 10)
    return true;
  else if(shapeType == mir::Pyramid && vID == 13)
    return true;
  else if(shapeType == mir::Triangular_Prism && vID == 15)
    return true;
  else if(shapeType == mir::Hexahedron && vID == 20)
    return true;
  else
    return false;
}

//--------------------------------------------------------------------------------

/**
 * \brief Computes the average value of the float values given.
 * 
 * \param values  A vector of float values.
 * 
 * \return  The average value.
 */
inline axom::float64 computeAverageFloat(const std::vector<axom::float64>& values)
{
  axom::float64 sum = 0.0;
  for(unsigned long i = 0; i < values.size(); ++i)
  {
    sum += values[i];
  }
  return sum / (axom::float64)values.size();
}

//--------------------------------------------------------------------------------

/**
 * \brief  Computes the centroid for the points given.
 * 
 * \param  A vector of points.
 * 
 * \return  The centroid point.
 */
inline mir::Point2 computeAveragePoint(const std::vector<mir::Point2>& points)
{
  mir::Point2 centroid;
  for(auto i = 0u; i < points.size(); ++i)
  {
    centroid.array() += points[i].array();
  }
  centroid.array() /= (axom::float64)points.size();

  return centroid;
}

//--------------------------------------------------------------------------------

/**
 * \brief Determine the shape type of an element.
 * 
 * \param parentShapeType  The shape of the element from which the new element is generated.
 * \param numVerts  The number of vertices of the new element.
 * 
 * \note It is assumed that the given cell is one that results from splitting its parent cell.
 */
inline mir::Shape determineElementShapeType(const Shape parentShapeType, const int numVerts)
{
  mir::Shape newShapeType;
  if(parentShapeType == mir::Shape::Triangle || parentShapeType == mir::Shape::Quad)
  {
    // Handle the two-dimensional case
    switch(numVerts)
    {
    case 3:
      newShapeType = mir::Shape::Triangle;
      break;
    case 4:
      newShapeType = mir::Shape::Quad;
      break;
    default:
      newShapeType = mir::Shape::Triangle;
      printf(
        "2D Case: Invalid number of vertices in "
        "determineElementShapeType().\n");
      break;
    }
  }
  else
  {
    // Handle the three-dimensional case
    switch(numVerts)
    {
    case 4:
      newShapeType = mir::Shape::Tetrahedron;
      break;
    case 5:
      newShapeType = mir::Shape::Pyramid;
      break;
    case 6:
      newShapeType = mir::Shape::Triangular_Prism;
      break;
    case 8:
      newShapeType = mir::Shape::Hexahedron;
      break;
    default:
      newShapeType = mir::Shape::Tetrahedron;
      printf(
        "3D Case: Invalid number of vertices in "
        "determineElementShapeType().\n");
      break;
    }
  }
  return newShapeType;
}

//--------------------------------------------------------------------------------

}  // namespace utilities
}  // namespace mir
}  // namespace axom

#endif
