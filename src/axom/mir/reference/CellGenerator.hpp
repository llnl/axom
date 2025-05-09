// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file CellGenerator.hpp
 * 
 * \brief Contains the specification for the CellGenerator class.
 * 
 */

#ifndef __CELL_GENERATOR_H__
#define __CELL_GENERATOR_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"  // unified header for slam classes and functions

#include "axom/mir/reference/MIRMesh.hpp"
#include "axom/mir/reference/MIRUtilities.hpp"
#include "axom/mir/reference/MIRMeshTypes.hpp"
#include "axom/mir/reference/CellData.hpp"
#include "axom/mir/reference/ZooClippingTables.hpp"
#include "axom/mir/reference/MIRUtilities.hpp"
#include "axom/mir/reference/CellClipper.hpp"

//--------------------------------------------------------------------------------

namespace axom
{
namespace mir
{
//--------------------------------------------------------------------------------

/**
   * \class CellGenerator
   * 
   * \brief A class that generates that uses the clipping information to generate
   *        the data needed in order to produce a clean, reconstructed mesh.
   */
class CellGenerator
{
public:
  /**
       * \brief Default constructor.
       */
  CellGenerator();

  /**
       * \brief Default destructor.
       */
  ~CellGenerator();

  /**
       * \brief Generates the topology of the new elements resulting from a split.
       * 
       * \param newElements  An ordered map of the generated elements' IDs to a list of the vertex IDs in the local frame of the output element.
       * \param newVertices  An order map of the vertex IDs in the local frame of the output element to the generated elements' IDs it is associated with.
       * \param out_cellData  Container to store the topology data of the generated elements.
       */
  void generateTopologyData(const std::map<int, std::vector<int>>& newElements,
                            const std::map<int, std::vector<int>>& newVertices,
                            CellData& out_cellData);

  /**
       * \brief Generates the vertex positions values for each of the new vertices of the generated element.
       * 
       * \param  shapeType  The shape type of the element.
       * \param  newVertices  An ordered map of vertices that compose the newly generated elements.
       * \param  vertexPositions  A vector of positions for each of the original element's vertices.
       * \param  tValues  An array of t values where each of the midpoint vertices are for the newly generated elements.
       * \param  out_cellData  Container to store the vertex position data of the generated elements.
       * 
       * \note  New vertex positions are interpolated from the original vertex positions.
       */
  void generateVertexPositions(const mir::Shape shapeType,
                               const std::map<int, std::vector<int>>& newVertices,
                               const std::vector<mir::Point2>& vertexPositions,
                               axom::float64* tValues,
                               CellData& out_cellData);

  /**
       * \brief Generates the vertex volume fractions for each of the new vertices of the generated element.
       * 
       * \param  shapeType  The shape type of the element.
       * \param  newVertices  An ordered map of vertices that compose the newly generated elements.
       * \param  vertexVF  A vector of positions for each of the original element's vertices.
       * \param  tValues  An array of t values where each of the midpoint vertices are for the newly generated elements.
       * \param  out_cellData  Container to store the vertex position data of the generated elements.
       * 
       * \note  New vertex positions are interpolated from the original vertex volume fractions.
       */
  void generateVertexVolumeFractions(const mir::Shape shapeType,
                                     const std::map<int, std::vector<int>>& newVertices,
                                     const std::vector<std::vector<axom::float64>>& vertexVF,
                                     axom::float64* tValues,
                                     CellData& out_cellData);

  /**
       * \brief Determines the more dominant material of the two given for the given element.
       * 
       * \param shapeType  An enumerator denoting the element's shape.
       * \param vertexIDs  A list of vertex IDs into the vertexVF param.
       * \param matOne  The ID of the first material.
       * \param matTwo  The ID of the second material.
       * \param vertexVF  The list of volume fractions associated with the given vertices in the vertexIDs param.
       * 
       * \return  The ID of the dominant material of the element.
       * 
       * \note The dominant element for the 2D/3D cases will be the same as the material present at one of the 
       *       original vertices that existed prior to the split. So, if you can find this vertex and its dominant
       *       material, then you know the dominant material of this new element.
       * 
       * \note It is assumed that the given cell is one that results from splitting its parent cell.
       */
  int determineCleanCellMaterial(const Shape shapeType,
                                 const std::vector<int>& vertexIDs,
                                 const int matOne,
                                 const int matTwo,
                                 const std::vector<std::vector<axom::float64>>& vertexVF);
};

//--------------------------------------------------------------------------------

}  // namespace mir
}  // namespace axom

#endif
