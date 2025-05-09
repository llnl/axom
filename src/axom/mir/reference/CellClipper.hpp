// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file CellClipper.hpp
 * 
 * \brief Contains the specification for the CellClipper class.
 * 
 */

#ifndef __CELL_CLIPPER_H
#define __CELL_CLIPPER_H

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"  // unified header for slam classes and functions

#include "axom/mir/reference/MIRMesh.hpp"
#include "axom/mir/reference/MIRUtilities.hpp"
#include "axom/mir/reference/MIRMeshTypes.hpp"
#include "axom/mir/reference/CellData.hpp"
#include "axom/mir/reference/ZooClippingTables.hpp"

//--------------------------------------------------------------------------------

namespace axom
{
namespace mir
{
//--------------------------------------------------------------------------------

/**
   * \class CellClipper
   * 
   * \brief A class that contains the functionality for taking an input cell
   *        and determining how it should be clipped.
   * 
   */
class CellClipper
{
public:
  /**
       * \brief Default constructor.
       */
  CellClipper();

  /**
       * \brief Default destructor.
       */
  ~CellClipper();

  /**
       * \brief Computes the elements, vertices, and t values resulting from splitting the element with the given vertex volume fractions.
       * 
       * \param shapeType  The shape type of the element.
       * \param vertexVF  The vertex volume fractions of the two materials with which to clip the cell.
       * \param newElements  An ordered map of the generated elements' IDs to a list of the vertex IDs in the local frame of the output element.
       * \param newVertices  An order map of the vertex IDs in the local frame of the output element to the generated elements' IDs it is associated with.
       * \param  tValues  An array of t values where each of the midpoint vertices are for the newly generated elements.
       * 
       */
  void computeClippingPoints(const mir::Shape shapeType,
                             const std::vector<std::vector<axom::float64>>& vertexVF,
                             std::map<int, std::vector<int>>& newElements,
                             std::map<int, std::vector<int>>& newVertices,
                             axom::float64* tValues);

  /**
       * \brief  Determines the index into the clipping table of the given shape based on the volume fractions at its vertices.
       * 
       * \param shapeType  The shape type of the element.
       * \param matOneVF  The volume fractions at the vertices of the element for the first material.
       * \param matTwoVF  The volume fractions at the vertices of the element for the second material.
       * 
       * \return The index into the clipping table.
       */
  unsigned int determineClippingCase(const mir::Shape shapeType,
                                     const std::vector<axom::float64>& matOneVF,
                                     const std::vector<axom::float64>& matTwoVF);

  /**
       * \brief Computes the t value as a percent from vertex one to vertex two based on the materials given.
       * 
       * \param vfMatOneVertexOne  The volume fraction of material one present at vertex one.
       * \param vfMatTwoVertexOne  The volume fraction of material two present at vertex one.
       * \param vfMatOneVertexTwo  The volume fraction of material one present at vertex two.
       * \param vfMatTwoVertexTwo  The volume fraction of material two present at vertex two.
       * 
       * \return The percent of the distance from vertex one to vertex two where the edge should be clipped.
       * 
       * \note  When one material's volume fractions dominates the other material's, then the edge should not be clipped and this function will return 0.0.
       * \note  When one of the materials is the NULL_MAT (and has vf = -1.0), these values are set to 0 in order to interpolate properly.
       */
  axom::float64 computeTValueOnEdge(axom::float64 vfMatOneVertexOne,
                                    axom::float64 vfMatTwoVertexOne,
                                    axom::float64 vfMatOneVertexTwo,
                                    axom::float64 vfMatTwoVertexTwo);

private:
  /**
      * \brief Returns a reference to the appropriate clipping table to use for the shape type.
      * 
      * \param The shape type of the element.
      * 
      * \return A reference to the clipping table. 
      */
  const std::vector<std::vector<int>>& getClipTable(const mir::Shape shapeType);
};

//--------------------------------------------------------------------------------

}  // namespace mir
}  // namespace axom

#endif
