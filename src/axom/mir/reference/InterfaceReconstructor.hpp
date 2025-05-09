// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InterfaceReconstructor.hpp
 * 
 * \brief Contains the specification for the InterfaceReconstructor class.
 * 
 */

#ifndef __INTERFACE_RECONSTRUCTOR_H__
#define __INTERFACE_RECONSTRUCTOR_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"

#include "axom/mir/reference/MIRMesh.hpp"
#include "axom/mir/reference/CellData.hpp"
#include "axom/mir/reference/ZooClippingTables.hpp"
#include "axom/mir/reference/MIRUtilities.hpp"
#include "axom/mir/reference/CellClipper.hpp"
#include "axom/mir/reference/CellGenerator.hpp"

#include <map>

namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{
/**
   * \class InterfaceReconstructor
   * 
   * \brief A class that contains the functionality for taking an input mesh
   *        with mixed cells and outputs a mesh with clean cells.
   * 
   * \detail This class requires that the user create a mesh of type MIRMesh 
   *         and pass that into one of the reconstruction methods. There are 
   *         currently two reconstructions methods implemented: one is the 
   *         a zoo-based algorithm described in Meredith 2004, and the other 
   *         is the iterative version described in Meredith and Childs 2010.
   */
class InterfaceReconstructor
{
public:
  /**
       * \brief Default constructor.
       */
  InterfaceReconstructor();

  /**
       * \brief Default destructor.
       */
  ~InterfaceReconstructor();

  /**
       * \brief  Performs material interface reconstruction using the zoo-based algorithm.
       * 
       * \param inputMesh  The mesh composed of mixed cells.
       * \param outputMesh The mesh composed of clean cells. 
       */
  void computeReconstructedInterface(mir::MIRMesh& inputMesh, mir::MIRMesh& outputMesh);

  /**
       * \brief Performs material interface reconstruction using an iterative version of the zoo-based algorithm.
       * 
       * \param inputMesh  The mesh made up of mixed cells.
       * \param numIterations  The number of iterations for which to run the algorithm.
       * \param percent  The percent of the difference to use when modifying the original element volume fractions.
       * \param outputMesh  The mesh composed of clean cells.
       */
  void computeReconstructedInterfaceIterative(mir::MIRMesh& inputMesh,
                                              const int numIterations,
                                              const axom::float64 percent,
                                              mir::MIRMesh& outputMesh);

  /**
       * \brief Generates a set of clean cells by splitting the element with the two given materials.
       * 
       * \param  shapeType  The shape type of the element to be split.
       * \param  parentElementID  The ID parent element.
       * \param  matOne  The first material to split with.
       * \param  matTwo  The second material to split with.
       * \param  elementVertices  The vertices of the element to be split.
       * \param  originalElementVertexVF  The original vertex volume fractions associated with the vertices of the element to be split.
       * \param  originalElementVertexPositions  The original vertex positions associated with the vertices of the element to be split.
       * \param out_cellData  Container to store the data of the generated elements.
       */
  void generateCleanCells(mir::Shape shapeType,
                          const int parentElementID,
                          const int matOne,
                          const int matTwo,
                          const std::vector<int>& elementVertices,
                          const std::vector<std::vector<axom::float64>>& originalElementVertexVF,
                          const std::vector<mir::Point2>& originalElementVertexPositions,
                          CellData& out_cellData);

private:
  mir::MIRMesh m_originalMesh;
};
}  // namespace mir
}  // namespace axom
#endif
