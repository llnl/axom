// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MeshTester.hpp
 * 
 * \brief Contains the specification for the MeshTester class.
 * 
 */

#ifndef __AXOM_BUMP_MESH_TESTER_HPP__
#define __AXOM_BUMP_MESH_TESTER_HPP__

#include "axom/core.hpp"  // for axom macros
#include "axom/primal.hpp"

#include <algorithm>
#include <vector>

#include <conduit.hpp>

namespace axom
{
namespace bump
{
namespace data
{
/*!
 * \class MeshTester
 * 
 * \brief A class used to generate meshes with specific properties so
 *        that the reconstruction output can be validated visually.
 * 
 */
class MeshTester
{
public:
  template <typename T>
  using Vec = std::vector<T>;

  using IndexVec = Vec<int>;
  using VolFracVec = Vec<axom::float64>;
  using VolumeFractions = Vec<VolFracVec>;
  using Point2 = axom::primal::Point<float, 2>;

public:
  /*!
   * \brief Default constructor. 
   */
  MeshTester() = default;

  /*!
   * \brief Default destructor.
   */
  ~MeshTester() = default;

  /*!
   * \brief Make a structured mesh instead of an unstructured mesh.
   *
   * \param structured True for structured; false otherwise.
   */
  void setStructured(bool structured);

public:
  /*!
   * \brief Initializes an MIRMesh based on the example from Meredith 2004 paper.
   * 
   * \note The mesh is a 3x3 uniform grid of quads with 2 materials.
   * 
   * \return  The generated mesh.
   */
  void initTestCaseOne(conduit::Node &mesh);

  /*!
   * \brief Initializes an MIRMesh based on the example from Meredith and Childs 2010 paper.
   * 
   * \note The mesh is a 3x3 uniform grid of quads with 3 materials.
   * 
   * \return  The generated mesh.
   */
  void initTestCaseTwo(conduit::Node &mesh);

  /*!
   * \brief Initializes an MIRMesh used for testing triangle clipping cases.
   * 
   * \note The mesh is a set of four triangles with 2 materials
   * 
   * \return  The generated mesh.
   */
  void initTestCaseThree(conduit::Node &mesh);

  /*!
   * \brief Intializes a mesh used for testing a single circle of one materials surrounded by another.
   * 
   * \note The mesh is a 3x3 uniform grid with 2 materials and has a single circle in the center composed
   *       of one material and is surrounded by a second material.
   * 
   * \return  The generated mesh.
   */
  void initTestCaseFour(conduit::Node &mesh);

  /*!
   * \brief Initializes a mesh to be used for testing a set of concentric circles centered in a uniform 2D grid.
   * 
   * \param gridSize  The number of elements in the width and the height of the uniform grid.
   * \param numCircles  The number of concentric circles that are centered in the grid.
   * 
   * \note  Each circle is composed of a different material.
   * 
   * \return  The generated mesh.
   */
  void initTestCaseFive(int gridSize, int numCircles, conduit::Node &mesh);

  /*!
   * \brief Initializes a mesh to be used for testing a set of concentric spheres centered in a uniform 3D grid.
   * 
   * \param gridSize  The number of elements in the width and the height of the uniform grid.
   * \param numSpheres  The number of concentric spheres that are centered in the grid.
   * 
   * \note  Each sphere is composed of a different material.
   * 
   * \return  The generated mesh.
   */
  void initTestCaseSix(int gridSize, int numSpheres, conduit::Node &mesh);

  /*!
   * \brief Initializes a mesh composed of a uniform grid with a circle of material in it.
   * 
   * \param gridSize  The number of elements in the width and height of the uniform grid.
   * \param circleCenter  The center point of the circle.
   * \param circleRadius  The radius of the circle.
   * 
   * \return  The generated mesh.
   */
  void createUniformGridTestCaseMesh(int gridSize,
                                     const Point2 &circleCenter,
                                     axom::float64 circleRadius,
                                     conduit::Node &mesh);
  /*!
   * \brief Initializes a mesh to be used for validating the results of quad clipping.
   * 
   * \note The mesh is a 3x3 uniform grid with 2 materials and element volume fraction such
   *       that the mesh would be split horizontally through the middle.
   * 
   * \return  The generated mesh.
   */
  void initQuadClippingTestMesh(conduit::Node &mesh);

private:
  /*!
   * \brief make a 3x3 mesh of quads.
   * \param mesh A conduit node that will contain the new mesh.
   */
  void mesh3x3(conduit::Node &mesh);

  /*!
   * \brief Generates a 2D uniform grid of n x n elements.
   * 
   * \param gridSize  The number of elements in the width and height of the uniform grid.
   */
  void generateGrid(int gridSize, conduit::Node &mesh);

  /*!
   * \brief Generates a 3D uniform grid of n x n x n elements.
   * 
   * \param gridSize  The number of elements in the width, height, and depth of the uniform grid.
   */
  void generateGrid3D(int gridSize, conduit::Node &mesh);

  /*!
   * \brief Calculates the number of corners of the quad that are within the circle.
   * 
   * \param circleCenter The center point of the circle.
   * \param circleRadius The radius of the circle.
   * \param quadP0  The upper left vertex of the quad.
   * \param quadP1  The lower left vertex of the quad.
   * \param quadP2  The lower right vertex of the quad.
   * \param quadP3  The upper right vertex of the quad.
   * 
   * \return The number of corners of the quad that are within the circle.
   */
  int circleQuadCornersOverlaps(const Point2 &circleCenter,
                                axom::float64 circleRadius,
                                const Point2 &quadP0,
                                const Point2 &quadP1,
                                const Point2 &quadP2,
                                const Point2 &quadP3);

  bool m_structured {false};
};

}  // namespace data
}  // namespace bump
}  // namespace axom

#endif
