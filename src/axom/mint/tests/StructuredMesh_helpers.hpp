// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef MINT_STRUCTURED_MESH_TEST_HELPERS_HPP_
#define MINT_STRUCTURED_MESH_TEST_HELPERS_HPP_

#include "axom/mint/config.hpp"  // for compile-time definitions

#include "axom/mint/mesh/blueprint.hpp"       // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"       // for CellTypes enum definition
#include "axom/mint/mesh/UniformMesh.hpp"     // for UniformMesh
#include "axom/mint/mesh/StructuredMesh.hpp"  // for StructuredMesh

#include "axom/slic/interface/slic.hpp"  // for slic macros

#include "gtest/gtest.h"  // for gtest macros

/*!
 * \file StructuredMesh_helpers.hpp
 * The functions declared in this header file test the methods of the
 * StructuredMesh class and are used in mint_mesh_uniform_mesh.cpp,
 * mint_mesh_rectilinear_mesh.cpp, and mint_mesh_curvilinear_mesh.cpp.
 */

namespace axom
{
namespace mint
{
namespace internal
{
constexpr double PI = 3.14;

// Always define these sizes for 3D because some of the mesh methods in mint
// meshes do not know their dimension until runtime and may contain 1D, 2D, 3D
// implementations so the passed buffers better be large enough - and to silence
// potential compiler warnings to that effect.
constexpr int MAX_FACE_NODES = 6;
constexpr int MAX_CELL_NODES = 8;

namespace
{
/*!
 * \brief Compare the values of two arrays.
 *
 * \param [in] p1 the first pointer to compare.
 * \param [in] p2 the second pointer to compare.
 * \param [in] n the number of elements to compare.
 */
template <typename T>
void compare_pointers(const T* p1, const T* p2, IndexType n)
{
  for(IndexType i = 0; i < n; ++i)
  {
    EXPECT_EQ(p1[i], p2[i]);
  }
}

/*!
 * \brief Check the indexing methods of a 2D StructuredMesh.
 *
 * \param [in] m the mesh to check.
 */
inline void check_indexing2D(const StructuredMesh* m)
{
  SLIC_ASSERT(m != nullptr);
  SLIC_ASSERT(m->getDimension() == 2);

  /* Check getNode Linear/Grid Index */
  IndexType idx = 0;
  for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
    {
      const IndexType nodeID = m->getNodeLinearIndex(i, j);
      EXPECT_EQ(nodeID, idx++);

      IndexType ii, jj, kk;
      m->getNodeGridIndex(nodeID, ii, jj);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);

      m->getNodeGridIndex(nodeID, ii, jj, kk);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);
      EXPECT_EQ(0, kk);
    }
  }

  /* Check getCell Linear/Grid Index */
  idx = 0;
  for(IndexType j = 0; j < m->getCellResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getCellResolution(0); ++i)
    {
      const IndexType cellID = m->getCellLinearIndex(i, j);
      EXPECT_EQ(cellID, idx++);

      IndexType ii, jj, kk;
      m->getCellGridIndex(cellID, ii, jj);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);

      m->getCellGridIndex(cellID, ii, jj, kk);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);
      EXPECT_EQ(0, kk);
    }
  }

  /* Check the getIFace Linear/Grid Index. */
  idx = 0;
  for(IndexType j = 0; j < m->getCellResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
    {
      const IndexType faceID = m->getIFaceLinearIndex(i, j);
      const IndexType faceID2 = m->getFaceLinearIndex(I_DIRECTION, i, j);
      EXPECT_EQ(faceID, idx++);
      EXPECT_EQ(faceID, faceID2);

      IndexType ii, jj, kk;
      m->getIFaceGridIndex(faceID, ii, jj);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);

      m->getIFaceGridIndex(faceID, ii, jj, kk);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);
      EXPECT_EQ(0, kk);
    }
  }

  /* Check the getJFace Linear/Grid Index. */
  for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getCellResolution(0); ++i)
    {
      const IndexType faceID = m->getJFaceLinearIndex(i, j);
      const IndexType faceID2 = m->getFaceLinearIndex(J_DIRECTION, i, j);
      EXPECT_EQ(faceID, idx++);
      EXPECT_EQ(faceID, faceID2);

      IndexType ii, jj, kk;
      m->getJFaceGridIndex(faceID, ii, jj);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);

      m->getJFaceGridIndex(faceID, ii, jj, kk);
      EXPECT_EQ(i, ii);
      EXPECT_EQ(j, jj);
      EXPECT_EQ(0, kk);
    }
  }
}

/*!
 * \brief Check the indexing methods of a 3D StructuredMesh.
 *
 * \param [in] m the mesh to check.
 */
inline void check_indexing3D(const StructuredMesh* m)
{
  SLIC_ASSERT(m != nullptr);
  SLIC_ASSERT(m->getDimension() == 3);

  /* Check getNode Linear/Grid Index */
  IndexType idx = 0;
  for(IndexType k = 0; k < m->getNodeResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
      {
        const IndexType nodeID = m->getNodeLinearIndex(i, j, k);
        EXPECT_EQ(nodeID, idx++);

        IndexType ii, jj, kk;
        m->getNodeGridIndex(nodeID, ii, jj, kk);
        EXPECT_EQ(i, ii);
        EXPECT_EQ(j, jj);
        EXPECT_EQ(k, kk);
      }
    }
  }

  /* Check getCell Linear/Grid Index */
  idx = 0;
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType cellID = m->getCellLinearIndex(i, j, k);
        EXPECT_EQ(cellID, idx++);

        IndexType ii, jj, kk;
        m->getCellGridIndex(cellID, ii, jj, kk);
        EXPECT_EQ(i, ii);
        EXPECT_EQ(j, jj);
        EXPECT_EQ(k, kk);
      }
    }
  }

  /* Check the getIFace Linear/Grid Index. */
  idx = 0;
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
      {
        const IndexType faceID = m->getIFaceLinearIndex(i, j, k);
        const IndexType faceID2 = m->getFaceLinearIndex(I_DIRECTION, i, j, k);
        EXPECT_EQ(faceID, idx++);
        EXPECT_EQ(faceID, faceID2);

        IndexType ii, jj, kk;
        m->getIFaceGridIndex(faceID, ii, jj, kk);
        EXPECT_EQ(i, ii);
        EXPECT_EQ(j, jj);
        EXPECT_EQ(k, kk);
      }
    }
  }

  /* Check the getJFace Linear/Grid Index. */
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType faceID = m->getJFaceLinearIndex(i, j, k);
        const IndexType faceID2 = m->getFaceLinearIndex(J_DIRECTION, i, j, k);
        EXPECT_EQ(faceID, idx++);
        EXPECT_EQ(faceID, faceID2);

        IndexType ii, jj, kk;
        m->getJFaceGridIndex(faceID, ii, jj, kk);
        EXPECT_EQ(i, ii);
        EXPECT_EQ(j, jj);
        EXPECT_EQ(k, kk);
      }
    }
  }

  /* Check the getKFace Linear/Grid Index. */
  for(IndexType k = 0; k < m->getNodeResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType faceID = m->getKFaceLinearIndex(i, j, k);
        const IndexType faceID2 = m->getFaceLinearIndex(K_DIRECTION, i, j, k);
        EXPECT_EQ(faceID, idx++);
        EXPECT_EQ(faceID, faceID2);

        IndexType ii, jj, kk;
        m->getKFaceGridIndex(faceID, ii, jj, kk);
        EXPECT_EQ(i, ii);
        EXPECT_EQ(j, jj);
        EXPECT_EQ(k, kk);
      }
    }
  }
}

/*!
 * \brief Check the topology related methods of a 1D StructuredMesh.
 *
 * \param [in] m the mesh to check.
 *
 * \note Since 1D meshes don't have any faces the only method tested here
 *  is getCellNodeIDs.
 */
inline void check_topology1D(const StructuredMesh* m)
{
  SLIC_ASSERT(m != nullptr);
  SLIC_ASSERT(m->getDimension() == 1);

  EXPECT_EQ(m->getNumberOfCellNodes(), 2);
  EXPECT_EQ(m->getNumberOfCellFaces(), 0);
  EXPECT_EQ(m->getNumberOfFaces(), 0);
  EXPECT_EQ(m->getNumberOfFaceNodes(), 0);

  IndexType cellNodes[2];
  const IndexType nCells = m->getNumberOfCells();
  for(IndexType cellID = 0; cellID < nCells; ++cellID)
  {
    m->getCellNodeIDs(cellID, cellNodes);
    EXPECT_EQ(cellNodes[0], cellID);
    EXPECT_EQ(cellNodes[1], cellID + 1);
  }
}

/*!
 * \brief Check the topology related methods of a 2D StructuredMesh.
 *
 * \param [in] m the mesh to check.
 */
inline void check_topology2D(const StructuredMesh* m)
{
  SLIC_ASSERT(m != nullptr);
  SLIC_ASSERT(m->getDimension() == 2);

  EXPECT_EQ(m->getNumberOfCellNodes(), 4);

  /* Check the cell to node map. */
  IndexType cellNodes[MAX_CELL_NODES], cellNodesIJ[MAX_CELL_NODES];
  const IndexType nodeJp = m->nodeJp();
  for(IndexType j = 0; j < m->getCellResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getCellResolution(0); ++i)
    {
      const IndexType cellID = m->getCellLinearIndex(i, j);
      m->getCellNodeIDs(cellID, cellNodes);
      m->getCellNodeIDs(i, j, cellNodesIJ);
      compare_pointers(cellNodes, cellNodesIJ, 4);

      const IndexType nodeID = m->getNodeLinearIndex(i, j);
      EXPECT_EQ(cellNodes[0], nodeID);
      EXPECT_EQ(cellNodes[1], nodeID + 1);
      EXPECT_EQ(cellNodes[2], nodeID + nodeJp + 1);
      EXPECT_EQ(cellNodes[3], nodeID + nodeJp);
    }
  }

  EXPECT_EQ(m->getNumberOfCellFaces(), 4);
  EXPECT_EQ(m->getNumberOfFaceNodes(), 2);

  const IndexType nIFaces = m->getNodeResolution(0) * m->getCellResolution(1);
  const IndexType nJFaces = m->getCellResolution(0) * m->getNodeResolution(1);
  EXPECT_EQ(nIFaces, m->getTotalNumFaces(0));
  EXPECT_EQ(nJFaces, m->getTotalNumFaces(1));

  const IndexType nFaces = nIFaces + nJFaces;
  EXPECT_EQ(m->getNumberOfFaces(), nFaces);

  /* Check the cell to face map */
  const IndexType cellJp = m->cellJp();
  IndexType cellFaces[MAX_FACE_NODES], cellFacesIJ[MAX_FACE_NODES];
  for(IndexType j = 0; j < m->getCellResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getCellResolution(0); ++i)
    {
      const IndexType cellID = m->getCellLinearIndex(i, j);
      m->getCellFaceIDs(cellID, cellFaces);
      m->getCellFaceIDs(i, j, cellFacesIJ);
      compare_pointers(cellFaces, cellFacesIJ, 4);

      EXPECT_EQ(cellFaces[0], cellID + j);
      EXPECT_EQ(cellFaces[1], cellID + j + 1);
      EXPECT_EQ(cellFaces[2], cellID + nIFaces);
      EXPECT_EQ(cellFaces[3], cellID + nIFaces + cellJp);
    }
  }

  /* Check the face to node map. */
  IndexType faceNodes[MAX_FACE_NODES];
  IndexType faceNodesDir[MAX_FACE_NODES];

  /* I-faces */
  for(IndexType j = 0; j < m->getCellResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
    {
      const IndexType faceID = m->getIFaceLinearIndex(i, j);
      m->getFaceNodeIDs(faceID, faceNodes);
      m->getIFaceNodeIDs(faceID, faceNodesDir);
      compare_pointers(faceNodes, faceNodesDir, 2);

      EXPECT_EQ(faceNodes[0], m->getNodeLinearIndex(i, j));
      EXPECT_EQ(faceNodes[1], m->getNodeLinearIndex(i, j + 1));
    }
  }

  /* J-faces */
  for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getCellResolution(0); ++i)
    {
      const IndexType faceID = m->getJFaceLinearIndex(i, j);
      m->getFaceNodeIDs(faceID, faceNodes);
      m->getJFaceNodeIDs(faceID, faceNodesDir);
      compare_pointers(faceNodes, faceNodesDir, 2);

      EXPECT_EQ(faceNodes[0], m->getNodeLinearIndex(i, j));
      EXPECT_EQ(faceNodes[1], m->getNodeLinearIndex(i + 1, j));
    }
  }

  /* Check the face to cell map */
  IndexType cell0, cell1, cell2, cell3;

  /* I-faces */
  for(IndexType j = 0; j < m->getCellResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
    {
      const IndexType faceID = m->getIFaceLinearIndex(i, j);
      m->getFaceCellIDs(faceID, cell0, cell1);
      m->getIFaceCellIDs(faceID, cell2, cell3);
      EXPECT_EQ(cell0, cell2);
      EXPECT_EQ(cell1, cell3);

      /* Check the boundaries. */
      if(i == 0)
      {
        EXPECT_EQ(cell0, m->getCellLinearIndex(i, j));
        EXPECT_EQ(cell1, -1);
      }
      else if(i == m->getCellResolution(0))
      {
        EXPECT_EQ(cell0, m->getCellLinearIndex(i - 1, j));
        EXPECT_EQ(cell1, -1);
      }
      else
      {
        EXPECT_EQ(cell0, m->getCellLinearIndex(i - 1, j));
        EXPECT_EQ(cell1, m->getCellLinearIndex(i, j));
      }
    }
  }

  /* J-faces */
  for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
  {
    for(IndexType i = 0; i < m->getCellResolution(0); ++i)
    {
      const IndexType faceID = m->getJFaceLinearIndex(i, j);
      m->getFaceCellIDs(faceID, cell0, cell1);
      m->getJFaceCellIDs(faceID, cell2, cell3);
      EXPECT_EQ(cell0, cell2);
      EXPECT_EQ(cell1, cell3);

      /* Check the boundaries. */
      if(j == 0)
      {
        EXPECT_EQ(cell0, m->getCellLinearIndex(i, j));
        EXPECT_EQ(cell1, -1);
      }
      else if(j == m->getCellResolution(1))
      {
        EXPECT_EQ(cell0, m->getCellLinearIndex(i, j - 1));
        EXPECT_EQ(cell1, -1);
      }
      else
      {
        EXPECT_EQ(cell0, m->getCellLinearIndex(i, j - 1));
        EXPECT_EQ(cell1, m->getCellLinearIndex(i, j));
      }
    }
  }
}

/*!
 * \brief Check the topology related methods of a 3D StructuredMesh.
 *
 * \param [in] m the mesh to check.
 */
inline void check_topology3D(const StructuredMesh* m)
{
  SLIC_ASSERT(m != nullptr);
  SLIC_ASSERT(m->getDimension() == 3);

  EXPECT_EQ(m->getNumberOfCellNodes(), 8);

  /* Check the cell to node map. */
  IndexType cellNodes[MAX_CELL_NODES], cellNodesIJK[MAX_CELL_NODES];
  const IndexType nodeJp = m->nodeJp();
  const IndexType nodeKp = m->nodeKp();
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType cellID = m->getCellLinearIndex(i, j, k);
        m->getCellNodeIDs(cellID, cellNodes);
        m->getCellNodeIDs(i, j, k, cellNodesIJK);
        compare_pointers(cellNodes, cellNodesIJK, 8);

        const IndexType nodeID = m->getNodeLinearIndex(i, j, k);
        EXPECT_EQ(cellNodes[0], nodeID);
        EXPECT_EQ(cellNodes[1], nodeID + 1);
        EXPECT_EQ(cellNodes[2], nodeID + nodeJp + 1);
        EXPECT_EQ(cellNodes[3], nodeID + nodeJp);
        EXPECT_EQ(cellNodes[4], nodeID + nodeKp);
        EXPECT_EQ(cellNodes[5], nodeID + nodeKp + 1);
        EXPECT_EQ(cellNodes[6], nodeID + nodeKp + nodeJp + 1);
        EXPECT_EQ(cellNodes[7], nodeID + nodeKp + nodeJp);
      }
    }
  }

  EXPECT_EQ(m->getNumberOfCellFaces(), 6);
  EXPECT_EQ(m->getNumberOfFaceNodes(), 4);

  const IndexType nIFaces =
    m->getNodeResolution(0) * m->getCellResolution(1) * m->getCellResolution(2);
  const IndexType nJFaces =
    m->getCellResolution(0) * m->getNodeResolution(1) * m->getCellResolution(2);
  const IndexType nKFaces =
    m->getCellResolution(0) * m->getCellResolution(1) * m->getNodeResolution(2);
  EXPECT_EQ(nIFaces, m->getTotalNumFaces(0));
  EXPECT_EQ(nJFaces, m->getTotalNumFaces(1));
  EXPECT_EQ(nKFaces, m->getTotalNumFaces(2));

  const IndexType nFaces = nIFaces + nJFaces + nKFaces;
  EXPECT_EQ(m->getNumberOfFaces(), nFaces);

  /* Check the cell to face map */
  IndexType cellFaces[MAX_FACE_NODES], cellFacesIJK[MAX_FACE_NODES];
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType cellID = m->getCellLinearIndex(i, j, k);
        m->getCellFaceIDs(cellID, cellFaces);
        m->getCellFaceIDs(i, j, k, cellFacesIJK);
        compare_pointers(cellFaces, cellFacesIJK, 6);

        EXPECT_EQ(cellFaces[0], m->getIFaceLinearIndex(i, j, k));
        EXPECT_EQ(cellFaces[1], m->getIFaceLinearIndex(i + 1, j, k));
        EXPECT_EQ(cellFaces[2], m->getJFaceLinearIndex(i, j, k));
        EXPECT_EQ(cellFaces[3], m->getJFaceLinearIndex(i, j + 1, k));
        EXPECT_EQ(cellFaces[4], m->getKFaceLinearIndex(i, j, k));
        EXPECT_EQ(cellFaces[5], m->getKFaceLinearIndex(i, j, k + 1));
      }
    }
  }

  /* Check the face to node map */
  IndexType faceNodes[MAX_FACE_NODES];
  IndexType faceNodesDir[MAX_FACE_NODES];

  /* I-faces */
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
      {
        const IndexType faceID = m->getIFaceLinearIndex(i, j, k);
        m->getFaceNodeIDs(faceID, faceNodes);
        m->getIFaceNodeIDs(faceID, faceNodesDir);
        compare_pointers(faceNodes, faceNodesDir, 4);

        EXPECT_EQ(faceNodes[0], m->getNodeLinearIndex(i, j, k));
        EXPECT_EQ(faceNodes[1], m->getNodeLinearIndex(i, j, k + 1));
        EXPECT_EQ(faceNodes[2], m->getNodeLinearIndex(i, j + 1, k + 1));
        EXPECT_EQ(faceNodes[3], m->getNodeLinearIndex(i, j + 1, k));
      }
    }
  }

  /* J-faces */
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType faceID = m->getJFaceLinearIndex(i, j, k);
        m->getFaceNodeIDs(faceID, faceNodes);
        m->getJFaceNodeIDs(faceID, faceNodesDir);
        compare_pointers(faceNodes, faceNodesDir, 4);

        EXPECT_EQ(faceNodes[0], m->getNodeLinearIndex(i, j, k));
        EXPECT_EQ(faceNodes[1], m->getNodeLinearIndex(i + 1, j, k));
        EXPECT_EQ(faceNodes[2], m->getNodeLinearIndex(i + 1, j, k + 1));
        EXPECT_EQ(faceNodes[3], m->getNodeLinearIndex(i, j, k + 1));
      }
    }
  }

  /* K-faces */
  for(IndexType k = 0; k < m->getNodeResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType faceID = m->getKFaceLinearIndex(i, j, k);
        m->getFaceNodeIDs(faceID, faceNodes);
        m->getKFaceNodeIDs(faceID, faceNodesDir);
        compare_pointers(faceNodes, faceNodesDir, 4);

        EXPECT_EQ(faceNodes[0], m->getNodeLinearIndex(i, j, k));
        EXPECT_EQ(faceNodes[1], m->getNodeLinearIndex(i + 1, j, k));
        EXPECT_EQ(faceNodes[2], m->getNodeLinearIndex(i + 1, j + 1, k));
        EXPECT_EQ(faceNodes[3], m->getNodeLinearIndex(i, j + 1, k));
      }
    }
  }

  /* Check the face to cell map */
  IndexType cell0, cell1, cell2, cell3;

  /* I-faces */
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getNodeResolution(0); ++i)
      {
        const IndexType faceID = m->getIFaceLinearIndex(i, j, k);
        m->getFaceCellIDs(faceID, cell0, cell1);
        m->getIFaceCellIDs(faceID, cell2, cell3);
        EXPECT_EQ(cell0, cell2);
        EXPECT_EQ(cell1, cell3);

        /* Check the boundaries. */
        if(i == 0)
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i, j, k));
          EXPECT_EQ(cell1, -1);
        }
        else if(i == m->getCellResolution(0))
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i - 1, j, k));
          EXPECT_EQ(cell1, -1);
        }
        else
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i - 1, j, k));
          EXPECT_EQ(cell1, m->getCellLinearIndex(i, j, k));
        }
      }
    }
  }

  /* J-faces */
  for(IndexType k = 0; k < m->getCellResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getNodeResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType faceID = m->getJFaceLinearIndex(i, j, k);
        m->getFaceCellIDs(faceID, cell0, cell1);
        m->getJFaceCellIDs(faceID, cell2, cell3);
        EXPECT_EQ(cell0, cell2);
        EXPECT_EQ(cell1, cell3);

        /* Check the boundaries. */
        if(j == 0)
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i, j, k));
          EXPECT_EQ(cell1, -1);
        }
        else if(j == m->getCellResolution(1))
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i, j - 1, k));
          EXPECT_EQ(cell1, -1);
        }
        else
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i, j - 1, k));
          EXPECT_EQ(cell1, m->getCellLinearIndex(i, j, k));
        }
      }
    }
  }

  /* K-faces */
  for(IndexType k = 0; k < m->getNodeResolution(2); ++k)
  {
    for(IndexType j = 0; j < m->getCellResolution(1); ++j)
    {
      for(IndexType i = 0; i < m->getCellResolution(0); ++i)
      {
        const IndexType faceID = m->getKFaceLinearIndex(i, j, k);
        m->getFaceCellIDs(faceID, cell0, cell1);
        m->getKFaceCellIDs(faceID, cell2, cell3);
        EXPECT_EQ(cell0, cell2);
        EXPECT_EQ(cell1, cell3);

        /* Check the boundaries. */
        if(k == 0)
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i, j, k));
          EXPECT_EQ(cell1, -1);
        }
        else if(k == m->getCellResolution(2))
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i, j, k - 1));
          EXPECT_EQ(cell1, -1);
        }
        else
        {
          EXPECT_EQ(cell0, m->getCellLinearIndex(i, j, k - 1));
          EXPECT_EQ(cell1, m->getCellLinearIndex(i, j, k));
        }
      }
    }
  }
}

/*!
 * \brief Check the indexing and topology related methods of a StructuredMesh.
 *
 * \param [in] m the mesh to check.
 */
inline void check_topology(const StructuredMesh* m)
{
  SLIC_ASSERT(m != nullptr);

  const int ndims = m->getDimension();
  IndexType nCells = 1;
  IndexType nNodes = 1;
  for(int dim = 0; dim < ndims; ++dim)
  {
    EXPECT_EQ(m->getNodeResolution(dim), m->getCellResolution(dim) + 1);
    nNodes *= m->getNodeResolution(dim);
    nCells *= m->getCellResolution(dim);
  }

  EXPECT_EQ(nNodes, m->getNumberOfNodes());
  EXPECT_EQ(nCells, m->getNumberOfCells());

  if(ndims == 1)
  {
    check_topology1D(m);
  }
  else if(ndims == 2)
  {
    check_indexing2D(m);
    check_topology2D(m);
  }
  else
  {
    SLIC_ASSERT(ndims == 3);
    check_indexing3D(m);
    check_topology3D(m);
  }
}

} /* namespace */

/*!
 * \brief Create a field with the given parameters and check for errors.
 *
 * \param [in/out] m the mesh to check.
 * \param [in] association the association of the field.
 * \param [in] name the name of the field.
 * \param [in] numComponents the number of components (optional).
 */
inline void check_create_field(StructuredMesh* m,
                               int association,
                               const std::string& name,
                               int numComponents = 1)
{
  EXPECT_TRUE(m != nullptr);
  EXPECT_FALSE(m->hasField(name, association));

  double* f = m->createField<double>(name, association, numComponents);
  EXPECT_TRUE(m->hasField(name, association));

  IndexType expected_num_tuples;
  if(association == NODE_CENTERED)
  {
    expected_num_tuples = m->getNumberOfNodes();
  }
  else if(association == CELL_CENTERED)
  {
    expected_num_tuples = m->getNumberOfCells();
  }
  else if(association == FACE_CENTERED)
  {
    expected_num_tuples = m->getNumberOfFaces();
  }
  else
  {
    expected_num_tuples = m->getNumberOfEdges();
  }

  const Field* field = m->getFieldData(association)->getField(name);
  EXPECT_TRUE(field != nullptr);
  EXPECT_EQ(f, Field::getDataPtr<double>(field));
  EXPECT_EQ(numComponents, field->getNumComponents());
  EXPECT_EQ(expected_num_tuples, field->getNumTuples());

  for(IndexType i = 0; i < expected_num_tuples * numComponents; ++i)
  {
    f[i] = PI * i;
  }
}

/*!
 * \brief Check that the field has the proper values.
 *
 * \param [in] field the field to check.
 */
inline void check_field_values(const Field* field)
{
  const IndexType num_values = field->getNumTuples() * field->getNumComponents();
  const double* values = Field::getDataPtr<double>(field);
  for(IndexType i = 0; i < num_values; ++i)
  {
    EXPECT_EQ(values[i], PI * i);
  }
}

/*!
 * \brief Create a field on each of the topological objects.
 *
 * \param [in] m the StructuredMesh to create the fields on.
 */
inline void check_create_fields(StructuredMesh* m)
{
  check_create_field(m, NODE_CENTERED, "n1", 1);
  check_create_field(m, CELL_CENTERED, "c1", 2);
  check_create_field(m, FACE_CENTERED, "f1", 3);
  check_create_field(m, EDGE_CENTERED, "e1", 4);
}

/*!
 * \brief Check that the fields have the proper values.
 *
 * \param [in] m the StructuredMesh to check.
 * \param [in] isInSidre if the fields should be in sidre.
 */
inline void check_fields(const StructuredMesh* m, bool isInSidre)
{
  EXPECT_TRUE(m->hasField("n1", NODE_CENTERED));
  EXPECT_TRUE(m->hasField("c1", CELL_CENTERED));
  EXPECT_TRUE(m->hasField("f1", FACE_CENTERED));
  EXPECT_TRUE(m->hasField("e1", EDGE_CENTERED));

  /* check node-centered field on m. */
  const FieldData* fd = m->getFieldData(NODE_CENTERED);
  const Field* field = fd->getField("n1");
  EXPECT_EQ(field->getNumTuples(), m->getNumberOfNodes());
  EXPECT_EQ(field->getNumComponents(), 1);
  EXPECT_EQ(field->isInSidre(), isInSidre);
  EXPECT_FALSE(field->isExternal());
  internal::check_field_values(field);

  /* check cell-centered field on m. */
  fd = m->getFieldData(CELL_CENTERED);
  field = fd->getField("c1");
  EXPECT_EQ(field->getNumTuples(), m->getNumberOfCells());
  EXPECT_EQ(field->getNumComponents(), 2);
  EXPECT_EQ(field->isInSidre(), isInSidre);
  EXPECT_FALSE(field->isExternal());
  internal::check_field_values(field);

  /* check face-centered field on m. */
  fd = m->getFieldData(FACE_CENTERED);
  field = fd->getField("f1");
  EXPECT_EQ(field->getNumTuples(), m->getNumberOfFaces());
  EXPECT_EQ(field->getNumComponents(), 3);
  EXPECT_EQ(field->isInSidre(), isInSidre);
  EXPECT_FALSE(field->isExternal());
  internal::check_field_values(field);

  /* check edge-centered field on m. */
  fd = m->getFieldData(EDGE_CENTERED);
  field = fd->getField("e1");
  EXPECT_EQ(field->getNumTuples(), m->getNumberOfEdges());
  EXPECT_EQ(field->getNumComponents(), 4);
  EXPECT_EQ(field->isInSidre(), isInSidre);
  EXPECT_FALSE(field->isExternal());
  internal::check_field_values(field);
}

/*!
 * \brief Check that the StructuredMesh was constructed correctly.
 *
 * \param [in] m the StructuredMesh to check.
 * \param [in] mesh_type the expected type of the StructuredMesh.
 * \param [in] ndims the expected number of dimensions.
 * \param [in] node_dims the expected nodal dimensions of the mesh.
 */
inline void check_constructor(const StructuredMesh* m,
                              int mesh_type,
                              int ndims,
                              const IndexType* node_dims)
{
  SLIC_ASSERT(mesh_type == STRUCTURED_CURVILINEAR_MESH ||
              mesh_type == STRUCTURED_RECTILINEAR_MESH || mesh_type == STRUCTURED_UNIFORM_MESH);

  EXPECT_TRUE(m != nullptr);
  EXPECT_EQ(m->getMeshType(), mesh_type);
  EXPECT_TRUE(ndims >= 1 && ndims <= 3);

  const int mesh_dimension = m->getDimension();

  EXPECT_EQ(mesh_dimension, ndims);
  EXPECT_EQ(m->getMeshType(), mesh_type);
  EXPECT_FALSE(m->hasExplicitConnectivity());
  EXPECT_FALSE(m->hasMixedCellTypes());

  if(mesh_type == STRUCTURED_UNIFORM_MESH)
  {
    EXPECT_FALSE(m->hasExplicitCoordinates());
  }
  else
  {
    EXPECT_TRUE(m->hasExplicitCoordinates());
  }

  CellType cell_type = (mesh_dimension == 3) ? HEX : (mesh_dimension == 2) ? QUAD : SEGMENT;
  EXPECT_EQ(m->getCellType(), cell_type);
  EXPECT_EQ(m->getNumberOfCellNodes(), getCellInfo(cell_type).num_nodes);
  EXPECT_EQ(m->getNumberOfCellFaces(), getCellInfo(cell_type).num_faces);

  for(int i = 0; i < mesh_dimension; ++i)
  {
    EXPECT_EQ(m->getNodeResolution(i), node_dims[i]);

    if(mesh_type != STRUCTURED_UNIFORM_MESH)
    {
      EXPECT_TRUE(m->getCoordinateArray(i) != nullptr);
    }
  }

  check_topology(m);
}

/*!
 * \brief Check that the UniformMesh was constructed correctly.
 *
 * \param [in] m the UniformMesh to check.
 * \param [in] ndims the expected number of dimensions.
 * \param [in] origin the origin of the UniformMesh.
 * \param [in] spacing the spacing between the nodes.
 * \param [in] node_dims the expected nodal dimensions of the mesh.
 */
inline void check_constructor(const UniformMesh* m,
                              int ndims,
                              const double* origin,
                              const double* spacing,
                              const IndexType* node_dims)
{
  EXPECT_TRUE(m != nullptr);
  EXPECT_TRUE(ndims >= 1 && ndims <= 3);

  for(int i = 0; i < ndims; ++i)
  {
    EXPECT_DOUBLE_EQ(origin[i], m->getOrigin()[i]);
    EXPECT_DOUBLE_EQ(spacing[i], m->getSpacing()[i]);
  }

  check_constructor(m, STRUCTURED_UNIFORM_MESH, ndims, node_dims);
}

/*!
 * \brief Check that the nodal extent of the StructuredMesh was set correctly.
 *
 * \param [in] m the StructuredMesh to check.
 * \param [in] extent the expected extent of the mesh.
 */
inline void check_node_extent(const StructuredMesh* m, const int64* extent)
{
  SLIC_ASSERT(m != nullptr);
  SLIC_ASSERT(extent != nullptr);

  int64 m_extent[6];
  m->getExtent(m_extent);

  for(int dim = 0; dim < m->getDimension(); ++dim)
  {
    EXPECT_EQ(m_extent[2 * dim], extent[2 * dim]);
    EXPECT_EQ(m_extent[2 * dim + 1], extent[2 * dim + 1]);
  }

  for(int dim = m->getDimension(); dim < 3; ++dim)
  {
    EXPECT_EQ(m_extent[2 * dim], m_extent[2 * dim + 1]);
  }
}

} /* namespace internal */
} /* namespace mint */
} /* namespace axom */

#endif /* MINT_STRUCTURED_MESH_TEST_HELPERS_HPP_ */
