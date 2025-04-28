// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 *
 * \brief Illustrates how to construct and use a mixed topology
 *  UnstructuredMesh by creating a 2D structured grid with both quads
 *  and triangles.
 */

// Axom includes
#include "axom/mint.hpp"
#include "axom/slic.hpp"

// C/C++ includes
#include <random>

using namespace axom;

inline bool appendQuad(axom::IndexType i, axom::IndexType j) { return (i + j) % 2 == 0; }

//------------------------------------------------------------------------------
int main(int AXOM_UNUSED_PARAM(argc), char** AXOM_UNUSED_PARAM(argv))
{
  axom::slic::SimpleLogger logger;  // create & initialize test logger

  constexpr int DIMENSION = 2;
  constexpr axom::IndexType X_EXTENT = 11;
  constexpr axom::IndexType Y_EXTENT = 11;
  constexpr double SPACING = 1.0;
  constexpr axom::IndexType NUM_NODES = X_EXTENT * Y_EXTENT;
  constexpr axom::IndexType NUM_CELLS = (X_EXTENT - 1) * (Y_EXTENT - 1) * 3 / 2;

  constexpr double HI = -10.0;
  constexpr double LO = 10.0;
  constexpr double PLO = 0.0;
  constexpr double PHI = 1.0;

  /* STEP 0: create the UnstructuredMesh */
  mint::UnstructuredMesh<mint::MIXED_SHAPE> mesh(DIMENSION, NUM_NODES, NUM_CELLS);

  /* STEP 1: Add fields to the nodes and cells.
   * Note that we can only use the pointers below because we specified
   * the node and cell capacities in the constructor. In general calling
   * append or insert can invalidate any associated pointers. */
  double* vx = mesh.createField<double>("vx", mint::NODE_CENTERED);
  double* vy = mesh.createField<double>("vy", mint::NODE_CENTERED);
  double* vz = mesh.createField<double>("vz", mint::NODE_CENTERED);
  double* p = mesh.createField<double>("pressure", mint::CELL_CENTERED);
  double* p_avg = mesh.createField<double>("average_pressure", mint::NODE_CENTERED);
  int* cells_per_node = mesh.createField<int>("cells_per_node", mint::NODE_CENTERED);

  /* STEP 3: Add the nodes */
  axom::IndexType node_ID = 0;
  for(axom::IndexType j = 0; j < Y_EXTENT; ++j)
  {
    for(axom::IndexType i = 0; i < X_EXTENT; ++i)
    {
      mesh.appendNode(i * SPACING, j * SPACING);

      vx[node_ID] = utilities::random_real(LO, HI);
      vy[node_ID] = utilities::random_real(LO, HI);
      vz[node_ID] = utilities::random_real(LO, HI);
      p_avg[node_ID] = 0.0;
      cells_per_node[node_ID] = 0;
      node_ID++;
    }
  }

  /* STEP 4: Add the cells */
  axom::IndexType cell_ID = 0;
  for(axom::IndexType j = 0; j < Y_EXTENT - 1; ++j)
  {
    for(axom::IndexType i = 0; i < X_EXTENT - 1; ++i)
    {
      const axom::IndexType bottom_left = j * X_EXTENT + i;
      const axom::IndexType bottom_right = bottom_left + 1;
      const axom::IndexType top_right = bottom_right + X_EXTENT;
      const axom::IndexType top_left = bottom_left + X_EXTENT;

      if(appendQuad(i, j))
      {
        /* Append a quad. */
        const axom::IndexType quad[4] = {bottom_left, bottom_right, top_right, top_left};
        mesh.appendCell(quad, mint::QUAD);

        p[cell_ID] = utilities::random_real(PLO, PHI);
        cell_ID++;
      }
      else
      {
        /* Append two triangles. */
        const axom::IndexType tri0[3] = {bottom_left, bottom_right, top_right};
        mesh.appendCell(tri0, mint::TRIANGLE);
        p[cell_ID] = utilities::random_real(PLO, PHI);
        cell_ID++;

        const axom::IndexType tri1[3] = {top_right, top_left, bottom_left};
        mesh.appendCell(tri1, mint::TRIANGLE);
        p[cell_ID] = utilities::random_real(PLO, PHI);
        cell_ID++;
      }
    }
  }

  /* STEP 5: Calculate the average pressure at each of the nodes. */
  const axom::IndexType n_cells = mesh.getNumberOfCells();
  for(cell_ID = 0; cell_ID < n_cells; ++cell_ID)
  {
    const double cell_pressure = p[cell_ID];
    const int n_cell_nodes = mesh.getNumberOfCellNodes(cell_ID);
    const axom::IndexType* connec = mesh.getCellNodeIDs(cell_ID);
    for(int i = 0; i < n_cell_nodes; ++i)
    {
      node_ID = connec[i];
      p_avg[node_ID] += cell_pressure;
      cells_per_node[node_ID]++;
    }
  }

  const axom::IndexType n_nodes = mesh.getNumberOfNodes();
  for(node_ID = 0; node_ID < n_nodes; ++node_ID)
  {
    p_avg[node_ID] /= cells_per_node[node_ID];
  }

  /* STEP 6: write the particle mesh in VTK format for visualization. */
  mint::write_vtk(&mesh, "UnstructuredMesh_mixed_topology.vtk");

  return 0;
}
