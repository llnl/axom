// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 *
 * \brief Illustrates how to construct and use a mint::UniformMesh object.
 */

// sphinx_tutorial_basic_example_start

// sphinx_tutorial_walkthrough_includes_start

#include "axom/config.hpp"                          // compile time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits

#include "axom/mint.hpp"                  // for mint classes and functions
#include "axom/core/numerics/Matrix.hpp"  // for numerics::Matrix

// sphinx_tutorial_walkthrough_includes_end

// namespace aliases
namespace mint = axom::mint;
namespace numerics = axom::numerics;
namespace xargs = mint::xargs;

using IndexType = axom::IndexType;

// compile-time switch for execution policy
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
constexpr int NUM_BLOCKS = 512;
using ExecPolicy = axom::CUDA_EXEC<NUM_BLOCKS>;
#elif defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
constexpr int NUM_BLOCKS = 512;
using ExecPolicy = axom::HIP_EXEC<NUM_BLOCKS>;
#elif defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
using ExecPolicy = axom::OMP_EXEC;
#else
using ExecPolicy = axom::SEQ_EXEC;
#endif

constexpr IndexType NUM_NODES_PER_CELL = 4;
constexpr double ONE_OVER_4 = 1. / static_cast<double>(NUM_NODES_PER_CELL);

/*!
 * \brief Holds command-line arguments
 */
static struct
{
  int res;
  bool useUnstructured;
} Arguments;

//------------------------------------------------------------------------------
// FUNCTION PROTOTYPES
//------------------------------------------------------------------------------
void parse_args(int argc, char** argv);
mint::Mesh* getUniformMesh();
mint::Mesh* getUnstructuredMesh();

//------------------------------------------------------------------------------
// PROGRAM MAIN
//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  parse_args(argc, argv);

  // sphinx_tutorial_walkthrough_set_memory_start
  // NOTE: use unified memory if we are using CUDA or HIP
  const int allocID = axom::execution_space<ExecPolicy>::allocatorID();
  axom::setDefaultAllocator(allocID);
  // sphinx_tutorial_walkthrough_set_memory_end

  // sphinx_tutorial_walkthrough_construct_mesh_start

  mint::Mesh* mesh = (Arguments.useUnstructured) ? getUnstructuredMesh() : getUniformMesh();

  // sphinx_tutorial_walkthrough_construct_mesh_end

  // sphinx_tutorial_walkthrough_add_fields_start

  // add a cell-centered and a node-centered field
  double* phi = mesh->createField<double>("phi", mint::NODE_CENTERED);
  double* hc = mesh->createField<double>("hc", mint::CELL_CENTERED);

  constexpr int NUM_COMPONENTS = 2;
  double* xc = mesh->createField<double>("xc", mint::CELL_CENTERED, NUM_COMPONENTS);

  // sphinx_tutorial_walkthrough_add_fields_end

  // sphinx_tutorial_walkthrough_compute_hf_start

  // loop over the nodes and evaluate Himmelblaus Function
  mint::for_all_nodes<ExecPolicy, xargs::xy>(
    mesh,
    AXOM_LAMBDA(IndexType nodeIdx, double x, double y) {
      const double x_2 = x * x;
      const double y_2 = y * y;
      const double A = x_2 + y - 11.0;
      const double B = x + y_2 - 7.0;

      phi[nodeIdx] = A * A + B * B;
    });

  // sphinx_tutorial_walkthrough_compute_hf_end

  // sphinx_tutorial_walkthrough_cell_centers_start

  // loop over cells and compute cell centers
  mint::for_all_cells<ExecPolicy, xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx, const numerics::Matrix<double>& coords, const IndexType* nodeIds) {
      // NOTE: A column vector of the coords matrix corresponds to a nodes coords

      // Sum the cell's nodal coordinates
      double xsum = 0.0;
      double ysum = 0.0;
      double hsum = 0.0;

      const IndexType numNodes = coords.getNumColumns();
      for(IndexType inode = 0; inode < numNodes; ++inode)
      {
        const double* node = coords.getColumn(inode);
        xsum += node[mint::X_COORDINATE];
        ysum += node[mint::Y_COORDINATE];

        hsum += phi[nodeIds[inode]];
      }  // END for all cell nodes

      // compute cell centroid by averaging the nodal coordinate sums
      const IndexType offset = cellIdx * NUM_COMPONENTS;
      const double invnnodes = 1.f / static_cast<double>(numNodes);
      xc[offset] = xsum * invnnodes;
      xc[offset + 1] = ysum * invnnodes;

      hc[cellIdx] = hsum * invnnodes;
    });

  // sphinx_tutorial_walkthrough_cell_centers_end

  // sphinx_tutorial_walkthrough_vtk_start

  // write the mesh in a VTK file for visualization
  std::string vtkfile = (Arguments.useUnstructured) ? "unstructured_mesh.vtk" : "uniform_mesh.vtk";
  mint::write_vtk(mesh, vtkfile);

  // sphinx_tutorial_walkthrough_vtk_end

  delete mesh;
  mesh = nullptr;

  return 0;
}

//------------------------------------------------------------------------------
//  FUNCTION PROTOTYPE IMPLEMENTATION
//------------------------------------------------------------------------------
void parse_args(int argc, char** argv)
{
  Arguments.res = 25;
  Arguments.useUnstructured = false;

  for(int i = 1; i < argc; ++i)
  {
    if(strcmp(argv[i], "--unstructured") == 0)
    {
      Arguments.useUnstructured = true;
    }

    else if(strcmp(argv[i], "--resolution") == 0)
    {
      Arguments.res = std::atoi(argv[++i]);
    }

  }  // END for all arguments

  SLIC_ERROR_IF(Arguments.res < 2, "invalid mesh resolution! Please, pick a value greater than 2.");
}

//------------------------------------------------------------------------------
// sphinx_tutorial_walkthrough_construct_umesh_start
mint::Mesh* getUniformMesh()
{
  // construct a N x N grid within a domain defined in [-5.0, 5.0]
  const double lo[] = {-5.0, -5.0};
  const double hi[] = {5.0, 5.0};
  mint::Mesh* m = new mint::UniformMesh(lo, hi, Arguments.res, Arguments.res);
  return m;
}
// sphinx_tutorial_walkthrough_construct_umesh_end

//------------------------------------------------------------------------------
mint::Mesh* getUnstructuredMesh()
{
  mint::Mesh* umesh = getUniformMesh();
  const IndexType umesh_ncells = umesh->getNumberOfCells();
  const IndexType umesh_nnodes = umesh->getNumberOfNodes();

  const IndexType ncells = umesh_ncells * 4;  // split each quad into 4 triangles
  const IndexType nnodes = umesh_nnodes + umesh_ncells;

  constexpr int DIMENSION = 2;
  using MeshType = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  MeshType* m = new MeshType(DIMENSION, mint::TRIANGLE, nnodes, ncells);
  m->resize(nnodes, ncells);

  double* x = m->getCoordinateArray(mint::X_COORDINATE);
  double* y = m->getCoordinateArray(mint::Y_COORDINATE);
  IndexType* cells = m->getCellNodesArray();

  // fill coordinates from uniform mesh
  mint::for_all_nodes<ExecPolicy, xargs::xy>(
    umesh,
    AXOM_LAMBDA(IndexType nodeIdx, double nx, double ny) {
      x[nodeIdx] = nx;
      y[nodeIdx] = ny;
    });

  // loop over cells, compute cell centers and fill connectivity
  mint::for_all_cells<ExecPolicy, xargs::coords>(
    umesh,
    AXOM_LAMBDA(IndexType cellIdx, const numerics::Matrix<double>& coords, const IndexType* nodeIds) {
      // NOTE: A column vector of the coords matrix corresponds to a nodes coords

      // Sum the cell's nodal coordinates
      double xsum = 0.0;
      double ysum = 0.0;
      for(IndexType inode = 0; inode < NUM_NODES_PER_CELL; ++inode)
      {
        const double* node = coords.getColumn(inode);
        xsum += node[mint::X_COORDINATE];
        ysum += node[mint::Y_COORDINATE];
      }  // END for all cell nodes

      // compute cell centroid by averaging the nodal coordinate sums
      const IndexType nc = umesh_nnodes + cellIdx; /* centroid index */
      x[nc] = xsum * ONE_OVER_4;
      y[nc] = ysum * ONE_OVER_4;

      // triangulate
      const IndexType& n0 = nodeIds[0];
      const IndexType& n1 = nodeIds[1];
      const IndexType& n2 = nodeIds[2];
      const IndexType& n3 = nodeIds[3];

      const IndexType offset = cellIdx * 12;

      cells[offset] = n0;
      cells[offset + 1] = nc;
      cells[offset + 2] = n3;

      cells[offset + 3] = n0;
      cells[offset + 4] = n1;
      cells[offset + 5] = nc;

      cells[offset + 6] = n1;
      cells[offset + 7] = n2;
      cells[offset + 8] = nc;

      cells[offset + 9] = n2;
      cells[offset + 10] = n3;
      cells[offset + 11] = nc;
    });

  // delete uniform mesh
  delete umesh;
  umesh = nullptr;

  return m;
}

// sphinx_tutorial_basic_example_end
