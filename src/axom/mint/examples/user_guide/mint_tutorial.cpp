// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"  // for axom macros
#include "axom/mint.hpp"  // for Mint classes & functions

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
namespace sidre = axom::sidre;
#endif

// C/C++ includes
#include <cmath>  // for definition of M_PI, exp()

// namespace aliases
namespace mint = axom::mint;
namespace numerics = axom::numerics;

/*!
 * \file
 *
 * \brief Various code snippets/examples used for the Mint tutorial section.
 *
 * \note These examples are designed to illustrate specific Mint
 *  concepts and capabilities. Consult the Tutorial section of Mint's
 *  User Guide for more details.
 *
 */

//------------------------------------------------------------------------------
//    TUTORIAL FUNCTIONS
//------------------------------------------------------------------------------
void vtk_output(mint::Mesh* mesh, const std::string fileName)
{
  // sphinx_tutorial_vtk_output_start

  mint::write_vtk(mesh, fileName);

  // sphinx_tutorial_vtk_output_end
}

//------------------------------------------------------------------------------
void node_traversals()
{
  using exec_policy = axom::SEQ_EXEC;
  using IndexType = axom::IndexType;

  const double lo[] = {-5.0, -5.0};
  const double hi[] = {5.0, 5.0};
  mint::UniformMesh mesh(lo, hi, 50, 50);

  mesh.createField<double>("xold", mint::NODE_CENTERED);
  mesh.createField<double>("yold", mint::NODE_CENTERED);
  mesh.createField<double>("vx", mint::NODE_CENTERED);
  mesh.createField<double>("vy", mint::NODE_CENTERED);
  mesh.createField<double>("vmag", mint::NODE_CENTERED);
  mesh.createField<IndexType>("ID", mint::NODE_CENTERED);

  // Index example
  {
    // sphinx_tutorial_for_all_nodes_index_start

    const double* vx = mesh.getFieldPtr<double>("vx", mint::NODE_CENTERED);
    const double* vy = mesh.getFieldPtr<double>("vy", mint::NODE_CENTERED);

    double* vmag = mesh.getFieldPtr<double>("vmag", mint::NODE_CENTERED);

    mint::for_all_nodes<exec_policy>(
      &mesh,
      AXOM_LAMBDA(IndexType nodeIdx) {
        const double vx2 = vx[nodeIdx] * vx[nodeIdx];
        const double vy2 = vy[nodeIdx] * vy[nodeIdx];
        vmag[nodeIdx] = sqrt(vx2 + vy2);
      });

    // sphinx_tutorial_for_all_nodes_index_end
  }

  // IJ examples
  {
    // sphinx_tutorial_for_all_nodes_ij_start

    const IndexType jp = mesh.nodeJp();

    IndexType* ID = mesh.getFieldPtr<IndexType>("ID", mint::NODE_CENTERED);
    mint::for_all_nodes<exec_policy, mint::xargs::ij>(
      &mesh,
      AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j) { ID[nodeIdx] = i + j * jp; });

    // sphinx_tutorial_for_all_nodes_ij_end
  }

  // XY example
  {
    // sphinx_tutorial_for_all_nodes_xy_start

    double invdt = 0.5;

    const double* xold = mesh.getFieldPtr<double>("xold", mint::NODE_CENTERED);
    const double* yold = mesh.getFieldPtr<double>("yold", mint::NODE_CENTERED);

    double* vx = mesh.getFieldPtr<double>("vx", mint::NODE_CENTERED);
    double* vy = mesh.getFieldPtr<double>("vy", mint::NODE_CENTERED);

    mint::for_all_nodes<exec_policy, mint::xargs::xy>(
      &mesh,
      AXOM_LAMBDA(IndexType nodeIdx, double x, double y) {
        vx[nodeIdx] = (x - xold[nodeIdx]) * invdt;
        vy[nodeIdx] = (y - yold[nodeIdx]) * invdt;
      });

    // sphinx_tutorial_for_all_nodes_xy_end
  }
}

//------------------------------------------------------------------------------
void cell_traversals()
{
  using exec_policy = axom::SEQ_EXEC;
  using IndexType = axom::IndexType;

  const double lo[] = {-5.0, -5.0};
  const double hi[] = {5.0, 5.0};
  mint::UniformMesh mesh(lo, hi, 50, 50);

  mesh.createField<double>("vx", mint::NODE_CENTERED);
  mesh.createField<double>("vy", mint::NODE_CENTERED);

  mesh.createField<double>("cell_vx", mint::CELL_CENTERED);
  mesh.createField<double>("cell_vy", mint::CELL_CENTERED);
  mesh.createField<double>("den", mint::CELL_CENTERED);
  mesh.createField<double>("mass", mint::CELL_CENTERED);
  mesh.createField<double>("vol", mint::CELL_CENTERED);
  mesh.createField<double>("xc", mint::CELL_CENTERED);
  mesh.createField<double>("yc", mint::CELL_CENTERED);
  mesh.createField<double>("perimeter", mint::CELL_CENTERED);
  mesh.createField<IndexType>("ID", mint::CELL_CENTERED);

  mesh.createField<double>("area", mint::FACE_CENTERED);
  // Index example
  {
    // sphinx_tutorial_for_all_cells_index_start

    const double* mass = mesh.getFieldPtr<double>("mass", mint::CELL_CENTERED);
    const double* vol = mesh.getFieldPtr<double>("vol", mint::CELL_CENTERED);

    double* den = mesh.getFieldPtr<double>("den", mint::CELL_CENTERED);

    mint::for_all_cells<exec_policy>(
      &mesh,
      AXOM_LAMBDA(IndexType cellIdx) { den[cellIdx] = mass[cellIdx] / vol[cellIdx]; });

    // sphinx_tutorial_for_all_cells_index_end
  }

  // IJ example
  {
    // sphinx_tutorial_for_all_cells_ij_start

    const IndexType jp = mesh.cellJp();

    IndexType* ID = mesh.getFieldPtr<IndexType>("ID", mint::CELL_CENTERED);
    mint::for_all_cells<exec_policy, mint::xargs::ij>(
      &mesh,
      AXOM_LAMBDA(IndexType cellIdx, IndexType i, IndexType j) { ID[cellIdx] = i + j * jp; });

    // sphinx_tutorial_for_all_cells_ij_end
  }

  // nodeIds example
  {
    // sphinx_tutorial_for_all_cells_nodeids_start

    const double* vx = mesh.getFieldPtr<double>("vx", mint::NODE_CENTERED);
    const double* vy = mesh.getFieldPtr<double>("vy", mint::NODE_CENTERED);

    double* cell_vx = mesh.getFieldPtr<double>("cell_vx", mint::CELL_CENTERED);
    double* cell_vy = mesh.getFieldPtr<double>("cell_vy", mint::CELL_CENTERED);

    mint::for_all_cells<exec_policy, mint::xargs::nodeids>(
      &mesh,
      AXOM_LAMBDA(IndexType cellIdx, const IndexType* nodeIDs, IndexType N) {
        // sum nodal contributions
        cell_vx[cellIdx] = 0.0;
        cell_vy[cellIdx] = 0.0;
        for(IndexType inode = 0; inode < N; ++inode)
        {
          cell_vx[cellIdx] += vx[nodeIDs[inode]];
          cell_vy[cellIdx] += vy[nodeIDs[inode]];
        }  // END for all cell nodes

        // average at the cell center
        const double invf = 1.0 / static_cast<double>(N);
        cell_vx[cellIdx] *= invf;
        cell_vy[cellIdx] *= invf;
      });

    // sphinx_tutorial_for_all_cells_nodeids_end
  }

  // coords example
  {
    // sphinx_tutorial_for_all_cells_coords_start

    double* xc = mesh.getFieldPtr<double>("xc", mint::CELL_CENTERED);
    double* yc = mesh.getFieldPtr<double>("yc", mint::CELL_CENTERED);

    mint::for_all_cells<exec_policy, mint::xargs::coords>(
      &mesh,
      AXOM_LAMBDA(IndexType cellIdx,
                  const numerics::Matrix<double>& coords,
                  const IndexType* AXOM_UNUSED_PARAM(nodeIdx)) {
        // sum nodal coordinates
        double xsum = 0.0;
        double ysum = 0.0;
        const int numNodes = coords.getNumColumns();
        for(int inode = 0; inode < numNodes; ++inode)
        {
          const double* node = coords.getColumn(inode);
          xsum += node[mint::X_COORDINATE];
          ysum += node[mint::Y_COORDINATE];
        }  // end for all cell nodes

        // compute centroid by averaging nodal coordinates
        const double invf = 1.0 / static_cast<double>(numNodes);
        xc[cellIdx] = xsum * invf;
        yc[cellIdx] = ysum * invf;
      });

    // sphinx_tutorial_for_all_cells_coords_end
  }

  // faceids example
  {
    // sphinx_tutorial_for_all_cells_faceids_start

    const double* area = mesh.getFieldPtr<double>("area", mint::FACE_CENTERED);
    double* perimeter = mesh.getFieldPtr<double>("perimeter", mint::CELL_CENTERED);

    mint::for_all_cells<exec_policy, mint::xargs::faceids>(
      &mesh,
      AXOM_LAMBDA(IndexType cellIdx, const IndexType* faceIDs, IndexType N) {
        perimeter[cellIdx] = 0.0;
        for(IndexType iface = 0; iface < N; ++iface)
        {
          perimeter[cellIdx] += area[faceIDs[iface]];
        }
      });

    // sphinx_tutorial_for_all_cells_faceids_end
  }
}

//------------------------------------------------------------------------------
void face_traversals()
{
  using exec_policy = axom::SEQ_EXEC;
  using IndexType = axom::IndexType;

  const double lo[] = {-5.0, -5.0};
  const double hi[] = {5.0, 5.0};
  mint::UniformMesh mesh(lo, hi, 50, 50);

  mesh.createField<double>("vx", mint::NODE_CENTERED);
  mesh.createField<double>("vy", mint::NODE_CENTERED);

  mesh.createField<double>("t1", mint::FACE_CENTERED);
  mesh.createField<double>("t2", mint::FACE_CENTERED);
  mesh.createField<double>("w", mint::FACE_CENTERED);
  ;
  mesh.createField<double>("temp", mint::FACE_CENTERED);
  mesh.createField<double>("face_vx", mint::FACE_CENTERED);
  mesh.createField<double>("face_vy", mint::FACE_CENTERED);
  mesh.createField<double>("fx", mint::FACE_CENTERED);
  mesh.createField<double>("fy", mint::FACE_CENTERED);
  mesh.createField<IndexType>("boundary", mint::FACE_CENTERED);

  // Index example
  {
    // sphinx_tutorial_for_all_faces_index_start

    const double* t1 = mesh.getFieldPtr<double>("t1", mint::FACE_CENTERED);
    const double* t2 = mesh.getFieldPtr<double>("t2", mint::FACE_CENTERED);
    const double* w = mesh.getFieldPtr<double>("w", mint::FACE_CENTERED);

    double* temp = mesh.getFieldPtr<double>("temp", mint::FACE_CENTERED);
    mint::for_all_faces<exec_policy>(
      &mesh,
      AXOM_LAMBDA(IndexType faceIdx) {
        const double wf = w[faceIdx];
        const double a = t1[faceIdx];
        const double b = t2[faceIdx];

        temp[faceIdx] = wf * a + (1. - wf) * b;
      });

    // sphinx_tutorial_for_all_faces_index_end
  }

  // nodeids example
  {
    // sphinx_tutorial_for_all_faces_nodeids_start

    const double* vx = mesh.getFieldPtr<double>("vx", mint::NODE_CENTERED);
    const double* vy = mesh.getFieldPtr<double>("vy", mint::NODE_CENTERED);

    double* face_vx = mesh.getFieldPtr<double>("face_vx", mint::FACE_CENTERED);
    double* face_vy = mesh.getFieldPtr<double>("face_vy", mint::FACE_CENTERED);

    mint::for_all_faces<exec_policy, mint::xargs::nodeids>(
      &mesh,
      AXOM_LAMBDA(IndexType faceIdx, const IndexType* nodeIDs, IndexType N) {
        // sum constituent face node contributions
        face_vx[faceIdx] = 0.0;
        face_vy[faceIdx] = 0.0;
        for(int inode = 0; inode < N; ++inode)
        {
          face_vx[faceIdx] += vx[nodeIDs[inode]];
          face_vy[faceIdx] += vy[nodeIDs[inode]];
        }  // END for all face nodes

        // average
        const double invf = 1.0 / static_cast<double>(N);
        face_vx[faceIdx] *= invf;
        face_vy[faceIdx] *= invf;
      });

    // sphinx_tutorial_for_all_faces_nodeids_end
  }

  // coords example
  {
    // sphinx_tutorial_for_all_faces_coords_start

    double* fx = mesh.getFieldPtr<double>("fx", mint::FACE_CENTERED);
    double* fy = mesh.getFieldPtr<double>("fy", mint::FACE_CENTERED);

    mint::for_all_faces<exec_policy, mint::xargs::coords>(
      &mesh,
      AXOM_LAMBDA(IndexType faceIdx,
                  const numerics::Matrix<double>& coords,
                  const IndexType* AXOM_UNUSED_PARAM(nodeIdx)) {
        // sum nodal coordinates
        double xsum = 0.0;
        double ysum = 0.0;
        const int numNodes = coords.getNumColumns();
        for(int inode = 0; inode < numNodes; ++inode)
        {
          const double* node = coords.getColumn(inode);
          xsum += node[mint::X_COORDINATE];
          ysum += node[mint::Y_COORDINATE];
        }  // end for all face nodes

        // compute centroid by averaging nodal coordinates
        const double invf = 1.0 / static_cast<double>(numNodes);
        fx[faceIdx] = xsum * invf;
        fy[faceIdx] = ysum * invf;
      });

    // sphinx_tutorial_for_all_faces_coords_end
  }

  // cellids example
  {
    // sphinx_tutorial_for_all_faces_cellids_start

    constexpr IndexType ON_BOUNDARY = 1;
    constexpr IndexType INTERIOR = 0;

    IndexType* boundary = mesh.getFieldPtr<IndexType>("boundary", mint::FACE_CENTERED);

    mint::for_all_faces<exec_policy, mint::xargs::cellids>(
      &mesh,
      AXOM_LAMBDA(IndexType faceIdx, IndexType AXOM_UNUSED_PARAM(c1), IndexType c2) {
        boundary[faceIdx] = (c2 == -1) ? ON_BOUNDARY : INTERIOR;
      });

    // sphinx_tutorial_for_all_faces_cellids_end
  }
}

//------------------------------------------------------------------------------
void fem_tutorial()
{
  // sphinx_tutorial_create_fe_start

  constexpr bool ZERO_COPY = true;

  double coords[] = {
    0.0,
    0.0,  // x1, y1
    5.0,
    0.0,  // x2, y2
    5.0,
    5.0,  // x3, y3,
    0.0,
    5.0  // x4, y4
  };

  numerics::Matrix<double> nodes_matrix(2, 4, coords, ZERO_COPY);
  mint::FiniteElement fe(nodes_matrix, mint::QUAD);

  // bind to FE basis, wires the shape function pointers
  mint::bind_basis<MINT_LAGRANGE_BASIS, mint::QUAD>(fe);

  // sphinx_tutorial_create_fe_end

  // sphinx_tutorial_evaluate_shape_functions_start

  // isoparametric center
  double xi[] = {0.5, 0.5};
  double N[4];
  fe.evaluateShapeFunctions(xi, N);

  // sphinx_tutorial_evaluate_shape_functions_end

  // sphinx_tutorial_jacobian_start

  numerics::Matrix<double> J(2, 2);
  fe.jacobian(xi, J);

  const double jdet = numerics::determinant(J);
  std::cout << "jacobian determinant: " << jdet << std::endl;

  // sphinx_tutorial_jacobian_end

  // sphinx_tutorial_forward_map_start

  double xc[2];
  fe.computePhysicalCoords(xi, xc);
  std::cout << "xc: ( ";
  std::cout << xc[0] << ", " << xc[1];
  std::cout << " )\n";

  // sphinx_tutorial_forward_map_end

  // sphinx_tutorial_inverse_map_start

  double xr[2];
  int status = fe.computeReferenceCoords(xc, xr);

  switch(status)
  {
  case mint::INVERSE_MAP_FAILED:
    std::cout << "Newton-Raphson failed!";
    break;
  case mint::OUTSIDE_ELEMENT:
    std::cout << "point is outside!\n";
    break;
  default:
    // found the reference coordinates!
    std::cout << "xr: ( ";
    std::cout << xr[0] << ", " << xr[1];
    std::cout << " )\n";
  }

  // sphinx_tutorial_inverse_map_end
}

//------------------------------------------------------------------------------
void working_with_fields(mint::Mesh* mesh)
{
  constexpr int FIELD_ASSOCIATION = mint::NODE_CENTERED;

  {  // START unnamed namespace

    // sphinx_tutorial_add_fields_start

    double* den = mesh->createField<double>("den", mint::CELL_CENTERED);
    double* vel = mesh->createField<double>("vel", mint::NODE_CENTERED, 3);

    // sphinx_tutorial_add_fields_end

    AXOM_UNUSED_VAR(den);
    AXOM_UNUSED_VAR(vel);

  }  // END unnamed namespace

  // sphinx_tutorial_get_field_by_name_start

  double* den = mesh->getFieldPtr<double>("den", mint::CELL_CENTERED);

  axom::IndexType nc = -1;  // number of components
  double* vel = mesh->getFieldPtr<double>("vel", mint::NODE_CENTERED, nc);

  // sphinx_tutorial_get_field_by_name_end

  // sphinx_tutorial_check_fields_start

  const bool hasDen = mesh->hasField("den", mint::CELL_CENTERED);
  const bool hasVel = mesh->hasField("vel", mint::NODE_CENTERED);

  // sphinx_tutorial_check_fields_end

  // sphinx_tutorial_remove_fields_start

  bool isRemoved = mesh->removeField("den", mint::CELL_CENTERED);

  // sphinx_tutorial_remove_fields_end

  // sphinx_tutorial_query_fields_start

  const mint::FieldData* fieldData = mesh->getFieldData(FIELD_ASSOCIATION);

  const int numFields = fieldData->getNumFields();
  for(int ifield = 0; ifield < numFields; ++ifield)
  {
    const mint::Field* field = fieldData->getField(ifield);

    const std::string& fieldName = field->getName();
    axom::IndexType numTuples = field->getNumTuples();
    axom::IndexType numComponents = field->getNumComponents();

    std::cout << "field name: " << fieldName << std::endl;
    std::cout << "numTuples: " << numTuples << std::endl;
    std::cout << "numComponents: " << numComponents << std::endl;

    if(field->getType() == mint::DOUBLE_FIELD_TYPE)
    {
      double* data = mesh->getFieldPtr<double>(fieldName, FIELD_ASSOCIATION);
      data[0] = 42.0;
      // process double precision floating point data
      // ...
    }
    else if(field->getType() == mint::INT32_FIELD_TYPE)
    {
      int* data = mesh->getFieldPtr<int>(fieldName, FIELD_ASSOCIATION);
      data[0] = 42;
      // process integral data
      // ...
    }
    // ...

  }  // END for all fields

  // sphinx_tutorial_query_fields_end

  AXOM_UNUSED_VAR(hasDen);
  AXOM_UNUSED_VAR(hasVel);
  AXOM_UNUSED_VAR(den);
  AXOM_UNUSED_VAR(vel);
  AXOM_UNUSED_VAR(isRemoved);
}

//------------------------------------------------------------------------------
void construct_uniform()
{
  // sphinx_tutorial_construct_uniform_start

  const double lo[] = {-5.0, -5.0};
  const double hi[] = {5.0, 5.0};
  mint::UniformMesh mesh(lo, hi, 50, 50);

  // sphinx_tutorial_construct_uniform_end

  mint::write_vtk(&mesh, "tutorial_uniform_mesh.vtk");
}

//------------------------------------------------------------------------------
void construct_rectilinear()
{
  // sphinx_tutorial_construct_rectilinear_start

  constexpr double beta = 0.1;
  const double expbeta = exp(beta);
  const double invf = 1 / (expbeta - 1.0);

  // construct a N x N rectilinear mesh
  constexpr axom::IndexType N = 25;
  mint::RectilinearMesh mesh(N, N);
  double* x = mesh.getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh.getCoordinateArray(mint::Y_COORDINATE);

  // fill the coordinates along each axis
  x[0] = y[0] = 0.0;
  for(int i = 1; i < N; ++i)
  {
    const double delta = (exp(i * beta) - 1.0) * invf;
    x[i] = x[i - 1] + delta;
    y[i] = y[i - 1] + delta;
  }

  // sphinx_tutorial_construct_rectilinear_end

  mint::write_vtk(&mesh, "tutorial_rectilinear_mesh.vtk");
}

//------------------------------------------------------------------------------
void construct_curvilinear()
{
  // sphinx_tutorial_construct_curvilinear_start

  constexpr double R = 2.5;
  constexpr double M = (2 * M_PI) / 50.0;
  constexpr double h = 0.5;
  constexpr axom::IndexType N = 25;

  // construct the curvilinear mesh object
  mint::CurvilinearMesh mesh(N, N);

  // get a handle on the coordinate arrays
  double* x = mesh.getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh.getCoordinateArray(mint::Y_COORDINATE);

  // fill the coordinates of the curvilinear mesh
  const axom::IndexType jp = mesh.nodeJp();
  for(axom::IndexType j = 0; j < N; ++j)
  {
    const axom::IndexType j_offset = j * jp;

    for(axom::IndexType i = 0; i < N; ++i)
    {
      const axom::IndexType nodeIdx = i + j_offset;

      const double xx = h * i;
      const double yy = h * j;
      const double alpha = yy + R;
      const double beta = xx * M;

      x[nodeIdx] = alpha * cos(beta);
      y[nodeIdx] = alpha * sin(beta);
    }  // END for all i

  }  // END for all j

  // sphinx_tutorial_construct_curvilinear_end

  mint::write_vtk(&mesh, "tutorial_curvilinear_mesh.vtk");
}

//------------------------------------------------------------------------------
void construct_single_cell_type_unstructured()
{
  // sphinx_tutorial_construct_unstructured_start

  constexpr int DIMENSION = 2;
  constexpr mint::CellType CELL_TYPE = mint::TRIANGLE;

  // Construct the mesh object
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(DIMENSION, CELL_TYPE);

  // Append the mesh nodes
  const axom::IndexType n0 = mesh.appendNode(0.0, 0.0);
  const axom::IndexType n1 = mesh.appendNode(2.0, 0.0);
  const axom::IndexType n2 = mesh.appendNode(1.0, 1.0);
  const axom::IndexType n3 = mesh.appendNode(3.5, 1.0);
  const axom::IndexType n4 = mesh.appendNode(2.5, 2.0);
  const axom::IndexType n5 = mesh.appendNode(5.0, 0.0);

  // Append mesh cells
  const axom::IndexType c0[] = {n1, n3, n2};
  const axom::IndexType c1[] = {n2, n0, n1};
  const axom::IndexType c2[] = {n3, n4, n2};
  const axom::IndexType c3[] = {n1, n5, n3};

  mesh.appendCell(c0);
  mesh.appendCell(c1);
  mesh.appendCell(c2);
  mesh.appendCell(c3);

  // sphinx_tutorial_construct_unstructured_end

  // STEP 3: write output mesh
  mint::write_vtk(&mesh, "tutorial_triangle_mesh.vtk");
}

//------------------------------------------------------------------------------
void construct_mixed_cell_type_unstructured()
{
  // sphinx_tutorial_construct_unstructured_mixed_start

  constexpr int DIMENSION = 2;

  // Construct the mesh object
  mint::UnstructuredMesh<mint::MIXED_SHAPE> mesh(DIMENSION);

  // Append the mesh nodes
  const axom::IndexType n0 = mesh.appendNode(0.0, 0.0);
  const axom::IndexType n1 = mesh.appendNode(2.0, 0.0);
  const axom::IndexType n2 = mesh.appendNode(1.0, 1.0);
  const axom::IndexType n3 = mesh.appendNode(3.5, 1.0);
  const axom::IndexType n4 = mesh.appendNode(2.5, 2.0);
  const axom::IndexType n5 = mesh.appendNode(5.0, 0.0);

  // Append mesh cells
  const axom::IndexType c0[] = {n0, n1, n2};
  const axom::IndexType c1[] = {n1, n5, n3, n2};
  const axom::IndexType c2[] = {n3, n4, n2};

  mesh.appendCell(c0, mint::TRIANGLE);
  mesh.appendCell(c1, mint::QUAD);
  mesh.appendCell(c2, mint::TRIANGLE);

  // sphinx_tutorial_construct_unstructured_mixed_end

  // STEP 3: write output mesh
  mint::write_vtk(&mesh, "tutorial_mixed_unstructured_mesh.vtk");
}

//------------------------------------------------------------------------------
void using_external_storage()
{
  // sphinx_tutorial_using_external_storage_start

  constexpr axom::IndexType NUM_NODES = 6;
  constexpr axom::IndexType NUM_CELLS = 4;

  // application buffers
  double x[] = {0.0, 2.0, 1.0, 3.5, 2.5, 5.0};
  double y[] = {0.0, 0.0, 1.0, 1.0, 2.0, 0.0};

  axom::IndexType cell_connectivity[] = {
    1,
    3,
    2,  // c0
    2,
    0,
    1,  // c1
    3,
    4,
    2,  // c2
    1,
    5,
    3  // c3
  };

  // cell-centered density field values
  double den[] = {0.5, 1.2, 2.5, 0.9};

  // construct mesh object with external buffers
  using MeshType = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  MeshType* mesh = new MeshType(mint::TRIANGLE, NUM_CELLS, cell_connectivity, NUM_NODES, x, y);

  // register external field
  mesh->createField<double>("den", mint::CELL_CENTERED, den);

  // output external mesh to vtk
  mint::write_vtk(mesh, "tutorial_external_mesh.vtk");

  // delete the mesh, doesn't touch application buffers
  delete mesh;
  mesh = nullptr;

  // sphinx_tutorial_using_external_storage_end
}

//------------------------------------------------------------------------------
void using_sidre()
{
#ifdef AXOM_MINT_USE_SIDRE

  // sphinx_tutorial_using_sidre_create_mesh_start
  // create a Sidre Datastore to store the mesh
  sidre::DataStore ds;
  sidre::Group* group = ds.getRoot();

  // Construct the mesh object and populate the supplied Sidre Group
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(2, mint::TRIANGLE, group);

  // Append the mesh nodes
  const axom::IndexType n0 = mesh.appendNode(0.0, 0.0);
  const axom::IndexType n1 = mesh.appendNode(2.0, 0.0);
  const axom::IndexType n2 = mesh.appendNode(1.0, 1.0);
  const axom::IndexType n3 = mesh.appendNode(3.5, 1.0);
  const axom::IndexType n4 = mesh.appendNode(2.5, 2.0);
  const axom::IndexType n5 = mesh.appendNode(5.0, 0.0);

  // Append mesh cells
  const axom::IndexType c0[] = {n1, n3, n2};
  const axom::IndexType c1[] = {n2, n0, n1};
  const axom::IndexType c2[] = {n3, n4, n2};
  const axom::IndexType c3[] = {n1, n5, n3};

  mesh.appendCell(c0);
  mesh.appendCell(c1);
  mesh.appendCell(c2);
  mesh.appendCell(c3);

  // create a cell-centered field
  double* den = mesh.createField<double>("den", mint::CELL_CENTERED);

  // set density values at each cell
  den[0] = 0.5;  // c0
  den[1] = 1.2;  // c1
  den[2] = 2.5;  // c2
  den[3] = 0.9;  // c3

  // sphinx_tutorial_using_sidre_create_mesh_end

  std::ofstream ofs;
  ofs.open("sphinx_tutorial_sidre_tree.txt");
  group->print(ofs);
  ofs.close();

  // sphinx_tutorial_sidre_import_start

  mint::Mesh* imported_mesh = mint::getMesh(group);
  std::cout << "Mesh Type: " << imported_mesh->getMeshType() << std::endl;
  std::cout << "hasSidre: " << imported_mesh->hasSidreGroup() << std::endl;

  mint::write_vtk(imported_mesh, "tutorial_imported_mesh.vtk");

  delete imported_mesh;
  imported_mesh = nullptr;

  // sphinx_tutorial_sidre_import_end

#endif /* AXOM_MINT_USE_SIDRE */
}

/*!
 * \brief Tutorial main
 */
int main(int AXOM_UNUSED_PARAM(argc), char** AXOM_UNUSED_PARAM(argv))
{
  // Native construction of various mesh types
  construct_uniform();
  construct_rectilinear();
  construct_curvilinear();
  construct_single_cell_type_unstructured();
  construct_mixed_cell_type_unstructured();

  fem_tutorial();
  using_sidre();
  using_external_storage();

  node_traversals();
  cell_traversals();
  face_traversals();

  const double lo[] = {-5.0, -5.0};
  const double hi[] = {5.0, 5.0};
  mint::UniformMesh mesh(lo, hi, 50, 50);

  working_with_fields(&mesh);
  vtk_output(&mesh, "tutorial_mesh.vtk");

  return 0;
}
