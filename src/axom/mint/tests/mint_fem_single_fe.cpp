// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic.hpp"

#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/Determinants.hpp"

// Mint includes
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/config.hpp"
#include "axom/mint/fem/FEBasis.hpp"
#include "axom/mint/mesh/Field.hpp"
#include "axom/mint/mesh/FieldData.hpp"
#include "axom/mint/mesh/FieldVariable.hpp"
#include "axom/mint/fem/FiniteElement.hpp"
#include "axom/mint/fem/shape_functions/Lagrange.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/fem/shape_functions/ShapeFunction.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/utils/vtk_utils.hpp"

// C/C++ includes
#include <cmath>  // for sin(), cos(), sqrt(), etc.

//#define MINT_FEM_DEBUG

using namespace axom;

//------------------------------------------------------------------------------
//  INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
/*!
 * \brief Helper method to calculate the centroid of Finite Element instance.
 * \param [in] fe pointer to a finite element instance.
 * \param [out] centroid buffer to store the centroid
 */
void compute_centroid(mint::FiniteElement* fe, double* centroid)
{
  EXPECT_TRUE(fe != nullptr);
  EXPECT_TRUE(centroid != nullptr);

  const int ndims = fe->getPhysicalDimension();
  const int nnodes = fe->getNumNodes();
  numerics::Matrix<double> physical_nodes(ndims, nnodes, fe->getPhysicalNodes(), true);

  if(fe->getCellType() == mint::PYRAMID)
  {
    // Hard-code expected centroid from VTK, since averaging the nodes
    // doesn't yield the centroid for a pyramid.
    //
    // NOTE: these values are hard-coded for the specific pyramid in this test.

    centroid[0] = 4.00001;
    centroid[1] = 4.28284;
    centroid[2] = 2.0;
  }
  else
  {
    for(int i = 0; i < ndims; ++i)
    {
      IndexType p = 0;
      IndexType N = 0;
      const double* xp_i = physical_nodes.getRow(i, p, N);

      centroid[i] = 0.0;
      for(int j = 0; j < N; j += p)
      {
        centroid[i] += xp_i[j];
      }  // END for all nodes

      centroid[i] /= nnodes;

    }  // END for all dimensions

  }  // END else
}

/*!
 * \brief Helper method to construct a single element of the specified type.
 *
 * \note Constructs the mesh element by scaling and rotating the reference
 *  element.
 *
 * \note Ownership of the returned pointer is propagated to the caller. The
 *  calling method is therefore responsible for managing and properly
 *  deallocating the returned mesh object.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CELLTYPE the corresponding cell type, e.g., MINT_QUAD
 */
template <int BasisType, mint::CellType CELLTYPE>
void get_single_fe(mint::FiniteElement*& fe)
{
  EXPECT_TRUE(fe == nullptr);

  using FEMType = typename mint::FEBasis<BasisType, CELLTYPE>;
  using ShapeFunctionType = typename FEMType::ShapeFunctionType;

  const bool zero_copy = true;
  const double SCALE = 10.0;
  const double ANGLE = 0.785398;  // 45 degrees
  const double SINT = sin(ANGLE);
  const double COST = cos(ANGLE);

  const int ndims = ShapeFunctionType::dimension();
  EXPECT_TRUE(ndims == 2 || ndims == 3);

  const int ndofs = ShapeFunctionType::numDofs();

  axom::IndexType* cell = new axom::IndexType[ndofs];

  double* center = new double[ndims];
  ShapeFunctionType::center(center);

  double* nodes = new double[ndofs * ndims];
  ShapeFunctionType::coords(nodes);

  double centroid[] = {0.0, 0.0, 0.0};  // used to compute the pyramid apex

  for(int i = 0; i < ndofs; ++i)
  {
    double* node = &nodes[i * ndims];

    // scale & rotate cell from the reference space to get a  test cell
    const double dx = node[0] - center[0];
    const double dy = node[1] - center[1];

    node[0] = SCALE * ((dx * COST - dy * SINT) + center[0]);
    node[1] = SCALE * ((dx * SINT + dy * COST) + center[1]);
    if(ndims == 3)
    {
      node[2] *= SCALE;
    }

    if(CELLTYPE == mint::PYRAMID && i < 4)
    {
      centroid[0] += node[0];
      centroid[1] += node[1];
      centroid[2] += node[2];
    }

    if(CELLTYPE == mint::PYRAMID && i == 4)
    {
      // generate right pyramid, ensure the apex is prependicular to the
      // base of the pyramid to facilitate testing.
      node[0] = 0.25 * centroid[0];
      node[1] = 0.25 * centroid[1];
      node[2] = 0.25 * centroid[2] + SCALE;
    }
  }

  numerics::Matrix<double> m(ndims, ndofs, nodes, zero_copy);
  fe = new mint::FiniteElement(m, CELLTYPE);
  mint::bind_basis<BasisType, CELLTYPE>(*fe);
  EXPECT_FALSE(fe->getBasisType() == MINT_UNDEFINED_BASIS);

#ifdef MINT_FEM_DEBUG
  std::string vtkFile = std::string(mint::getCellInfo(CELLTYPE).name) + ".vtk";
  mint::write_vtk(*fe, vtkFile);
#endif

  // clean up
  delete[] center;
  delete[] nodes;
  delete[] cell;
}

/*!
 * \brief Ensures the FiniteElement instance is properly constructed by
 *  checking the underlying ShapeFunction.
 *
 * \param [in] fe pointer to the finite element instance to test
 *
 * \pre fe != nullptr
 */
template <typename ShapeFunctionType>
void check_reference_element(mint::FiniteElement* fe)
{
  EXPECT_TRUE(fe != nullptr);
  EXPECT_FALSE(fe->getBasisType() == MINT_UNDEFINED_BASIS);

  // STEP 0: check reference element attributes, eg, dofs, dimension, etc.
  ShapeFunctionType sf;
  EXPECT_EQ(sf.numDofs(), fe->getNumDofs());
  EXPECT_EQ(sf.dimension(), fe->getReferenceDimension());

  const int fe_dim = fe->getPhysicalDimension();
  EXPECT_TRUE((fe_dim >= 1) && (fe_dim <= 3));

  const int ref_dim = fe->getReferenceDimension();
  EXPECT_TRUE(ref_dim <= fe_dim);

  // STEP 1: check reference coordinates
  const int ndofs = sf.numDofs();
  const int ndims = sf.dimension();

  // populate a matrix with shape function coordinates (expected)
  numerics::Matrix<double> shape_coords(ndims, ndofs);
  sf.coords(shape_coords.data());

  // populate a matrix with the element's reference coordinates to test
  numerics::Matrix<double> ref_coords(ndims, ndofs, fe->getReferenceNodes(), true);

  for(int i = 0; i < ndofs; ++i)
  {
    const double* xr_expected = shape_coords.getColumn(i);
    const double* xr = ref_coords.getColumn(i);

    for(int j = 0; j < ndims; ++j)
    {
      EXPECT_DOUBLE_EQ(xr_expected[j], xr[j]);
    }

  }  // END for all dofs

  // STEP 2: check reference center
  double* shape_center = new double[ndims];
  sf.center(shape_center);

  for(int j = 0; j < ref_dim; ++j)
  {
    const double expected = shape_center[j];
    const double actual = fe->getReferenceCenter()[j];
    EXPECT_DOUBLE_EQ(expected, actual);
  }

  // STEP 3: clean up
  delete[] shape_center;
}

/*!
 * \brief Checks the mapping from reference space to physical space.
 *
 *  The test loops over all the nodes of the reference element, maps them to
 *  physical space and makes sure they match with the corresponding physical
 *  coordinates of the element.
 *
 * \param [in] fe pointer to the finite element instance to test
 *
 * \pre fe != nullptr
 *
 * \see check_forward_map
 */
void test_forward_map(mint::FiniteElement* fe, double TOL = 1.e-9)
{
  EXPECT_TRUE(fe != nullptr);
  EXPECT_FALSE(fe->getBasisType() == MINT_UNDEFINED_BASIS);
  EXPECT_TRUE(fe->getReferenceDimension() == fe->getPhysicalDimension());
  EXPECT_TRUE(fe->getNumNodes() == fe->getNumDofs());

  int numdofs = fe->getNumDofs();
  int nnodes = fe->getNumNodes();
  int ndims = fe->getReferenceDimension();

  // STEP 0: get the reference coordinates
  numerics::Matrix<double> reference_coords(ndims, numdofs, fe->getReferenceNodes(), true);

  // STEP 1: get the physical coordinates of the element
  numerics::Matrix<double> element_coords(ndims, nnodes, fe->getPhysicalNodes(), true);

  // STEP 2: buffer to store the mapped physical coordinates
  numerics::Matrix<double> physical_coords(ndims, nnodes);

  // STEP 3: map each reference coordinate to physical
  for(int i = 0; i < numdofs; ++i)
  {
    // forward map
    const double* xr = reference_coords.getColumn(i);
    double* xp = physical_coords.getColumn(i);
    fe->computePhysicalCoords(xr, xp);

    double* expected_xp = element_coords.getColumn(i);

    // check mapping
    for(int j = 0; j < ndims; ++j)
    {
      EXPECT_NEAR(expected_xp[j], xp[j], TOL);
    }

  }  // END for all dofs

  // STEP 3: map center
  const double* xr_center = fe->getReferenceCenter();
  double* xp_center = new double[ndims];
  fe->computePhysicalCoords(xr_center, xp_center);

  // STEP 4: compute expected centroid by averaging coordinates
  double* expected_centroid = new double[ndims];
  compute_centroid(fe, expected_centroid);

  // STEP 5: check mapping of centroid
  for(int i = 0; i < ndims; ++i)
  {
    EXPECT_NEAR(expected_centroid[i], xp_center[i], TOL);
  }

  // STEP 6: cleanup
  delete[] xp_center;
  delete[] expected_centroid;
}

/*!
 * \brief Checks the mapping from physical space to reference space.
 *
 *  The test loops over all the nodes of the element in physical space, maps
 *  them to the reference space and makes sure they match with the corresponding
 *  reference coordinates of the element.
 *
 * \param [in] fe pointer to the finite element instance to test
 *
 * \pre fe != nullptr
 *
 * \see check_inverse_map
 */
void test_inverse_map(mint::FiniteElement* fe, double TOL = 1.e-9)
{
  EXPECT_TRUE(fe != nullptr);
  EXPECT_FALSE(fe->getBasisType() == MINT_UNDEFINED_BASIS);
  EXPECT_TRUE(fe->getReferenceDimension() == fe->getPhysicalDimension());
  EXPECT_TRUE(fe->getNumNodes() == fe->getNumDofs());

  // STEP 0: get Matrix of physical nodes; nodal coordinates stored in columns.
  const int ndims = fe->getPhysicalDimension();
  const int nnodes = fe->getNumNodes();
  numerics::Matrix<double> physical_nodes(ndims, nnodes, fe->getPhysicalNodes(), true);

  // STEP 1: get reference coordinates
  numerics::Matrix<double> reference_nodes(ndims, nnodes, fe->getReferenceNodes(), true);

  // STEP 2: allocate buffer to store the computed reference coordinates
  double* xr = new double[ndims];

  // STEP 3: loop over physical nodes and map them to reference space
  for(int i = 0; i < nnodes; ++i)
  {
    if(fe->getCellType() == mint::PYRAMID && i == 4)
    {
      // skip inverse map at the apex of the pyramid, system is singular
      continue;
    }

    const double* expected_xr = reference_nodes.getColumn(i);
    const double* xp = physical_nodes.getColumn(i);
    int rc = fe->computeReferenceCoords(xp, xr, TOL);
    EXPECT_TRUE(rc == mint::INSIDE_ELEMENT);

    // check mapping
    for(int j = 0; j < ndims; ++j)
    {
      EXPECT_NEAR(expected_xr[j], xr[j], TOL);
    }

  }  // END for all physical nodes

  // STEP 4: calculate centroid in physical space by averaging coordinates
  double* centroid = new double[ndims];
  compute_centroid(fe, centroid);

  // STEP 4: map centroid
  int rc = fe->computeReferenceCoords(centroid, xr, TOL);
  EXPECT_TRUE(rc == mint::INSIDE_ELEMENT);

  // STEP 5: check centroid mapping
  const double* refcenter = fe->getReferenceCenter();
  for(int i = 0; i < ndims; ++i)
  {
    EXPECT_NEAR(refcenter[i], xr[i], TOL);
  }

  // STEP 6: clean up
  delete[] xr;
  delete[] centroid;
}

/*!
 * \brief Performs basic checks on a FiniteElement instance bound to a
 *  corresponding basis function.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 *
 * \see check_reference_element()
 */
template <int BasisType, mint::CellType CELLTYPE>
void check_shape()
{
  const int cell_value = mint::cellTypeToInt(CELLTYPE);
  EXPECT_TRUE((cell_value >= 0) && (cell_value < mint::NUM_CELL_TYPES));
  EXPECT_TRUE((BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES));

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  using FEMType = typename mint::FEBasis<BasisType, CELLTYPE>;
  using ShapeFunctionType = typename FEMType::ShapeFunctionType;

  // STEP 0: construct finite element mesh
  mint::FiniteElement* fe = nullptr;
  get_single_fe<BasisType, CELLTYPE>(fe);

  EXPECT_TRUE(fe != nullptr);
  EXPECT_EQ(mint::getCellInfo(CELLTYPE).num_nodes, fe->getNumNodes());

  // STEP 1: test FE instance
  check_reference_element<ShapeFunctionType>(fe);

  // STEP 2: clean up
  delete fe;
}

/*!
 * \brief Ensures the jacobian is positive within the element by evaluating
 *  the determinant of the jacobian at the element nodes, center and computed
 *  interior points.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CELLTYPE the corresponding cell type, e.g., MINT_QUAD
 */
template <int BasisType, mint::CellType CELLTYPE>
void check_jacobian(double TOL = 1.e-9)
{
  const int cell_value = mint::cellTypeToInt(CELLTYPE);
  EXPECT_TRUE((cell_value >= 0) && (cell_value < mint::NUM_CELL_TYPES));
  EXPECT_TRUE((BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES));

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  const double LTOL = 0.0 - TOL;
  double det = 0.0;

  // STEP 0: construct a single element
  mint::FiniteElement* fe = nullptr;
  get_single_fe<BasisType, CELLTYPE>(fe);

  EXPECT_TRUE(fe != nullptr);
  EXPECT_EQ(mint::getCellInfo(CELLTYPE).num_nodes, fe->getNumNodes());

  // STEP 1: construct a Matrix object to store the jacobian
  const int ndims = fe->getPhysicalDimension();
  numerics::Matrix<double> J(ndims, ndims);

  // STEP 2: test jacobian at the reference nodes
  const int ndofs = fe->getNumDofs();
  const int rdim = fe->getReferenceDimension();
  numerics::Matrix<double> n(rdim, ndofs, fe->getReferenceNodes(), true);
  for(int i = 0; i < ndofs; ++i)
  {
    const double* xi = n.getColumn(i);
    fe->jacobian(xi, J);

    det = numerics::determinant(J);
    EXPECT_GT(det, LTOL);
  }

  // STEP 3: test jacobian at the reference center
  const double* xi_c = fe->getReferenceCenter();
  fe->jacobian(xi_c, J);
  det = numerics::determinant(J);
  EXPECT_GT(det, LTOL);

  // STEP 4: test jacobian at interior points; the interior points are
  // computed by taking the midpoint of a reference node and the centroid
  double* rp = new double[ndims];
  for(int i = 0; i < ndofs; ++i)
  {
    const double* xi = n.getColumn(i);
    for(int j = 0; j < ndims; ++j)
    {
      rp[j] = 0.5 * (xi[j] + xi_c[j]);
    }

    fe->jacobian(rp, J);
    det = numerics::determinant(J);
    EXPECT_GT(det, LTOL);
  }

  // STEP 5: clean up
  delete[] rp;
  delete fe;
}

/*!
 * \brief Checks the forward map of a FiniteElement object.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CELLTYPE the corresponding cell type, e.g., MINT_QUAD
 *
 * \see test_forward_map()
 */
template <int BasisType, mint::CellType CELLTYPE>
void check_forward_map(double TOL = 1.e-9)
{
  const int cell_value = mint::cellTypeToInt(CELLTYPE);
  EXPECT_TRUE((cell_value >= 0) && (cell_value < mint::NUM_CELL_TYPES));
  EXPECT_TRUE((BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES));

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  // STEP 0: construct a mesh with a single element
  mint::FiniteElement* fe = nullptr;

  get_single_fe<BasisType, CELLTYPE>(fe);
  EXPECT_TRUE(fe != nullptr);
  EXPECT_EQ(mint::getCellInfo(CELLTYPE).num_nodes, fe->getNumNodes());

  // STEP 1: check forward mapping
  test_forward_map(fe, TOL);

  // STEP 2: clean up
  delete fe;
}

/*!
 * \brief Checks the inverse map of a FiniteElement object.
 *
 * \param [in] TOL optional user-supplied tolerange. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CELLTYPE the corresponding cell type, e.g., MINT_QUAD
 *
 * \see test_inverse_map()
 */
template <int BasisType, mint::CellType CELLTYPE>
void check_inverse_map(double TOL = 1.e-9)
{
  const int cell_value = mint::cellTypeToInt(CELLTYPE);
  EXPECT_TRUE((cell_value >= 0) && (cell_value < mint::NUM_CELL_TYPES));
  EXPECT_TRUE((BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES));

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  // STEP 0: construct a mesh with a single element
  mint::FiniteElement* fe = nullptr;

  get_single_fe<BasisType, CELLTYPE>(fe);
  EXPECT_TRUE(fe != nullptr);
  EXPECT_EQ(mint::getCellInfo(CELLTYPE).num_nodes, fe->getNumNodes());

  // STEP 1: check inverse map
  test_inverse_map(fe, TOL);

  // STEP 2: clean up
  delete fe;
}

/*!
 * \brief Checks correctness of the inverse map for point-in-cell queries.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CELLTYPE the corresponding cell type, e.g., MINT_QUAD
 */
template <int BasisType, mint::CellType CELLTYPE>
void point_in_cell(double TOL = 1.e-9)
{
  const int cell_value = mint::cellTypeToInt(CELLTYPE);
  EXPECT_TRUE((cell_value >= 0) && (cell_value < mint::NUM_CELL_TYPES));
  EXPECT_TRUE((BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES));

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  // STEP 0: construct a mesh with a single element
  mint::FiniteElement* fe = nullptr;

  get_single_fe<BasisType, CELLTYPE>(fe);
  EXPECT_TRUE(fe != nullptr);
  EXPECT_EQ(mint::getCellInfo(CELLTYPE).num_nodes, fe->getNumNodes());

  // STEP 0: test variables
  const int nnodes = fe->getNumNodes();
  const int ndims = fe->getPhysicalDimension();
  double* xi = new double[ndims];
  double* xc = new double[ndims];
  double* rp = new double[ndims];
  int status = 0;

  numerics::Matrix<double> nodes(ndims, nnodes, fe->getPhysicalNodes(), true);

  // STEP 1: Ensure element nodes are inside
  for(int i = 0; i < nnodes; ++i)
  {
    if(fe->getCellType() == mint::PYRAMID && i == 4)
    {
      // skip inverse map at the apex of the pyramid, system is singular
      continue;
    }

    const double* xp = nodes.getColumn(i);

    status = fe->computeReferenceCoords(xp, xi, TOL);
    EXPECT_TRUE(status == mint::INSIDE_ELEMENT);
  }

  // STEP 2: Ensure center is inside
  fe->computePhysicalCoords(fe->getReferenceCenter(), xc);

  status = fe->computeReferenceCoords(xc, xi, TOL);
  EXPECT_TRUE(status == mint::INSIDE_ELEMENT);

  // STEP 3: Ensure other interior nodes are inside
  numerics::Matrix<double> dofs(ndims, nnodes, fe->getReferenceNodes(), true);

  for(int i = 0; i < nnodes; ++i)
  {
    const double* dof = dofs.getColumn(i);

    for(int j = 0; j < ndims; ++j)
    {
      rp[j] = 0.5 * (dof[j] + fe->getReferenceCenter()[j]);
    }

    fe->computePhysicalCoords(rp, xc);

    status = fe->computeReferenceCoords(xc, xi, TOL);
    EXPECT_TRUE(status == mint::INSIDE_ELEMENT);

  }  // END for

  // STEP 4: test outside points by shifting reference element
  const double SHIFT = 10;
  for(int i = 0; i < nnodes; ++i)
  {
    const double* dof = dofs.getColumn(i);

    // Test +SHIFT
    for(int j = 0; j < ndims; ++j)
    {
      rp[j] = dof[j] + SHIFT;
    }

    fe->computePhysicalCoords(rp, xc);

    status = fe->computeReferenceCoords(xc, xi, TOL);
    EXPECT_TRUE(status == mint::OUTSIDE_ELEMENT);

    // Test -SHIFT
    for(int j = 0; j < ndims; ++j)
    {
      rp[j] = dof[j] - SHIFT;
    }

    fe->computePhysicalCoords(rp, xc);

    status = fe->computeReferenceCoords(xc, xi, TOL);
    EXPECT_TRUE(status == mint::OUTSIDE_ELEMENT);
  }

  // STEP 5: clean up
  delete[] xi;
  delete[] xc;
  delete[] rp;
  delete fe;
}

/*!
 * \brief Computes a fake node-centered field, using the nodal coordinates, to
 *  test interpolation within the element.
 *
 * \param [in] x pointer to the coordinates of the node.
 * \param [in] N the number of dimensions
 *
 * \return f the computed field.
 *
 * \see check_interp()
 */
double analytic_function(const double* x, int N)
{
  double f = 0.0;
  for(int i = 0; i < N; ++i)
  {
    f += x[i];
  }
  return (f);
}

/*!
 * \brief Given field values and associated weights at the nodes of an element,
 *  this method interpolates the field.
 *
 * \param [in] f user-supplied buffer where the nodal field values are stored
 * \param [in] wgts user-supplied buffer of corresponding interpolation weigths
 * \param [in] N the number of nodes, i.e., the length of f and wgts
 *
 * \pre f != nullptr
 * \pre wgts != nullptr
 *
 * \return finterp the interpolated field value
 *
 * \see check_interp()
 */
double interp(const double* f, const double* wgts, int N)
{
  EXPECT_TRUE(f != nullptr);
  EXPECT_TRUE(wgts != nullptr);

  double finterp = 0.0;
  for(int i = 0; i < N; ++i)
  {
    finterp += wgts[i] * f[i];
  }
  return (finterp);
}

/*!
 * \brief Checks the interpolation within the element against a field that
 *  is defined analytically.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CELLTYPE the corresponding cell type, e.g., MINT_QUAD
 */
template <int BasisType, mint::CellType CELLTYPE>
void check_interp(double TOL = 1.e-9)
{
  const int cell_value = mint::cellTypeToInt(CELLTYPE);
  EXPECT_TRUE((cell_value >= 0) && (cell_value < mint::NUM_CELL_TYPES));
  EXPECT_TRUE((BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES));

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  // STEP 0: construct a mesh with a single element
  mint::FiniteElement* fe = nullptr;

  get_single_fe<BasisType, CELLTYPE>(fe);
  EXPECT_TRUE(fe != nullptr);
  EXPECT_EQ(mint::getCellInfo(CELLTYPE).num_nodes, fe->getNumNodes());

  const int ndims = fe->getPhysicalDimension();
  const int nnodes = fe->getNumNodes();
  double* wgts = new double[nnodes];

  // STEP 1: setup a nodal field to interpolate
  double* f = new double[nnodes];

  numerics::Matrix<double> nodes(ndims, nnodes, fe->getPhysicalNodes(), true);

  for(int i = 0; i < nnodes; ++i)
  {
    const double* x = nodes.getColumn(i);
    f[i] = analytic_function(x, ndims);
  }

  // STEP 2: check interpolation at the nodes
  const int ndofs = fe->getNumDofs();
  numerics::Matrix<double> dofs(ndims, ndofs, fe->getReferenceNodes(), true);
  for(int i = 0; i < ndofs; ++i)
  {
    const double* xi = dofs.getColumn(i);
    fe->evaluateShapeFunctions(xi, wgts);
    double finterp = interp(f, wgts, ndofs);
    EXPECT_NEAR(f[i], finterp, TOL);
  }

  // STEP 3: compute centroid
  double* xc = new double[ndims];
  fe->computePhysicalCoords(fe->getReferenceCenter(), xc);

  // STEP 4: interpolate
  fe->evaluateShapeFunctions(fe->getReferenceCenter(), wgts);
  double finterp = interp(f, wgts, nnodes);

  // STEP 5: check interpolation
  const double fexpected = analytic_function(xc, ndims);
  EXPECT_NEAR(fexpected, finterp, TOL);

  // STEP 6: clean up
  delete fe;
  delete[] xc;
  delete[] f;
  delete[] wgts;
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_single_fe, check_override_max_newton)
{
  constexpr int MAX_NEWTON = 42;  // test value to override max newton

  // STEP 0: construct a mesh with a single element
  mint::FiniteElement* fe = nullptr;
  get_single_fe<MINT_LAGRANGE_BASIS, mint::QUAD>(fe);

  EXPECT_FALSE(MAX_NEWTON == fe->getMaxSolverIterations());

  // STEP 1: override max newton iterations
  fe->setMaxSolverIterations(MAX_NEWTON);
  EXPECT_EQ(MAX_NEWTON, fe->getMaxSolverIterations());

  // clean up
  delete fe;
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, matrix_constructor_deepcopy)
{
  // STEP 0: constants used in the test
  constexpr int NROWS = 2;
  constexpr int NCOLS = 3;
  constexpr int NSIZE = NROWS * NCOLS;

  // STEP 1: setup element coordinates matrix
  numerics::Matrix<double> M(NROWS, NCOLS);
  M.fillColumn(0, 1.0);
  M.fillColumn(1, 2.0);
  M.fillColumn(2, 3.0);

  // STEP 2: construct FE object by making a deep-copy
  mint::FiniteElement fe(M, mint::TRIANGLE);
  EXPECT_FALSE(fe.usesExternalBuffer());

  // STEP 3: ensure FE and Matrix objects are pointing to the same buffer
  double* physical_nodes = fe.getPhysicalNodes();
  EXPECT_FALSE(physical_nodes == M.data());

  // STEP 4: ensure the contents are the same
  const double* expected = M.data();
  for(int i = 0; i < NSIZE; ++i)
  {
    EXPECT_DOUBLE_EQ(expected[i], physical_nodes[i]);
  }

  // STEP 5: change the matrix data
  M.swapColumns(0, 2);
  M.swapColumns(1, 2);

  // STEP 6: ensure the contents *are not* the same
  for(int i = 0; i < NSIZE; ++i)
  {
    EXPECT_FALSE(utilities::isNearlyEqual(M.data()[i], physical_nodes[i]));
  }
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, matrix_constructor_shallowcopy)
{
  // STEP 0: constants used in the test
  constexpr int NROWS = 2;
  constexpr int NCOLS = 3;
  constexpr int NSIZE = NROWS * NCOLS;

  // STEP 1: setup element coordinates matrix
  numerics::Matrix<double> M(2, 3);
  M.fillColumn(0, 1.0);
  M.fillColumn(1, 2.0);
  M.fillColumn(2, 3.0);

  // STEP 2: construct FE object by making a shallow-copy
  mint::FiniteElement* fe = new mint::FiniteElement(M, mint::TRIANGLE, true);
  EXPECT_TRUE(fe->usesExternalBuffer());

  // STEP 3: ensure FE and Matrix objects are pointing to the same buffer
  double* physical_nodes = fe->getPhysicalNodes();
  EXPECT_TRUE(physical_nodes == M.data());

  // STEP 4: ensure the contents are the same
  const double* expected = M.data();
  for(int i = 0; i < NSIZE; ++i)
  {
    EXPECT_DOUBLE_EQ(expected[i], physical_nodes[i]);
  }

  // STEP 5: change the matrix data
  M.swapColumns(0, 2);
  M.swapColumns(1, 2);

  // STEP 6: ensure the contents *still* the same
  for(int i = 0; i < NSIZE; ++i)
  {
    EXPECT_DOUBLE_EQ(expected[i], physical_nodes[i]);
  }

  // STEP 7: delete the FE object, ensure Matrix buffer is not corrupted
  delete fe;
  fe = nullptr;
  EXPECT_FALSE(M.data() == nullptr);
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, check_fe_shape_function)
{
  check_shape<MINT_LAGRANGE_BASIS, mint::QUAD>();
  check_shape<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  check_shape<MINT_LAGRANGE_BASIS, mint::TET>();
  check_shape<MINT_LAGRANGE_BASIS, mint::HEX>();
  check_shape<MINT_LAGRANGE_BASIS, mint::PRISM>();
  check_shape<MINT_LAGRANGE_BASIS, mint::PYRAMID>();

  check_shape<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  check_shape<MINT_LAGRANGE_BASIS, mint::HEX27>();
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, check_fe_jacobian)
{
  check_jacobian<MINT_LAGRANGE_BASIS, mint::QUAD>();
  check_jacobian<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  check_jacobian<MINT_LAGRANGE_BASIS, mint::TET>();
  check_jacobian<MINT_LAGRANGE_BASIS, mint::HEX>();
  check_jacobian<MINT_LAGRANGE_BASIS, mint::PRISM>();
  check_jacobian<MINT_LAGRANGE_BASIS, mint::PYRAMID>();

  check_jacobian<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  check_jacobian<MINT_LAGRANGE_BASIS, mint::HEX27>();
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, check_fe_forward_map)
{
  check_forward_map<MINT_LAGRANGE_BASIS, mint::QUAD>();
  check_forward_map<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  check_forward_map<MINT_LAGRANGE_BASIS, mint::TET>();
  check_forward_map<MINT_LAGRANGE_BASIS, mint::HEX>();
  check_forward_map<MINT_LAGRANGE_BASIS, mint::PRISM>();
  check_forward_map<MINT_LAGRANGE_BASIS, mint::PYRAMID>(1.e-5);

  check_forward_map<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  check_forward_map<MINT_LAGRANGE_BASIS, mint::HEX27>();
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, check_fe_inverse_map)
{
  check_inverse_map<MINT_LAGRANGE_BASIS, mint::QUAD>();
  check_inverse_map<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  check_inverse_map<MINT_LAGRANGE_BASIS, mint::TET>();
  check_inverse_map<MINT_LAGRANGE_BASIS, mint::HEX>();
  check_inverse_map<MINT_LAGRANGE_BASIS, mint::PRISM>();
  check_inverse_map<MINT_LAGRANGE_BASIS, mint::PYRAMID>(1.e-5);

  check_inverse_map<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  check_inverse_map<MINT_LAGRANGE_BASIS, mint::HEX27>();
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, check_fe_point_in_cell)
{
  point_in_cell<MINT_LAGRANGE_BASIS, mint::QUAD>();
  point_in_cell<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  point_in_cell<MINT_LAGRANGE_BASIS, mint::TET>();
  point_in_cell<MINT_LAGRANGE_BASIS, mint::HEX>();
  point_in_cell<MINT_LAGRANGE_BASIS, mint::PRISM>();
  point_in_cell<MINT_LAGRANGE_BASIS, mint::PYRAMID>();

  point_in_cell<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  point_in_cell<MINT_LAGRANGE_BASIS, mint::HEX27>(1.e-7);
}

//------------------------------------------------------------------------------
TEST(mint_fem_single_fe, check_fe_interp)
{
  check_interp<MINT_LAGRANGE_BASIS, mint::QUAD>(1.e-12);
  check_interp<MINT_LAGRANGE_BASIS, mint::TRIANGLE>(1.e-12);
  check_interp<MINT_LAGRANGE_BASIS, mint::TET>(1.e-12);
  check_interp<MINT_LAGRANGE_BASIS, mint::HEX>(1.e-24);
  check_interp<MINT_LAGRANGE_BASIS, mint::PRISM>(1.e-12);
  check_interp<MINT_LAGRANGE_BASIS, mint::PYRAMID>(1.e-12);

  check_interp<MINT_LAGRANGE_BASIS, mint::QUAD9>(1.e-24);
  check_interp<MINT_LAGRANGE_BASIS, mint::HEX27>(1.e-24);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
