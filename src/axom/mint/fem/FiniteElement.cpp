// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "FiniteElement.hpp"

// Axom includes
#include "axom/core/Macros.hpp"  // for AXOM macros
#include "axom/core/Types.hpp"   // for nullptr

// Axom Utils includes
#include "axom/core/utilities/Utilities.hpp"    // for abs()
#include "axom/core/numerics/linear_solve.hpp"  // for linear_solve()
#include "axom/core/numerics/matvecops.hpp"     // for matrix/vector operators

// Mint includes
#include "axom/mint/mesh/CellTypes.hpp"                // for cell type definitions
#include "axom/mint/fem/shape_functions/Lagrange.hpp"  //Lagrange ShapeFunctions
#include "axom/mint/mesh/Mesh.hpp"                     // for mesh data structure
#include "axom/mint/fem/shape_functions/ShapeFunction.hpp"

// Slic includes
#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace detail
{
//------------------------------------------------------------------------------
bool diverged(const double* xi, int N)
{
  const double DIVERGED = 1.e6;

  bool divergence_detected = false;
  for(int i = 0; !divergence_detected && (i < N); ++i)
  {
    divergence_detected = divergence_detected || (xi[i] > DIVERGED);
  }

  return divergence_detected;
}

}  // namespace detail

//------------------------------------------------------------------------------
// FINITE ELEMENT CLASS IMPLEMENTATION
//------------------------------------------------------------------------------

FiniteElement::FiniteElement(numerics::Matrix<double>& M, CellType cellType, bool useExternal)
  : m_dim(M.getNumRows())
  , m_ctype(cellType)
  , m_shape_func_type(MINT_UNDEFINED_BASIS)
  , m_maxNewtonIterations(-1)
  , m_numnodes(M.getNumColumns())
  , m_jac(nullptr)
  , m_xyz(nullptr)
  , m_phi(nullptr)
  , m_phidot(nullptr)
  , m_usingExternal(useExternal)
  , m_shapeFunction(nullptr)
  , m_shapeFunctionDerivatives(nullptr)
  , m_reference_min(-1)
  , m_reference_max(-1)
  , m_reference_dim(-1)
  , m_numdofs(-1)
  , m_reference_coords(nullptr)
  , m_reference_center(nullptr)
{
  this->setUp();

  if(m_usingExternal)
  {
    // shallow copy the data
    m_xyz = M.data();
  }
  else
  {
    // make a deep copy of the data
    const int N = m_dim * m_numnodes;
    for(int i = 0; i < N; ++i)
    {
      m_xyz[i] = M.data()[i];
    }
  }
}

//------------------------------------------------------------------------------
FiniteElement::~FiniteElement() { this->tearDown(); }

//------------------------------------------------------------------------------
int FiniteElement::computeReferenceCoords(const double* xp, double* xr, double TOL)
{
  SLIC_ASSERT(xp != nullptr);
  SLIC_ASSERT(xr != nullptr);
  SLIC_ASSERT(this->getBasisType() != MINT_UNDEFINED_BASIS);

  if(this->getBasisType() == MINT_UNDEFINED_BASIS)
  {
    SLIC_WARNING("No associated FiniteElement basis!");
    return INVERSE_MAP_FAILED;
  }

  const int MAX_DIM = 3;
  double l1norm = 0.0;
  bool converged = false;

  double x[MAX_DIM];    // residual vector
  double psi[MAX_DIM];  // rhs

  // STEP 0: matrix instance to store the jacobian, J
  numerics::Matrix<double> J(m_dim, m_dim, m_jac, true);
  numerics::Matrix<double> coord_matrix(m_dim, m_numnodes, m_xyz, true);

  // STEP 1: set initial guess for Newton-Raphson at the parametric center
  const int ref_dim = this->getReferenceDimension();
  for(int i = 0; i < m_dim; ++i)
  {
    xr[i] = (i < ref_dim) ? this->getReferenceCenter()[i] : 0.0;
  }

  // STEP 2: start Newton-Raphson iteration
  for(int iter = 0; !converged && (iter < m_maxNewtonIterations); ++iter)
  {
    // evaluate the shape functions and jacobian @ \xi
    this->evaluateShapeFunctions(xr, m_phi);
    this->jacobian(xr, J);

    // compute the right-hand side term, -\psi
    for(int i = 0; i < m_dim; ++i)
    {
      psi[i] = (-1.0) * xp[i];
      for(int j = 0; j < m_numnodes; ++j)
      {
        psi[i] += m_phi[j] * coord_matrix(i, j);
      }

      psi[i] = (-1.0) * psi[i];
    }

    // calculate residual
    int rc = numerics::linear_solve(J, psi, x);

    // compute l1norm and generate improvements
    l1norm = 0.0;
    for(int i = 0; i < m_dim; ++i)
    {
      l1norm += utilities::abs(x[i]);
      xr[i] += x[i];
    }

    if(rc != 0)
    {
      double det = numerics::determinant(J);
      SLIC_WARNING("Newton-Raphson failed, system appears singular!");
      SLIC_INFO("singular system => det(J)=" << det);
      SLIC_INFO("l1norm=" << l1norm);
      SLIC_INFO("Newton iteration=" << iter);

      double xr1 = xr[0];
      double xr2 = xr[1];
      double xr3 = (m_dim == 3) ? xr[2] : 0.0;
      SLIC_INFO("xr=[" << xr1 << " " << xr2 << " " << xr3 << "]");
      break;
    }

    // convergence check
    converged = (l1norm < TOL) ? true : false;

    // divergence check
    if(!converged && detail::diverged(xr, m_dim))
    {
      SLIC_WARNING("Newton-Raphson divergence detected!");
      SLIC_INFO("l1norm=" << l1norm << " iter= " << iter);
      break;
    }

  }  // END newton iterations

  if(!converged)
  {
    SLIC_WARNING("Newton-Raphson did not converge!");
    SLIC_INFO("l1norm=" << l1norm);
    return INVERSE_MAP_FAILED;
  }

  int rc = (this->inReferenceElement(xr, TOL)) ? INSIDE_ELEMENT : OUTSIDE_ELEMENT;

  return rc;
}

//------------------------------------------------------------------------------
void FiniteElement::computePhysicalCoords(const double* xr, double* xp)
{
  SLIC_ASSERT(xr != nullptr);
  SLIC_ASSERT(xp != nullptr);
  SLIC_ASSERT(this->getBasisType() != MINT_UNDEFINED_BASIS);

  if(this->getBasisType() == MINT_UNDEFINED_BASIS)
  {
    SLIC_WARNING("No associated FiniteElement basis!");
    return;
  }

  // STEP 0: evaluate the shape functions, update m_phi
  this->evaluateShapeFunctions(xr, m_phi);

  // STEP 1: get a coordinates matrix instance
  numerics::Matrix<double> coords_matrix(m_dim, m_numnodes, m_xyz, true);

  // STEP 2: compute the global (physical) coordinates
  numerics::matrix_vector_multiply(coords_matrix, m_phi, xp);
}

//------------------------------------------------------------------------------
void FiniteElement::jacobian(const double* lc, numerics::Matrix<double>& J)
{
  SLIC_ASSERT(lc != nullptr);
  SLIC_ASSERT(J.getNumColumns() == this->getReferenceDimension());
  SLIC_ASSERT(J.getNumRows() == this->getPhysicalDimension());
  SLIC_ASSERT(this->getBasisType() != MINT_UNDEFINED_BASIS);

  if(this->getBasisType() == MINT_UNDEFINED_BASIS)
  {
    SLIC_WARNING("No associated FiniteElement basis!");
    return;
  }

  // STEP 0: evaluate the derivatives
  this->evaluateDerivatives(lc, m_phidot);
  numerics::Matrix<double> derivs_matrix(m_numdofs, this->getReferenceDimension(), m_phidot, true);

  // STEP 1: get the coordinates
  numerics::Matrix<double> coords_matrix(m_dim, m_numnodes, m_xyz, true);

  // STEP 2: compute the jacobian
  numerics::matrix_multiply(coords_matrix, derivs_matrix, J);
}

//------------------------------------------------------------------------------
void FiniteElement::evaluateShapeFunctions(const double* xr, double* phi)
{
  SLIC_ASSERT(xr != nullptr);
  SLIC_ASSERT(phi != nullptr);
  SLIC_ASSERT(m_shapeFunction != nullptr);
  SLIC_ASSERT(this->getBasisType() != MINT_UNDEFINED_BASIS);

  if(this->getBasisType() == MINT_UNDEFINED_BASIS)
  {
    SLIC_WARNING("No associated FiniteElement basis!");
    return;
  }

  m_shapeFunction(xr, phi);
}

//------------------------------------------------------------------------------
void FiniteElement::evaluateDerivatives(const double* xr, double* phidot)
{
  SLIC_ASSERT(xr != nullptr);
  SLIC_ASSERT(phidot != nullptr);
  SLIC_ASSERT(m_shapeFunctionDerivatives != nullptr);
  SLIC_ASSERT(this->getBasisType() != MINT_UNDEFINED_BASIS);

  if(this->getBasisType() == MINT_UNDEFINED_BASIS)
  {
    SLIC_WARNING("No associated FiniteElement basis!");
    return;
  }

  m_shapeFunctionDerivatives(xr, phidot);
}

//------------------------------------------------------------------------------
//  PRIVATE METHODS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void FiniteElement::setUp()
{
  m_jac = new double[m_dim * m_dim];

  if(!m_usingExternal)
  {
    m_xyz = new double[m_numnodes * m_dim];
  }

  m_phi = new double[m_numnodes];
  m_phidot = new double[m_numnodes * m_dim];
  m_reference_coords = new double[m_numnodes * m_dim];
  m_reference_center = new double[m_dim];
}

//------------------------------------------------------------------------------
void FiniteElement::tearDown()
{
  delete[] m_jac;

  if(!m_usingExternal)
  {
    delete[] m_xyz;
  }

  delete[] m_phi;
  delete[] m_phidot;
  delete[] m_reference_coords;
  delete[] m_reference_center;
}

//------------------------------------------------------------------------------
bool FiniteElement::inReferenceElement(const double* xi, double TOL)
{
  SLIC_ASSERT(xi != nullptr);

  const double LTOL = m_reference_min - TOL;
  const double HTOL = m_reference_max + TOL;

  bool is_inside = true;

  switch(m_ctype)
  {
  case mint::TRIANGLE:
  case mint::TET:
  case mint::PRISM:
  case mint::PYRAMID:
    this->evaluateShapeFunctions(xi, m_phi);
    for(int i = 0; is_inside && (i < m_numdofs); ++i)
    {
      is_inside = is_inside && (m_phi[i] > LTOL) && (m_phi[i] < HTOL);
    }
    break;
  default:
    for(int i = 0; is_inside && i < m_dim; ++i)
    {
      is_inside = is_inside && (xi[i] > LTOL) && (xi[i] < HTOL);
    }
  }  // END switch

  return is_inside;
}

} /* namespace mint */
} /* namespace axom */
