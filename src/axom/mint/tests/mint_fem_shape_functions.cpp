// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/NumericLimits.hpp"
#include "axom/mint/fem/shape_functions/Lagrange.hpp"
#include "axom/mint/fem/shape_functions/ShapeFunction.hpp"
#include "axom/mint/fem/FEBasis.hpp"
#include "axom/mint/fem/FEBasisTypes.hpp"
#include "axom/mint/mesh/CellTypes.hpp"

#include "axom/slic.hpp"

using namespace axom;
using mint::Lagrange;
using mint::ShapeFunction;

//------------------------------------------------------------------------------
//  INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
/*!
 * \brief Tests basic attributes of the shape function
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CELLTYPE the corresponding cell type, e.g., MINT_QUAD
 */
template <int BasisType, mint::CellType CELLTYPE>
void reference_element(double TOL = axom::numeric_limits<double>::epsilon())
{
  using FEMType = typename mint::FEBasis<BasisType, CELLTYPE>;
  using ShapeFunctionType = typename FEMType::ShapeFunctionType;
  ShapeFunctionType sf;

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  const mint::CellType ctype = sf.cellType();
  int ctype_val = mint::cellTypeToInt(ctype);
  EXPECT_TRUE((ctype_val >= 0) && (ctype_val < mint::NUM_CELL_TYPES));

  const int type = sf.type();
  EXPECT_TRUE((type >= 0) && (type < MINT_NUM_BASIS_TYPES));

  const int ndofs = sf.numDofs();
  EXPECT_TRUE(ndofs > 0);
  EXPECT_TRUE(sf.maxNewtonIters() >= 16);

  const int ndims = sf.dimension();
  EXPECT_TRUE(ndims > 0 && ndims < 4);

  const double LO = sf.min() - TOL;
  const double HI = sf.max() + TOL;

  double* xi = new double[ndims];
  sf.center(xi);
  for(int i = 0; i < ndims; ++i)
  {
    EXPECT_TRUE(xi[i] > LO);
    EXPECT_TRUE(xi[i] < HI);
  }
  delete[] xi;

  const int N = ndofs * ndims;
  double* coords = new double[N];
  sf.coords(coords);
  for(int i = 0; i < N; ++i)
  {
    EXPECT_TRUE(coords[i] > LO);
    EXPECT_TRUE(coords[i] < HI);
  }
  delete[] coords;
}

/*!
 * \brief Ensures shape functions satisfy the kronecker delta property.
 */
template <int BasisType, mint::CellType CELLTYPE>
void kronecker_delta()
{
  using FEMType = typename mint::FEBasis<BasisType, CELLTYPE>;
  using ShapeFunctionType = typename FEMType::ShapeFunctionType;
  ShapeFunctionType sf;

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  int ndims = sf.dimension();
  int ndofs = sf.numDofs();

  const int N = ndofs * ndims;
  double* phi = new double[ndofs];
  double* coords = new double[N];
  sf.coords(coords);

  for(int i = 0; i < ndofs; ++i)
  {
    double* nc = &coords[i * ndims];
    sf.evaluate(nc, phi);

    for(int j = 0; j < ndofs; ++j)
    {
      double expected = (i == j) ? 1.0 : 0.0;
      EXPECT_DOUBLE_EQ(expected, phi[j]);
    }

  }  // END for all dofs

  delete[] phi;
  delete[] coords;
}

/*!
 * \brief Ensures shape functions satisfy the partition of unity property.
 */
template <int BasisType, mint::CellType CELLTYPE>
void partition_of_unity()
{
  using FEMType = typename mint::FEBasis<BasisType, CELLTYPE>;
  using ShapeFunctionType = typename FEMType::ShapeFunctionType;
  ShapeFunctionType sf;

  SLIC_INFO("checking " << mint::basis_name[BasisType] << " / " << mint::getCellInfo(CELLTYPE).name);

  int ndims = sf.dimension();
  int ndofs = sf.numDofs();
  const int N = ndofs * ndims;
  double* coords = new double[N];
  double* phi = new double[ndofs];
  double* center = new double[ndims];

  sf.coords(coords);
  sf.center(center);

  for(int i = 0; i < ndofs; ++i)
  {
    double* nc = &coords[i * ndims];
    sf.evaluate(nc, phi);

    double sum = 0.0;
    for(int j = 0; j < ndofs; ++j)
    {
      sum += phi[j];
    }  // END

    // the sum of the weights should be unity
    EXPECT_DOUBLE_EQ(1.0, sum);

  }  // END for all dofs

  // test center
  sf.evaluate(center, phi);
  double sum = 0.0;
  for(int i = 0; i < ndofs; ++i)
  {
    sum += phi[i];
  }
  EXPECT_DOUBLE_EQ(1.0, sum);

  // test other interior nodes
  double* nc = new double[ndims];
  for(int i = 0; i < ndofs; ++i)
  {
    double* xi = &coords[i * ndims];
    for(int j = 0; j < ndims; ++j)
    {
      nc[j] = (xi[j] + center[j]) * 0.5;
    }

    sf.evaluate(nc, phi);

    sum = 0.0;
    for(int j = 0; j < ndofs; ++j)
    {
      sum += phi[j];
    }
    EXPECT_DOUBLE_EQ(1.0, sum);
  }

  delete[] phi;
  delete[] coords;
  delete[] center;
  delete[] nc;
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(mint_fem_shape_functions, check_reference_element)
{
  reference_element<MINT_LAGRANGE_BASIS, mint::QUAD>();
  reference_element<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  reference_element<MINT_LAGRANGE_BASIS, mint::TET>();
  reference_element<MINT_LAGRANGE_BASIS, mint::HEX>();
  reference_element<MINT_LAGRANGE_BASIS, mint::PRISM>();
  reference_element<MINT_LAGRANGE_BASIS, mint::PYRAMID>();

  reference_element<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  reference_element<MINT_LAGRANGE_BASIS, mint::HEX27>();
}

//------------------------------------------------------------------------------
TEST(mint_fem_shape_functions, check_kronecker_delta)
{
  kronecker_delta<MINT_LAGRANGE_BASIS, mint::QUAD>();
  kronecker_delta<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  kronecker_delta<MINT_LAGRANGE_BASIS, mint::TET>();
  kronecker_delta<MINT_LAGRANGE_BASIS, mint::HEX>();
  kronecker_delta<MINT_LAGRANGE_BASIS, mint::PRISM>();
  kronecker_delta<MINT_LAGRANGE_BASIS, mint::PYRAMID>();

  kronecker_delta<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  kronecker_delta<MINT_LAGRANGE_BASIS, mint::HEX27>();
}

//------------------------------------------------------------------------------
TEST(mint_fem_shape_functions, check_partition_of_unity)
{
  partition_of_unity<MINT_LAGRANGE_BASIS, mint::QUAD>();
  partition_of_unity<MINT_LAGRANGE_BASIS, mint::TRIANGLE>();
  partition_of_unity<MINT_LAGRANGE_BASIS, mint::TET>();
  partition_of_unity<MINT_LAGRANGE_BASIS, mint::HEX>();
  partition_of_unity<MINT_LAGRANGE_BASIS, mint::PRISM>();
  partition_of_unity<MINT_LAGRANGE_BASIS, mint::PYRAMID>();

  partition_of_unity<MINT_LAGRANGE_BASIS, mint::QUAD9>();
  partition_of_unity<MINT_LAGRANGE_BASIS, mint::HEX27>();
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
