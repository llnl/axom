// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core.hpp"
#include "axom/mint.hpp"

// C/C++ includes
#include <cmath>

// aliases
namespace mint = axom::mint;
namespace utilities = axom::utilities;
namespace xargs = mint::xargs;
using IndexType = axom::IndexType;

constexpr double R = 2.5;
constexpr double M = (2 * M_PI) / 50.0;
/*!
 * \file
 *
 * \brief Illustrates how to construct and use a mint::CurvilinearMesh object.
 */

int main(int AXOM_UNUSED_PARAM(argc), char** AXOM_UNUSED_PARAM(argv))
{
  constexpr int N = 100;
  constexpr double h = 0.5;

  // STEP 0: construct a N x N curvilinear mesh
  mint::CurvilinearMesh mesh(N, N);
  double* x = mesh.getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh.getCoordinateArray(mint::Y_COORDINATE);

  // STEP 1: add fiducial fields
  double* dx = mesh.createField<double>("dx", mint::NODE_CENTERED);
  double* dy = mesh.createField<double>("dy", mint::NODE_CENTERED);

  // STEP 2: fill in the coordinates
  mint::for_all_nodes<axom::SEQ_EXEC, xargs::ij>(
    &mesh,
    AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j) {
      const double xx = h * i;
      const double yy = h * j;
      const double alpha = yy + R;
      const double beta = xx * M;

      x[nodeIdx] = alpha * cos(beta);
      y[nodeIdx] = alpha * sin(beta);

      dx[nodeIdx] = x[nodeIdx] + y[nodeIdx];
      dy[nodeIdx] = x[nodeIdx] - y[nodeIdx];
    });

  // STEP 3: dump file
  mint::write_vtk(&mesh, "curvilinear_mesh.vtk");

  return 0;
}
