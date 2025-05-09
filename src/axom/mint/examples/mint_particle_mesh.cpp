// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 *
 * \brief Illustrates how to construct and use a ParticleMesh to perform
 *  operations on a set of particles.
 */

// Axom utilities
#include "axom/core.hpp"
#include "axom/mint.hpp"

// namespace aliases
namespace mint = axom::mint;
namespace utilities = axom::utilities;

//------------------------------------------------------------------------------
int main(int AXOM_UNUSED_PARAM(argc), char** AXOM_UNUSED_PARAM(argv))
{
  using int64 = axom::IndexType;
  const axom::IndexType NUM_PARTICLES = 100;
  const int DIMENSION = 3;

  const double HI = 10.0;
  const double LO = -10.0;
  const double VLO = 0.0;
  const double VHI = 1.0;

  // STEP 0: create the ParticleMesh
  mint::ParticleMesh particles(DIMENSION, NUM_PARTICLES);

  // STEP 1: Add fields to the Particles
  double* vx = particles.createField<double>("vx", mint::NODE_CENTERED);
  double* vy = particles.createField<double>("vy", mint::NODE_CENTERED);
  double* vz = particles.createField<double>("vz", mint::NODE_CENTERED);
  int64* id = particles.createField<int64>("id", mint::NODE_CENTERED);

  // STEP 2: grab handle to the particle position arrays
  double* px = particles.getCoordinateArray(mint::X_COORDINATE);
  double* py = particles.getCoordinateArray(mint::Y_COORDINATE);
  double* pz = particles.getCoordinateArray(mint::Z_COORDINATE);

  // STEP 3: loop over the particle data
  const int64 numParticles = particles.getNumberOfNodes();
  for(int64 i = 0; i < numParticles; ++i)
  {
    px[i] = utilities::random_real(LO, HI);
    py[i] = utilities::random_real(LO, HI);
    pz[i] = utilities::random_real(LO, HI);

    vx[i] = utilities::random_real(VLO, VHI);
    vy[i] = utilities::random_real(VLO, VHI);
    vz[i] = utilities::random_real(VLO, VHI);

    id[i] = i;

  }  // END

  // STEP 4: write the particle mesh in VTK format for visualization
  mint::write_vtk(&particles, "particles.vtk");

  return 0;
}
