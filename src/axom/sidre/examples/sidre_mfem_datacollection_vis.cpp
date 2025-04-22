// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file sidre_mfem_datacollection_vis.cpp
 * @brief This example code is a basic demonstration of Sidre's
 * MFEMSidreDataCollection class for visualizing a time-marching
 * simulation.  It is a more thorough version of the snippet
 * provided in Sidre's Sphinx documentation.
 */

#include "axom/config.hpp"

// Datacollection header
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"

// MFEM includes - needed to set up simulation
#include "mfem.hpp"

// Create a simple compiler define for whether we're using MPI
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  #define EXAMPLE_USES_MPI
  #include "mpi.h"
#else
  #undef EXAMPLE_USES_MPI
#endif

int main(int argc, char* argv[])
{
#ifdef EXAMPLE_USES_MPI
  MPI_Init(&argc, &argv);
#else
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);
#endif

  mfem::Mesh* mesh = nullptr;

  // Built a 2D mesh with 100 square elements
  auto serial_mesh = mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::QUADRILATERAL);
  mesh = &serial_mesh;

#ifdef EXAMPLE_USES_MPI
  mfem::ParMesh parallel_mesh(MPI_COMM_WORLD, serial_mesh);
  mesh = &parallel_mesh;
#endif

  // Set up the FiniteElementSpace - needed for the grid functions
  // Initialize with H1 elements of order 1
  mfem::H1_FECollection fec(/*order=*/1, /*dimension=*/2);
  mfem::FiniteElementSpace fes(mesh, &fec);

  // _sidredc_vis_construct_start
  // Initialize the datacollection with the mesh
  // Note: all fields (added with RegisterField) must be on this mesh
  axom::sidre::MFEMSidreDataCollection dc("sidre_mfem_datacoll_vis_ex", mesh);
  // _sidredc_vis_construct_end

// This is where the time-dependent operator would be set up...

// Initialize the solution field
#ifdef EXAMPLE_USES_MPI
  mfem::ParFiniteElementSpace par_fes(fes, parallel_mesh);
  mfem::ParGridFunction soln(&par_fes);
#else
  mfem::GridFunction soln(&fes);
#endif
  soln = 0.0;

  // _sidredc_vis_field_start
  // Note: Any number of fields can be registered
  dc.RegisterField("solution", &soln);
  // _sidredc_vis_field_end

  double dt = 0.05;
  // _sidredc_vis_state_start
  // Save the initial state
  dc.SetCycle(0);      // Iteration counter
  dc.SetTime(0.0);     // Simulation time
  dc.SetTimeStep(dt);  // Time step
  // _sidredc_vis_state_end

#ifdef AXOM_USE_HDF5
  std::string sidre_protocol = "sidre_hdf5";
#else
  std::string sidre_protocol = "sidre_conduit_json";
#endif

  // _sidredc_vis_save_start
  // Filename and protocol, both of which are optional
  dc.Save("sidre_mfem_datacoll_vis_ex", sidre_protocol);
  // _sidredc_vis_save_end

  // Sample time parameters
  int n_iter = 10;

  for(int i = 0; i < n_iter; i++)
  {
    // Calculate the next iteration of the solution field...
    // For simplicity, every element in the field is set to the current time
    soln = dt * i;

    // then save it after updating the time information...
    dc.SetCycle(i);
    dc.SetTime(dt * i);
    dc.Save("sidre_mfem_datacoll_vis_ex", sidre_protocol);
  }

#ifdef EXAMPLE_USES_MPI
  MPI_Finalize();
#endif
}
