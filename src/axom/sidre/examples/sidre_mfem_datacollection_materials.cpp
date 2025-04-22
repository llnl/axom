// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file sidre_mfem_datacollection_materials.cpp
 * @brief This example code is a basic demonstration of Sidre's
 * MFEMSidreDataCollection class for material-based fields
 * and species-based fields.
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

/**
 * Utility function to construct a serial or parallel grid function
 * and set all values to the provided \a value
 * The data collection will deallocate these during its destuction
 */
mfem::GridFunction* createAndInitGridFunction(mfem::FiniteElementSpace* fes, double value)
{
#ifdef EXAMPLE_USES_MPI
  auto* pfes = dynamic_cast<mfem::ParFiniteElementSpace*>(fes);
  SLIC_ASSERT(pfes != nullptr);
  mfem::GridFunction* gf = new mfem::ParGridFunction(pfes);
#else
  mfem::GridFunction* gf = new mfem::GridFunction(fes);
#endif

  (*gf) = value;

  return gf;
}

int main(int argc, char* argv[])
{
#ifdef EXAMPLE_USES_MPI
  MPI_Init(&argc, &argv);
#else
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);
#endif

  // Create a mesh; the data collection will delete it
  mfem::Mesh* mesh = nullptr;

  // Built a 2D mesh with 100 square elements
  auto* serial_mesh =
    new mfem::Mesh(mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::QUADRILATERAL));
  mesh = serial_mesh;

#ifdef EXAMPLE_USES_MPI
  auto* parallel_mesh = new mfem::ParMesh(MPI_COMM_WORLD, *serial_mesh);
  mesh = parallel_mesh;
  delete serial_mesh;
#endif

  // Set up the FiniteElementSpace - needed for the grid functions
  // Initialize with H1 elements of order 1
  mfem::H1_FECollection fec(/*order=*/1, /*dimension=*/2);

#ifdef EXAMPLE_USES_MPI
  mfem::FiniteElementSpace serial_fes(mesh, &fec);
  mfem::FiniteElementSpace* fes = new mfem::ParFiniteElementSpace(serial_fes, *parallel_mesh);
#else
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, &fec);
#endif

  // Initialize the datacollection with the mesh
  // Note: all fields (added with RegisterField) must be on this mesh
  // Note: the dc is responsible for deallocating the mesh and fields
  constexpr bool own_mesh_data = true;
  axom::sidre::MFEMSidreDataCollection dc("sidre_mfem_datacoll_materials_ex", mesh, own_mesh_data);

  // _sidredc_material_matset_start
  // Inform the DataCollection that volume fraction information for a material set
  // called "matset" will be in fields called "volume_fraction"
  dc.AssociateMaterialSet("volume_fraction", "matset");
  // _sidredc_material_matset_end

  // _sidredc_material_depfield_start
  // Inform the DataCollection of a material-dependent field "density" associated
  // with a material set called "matset"
  dc.AssociateMaterialDependentField("density", "matset");
  // _sidredc_material_depfield_end

  // _sidredc_material_specset_start
  // Inform the DataCollection that species set data for the species set "specset"
  // associated with material set "matset" will be in fields called "partial_density"
  constexpr bool volume_dependent = false;
  dc.AssociateSpeciesSet("partial_density", "specset", "matset", volume_dependent);
  // _sidredc_material_specset_end

  mfem::GridFunction* vol_frac_1 = createAndInitGridFunction(fes, 0.65);
  mfem::GridFunction* vol_frac_2 = createAndInitGridFunction(fes, 0.35);
  // _sidredc_material_matset_register_start
  dc.RegisterField("volume_fraction_001", vol_frac_1);
  dc.RegisterField("volume_fraction_002", vol_frac_2);
  // _sidredc_material_matset_register_end

  mfem::GridFunction* density_independent = createAndInitGridFunction(fes, 1.0);
  mfem::GridFunction* density_dependent_1 = createAndInitGridFunction(fes, 0.2);
  mfem::GridFunction* density_dependent_2 = createAndInitGridFunction(fes, 0.3);
  // _sidredc_material_depfield_register_start
  dc.RegisterField("density", density_independent);
  dc.RegisterField("density_001", density_dependent_1);
  dc.RegisterField("density_002", density_dependent_2);
  // _sidredc_material_depfield_register_end

  mfem::GridFunction* partial_density_1_1 = createAndInitGridFunction(fes, 0.4);
  mfem::GridFunction* partial_density_1_2 = createAndInitGridFunction(fes, 0.3);
  mfem::GridFunction* partial_density_2_1 = createAndInitGridFunction(fes, 0.2);
  mfem::GridFunction* partial_density_2_2 = createAndInitGridFunction(fes, 0.1);
  // _sidredc_material_specset_register_start
  dc.RegisterField("partial_density_001_001", partial_density_1_1);
  dc.RegisterField("partial_density_001_002", partial_density_1_2);
  dc.RegisterField("partial_density_002_001", partial_density_2_1);
  dc.RegisterField("partial_density_002_002", partial_density_2_2);
  // _sidredc_material_specset_register_end

  // This is where the time-dependent operator would be set up...

  // Save the initial state
  dc.SetCycle(0);   // Iteration counter
  dc.SetTime(0.0);  // Simulation time

  // Filename and protocol, both of which are optional
#ifdef AXOM_USE_HDF5
  std::string protocol = "sidre_hdf5";
#else
  std::string protocol = "sidre_conduit_json";
#endif
  dc.Save("sidre_mfem_datacoll_materials_ex", protocol);

  delete fes;
#ifdef EXAMPLE_USES_MPI
  MPI_Finalize();
#endif
}
