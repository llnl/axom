// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file sidre_mfem_datacollection_restart.cpp
 * @brief This example code is a basic demonstration of Sidre's
 * MFEMSidreDataCollection class for restarting a simulation.
 * 
 * To reload from a file, the cycle number to reload can be
 * specified as a command-line argument.
 * 
 * For example, run
 * @code{.sh}
 * ./sidre_mfem_datacollection_restart # generates cycles 0-9
 * # then...
 * ./sidre_mfem_datacollection_restart 9 # loads cycle 9 and continues
 * @endcode
 */

#include "axom/config.hpp"

// Datacollection header
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"

#include <memory>  // for unique_ptr

// Create a simple compiler define for whether we're using MPI
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  #define EXAMPLE_USES_MPI
  #include "mpi.h"
#else
  #undef EXAMPLE_USES_MPI
#endif

// MFEM includes - needed to set up simulation
#include "mfem.hpp"
#include "axom/CLI11.hpp"

// Stores the state of the simulation - a mesh, fields, and associated objects
class SimulationState
{
public:
  // Initializes the simulation, using an existing file with specified cycle if
  // cycle_to_load is specified (>= 0)
  SimulationState(axom::sidre::MFEMSidreDataCollection& dc, const int cycle_to_load)
    : m_datacoll(dc)
  {
#ifdef EXAMPLE_USES_MPI
    MPI_Comm_rank(m_datacoll.GetComm(), &m_rank);
#endif
    // Check if this is a restart run
    if(cycle_to_load >= 0)
    {
      // If it is, we can load everything in and "unwrap" to fill in the state
      reloadSim(cycle_to_load);
    }
    // Otherwise it's a nominal run so we have to create everything
    // In a realistic simulation this is where an input file might be used
    else
    {
      setupNewSim();
    }
  }

  ~SimulationState()
  {
    if(m_owns_data)
    {
      delete m_mesh;
      delete m_fecoll;
      delete m_fespace;
      delete m_soln_field;
      delete m_qspace;
      delete m_qfunc;
    }
  }

  // A simulated "step" of the simulation
  void step(double dt)
  {
    // Update simulation state variables
    double t = m_datacoll.GetTime();
    t += dt;
    m_datacoll.SetTime(t);

    const int cycle = m_datacoll.GetCycle();
    m_datacoll.SetCycle(cycle + 1);

    // Calculate the next iteration of the solution field...
    // For simplicity, every element in the field is set to the current time
    *m_soln_field = t;
  }

private:
  // Simulation state setup
  void setupNewSim()
  {
    SLIC_INFO_IF(m_rank == 0, "Starting a new simulation");
    // Everything here is managed by the SimulationState object
    m_owns_data = true;

    // Build a 2D mesh with 100 square elements
    m_mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::QUADRILATERAL));

#ifdef EXAMPLE_USES_MPI
    mfem::Mesh* tmp_mesh = m_mesh;
    m_mesh = new mfem::ParMesh(MPI_COMM_WORLD, *tmp_mesh);
    delete tmp_mesh;
#endif
    // Set up the DataCollection with the newly created mesh
    m_datacoll.SetMesh(m_mesh);

    // Set up the FiniteElementSpace - needed for the grid functions
    // Initialize with H1 elements of order 1
    m_fecoll = new mfem::H1_FECollection(/*order=*/1, m_mesh->Dimension());
#ifdef EXAMPLE_USES_MPI
    auto par_mesh = dynamic_cast<mfem::ParMesh*>(m_mesh);
    m_fespace = new mfem::ParFiniteElementSpace(par_mesh, m_fecoll);
#else
    m_fespace = new mfem::FiniteElementSpace(m_mesh, m_fecoll);
#endif

    // Initialize the solution field

    // Set the data to nullptr so the datacollection will initialize it with
    // its own managed data (needed for a restart)
#ifdef EXAMPLE_USES_MPI
    auto par_fespace = dynamic_cast<mfem::ParFiniteElementSpace*>(m_fespace);
    m_soln_field = new mfem::ParGridFunction(par_fespace, static_cast<double*>(nullptr));
#else
    m_soln_field = new mfem::GridFunction(m_fespace, nullptr);
#endif
    m_datacoll.RegisterField("solution", m_soln_field);

    // Intialize to zero as our "initial conditions"
    *m_soln_field = 0.0;

    // Initialize a quadrature function (per-qpt data)
    m_qspace = new mfem::QuadratureSpace(m_mesh, /*order=*/1);
    // Set the data to nullptr so the datacollection will initialize it with
    // its own managed data (needed for a restart)
    m_qfunc = new mfem::QuadratureFunction(m_qspace);
    m_qfunc->NewDataAndSize(nullptr, m_qfunc->GetSpace()->GetSize());
    m_datacoll.RegisterQField("qpt_data", m_qfunc);
    *m_qfunc = 0.0;

    // Set t = 0 state info
    m_datacoll.SetCycle(0);   // Iteration counter
    m_datacoll.SetTime(0.0);  // Simulation time
  }

  // Sets up the state with non-owning pointers
  void reloadSim(const int cycle_to_load)
  {
    m_datacoll.Load(cycle_to_load);
    SLIC_INFO_IF(m_rank == 0,
                 "Reading in existing data and restarting from iteration "
                   << m_datacoll.GetCycle() << " at time " << m_datacoll.GetTime());
    // The Mesh, GridFunction, etc, objects already exist and can be accessed
    m_mesh = m_datacoll.GetMesh();
    m_soln_field = m_datacoll.GetField("solution");
    m_fespace = m_soln_field->FESpace();
    m_fecoll = m_fespace->FEColl();
    m_qfunc = m_datacoll.GetQField("qpt_data");
    m_qspace = dynamic_cast<mfem::QuadratureSpace*>(m_qfunc->GetSpace());
  }

  // FEM-related objects needed as part of a simulation
  // In a real simulation these would be exposed via accessors
  mfem::Mesh* m_mesh {nullptr};
  const mfem::FiniteElementCollection* m_fecoll {nullptr};
  mfem::FiniteElementSpace* m_fespace {nullptr};
  mfem::GridFunction* m_soln_field {nullptr};
  mfem::QuadratureSpace* m_qspace {nullptr};
  mfem::QuadratureFunction* m_qfunc {nullptr};
  bool m_owns_data {false};

  // A reference to the datacollection so it can be updated with time/cycle
  // information on each time step
  axom::sidre::MFEMSidreDataCollection& m_datacoll;

  // The MPI rank, used to display messages on just one node
  int m_rank {0};
};

int main(int argc, char* argv[])
{
  int num_procs = 1;
#ifdef EXAMPLE_USES_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

  // Initialize the datacollection
  const bool owns_mesh_data = false;
  axom::sidre::MFEMSidreDataCollection dc("sidre_mfem_datacoll_restart_ex", nullptr, owns_mesh_data);
#ifdef EXAMPLE_USES_MPI
  dc.SetComm(MPI_COMM_WORLD);
#endif

  // Command-line argument to load in a specific cycle - optional
  axom::CLI::App app {"Example of Axom's MFEMSidreDataCollection for restarts"};
  int cycle_to_load = -1;

  std::vector<std::string> sidre_protocols = {"sidre_conduit_json",
                                              "sidre_json",
                                              "conduit_bin",
                                              "conduit_json",
                                              "json"};
#ifdef AXOM_USE_HDF5
  std::string protocol = "sidre_hdf5";
  sidre_protocols.push_back("sidre_hdf5");
  sidre_protocols.push_back("conduit_hdf5");
#else
  std::string protocol = "sidre_conduit_json";
#endif

  int num_files = -1;
  app.add_option("--cycle", cycle_to_load)
    ->description("Optional simulation cycle to load")
    ->capture_default_str();
  app.add_option("--protocol", protocol)
    ->description("Optional sidre protocol to use for checkpoints and restarts")
    ->check(axom::CLI::IsMember(sidre_protocols))
    ->capture_default_str();
  app.add_option("--num_files", num_files)
    ->description(
      "Optional flag to set the number of output files for parallel "
      "simulations (default one output file per rank)")
    ->capture_default_str()
    ->check(axom::CLI::Range(1, num_procs).description("Range [1,num_procs]"));
  CLI11_PARSE(app, argc, argv);

#ifdef EXAMPLE_USES_MPI
  if(num_files > 0)
  {
    dc.SetNumFiles(num_files);
  }
#endif

  // Initialize the simulation data structures
  SimulationState sim_state(dc, cycle_to_load);

  // This is where the time-dependent operator would be set up...

  // Save initial state of simulation
  dc.Save("sidre_mfem_datacoll_restart_ex", protocol);

  // Sample time parameters
  const int n_iter = 10;
  const int n_checkpoint = 5;
  const double dt = 0.05;

  for(int i = 0; i < n_iter; i++)
  {
    sim_state.step(dt);
    if(i % n_checkpoint == 0)
    {
      // then save it at each checkpoint
      dc.Save("sidre_mfem_datacoll_restart_ex", protocol);
    }
  }

#ifdef EXAMPLE_USES_MPI
  MPI_Finalize();
#endif
}
