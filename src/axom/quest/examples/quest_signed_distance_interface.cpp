// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
// _quest_distance_interface_include_start
#include "axom/quest.hpp"
// _quest_distance_interface_include_end
#include "axom/slic.hpp"

#include "axom/CLI11.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#else
using MPI_Comm = int;
#endif

// C/C++ includes
#include <cstring>

/*!
 * \file quest_signed_distance_interfaces.cpp
 *
 * Example that illustrates the use of Quest's C-style Signed Distance 
 * interface from within a C++ application. Supports sequential execution ("seq")
 * as well as accelerated execution via RAJA using openmp on the host ("omp") 
 * and cuda or hip on the device ("gpu")
 */

// namespace aliases
namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace utilities = axom::utilities;

// Predeclare types
struct Arguments;

// Function prototypes
void initialize_logger();
void finalize_logger();
void generate_uniform_box_mesh(mint::UniformMesh*& mesh, Arguments& args);
void run_batched_query(mint::UniformMesh*& mesh);

//------------------------------------------------------------------------------
// GLOBALS
//------------------------------------------------------------------------------
MPI_Comm global_comm;
int mpirank;
int numranks;

const std::map<std::string, quest::SignedDistExec> validExecPolicies {
  {"seq", quest::SignedDistExec::CPU},
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  {"omp", quest::SignedDistExec::OpenMP},
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  {"gpu", quest::SignedDistExec::GPU}
#endif
};

/*!
 * \brief Holds command-line arguments
 */
struct Arguments
{
  std::string fileName;
  int ndims {3};
  std::vector<axom::IndexType> box_dims {32, 32, 32};
  std::vector<double> box_min;
  std::vector<double> box_max;
  bool is_water_tight {true};
  bool dump_vtk {true};
  bool use_shared {false};
  bool use_batched_query {false};
  bool ignore_signs {false};
  quest::SignedDistExec exec_space {quest::SignedDistExec::CPU};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-f,--file", this->fileName, "specifies the input mesh file")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_option("--dimension", this->ndims, "the problem dimension")->capture_default_str();

    app.add_option("--box-dims", box_dims, "the dimensions of the box mesh")->expected(3);

    // If either box-min or box-max are provided, they must both be present
    auto* minbb =
      app.add_option("--box-min", box_min, "the lower corner of the box mesh")->expected(3);
    auto* maxbb =
      app.add_option("--box-max", box_max, "the upper corner of the box mesh")->expected(3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    app.add_flag("!--no-vtk", this->dump_vtk, "disables VTK output")->capture_default_str();

    app
      .add_flag("!--not-watertight", this->is_water_tight, "indicates that input is not water-tight")
      ->capture_default_str();

    app.add_flag("--use-shared", this->use_shared, "stores the surface using MPI-3 shared memory")
      ->capture_default_str();

    app
      .add_flag("--batched",
                this->use_batched_query,
                "uses a single batched query on all points instead of many "
                "individual queries")
      ->capture_default_str();

    app.add_flag("--ignore-signs", this->ignore_signs, "distance query should ignore signs")
      ->capture_default_str();

    std::string pol_info = "Sets execution space of the SignedDistance query.\n";
    pol_info += "Set to \'seq\' to use sequential execution policy.";
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
    pol_info += "\nSet to \'omp\' to use an OpenMP execution policy.";
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
    pol_info += "\nSet to \'gpu\' to use a GPU execution policy.";
#endif
    app.add_option("-e, --exec_space", this->exec_space, pol_info)
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validExecPolicies));

    app.get_formatter()->column_width(40);

    // could throw an exception
    app.parse(argc, argv);

    SLIC_ERROR_IF((this->ndims != 3), "The signed distance is currently only supported in 3-D");

    slic::flushStreams();
  }
};

//------------------------------------------------------------------------------
// PROGRAM MAIN
//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // STEP 0: initialize MPI
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
  global_comm = MPI_COMM_WORLD;

  MPI_Comm_rank(global_comm, &mpirank);
  MPI_Comm_size(global_comm, &numranks);
#else
  mpirank = 0;
  numranks = 1;
#endif

  utilities::Timer timer;

  // STEP 1: initialize the logger
  initialize_logger();

  // STEP 2: parse command line arguments
  Arguments args;
  axom::CLI::App app {"Driver for signed distance query"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(mpirank == 0)
    {
      retval = app.exit(e);
    }
    finalize_logger();

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(retval);
  }
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  if(args.exec_space == quest::SignedDistExec::GPU)
  {
  #if defined(__CUDACC__)
    using GPUExec = axom::CUDA_EXEC<256>;
  #elif defined(__HIPCC__)
    using GPUExec = axom::HIP_EXEC<256>;
  #else
    using GPUExec = axom::SEQ_EXEC;
  #endif
    axom::setDefaultAllocator(axom::execution_space<GPUExec>::allocatorID());
  }
#endif

  // STEP 3: initialize the signed distance interface
  SLIC_INFO("initializing signed distance function...");
  SLIC_INFO("input file: " << args.fileName);
  SLIC_INFO("exec_space=" << (int)args.exec_space);
  slic::flushStreams();

  timer.start();
  quest::signed_distance_use_shared_memory(args.use_shared);
  quest::signed_distance_set_closed_surface(args.is_water_tight);
  quest::signed_distance_set_compute_signs(!args.ignore_signs);
  quest::signed_distance_set_execution_space(args.exec_space);
  // _quest_distance_interface_init_start
  int rc = quest::signed_distance_init(args.fileName, global_comm);
  // _quest_distance_interface_init_end
  timer.stop();

  SLIC_ERROR_IF((rc != 0), "Signed Distance query initialization failed!");
  SLIC_INFO("time to initialize: " << timer.elapsed() << "s");
  slic::flushStreams();

  // STEP 5: Generate computational mesh
  mint::UniformMesh* mesh = nullptr;
  generate_uniform_box_mesh(mesh, args);
  SLIC_ERROR_IF(mesh == nullptr, "box mesh is null!");
  double* phi = mesh->createField<double>("phi", mint::NODE_CENTERED);

  // STEP 6: evaluate the signed distance field on the given mesh
  SLIC_INFO("evaluating signed distance field on specified box mesh...");
  slic::flushStreams();

  if(!args.use_batched_query)
  {
    const axom::IndexType nnodes = mesh->getNumberOfNodes();

    timer.reset();
    timer.start();
    for(axom::IndexType inode = 0; inode < nnodes; ++inode)
    {
      double pt[3];
      mesh->getNode(inode, pt);
      // _quest_distance_interface_test_start
      phi[inode] = quest::signed_distance_evaluate(pt[0], pt[1], pt[2]);
      // _quest_distance_interface_test_end
    }  // END for all nodes
    timer.stop();
    SLIC_INFO("time to evaluate: " << timer.elapsed() << "s");
    slic::flushStreams();
  }
  else
  {
    run_batched_query(mesh);
  }

  // STEP 7: vtk output
  if(args.dump_vtk)
  {
    SLIC_INFO("writing vtk output");
    slic::flushStreams();

    std::ostringstream oss;
    oss << "uniform_mesh_" << mpirank << ".vtk";
    mint::write_vtk(mesh, oss.str());
  }

  // STEP 8: finalize
  delete mesh;
  mesh = nullptr;

  // _quest_distance_interface_finalize_start
  quest::signed_distance_finalize();
  // _quest_distance_interface_finalize_end
  finalize_logger();
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}

//------------------------------------------------------------------------------
//  FUNCTION PROTOTYPE IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void generate_uniform_box_mesh(mint::UniformMesh*& mesh, Arguments& args)
{
  SLIC_ASSERT(mesh == nullptr);
  SLIC_ASSERT(quest::signed_distance_initialized());

  double mesh_box_min[3];
  double mesh_box_max[3];

  double* lo = nullptr;
  double* hi = nullptr;

  if(!args.box_min.empty() && !args.box_max.empty())
  {
    SLIC_INFO("using specified box bounds");
    lo = args.box_min.data();
    hi = args.box_max.data();
  }
  else
  {
    SLIC_INFO("computing mesh bounds...");
    quest::signed_distance_get_mesh_bounds(mesh_box_min, mesh_box_max);
    lo = mesh_box_min;
    hi = mesh_box_max;
  }

  //output some information
  {
    const primal::Point<double, 3> lowerPoint(lo, 3);
    const primal::Point<double, 3> upperPoint(hi, 3);
    auto bbox = primal::BoundingBox<double, 3>(lowerPoint, upperPoint);
    SLIC_INFO("bounding box " << bbox);

    const primal::Point<axom::IndexType, 3> bdims(args.box_dims.data(), 3);
    SLIC_INFO("constructing Uniform Mesh of resolution " << bdims);
  }

  axom::IndexType Ni = static_cast<axom::IndexType>(args.box_dims[0]);
  axom::IndexType Nj = static_cast<axom::IndexType>(args.box_dims[1]);
  axom::IndexType Nk = static_cast<axom::IndexType>(args.box_dims[2]);
  mesh = new mint::UniformMesh(lo, hi, Ni, Nj, Nk);

  SLIC_ASSERT(mesh != nullptr);
  slic::flushStreams();
}

//------------------------------------------------------------------------------
void initialize_logger()
{
  // initialize logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  // setup the logstreams
  std::string fmt = "";
  slic::LogStream* logStream = nullptr;

#ifdef AXOM_USE_MPI
  fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::SynchronizedStream(&std::cout, global_comm, fmt);
#else
  fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif

  // register stream objects with the logger
  slic::addStreamToAllMsgLevels(logStream);
}

//------------------------------------------------------------------------------
void finalize_logger()
{
  slic::flushStreams();
  slic::finalize();
}

//------------------------------------------------------------------------------
void run_batched_query(mint::UniformMesh*& mesh)
{
  double* phi = mesh->getFieldPtr<double>("phi", mint::NODE_CENTERED);

  const axom::IndexType nnodes = mesh->getNumberOfNodes();

  utilities::Timer timer;
  timer.start();

  // Allocate space for the coordinate arrays
  double* x = axom::allocate<double>(nnodes);
  double* y = axom::allocate<double>(nnodes);
  double* z = axom::allocate<double>(nnodes);

  // Determine an appropriate policy
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)
  using ExecPolicy = axom::OMP_EXEC;
#else
  using ExecPolicy = axom::SEQ_EXEC;
#endif

  // Fill the coordinate arrays
  mint::for_all_nodes<ExecPolicy, mint::xargs::xyz>(
    mesh,
    AXOM_LAMBDA(axom::IndexType idx, double xx, double yy, double zz) {
      x[idx] = xx;
      y[idx] = yy;
      z[idx] = zz;
    });

  // Call the vectorized version of the signed distance query
  quest::signed_distance_evaluate(x, y, z, nnodes, phi);

  // Deallocate the coordinate arrays
  axom::deallocate(x);
  axom::deallocate(y);
  axom::deallocate(z);

  timer.stop();
  SLIC_INFO("time to evaluate batched query: " << timer.elapsed() << "s");
  slic::flushStreams();
}
