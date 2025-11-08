// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions
#include "runMIR.hpp"
#include "MIRApplication.hpp"

#include <conduit.hpp>
#include <conduit_relay.hpp>

#include <string>

using RuntimePolicy = axom::runtime_policy::Policy;

//--------------------------------------------------------------------------------
MIRApplication::MIRApplication()
  : handler(true)
  , gridSize(5)
  , numCircles(2)
  , dimension(2)
  , numTrials(1)
  , writeFiles(true)
  , outputFilePath("output")
  , method("equiz")
  , policy(RuntimePolicy::seq)
  , annotationMode("report")
{ }

//--------------------------------------------------------------------------------
int MIRApplication::initialize(int argc, char **argv)
{
  axom::CLI::App app;
  app.add_flag("--handler", handler)
    ->description("Install a custom error handler that loops forever.")
    ->capture_default_str();
  app.add_option("--gridsize", gridSize)
    ->check(axom::CLI::PositiveNumber)
    ->description("The number of zones along an axis.");
  app.add_option("--method", method)
    ->description("The MIR method (or operation) name (equiz, elvira, traversal)");
  app.add_option("--numcircles", numCircles)
    ->check(axom::CLI::PositiveNumber)
    ->description("The number of circles to use for material creation.");
  app.add_option("--output", outputFilePath)
    ->description("The file path for HDF5/YAML output files");
  app.add_option("--dimension", dimension, "Specify the dimension (2 or 3)")
    ->check(axom::CLI::Range(2, 3));  // Restrict the value to the range [2, 3]
  bool disable_write = !writeFiles;
  app.add_flag("--disable-write", disable_write)->description("Disable writing data files");
  app.add_option("--trials", numTrials)
    ->check(axom::CLI::PositiveNumber)
    ->description("The number of MIR trials to run on the mesh.");

#if defined(AXOM_USE_CALIPER)
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

  std::stringstream pol_sstr;
  pol_sstr << "Set MIR runtime policy method.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  pol_sstr << "\nSet to 'seq' or 0 to use the RAJA sequential policy.";
  #ifdef AXOM_USE_OPENMP
  pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
  #endif
  #ifdef AXOM_USE_CUDA
  pol_sstr << "\nSet to 'cuda' or 2 to use the RAJA CUDA policy.";
  #endif
  #ifdef AXOM_USE_HIP
  pol_sstr << "\nSet to 'hip' or 3 to use the RAJA HIP policy.";
  #endif
#endif
  app.add_option("-p, --policy", policy, pol_sstr.str())
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

  // Parse command line options.
  int retval = 0;
  try
  {
    app.parse(argc, argv);
    writeFiles = !disable_write;
  }
  catch(axom::CLI::CallForHelp &e)
  {
    std::cout << app.help() << std::endl;
    retval = -1;
  }
  catch(axom::CLI::ParseError &e)
  {
    // Handle other parsing errors
    std::cerr << e.what() << std::endl;
    retval = -2;
  }
  return retval;
}

//--------------------------------------------------------------------------------
int MIRApplication::execute()
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);
  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);

  if(handler)
  {
    conduit::utils::set_error_handler(conduit_debug_err_handler);
  }
#if defined(AXOM_USE_CALIPER)
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(annotationMode);
#endif
  int retval = 0;
  try
  {
    retval = runMIR();
  }
  catch(std::invalid_argument const &e)
  {
    SLIC_WARNING("Bad input. " << e.what());
    retval = -2;
  }
  catch(std::out_of_range const &e)
  {
    SLIC_WARNING("Integer overflow. " << e.what());
    retval = -3;
  }
  return retval;
}

//--------------------------------------------------------------------------------
bool MIRApplication::requiresStructuredMesh(const std::string &method) const
{
  return method == "elvira";
}

//--------------------------------------------------------------------------------
int MIRApplication::runMIR()
{
  // Initialize a mesh for testing MIR
  auto timer = axom::utilities::Timer(true);
  axom::bump::data::MeshTester tester;
  conduit::Node mesh;
  {
    AXOM_ANNOTATE_SCOPE("generate");
    tester.setStructured(requiresStructuredMesh(method));
    if(dimension == 3)
    {
      SLIC_INFO("Generating 3D mesh");
      tester.initTestCaseSix(gridSize, numCircles, mesh);
    }
    else
    {
      SLIC_INFO("Generating 2D mesh");
      tester.initTestCaseFive(gridSize, numCircles, mesh);
    }
    adjustMesh(mesh);
  }
  timer.stop();
  SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Output initial mesh.
  if(writeFiles)
  {
    saveMesh(mesh, "concentric_circles");
  }

  // Begin material interface reconstruction
  timer.reset();
  timer.start();
  conduit::Node options, resultMesh;
  options["matset"] = "mat";
  options["method"] = method;  // pass method via options.
  options["trials"] = numTrials;

  int retval = 0;
  if(policy == RuntimePolicy::seq)
  {
    retval = runMIR_seq(dimension, mesh, options, resultMesh);
  }
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_OPENMP)
  else if(policy == RuntimePolicy::omp)
  {
    retval = runMIR_omp(dimension, mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_CUDA)
  else if(policy == RuntimePolicy::cuda)
  {
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
    retval = runMIR_cuda(dimension, mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_HIP)
  else if(policy == RuntimePolicy::hip)
  {
    retval = runMIR_hip(dimension, mesh, options, resultMesh);
  }
  #endif
#endif
  else
  {
    retval = -1;
    SLIC_ERROR("Unhandled policy.");
  }
  timer.stop();
  SLIC_INFO("Material interface reconstruction time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Output results
  if(writeFiles)
  {
    AXOM_ANNOTATE_SCOPE("save_output");
    saveMesh(resultMesh, outputFilePath);
  }

  return retval;
}

//--------------------------------------------------------------------------------
void MIRApplication::adjustMesh(conduit::Node &) { }

//--------------------------------------------------------------------------------
void MIRApplication::saveMesh(const conduit::Node &n_mesh, const std::string &path)
{
#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
  std::string protocol("hdf5");
#else
  std::string protocol("yaml");
#endif
  conduit::relay::io::blueprint::save_mesh(n_mesh, path, protocol);
}

//--------------------------------------------------------------------------------
void MIRApplication::conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  SLIC_ERROR(axom::fmt::format("Error from Conduit: s1={}, s2={}, i1={}", s1, s2, i1));
  // This is on purpose.
  while(1);
}
