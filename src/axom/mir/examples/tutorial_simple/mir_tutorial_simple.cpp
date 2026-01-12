// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"

#include "runMIR.hpp"

#include <conduit.hpp>
#include <conduit_relay.hpp>

#include <string>

// namespace aliases
namespace numerics = axom::numerics;
namespace slam = axom::slam;
namespace mir = axom::mir;
namespace fs = axom::utilities::filesystem;

using RuntimePolicy = axom::runtime_policy::Policy;

// Enable when EquiZ supports iteration.
// #define AXOM_EQUIZ_SUPPORTS_ITERATION

//--------------------------------------------------------------------------------
/// Contain program options.
struct Input
{
  int m_test_case {1};  // valid values 1,2,3,4,5
  bool m_should_iterate {false};
  int m_iter_count {0};
  double m_iter_percent {0.};
  bool m_verbose {false};
  bool m_disable_write {false};
  std::string m_output_dir {};
  RuntimePolicy m_policy {RuntimePolicy::seq};
  std::string m_annotationMode {"report"};
  axom::CLI::App m_app {};

  /// Parse command line.
  void parse(int argc, char **argv)
  {
    m_app.add_option("--test-case", m_test_case)
      ->check(axom::CLI::Range(1, 5))
      ->description("Select the test case.")
      ->capture_default_str();

    m_output_dir = axom::utilities::filesystem::getCWD();
    m_app.add_option("--output-dir", m_output_dir)
      ->check(axom::CLI::ExistingDirectory)
      ->description("The directory for HDF5/YAML output files");

#if defined(AXOM_EQUIZ_SUPPORTS_ITERATION)
    m_app.add_option("--iter-count", m_iter_count)
      ->check(axom::CLI::Range(1, 100))
      ->description("The number of iterations for MIR");

    m_app.add_option("--iter-percent", m_iter_percent)
      ->check(axom::CLI::Bound(1.e-6, 10.))
      ->description("The percent error for iterative MIR");
#endif

    m_app.add_flag("-v,--verbose", m_verbose)->description("Verbose output")->capture_default_str();
    m_app.add_flag("--disable-write", m_disable_write)
      ->description("Disable writing data files")
      ->capture_default_str();

#if defined(AXOM_USE_CALIPER)
    m_app.add_option("--caliper", m_annotationMode)
      ->description("caliper annotation mode. Valid options include 'none' and 'report'. ")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    std::stringstream pol_sstr;
    pol_sstr << "Set MIR runtime policy.";
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
    m_app.add_option("-p, --policy", m_policy, pol_sstr.str())
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    m_app.parse(argc, argv);
  }

  bool shouldIterate() const { return m_should_iterate; }
  int numIterations() const { return m_iter_count; }
  int iterPercentage() const { return m_iter_percent; }
  bool writeFiles() const { return !m_disable_write; }
};

//--------------------------------------------------------------------------------
/// Print a Conduit node.
void printNode(const conduit::Node &n)
{
  conduit::Node options;
  options["num_children_threshold"] = 10000;
  options["num_elements_threshold"] = 10000;
  n.to_summary_string_stream(std::cout, options);
}

//--------------------------------------------------------------------------------
/*!
 * \brief Tutorial main showing how to initialize test cases and perform mir.
 */
int main(int argc, char **argv)
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  // Parse arguments
  Input params;
  try
  {
    params.parse(argc, argv);
  }
  catch(const axom::CLI::ParseError &e)
  {
    return params.m_app.exit(e);
  }

#if defined(AXOM_USE_CALIPER)
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(params.m_annotationMode);
#endif

  // Make the mesh
  conduit::Node mesh;
  axom::bump::data::MeshTester tester;
  auto timer = axom::utilities::Timer(true);
  switch(params.m_test_case)
  {
  case 1:
    tester.initTestCaseOne(mesh);
    break;
  case 2:
    tester.initTestCaseTwo(mesh);
    break;
  case 3:
    tester.initTestCaseThree(mesh);
    break;
  case 4:
    tester.initTestCaseFour(mesh);
    break;
  case 5:
  {
    constexpr int GRIDSIZE = 25;
    constexpr int MAXMATERIALS = 12;
    tester.initTestCaseFive(GRIDSIZE, MAXMATERIALS, mesh);
  }
  break;
  case 6:
  {
    constexpr int GRIDSIZE = 15;
    constexpr int MAXMATERIALS = 3;
    tester.initTestCaseSix(GRIDSIZE, MAXMATERIALS, mesh);
  }
  break;
  }
  timer.stop();
  SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
  std::string protocol("hdf5");
#else
  std::string protocol("yaml");
#endif

  // Save input mesh
  if(params.writeFiles())
  {
    std::string filepath, filename("inputMesh");
    filepath = axom::utilities::filesystem::joinPath(params.m_output_dir, filename);
    conduit::relay::io::blueprint::save_mesh(mesh, filepath, protocol);
  }
  if(params.m_verbose)
  {
    SLIC_INFO("Initial mesh:");
    printNode(mesh);
  }

  // Begin material interface reconstruction
  timer.start();

  // Set up options.
  conduit::Node options;
  options["matset"] = "mat";
#if defined(AXOM_EQUIZ_SUPPORTS_ITERATION)
  // Future options
  options["iterate"] = params.shouldIterate() ? 1 : 0;
  options["iterate_percentage"] = params.iterPercentage();
#endif

  int retval {};

  // Run MIR. Note - the runMIR_xxx functions currently handle just the
  // topology types that are created by MeshTester: unstructured (tet, quad, hex).
  conduit::Node resultMesh;
  if(params.m_policy == RuntimePolicy::seq)
  {
    retval = runMIR_seq(mesh, options, resultMesh);
  }
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_OPENMP)
  else if(params.m_policy == RuntimePolicy::omp)
  {
    retval = runMIR_omp(mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_CUDA)
  else if(params.m_policy == RuntimePolicy::cuda)
  {
    retval = runMIR_cuda(mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_HIP)
  else if(params.m_policy == RuntimePolicy::hip)
  {
    retval = runMIR_hip(mesh, options, resultMesh);
  }
  #endif
#endif
  else
  {
    retval = -1;
    SLIC_ERROR("Unhandled policy.");
  }
  timer.stop();
  SLIC_INFO("Reconstruction time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Save output.
  if(retval == 0 && params.writeFiles())
  {
    std::string filepath, filename("processedMesh");
    filepath = axom::utilities::filesystem::joinPath(params.m_output_dir, filename);
    conduit::relay::io::blueprint::save_mesh(resultMesh, filepath, protocol);
  }

  if(params.m_verbose)
  {
    SLIC_INFO("Final mesh:");
    printNode(resultMesh);
  }

  return retval;
}

//--------------------------------------------------------------------------------
