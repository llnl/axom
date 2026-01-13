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
#include "HMApplication.hpp"

#include <conduit.hpp>
#include <conduit_relay.hpp>

#include <string>

using RuntimePolicy = axom::runtime_policy::Policy;

namespace detail
{
/*!
 * \brief Turn a field on a fine mesh into a matset on the coarse mesh.
 *
 * \param topoName The name of the topology.
 * \param dims The number of cells in each logical dimension.
 * \param refinement The amount of refinement between coarse/fine meshes.
 * \param n_coarse The coarse mesh (it gets the matset).
 * \param n_field The field used for matset creation.
 * \param nmats The number of materials to make.
 */
void heavily_mixed_matset(const std::string &topoName,
                          int dims[3],
                          int refinement,
                          conduit::Node &n_coarse,
                          const conduit::Node &n_field,
                          int nmats)
{
  const auto fine = n_field.as_float64_accessor();
  int nzones = dims[0] * dims[1] * dims[2];
  int nslots = nzones * nmats;
  std::vector<double> vfs(nslots, 0.);

  // break the data range into nmats parts.
  const double matSize = 1000. / nmats;  //fine.max() / nmats;

  const int rdims[] = {dims[0] * refinement, dims[1] * refinement, dims[2] * refinement};

  for(int k = 0; k < dims[2]; k++)
  {
    for(int j = 0; j < dims[1]; j++)
    {
      for(int i = 0; i < dims[0]; i++)
      {
        const int zoneIndex = k * dims[1] * dims[0] + j * dims[0] + i;
        const int kr = k * refinement;
        const int jr = j * refinement;
        const int ir = i * refinement;
        for(int jj = 0; jj < refinement; jj++)
          for(int ii = 0; ii < refinement; ii++)
          {
            const int fine_index = (kr * rdims[0] * rdims[1]) + ((jr + jj) * rdims[0]) + (ir + ii);
            const int matid =
              axom::utilities::clampVal(static_cast<int>(fine[fine_index] / matSize), 0, nmats - 1);
            const int matslot = zoneIndex * nmats + matid;
            vfs[matslot] += 1. / (refinement * refinement);
          }
      }
    }
  }

  std::vector<int> material_ids;
  std::vector<double> volume_fractions;
  std::vector<int> indices;
  std::vector<int> sizes;
  std::vector<int> offsets;
  for(int k = 0; k < dims[2]; k++)
  {
    for(int j = 0; j < dims[1]; j++)
    {
      for(int i = 0; i < dims[0]; i++)
      {
        const int zoneIndex = k * dims[0] * dims[1] + j * dims[0] + i;

        int size = 0;
        offsets.push_back(indices.size());
        for(int m = 0; m < nmats; m++)
        {
          int matslot = zoneIndex * nmats + m;
          if(vfs[matslot] > 0)
          {
            indices.push_back(material_ids.size());
            material_ids.push_back(m + 1);
            volume_fractions.push_back(vfs[matslot]);
            size++;
          }
        }
        sizes.push_back(size);
      }
    }
  }
  conduit::Node &n_matset = n_coarse["matsets/mat"];
  n_matset["topology"] = topoName;
  conduit::Node &n_material_map = n_matset["material_map"];
  for(int i = 0; i < nmats; i++)
  {
    int matno = i + 1;
    const std::string name = axom::fmt::format("mat{:02d}", matno);
    n_material_map[name] = matno;
  }
  n_matset["material_ids"].set(material_ids);
  n_matset["volume_fractions"].set(volume_fractions);
  n_matset["indices"].set(indices);
  n_matset["sizes"].set(sizes);
  n_matset["offsets"].set(offsets);

  n_coarse["fields/nmats/association"] = "element";
  n_coarse["fields/nmats/topology"] = topoName;
  n_coarse["fields/nmats/values"].set(sizes);
}

template <typename CPUExecSpace>
void heavily_mixed(conduit::Node &n_mesh, int dims[3], int refinement, int nmats)
{
  const int rdims[] = {refinement * dims[0], refinement * dims[1], refinement * dims[2]};

  // Default window
  const conduit::float64 x_min = -0.6;
  const conduit::float64 x_max = 0.6;
  const conduit::float64 y_min = -0.5;
  const conduit::float64 y_max = 0.5;
  const conduit::float64 c_re = -0.5125;
  const conduit::float64 c_im = 0.5213;

  conduit::blueprint::mesh::examples::julia(dims[0], dims[1], x_min, x_max, y_min, y_max, c_re, c_im, n_mesh);
  if(dims[2] > 1)
  {
    // Add another dimension to the coordset.
    const conduit::float64 z_min = 0.;
    const conduit::float64 z_max = x_max - x_min;
    std::vector<conduit::float64> z;
    for(int i = 0; i <= dims[2]; i++)
    {
      const auto t = static_cast<conduit::float64>(i) / dims[2];
      const auto zc = axom::utilities::lerp(z_min, z_max, t);
      z.push_back(zc);
    }
    n_mesh["coordsets/coords/values/z"].set(z);

    // Destination window
    const conduit::float64 s = 0.9;
    const conduit::float64 x1_min = x_min * s;
    const conduit::float64 x1_max = x_max * s;
    const conduit::float64 y1_min = y_min * s;
    const conduit::float64 y1_max = y_max * s;

    conduit::Node n_field;
    n_field.set(conduit::DataType::int32(rdims[0] * rdims[1] * rdims[2]));
    conduit::int32 *destPtr = n_field.as_int32_ptr();
    axom::for_all<CPUExecSpace>(
      rdims[2],
      AXOM_LAMBDA(int k) {
        const auto t = static_cast<conduit::float64>(k) / (dims[2] - 1);
        // Interpolate the window
        const conduit::float64 x0 = axom::utilities::lerp(x_min, x1_min, t);
        const conduit::float64 x1 = axom::utilities::lerp(x_max, x1_max, t);
        const conduit::float64 y0 = axom::utilities::lerp(y_min, y1_min, t);
        const conduit::float64 y1 = axom::utilities::lerp(y_max, y1_max, t);
        conduit::Node n_rmesh;
        conduit::blueprint::mesh::examples::julia(rdims[0], rdims[1], x0, x1, y0, y1, c_re, c_im, n_rmesh);
        const conduit::Node &n_src_field = n_rmesh["fields/iters/values"];
        const conduit::int32 *srcPtr = n_src_field.as_int32_ptr();
        conduit::int32 *currentDestPtr = destPtr + k * rdims[0] * rdims[1];
        axom::copy(currentDestPtr, srcPtr, rdims[0] * rdims[1] * sizeof(conduit::int32));
#ifndef AXOM_DEVICE_CODE
        SLIC_INFO(axom::fmt::format("Made slice {}/{}", k + 1, rdims[2]));
#endif
      });

    // Make a matset based on the higher resolution julia field.
    heavily_mixed_matset("topo", dims, refinement, n_mesh, n_field, nmats);
  }
  else
  {
    // Generate the same julia set at higher resolution to use as materials.
    conduit::Node n_rmesh;
    conduit::blueprint::mesh::examples::julia(rdims[0],
                                              rdims[1],
                                              x_min,
                                              x_max,
                                              y_min,
                                              y_max,
                                              c_re,
                                              c_im,
                                              n_rmesh);

    // Make a matset based on the higher resolution julia field.
    const conduit::Node &n_field = n_rmesh["fields/iters/values"];
    heavily_mixed_matset("topo", dims, refinement, n_mesh, n_field, nmats);
  }
}

}  // end namespace detail

//--------------------------------------------------------------------------------
HMApplication::HMApplication()
  : m_handler(true)
  , m_dims {100, 100, 1}
  , m_numMaterials(40)
  , m_refinement(40)
  , m_numTrials(1)
  , m_writeFiles(true)
  , m_outputFilePath("output")
  , m_method("elvira")
  , m_policy(RuntimePolicy::seq)
  , m_annotationMode("report")
{ }

//--------------------------------------------------------------------------------
int HMApplication::initialize(int argc, char **argv)
{
  axom::CLI::App app;
  app.add_flag("--handler", m_handler)
    ->description("Install a custom error handler that loops forever.")
    ->capture_default_str();

  std::vector<int> dims;
  app.add_option("--dims", dims, "Dimensions in x,y,z (z optional)")
    ->expected(2, 3)
    ->check(axom::CLI::PositiveNumber);

  app.add_option("--method", m_method)
    ->description("The MIR method (or operation) name (equiz, elvira, traversal)");
  app.add_option("--materials", m_numMaterials)
    ->check(axom::CLI::PositiveNumber)
    ->description("The number of materials to create.");
  app.add_option("--refinement", m_refinement)
    ->check(axom::CLI::PositiveNumber)
    ->description("The refinement within a zone when creating materials.");
  app.add_option("--output", m_outputFilePath)
    ->description("The file path for HDF5/YAML output files");
  bool disable_write = !m_writeFiles;
  app.add_flag("--disable-write", disable_write)->description("Disable writing data files");
  app.add_option("--trials", m_numTrials)
    ->check(axom::CLI::PositiveNumber)
    ->description("The number of MIR trials to run on the mesh.");

#if defined(AXOM_USE_CALIPER)
  app.add_option("--caliper", m_annotationMode)
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
  app.add_option("-p, --policy", m_policy, pol_sstr.str())
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

  // Parse command line options.
  int retval = 0;
  try
  {
    app.parse(argc, argv);
    m_writeFiles = !disable_write;
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

  for(size_t i = 0; i < axom::utilities::min(dims.size(), static_cast<size_t>(3)); i++)
  {
    m_dims[i] = dims[i];
  }

  return retval;
}

//--------------------------------------------------------------------------------
int HMApplication::execute()
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);
  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);

  if(m_handler)
  {
    conduit::utils::set_error_handler(conduit_debug_err_handler);
  }
#if defined(AXOM_USE_CALIPER)
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(m_annotationMode);
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
int HMApplication::runMIR()
{
  // Initialize a mesh for testing MIR
  auto timer = axom::utilities::Timer(true);
  conduit::Node mesh;
  {
    AXOM_ANNOTATE_SCOPE("generate");
    int dims[3];
    // NOTE: Use axom::copy to copy m_dims into dims. This way, the Umpire resource
    //       manager gets created now so by the time we have multiple threads using
    //       axom::copy, there is no race condition.
    axom::copy(dims, m_dims, sizeof(int) * 3);

    SLIC_INFO(axom::fmt::format("dims: {},{},{}", dims[0], dims[1], dims[2]));
    SLIC_INFO(axom::fmt::format("refinement: {}", m_refinement));
    SLIC_INFO(axom::fmt::format("numMaterials: {}", m_numMaterials));
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_OPENMP)
    using CPUExecSpace = axom::OMP_EXEC;
#else
    using CPUExecSpace = axom::SEQ_EXEC;
#endif
    detail::heavily_mixed<CPUExecSpace>(mesh, dims, m_refinement, m_numMaterials);
  }
  timer.stop();
  SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Output initial mesh.
  if(m_writeFiles)
  {
    saveMesh(mesh, "heavily_mixed");
  }

  // Begin material interface reconstruction
  timer.reset();
  timer.start();
  conduit::Node options, resultMesh;
  options["matset"] = "mat";
  options["method"] = m_method;  // pass method via options.
  options["trials"] = m_numTrials;

  const int dimension = (m_dims[2] > 1) ? 3 : 2;
  int retval = 0;
  if(m_policy == RuntimePolicy::seq)
  {
    retval = runMIR_seq(dimension, mesh, options, resultMesh);
  }
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_OPENMP)
  else if(m_policy == RuntimePolicy::omp)
  {
    retval = runMIR_omp(dimension, mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_CUDA)
  else if(m_policy == RuntimePolicy::cuda)
  {
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
    retval = runMIR_cuda(dimension, mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_HIP)
  else if(m_policy == RuntimePolicy::hip)
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
  if(m_writeFiles)
  {
    AXOM_ANNOTATE_SCOPE("save_output");
    saveMesh(resultMesh, m_outputFilePath);
  }

  return retval;
}

//--------------------------------------------------------------------------------
void HMApplication::adjustMesh(conduit::Node &) { }

//--------------------------------------------------------------------------------
void HMApplication::saveMesh(const conduit::Node &n_mesh, const std::string &path)
{
#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
  std::string protocol("hdf5");
#else
  std::string protocol("yaml");
#endif
  conduit::relay::io::blueprint::save_mesh(n_mesh, path, protocol);
}

//--------------------------------------------------------------------------------
void HMApplication::conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  SLIC_ERROR(axom::fmt::format("Error from Conduit: s1={}, s2={}, i1={}", s1, s2, i1));
  // This is on purpose.
  while(1);
}
