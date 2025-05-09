// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/ClipFieldFilter.hpp"

// clang-format off
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  using seq_exec = axom::SEQ_EXEC;

  #if defined(AXOM_USE_OPENMP)
    using omp_exec = axom::OMP_EXEC;
  #endif

  #if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
  #endif

  #if defined(AXOM_USE_HIP)
    constexpr int HIP_BLOCK_SIZE = 64;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
  #endif
#endif
// clang-format on

namespace axom
{
namespace mir
{
namespace clipping
{
void ClipFieldFilter::execute(const conduit::Node &n_input,
                              const conduit::Node &n_options,
                              conduit::Node &n_output)
{
  ClipOptions opts(n_options);
  const std::string clipFieldName = opts.clipField();

  const conduit::Node &n_fields = n_input.fetch_existing("fields");
  const conduit::Node &n_clipField = n_fields.fetch_existing(clipFieldName);
  const std::string &topoName = n_clipField["topology"].as_string();
  const conduit::Node &n_topo = n_input.fetch_existing("topologies/" + topoName);
  const std::string &coordsetName = n_topo["coordset"].as_string();
  const conduit::Node &n_coordset = n_input.fetch_existing("coordsets/" + coordsetName);

  execute(n_topo,
          n_coordset,
          n_fields,
          n_options,
          n_output["topologies/" + opts.topologyName(topoName)],
          n_output["coordsets/" + opts.coordsetName(coordsetName)],
          n_output["fields"]);
}

void ClipFieldFilter::execute(const conduit::Node &n_topo,
                              const conduit::Node &n_coordset,
                              const conduit::Node &n_fields,
                              const conduit::Node &n_options,
                              conduit::Node &n_newTopo,
                              conduit::Node &n_newCoordset,
                              conduit::Node &n_newFields)
{
  // Instantiate the algorithm for the right device and invoke it.
  if(m_runtime == axom::runtime_policy::Policy::seq)
  {
    ClipFieldFilterDevice<seq_exec> clipper;
    clipper.execute(n_topo, n_coordset, n_fields, n_options, n_newTopo, n_newCoordset, n_newFields);
  }
#if 0
  #if defined(AXOM_USE_OPENMP)
  else if(m_runtime == axom::runtime_policy::Policy::omp)
  {
    ClipFieldFilterDevice<omp_exec> clipper;
    clipper.execute(n_topo, n_coordset, n_fields, n_options, n_newTopo, n_newCoordset, n_newFields);
  }
  #endif
  #if defined(AXOM_USE_CUDA)
  else if(m_runtime == axom::runtime_policy::Policy::cuda)
  {
    ClipFieldFilterDevice<cuda_exec> clipper;
    clipper.execute(n_topo, n_coordset, n_fields, n_options, n_newTopo, n_newCoordset, n_newFields);
  }
  #endif
  #if defined(AXOM_USE_HIP)
  else if(m_runtime == axom::runtime_policy::Policy::hip)
  {
    ClipFieldFilterDevice<hip_exec> clipper;
    clipper.execute(n_topo, n_coordset, n_fields, n_options, n_newTopo, n_newCoordset, n_newFields);
  }
  #endif
#endif
}

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom
