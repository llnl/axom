// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EXAMPLES_HEAVILY_MIXED_RUNMIR_HPP
#define AXOM_MIR_EXAMPLES_HEAVILY_MIXED_RUNMIR_HPP
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"

template <typename ExecSpace, int NDIMS>
int runMIR(const conduit::Node &hostMesh, const conduit::Node &options, conduit::Node &hostResult)
{
  AXOM_ANNOTATE_SCOPE("runMIR");

  namespace utils = axom::bump::utilities;

  // Pick the method out of the options.
  std::string method("equiz");
  if(options.has_child("method"))
  {
    method = options["method"].as_string();
  }
  SLIC_INFO(axom::fmt::format("Using policy {} for {} {}D",
                              axom::execution_space<ExecSpace>::name(),
                              method,
                              NDIMS));

  // Get the number of times we want to run MIR.
  int trials = 1;
  if(options.has_child("trials"))
  {
    trials = std::max(1, options["trials"].to_int());
  }

  // Check materials.
  constexpr int MAXMATERIALS = 100;
  auto materialInfo = axom::bump::views::materials(hostMesh["matsets/mat"]);
  if(materialInfo.size() >= MAXMATERIALS)
  {
    SLIC_WARNING(
      axom::fmt::format("To use more than {} materials, recompile with "
                        "larger MAXMATERIALS value.",
                        MAXMATERIALS));
    return -4;
  }

  conduit::Node deviceMesh;
  {
    AXOM_ANNOTATE_SCOPE("host->device");
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
  }

  const conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  const conduit::Node &n_topology = deviceMesh["topologies/topo"];
  const conduit::Node &n_matset = deviceMesh["matsets/mat"];
  conduit::Node deviceResult;
  for(int trial = 0; trial < trials; trial++)
  {
    deviceResult.reset();

    // Make views
    using namespace axom::bump::views;
    auto coordsetView = make_rectilinear_coordset<double, NDIMS>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    auto topologyView = make_rectilinear_topology<NDIMS>::view(n_topology);
    using TopologyView = decltype(topologyView);

    auto matsetView = make_unibuffer_matset<int, double, MAXMATERIALS>::view(n_matset);
    using MatsetView = decltype(matsetView);

    if(method == "equiz")
    {
      using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
      MIR m(topologyView, coordsetView, matsetView);
      m.execute(deviceMesh, options, deviceResult);
    }
    else if(method == "elvira")
    {
      using IndexingPolicy = typename TopologyView::IndexingPolicy;
      using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
      MIR m(topologyView, coordsetView, matsetView);
      m.execute(deviceMesh, options, deviceResult);
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Unsupported MIR method {}", method));
    }
  }

  {
    AXOM_ANNOTATE_SCOPE("device->host");
    utils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);
  }

  return 0;
}

// Prototypes.
int runMIR_seq(int dimension,
               const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result);
int runMIR_omp(int dimension,
               const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result);
int runMIR_cuda(int dimension,
                const conduit::Node &mesh,
                const conduit::Node &options,
                conduit::Node &result);
int runMIR_hip(int dimension,
               const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result);

#endif
