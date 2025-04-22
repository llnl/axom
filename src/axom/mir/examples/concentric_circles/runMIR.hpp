// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EXAMPLES_CONCENTRIC_CIRCLES_RUNMIR_HPP
#define AXOM_MIR_EXAMPLES_CONCENTRIC_CIRCLES_RUNMIR_HPP
#include "axom/config.hpp"
#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions

template <typename ExecSpace>
int runMIR(const conduit::Node &hostMesh, const conduit::Node &options, conduit::Node &hostResult)
{
  AXOM_ANNOTATE_BEGIN("runMIR");

  namespace bputils = axom::mir::utilities::blueprint;
  using namespace axom::mir::views;
  SLIC_INFO(axom::fmt::format("Using policy {}", axom::execution_space<ExecSpace>::name()));

  // Check materials.
  constexpr int MAXMATERIALS = 20;
  auto materialInfo = materials(hostMesh["matsets/mat"]);
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
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);
  }

  const conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  const conduit::Node &n_topology = deviceMesh["topologies/mesh"];
  const conduit::Node &n_matset = deviceMesh["matsets/mat"];
  conduit::Node deviceResult;

  // Pick the method out of the options.
  std::string method("equiz");
  if(options.has_child("method"))
  {
    method = options["method"].as_string();
  }

  if(method == "equiz")
  {
    // _equiz_mir_start
    // Make views (we know beforehand which types to make)
    auto coordsetView = make_explicit_coordset<float, 2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    using ShapeType = QuadShape<int>;
    auto topologyView = make_unstructured_single_shape<ShapeType>::view(n_topology);
    using TopologyView = decltype(topologyView);

    auto matsetView = make_unibuffer_matset<int, float, MAXMATERIALS>::view(n_matset);
    using MatsetView = decltype(matsetView);

    using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    m.execute(deviceMesh, options, deviceResult);
    // _equiz_mir_end
  }
  else if(method == "elvira")
  {
    // Make views (we know beforehand which types to make)
    auto coordsetView = make_explicit_coordset<float, 2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    auto topologyView = make_structured<2>::view(n_topology);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    auto matsetView = make_unibuffer_matset<int, float, MAXMATERIALS>::view(n_matset);
    using MatsetView = decltype(matsetView);

    using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    m.execute(deviceMesh, options, deviceResult);
  }
  else
  {
    SLIC_ERROR(axom::fmt::format("Unsupported MIR method {}", method));
  }

  {
    AXOM_ANNOTATE_SCOPE("device->host");
    bputils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);
  }

  return 0;
}

// Prototypes.
int runMIR_seq(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_omp(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_cuda(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_hip(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);

#endif
