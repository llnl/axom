// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EXAMPLES_HEAVILY_MIXED_RUNMIR_HPP
#define AXOM_MIR_EXAMPLES_HEAVILY_MIXED_RUNMIR_HPP
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"

//--------------------------------------------------------------------------------
template <typename ExecSpace>
int installAllocator(size_t initialPoolSizeBytes)
{
  int allocator_id = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_UMPIRE)
  auto &rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(allocator_id);

  const std::string newName = allocator.getName() + "_POOL";
  SLIC_INFO(
    axom::fmt::format("Creating pool allocator {} with {} bytes.", newName, initialPoolSizeBytes));

  // Create a pool on top of the allocator.
  auto pooled = rm.makeAllocator<umpire::strategy::QuickPool>(
    newName,
    allocator,
    initialPoolSizeBytes,  // first_minimum_pool_allocation_size
    1 << 20,               // next_minimum_pool_allocation_size = 1 MiB chunks
    256                    // alignment
  );

  allocator_id = pooled.getId();
#endif
  return allocator_id;
}

//--------------------------------------------------------------------------------
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

  // See whether we were directed to make a memory pool.
  int allocator_id = axom::execution_space<ExecSpace>::allocatorID();
  if(options.has_path("pool_size"))
  {
    const auto pool_size = options["pool_size"].to_uint64();
    allocator_id = installAllocator<ExecSpace>(pool_size);
    SLIC_INFO(axom::fmt::format("Using custom allocator {}", allocator_id));
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

  if(options.has_path("pool_size"))
  {
#if defined(AXOM_USE_UMPIRE)
    try
    {
      auto &rm = umpire::ResourceManager::getInstance();
      umpire::Allocator allocator = rm.getAllocator(allocator_id);
      SLIC_INFO("Allocator Information:");
      SLIC_INFO(axom::fmt::format("\tname: {}", allocator.getName()));
      SLIC_INFO(axom::fmt::format("\thighwatermark: {}", allocator.getHighWatermark()));
      SLIC_INFO(axom::fmt::format("\tcurrentsize: {}", allocator.getCurrentSize()));
      SLIC_INFO(axom::fmt::format("\tactualsize: {}", allocator.getActualSize()));
      SLIC_INFO(axom::fmt::format("\tallocationcount: {}", allocator.getAllocationCount()));
    }
    catch(...)
    {
      SLIC_ERROR("Allocator information could not be retrieved.");
    }
#endif
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
