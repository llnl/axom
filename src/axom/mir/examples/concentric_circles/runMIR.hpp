// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EXAMPLES_CONCENTRIC_CIRCLES_RUNMIR_HPP
#define AXOM_MIR_EXAMPLES_CONCENTRIC_CIRCLES_RUNMIR_HPP
#include "axom/config.hpp"
#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions

template <typename ExecSpace, typename MatsetView>
void test_matset_traversal(MatsetView matsetView)
{
  AXOM_ANNOTATE_SCOPE("test_matset_traversal");
  double vf1, vf2;
  {
    AXOM_ANNOTATE_SCOPE("zoneMaterials");
    axom::ReduceSum<ExecSpace, double> vfSum(0.);
    axom::for_all<ExecSpace>(
      matsetView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        typename MatsetView::IDList ids;
        typename MatsetView::VFList vfs;
        matsetView.zoneMaterials(zoneIndex, ids, vfs);
        double sum = 0.;
        for(axom::IndexType i = 0; i < vfs.size(); i++)
        {
          sum += vfs[i];
        }
        vfSum += sum;
      });
    vf1 = vfSum.get();
  }
  {
    AXOM_ANNOTATE_SCOPE("iterators");
    axom::ReduceSum<ExecSpace, double> vfSum(0.);
    axom::for_all<ExecSpace>(
      matsetView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto end = matsetView.endZone(zoneIndex);
        double sum = 0.;
        for(auto it = matsetView.beginZone(zoneIndex); it != end; it++)
        {
          sum += it.volume_fraction();
        }
        vfSum += sum;
      });
    vf2 = vfSum.get();
  }

  const double eps = 1.e-10;
  SLIC_INFO(axom::fmt::format("test_matset_traversal: vf1={}, vf2={}, nzones={}, result={}",
                              vf1,
                              vf2,
                              matsetView.numberOfZones(),
                              (axom::utilities::abs(vf1 - vf2) < eps) ? "pass" : "fail"));
}

//--------------------------------------------------------------------------------
template <typename ExecSpace>
int installAllocator([[maybe_unused]] size_t initialPoolSizeBytes)
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
  constexpr int MAXMATERIALS = 20;
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
    // Pass allocator_id since it might be custom for device.
    utils::copy<ExecSpace>(deviceMesh, hostMesh, allocator_id);
  }

  const conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  const conduit::Node &n_topology = deviceMesh["topologies/mesh"];
  const conduit::Node &n_matset = deviceMesh["matsets/mat"];
  conduit::Node deviceResult;
  for(int trial = 0; trial < trials; trial++)
  {
    deviceResult.reset();
    if(method == "equiz")
    {
      // _equiz_mir_start
      using namespace axom::bump::views;
      // Make views (we know beforehand which types to make)
      auto coordsetView = make_explicit_coordset<float, NDIMS>::view(n_coordset);
      using CoordsetView = decltype(coordsetView);

      using ShapeType = typename std::conditional<NDIMS == 3, HexShape<int>, QuadShape<int>>::type;
      auto topologyView = make_unstructured_single_shape_topology<ShapeType>::view(n_topology);
      using TopologyView = decltype(topologyView);

      auto matsetView = make_unibuffer_matset<int, float, MAXMATERIALS>::view(n_matset);
      using MatsetView = decltype(matsetView);

      using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
      MIR m(topologyView, coordsetView, matsetView);
      m.setAllocatorID(allocator_id);
      m.execute(deviceMesh, options, deviceResult);
      // _equiz_mir_end
    }
    else if(method == "elvira")
    {
      // Make views (we know beforehand which types to make)
      using namespace axom::bump::views;
      auto coordsetView = make_explicit_coordset<float, NDIMS>::view(n_coordset);
      using CoordsetView = decltype(coordsetView);

      auto topologyView = make_structured_topology<NDIMS>::view(n_topology);
      using TopologyView = decltype(topologyView);
      using IndexingPolicy = typename TopologyView::IndexingPolicy;

      auto matsetView = make_unibuffer_matset<int, float, MAXMATERIALS>::view(n_matset);
      using MatsetView = decltype(matsetView);

      using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
      MIR m(topologyView, coordsetView, matsetView);
      m.setAllocatorID(allocator_id);
      m.execute(deviceMesh, options, deviceResult);
    }
    else if(method == "traversal")
    {
      using namespace axom::bump::views;
      auto matsetView = make_unibuffer_matset<int, float, MAXMATERIALS>::view(n_matset);
      test_matset_traversal<ExecSpace>(matsetView);
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
