// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_NODE_TO_ZONE_RELATION_BUILDER_HPP_
#define AXOM_BUMP_NODE_TO_ZONE_RELATION_BUILDER_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/dispatch_unstructured_topology.hpp"
#include "axom/bump/MakeUnstructured.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

#include <map>
#include <vector>

namespace axom
{
namespace bump
{
namespace details
{
/*!
 * \brief Build the node to zone relation.
 */
template <typename ExecSpace, typename ViewType>
struct BuildRelation
{
  /*!
   * \brief Given views that contain the nodes and zones, sort the zones using the
   *        node numbers to produce a list of zones for each node and an offsets array
   *        that points to the start of each list of zones.
   * 
   * \param[in]    nodesView   A view that contains the set of all of the nodes in the topology (the connectivity)
   * \param[inout] zonesView   A view (same size as \a nodesView) that contains the zone number of each node.
   * \param[out]   sizesView A view that we fill with sizes.
   * \param[out]   offsetsView A view that we fill with offsets so offsetsView[i] points to the start of the i'th list in \a zonesView.
   *
   * \note axom::sort_pairs can be slow if there are a lot of nodes (depends on ExecSpace too).
   */
  static void execute(ViewType nodesView, ViewType zonesView, ViewType sizesView, ViewType offsetsView)
  {
    AXOM_ANNOTATE_SCOPE("FillZonesAndOffsets");
    SLIC_ASSERT(nodesView.size() == zonesView.size());

    using value_type = typename ViewType::value_type;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Make a copy of the nodes that we'll use as keys.
    const auto n = nodesView.size();
    axom::Array<value_type> keys(n, n, allocatorID);
    auto keysView = keys.view();
    axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i) { keysView[i] = nodesView[i]; });

    // Sort the keys, zones in place. This sorts the zonesView which we want for output.
    axom::sort_pairs<ExecSpace>(keysView, zonesView);

    // Make a mask array for where differences occur.
    axom::Array<axom::IndexType> mask(n, n, allocatorID);
    auto maskView = mask.view();
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) {
        maskView[i] = (i >= 1) ? ((keysView[i] != keysView[i - 1]) ? 1 : 0) : 1;
      });

    // Do a scan on the mask array to build an offset array.
    axom::Array<axom::IndexType> dest_offsets(n, n, allocatorID);
    auto dest_offsetsView = dest_offsets.view();
    axom::exclusive_scan<ExecSpace>(maskView, dest_offsetsView);

    // Build the offsets to each node's zone ids.
    axom::for_all<ExecSpace>(
      offsetsView.size(),
      AXOM_LAMBDA(axom::IndexType i) { offsetsView[i] = 0; });
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) {
        if(maskView[i])
        {
          offsetsView[dest_offsetsView[i]] = i;
        }
      });

    // Compute sizes from offsets.
    const value_type totalSize = nodesView.size();
    axom::for_all<ExecSpace>(
      offsetsView.size(),
      AXOM_LAMBDA(axom::IndexType i) {
        sizesView[i] = (i < offsetsView.size() - 1) ? (offsetsView[i + 1] - offsetsView[i])
                                                    : (totalSize - offsetsView[i]);
      });
  }
};

/// Partial specialization for axom::SEQ_EXEC.
template <typename ViewType>
struct BuildRelation<axom::SEQ_EXEC, ViewType>
{
  static void execute(ViewType nodesView, ViewType zonesView, ViewType sizesView, ViewType offsetsView)
  {
    AXOM_ANNOTATE_SCOPE("FillZonesAndOffsets");
    SLIC_ASSERT(nodesView.size() == zonesView.size());
    using value_type = typename ViewType::value_type;
    using ExecSpace = axom::SEQ_EXEC;
    const int allocatorID = execution_space<ExecSpace>::allocatorID();

    // Count how many times a node is used.
    axom::for_all<ExecSpace>(
      sizesView.size(),
      AXOM_LAMBDA(axom::IndexType index) { sizesView[index] = 0; });
    axom::for_all<ExecSpace>(
      nodesView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        // Works because ExecSpace=SEQ_EXEC.
        sizesView[nodesView[index]]++;
      });
    // Make offsets
    axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);

    axom::for_all<ExecSpace>(
      sizesView.size(),
      AXOM_LAMBDA(axom::IndexType index) { sizesView[index] = 0; });

    // Make a copy of zonesView so we can reorganize zonesView.
    axom::Array<value_type> zcopy(zonesView.size(), zonesView.size(), allocatorID);
    axom::copy(zcopy.data(), zonesView.data(), zonesView.size() * sizeof(value_type));
    auto zcopyView = zcopy.view();

    // Fill in zonesView, sizesView with each node's zones.
    axom::for_all<ExecSpace>(
      nodesView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const auto ni = nodesView[index];
        const auto destOffset = offsetsView[ni] + sizesView[ni];
        zonesView[destOffset] = zcopyView[index];
        // Works because ExecSpace=SEQ_EXEC.
        sizesView[ni]++;
      });
  }
};

}  // end namespace details

/*!
 * \brief Build an o2m relation that lets us look up the zones for a node.
 *
 * \note The zone list for each point is not sorted.
 */
template <typename ExecSpace>
class NodeToZoneRelationBuilder
{
public:
  /*!
   * \brief Build a node to zone relation and store the resulting O2M relation in the \a relation conduit node.
   *
   * \param topo The topology for which we're building the O2M relation.
   * \param coordset The topology's coordset.
   * \param[out] The node that will contain the O2M relation.
   */
  void execute(const conduit::Node &topo, const conduit::Node &coordset, conduit::Node &relation)
  {
    namespace utils = axom::bump::utilities;
    const std::string type = topo.fetch_existing("type").as_string();

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const int conduitAllocatorID = c2a.getConduitAllocatorID();

    conduit::Node &n_zones = relation["zones"];
    conduit::Node &n_sizes = relation["sizes"];
    conduit::Node &n_offsets = relation["offsets"];
    n_zones.set_allocator(conduitAllocatorID);
    n_sizes.set_allocator(conduitAllocatorID);
    n_offsets.set_allocator(conduitAllocatorID);

    if(type == "unstructured")
    {
      conduit::blueprint::mesh::utils::ShapeType shape(topo);
      const conduit::Node &n_connectivity = topo["elements/connectivity"];
      const std::string shapeType = topo["elements/shape"].as_string();
      const auto intTypeId = n_connectivity.dtype().id();
      const auto connSize = n_connectivity.dtype().number_of_elements();

      // Use the coordset to get the number of nodes. Conduit should be able to do this using only metadata.
      const auto nnodes = conduit::blueprint::mesh::utils::coordset::length(coordset);

      if(shape.is_polyhedral())
      {
        views::dispatch_unstructured_polyhedral_topology(
          topo,
          [&](auto AXOM_UNUSED_PARAM(shape), auto topoView) {
            handlePolyhedralView(topoView, n_zones, n_sizes, n_offsets, nnodes, intTypeId);
          });
      }
      else if(shape.is_polygonal() || shapeType == "mixed")
      {
        const conduit::Node &n_topo_sizes = topo["elements/sizes"];
        const conduit::Node &n_topo_offsets = topo["elements/offsets"];

        const auto nzones = n_topo_sizes.dtype().number_of_elements();

        // Allocate Conduit arrays on the device in a data type that matches the connectivity.
        n_zones.set(conduit::DataType(intTypeId, connSize));
        n_sizes.set(conduit::DataType(intTypeId, nnodes));
        n_offsets.set(conduit::DataType(intTypeId, nnodes));

        // Make zones for each node
        views::IndexNode_to_ArrayView_same(
          n_zones,
          n_topo_sizes,
          n_topo_offsets,
          [&](auto zonesView, auto sizesView, auto offsetsView) {
            fillZonesMixed(nzones, zonesView, sizesView, offsetsView);
          });

        views::IndexNode_to_ArrayView_same(
          n_connectivity,
          n_zones,
          n_sizes,
          n_offsets,
          [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView) {
            // Make the relation.
            using ViewType = decltype(connectivityView);
            details::BuildRelation<ExecSpace, ViewType>::execute(connectivityView,
                                                                 zonesView,
                                                                 sizesView,
                                                                 offsetsView);
          });
      }
      else
      {
        // Shapes are all the same size.
        const auto nodesPerShape = shape.indices;

        // Allocate Conduit arrays on the device in a data type that matches the connectivity.
        n_zones.set(conduit::DataType(intTypeId, connSize));
        n_sizes.set(conduit::DataType(intTypeId, nnodes));
        n_offsets.set(conduit::DataType(intTypeId, nnodes));

        views::IndexNode_to_ArrayView_same(
          n_connectivity,
          n_zones,
          n_sizes,
          n_offsets,
          [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView) {
            // Make zones for each node
            fillZones(zonesView, connSize, nodesPerShape);

            // Make the relation.
            using ViewType = decltype(connectivityView);
            details::BuildRelation<ExecSpace, ViewType>::execute(connectivityView,
                                                                 zonesView,
                                                                 sizesView,
                                                                 offsetsView);
          });
      }
    }
    else
    {
      // These are all structured topos of some sort. Make an unstructured representation and recurse.

      conduit::Node mesh;
      MakeUnstructured<ExecSpace>::execute(topo, coordset, "newtopo", mesh);

      // Recurse using the unstructured mesh.
      execute(mesh.fetch_existing("topologies/newtopo"), coordset, relation);
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Handle a polyhedral view.
   *
   * \param topoView A polyhedral topology view.
   * \param[out] n_zones The new zones node for the relation.
   * \param[out] n_sizes The new sizes node for the relation.
   * \param[out] n_offsets The new offsets node for the relation.
   * \param nnodes The number of nodes in the mesh's coordset.
   * \param intTypeId The dtype id for the connectivity.
   * \param connSize The length of the connectivity.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename PHView>
  void handlePolyhedralView(PHView topoView,
                            conduit::Node &n_zones,
                            conduit::Node &n_sizes,
                            conduit::Node &n_offsets,
                            axom::IndexType nnodes,
                            int intTypeId) const
  {
    namespace utils = axom::bump::utilities;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const int conduitAllocatorID = c2a.getConduitAllocatorID();
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    const auto nzones = topoView.numberOfZones();
    axom::Array<axom::IndexType> sizes(nzones, nzones, allocatorID);
    auto sizes_view = sizes.view();

    // Run through the topology once to do a count of each zone's unique node ids.
    axom::ReduceSum<ExecSpace, axom::IndexType> count(0);
    const PHView deviceTopologyView(topoView);
    axom::for_all<ExecSpace>(
      topoView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = deviceTopologyView.zone(zoneIndex);
        const auto uniqueIds = zone.getUniqueIds();
        sizes_view[zoneIndex] = uniqueIds.size();
        count += uniqueIds.size();
      });
    const auto connSize = count.get();

    // Do a scan on the size array to build an offset array.
    axom::Array<axom::IndexType> offsets(nzones, nzones, allocatorID);
    auto offsets_view = offsets.view();
    axom::exclusive_scan<ExecSpace>(sizes_view, offsets_view);
    sizes.clear();

    // Allocate Conduit arrays on the device in a data type that matches the connectivity.
    conduit::Node n_conn;
    n_conn.set_allocator(conduitAllocatorID);
    n_conn.set(conduit::DataType(intTypeId, connSize));

    n_zones.set(conduit::DataType(intTypeId, connSize));
    n_sizes.set(conduit::DataType(intTypeId, nnodes));
    n_offsets.set(conduit::DataType(intTypeId, nnodes));

    views::IndexNode_to_ArrayView_same(
      n_conn,
      n_zones,
      n_sizes,
      n_offsets,
      [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView) {
        fillZonesPH(topoView, connectivityView, zonesView, offsets_view);

        // Make the relation.
        using ViewType = decltype(connectivityView);
        details::BuildRelation<ExecSpace, ViewType>::execute(connectivityView,
                                                             zonesView,
                                                             sizesView,
                                                             offsetsView);
      });
  }

  /*!
   * \brief Fill in the zone numbers for each mixed-sized zone.
   *
   * \param topoView The topology view for the PH mesh.
   * \param connectivityView The view that contains the connectivity.
   * \param zonesView The view that will contain the zone ids.
   * \param offsetsView The view that contains the offsets.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename TopologyView, typename IntegerView, typename OffsetsView>
  void fillZonesPH(const TopologyView &topoView,
                   IntegerView connectivityView,
                   IntegerView zonesView,
                   OffsetsView offsets_view) const
  {
    // Run through the data one more time to build the nodes and zones arrays.
    const TopologyView deviceTopologyView(topoView);
    axom::for_all<ExecSpace>(
      topoView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = deviceTopologyView.zone(zoneIndex);
        const auto uniqueIds = zone.getUniqueIds();
        auto destIdx = offsets_view[zoneIndex];
        for(axom::IndexType i = 0; i < uniqueIds.size(); i++, destIdx++)
        {
          connectivityView[destIdx] = uniqueIds[i];
          zonesView[destIdx] = zoneIndex;
        }
      });
  }

  /*!
   * \brief Fill in the zone numbers for each mixed-sized zone.
   *
   * \param nzones The number of zones.
   * \param zonesView The view that will contain the zone ids.
   * \param sizesView The view that contains the sizes.
   * \param offsetsView The view that contains the offsets.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename IntegerView>
  void fillZonesMixed(axom::IndexType nzones,
                      IntegerView zonesView,
                      IntegerView sizesView,
                      IntegerView offsetsView) const
  {
    using DataType = typename decltype(zonesView)::value_type;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        for(DataType i = 0; i < sizesView[zoneIndex]; i++)
          zonesView[offsetsView[zoneIndex] + i] = zoneIndex;
      });
  }

  /*!
   * \brief Fill in the zone numbers for each node in the connectivity.
   *
   * \param zonesView The view that will contain the zone ids.
   * \param connSize The length of the connectivity.
   * \param nodesPerShape The number of nodes per shape.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename IntegerView>
  void fillZones(IntegerView zonesView, axom::IndexType connSize, axom::IndexType nodesPerShape) const
  {
    axom::for_all<ExecSpace>(
      connSize,
      AXOM_LAMBDA(axom::IndexType index) { zonesView[index] = index / nodesPerShape; });
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
