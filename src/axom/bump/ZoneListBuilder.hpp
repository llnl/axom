// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_ZONELIST_BUILDER_HPP
#define AXOM_BUMP_ZONELIST_BUILDER_HPP

#include "axom/core.hpp"

#include <conduit.hpp>

namespace axom
{
namespace bump
{

/*!
 * \brief This struct builds lists of clean and mixed zones using the input topology and matset views.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The topology view type on which the algorithm will run.
 * \tparam MatsetView The matset view type on which the algorithm will run.
 */
template <typename ExecSpace, typename TopologyView, typename MatsetView>
class ZoneListBuilder
{
  using MaskType = typename axom::bump::utilities::mask_traits<ExecSpace, int>::type;

public:
  using SelectedZonesView = axom::ArrayView<axom::IndexType>;
  using ZoneType = typename TopologyView::ShapeType;

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for creating the zone lists.
   * \param matsetView The matset view to use for creating the zone lists.
   */
  ZoneListBuilder(const TopologyView &topoView, const MatsetView &matsetView)
    : m_topologyView(topoView)
    , m_matsetView(matsetView)
  { }

  /*!
   * \brief Build the list of clean and mixed zones using the number of materials
   *        per zone, maxed to the nodes.
   *
   * \param nnodes The number of nodes in the topology's coordset.
   * \param[out] cleanIndices An array that will contain the list of clean material zone ids.
   * \param[out] mixedIndices An array that will contain the list of mixed material zone ids.
   *
   * \note The clean/mixed index arrays are not strictly what could be determined by the matset alone.
   *       We figure out which nodes touch multiple materials. Then we iterate over the zones and
   *       those that touch only nodes that have 1 material are marked clean, otherwise they are
   *       considered mixed as we might have to split those zones.
   */
  void execute(axom::IndexType nnodes,
               axom::Array<axom::IndexType> &cleanIndices,
               axom::Array<axom::IndexType> &mixedIndices) const
  {
    AXOM_ANNOTATE_SCOPE("ZoneListBuilder.1");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_BEGIN("nMatsPerNode");
    axom::Array<int> nMatsPerNode(axom::ArrayOptions::Uninitialized(), nnodes, nnodes, allocatorID);
    auto nMatsPerNodeView = nMatsPerNode.view();
    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType nodeIndex) { nMatsPerNodeView[nodeIndex] = 1; });

    // Determine max number of materials a node might touch.
    MatsetView deviceMatsetView(m_matsetView);
    const TopologyView deviceTopologyView(m_topologyView);
    axom::for_all<ExecSpace>(
      m_topologyView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const int nmats = deviceMatsetView.numberOfMaterials(zoneIndex);
        if(nmats > 1)
        {
          const auto zone = deviceTopologyView.zone(zoneIndex);
          const auto nnodesThisZone = zone.numberOfNodes();
          int *nodeData = nMatsPerNodeView.data();
          for(axom::IndexType i = 0; i < nnodesThisZone; i++)
          {
            const auto nodeId = zone.getId(i);
            int *nodePtr = nodeData + nodeId;
            axom::atomicMax<ExecSpace>(nodePtr, nmats);
          }
        }
      });
    axom::synchronize<ExecSpace>();
    AXOM_ANNOTATE_END("nMatsPerNode");

    // Now, mark all zones that have 1 mat per node as clean.
    AXOM_ANNOTATE_BEGIN("mask");
    const auto nzones = m_topologyView.numberOfZones();
    axom::Array<MaskType> mask(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
    auto maskView = mask.view();
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = deviceTopologyView.zone(zoneIndex);

        MaskType clean {1};
        const axom::IndexType nnodesThisZone = zone.numberOfNodes();
        const auto &zoneNodeIds = zone.getIdsStorage();
        for(axom::IndexType i = 0; i < nnodesThisZone; i++)
        {
          const auto nodeId = zoneNodeIds[i];
          clean &= (nMatsPerNodeView[nodeId] == 1) ? MaskType{1} : MaskType{0};
        }

        maskView[zoneIndex] = clean;
      });
    AXOM_ANNOTATE_END("mask");

    axom::IndexType nClean = 0;
    axom::Array<int> maskOffsets;
    axom::ArrayView<int> maskOffsetsView;
    if(nzones > 0)
    {
      AXOM_ANNOTATE_SCOPE("numClean");
      if constexpr(axom::execution_space<ExecSpace>::onDevice())
      {
        // On device, use a reduction on maskView to count clean zones.
        axom::ReduceSum<ExecSpace, int> mask_reduce(0);
        axom::for_all<ExecSpace>(
          nzones,
          AXOM_LAMBDA(axom::IndexType zoneIndex) {
            mask_reduce += static_cast<int>(maskView[zoneIndex]);
          });
        nClean = mask_reduce.get();
      }
      else
      {
        // Off device, do the offset scan early and use the results to compute
        // nClean instead of a reduction.
        AXOM_ANNOTATE_SCOPE("offsets");
        maskOffsets =
          axom::Array<int>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
        maskOffsetsView = maskOffsets.view();
        axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

        nClean = maskOffsetsView[nzones - 1] + static_cast<int>(maskView[nzones - 1]);
      }
    }

    if(nClean > 0)
    {
      // Compute maskOffsets if we did not do it yet.
      if(maskOffsets.empty())
      {
        AXOM_ANNOTATE_SCOPE("offsets");
        maskOffsets =
          axom::Array<int>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
        maskOffsetsView = maskOffsets.view();
        axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      }

      // Make the output cleanIndices array.
      AXOM_ANNOTATE_BEGIN("cleanIndices");
      cleanIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nClean, nClean, allocatorID);
      auto cleanIndicesView = cleanIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            cleanIndicesView[maskOffsetsView[index]] = index;
          }
        });
      AXOM_ANNOTATE_END("cleanIndices");

      // Make the mixedIndices array.
      AXOM_ANNOTATE_BEGIN("mixedIndices");
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          maskView[index] = (maskView[index] == MaskType {1}) ? MaskType {0} : MaskType {1};
        });
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      const int nMixed = nzones - nClean;
      mixedIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nMixed, nMixed, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            mixedIndicesView[maskOffsetsView[index]] = index;
          }
        });
      AXOM_ANNOTATE_END("mixedIndices");
    }
    else
    {
      AXOM_ANNOTATE_SCOPE("mixedIndices");
      cleanIndices = axom::Array<axom::IndexType>();

      // There were no clean, so it must all be mixed.
      mixedIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) { mixedIndicesView[index] = index; });
    }
  }

  /*!
   * \brief Build the list of clean and mixed zones using the number of materials
   *        per zone, maxed to the nodes. Limit the number of zones.
   *
   * \param nnodes The number of nodes in the topology's coordset.
   * \param selectedZonesView A view containing the zone indices we're considering. We pass
   *                          it by value so it can be captured into lambdas.
   * \param[out] cleanIndices An array that will contain the list of clean material zone ids.
   * \param[out] mixedIndices An array that will contain the list of mixed material zone ids.
   *
   * \note The clean/mixed index arrays are not strictly what could be determined by the matset alone.
   *       We figure out which nodes touch multiple materials. Then we iterate over the zones and
   *       those that touch only nodes that have 1 material are marked clean, otherwise they are
   *       considered mixed as we might have to split those zones.
   */
  void execute(axom::IndexType nnodes,
               const SelectedZonesView selectedZonesView,
               axom::Array<axom::IndexType> &cleanIndices,
               axom::Array<axom::IndexType> &mixedIndices) const
  {
    AXOM_ANNOTATE_SCOPE("ZoneListBuilder.2");
    SLIC_ASSERT(selectedZonesView.size() > 0);

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_BEGIN("nMatsPerNode");
    axom::Array<int> nMatsPerNode(axom::ArrayOptions::Uninitialized(), nnodes, nnodes, allocatorID);
    auto nMatsPerNodeView = nMatsPerNode.view();
    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType nodeIndex) { nMatsPerNodeView[nodeIndex] = 1; });

    // Determine max number of materials a node might touch.
    MatsetView deviceMatsetView(m_matsetView);
    const TopologyView deviceTopologyView(m_topologyView);
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const int nmats = deviceMatsetView.numberOfMaterials(zoneIndex);
        if(nmats > 1)
        {
          const auto zone = deviceTopologyView.zone(zoneIndex);
          const auto nnodesThisZone = zone.numberOfNodes();
          int *nodeData = nMatsPerNodeView.data();
          for(axom::IndexType i = 0; i < nnodesThisZone; i++)
          {
            const auto nodeId = zone.getId(i);
            int *nodePtr = nodeData + nodeId;
            axom::atomicMax<ExecSpace>(nodePtr, nmats);
          }
        }
      });
    AXOM_ANNOTATE_END("nMatsPerNode");

    // Now, mark all selected zones that have 1 mat per node as clean.
    AXOM_ANNOTATE_BEGIN("mask");
    const auto nzones = selectedZonesView.size();
    axom::Array<MaskType> mask(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
    auto maskView = mask.view();
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto zone = deviceTopologyView.zone(zoneIndex);

        MaskType clean {1};
        const axom::IndexType nnodesThisZone = zone.numberOfNodes();
        for(axom::IndexType i = 0; i < nnodesThisZone; i++)
        {
          const auto nodeId = zone.getId(i);
          clean &= (nMatsPerNodeView[nodeId] == 1) ? MaskType {1} : MaskType {0};
        }

        maskView[szIndex] = clean;
      });
    AXOM_ANNOTATE_END("mask");

    axom::IndexType nClean = 0;
    axom::Array<int> maskOffsets;
    axom::ArrayView<int> maskOffsetsView;
    if(nzones > 0)
    {
      AXOM_ANNOTATE_SCOPE("numClean");
      if constexpr(axom::execution_space<ExecSpace>::onDevice())
      {
        // On device, use a reduction on maskView to count clean zones.
        axom::ReduceSum<ExecSpace, int> mask_reduce(0);
        axom::for_all<ExecSpace>(
          nzones,
          AXOM_LAMBDA(axom::IndexType szIndex) {
            mask_reduce += static_cast<int>(maskView[szIndex]);
          });
        nClean = mask_reduce.get();
      }
      else
      {
        // Off device, do the offset scan early and use the results to compute
        // nClean instead of a reduction.
        AXOM_ANNOTATE_SCOPE("offsets");
        maskOffsets =
          axom::Array<int>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
        maskOffsetsView = maskOffsets.view();
        axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

        nClean = maskOffsetsView[nzones - 1] + static_cast<int>(maskView[nzones - 1]);
      }
    }

    if(nClean > 0)
    {
      // Compute maskOffsets if we did not do it yet.
      if(maskOffsets.empty())
      {
        AXOM_ANNOTATE_SCOPE("offsets");
        maskOffsets =
          axom::Array<int>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
        maskOffsetsView = maskOffsets.view();
        axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      }

      // Make the output cleanIndices array.
      AXOM_ANNOTATE_BEGIN("cleanIndices");
      cleanIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nClean, nClean, allocatorID);
      auto cleanIndicesView = cleanIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            cleanIndicesView[maskOffsetsView[index]] = selectedZonesView[index];
          }
        });
      AXOM_ANNOTATE_END("cleanIndices");

      // Make the mixedIndices array.
      AXOM_ANNOTATE_BEGIN("mixedIndices");
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          maskView[index] = (maskView[index] == MaskType {1}) ? MaskType {0} : MaskType {1};
        });
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      const int nMixed = nzones - nClean;
      mixedIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nMixed, nMixed, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            mixedIndicesView[maskOffsetsView[index]] = selectedZonesView[index];
          }
        });
      AXOM_ANNOTATE_END("mixedIndices");
    }
    else
    {
      AXOM_ANNOTATE_SCOPE("mixedIndices");
      cleanIndices = axom::Array<axom::IndexType>();

      // There were no clean, so it must all be mixed.
      mixedIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) { mixedIndicesView[index] = selectedZonesView[index]; });
    }
  }

  /*!
   * \brief Build the list of clean and mixed zones using the number of materials
   *        per zone. This method essentially partitions the input selectedZonesView
   *        into clean and mixed lists.
   *
   * \param selectedZonesView A view containing the zone indices we're considering. We pass
   *                          it by value so it can be captured into lambdas.
   * \param[out] cleanIndices An array that will contain the list of clean material zone ids.
   * \param[out] mixedIndices An array that will contain the list of mixed material zone ids.
   *
   */
  void execute(const SelectedZonesView selectedZonesView,
               axom::Array<axom::IndexType> &cleanIndices,
               axom::Array<axom::IndexType> &mixedIndices) const
  {
    AXOM_ANNOTATE_SCOPE("ZoneListBuilder.3");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_BEGIN("mask");
    const auto nzones = selectedZonesView.size();
    axom::Array<MaskType> mask(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
    auto maskView = mask.view();
    axom::ReduceSum<ExecSpace, int> mask_reduce(0);
    const MatsetView deviceMatsetView(m_matsetView);
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto matZoneIndex = zoneIndex;

        // clean zone == 1, mixed zone = 0
        const int ival = (deviceMatsetView.numberOfMaterials(matZoneIndex) == 1) ? 1 : 0;
        maskView[szIndex] = static_cast<MaskType>(ival);
        mask_reduce += ival;
      });
    AXOM_ANNOTATE_END("mask");

    const axom::IndexType numCleanZones = mask_reduce.get();
    const axom::IndexType numMixedZones = nzones - numCleanZones;
    if(numCleanZones > 0 && numMixedZones > 0)
    {
      AXOM_ANNOTATE_SCOPE("mixedIndices");

      axom::Array<int> maskOffsets(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
      auto maskOffsetsView = maskOffsets.view();
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

      // Fill in clean zone ids.
      cleanIndices = axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(),
                                                  numCleanZones,
                                                  numCleanZones,
                                                  allocatorID);
      auto cleanIndicesView = cleanIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType szIndex) {
          if(maskView[szIndex] > 0)
          {
            cleanIndicesView[maskOffsetsView[szIndex]] = selectedZonesView[szIndex];
          }
          maskView[szIndex] = (maskView[szIndex] > 0) ? MaskType {0} : MaskType {1};
        });

      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

      // Fill in mixed zone ids.
      mixedIndices = axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(),
                                                  numMixedZones,
                                                  numMixedZones,
                                                  allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType szIndex) {
          if(maskView[szIndex] > 0)
          {
            mixedIndicesView[maskOffsetsView[szIndex]] = selectedZonesView[szIndex];
          }
        });
    }
    else if(numCleanZones > 0)
    {
      AXOM_ANNOTATE_SCOPE("cleanIndices");

      // There were no mixed, so it must all be clean.
      cleanIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
      auto cleanIndicesView = cleanIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) { cleanIndicesView[index] = selectedZonesView[index]; });

      mixedIndices = axom::Array<axom::IndexType>();
    }
    else if(numMixedZones > 0)
    {
      AXOM_ANNOTATE_SCOPE("mixedIndices");

      cleanIndices = axom::Array<axom::IndexType>();

      // There were no clean, so it must all be mixed.
      mixedIndices =
        axom::Array<axom::IndexType>(axom::ArrayOptions::Uninitialized(), nzones, nzones, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) { mixedIndicesView[index] = selectedZonesView[index]; });
    }
  }

private:
  TopologyView m_topologyView;
  MatsetView m_matsetView;
};

}  // end namespace bump
}  // end namespace axom

#endif
