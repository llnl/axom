// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_EXTRACT_ZONES_AND_MATSET_POLYHEDRAL_HPP
#define AXOM_BUMP_EXTRACT_ZONES_AND_MATSET_POLYHEDRAL_HPP

#include "axom/core.hpp"
#include "axom/bump/ExtractZones.hpp"
#include "axom/bump/CoordsetBlender.hpp"
#include "axom/bump/CoordsetSlicer.hpp"
#include "axom/bump/FieldSlicer.hpp"
#include "axom/bump/MatsetSlicer.hpp"

namespace axom
{
namespace bump
{

/*!
 * \brief Make a new polyhedral topology, coordset, and matset by extracting selected
 *        zones from the input 3D structured mesh. We make assumptions for structured 3D
 *        meshes that let us make unique faces by construction.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam IndexPolicy The structured indexing policy used in the topology.
 * \tparam CoordsetView The coordset view type.
 * \tparam MatsetView The matset view type.
 */
template <typename ExecSpace, typename IndexPolicy, typename CoordsetView, typename MatsetView>
class ExtractZonesAndMatsetPolyhedral
  : public ExtractZonesAndMatset<ExecSpace, axom::bump::views::StructuredTopologyView<IndexPolicy>, CoordsetView, MatsetView>
{
public:
  using ParentClass =
    ExtractZonesAndMatset<ExecSpace, axom::bump::views::StructuredTopologyView<IndexPolicy>, CoordsetView, MatsetView>;
  using Sizes = typename ParentClass::Sizes;
  using SelectedZonesView = typename ParentClass::SelectedZonesView;
  using TopologyView = axom::bump::views::StructuredTopologyView<IndexPolicy>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

  static_assert(TopologyView::dimension() == 3, "This class requires 3D structured topology views");

  /*!
   * \brief Constructor
   *
   * \param topoView The input topology view.
   * \param coordsetView The input coordset view.
   * \param matsetView The input matset view.
   */
  ExtractZonesAndMatsetPolyhedral(const TopologyView &topoView,
                                  const CoordsetView &coordsetView,
                                  const MatsetView &matsetView)
    : ParentClass(topoView, coordsetView, matsetView)
  { }

  /*!
   * \brief Destructor
   */
  virtual ~ExtractZonesAndMatsetPolyhedral() = default;

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  /*!
   * \brief Make the output topology for just the selected zones. This is an
   *        override of the base class' behavior and we make assumptions for
   *        3D structured topologies to help make faces more easily.
   *
   * \param selectedZonesView A view that contains the ids of the zones to extract.
   * \param dataSizes Array sizes for connectivity, sizes, etc.
   * \param extra Extra sizes for connectivity, sizes, etc.
   * \param old2newView A view that lets us map old node numbers to new node numbers.
   * \param n_topo The input topology.
   * \param n_options A node containing the options (ignored).
   * \param n_newTopo A node to contain the new topology.
   */
  virtual void makeTopology(const SelectedZonesView selectedZonesView,
                            const Sizes &AXOM_UNUSED_PARAM(dataSizes),
                            const Sizes &AXOM_UNUSED_PARAM(extra),
                            const axom::ArrayView<ConnectivityType> &old2newView,
                            const conduit::Node &n_topo,
                            const conduit::Node &AXOM_UNUSED_PARAM(n_options),
                            conduit::Node &n_newTopo) const override
  {
    AXOM_ANNOTATE_SCOPE("makeTopology(polyhedral)");
    namespace utils = axom::bump::utilities;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // We know that we have a 3D structured mesh for which we need to make a
    // polyhedral output topology.

    constexpr int FacesPerHex = 6;
    // Note: MSVC needs constexpr lambda capture to be marked `static` even though constexpr should suffice
    static constexpr int PointsPerQuad = 4;
    constexpr std::uint8_t ZoneEmpty = 0;
    constexpr std::uint8_t ZoneSelected = 1 << 7;

    const auto numSelectedZones = selectedZonesView.size();
    const auto nzones = ParentClass::m_topologyView.numberOfZones();

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("select");
    // Expand the selected zones over all the zones in the mesh.
    axom::Array<std::uint8_t> zoneFlags(nzones, nzones, allocatorID);
    auto zoneFlagsView = zoneFlags.view();
    zoneFlags.fill(ZoneEmpty);
    axom::for_all<ExecSpace>(
      numSelectedZones,
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        zoneFlagsView[zoneIndex] = ZoneSelected;
      });
    AXOM_ANNOTATE_END("select");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("identify");

    // Figure out which zones need to make faces. We know the input mesh was structured
    // so we can iterate over the selected zones and check their face neighbors to see
    // which faces we need to make.

    axom::Array<int> zoneFaceSizes(nzones, nzones, allocatorID);
    zoneFaceSizes.fill(0);
    auto zoneFaceSizesView = zoneFaceSizes.view();

    const auto deviceTopologyView(ParentClass::m_topologyView);
    axom::ReduceSum<ExecSpace, int> faceCount_reduce(0);
    axom::for_all<ExecSpace>(
      numSelectedZones,
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto logicalIndex = deviceTopologyView.indexing().IndexToLogicalIndex(zoneIndex);

        int faceCount = 0;
        std::uint8_t zf = ZoneEmpty;
        const int faceOrder[] = {0, 2, 4, 1, 3, 5};
        for(int fi = 0; fi < FacesPerHex; fi++)
        {
          const int queryFace = faceOrder[fi];

          // x y z x y z
          const int dim = queryFace >> 1;
          const int delta = (fi < 3) ? -1 : 1;

          bool makeFace = false;
          if(fi < 3)
          {
            // The zone always makes the lower faces.
            makeFace = true;
          }
          else if(logicalIndex[dim] == (deviceTopologyView.indexing().logicalDimensions()[dim] - 1))
          {
            // The zone makes the high faces if on the edge of the mesh.
            makeFace = true;
          }
          else
          {
            // Get the index of the face neighbor.
            auto logicalNeighbor(logicalIndex);
            logicalNeighbor[dim] += delta;
            const auto neighbor = deviceTopologyView.indexing().LogicalIndexToIndex(logicalNeighbor);

            // If that neighbor is NOT selected, we need to make the face.
            makeFace = (zoneFlagsView[neighbor] & ZoneSelected) == 0;
          }

          zf |= makeFace ? (1 << queryFace) : 0;
          faceCount += makeFace ? 1 : 0;
        }

        zoneFlagsView[zoneIndex] |= zf;
        zoneFaceSizesView[zoneIndex] = faceCount;
        faceCount_reduce += faceCount;
      });
    const int faceCount = faceCount_reduce.get();
    AXOM_ANNOTATE_END("identify");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("offsets");
    axom::Array<int> zoneFaceOffsets(nzones, nzones, allocatorID);
    auto zoneFaceOffsetsView = zoneFaceOffsets.view();
    axom::exclusive_scan<ExecSpace>(zoneFaceSizesView, zoneFaceOffsetsView);
    zoneFaceSizes.clear();
    AXOM_ANNOTATE_END("offsets");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("allocate");

    n_newTopo["type"] = "unstructured";
    n_newTopo["coordset"] = n_topo["coordset"].as_string();
    n_newTopo["elements/shape"] = "polyhedral";
    n_newTopo["subelements/shape"] = "polygonal";

    conduit::Node &n_conn = n_newTopo["elements/connectivity"];
    n_conn.set_allocator(c2a.getConduitAllocatorID());
    n_conn.set(
      conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numSelectedZones * FacesPerHex));
    auto connView = utils::make_array_view<ConnectivityType>(n_conn);

    conduit::Node &n_sizes = n_newTopo["elements/sizes"];
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numSelectedZones));
    auto sizesView = utils::make_array_view<ConnectivityType>(n_sizes);

    conduit::Node &n_offsets = n_newTopo["elements/offsets"];
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numSelectedZones));
    auto offsetsView = utils::make_array_view<ConnectivityType>(n_offsets);

    conduit::Node &n_se_conn = n_newTopo["subelements/connectivity"];
    n_se_conn.set_allocator(c2a.getConduitAllocatorID());
    n_se_conn.set(
      conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, faceCount * PointsPerQuad));
    auto seConnView = utils::make_array_view<ConnectivityType>(n_se_conn);

    conduit::Node &n_se_sizes = n_newTopo["subelements/sizes"];
    n_se_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_se_sizes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, faceCount));
    auto seSizesView = utils::make_array_view<ConnectivityType>(n_se_sizes);

    conduit::Node &n_se_offsets = n_newTopo["subelements/offsets"];
    n_se_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_se_offsets.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, faceCount));
    auto seOffsetsView = utils::make_array_view<ConnectivityType>(n_se_offsets);

    AXOM_ANNOTATE_END("allocate");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("subelements");

    // Fill subelement sizes, offsets
    axom::for_all<ExecSpace>(
      faceCount,
      AXOM_LAMBDA(axom::IndexType index) {
        seSizesView[index] = PointsPerQuad;
        seOffsetsView[index] = PointsPerQuad * index;
      });

    // Fill subelement connectivity to define the faces.
    axom::for_all<ExecSpace>(
      numSelectedZones,
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto zone = deviceTopologyView.zone(zoneIndex);
        const auto zf = zoneFlagsView[zoneIndex];
        int localFace = 0;
        // Query faces in this order so we do all the low faces then all the high faces.
        const int faceOrder[] = {0, 2, 4, 1, 3, 5};
        for(int fi = 0; fi < FacesPerHex; fi++)
        {
          const int queryFace = faceOrder[fi];

          // Check whether this zone defines the face.
          const bool zoneDefinesFace = (zf & (1 << queryFace)) > 0;
          if(zoneDefinesFace)
          {
            // Get the node ids from the zone.
            ConnectivityType ids[PointsPerQuad];
            axom::IndexType numIds = 0;
            zone.getFace(queryFace, ids, numIds);
            SLIC_ASSERT(numIds == PointsPerQuad);

            // Map the face's node ids to the new coordset.
            const auto faceId = zoneFaceOffsetsView[zoneIndex] + localFace;
            const auto offset = faceId * PointsPerQuad;
            for(int idx = 0; idx < numIds; idx++)
            {
              seConnView[offset + idx] = old2newView[ids[idx]];
            }
            localFace++;
          }
        }
      });
    AXOM_ANNOTATE_END("subelements");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("elements");

    // Fill element sizes, offsets
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        sizesView[szIndex] = FacesPerHex;
        offsetsView[szIndex] = FacesPerHex * szIndex;
      });

    // Fill element connectivity
    axom::for_all<ExecSpace>(
      numSelectedZones,
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto zf = zoneFlagsView[zoneIndex];

        int localFace = 0, connOffset = 0;
        // Query faces in this order so we do all the low faces then all the high faces.
        const int faceOrder[] = {0, 2, 4, 1, 3, 5};
        for(int fi = 0; fi < FacesPerHex; fi++)
        {
          const int queryFace = faceOrder[fi];

          const bool zoneDefinesFace = (zf & (1 << queryFace)) > 0;
          if(zoneDefinesFace)
          {
            const ConnectivityType faceId = zoneFaceOffsetsView[zoneIndex] + localFace;
            connView[szIndex * FacesPerHex + connOffset] = faceId;
            localFace++;
            connOffset++;
          }
          else
          {
            // Get the index of the relevant neighbor.
            const int dim = queryFace >> 1;  // x y z x y z
            const int delta = (fi < 3) ? -1 : 1;
            auto logicalNeighbor = deviceTopologyView.indexing().IndexToLogicalIndex(zoneIndex);
            logicalNeighbor[dim] += delta;
            const auto neighbor = deviceTopologyView.indexing().LogicalIndexToIndex(logicalNeighbor);

            // Get whether the neighbor defines the face companion.
            const auto nzf = zoneFlagsView[neighbor];
            const int neighborQueryFace = (fi < 3) ? (queryFace + 1) : (queryFace - 1);
            const bool neighborDefinesFace = (nzf & (1 << neighborQueryFace)) > 0;

            if(neighborDefinesFace)
            {
              const ConnectivityType faceId = zoneFaceOffsetsView[neighbor] + dim;
              connView[szIndex * FacesPerHex + connOffset] = faceId;
              connOffset++;
            }
          }
        }
      });
    AXOM_ANNOTATE_END("elements");
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
