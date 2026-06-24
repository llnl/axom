// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file MarchingCubesBumpAdaptor.hpp
 *
 * @brief Adapts a bump::extraction::CutField Blueprint output mesh into the
 * legacy quest::MarchingCubes fixed-stride output buffers.
 *
 * bump's CutField output is a *welded*, *mixed-shape* unstructured Blueprint
 * topology:
 *   topologies/<t>/type                == "unstructured"
 *   topologies/<t>/elements/connectivity (flat, ConnectivityType)
 *   topologies/<t>/elements/sizes        (per-zone corner count)
 *   topologies/<t>/elements/offsets      (per-zone start into connectivity)
 *   topologies/<t>/elements/shapes       (per-zone Blueprint ShapeID)
 *   coordsets/<c>/values/{x,y[,z]}       (explicit, blended/welded points)
 *   fields/originalElements/values       (element-assoc, input zone per fragment)
 *
 * The legacy MarchingCubes output is an *un-welded* fixed-stride representation:
 *   m_facetNodeCoords : (facetCount*DIM, DIM)   one row per facet-corner
 *   m_facetNodeIds    : (facetCount, DIM)        DIM corner ids per facet, where
 *                       the ids index into m_facetNodeCoords and are offset by
 *                       m_facetIndexOffset*DIM (the parent concatenates domains)
 *   m_facetParentIds  : (facetCount)             parent-cell id per facet
 *
 * Conversion, all in ExecSpace memory:
 *   1. For DIM==2: each welded segment (Line_ShapeID, size 2) is one facet.
 *      For DIM==3: each welded polygon of p corners fan-triangulates into
 *      (p-2) triangles (corners {0,k,k+1}).
 *   2. Re-expand: every output facet writes DIM fresh rows into
 *      m_facetNodeCoords (so the legacy facetCount*DIM node-count invariant
 *      holds) and its m_facetNodeIds row is the consecutive ids of those rows.
 *   3. Parent id per facet := originalElements[srcZone], optionally remapped to
 *      the legacy field-stride flat order (structured input + legacyFieldOrder).
 *
 * This file performs no triangle *welding* of its own; it intentionally expands
 * back to the legacy soup so existing users are byte-compatible.  Callers who
 * want bump's richer welded mesh use the additive Blueprint accessors (Phase 3).
 */

#ifndef AXOM_QUEST_MARCHINGCUBESBUMPADAPTOR_H_
#define AXOM_QUEST_MARCHINGCUBESBUMPADAPTOR_H_

#include "axom/config.hpp"

#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)

  #include "axom/core/execution/execution_space.hpp"
  #include "axom/core/execution/for_all.hpp"
  #include "axom/core/memory_management.hpp"
  #include "axom/core/Array.hpp"
  #include "axom/core/ArrayView.hpp"
  #include "axom/core/MDMapping.hpp"
  #include "axom/core/numerics/floating_point_limits.hpp"
  #include "axom/slic/interface/slic_macros.hpp"

  #include "axom/bump/utilities/blueprint_utilities.hpp"
  #include "axom/bump/views/NodeArrayView.hpp"
  #include "axom/bump/views/Shapes.hpp"

  #include "conduit_node.hpp"

  #include <string>

namespace axom
{
namespace quest
{
namespace detail
{
namespace marching_cubes
{
/*!
 * @brief Number of legacy facets a bump zone of \a nCorners contributes.
 *
 * 2D: segment -> 1.  3D: p-gon fans into (p-2) triangles (0 if degenerate).
 */
template <int DIM>
AXOM_HOST_DEVICE inline axom::IndexType facetsPerZone(axom::IndexType nCorners)
{
  if(DIM == 3)
  {
    return nCorners >= 3 ? (nCorners - 2) : 0;
  }
  // DIM == 2: a line segment.
  return nCorners >= 2 ? 1 : 0;
}

/*!
 * @brief Convert one bump CutField output (single domain) into the legacy
 * fixed-stride output buffers supplied by the parent MarchingCubes.
 *
 * @tparam DIM Spatial dimension (2 or 3).
 * @tparam ExecSpace Axom execution space.
 *
 * @param n_output The bump CutField output Blueprint mesh (in ExecSpace memory).
 * @param facetNodeIds [out] view, shape (totalFacetCount, DIM).
 * @param facetNodeCoords [out] view, shape (totalFacetCount*DIM, DIM).
 * @param facetParentIds [out] view, shape (totalFacetCount).
 * @param facetIndexOffset This domain's first facet index in the concatenated
 *   output (the parent's m_facetIndexOffsets[d]).  Node-coord rows and node ids
 *   are written starting at facetIndexOffset (ids offset by *DIM).
 * @param thisDomainFacetCount Number of facets this domain produces (already
 *   computed by the caller; equals sum of facetsPerZone over the bump zones).
 * @param fieldStrideRemap If non-null, a precomputed per-input-zone map from
 *   bump's i-fastest Blueprint zone id to the legacy field-stride flat id.  When
 *   null, originalElements ids are written through unchanged.
 *
 * @pre All output views and \a n_output live in ExecSpace's memory space.
 */
template <int DIM, typename ExecSpace>
void adaptCutFieldOutput(const conduit::Node& n_output,
                         axom::ArrayView<axom::IndexType, 2> facetNodeIds,
                         axom::ArrayView<double, 2> facetNodeCoords,
                         axom::ArrayView<axom::IndexType, 1> facetParentIds,
                         axom::IndexType facetIndexOffset,
                         axom::IndexType thisDomainFacetCount,
                         axom::ArrayView<const axom::IndexType, 1> fieldStrideRemap)
{
  namespace bputils = axom::bump::utilities;
  namespace bpviews = axom::bump::views;

  if(thisDomainFacetCount == 0)
  {
    return;
  }

  // --- Locate the single output topology + coordset ------------------------
  const conduit::Node& n_topos = n_output.fetch_existing("topologies");
  SLIC_ASSERT(n_topos.number_of_children() == 1);
  const conduit::Node& n_topo = n_topos.child(0);
  const std::string coordsetName = n_topo.fetch_existing("coordset").as_string();
  const conduit::Node& n_coords =
    n_output.fetch_existing(axom::fmt::format("coordsets/{}", coordsetName));
  const conduit::Node& n_elems = n_topo.fetch_existing("elements");

  // bump always emits explicit sizes/offsets/connectivity for cut output.
  const conduit::Node& n_sizes = n_elems.fetch_existing("sizes");
  const conduit::Node& n_offsets = n_elems.fetch_existing("offsets");
  const conduit::Node& n_conn = n_elems.fetch_existing("connectivity");

  // originalElements: element-associated, one entry per output zone (fragment).
  const conduit::Node& n_orig = n_output.fetch_existing("fields/originalElements/values");

  SLIC_ERROR_IF(n_offsets.dtype().id() != n_sizes.dtype().id() ||
                  n_conn.dtype().id() != n_sizes.dtype().id() ||
                  n_orig.dtype().id() != n_sizes.dtype().id(),
                "MarchingCubes bump adaptor expects connectivity, sizes, "
                "offsets, and originalElements to use the same integer type.");

  auto adaptViews = [&](auto sizesView, auto offsetsView, auto connView, auto origView) {
    const conduit::Node& n_x = n_coords.fetch_existing("values/x");
    const conduit::Node& n_y = n_coords.fetch_existing("values/y");
    auto xView = bputils::make_array_view<double>(n_x);
    auto yView = bputils::make_array_view<double>(n_y);
    // z only in 3D.
    axom::ArrayView<double> zView;
    if(DIM == 3)
    {
      const conduit::Node& n_z = n_coords.fetch_existing("values/z");
      zView = bputils::make_array_view<double>(n_z);
    }

    const axom::IndexType numZones = static_cast<axom::IndexType>(sizesView.size());

    // --- Per-zone facet offset (exclusive scan of facetsPerZone) -----------
    // We need, for each bump zone, the index of its first facet within this
    // domain so kernels can write without atomics.
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<axom::IndexType> zoneFacetCounts(numZones, numZones, allocatorID);
    auto zoneFacetCountsView = zoneFacetCounts.view();
    axom::for_all<ExecSpace>(
      numZones,
      AXOM_LAMBDA(axom::IndexType z) {
        zoneFacetCountsView[z] = facetsPerZone<DIM>(static_cast<axom::IndexType>(sizesView[z]));
      });

    axom::Array<axom::IndexType> zoneFacetOffsets(numZones, numZones, allocatorID);
    auto zoneFacetOffsetsView = zoneFacetOffsets.view();
    axom::exclusive_scan<ExecSpace>(zoneFacetCountsView, zoneFacetOffsetsView);

    // Capture raw views for the kernel.
    const bool doRemap = !fieldStrideRemap.empty();

    // --- The fan-triangulation + re-expansion kernel -----------------------
    // One thread per bump zone.  Each zone writes facetsPerZone facets; for
    // each facet we emit DIM corner coords (re-expanded / un-welded) and DIM
    // ids.
    axom::for_all<ExecSpace>(
      numZones,
      AXOM_LAMBDA(axom::IndexType z) {
        const axom::IndexType nCorners = static_cast<axom::IndexType>(sizesView[z]);
        const axom::IndexType nFacets = facetsPerZone<DIM>(nCorners);
        if(nFacets == 0)
        {
          return;
        }
        const axom::IndexType connStart = static_cast<axom::IndexType>(offsetsView[z]);

        // Parent-cell id for every facet of this zone.
        axom::IndexType parentId = static_cast<axom::IndexType>(origView[z]);
        if(doRemap)
        {
          parentId = fieldStrideRemap[parentId];
        }

        // This zone's first facet within the whole concatenated output.
        const axom::IndexType facetBase = facetIndexOffset + zoneFacetOffsetsView[z];

        for(axom::IndexType f = 0; f < nFacets; ++f)
        {
          const axom::IndexType facetIdx = facetBase + f;

          // Local corner indices of this facet within the zone.
          //   DIM==2: the segment endpoints {0,1}
          //   DIM==3: fan triangle {0, f+1, f+2}
          axom::IndexType local[DIM];
          if constexpr(DIM == 3)
          {
            local[0] = 0;
            local[1] = f + 1;
            local[2] = f + 2;
          }
          else
          {
            local[0] = 0;
            local[1] = 1;
          }

          // Re-expanded node rows for this facet are contiguous: each facet
          // owns exactly DIM rows in m_facetNodeCoords at facetIdx*DIM .. +DIM.
          const axom::IndexType nodeRowBase = facetIdx * DIM;

          for(int c = 0; c < DIM; ++c)
          {
            const axom::IndexType weldedNode =
              static_cast<axom::IndexType>(connView[connStart + local[c]]);
            const axom::IndexType outRow = nodeRowBase + c;

            facetNodeCoords(outRow, 0) = xView[weldedNode];
            facetNodeCoords(outRow, 1) = yView[weldedNode];
            if(DIM == 3)
            {
              facetNodeCoords(outRow, 2) = zView[weldedNode];
            }

            // Legacy ids index into m_facetNodeCoords directly.
            facetNodeIds(facetIdx, c) = outRow;
          }

          facetParentIds[facetIdx] = parentId;
        }
      });
  };

  #if defined(_WIN32)
  adaptViews(bputils::make_array_view<axom::IndexType>(n_sizes),
             bputils::make_array_view<axom::IndexType>(n_offsets),
             bputils::make_array_view<axom::IndexType>(n_conn),
             bputils::make_array_view<axom::IndexType>(n_orig));
  #else
  bpviews::indexNodeToArrayViewSame(n_sizes, n_offsets, n_conn, n_orig, adaptViews);
  #endif
}

/*!
 * @brief Build the per-input-zone remap from bump's i-fastest Blueprint zone id
 * to the legacy field-stride flat id, for structured input.
 *
 * bump's StructuredIndexing numbers zones flat = i + j*nx + k*nx*ny (i-fastest),
 * independent of memory layout.  The legacy parent-cell id is the flat index in
 * the function field's stride order.  This routine, given the per-dimension cell
 * counts \a cellDims (logical, in i,j,k order) and the function field's
 * \a fieldSlowestDirs (the slowest-to-fastest permutation from the field's
 * MDMapping), produces remap[bumpZoneId] = legacyZoneId.
 *
 * Returns an empty Array when the two orderings coincide (i-fastest field), so
 * callers can skip remapping entirely.
 *
 * @note Built in host memory then copied to ExecSpace memory by the caller.
 */
template <int DIM>
axom::Array<axom::IndexType> buildFieldStrideRemap(
  const axom::StackArray<axom::IndexType, DIM>& cellDims,
  const axom::StackArray<std::uint16_t, DIM>& fieldSlowestDirs)
{
  // Identity stride order is "i fastest" == slowestDirs {DIM-1, ..., 1, 0}.
  bool isIFastest = true;
  for(int d = 0; d < DIM; ++d)
  {
    if(fieldSlowestDirs[d] != static_cast<std::uint16_t>(DIM - 1 - d))
    {
      isIFastest = false;
      break;
    }
  }
  if(isIFastest)
  {
    return axom::Array<axom::IndexType>(0, 0);  // no remap needed
  }

  axom::IndexType numZones = 1;
  for(int d = 0; d < DIM; ++d)
  {
    numZones *= cellDims[d];
  }

  // bump mapping: i-fastest.
  axom::MDMapping<DIM> bumpMap(cellDims, axom::ArrayStrideOrder::COLUMN);
  // legacy mapping: field stride order via slowestDirs.
  axom::MDMapping<DIM> legacyMap;
  legacyMap.initializeShape(cellDims, fieldSlowestDirs);

  axom::Array<axom::IndexType> remap(numZones, numZones);
  auto remapView = remap.view();
  // Host loop: enumerate logical multi-indices, map each to both flat ids.
  for(axom::IndexType bumpId = 0; bumpId < numZones; ++bumpId)
  {
    const auto multi = bumpMap.toMultiIndex(bumpId);
    remapView[bumpId] = legacyMap.toFlatIndex(multi);
  }
  return remap;
}

}  // namespace marching_cubes
}  // namespace detail
}  // namespace quest
}  // namespace axom

#endif  // AXOM_USE_CONDUIT && AXOM_USE_BUMP
#endif  // AXOM_QUEST_MARCHINGCUBESBUMPADAPTOR_H_
