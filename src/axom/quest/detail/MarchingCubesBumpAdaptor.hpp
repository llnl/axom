// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file MarchingCubesBumpAdaptor.hpp
 *
 * @brief Adapts a bump::extraction::CutField Blueprint output mesh into
 * the legacy quest::MarchingCubes fixed-stride output buffers.
 *
 * bump's CutField output is a welded, mixed-shape unstructured Blueprint topology:
 *   <bp_root>
 *    ├── topologies
 *    │   └── <t>
 *    │        ├─•  type                == "unstructured"
 *    │        └──  elements
 *    │             ├─• connectivity   (flat, ConnectivityType)
 *    │             ├─• sizes          (per-zone corner count)
 *    │             ├─• offsets        (per-zone start into connectivity)
 *    │             └─• shapes         (per-zone Blueprint ShapeID)
 *    ├── coordsets
 *    │   └── <c>
 *    │        └──  values             (explicit, blended/welded points)
 *    │             ├─•  x
 *    │             ├─•  y
 *    │             └─• [z]
 *    └── fields
 *        └──  originalElements
 *             └─• values              (element-assoc, input zone per fragment)
 *
 * The legacy MarchingCubes output is composed of triangles (3D) or segments (2D):
 *   m_facetNodeCoords : (nodeCount, DIM)         vertices of the mesh
 *   m_facetNodeIds    : (facetCount, DIM)        indices of each facet, where
 *                       the ids index into m_facetNodeCoords and are offset by
 *                       the domain's nodeIndexOffset (the parent concatenates domains)
 *   m_facetParentIds  : (facetCount)             parent-cell id per facet
 *
 * Conversion, all in ExecSpace memory:
 *   1. For DIM==2: each welded segment (Line_ShapeID, size 2) is one facet.
 *      For DIM==3: each welded polygon of p corners fan-triangulates into (p-2) triangles (corners {0,k,k+1}).
 *   2. Reuse bump's welded vertex coordinates and write only triangle/segment connectivity.
 *   3. Parent id per facet := originalElements[srcZone].
 */

#pragma once

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesBumpAdaptor.hpp requires conduit"
#endif
#ifndef AXOM_USE_BUMP
  #error "MarchingCubesBumpAdaptor.hpp requires bump"
#endif

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/reductions.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/numerics/floating_point_limits.hpp"
#include "axom/slic/interface/slic_macros.hpp"

#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/views/Shapes.hpp"

#include "conduit_node.hpp"

#include <string>
#include <type_traits>
#include <utility>

namespace axom::quest::detail::marching_cubes
{
/*!
 * @brief Private name for the parent-zone field requested from bump.
 *
 * Not bump's default: TableBasedExtractor::makeOriginalElements branches on
 * whether the input mesh already carries a field of the configured name
 * and, if so, maps those values forward instead of writing zone indices.
 * Any mesh produced by a prior bump operation contains that field,
 * and the empty "fields" option does not suppress the branch,
 * so a plausible input silently redefines what a parent cell id means.
 */
constexpr const char* kOriginalElementsField = "__axom_mc_originalElements";

/*!
 * @brief Name the parent-zone field carries on the PUBLIC Blueprint output.
 *
 * populateContourMeshBlueprint() hands the mesh to the caller and
 * quest_marching_cubes_bump.cpp asserts on this name, so it is an API contract.
 * bump's output is renamed from the private request name to this immediately
 * after extraction, keeping the private name confined to the request.
 */
constexpr const char* kPublicOriginalElementsField = "originalElements";

/*!
 * @brief Number of legacy facets a bump zone of \a nCorners contributes.
 *
 * 2D: segment -> 1.
 * 3D: p-gon fans into (p-2) triangles (0 if degenerate).
 */
template <int DIM>
AXOM_HOST_DEVICE inline axom::IndexType facetsPerZone(axom::IndexType nCorners)
{
  if constexpr(DIM == 3)
  {
    return nCorners >= 3 ? (nCorners - 2) : 0;
  }
  // DIM == 2: a line segment.
  return nCorners >= 2 ? 1 : 0;
}

template <typename ExecSpace, typename InValuesView, typename OutValuesView, typename CountsView, typename OffsetsView>
void duplicateElementValuesForTriangulationViews(InValuesView inValues,
                                                 OutValuesView outValues,
                                                 CountsView zoneFacetCounts,
                                                 OffsetsView zoneFacetOffsets)
{
  const axom::IndexType inputZoneCount = static_cast<axom::IndexType>(inValues.size());
  axom::for_all<ExecSpace>(
    inputZoneCount,
    AXOM_LAMBDA(axom::IndexType z) {
      const axom::IndexType outBegin = zoneFacetOffsets[z];
      const axom::IndexType nOut = zoneFacetCounts[z];
      for(axom::IndexType f = 0; f < nOut; ++f)
      {
        outValues[outBegin + f] = inValues[z];
      }
    });
}

template <typename ExecSpace, typename CountsView, typename OffsetsView>
void duplicateElementValuesForTriangulation(conduit::Node& n_values,
                                            axom::IndexType outputZoneCount,
                                            CountsView zoneFacetCounts,
                                            OffsetsView zoneFacetOffsets,
                                            int allocatorID)
{
  namespace bpviews = axom::bump::views;

  if(n_values.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < n_values.number_of_children(); ++i)
    {
      duplicateElementValuesForTriangulation<ExecSpace>(n_values[i],
                                                        outputZoneCount,
                                                        zoneFacetCounts,
                                                        zoneFacetOffsets,
                                                        allocatorID);
    }
    return;
  }

  conduit::Node newValues;
  newValues.set_allocator(axom::sidre::ConduitMemory::axomAllocIdToConduit(allocatorID));
  newValues.set(conduit::DataType(n_values.dtype().id(), outputZoneCount));

  bpviews::nodeToArrayViewSame(n_values, newValues, [&](auto inValues, auto outValues) {
    duplicateElementValuesForTriangulationViews<ExecSpace>(inValues,
                                                           outValues,
                                                           zoneFacetCounts,
                                                           zoneFacetOffsets);
  });

  n_values.move(newValues);
}

template <int DIM, typename ExecSpace, typename SizesView, typename OffsetsView, typename ConnView>
void triangulateBlueprintMeshViews(conduit::Node& n_output,
                                   conduit::Node& n_conn,
                                   conduit::Node& n_sizes,
                                   conduit::Node& n_offsets,
                                   const std::string& topologyName,
                                   SizesView sizesView,
                                   OffsetsView offsetsView,
                                   ConnView connView,
                                   int allocatorID)
{
  namespace bputils = axom::bump::utilities;
  namespace bpviews = axom::bump::views;

  using ConnectivityType = typename std::decay_t<ConnView>::value_type;

  const axom::IndexType inputZoneCount = static_cast<axom::IndexType>(sizesView.size());
  if(inputZoneCount == 0)
  {
    return;
  }

  axom::Array<axom::IndexType> zoneFacetCounts(inputZoneCount, inputZoneCount, allocatorID);
  auto zoneFacetCountsView = zoneFacetCounts.view();

  axom::ReduceSum<ExecSpace, axom::IndexType> totalFacetsReduce(0);
  axom::ReduceSum<ExecSpace, axom::IndexType> nonTriReduce(0);
  axom::for_all<ExecSpace>(
    inputZoneCount,
    AXOM_LAMBDA(axom::IndexType z) {
      const auto nCorners = static_cast<axom::IndexType>(sizesView[z]);
      const auto nFacets = facetsPerZone<3>(nCorners);
      zoneFacetCountsView[z] = nFacets;
      totalFacetsReduce += nFacets;
      nonTriReduce += (nCorners == 3) ? 0 : 1;
    });

  const axom::IndexType outputZoneCount = totalFacetsReduce.get();
  if(nonTriReduce.get() == 0)
  {
    return;
  }

  axom::Array<axom::IndexType> zoneFacetOffsets(inputZoneCount, inputZoneCount, allocatorID);
  auto zoneFacetOffsetsView = zoneFacetOffsets.view();
  axom::exclusive_scan<ExecSpace>(zoneFacetCountsView, zoneFacetOffsetsView);

  const auto conduitAllocatorID = axom::sidre::ConduitMemory::axomAllocIdToConduit(allocatorID);
  conduit::Node newConn;
  conduit::Node newSizes;
  conduit::Node newOffsets;
  conduit::Node newShapes;
  newConn.set_allocator(conduitAllocatorID);
  newSizes.set_allocator(conduitAllocatorID);
  newOffsets.set_allocator(conduitAllocatorID);
  newShapes.set_allocator(conduitAllocatorID);
  newConn.set(conduit::DataType(n_conn.dtype().id(), outputZoneCount * 3));
  newSizes.set(conduit::DataType(n_sizes.dtype().id(), outputZoneCount));
  newOffsets.set(conduit::DataType(n_offsets.dtype().id(), outputZoneCount));
  newShapes.set(conduit::DataType(n_sizes.dtype().id(), outputZoneCount));

  auto newConnView = bputils::make_array_view<ConnectivityType>(newConn);
  auto newSizesView = bputils::make_array_view<ConnectivityType>(newSizes);
  auto newOffsetsView = bputils::make_array_view<ConnectivityType>(newOffsets);
  auto newShapesView = bputils::make_array_view<ConnectivityType>(newShapes);

  axom::for_all<ExecSpace>(
    inputZoneCount,
    AXOM_LAMBDA(axom::IndexType z) {
      const axom::IndexType nFacets = zoneFacetCountsView[z];
      const axom::IndexType connStart = static_cast<axom::IndexType>(offsetsView[z]);
      const axom::IndexType triStart = zoneFacetOffsetsView[z];

      for(axom::IndexType f = 0; f < nFacets; ++f)
      {
        const axom::IndexType tri = triStart + f;
        const axom::IndexType outConn = tri * 3;
        newConnView[outConn + 0] = connView[connStart + 0];
        newConnView[outConn + 1] = connView[connStart + f + 1];
        newConnView[outConn + 2] = connView[connStart + f + 2];
        newSizesView[tri] = static_cast<ConnectivityType>(3);
        newOffsetsView[tri] = static_cast<ConnectivityType>(outConn);
        newShapesView[tri] = static_cast<ConnectivityType>(bpviews::Tri_ShapeID);
      }
    });

  if(n_output.has_child("fields"))
  {
    conduit::Node& n_fields = n_output["fields"];
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); ++i)
    {
      conduit::Node& n_field = n_fields[i];
      if(n_field.has_path("association") && n_field["association"].as_string() == "element" &&
         n_field.has_path("topology") && n_field["topology"].as_string() == topologyName &&
         n_field.has_child("values"))
      {
        duplicateElementValuesForTriangulation<ExecSpace>(n_field["values"],
                                                          outputZoneCount,
                                                          zoneFacetCountsView,
                                                          zoneFacetOffsetsView,
                                                          allocatorID);
      }
    }
  }

  n_conn.move(newConn);
  n_sizes.move(newSizes);
  n_offsets.move(newOffsets);
  conduit::Node& n_topo = n_output["topologies"][topologyName];
  conduit::Node& n_elems = n_topo.fetch_existing("elements");
  n_elems["shapes"].move(newShapes);
  n_elems["shape_map"].reset();
  n_elems["shape_map"][bpviews::TriTraits::name()] = bpviews::Tri_ShapeID;
}

/*!
 * @brief Convert a bump CutField Blueprint domain from polygonal surface
 * elements to a welded triangle mesh in place.
 *
 * This rewrites only topology connectivity and element-associated fields.  The
 * coordset is left untouched, and generated triangles reference the existing
 * welded vertex ids.
 */
template <int DIM, typename ExecSpace>
void triangulateBlueprintMesh(conduit::Node& n_output, int allocatorID)
{
  if constexpr(DIM != 3)
  {
    return;
  }

  namespace bputils = axom::bump::utilities;
  namespace bpviews = axom::bump::views;

  if(!n_output.has_child("topologies"))
  {
    return;  // empty contour (isovalue outside the data range): nothing to triangulate
  }
  const conduit::Node& n_topos = n_output.fetch_existing("topologies");
  SLIC_ASSERT(n_topos.number_of_children() == 1);
  const std::string topologyName = n_topos.child(0).name();
  conduit::Node& n_topo = n_output["topologies"][topologyName];
  conduit::Node& n_elems = n_topo.fetch_existing("elements");

  conduit::Node& n_conn = n_elems.fetch_existing("connectivity");
  conduit::Node& n_sizes = n_elems.fetch_existing("sizes");
  conduit::Node& n_offsets = n_elems.fetch_existing("offsets");

  SLIC_ERROR_IF(
    n_offsets.dtype().id() != n_sizes.dtype().id() || n_conn.dtype().id() != n_sizes.dtype().id(),
    "MarchingCubes bump Blueprint triangulation expects connectivity, "
    "sizes, and offsets to use the same integer type.");

  auto triangulateViews = [&](auto sizesView, auto offsetsView, auto connView) {
    triangulateBlueprintMeshViews<DIM, ExecSpace>(n_output,
                                                  n_conn,
                                                  n_sizes,
                                                  n_offsets,
                                                  topologyName,
                                                  sizesView,
                                                  offsetsView,
                                                  connView,
                                                  allocatorID);
  };

#if defined(_WIN32)
  triangulateViews(bputils::make_array_view<axom::IndexType>(n_sizes),
                   bputils::make_array_view<axom::IndexType>(n_offsets),
                   bputils::make_array_view<axom::IndexType>(n_conn));
#else
  bpviews::indexNodeToArrayViewSame(n_sizes, n_offsets, n_conn, std::move(triangulateViews));
#endif
}

template <int DIM, typename ExecSpace, typename SizesView, typename OffsetsView, typename ConnView, typename OrigView>
void adaptCutFieldOutputViews(const conduit::Node& n_coords,
                              SizesView sizesView,
                              OffsetsView offsetsView,
                              ConnView connView,
                              OrigView origView,
                              axom::ArrayView<axom::IndexType, 2> facetNodeIds,
                              axom::ArrayView<double, 2> facetNodeCoords,
                              axom::ArrayView<axom::IndexType, 1> facetParentIds,
                              axom::IndexType facetIndexOffset,
                              axom::IndexType nodeIndexOffset,
                              int objectAllocatorID)
{
  namespace bputils = axom::bump::utilities;

  const conduit::Node& n_x = n_coords.fetch_existing("values/x");
  const conduit::Node& n_y = n_coords.fetch_existing("values/y");
  auto xView = bputils::make_array_view<double>(n_x);
  auto yView = bputils::make_array_view<double>(n_y);
  // z only in 3D.
  axom::ArrayView<double> zView;
  if constexpr(DIM == 3)
  {
    const conduit::Node& n_z = n_coords.fetch_existing("values/z");
    zView = bputils::make_array_view<double>(n_z);
  }

  const axom::IndexType numZones = static_cast<axom::IndexType>(sizesView.size());
  const axom::IndexType numNodes = static_cast<axom::IndexType>(xView.size());

  axom::for_all<ExecSpace>(
    numNodes,
    AXOM_LAMBDA(axom::IndexType n) {
      facetNodeCoords(nodeIndexOffset + n, 0) = xView[n];
      facetNodeCoords(nodeIndexOffset + n, 1) = yView[n];
      if constexpr(DIM == 3)
      {
        facetNodeCoords(nodeIndexOffset + n, 2) = zView[n];
      }
    });

  // --- Per-zone facet offset (exclusive scan of facetsPerZone) -----------
  // We need, for each bump zone, the index of its first facet within this
  // domain so kernels can write without atomics.
  // Use the object's allocator, not the execution space default.
  const int allocatorID = objectAllocatorID;
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

  // --- The fan-triangulation kernel -------------------------------------
  // One thread per bump zone.  Each zone writes facetsPerZone facets;
  // each facet reuses bump's welded coordset vertex ids.
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
      const axom::IndexType parentId = static_cast<axom::IndexType>(origView[z]);

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

        for(int c = 0; c < DIM; ++c)
        {
          const axom::IndexType weldedNode =
            static_cast<axom::IndexType>(connView[connStart + local[c]]);
          facetNodeIds(facetIdx, c) = nodeIndexOffset + weldedNode;
        }

        facetParentIds[facetIdx] = parentId;
      }
    });
}

/*!
 * @brief Convert one bump CutField output (single domain) into the
 * fixed-corners-per-facet output buffers supplied by the parent MarchingCubes.
 *
 * @tparam DIM Spatial dimension (2 or 3).
 * @tparam ExecSpace Axom execution space.
 *
 * @param n_output The bump CutField output Blueprint mesh (in ExecSpace memory).
 * @param facetNodeIds [out] view, shape (totalFacetCount, DIM).
 * @param facetNodeCoords [out] view, shape (totalNodeCount, DIM).
 * @param facetParentIds [out] view, shape (totalFacetCount).
 * @param facetIndexOffset This domain's first facet index in the concatenated
 *   output (the parent's m_facetIndexOffsets[d]).
 * @param nodeIndexOffset This domain's first node index in the concatenated output.
 * @param thisDomainFacetCount Number of facets this domain produces (already
 *   computed by the caller; equals sum of facetsPerZone over the bump zones).
 * @param objectAllocatorID Allocator used for temporary arrays.
 *
 * @pre All output views and \a n_output live in ExecSpace's memory space.
 */
template <int DIM, typename ExecSpace>
void adaptCutFieldOutput(const conduit::Node& n_output,
                         axom::ArrayView<axom::IndexType, 2> facetNodeIds,
                         axom::ArrayView<double, 2> facetNodeCoords,
                         axom::ArrayView<axom::IndexType, 1> facetParentIds,
                         axom::IndexType facetIndexOffset,
                         axom::IndexType nodeIndexOffset,
                         axom::IndexType thisDomainFacetCount,
                         int objectAllocatorID)
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

  // bump always emits explicit sizes/offsets/connectivity for cut output.
  const conduit::Node& n_elems = n_topo.fetch_existing("elements");
  const conduit::Node& n_conn = n_elems.fetch_existing("connectivity");
  const conduit::Node& n_sizes = n_elems.fetch_existing("sizes");
  const conduit::Node& n_offsets = n_elems.fetch_existing("offsets");

  const std::string coordsetName = n_topo.fetch_existing("coordset").as_string();
  const conduit::Node& n_coords =
    n_output.fetch_existing(axom::fmt::format("coordsets/{}", coordsetName));

  // originalElements: element-associated, one entry per output zone (fragment).
  const conduit::Node& n_orig =
    n_output.fetch_existing(axom::fmt::format("fields/{}/values", kPublicOriginalElementsField));

  SLIC_ERROR_IF(n_offsets.dtype().id() != n_sizes.dtype().id() ||
                  n_conn.dtype().id() != n_sizes.dtype().id() ||
                  n_orig.dtype().id() != n_sizes.dtype().id(),
                "MarchingCubes bump adaptor expects connectivity, sizes, "
                "offsets, and originalElements to use the same integer type.");

  auto adaptViews = [&](auto sizesView, auto offsetsView, auto connView, auto origView) {
    adaptCutFieldOutputViews<DIM, ExecSpace>(n_coords,
                                             sizesView,
                                             offsetsView,
                                             connView,
                                             origView,
                                             facetNodeIds,
                                             facetNodeCoords,
                                             facetParentIds,
                                             facetIndexOffset,
                                             nodeIndexOffset,
                                             objectAllocatorID);
  };

#if defined(_WIN32)
  adaptViews(bputils::make_array_view<axom::IndexType>(n_sizes),
             bputils::make_array_view<axom::IndexType>(n_offsets),
             bputils::make_array_view<axom::IndexType>(n_conn),
             bputils::make_array_view<axom::IndexType>(n_orig));
#else
  bpviews::indexNodeToArrayViewSame(n_sizes, n_offsets, n_conn, n_orig, std::move(adaptViews));
#endif
}

}  // namespace axom::quest::detail::marching_cubes
