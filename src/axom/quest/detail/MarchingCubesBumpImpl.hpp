// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file MarchingCubesBumpImpl.hpp
 *
 * @brief Implements single-domain isocontouring with \c axom::bump::extraction::CutField.
 *
 * Bump accepts uniform, rectilinear, and explicit structured meshes,
 * plus single-shape quad meshes in 2D and hex meshes in 3D.
 * The selected execution space may be sequential, OpenMP, CUDA, or HIP.
 *
 * CutField performs extraction in one call, while ImplBase separates marking,
 * scanning, and facet generation. scanCrossings() therefore runs CutField and
 * caches its Blueprint output. computeFacets() copies that output into the
 * parent MarchingCubes buffers. The adaptor can fan-triangulate Bump's welded
 * polygonal faces when the caller requests legacy triangle output.
 */

#pragma once

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesBumpImpl.hpp requires conduit"
#endif
#ifndef AXOM_USE_BUMP
  #error "MarchingCubesBumpImpl.hpp requires bump"
#endif

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/reductions.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/quest/detail/MarchingCubesBumpAdaptor.hpp"

// Bump extraction and views.
#include "axom/bump/extraction/CutField.hpp"
#include "axom/bump/extraction/FieldIntersector.hpp"
#include "axom/bump/IndexingPolicies.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/views/dispatch_coordset.hpp"
#include "axom/bump/views/dispatch_topology.hpp"
#include "axom/bump/views/Shapes.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"

#include "conduit_node.hpp"
#include "conduit_blueprint.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace axom::quest::detail::marching_cubes
{
template <int DIM, typename ExecSpace, typename SizesView>
axom::IndexType computeTriangulatedFacetCountView(SizesView sizes)
{
  const axom::IndexType n = static_cast<axom::IndexType>(sizes.size());
  axom::ReduceSum<ExecSpace, axom::IndexType> facetCount(0);
  axom::for_all<ExecSpace>(n, [=] AXOM_HOST_DEVICE(axom::IndexType i) {
    const auto p = static_cast<axom::IndexType>(sizes[i]);
    facetCount += facetsPerZone<DIM>(p);
  });
  return facetCount.get();
}

/*!
 * @brief Bump-backed single-domain marching cubes implementation.
 *
 * @tparam DIM Spatial dimension (2 or 3).
 * @tparam ExecSpace Axom execution space (SEQ_EXEC, OMP_EXEC, CUDA_EXEC<>, HIP_EXEC<>).
 *
 * This object retains a reference to one Blueprint domain. scanCrossings()
 * runs CutField and caches its Blueprint output. computeFacets() then copies
 * that output into the buffers supplied by MarchingCubes.
 */
template <int DIM, typename ExecSpace>
class MarchingCubesBumpImpl : public MarchingCubesSingleDomain::ImplBase
{
public:
  static constexpr int SelectedDimensions = axom::bump::views::select_dimensions(DIM);
  static constexpr int ShapeTypes = (DIM == 3)
    ? axom::bump::views::select_shapes(axom::bump::views::Hex_ShapeID)
    : axom::bump::views::select_shapes(axom::bump::views::Quad_ShapeID);

  MarchingCubesBumpImpl(int allocatorID) : m_allocatorID(allocatorID) { }

  /*!
   * @brief Cache and validate the input domain.
   *
   * Extraction waits until scanCrossings(), after the field and isovalue are set.
   */
  void setDomain(const conduit::Node& dom,
                 const std::string& topologyName,
                 const std::string& maskFieldName) override
  {
    m_dom = &dom;
    m_topologyName = topologyName;
    m_maskFieldName = maskFieldName;

    if(!m_maskFieldName.empty())
    {
      const conduit::Node& n_mask =
        dom.fetch_existing(axom::fmt::format("fields/{}", m_maskFieldName));
      SLIC_ERROR_IF(n_mask.fetch_existing("association").as_string() != "element",
                    "MarchingCubes mask fields must be element-associated.");
      SLIC_ERROR_IF(!n_mask.has_path("values"),
                    "MarchingCubes mask field is missing a values node.");
      SLIC_ERROR_IF(!n_mask.fetch_existing("values").dtype().is_int32(),
                    "MarchingCubes mask field values must be int32.");
    }

    // Accept structured topologies and single-shape quads or hexes.
    const conduit::Node& n_topo =
      dom.fetch_existing(axom::fmt::format("topologies/{}", topologyName));
    const std::string topoType = n_topo.fetch_existing("type").as_string();

    // The output adaptor reads coordinates as double. Validate the input here
    // so an error names the offending path.
    const std::string csPath =
      axom::fmt::format("coordsets/{}/values", n_topo.fetch_existing("coordset").as_string());
    for(const char* comp : {"x", "y", "z"})
    {
      validateFieldIsFloat64(axom::fmt::format("{}/{}", csPath, comp), "coordset component");
    }
    if(!m_fcnFieldName.empty())
    {
      validateFieldIsFloat64(axom::fmt::format("fields/{}/values", m_fcnFieldName),
                             "function field");
    }
    if(topoType == "unstructured")
    {
      const std::string shape = n_topo.fetch_existing("elements/shape").as_string();
      const char* expected = (DIM == 3) ? "hex" : "quad";
      SLIC_ERROR_IF(shape != expected,
                    axom::fmt::format("MarchingCubes (Bump backend) supports unstructured "
                                      "single-shape '{}' in {}D, but got shape '{}'.",
                                      expected,
                                      DIM,
                                      shape));
    }
    else
    {
      SLIC_ERROR_IF(
        topoType != "uniform" && topoType != "rectilinear" && topoType != "structured",
        axom::fmt::format("MarchingCubes (Bump backend) does not support topology type '{}'.",
                          topoType));
    }
  }

  /*!
   * @brief Set the nodal function field, validating its type.
   *
   * MarchingCubes requires a float64 function field.
   */
  void setFunctionField(const std::string& fcnFieldName) override
  {
    m_fcnFieldName = fcnFieldName;
    if(m_dom != nullptr && !m_fcnFieldName.empty())
    {
      validateFieldIsFloat64(axom::fmt::format("fields/{}/values", m_fcnFieldName),
                             "function field");
    }
  }

  //! @brief Validate that an existing Blueprint array is float64.
  void validateFieldIsFloat64(const std::string& path, const std::string& what) const
  {
    if(m_dom == nullptr || !m_dom->has_path(path))
    {
      return;  // absence is reported elsewhere, with a better message
    }
    const conduit::Node& n = m_dom->fetch_existing(path);
    SLIC_ERROR_IF(!n.dtype().is_float64(),
                  axom::fmt::format("MarchingCubes (Bump backend) requires a float64 {} at '{}', "
                                    "but found '{}'.",
                                    what,
                                    path,
                                    n.dtype().name()));
  }

  void setContourValue(double contourVal) override { m_contourVal = contourVal; }

  void setMaskValue(int maskVal) override { m_maskVal = maskVal; }

  //! Retain the value for the shared interface; only the legacy backend reads it.
  void setDataParallelism(MarchingCubesDataParallelism dataPar) override
  {
    m_dataParallelism = dataPar;
  }

  //! @brief No-op for the Bump backend. scanCrossings() performs extraction.
  void markCrossings() override { }

  /*!
   * @brief Run the Bump extraction so the facet count is known.
   *
   * MarchingCubes needs each domain's facet count before it allocates the
   * shared output buffers. CutField cannot provide that count without running,
   * so this phase performs and caches the extraction.
   */
  void scanCrossings() override
  {
    clearDomain();
    runExtraction();
    m_extractionRan = true;
  }

  //! @brief Copy cached Bump output into the parent-allocated output buffers.
  void computeFacets() override { fillLegacyOutputBuffers(); }

  axom::IndexType getContourCellCount() const override { return m_facetCount; }

  axom::IndexType getContourNodeCount() const override
  {
    if(m_facetCount == 0)
    {
      return 0;
    }

    const conduit::Node& n_topos = m_output.fetch_existing("topologies");
    SLIC_ASSERT(n_topos.number_of_children() == 1);
    const conduit::Node& n_topo = n_topos.child(0);
    const std::string coordsetName = n_topo.fetch_existing("coordset").as_string();
    const conduit::Node& n_coords =
      m_output.fetch_existing(axom::fmt::format("coordsets/{}", coordsetName));
    return static_cast<axom::IndexType>(
      n_coords.fetch_existing("values/x").dtype().number_of_elements());
  }

  /*!
   * @brief Whether a Blueprint contour can be produced.
   *
   * True once computeIsocontour() has run, even if the contour is empty.
   */
  bool hasContourMeshBlueprint() const override { return m_extractionRan; }

  void copyContourMeshBlueprint(conduit::Node& bpMesh, bool triangulate) const override
  {
    SLIC_ERROR_IF(!m_extractionRan,
                  "MarchingCubes Bump backend has no Blueprint contour output. "
                  "Call computeIsocontour() before requesting it.");
    axom::bump::utilities::copy<ExecSpace>(bpMesh, m_output, m_allocatorID);
    if(triangulate)
    {
      triangulateBlueprintMesh<DIM, ExecSpace>(bpMesh, m_allocatorID);
    }
  }

  void relinquishContourMeshBlueprint(conduit::Node& bpMesh) override
  {
    SLIC_ERROR_IF(!m_extractionRan,
                  "MarchingCubes Bump backend has no Blueprint contour output. "
                  "Call computeIsocontour() before requesting it.");
    bpMesh.reset();
    bpMesh.swap(m_output);
    m_facetCount = 0;
    m_extractionRan = false;
  }

  void clearDomain() override
  {
    m_output.reset();
    m_facetCount = 0;
    m_extractionRan = false;
  }

#if !defined(__CUDACC__)
private:
#endif
  /*!
   * @brief Convert the requested isovalue to Bump's field type while matching
   * the legacy backend's greater-than-or-equal corner classification.
   */
  template <typename FieldType>
  static FieldType isoValueForBump(double contourVal)
  {
    const auto value = static_cast<FieldType>(contourVal);
    return std::nextafter(value, -std::numeric_limits<FieldType>::infinity());
  }

  template <typename SelectionPredicate>
  void buildSelectedZones(axom::IndexType nZones,
                          SelectionPredicate isSelected,
                          conduit::Node& n_options,
                          axom::Array<axom::IndexType>& selectedZones) const
  {
    axom::Array<axom::IndexType> maskFlags(nZones, nZones, m_allocatorID);
    auto maskFlagsView = maskFlags.view();

    axom::ReduceSum<ExecSpace, axom::IndexType> selectedCountReduce(0);
    axom::for_all<ExecSpace>(nZones, [=] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
      const axom::IndexType selected = isSelected(zoneIndex) ? 1 : 0;
      maskFlagsView[zoneIndex] = selected;
      selectedCountReduce += selected;
    });

    const axom::IndexType selectedCount = selectedCountReduce.get();
    selectedZones = axom::Array<axom::IndexType>(selectedCount, selectedCount, m_allocatorID);
    if(selectedCount == 0)
    {
      return;
    }

    axom::Array<axom::IndexType> selectedOffsets(nZones, nZones, m_allocatorID);
    auto selectedOffsetsView = selectedOffsets.view();
    axom::exclusive_scan<ExecSpace>(maskFlagsView, selectedOffsetsView);

    auto selectedZonesView = selectedZones.view();
    axom::for_all<ExecSpace>(nZones, [=] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
      if(maskFlagsView[zoneIndex] != 0)
      {
        selectedZonesView[selectedOffsetsView[zoneIndex]] = zoneIndex;
      }
    });

    n_options["selectedZones"].set_external(selectedZones.data(), selectedZones.size());
  }

  /*!
   * @brief Select cells that cross the isovalue and match the optional mask.
   *
   * Use FieldIntersector to match the extractor's corner classification.
   */
  template <typename TopologyView, typename CoordsetView>
  bool attachCrossingSelectedZonesOption(const TopologyView& topologyView,
                                         const CoordsetView& coordsetView,
                                         const conduit::Node& n_topo,
                                         const conduit::Node& n_coords,
                                         const conduit::Node& n_fields,
                                         conduit::Node& n_options,
                                         axom::Array<axom::IndexType>& crossingZones) const
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesBumpImpl::attachCrossingSelectedZonesOption");
    namespace bumpx = axom::bump::extraction;

    bumpx::FieldIntersector<ExecSpace, TopologyView, CoordsetView> intersector;
    intersector.setAllocatorID(m_allocatorID);
    intersector.initialize(topologyView, coordsetView, n_options, n_topo, n_coords, n_fields);
    const auto intersectorView = intersector.view();

    const TopologyView deviceTopologyView(topologyView);
    const auto isCrossing = [=] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
      const auto zone = deviceTopologyView.zone(zoneIndex);
      const auto ids = zone.getIds();
      const auto caseNumber = intersectorView.determineTableCase(zoneIndex, ids);
      const auto allPositive = (axom::IndexType {1} << ids.size()) - axom::IndexType {1};
      return caseNumber != 0 && caseNumber != allPositive;
    };

    const axom::IndexType nZones = topologyView.numberOfZones();
    if(m_maskFieldName.empty())
    {
      buildSelectedZones(nZones, isCrossing, n_options, crossingZones);
    }
    else
    {
      const conduit::Node& n_mask =
        m_dom->fetch_existing(axom::fmt::format("fields/{}", m_maskFieldName));
      const auto maskIndexing = axom::bump::makeElementFieldIndexing(topologyView, n_mask);
      const auto maskView =
        axom::bump::utilities::make_array_view<int>(n_mask.fetch_existing("values"));
      const int maskVal = m_maskVal;
      buildSelectedZones(
        nZones,
        [=] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
          return maskView[maskIndexing[zoneIndex]] == maskVal && isCrossing(zoneIndex);
        },
        n_options,
        crossingZones);
    }
    return !crossingZones.empty();
  }

  /*!
   * @brief Run CutField and store its Blueprint output.
   *
   * Bump dispatch converts the Blueprint topology and coordset to the view
   * types required by CutField. Input arrays must be accessible from \c ExecSpace.
   */
  void runExtraction()
  {
    SLIC_ASSERT(m_dom != nullptr);
    SLIC_ASSERT(!m_fcnFieldName.empty());

    namespace bumpviews = axom::bump::views;
    namespace bumpx = axom::bump::extraction;

    const conduit::Node& n_topo =
      m_dom->fetch_existing(axom::fmt::format("topologies/{}", m_topologyName));
    const std::string coordsetName = n_topo.fetch_existing("coordset").as_string();
    const conduit::Node& n_coords =
      m_dom->fetch_existing(axom::fmt::format("coordsets/{}", coordsetName));

    // Options shared by all dispatch branches.
    conduit::Node n_options;
    n_options["field"] = m_fcnFieldName;
    // Ask Bump to record the input zone that produced each output element.
    std::string originalElementsField = kOriginalElementsField;
    const conduit::Node& n_fields = m_dom->fetch_existing("fields");
    while(n_fields.has_child(originalElementsField))
    {
      originalElementsField += "_";
    }
    n_options["originalElementsField"] = originalElementsField;
    // Do not interpolate other input fields into the contour.
    n_options["fields"].set(conduit::DataType::object());

    conduit::Node& n_out = m_output;

    // Dispatch only this dimension and the supported unstructured shapes.
    // A valid view pair sets dispatched. The callback sets extracted only after
    // CutField runs, distinguishing an empty contour from an unsupported mesh.
    bool dispatched = false;
    bool extracted = false;
    bumpviews::dispatch_coordset<SelectedDimensions>(n_coords, [&](auto coordsetView) {
      using CoordsetView = decltype(coordsetView);
      bumpviews::dispatch_topology<SelectedDimensions, ShapeTypes>(
        n_topo,
        [&](const std::string& AXOM_UNUSED_PARAM(shape), auto topologyView) {
          using TopologyView = decltype(topologyView);
          dispatched = true;

          using Intersector = bumpx::FieldIntersector<ExecSpace, TopologyView, CoordsetView>;
          using Cut = bumpx::CutField<ExecSpace, TopologyView, CoordsetView, Intersector>;

          Intersector intersector;
          intersector.setAllocatorID(m_allocatorID);
          Cut iso(topologyView, coordsetView, intersector);
          iso.setAllocatorID(m_allocatorID);

          // Shift the threshold so Bump's strict comparison matches the legacy
          // kernel's greater-than-or-equal comparison. The prefilter and
          // extractor must use the same shifted value.
          using IsoFieldType =
            typename bumpx::FieldIntersector<ExecSpace, TopologyView, CoordsetView>::FieldType;
          const IsoFieldType isoForBump = isoValueForBump<IsoFieldType>(m_contourVal);
          n_options["value"] = static_cast<double>(isoForBump);

          axom::Array<axom::IndexType> selectedZones;
          const bool hasCrossingZones = attachCrossingSelectedZonesOption(topologyView,
                                                                          coordsetView,
                                                                          n_topo,
                                                                          n_coords,
                                                                          n_fields,
                                                                          n_options,
                                                                          selectedZones);
          if(!hasCrossingZones)
          {
            m_facetCount = 0;
            return;
          }

          conduit::Node execOptions;
          axom::bump::utilities::copy<ExecSpace>(execOptions, n_options, m_allocatorID);
          {
            AXOM_ANNOTATE_SCOPE("MarchingCubesBumpImpl::CutField::execute");
            iso.execute(*m_dom, execOptions, n_out);

            // Rename the field before exposing the result. The private request
            // name prevents an input field from replacing the parent-zone ids.
            const std::string privateField = axom::fmt::format("fields/{}", originalElementsField);
            if(n_out.has_path(privateField))
            {
              n_out["fields"].rename_child(originalElementsField, kPublicOriginalElementsField);
            }
          }
          extracted = true;
        });
    });

    SLIC_ERROR_IF(
      !dispatched,
      axom::fmt::format("MarchingCubes (Bump backend) could not build views for topology '{}' "
                        "(type '{}') with coordset '{}' (type '{}') in {}D.",
                        m_topologyName,
                        n_topo.fetch_existing("type").as_string(),
                        coordsetName,
                        n_coords.fetch_existing("type").as_string(),
                        DIM));

    if(!extracted)
    {
      // No selected zone crosses the isovalue, so the contour is empty.
      m_output.reset();
      m_facetCount = 0;
      return;
    }

    // Count the facets in the fixed-stride output after polygon triangulation.
    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesBumpImpl::computeTriangulatedFacetCount");
      m_facetCount = computeTriangulatedFacetCount(n_out);
    }
  }

  /*!
   * @brief Count fixed-stride facets after fan-triangulating Bump output.
   *
   * In 2D, each two-node segment is one facet.
   * In 3D, a polygon with \c p corners produces \c p-2 triangles.
   */
  axom::IndexType computeTriangulatedFacetCount(const conduit::Node& n_out) const
  {
    const std::string newTopoName = onlyTopologyName(n_out);
    const conduit::Node& n_elems =
      n_out.fetch_existing(axom::fmt::format("topologies/{}/elements", newTopoName));

    // Polygonal/segment output carries an explicit "sizes" array.
    if(n_elems.has_child("sizes"))
    {
      const conduit::Node& n_sizes = n_elems.fetch_existing("sizes");
      axom::IndexType facets = 0;
      axom::bump::views::nodeToArrayView(n_sizes, [&](auto sizes) {
        facets = computeTriangulatedFacetCountView<DIM, ExecSpace>(sizes);
      });
      return facets;
    }

    SLIC_ERROR(axom::fmt::format(
      "MarchingCubes Bump backend: cut output topology '{}' has no 'sizes' array. "
      "The adaptor requires explicit sizes on Bump's cut output.",
      newTopoName));
    return 0;
  }

  /*!
   * @brief Fill the parent-allocated fixed-stride buffers from cached Bump output.
   *
   * Polygonal faces are fan-triangulated and reuse Bump's welded vertices.
   */
  void fillLegacyOutputBuffers()
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesBumpImpl::triangulateAndAdapt");
    if(m_facetCount == 0)
    {
      return;
    }

    adaptCutFieldOutput<DIM, ExecSpace>(m_output,
                                        m_facetNodeIds,
                                        m_facetNodeCoords,
                                        m_facetParentIds,
                                        m_facetIndexOffset,
                                        m_nodeIndexOffset,
                                        m_facetCount,
                                        m_allocatorID);
  }

  //! @brief Return the single topology name in a Bump output node.
  static std::string onlyTopologyName(const conduit::Node& n_out)
  {
    const conduit::Node& n_topos = n_out.fetch_existing("topologies");
    SLIC_ASSERT(n_topos.number_of_children() == 1);
    return n_topos.child(0).name();
  }

private:
  int m_allocatorID = axom::INVALID_ALLOCATOR_ID;

  const conduit::Node* m_dom {nullptr};
  std::string m_topologyName;
  std::string m_fcnFieldName;
  std::string m_maskFieldName;

  //! @brief Cached Bump CutField output (Blueprint mesh).
  conduit::Node m_output;

  //! @brief Fixed-stride facet count after fan triangulation.
  axom::IndexType m_facetCount {};

  //! @brief Whether extraction ran, distinguishing unavailable from empty output.
  bool m_extractionRan {false};
};

}  // namespace axom::quest::detail::marching_cubes
