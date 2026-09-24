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
#include "axom/core/MDMapping.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/quest/MeshViewUtil.hpp"
#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/quest/detail/MarchingCubesBumpAdaptor.hpp"

// Bump extraction and views.
#include "axom/bump/extraction/CutField.hpp"
#include "axom/bump/extraction/FieldIntersector.hpp"
#include "axom/bump/SelectedZones.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/views/dispatch_coordset.hpp"
#include "axom/bump/views/dispatch_topology.hpp"
#include "axom/bump/views/dispatch_unstructured_topology.hpp"
#include "axom/bump/views/Shapes.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"

#include "conduit_node.hpp"
#include "conduit_blueprint.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
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
    facetCount += (DIM == 3) ? (p >= 3 ? p - 2 : 0) : (p >= 2 ? 1 : 0);
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
  static constexpr auto MemorySpace = execution_space<ExecSpace>::memory_space;
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

    // MeshViewUtil supports only structured topologies with explicit
    // coordsets. Use it for the crossing prefilter only in that case.

    const std::string coordsetTypeForPath =
      dom
        .fetch_existing(
          axom::fmt::format("coordsets/{}", n_topo.fetch_existing("coordset").as_string()))
        .fetch_existing("type")
        .as_string();
    m_useMeshViewUtilPath = (topoType == "structured") && (coordsetTypeForPath == "explicit");

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

  // Retain the value for the shared interface; only the legacy backend reads it.
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
    m_extractionRan = true;
    runExtraction();
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

    SLIC_ASSERT(m_output != nullptr);
    const conduit::Node& n_topos = m_output->fetch_existing("topologies");
    SLIC_ASSERT(n_topos.number_of_children() == 1);
    const conduit::Node& n_topo = n_topos.child(0);
    const std::string coordsetName = n_topo.fetch_existing("coordset").as_string();
    const conduit::Node& n_coords =
      m_output->fetch_existing(axom::fmt::format("coordsets/{}", coordsetName));
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
    if(m_output == nullptr)
    {
      bpMesh.reset();
      return;
    }
    axom::bump::utilities::copy<ExecSpace>(bpMesh, *m_output, m_allocatorID);
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
    if(m_output != nullptr)
    {
      bpMesh.swap(*m_output);
      m_output.reset();
    }
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

  /*! @brief Dispatch a coordset view restricted to this implementation's DIM. */
  template <typename FuncType>
  static void dispatchCoordset(const conduit::Node& n_coords, FuncType&& func)
  {
    axom::bump::views::dispatch_coordset<SelectedDimensions>(n_coords, std::forward<FuncType>(func));
  }

  /*!
   * @brief Dispatch a topology in the template's spatial dimension.
   *
   * Unstructured meshes are limited to quads in 2D and hexes in 3D.
   */
  template <typename FuncType>
  static void dispatchTopology(const conduit::Node& n_topo, FuncType&& func)
  {
    axom::bump::views::dispatch_topology<SelectedDimensions, ShapeTypes>(
      n_topo,
      std::forward<FuncType>(func));
  }

  void attachSelectedZonesOption(conduit::Node& n_options,
                                 axom::Array<axom::IndexType>& selectedZones) const
  {
    conduit::Node& n_selectedZones = n_options["selectedZones"];
    if(selectedZones.empty())
    {
      n_selectedZones.set(
        conduit::DataType(axom::bump::utilities::cpp2conduit<axom::IndexType>::id, 0));
    }
    else
    {
      n_selectedZones.set_external(selectedZones.data(), selectedZones.size());
    }
  }

  template <typename MaskPredicate>
  void buildSelectedZonesFromMask(axom::IndexType nZones,
                                  MaskPredicate isSelected,
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

    attachSelectedZonesOption(n_options, selectedZones);
  }

  /*!
   * @brief Restrict extraction to zones with the requested mask value.
   *
   * The topology view maps compact zone indices to mask field indices.
   * This accounts for padding in strided structured fields.
   */
  template <typename TopologyView>
  void addMaskSelectedZonesOption(const TopologyView& topologyView,
                                  conduit::Node& n_options,
                                  axom::Array<axom::IndexType>& selectedZones) const
  {
    namespace bputils = axom::bump::utilities;

    if(m_maskFieldName.empty())
    {
      return;
    }

    const axom::IndexType nZones = topologyView.numberOfZones();
    const conduit::Node& n_mask =
      m_dom->fetch_existing(axom::fmt::format("fields/{}", m_maskFieldName));
    const conduit::Node& n_maskValues = n_mask.fetch_existing("values");

    // Copy the member so the device predicate does not dereference the host `this` pointer.
    const int maskVal = m_maskVal;

    auto maskView = bputils::make_array_view<int>(n_maskValues);
    const TopologyView deviceTopologyView(topologyView);
    buildSelectedZonesFromMask(
      nZones,
      [maskView, maskVal, deviceTopologyView] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
        return maskView[deviceTopologyView.zoneFieldIndex(zoneIndex)] == maskVal;
      },
      n_options,
      selectedZones);
  }

  /*!
   * @brief Filter the current zone selection to cells that cross the isovalue.
   *
   * This generic path uses FieldIntersector when MeshViewUtil cannot supply structured indexing.
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

    axom::bump::SelectedZones<ExecSpace> selectedZones(topologyView.numberOfZones(),
                                                       n_options,
                                                       "selectedZones",
                                                       m_allocatorID);
    const auto selectedZonesView = selectedZones.view();
    if(selectedZonesView.empty())
    {
      return false;
    }

    bumpx::FieldIntersector<ExecSpace, TopologyView, CoordsetView> intersector;
    intersector.setAllocatorID(m_allocatorID);
    intersector.initialize(topologyView, coordsetView, n_options, n_topo, n_coords, n_fields);
    const auto intersectorView = intersector.view();

    axom::ReduceSum<ExecSpace, axom::IndexType> crossingCount(0);
    axom::Array<axom::IndexType> crossingFlags(selectedZonesView.size(),
                                               selectedZonesView.size(),
                                               m_allocatorID);
    auto crossingFlagsView = crossingFlags.view();
    const TopologyView deviceTopologyView(topologyView);
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      [=] AXOM_HOST_DEVICE(axom::IndexType selectedIndex) {
        const auto zoneIndex = selectedZonesView[selectedIndex];
        const auto zone = deviceTopologyView.zone(zoneIndex);
        const auto ids = zone.getIds();
        const auto caseNumber = intersectorView.determineTableCase(zoneIndex, ids);
        const auto allPositive = (axom::IndexType {1} << ids.size()) - axom::IndexType {1};
        const axom::IndexType crosses = (caseNumber != 0 && caseNumber != allPositive) ? 1 : 0;
        crossingFlagsView[selectedIndex] = crosses;
        crossingCount += crosses;
      });

    const axom::IndexType crossingCountValue = crossingCount.get();
    crossingZones =
      axom::Array<axom::IndexType>(crossingCountValue, crossingCountValue, m_allocatorID);

    axom::Array<axom::IndexType> crossingOffsets(selectedZonesView.size(),
                                                 selectedZonesView.size(),
                                                 m_allocatorID);
    auto crossingOffsetsView = crossingOffsets.view();
    axom::exclusive_scan<ExecSpace>(crossingFlagsView, crossingOffsetsView);

    auto crossingZonesView = crossingZones.view();
    axom::for_all<ExecSpace>(selectedZonesView.size(),
                             [=] AXOM_HOST_DEVICE(axom::IndexType selectedIndex) {
                               if(crossingFlagsView[selectedIndex] != 0)
                               {
                                 crossingZonesView[crossingOffsetsView[selectedIndex]] =
                                   selectedZonesView[selectedIndex];
                               }
                             });

    attachSelectedZonesOption(n_options, crossingZones);
    return crossingCountValue > 0;
  }

  /*!
   * @brief Fast structured crossing pre-filter.
   *
   * @param isoForBump Threshold in the intersector's field type; see isoValueForBump().
   *   This prefilter must use the same corner test as Bump's FieldIntersector.
   *   Otherwise it can discard a zone that FieldIntersector would cut.
   */
  template <typename IsoFieldType>
  bool attachStructuredCrossingSelectedZonesOption(IsoFieldType isoForBump,
                                                   conduit::Node& n_options,
                                                   axom::Array<axom::IndexType>& crossingZones) const
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesBumpImpl::attachStructuredCrossingSelectedZonesOption");

    axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(*m_dom, m_topologyName);
    const auto fcnView = mvu.template getConstFieldView<double>(m_fcnFieldName, false);
    axom::ArrayView<const int, DIM, MemorySpace> maskView;
    if(!m_maskFieldName.empty())
    {
      maskView = mvu.template getConstFieldView<int>(m_maskFieldName, false);
    }

    const auto cellShape = mvu.getCellShape();
    const axom::MDMapping<DIM> topoMap(cellShape, axom::ArrayStrideOrder::COLUMN);
    const axom::IndexType nZones = mvu.getCellCount();

    if constexpr(std::is_same_v<ExecSpace, axom::SEQ_EXEC>)
    {
      /*
        Iterate over logical indices. Calling topoMap.toMultiIndex() for every zone
        costs DIM integer divisions, while nested loops update the flat index by addition.

        Cache each node's sign in a byte. Adjacent zones share corners, so this
        avoids repeated strided reads of the double field.
      */
      crossingZones = axom::Array<axom::IndexType>(0, 0, m_allocatorID);
      crossingZones.reserve(nZones);
      {
        // Node-sign plane cache: signs for logical k and k+1 (3D), or the single plane (2D).
        // Indexed [j * pi + i] over NODE counts.
        const axom::IndexType pi = cellShape[0] + 1;
        const axom::IndexType pj = cellShape[1] + 1;
        const axom::IndexType planeSize = pi * pj;
        axom::Array<std::uint8_t> signPlanes(2 * planeSize, 2 * planeSize);
        auto signs = signPlanes.view();

        auto fillPlane = [&](axom::IndexType which, axom::IndexType k) {
          std::uint8_t* dst = signs.data() + which * planeSize;
          for(axom::IndexType j = 0; j < pj; ++j)
          {
            for(axom::IndexType i = 0; i < pi; ++i)
            {
              if constexpr(DIM == 2)
              {
                AXOM_UNUSED_VAR(k);
                dst[j * pi + i] = static_cast<IsoFieldType>(fcnView(i, j)) > isoForBump ? 1 : 0;
              }
              else
              {
                dst[j * pi + i] = static_cast<IsoFieldType>(fcnView(i, j, k)) > isoForBump ? 1 : 0;
              }
            }
          }
        };

        const axom::IndexType nk = (DIM == 3) ? cellShape[DIM - 1] : 1;
        fillPlane(0, 0);

        for(axom::IndexType k = 0; k < nk; ++k)
        {
          if constexpr(DIM == 3)
          {
            // Plane k is already in slot (k % 2); fill k+1 into the other slot.
            fillPlane((k + 1) % 2, k + 1);
          }
          const std::uint8_t* lo = signs.data() + (DIM == 3 ? (k % 2) : 0) * planeSize;
          const std::uint8_t* hi = signs.data() + (DIM == 3 ? ((k + 1) % 2) : 0) * planeSize;

          for(axom::IndexType j = 0; j < cellShape[1]; ++j)
          {
            const axom::IndexType row = j * pi;
            const axom::IndexType rowUp = (j + 1) * pi;
            for(axom::IndexType i = 0; i < cellShape[0]; ++i)
            {
              bool useZone = maskView.empty();
              if(!useZone)
              {
                if constexpr(DIM == 2)
                {
                  useZone = (maskView(i, j) == m_maskVal);
                }
                else
                {
                  useZone = (maskView(i, j, k) == m_maskVal);
                }
              }
              if(!useZone)
              {
                continue;
              }

              int nPos = lo[row + i] + lo[row + i + 1] + lo[rowUp + i] + lo[rowUp + i + 1];
              int nCorners = 4;
              if constexpr(DIM == 3)
              {
                nPos += hi[row + i] + hi[row + i + 1] + hi[rowUp + i] + hi[rowUp + i + 1];
                nCorners = 8;
              }

              if(nPos != 0 && nPos != nCorners)
              {
                if constexpr(DIM == 2)
                {
                  crossingZones.push_back(i + j * cellShape[0]);
                }
                else
                {
                  crossingZones.push_back(i + cellShape[0] * (j + cellShape[1] * k));
                }
              }
            }
          }
        }
      }

      attachSelectedZonesOption(n_options, crossingZones);
      return !crossingZones.empty();
    }

    axom::IndexType crossingCountValue = 0;
    {
      AXOM_ANNOTATE_BEGIN("MarchingCubesBumpImpl::crossingFlagAllocation");
      axom::Array<axom::IndexType> crossingFlags(nZones, nZones, m_allocatorID);
      AXOM_ANNOTATE_END("MarchingCubesBumpImpl::crossingFlagAllocation");
      auto crossingFlagsView = crossingFlags.view();

      // Copy member values to locals so the device lambda captures values rather than `this`.
      const IsoFieldType isoVal = isoForBump;
      const int maskVal = m_maskVal;
      axom::ReduceSum<ExecSpace, axom::IndexType> crossingCount(0);
      AXOM_ANNOTATE_BEGIN("MarchingCubesBumpImpl::crossingClassification");
      axom::for_all<ExecSpace>(
        nZones,
        [topoMap, maskView, fcnView, isoVal, maskVal, crossingFlagsView, crossingCount] AXOM_HOST_DEVICE(
          axom::IndexType zoneIndex) {
          const auto idx = topoMap.toMultiIndex(zoneIndex);
          bool useZone = maskView.empty();
          if(!useZone)
          {
            if constexpr(DIM == 2)
            {
              useZone = (maskView(idx[0], idx[1]) == maskVal);
            }
            else
            {
              useZone = (maskView(idx[0], idx[1], idx[2]) == maskVal);
            }
          }

          bool hasPositive = false;
          bool hasNonPositive = false;
          if(useZone)
          {
            if constexpr(DIM == 2)
            {
              const bool p0 = static_cast<IsoFieldType>(fcnView(idx[0], idx[1])) > isoVal;
              const bool p1 = static_cast<IsoFieldType>(fcnView(idx[0] + 1, idx[1])) > isoVal;
              const bool p2 = static_cast<IsoFieldType>(fcnView(idx[0] + 1, idx[1] + 1)) > isoVal;
              const bool p3 = static_cast<IsoFieldType>(fcnView(idx[0], idx[1] + 1)) > isoVal;
              hasPositive = p0 || p1 || p2 || p3;
              hasNonPositive = !p0 || !p1 || !p2 || !p3;
            }
            else
            {
              const bool p0 = static_cast<IsoFieldType>(fcnView(idx[0], idx[1], idx[2])) > isoVal;
              const bool p1 = static_cast<IsoFieldType>(fcnView(idx[0] + 1, idx[1], idx[2])) > isoVal;
              const bool p2 = static_cast<IsoFieldType>(fcnView(idx[0], idx[1] + 1, idx[2])) > isoVal;
              const bool p3 =
                static_cast<IsoFieldType>(fcnView(idx[0] + 1, idx[1] + 1, idx[2])) > isoVal;
              const bool p4 = static_cast<IsoFieldType>(fcnView(idx[0], idx[1], idx[2] + 1)) > isoVal;
              const bool p5 =
                static_cast<IsoFieldType>(fcnView(idx[0] + 1, idx[1], idx[2] + 1)) > isoVal;
              const bool p6 =
                static_cast<IsoFieldType>(fcnView(idx[0], idx[1] + 1, idx[2] + 1)) > isoVal;
              const bool p7 =
                static_cast<IsoFieldType>(fcnView(idx[0] + 1, idx[1] + 1, idx[2] + 1)) > isoVal;
              hasPositive = p0 || p1 || p2 || p3 || p4 || p5 || p6 || p7;
              hasNonPositive = !p0 || !p1 || !p2 || !p3 || !p4 || !p5 || !p6 || !p7;
            }
          }

          const axom::IndexType crosses = (hasPositive && hasNonPositive) ? 1 : 0;
          crossingFlagsView[zoneIndex] = crosses;
          crossingCount += crosses;
        });

      crossingCountValue = crossingCount.get();
      AXOM_ANNOTATE_END("MarchingCubesBumpImpl::crossingClassification");

      AXOM_ANNOTATE_BEGIN("MarchingCubesBumpImpl::crossingZoneAllocation");
      crossingZones =
        axom::Array<axom::IndexType>(crossingCountValue, crossingCountValue, m_allocatorID);
      AXOM_ANNOTATE_END("MarchingCubesBumpImpl::crossingZoneAllocation");

      AXOM_ANNOTATE_BEGIN("MarchingCubesBumpImpl::crossingOffsetAllocation");
      axom::Array<axom::IndexType> crossingOffsets(nZones, nZones, m_allocatorID);
      AXOM_ANNOTATE_END("MarchingCubesBumpImpl::crossingOffsetAllocation");
      auto crossingOffsetsView = crossingOffsets.view();

      AXOM_ANNOTATE_BEGIN("MarchingCubesBumpImpl::crossingScan");
      axom::exclusive_scan<ExecSpace>(crossingFlagsView, crossingOffsetsView);
      AXOM_ANNOTATE_END("MarchingCubesBumpImpl::crossingScan");

      auto crossingZonesView = crossingZones.view();
      AXOM_ANNOTATE_BEGIN("MarchingCubesBumpImpl::crossingCompaction");
      axom::for_all<ExecSpace>(nZones, [=] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
        if(crossingFlagsView[zoneIndex] != 0)
        {
          crossingZonesView[crossingOffsetsView[zoneIndex]] = zoneIndex;
        }
      });
      AXOM_ANNOTATE_END("MarchingCubesBumpImpl::crossingCompaction");

      AXOM_ANNOTATE_BEGIN("MarchingCubesBumpImpl::crossingScratchRelease");
    }
    AXOM_ANNOTATE_END("MarchingCubesBumpImpl::crossingScratchRelease");

    attachSelectedZonesOption(n_options, crossingZones);
    return crossingCountValue > 0;
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
    n_options["value"] = m_contourVal;
    // Ask Bump to record the input zone that produced each output element.
    n_options["originalElementsField"] = kOriginalElementsField;
    // Do not interpolate other input fields into the contour.
    n_options["fields"].set(conduit::DataType::object());

    m_output = std::make_unique<conduit::Node>();
    conduit::Node& n_out = *m_output;

    // Dispatch only this dimension and the supported unstructured shapes.
    // A valid view pair sets dispatched. The callback sets extracted only after
    // CutField runs, distinguishing an empty contour from an unsupported mesh.
    bool dispatched = false;
    bool extracted = false;
    dispatchCoordset(n_coords, [&](auto coordsetView) {
      using CoordsetView = decltype(coordsetView);
      dispatchTopology(n_topo, [&](const std::string& AXOM_UNUSED_PARAM(shape), auto topologyView) {
        using TopologyView = decltype(topologyView);
        dispatched = true;

        using Cut = bumpx::CutField<ExecSpace, TopologyView, CoordsetView>;

        Cut iso(topologyView, coordsetView);
        iso.setAllocatorID(m_allocatorID);

        // Shift the threshold so Bump's strict comparison matches the legacy
        // kernel's greater-than-or-equal comparison. The prefilter and
        // extractor must use the same shifted value.
        using IsoFieldType =
          typename bumpx::FieldIntersector<ExecSpace, TopologyView, CoordsetView>::FieldType;
        const IsoFieldType isoForBump = isoValueForBump<IsoFieldType>(m_contourVal);
        n_options["value"] = static_cast<double>(isoForBump);

        axom::Array<axom::IndexType> selectedZones;
        const bool hasCrossingZones = m_useMeshViewUtilPath
          ? attachStructuredCrossingSelectedZonesOption<IsoFieldType>(isoForBump,
                                                                      n_options,
                                                                      selectedZones)
          : [&]() {
              addMaskSelectedZonesOption(topologyView, n_options, selectedZones);
              return attachCrossingSelectedZonesOption(topologyView,
                                                       coordsetView,
                                                       n_topo,
                                                       n_coords,
                                                       m_dom->fetch_existing("fields"),
                                                       n_options,
                                                       selectedZones);
            }();
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
          const std::string privateField = axom::fmt::format("fields/{}", kOriginalElementsField);
          if(n_out.has_path(privateField))
          {
            n_out["fields"].rename_child(kOriginalElementsField, kPublicOriginalElementsField);
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
    SLIC_ASSERT(m_output != nullptr);

    adaptCutFieldOutput<DIM, ExecSpace>(*m_output,
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

  //! @brief Whether the structured-explicit crossing prefilter is available.
  bool m_useMeshViewUtilPath {false};

  //! @brief Cached Bump CutField output (Blueprint mesh).
  std::unique_ptr<conduit::Node> m_output;

  //! @brief Fixed-stride facet count after fan triangulation.
  axom::IndexType m_facetCount {};

  //! @brief Whether extraction ran, distinguishing unavailable from empty output.
  bool m_extractionRan {false};
};

}  // namespace axom::quest::detail::marching_cubes
