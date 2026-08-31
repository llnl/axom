// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file MarchingCubesBumpImpl.hpp
 *
 * @brief A MarchingCubesSingleDomain::ImplBase implementation
 * that delegates isocontour extraction to axom::bump::extraction::CutField.
 *
 * Unlike the legacy MarchingCubesImpl (which contains a hand-written marching cubes kernel over structured data),
 * this implementation wraps the bump CutField extractor.
 * CutField is templated on <ExecSpace, TopologyView, CoordsetView>, so it transparently supports:
 *   - structured (uniform / rectilinear / explicit-structured) topologies, and
 *   - unstructured *single-shape* quad (2D) and hex (3D) topologies,
 * on all execution spaces (seq, omp, cuda, hip), via the bump view dispatch.
 *
 * Design notes:
 *  - CutField is a Blueprint-in / Blueprint-out operation.
 *    We run it once per domain in computeFacets();
 *    markCrossings()/scanCrossings() do the cheap bookkeeping the parent MarchingCubes orchestration expects.
 *  - The phased ImplBase interface (mark/scan/compute) was designed around the legacy kernel.
 *    bump does everything in one execute() call, so we run the extractor lazily and cache its result,
 *    then satisfy the count queries from the cached result.
 *  - bump produces a welded, topologically-connected surface (blend-group uniquification).
 *    The 3D output may be polygonal (tri/quad/poly5..8), and the 2D output is line segments.
 *    The adaptor can optionally triangulate the polygon.
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

// bump extraction + views
#include "axom/bump/extraction/CutField.hpp"
#include "axom/bump/extraction/FieldIntersector.hpp"
#include "axom/bump/SelectedZones.hpp"
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
  axom::for_all<ExecSpace>(
    n,
    AXOM_LAMBDA(axom::IndexType i) {
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
 * This object holds a reference to a single Blueprint domain and, on
 * computeFacets(), invokes bump's CutField extractor to produce the isocontour.
 * It then adapts the bump Blueprint output into the legacy output buffers
 * supplied by the parent MarchingCubes via setOutputBuffers().
 */
template <int DIM, typename ExecSpace>
class MarchingCubesBumpImpl : public MarchingCubesSingleDomain::ImplBase
{
public:
  static constexpr auto MemorySpace = execution_space<ExecSpace>::memory_space;
  static constexpr int SelectedDimensions = axom::bump::views::select_dimensions(DIM);
  static constexpr int ShapeTypes =
    (DIM == 3) ? (1 << axom::bump::views::Hex_ShapeID) : (1 << axom::bump::views::Quad_ShapeID);

  MarchingCubesBumpImpl(int allocatorID) : m_allocatorID(allocatorID) { }

  /*!
   * @brief Cache the domain and topology/mask names.
   *
   * We do not build views here because CutField needs the field name and
   * contour value too; we defer all heavy work to computeFacets().
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

    // Validate that this is a topology bump+MarchingCubes supports: a DIM-dimensional structured topology,
    // or an unstructured single-shape quad (DIM==2) / hex (DIM==3) topology.  Mixed/polyhedral are rejected here.
    const conduit::Node& n_topo =
      dom.fetch_existing(axom::fmt::format("topologies/{}", topologyName));
    const std::string topoType = n_topo.fetch_existing("type").as_string();

    // MeshViewUtil::isValid() requires both a "structured" topology and an explicit coordset
    // Gating those paths on m_isStructured meant a uniform/rectilinear mesh passed setDomain()'s validation,
    // and then hard-errored deep inside the crossing pre-filter.
    // bump handles uniform and rectilinear, so use bump's views rather than MeshViewUtil.

    const std::string coordsetTypeForPath =
      dom
        .fetch_existing(
          axom::fmt::format("coordsets/{}", n_topo.fetch_existing("coordset").as_string()))
        .fetch_existing("type")
        .as_string();
    m_useMeshViewUtilPath = (topoType == "structured") && (coordsetTypeForPath == "explicit");

    // Strided-structured (ghost-padded) input needs no special handling here.
    // adaptCutFieldOutputViews() reads bump's blended output coordset with make_array_view<double>,
    // which errors late and opaquely on a float32 coordset.
    // Catch it here instead, where the path can be named.
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
      validateFieldStrideOrder(axom::fmt::format("fields/{}", m_fcnFieldName));
    }
    if(topoType == "unstructured")
    {
      const std::string shape = n_topo.fetch_existing("elements/shape").as_string();
      const char* expected = (DIM == 3) ? "hex" : "quad";
      SLIC_ERROR_IF(shape != expected,
                    axom::fmt::format("MarchingCubes (bump backend) supports unstructured "
                                      "single-shape '{}' in {}D, but got shape '{}'.",
                                      expected,
                                      DIM,
                                      shape));
    }
    else
    {
      SLIC_ERROR_IF(
        topoType != "uniform" && topoType != "rectilinear" && topoType != "structured",
        axom::fmt::format("MarchingCubes (bump backend) does not support topology type '{}'.",
                          topoType));
    }
  }

  /*!
   * @brief Set the nodal function field, validating its type.
   *
   * The structured pre-filter reads the field through MeshViewUtil::getConstFieldView<double>(),
   * which assumes that the values are `double`. Check the type here, and provide an error
   * message with the actual type when necessary.
   */
  void setFunctionField(const std::string& fcnFieldName) override
  {
    m_fcnFieldName = fcnFieldName;
    if(m_dom != nullptr && !m_fcnFieldName.empty())
    {
      validateFieldIsFloat64(axom::fmt::format("fields/{}/values", m_fcnFieldName),
                             "function field");
      validateFieldStrideOrder(axom::fmt::format("fields/{}", m_fcnFieldName));
    }
  }

  /*!
   * @brief Validate the function field layout used by bump's flat field view.
   *
   * bump's FieldIntersector reads the field with a flat make_array_view,
   * ignoring the field's Blueprint offsets/strides. That is only correct for a compact i-fastest layout,
   * or a ghost-padded i-fastest layout whose offsets and strides match the structured topology.
   * Reject independently strided fields and topology arrays cases.
   *
   * i-fastest is characterised by strides[0] == 1 and non-decreasing strides.
   *
   * @note Uniform, rectilinear, and unstructured topologies use compact node
   *   numbering, so field offsets/strides are not supported on those paths.
   */
  void validateFieldStrideOrder(const std::string& fieldPath) const
  {
    if(m_dom == nullptr || !m_dom->has_path(fieldPath))
    {
      return;
    }
    const conduit::Node& n_field = m_dom->fetch_existing(fieldPath);
    const bool hasFieldOffsets = n_field.has_child("offsets");
    const bool hasFieldStrides = n_field.has_child("strides");

    const conduit::Node& n_topo =
      m_dom->fetch_existing(axom::fmt::format("topologies/{}", m_topologyName));
    const bool hasTopoOffsets = n_topo.has_path("elements/dims/offsets");
    const bool hasTopoStrides = n_topo.has_path("elements/dims/strides");

    if(!hasFieldOffsets && !hasFieldStrides && !hasTopoOffsets && !hasTopoStrides)
    {
      return;
    }
    if(!m_useMeshViewUtilPath)
    {
      SLIC_ERROR(axom::fmt::format(
        "MarchingCubes (bump backend) does not support function-field offsets/strides "
        "on topology '{}'; bump indexes this field with compact topology node ids.",
        m_topologyName));
      return;
    }

    axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(*m_dom, m_topologyName);
    const auto nodeShape = mvu.getNodeShape();
    axom::StackArray<axom::IndexType, DIM> fieldOffsets {};
    axom::StackArray<axom::IndexType, DIM> topoOffsets {};
    axom::StackArray<axom::IndexType, DIM> fieldStrides {};
    axom::StackArray<axom::IndexType, DIM> topoStrides {};

    axom::IndexType compactStride = 1;
    for(int d = 0; d < DIM; ++d)
    {
      fieldStrides[d] = compactStride;
      topoStrides[d] = compactStride;
      compactStride *= nodeShape[d];
    }

    auto readMetadata = [](const conduit::Node& node, const std::string& path, auto& values) {
      if(!node.has_path(path))
      {
        return true;
      }
      const conduit::Node& metadata = node.fetch_existing(path);
      if(metadata.dtype().number_of_elements() != DIM)
      {
        SLIC_ERROR(axom::fmt::format("MarchingCubes metadata '{}' has {} values; expected {}.",
                                     path,
                                     metadata.dtype().number_of_elements(),
                                     DIM));
        return false;
      }
      // The input mesh can reside in device memory. Copy its small metadata
      // array to the host before inspecting it.
      axom::bump::utilities::fillFromNode(node, path, values, true);
      return true;
    };

    if(!readMetadata(n_field, "offsets", fieldOffsets) ||
       !readMetadata(n_field, "strides", fieldStrides) ||
       !readMetadata(n_topo, "elements/dims/offsets", topoOffsets) ||
       !readMetadata(n_topo, "elements/dims/strides", topoStrides))
    {
      return;
    }

    bool iFastest = fieldStrides[0] == 1;
    for(int d = 1; d < DIM; ++d)
    {
      iFastest = iFastest && (fieldStrides[d] >= fieldStrides[d - 1]);
    }
    if(!iFastest)
    {
      SLIC_ERROR(
        axom::fmt::format("MarchingCubes (bump backend) requires an i-fastest function field: bump "
                          "reads field values as a flat array and does not honor Blueprint field "
                          "strides, so a permuted layout is silently transposed. Field '{}' has "
                          "strides that are not i-fastest.",
                          fieldPath));
      return;
    }

    SLIC_ERROR_IF(
      fieldOffsets != topoOffsets || fieldStrides != topoStrides,
      axom::fmt::format("MarchingCubes (bump backend) requires function field '{}' to use the "
                        "same offsets and strides as its structured topology. bump indexes "
                        "field values directly with topology node ids.",
                        fieldPath));
  }

  //! @brief Require a float64 Blueprint array, naming the offending type if not.
  void validateFieldIsFloat64(const std::string& path, const std::string& what) const
  {
    if(m_dom == nullptr || !m_dom->has_path(path))
    {
      return;  // absence is reported elsewhere, with a better message
    }
    const conduit::Node& n = m_dom->fetch_existing(path);
    SLIC_ERROR_IF(!n.dtype().is_float64(),
                  axom::fmt::format("MarchingCubes (bump backend) requires a float64 {} at '{}', "
                                    "but found '{}'.",
                                    what,
                                    path,
                                    n.dtype().name()));
  }

  void setContourValue(double contourVal) override { m_contourVal = contourVal; }

  void setMaskValue(int maskVal) override { m_maskVal = maskVal; }

  /*!
   * @brief Record the requested robustness policy (Phase 6 seam).
   *
   * Currently advisory: `standard` and `robust` both run bump's default intersector + tables.
   * When a robust (double-precision, +/-/0-aware) intersector is available,
   * runExtraction() will select it for `robust` with no change to the calling code.
   */
  void setRobustnessPolicy(MarchingCubesRobustnessPolicy policy) override
  {
    m_robustnessPolicy = policy;
  }

  // The data-parallelism knob is a legacy-kernel concept; bump manages its own
  // parallelism.  We accept and ignore it (kept for API compatibility).
  void setDataParallelism(MarchingCubesDataParallelism dataPar) override
  {
    m_dataParallelism = dataPar;
  }

  // ---- Phased interface (mark/scan/compute) -------------------------------
  // bump does extraction in a single execute() call.
  // We run it lazily in runExtraction() and have the phase methods drive/observe that.

  //! @brief No-op for the bump backend (extraction is deferred).
  void markCrossings() override { /* no-op: deferred to computeFacets */ }

  /*!
   * @brief Run the bump extraction so the facet count is known.
   *
   * The parent MarchingCubes allocates the shared output buffers after the scan phase
   * (it needs per-domain counts to size them) and before the compute phase.
   *
   * bump cannot give us a count without doing the full extraction, so we perform extraction here and cache the result.
   * The count then becomes available to the parent, and computeFacets() copies cached data into the buffers the parent allocated.
   */
  void scanCrossings() override
  {
    m_extractionRan = true;
    runExtraction();
  }

  //! @brief Copy cached bump output into the parent-allocated output buffers.
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
                  "MarchingCubes bump backend has no Blueprint contour output. "
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
                  "MarchingCubes bump backend has no Blueprint contour output. "
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
   * @brief Convert the requested isovalue to bump's field type while preserving
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
    namespace bumpviews = axom::bump::views;

    const std::string cstype = n_coords.fetch_existing("type").as_string();
    if(cstype == "uniform")
    {
      auto coordsetView = bumpviews::make_uniform_coordset<DIM>::view(n_coords);
      func(coordsetView);
    }
    else if(cstype == "rectilinear")
    {
      const conduit::Node& values = n_coords.fetch_existing("values");
      if constexpr(DIM == 2)
      {
        SLIC_ERROR_IF(values.number_of_children() != 2,
                      "2D rectilinear coordsets require 2 component arrays.");
        bumpviews::floatNodeToArrayViewSame(values[0], values[1], [&](auto xView, auto yView) {
          bumpviews::RectilinearCoordsetView2<typename decltype(xView)::value_type> coordsetView(
            xView,
            yView);
          func(coordsetView);
        });
      }
      else
      {
        SLIC_ERROR_IF(values.number_of_children() != 3,
                      "3D rectilinear coordsets require 3 component arrays.");
        bumpviews::floatNodeToArrayViewSame(
          values[0],
          values[1],
          values[2],
          [&](auto xView, auto yView, auto zView) {
            bumpviews::RectilinearCoordsetView3<typename decltype(xView)::value_type> coordsetView(
              xView,
              yView,
              zView);
            func(coordsetView);
          });
      }
    }
    else if(cstype == "explicit")
    {
      const conduit::Node& values = n_coords.fetch_existing("values");
      if constexpr(DIM == 2)
      {
        SLIC_ERROR_IF(values.number_of_children() != 2,
                      "2D explicit coordsets require 2 component arrays.");
        bumpviews::floatNodeToArrayViewSame(values[0], values[1], [&](auto xView, auto yView) {
          bumpviews::ExplicitCoordsetView<typename decltype(xView)::value_type, 2> coordsetView(
            xView,
            yView);
          func(coordsetView);
        });
      }
      else
      {
        SLIC_ERROR_IF(values.number_of_children() != 3,
                      "3D explicit coordsets require 3 component arrays.");
        bumpviews::floatNodeToArrayViewSame(
          values[0],
          values[1],
          values[2],
          [&](auto xView, auto yView, auto zView) {
            bumpviews::ExplicitCoordsetView<typename decltype(xView)::value_type, 3> coordsetView(
              xView,
              yView,
              zView);
            func(coordsetView);
          });
      }
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Unsupported coordset type '{}'.", cstype));
    }
  }

  /*! @brief Dispatch a topology view restricted to MarchingCubes-supported shapes. */
  template <typename FuncType>
  static void dispatchTopology(const conduit::Node& n_topo, FuncType&& func)
  {
    namespace bumpviews = axom::bump::views;

#if defined(_WIN32)
    // Windows shared-library builds auto-export template instantiations from
    // axom_quest.dll.  Keep this opt-in bump path narrow enough to link there,
    // while preserving the generic bump dispatcher on other platforms.
    const std::string topoType = n_topo.fetch_existing("type").as_string();
    if(topoType == "unstructured")
    {
      const std::string shape = n_topo.fetch_existing("elements/shape").as_string();
      if constexpr(DIM == 3)
      {
        SLIC_ERROR_IF(shape != "hex",
                      axom::fmt::format("MarchingCubes bump backend expected "
                                        "unstructured hex topology, but got '{}'.",
                                        shape));
        using ShapeType = bumpviews::HexShape<axom::IndexType>;
        auto topologyView =
          bumpviews::make_unstructured_single_shape_topology<ShapeType>::view(n_topo);
        func(shape, topologyView);
      }
      else
      {
        SLIC_ERROR_IF(shape != "quad",
                      axom::fmt::format("MarchingCubes bump backend expected "
                                        "unstructured quad topology, but got '{}'.",
                                        shape));
        using ShapeType = bumpviews::QuadShape<axom::IndexType>;
        auto topologyView =
          bumpviews::make_unstructured_single_shape_topology<ShapeType>::view(n_topo);
        func(shape, topologyView);
      }
    }
    else if(topoType == "uniform" || topoType == "rectilinear" || topoType == "structured")
    {
      SLIC_ERROR_IF(
        n_topo.has_path("elements/dims/offsets") || n_topo.has_path("elements/dims/strides"),
        "MarchingCubes bump backend does not support strided structured topology "
        "on Windows shared-library builds.");

      const std::string shape = (DIM == 3) ? "hex" : "quad";
      bumpviews::StructuredTopologyView<bumpviews::StructuredIndexing<axom::IndexType, DIM>> topologyView;
      if(topoType == "uniform")
      {
        topologyView = bumpviews::make_uniform_topology<DIM>::view(n_topo);
      }
      else if(topoType == "rectilinear")
      {
        topologyView = bumpviews::make_rectilinear_topology<DIM>::view(n_topo);
      }
      else
      {
        topologyView = bumpviews::make_structured_topology<DIM>::view(n_topo);
      }
      func(shape, topologyView);
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Unsupported topology type '{}'.", topoType));
    }
#else
    bumpviews::dispatch_topology<SelectedDimensions, ShapeTypes>(n_topo,
                                                                 std::forward<FuncType>(func));
#endif
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
    axom::for_all<ExecSpace>(
      nZones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
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
    axom::for_all<ExecSpace>(
      nZones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        if(maskFlagsView[zoneIndex] != 0)
        {
          selectedZonesView[selectedOffsetsView[zoneIndex]] = zoneIndex;
        }
      });

    attachSelectedZonesOption(n_options, selectedZones);
  }

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

    // Copy mask value to a local so the device predicates below capture it by value.
    // AXOM_LAMBDA is [=]; capturing the m_maskVal *member* would instead capture `this`,
    // and dereferencing a host `this` pointer inside a CUDA/HIP kernel is undefined behavior.
    // (Compiles and passes on seq/omp regardless, which is why this must be a local, not the member.)
    const int maskVal = m_maskVal;

    if(m_useMeshViewUtilPath)
    {
      // Structured + explicit: read the mask through MeshViewUtil so any
      // ghost offsets/strides on the field are honored.
      axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(*m_dom, m_topologyName);
      const auto maskView = mvu.template getConstFieldView<int>(m_maskFieldName, false);
      const axom::MDMapping<DIM> topoMap(mvu.getCellShape(), axom::ArrayStrideOrder::COLUMN);

      buildSelectedZonesFromMask(
        nZones,
        [topoMap, maskView, maskVal] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
          const auto zoneIdx = topoMap.toMultiIndex(zoneIndex);
          if constexpr(DIM == 2)
          {
            return maskView(zoneIdx[0], zoneIdx[1]) == maskVal;
          }
          else
          {
            return maskView(zoneIdx[0], zoneIdx[1], zoneIdx[2]) == maskVal;
          }
        },
        n_options,
        selectedZones);
    }
    else
    {
      // Everything else (unstructured, uniform, rectilinear): the mask values
      // are a flat array in the topology's zone order, which matches bump's
      // zone numbering.  That is only true without per-field offsets/strides,
      // so reject those explicitly rather than silently misindexing.
      SLIC_ERROR_IF(n_mask.has_child("offsets") || n_mask.has_child("strides"),
                    "MarchingCubes (bump backend) does not support a mask field with "
                    "Blueprint offsets/strides on a non-structured-explicit topology.");
      auto maskView = bputils::make_array_view<int>(n_maskValues);
      SLIC_ERROR_IF(maskView.size() < nZones,
                    "MarchingCubes mask field has fewer values than topology zones.");
      buildSelectedZonesFromMask(
        nZones,
        [maskView, maskVal] AXOM_HOST_DEVICE(axom::IndexType zoneIndex) {
          return maskView[zoneIndex] == maskVal;
        },
        n_options,
        selectedZones);
    }
  }

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
      AXOM_LAMBDA(axom::IndexType selectedIndex) {
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
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType selectedIndex) {
        if(crossingFlagsView[selectedIndex] != 0)
        {
          crossingZonesView[crossingOffsetsView[selectedIndex]] = selectedZonesView[selectedIndex];
        }
      });

    attachSelectedZonesOption(n_options, crossingZones);
    return crossingCountValue > 0;
  }

  /*!
   * @brief Fast structured crossing pre-filter.
   *
   * @param isoForBump Threshold in the intersector's field type; see isoValueForBump().
   *   The corner test below MUST be the same expression bump's FieldIntersector uses,
   *   or this pre-filter can exclude a zone that bump would have cut (silently dropping facets).
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
        Iterate the logical index space directly rather than deriving it from a flat zone index.

        topoMap.toMultiIndex(zoneIndex) costs DIM integer divisions per zone,
        and this loop runs over EVERY zone, not just crossing ones.
        Nested loops make the flat index incremental and the divisions disappear.

        Node signs are also hoisted: adjacent zones share four (2D) or eight (3D) corners,
        so classifying each NODE once into a byte array and combining bytes per zone
        replaces 2^DIM strided double reads per zone with 2^DIM byte reads.
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

    axom::Array<axom::IndexType> crossingFlags(nZones, nZones, m_allocatorID);
    auto crossingFlagsView = crossingFlags.view();

    // Local copies for device capture: AXOM_LAMBDA is [=], so capturing members
    // would capture `this` and dereference a host pointer in a device kernel.
    const IsoFieldType isoVal = isoForBump;
    const int maskVal = m_maskVal;
    axom::ReduceSum<ExecSpace, axom::IndexType> crossingCount(0);
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

    const axom::IndexType crossingCountValue = crossingCount.get();
    crossingZones =
      axom::Array<axom::IndexType>(crossingCountValue, crossingCountValue, m_allocatorID);

    axom::Array<axom::IndexType> crossingOffsets(nZones, nZones, m_allocatorID);
    auto crossingOffsetsView = crossingOffsets.view();
    axom::exclusive_scan<ExecSpace>(crossingFlagsView, crossingOffsetsView);

    auto crossingZonesView = crossingZones.view();
    axom::for_all<ExecSpace>(
      nZones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        if(crossingFlagsView[zoneIndex] != 0)
        {
          crossingZonesView[crossingOffsetsView[zoneIndex]] = zoneIndex;
        }
      });

    attachSelectedZonesOption(n_options, crossingZones);
    return crossingCountValue > 0;
  }

  /*!
   * @brief Instantiate CutField for (DIM, ExecSpace, this domain's view types)
   * and run it, storing the Blueprint output.
   *
   * Uses bump's dispatch_topology / dispatch_coordset to turn the runtime
   * Blueprint topology+coordset into compile-time view types, then instantiates
   * CutField<ExecSpace, TopoView, CoordView> and calls execute().
   *
   * NOTE: The input domain arrays must already be in a memory space compatible
   * with ExecSpace (the same precondition the legacy backend has).
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
    // Ask bump to record each output facet's originating input zone, which we
    // map onto the legacy "parent cell id" output.
    n_options["originalElementsField"] = kOriginalElementsField;
    // MarchingCubes only consumes the generated originalElements field from CutField.
    // An explicit empty fields map avoids blending/slicing all input fields by default.
    n_options["fields"].set(conduit::DataType::object());

    m_output = std::make_unique<conduit::Node>();
    conduit::Node& n_out = *m_output;

    // Restrict the unstructured shape set to {quad, hex} as requested, to bound template instantiation.
    // Structured dimensions restricted to DIM. Dispatch coordset, then topology, building the matching views
    // and running CutField.  The double dispatch yields the concrete (CoordView, TopoView) pair at compile time.
    bool extracted = false;
    dispatchCoordset(n_coords, [&](auto coordsetView) {
      using CoordsetView = decltype(coordsetView);
      dispatchTopology(n_topo, [&](const std::string& AXOM_UNUSED_PARAM(shape), auto topologyView) {
        using TopologyView = decltype(topologyView);

        // --- Phase 6 robustness seam --------------------------------------
        // The intersector policy is the single point that determines per-cell topology + crossing precision.
        // bump's default FieldIntersector is float-precision and two-label (no +/-/0).
        // When a robust intersector (double precision, +/-/0 / asymptotic-decider aware) is added to bump,
        // alias `Cut` to CutField<..., RobustIntersector> in the `robust` branch below;
        // no other quest code changes.  Until then both branches use the default, and selecting `robust` emits a one-time note.
        using StandardCut = bumpx::CutField<ExecSpace, TopologyView, CoordsetView>;
        // using RobustCut = bumpx::CutField<ExecSpace, TopologyView, CoordsetView,
        //   axom::bump::extraction::RobustFieldIntersector<ExecSpace, TopologyView, CoordsetView>>;
        using Cut = StandardCut;

        if(m_robustnessPolicy == MarchingCubesRobustnessPolicy::robust)
        {
          static bool warnedOnce = false;
          if(!warnedOnce)
          {
            warnedOnce = true;
            SLIC_INFO(
              "MarchingCubes: robust isosurface policy requested, but a robust "
              "bump intersector is not yet available; using the standard "
              "(single-precision, two-label) intersector.");
          }
        }

        Cut iso(topologyView, coordsetView);
        iso.setAllocatorID(m_allocatorID);

        // --- Corner-classification convention ------------------------------
        // Reconcile bump's strict corner test with the legacy kernel's `>=`.
        // Both the value handed to bump AND the structured pre-filter's own
        // test must use this, or the pre-filter and the extractor disagree.
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

          // Restore the public field name before the adaptor or any caller sees the output.
          // The private request name prevents bump from forwarding a same-named field from the input mesh.
          const std::string privateField = axom::fmt::format("fields/{}", kOriginalElementsField);
          if(n_out.has_path(privateField))
          {
            n_out["fields"].rename_child(kOriginalElementsField, kPublicOriginalElementsField);
          }
        }
        extracted = true;
      });
    });

    if(!extracted)
    {
      // An out-of-range isovalue or an empty mask is valid and produces an available, empty contour
      // rather than an empty Blueprint node.
      m_output.reset();
      m_facetCount = 0;
      return;
    }

    // Determine the facet count from the bump output.  After fan-triangulation (see fillLegacyOutputBuffers)
    // the legacy facet count is the number of triangles/segments, not the number of bump polygons;
    // compute it from the output element sizes.
    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesBumpImpl::computeTriangulatedFacetCount");
      m_facetCount = computeTriangulatedFacetCount(n_out);
    }
  }

  /*!
   * @brief Number of (DIM-cornered) facets after fan-triangulating bump output.
   *
   * For DIM==2 the bump output elements are 2-node segments -> 1 facet each.
   * For DIM==3 a p-gon fans into (p-2) triangles.
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
      "MarchingCubes bump backend: cut output topology '{}' has no 'sizes' array. "
      "The adaptor requires explicit sizes on bump's cut output.",
      newTopoName));
    return 0;
  }

  /*!
   * @brief Fill the parent-allocated legacy output buffers from cached bump output,
   * triangulating the polygons while reusing bump's welded vertex coordinates.
   */
  void fillLegacyOutputBuffers()
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesBumpImpl::fillLegacyOutputBuffers");
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

  //! @brief Return the (single) topology name present in a bump output node.
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

  MarchingCubesRobustnessPolicy m_robustnessPolicy {MarchingCubesRobustnessPolicy::standard};

  //! @brief Whether the MeshViewUtil fast paths apply (structured + explicit only).
  bool m_useMeshViewUtilPath {false};

  //! @brief Cached bump CutField output (Blueprint mesh).
  std::unique_ptr<conduit::Node> m_output;

  //! @brief Legacy facet count (post fan-triangulation).
  axom::IndexType m_facetCount {};

  //! @brief Whether extraction ran, distinguishing unavailable from empty output.
  bool m_extractionRan {false};
};

}  // namespace axom::quest::detail::marching_cubes
