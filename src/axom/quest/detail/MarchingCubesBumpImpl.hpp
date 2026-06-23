// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file MarchingCubesBumpImpl.hpp
 *
 * @brief A MarchingCubesSingleDomain::ImplBase implementation that delegates
 * isocontour extraction to axom::bump::extraction::CutField.
 *
 * Unlike the legacy MarchingCubesImpl (which contains a hand-written marching
 * cubes kernel over *structured* data), this implementation wraps the bump
 * CutField extractor.  CutField is templated on <ExecSpace, TopologyView,
 * CoordsetView>, so it transparently supports:
 *   - structured (uniform / rectilinear / explicit-structured) topologies, and
 *   - unstructured *single-shape* quad (2D) and hex (3D) topologies,
 * on all execution spaces (seq, omp, cuda, hip), via the bump view dispatch.
 *
 * Design notes (see MARCHINGCUBES_ON_BUMP.md):
 *  - CutField is a Blueprint-in / Blueprint-out operation.  We run it once per
 *    domain in computeFacets(); markCrossings()/scanCrossings() do the cheap
 *    bookkeeping the parent MarchingCubes orchestration expects.
 *  - The phased ImplBase interface (mark/scan/compute) was designed around the
 *    legacy kernel.  bump does everything in one execute() call, so we run the
 *    extractor lazily and cache its result, then satisfy the count queries from
 *    the cached result.
 *  - bump produces a WELDED, topologically-connected surface (blend-group
 *    uniquification).  The 3D output may be polygonal (tri/quad/poly5..8), and
 *    the 2D output is line segments.  To preserve the legacy fixed-stride
 *    (facetCount, DIM) "triangle soup" output contract, the adaptor
 *    fan-triangulates polygons and re-expands welded points into per-facet
 *    corners when filling the legacy output buffers.  (Richer, welded output is
 *    exposed via additive accessors on MarchingCubes; see Phase 3.)
 */

#ifndef AXOM_QUEST_MARCHINGCUBESBUMPIMPL_H_
#define AXOM_QUEST_MARCHINGCUBESBUMPIMPL_H_

#include "axom/config.hpp"

#if defined(AXOM_USE_CONDUIT) && defined(AXOM_ENABLE_BUMP)

  #include "axom/core/execution/execution_space.hpp"
  #include "axom/core/execution/for_all.hpp"
  #include "axom/core/MDMapping.hpp"
  #include "axom/slic/interface/slic_macros.hpp"
  #include "axom/quest/MeshViewUtil.hpp"
  #include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
  #include "axom/quest/detail/MarchingCubesBumpAdaptor.hpp"

  // bump extraction + views
  #include "axom/bump/extraction/CutField.hpp"
  #include "axom/bump/views/dispatch_coordset.hpp"
  #include "axom/bump/views/dispatch_topology.hpp"
  #include "axom/bump/views/Shapes.hpp"
  #include "axom/bump/utilities/conduit_memory.hpp"

  #include "conduit_node.hpp"
  #include "conduit_blueprint.hpp"

  #include <memory>
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

    // Validate that this is a topology bump+MarchingCubes supports: a DIM-
    // dimensional structured topology, or an unstructured single-shape quad
    // (DIM==2) / hex (DIM==3) topology.  Mixed/polyhedral are rejected here.
    const conduit::Node& n_topo =
      dom.fetch_existing(axom::fmt::format("topologies/{}", topologyName));
    const std::string topoType = n_topo.fetch_existing("type").as_string();
    m_isStructured = (topoType != "unstructured");
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
      SLIC_ERROR_IF(topoType != "uniform" && topoType != "rectilinear" &&
                      topoType != "structured",
                    axom::fmt::format("MarchingCubes (bump backend) does not support "
                                      "topology type '{}'.",
                                      topoType));
    }
  }

  void setFunctionField(const std::string& fcnFieldName) override
  {
    m_fcnFieldName = fcnFieldName;
  }

  void setContourValue(double contourVal) override { m_contourVal = contourVal; }

  void setMaskValue(int maskVal) override { m_maskVal = maskVal; }

  /*!
   * @brief Honor the requested parent-cell-id numbering.
   *
   * blueprintZoneId (default): use bump's originalElements directly.
   * legacyFieldOrder: for structured input, remap the Blueprint zone index
   *   (which bump orders i-fastest, independent of memory layout) to the legacy
   *   flat index in the function field's stride order.  For unstructured input
   *   this mode is ignored (no canonical "field stride order" exists) and the
   *   Blueprint zone id is used.
   */
  void setParentCellIdMode(MarchingCubesParentCellIdMode mode) override
  {
    m_parentCellIdMode = mode;
  }

  // The data-parallelism knob is a legacy-kernel concept; bump manages its own
  // parallelism.  We accept and ignore it (kept for API compatibility).
  void setDataParallelism(MarchingCubesDataParallelism dataPar) override
  {
    m_dataParallelism = dataPar;
  }

  // ---- Phased interface (mark/scan/compute) -------------------------------
  // bump does extraction in a single execute() call.  We run it lazily in
  // runExtraction() and have the phase methods drive/observe that.

  //! @brief No-op for the bump backend (extraction is deferred).
  void markCrossings() override { /* no-op: deferred to computeFacets */ }

  /*!
   * @brief Run the bump extraction so the facet count is known.
   *
   * The parent MarchingCubes allocates the shared output buffers *after* the
   * scan phase (it needs per-domain counts to size them) and *before* the
   * compute phase.  bump cannot give us a count without doing the full
   * extraction, so we perform extraction here and cache the result; the count
   * then becomes available to the parent, and computeFacets() simply copies
   * cached data into the buffers the parent allocated.
   */
  void scanCrossings() override { runExtraction(); }

  //! @brief Copy cached bump output into the parent-allocated output buffers.
  void computeFacets() override { fillLegacyOutputBuffers(); }

  axom::IndexType getContourCellCount() const override { return m_facetCount; }

  void clearDomain() override
  {
    m_output.reset();
    m_facetCount = 0;
  }

private:
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
    n_options["originalElementsField"] = "originalElements";
    // TODO(masking): plumb m_maskFieldName/m_maskVal to a selectedZones list or
    // an intersector that honors the mask (Phase 2 follow-up).  bump does not
    // take a mask field directly; it takes selectedZones.

    m_output = std::make_unique<conduit::Node>();
    conduit::Node& n_out = *m_output;

    // Restrict the unstructured shape set to {quad, hex} as requested, to bound
    // template instantiation.  Structured dimensions restricted to DIM.
    constexpr int shapeMask =
      (DIM == 3) ? (1 << bumpviews::Hex_ShapeID) : (1 << bumpviews::Quad_ShapeID);
    constexpr int dimMask = bumpviews::select_dimensions(DIM);

    // Dispatch coordset, then topology, building the matching views and running
    // CutField.  The double dispatch yields the concrete (CoordView, TopoView)
    // pair at compile time.
    bumpviews::dispatch_coordset(n_coords, [&](auto coordsetView) {
      using CoordsetView = decltype(coordsetView);
      bumpviews::dispatch_topology<dimMask, shapeMask>(
        n_topo,
        [&](const std::string& AXOM_UNUSED_PARAM(shape), auto topologyView) {
          using TopologyView = decltype(topologyView);
          using Cut = bumpx::CutField<ExecSpace, TopologyView, CoordsetView>;
          Cut iso(topologyView, coordsetView);
          iso.setAllocatorID(m_allocatorID);
          iso.execute(*m_dom, n_options, n_out);
        });
    });

    // Determine the facet count from the bump output.  After fan-triangulation
    // (see fillLegacyOutputBuffers) the legacy facet count is the number of
    // triangles/segments, not the number of bump polygons; compute it from the
    // output element sizes.
    m_facetCount = computeTriangulatedFacetCount(n_out);

    // For the opt-in legacyFieldOrder numbering on structured input, capture the
    // logical cell dims and the function field's stride order so the output
    // adaptor can remap bump's i-fastest zone ids back to the legacy ordering.
    if(m_parentCellIdMode == MarchingCubesParentCellIdMode::legacyFieldOrder &&
       m_isStructured)
    {
      captureStructuredMetadata();
    }
  }

  /*!
   * @brief Populate m_cellDims and m_fieldSlowestDirs from the structured
   * domain, for the legacyFieldOrder parent-id remap.
   *
   * Uses MeshViewUtil to read the logical cell shape and the function field's
   * strides, from which an MDMapping yields the slowest->fastest permutation.
   *
   * NOTE: This path is exercised only when a caller explicitly opts into
   * legacyFieldOrder on structured input; it is the part of the bump backend
   * most in need of build/test validation (MeshViewUtil templating x ExecSpace).
   */
  void captureStructuredMetadata()
  {
    axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(*m_dom, m_topologyName);
    m_cellDims = mvu.getCellShape();
    const auto fcnView = mvu.template getConstFieldView<double>(m_fcnFieldName, false);
    // Build an MDMapping from the field strides to extract the stride order.
    axom::MDMapping<DIM> fcnMap(fcnView.strides());
    m_fieldSlowestDirs = fcnMap.slowestDirs();
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
      const auto sizes = n_sizes.as_index_t_accessor();
      const conduit::index_t n = sizes.number_of_elements();
      axom::IndexType facets = 0;
      for(conduit::index_t i = 0; i < n; ++i)
      {
        const auto p = static_cast<axom::IndexType>(sizes[i]);
        facets += (DIM == 3) ? (p >= 3 ? p - 2 : 0) : 1;
      }
      return facets;
    }

    // Fixed-shape output (e.g. all-tri or all-segment): infer count from
    // connectivity length / corners-per-element.
    const conduit::Node& n_conn = n_elems.fetch_existing("connectivity");
    const auto cornersPerElem = (DIM == 3) ? 3 : 2;  // tri or segment
    return static_cast<axom::IndexType>(n_conn.dtype().number_of_elements() / cornersPerElem);
  }

  /*!
   * @brief Fill the parent-allocated legacy output buffers from cached bump
   * output, fan-triangulating polygons and re-expanding welded points so the
   * legacy (facetCount, DIM) un-welded contract is preserved exactly.
   *
   * Implemented in Phase 1.3 (host path) and Phase 1.4 (device path).  Left as
   * a clearly-marked stub here so the scaffolding compiles and the design is
   * reviewable before the data-movement code lands.
   */
  void fillLegacyOutputBuffers()
  {
    SLIC_ASSERT(m_output != nullptr);

    // Build the legacy field-stride remap only when the user asked for the
    // legacy numbering AND the input is structured (unstructured has no
    // canonical field stride order; we leave the remap empty -> pass-through).
    axom::Array<axom::IndexType> remapHost(0, 0);
    if(m_parentCellIdMode == MarchingCubesParentCellIdMode::legacyFieldOrder &&
       m_isStructured)
    {
      remapHost = buildFieldStrideRemap<DIM>(m_cellDims, m_fieldSlowestDirs);
    }

    // Move the (possibly empty) remap into ExecSpace memory for the kernel.
    axom::Array<axom::IndexType> remapDevice;
    axom::ArrayView<const axom::IndexType, 1> remapView;
    if(!remapHost.empty())
    {
      remapDevice = axom::Array<axom::IndexType>(remapHost, m_allocatorID);
      remapView = remapDevice.view();
    }

    adaptCutFieldOutput<DIM, ExecSpace>(*m_output,
                                        m_facetNodeIds,
                                        m_facetNodeCoords,
                                        m_facetParentIds,
                                        m_facetIndexOffset,
                                        m_facetCount,
                                        remapView);
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

  const conduit::Node* m_dom = nullptr;
  std::string m_topologyName;
  std::string m_fcnFieldName;
  std::string m_maskFieldName;

  //! @brief How to number parent-cell ids of generated facets.
  MarchingCubesParentCellIdMode m_parentCellIdMode = MarchingCubesParentCellIdMode::blueprintZoneId;

  //! @name Structured metadata, captured only for the legacyFieldOrder remap.
  //! @{
  bool m_isStructured = false;
  axom::StackArray<axom::IndexType, DIM> m_cellDims {};
  axom::StackArray<std::uint16_t, DIM> m_fieldSlowestDirs {};
  //! @}

  //! @brief Cached bump CutField output (Blueprint mesh).
  std::unique_ptr<conduit::Node> m_output;

  //! @brief Legacy facet count (post fan-triangulation).
  axom::IndexType m_facetCount = 0;
};

}  // namespace marching_cubes
}  // namespace detail
}  // namespace quest
}  // namespace axom

#endif  // AXOM_USE_CONDUIT && AXOM_ENABLE_BUMP
#endif  // AXOM_QUEST_MARCHINGCUBESBUMPIMPL_H_
