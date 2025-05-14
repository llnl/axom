// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_ELVIRA_ALGORITHM_HPP_
#define AXOM_MIR_ELVIRA_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/utilities/ExtractZones.hpp"
#include "axom/mir/utilities/MakeZoneCenters.hpp"
#include "axom/mir/utilities/MergeMeshes.hpp"
#include "axom/mir/utilities/PrimalAdaptor.hpp"
#include "axom/mir/utilities/SelectedZones.hpp"
#include "axom/mir/utilities/ZoneListBuilder.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"
#include "axom/mir/utilities/utilities.hpp"
#include "axom/mir/views/MaterialView.hpp"
#include "axom/mir/views/StructuredTopologyView.hpp"

#include "axom/primal/operators/compute_bounding_box.hpp"
#include "axom/primal/operators/clip.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <algorithm>
#include <string>

// Uncomment to save inputs and outputs.
//#define AXOM_ELVIRA_DEBUG

// Uncomment to debug make fragments.
//#define AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS

// Uncomment to gather ELVIRA data and save to YAML file.
//#define AXOM_ELVIRA_GATHER_INFO

#if defined(AXOM_ELVIRA_DEBUG)
  #include <conduit/conduit_relay_io.hpp>
  #include <conduit/conduit_relay_io_blueprint.hpp>
#endif

#include "axom/mir/detail/elvira_impl.hpp"
#include "axom/mir/detail/elvira_detail.hpp"

namespace axom
{
namespace mir
{

//------------------------------------------------------------------------------
/*!
 * \brief Implements Elvira algorithm for structured meshes.
 */
template <typename ExecSpace, typename IndexPolicy, typename CoordsetView, typename MatsetView>
class ElviraAlgorithm : public axom::mir::MIRAlgorithm
{
public:
  using TopologyView = axom::mir::views::StructuredTopologyView<IndexPolicy>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

protected:
  static constexpr int NDIMS = CoordsetView::dimension();
  static constexpr int StencilSize = elvira::getStencilSize(NDIMS);
  static constexpr int numVectorComponents = 3;

  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;

  // Determine the output type from the clip operations. Those are the shape
  // types that we're emitting into the MIR output. Create the builder.
  using CoordType = typename CoordsetView::value_type;
  using ClipResultType =
    typename std::conditional<NDIMS == 2,
                              axom::primal::Polygon<CoordType, 2, axom::primal::PolygonArray::Static>,
                              axom::primal::Polyhedron<CoordType, 3>>::type;

  using VectorType = axom::primal::Vector<CoordType, NDIMS>;
  using PointType = axom::primal::Point<CoordType, NDIMS>;
  using PlaneType = axom::primal::Plane<CoordType, NDIMS>;

  using ShapeView = axom::mir::utilities::blueprint::PrimalAdaptor<TopologyView, CoordsetView>;
  using Builder =
    detail::TopologyBuilder<ExecSpace, CoordsetView, TopologyView, MatsetView, ClipResultType, NDIMS>;
  using BuilderView = typename Builder::View;

public:
  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  ElviraAlgorithm(const TopologyView &topoView,
                  const CoordsetView &coordsetView,
                  const MatsetView &matsetView)
    : axom::mir::MIRAlgorithm()
    , m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~ElviraAlgorithm() = default;

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  /*!
   * \brief Perform material interface reconstruction on a single domain.
   *
   * \param[in] n_topo The Conduit node containing the topology that will be used for MIR.
   * \param[in] n_coordset The Conduit node containing the coordset.
   * \param[in] n_fields The Conduit node containing the fields.
   * \param[in] n_matset The Conduit node containing the matset.
   * \param[in] n_options The Conduit node containing the options that help govern MIR execution.
   *
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   * \param[out] n_newMatset A Conduit node that will contain the new matset.
   * 
   */
  virtual void executeDomain(const conduit::Node &n_topo,
                             const conduit::Node &n_coordset,
                             const conduit::Node &n_fields,
                             const conduit::Node &n_matset,
                             const conduit::Node &n_options,
                             conduit::Node &n_newTopo,
                             conduit::Node &n_newCoordset,
                             conduit::Node &n_newFields,
                             conduit::Node &n_newMatset) override
  {
    namespace bputils = axom::mir::utilities::blueprint;

    AXOM_ANNOTATE_SCOPE("ElviraAlgorithm");

    // Copy the options to make sure they are in the right memory space.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

    // _mir_utilities_selectedzones_begin
    // Get selected zones from the options.
    bputils::SelectedZones<ExecSpace> selectedZones(m_topologyView.numberOfZones(), n_options_copy);
    const auto selectedZonesView = selectedZones.view();
    // _mir_utilities_selectedzones_end

    // Partition the selected zones into clean, mixed lists.
    axom::Array<axom::IndexType> cleanZones, mixedZones;
    bputils::ZoneListBuilder<ExecSpace, TopologyView, MatsetView> zlb(m_topologyView, m_matsetView);
    zlb.execute(selectedZonesView, cleanZones, mixedZones);
    SLIC_ASSERT((cleanZones.size() + mixedZones.size()) == selectedZonesView.size());
    SLIC_INFO(
      axom::fmt::format("cleanZones: {}, mixedZones: {}", cleanZones.size(), mixedZones.size()));

    if(cleanZones.size() > 0 && mixedZones.size() > 0)
    {
      // Gather the inputs into a single root but replace the fields with
      // a new node to which we can add additional fields.
      conduit::Node n_root;
      n_root[n_coordset.path()].set_external(n_coordset);
      n_root[n_topo.path()].set_external(n_topo);
      n_root[n_matset.path()].set_external(n_matset);
      conduit::Node &n_root_coordset = n_root[n_coordset.path()];
      conduit::Node &n_root_topo = n_root[n_topo.path()];
      conduit::Node &n_root_matset = n_root[n_matset.path()];
      conduit::Node n_root_fields = n_root["fields"];

      // Make the clean mesh.
      conduit::Node n_cleanOutput;
      makeCleanOutput(n_root, n_topo.name(), n_options_copy, cleanZones.view(), n_cleanOutput);

      // Process the mixed part of the mesh.
      processMixedZones(mixedZones.view(),
                        n_root_topo,
                        n_root_coordset,
                        n_root_fields,
                        n_root_matset,
                        n_options_copy,
                        n_newTopo,
                        n_newCoordset,
                        n_newFields,
                        n_newMatset);

      // Gather the MIR output into a single node.
      conduit::Node n_mirOutput;
      n_mirOutput[n_newTopo.path()].set_external(n_newTopo);
      n_mirOutput[n_newCoordset.path()].set_external(n_newCoordset);
      n_mirOutput[n_newFields.path()].set_external(n_newFields);
      n_mirOutput[n_newMatset.path()].set_external(n_newMatset);
#if defined(AXOM_ELVIRA_DEBUG)
      saveMesh(n_mirOutput, "debug_elvira_mir");
      SLIC_DEBUG("--- clean ---");
      printNode(n_cleanOutput);
      SLIC_DEBUG("--- MIR ---");
      printNode(n_mirOutput);
#endif

      // Merge the clean zones and MIR output
      conduit::Node n_merged;
      merge(n_newTopo.name(), n_cleanOutput, n_mirOutput, n_merged);
#if defined(AXOM_ELVIRA_DEBUG)
      SLIC_DEBUG("--- merged ---");
      printNode(n_merged);

      // Save merged output.
      saveMesh(n_merged, "debug_elvira_merged");
#endif

      // Move the merged output into the output variables.
      n_newCoordset.move(n_merged[n_newCoordset.path()]);
      n_newTopo.move(n_merged[n_newTopo.path()]);
      n_newFields.move(n_merged[n_newFields.path()]);
      n_newMatset.move(n_merged[n_newMatset.path()]);
    }
    else if(cleanZones.size() == 0 && mixedZones.size() > 0)
    {
      // Only mixed zones.
      processMixedZones(mixedZones.view(),
                        n_topo,
                        n_coordset,
                        n_fields,
                        n_matset,
                        n_options_copy,
                        n_newTopo,
                        n_newCoordset,
                        n_newFields,
                        n_newMatset);
    }
    else if(cleanZones.size() > 0 && mixedZones.size() == 0)
    {
      // There were no mixed zones. We can copy the input to the output.
      {
        AXOM_ANNOTATE_SCOPE("copy");
        bputils::copy<ExecSpace>(n_newCoordset, n_coordset);
        bputils::copy<ExecSpace>(n_newTopo, n_topo);
        bputils::copy<ExecSpace>(n_newFields, n_fields);
        bputils::copy<ExecSpace>(n_newMatset, n_matset);
      }

      // Add an originalElements array.
      const std::string originalElementsField(MIROptions(n_options).originalElementsField());
      addOriginal(n_newFields[originalElementsField],
                  n_topo.name(),
                  "element",
                  m_topologyView.numberOfZones());
    }
  }

  /*!
   * \brief Merge meshes for clean and MIR outputs.
   *
   * \param topoName The name of the topology.
   * \param n_cleanOutput The mesh that contains the clean zones.
   * \param n_mirOutput The mesh that contains the MIR output.
   * \param[out] n_merged The output node for the merged mesh.
   */
  void merge(const std::string &topoName,
             conduit::Node &n_cleanOutput,
             conduit::Node &n_mirOutput,
             conduit::Node &n_merged) const
  {
    AXOM_ANNOTATE_SCOPE("merge");
    namespace bputils = axom::mir::utilities::blueprint;

    // Create a MergeMeshesAndMatsets type that will operate on the material
    // inputs, which at this point will be unibuffer with known types. We can
    // reduce code bloat and compile time by passing a MaterialDispatch policy.
    using IntElement = typename MatsetView::IndexType;
    using FloatElement = typename MatsetView::FloatType;
    constexpr size_t MAXMATERIALS = MatsetView::MaxMaterials;
    using DispatchPolicy =
      bputils::DispatchTypedUnibufferMatset<IntElement, FloatElement, MAXMATERIALS>;
    using MergeMeshes = bputils::MergeMeshesAndMatsets<ExecSpace, DispatchPolicy>;

    // Merge clean and MIR output.
    std::vector<bputils::MeshInput> inputs(2);
    inputs[0].m_input = &n_cleanOutput;
    inputs[0].topologyName = topoName;

    inputs[1].m_input = &n_mirOutput;
    inputs[1].topologyName = topoName;

    conduit::Node mmOpts;
    mmOpts["topologyName"] = topoName;
    MergeMeshes mm;
    mm.execute(inputs, mmOpts, n_merged);
  }

  /*!
   * \brief Adds original ids field to supplied fields node.
   *
   * \param n_field The new field node.
   * \param topoName The topology name for the field.
   * \param association The field association.
   * \param nvalues The number of nodes in the field.
   *
   */
  void addOriginal(conduit::Node &n_field,
                   const std::string &topoName,
                   const std::string &association,
                   axom::IndexType nvalues) const
  {
    AXOM_ANNOTATE_SCOPE("addOriginal");
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Add a new field for the original ids.
    n_field["topology"] = topoName;
    n_field["association"] = association;
    n_field["values"].set_allocator(c2a.getConduitAllocatorID());
    n_field["values"].set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, nvalues));
    auto view = bputils::make_array_view<ConnectivityType>(n_field["values"]);
    axom::for_all<ExecSpace>(
      nvalues,
      AXOM_LAMBDA(axom::IndexType index) { view[index] = static_cast<ConnectivityType>(index); });
  }

  /*!
   * \brief Take the mesh in n_root and extract the zones identified by the
   *        \a cleanZones array and store the results into the \a n_cleanOutput
   *        node.
   *
   * \param n_root The input mesh from which zones are being extracted.
   * \param topoName The name of the topology.
   * \param n_options The options to inherit.
   * \param cleanZones An array of clean zone ids.
   * \param[out] n_cleanOutput The node that will contain the clean mesh output.
   *
   * \return The number of nodes in the clean mesh output.
   */
  void makeCleanOutput(const conduit::Node &n_root,
                       const std::string &topoName,
                       const conduit::Node &n_options,
                       const axom::ArrayView<axom::IndexType> &cleanZones,
                       conduit::Node &n_cleanOutput) const
  {
    AXOM_ANNOTATE_SCOPE("makeCleanOutput");
    namespace bputils = axom::mir::utilities::blueprint;

    // Make the clean mesh.
    bputils::ExtractZonesAndMatset<ExecSpace, TopologyView, CoordsetView, MatsetView> ez(
      m_topologyView,
      m_coordsetView,
      m_matsetView);
    conduit::Node n_ezopts;
    n_ezopts["topology"] = topoName;
    n_ezopts["compact"] = 1;
    n_ezopts["originalElementsField"] = Options(n_options).originalElementsField();
    // Forward some options involved in naming the objects.
    const std::vector<std::string> keys {"topologyName", "coordsetName", "matsetName"};
    for(const auto &key : keys)
    {
      if(n_options.has_path(key))
      {
        n_ezopts[key].set(n_options[key]);
      }
    }
    ez.execute(cleanZones, n_root, n_ezopts, n_cleanOutput);
#if defined(AXOM_ELVIRA_DEBUG)
    {
      AXOM_ANNOTATE_SCOPE("saveClean");
      saveMesh(n_cleanOutput, "debug_elvira_clean");
    }
#endif
  }

  /*!
   * \brief Perform material interface reconstruction on mixed zones.
   *
   * \param[in] mixedZonesView A view that contains the list of mixed zones that we'll operate on.
   * \param[in] n_topo The Conduit node containing the topology that will be used for MIR.
   * \param[in] n_coordset The Conduit node containing the coordset.
   * \param[in] n_fields The Conduit node containing the fields.
   * \param[in] n_matset The Conduit node containing the matset.
   * \param[in] n_options The Conduit node containing the options that help govern MIR execution.
   *
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   * \param[out] n_newMatset A Conduit node that will contain the new matset.
   * 
   */
  void processMixedZones(const axom::ArrayView<axom::IndexType> mixedZonesView,
                         const conduit::Node &n_topo,
                         const conduit::Node &n_coordset,
                         const conduit::Node &AXOM_UNUSED_PARAM(n_fields),
                         const conduit::Node &n_matset,
                         const conduit::Node &n_options,
                         conduit::Node &n_newTopo,
                         conduit::Node &n_newCoordset,
                         conduit::Node &n_newFields,
                         conduit::Node &n_newMatset)
  {
    AXOM_ANNOTATE_SCOPE("processMixedZones");
    namespace bputils = axom::mir::utilities::blueprint;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    constexpr int NDIMS = TopologyView::dimension();

    // Handle options.
    // When coordinates have float value, we can't necessarily get beyond a
    // certain tolerance so set the tolerance accordingly.
    constexpr double DEFAULT_TOLERANCE = std::is_same<CoordType, float>::value
      ? (axom::numeric_limits<float>::epsilon() * 4.f)
      : 1.e-10;
    constexpr int DEFAULT_MAX_ITERATIONS = 50;
    double tolerance = DEFAULT_TOLERANCE;
    double point_tolerance = DEFAULT_TOLERANCE;
    if(n_options.has_child("tolerance") && n_options["tolerance"].dtype().is_number())
    {
      tolerance = axom::utilities::abs(n_options["tolerance"].to_double());
      if(tolerance < axom::numeric_limits<CoordType>::epsilon())
      {
        tolerance = axom::numeric_limits<CoordType>::epsilon();
      }
      SLIC_ASSERT(tolerance > 0.);
    }
    if(n_options.has_child("point_tolerance") && n_options["point_tolerance"].dtype().is_number())
    {
      point_tolerance = axom::utilities::abs(n_options["point_tolerance"].to_double());
      if(point_tolerance < axom::numeric_limits<CoordType>::epsilon())
      {
        point_tolerance = axom::numeric_limits<CoordType>::epsilon();
      }
      SLIC_ASSERT(point_tolerance > 0.);
    }
    int max_iterations = DEFAULT_MAX_ITERATIONS;
    if(n_options.has_child("max_iterations") && n_options["max_iterations"].dtype().is_number())
    {
      max_iterations = axom::utilities::abs(n_options["max_iterations"].to_int());
      max_iterations = axom::utilities::clampVal(max_iterations, 1, 100);
    }

#if defined(AXOM_ELVIRA_GATHER_INFO)
    // Let's output the normals
    conduit::Node *n_result = new conduit::Node;
#endif

    //--------------------------------------------------------------------------
    // Count the number of fragments we'll make for the mixed zones.
    AXOM_ANNOTATE_BEGIN("counting");

    const auto nzones = mixedZonesView.size();
    axom::Array<axom::IndexType> matCount(nzones, nzones, allocatorID);
    axom::Array<axom::IndexType> matZone(nzones, nzones, allocatorID);
    auto matCountView = matCount.view();
    auto matZoneView = matZone.view();

    // Get the material count per zone and the zone number (in case of strided structured)
    const TopologyView deviceTopologyView(m_topologyView);
    const MatsetView deviceMatsetView(m_matsetView);
    RAJA::ReduceSum<reduce_policy, axom::IndexType> num_reduce(0);
    RAJA::ReduceMax<reduce_policy, axom::IndexType> reduce_maxcuts(0);
    axom::for_all<ExecSpace>(
      mixedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // Get the material data for the zone.
        const auto zoneIndex = mixedZonesView[szIndex];
        const auto matZoneIndex = deviceTopologyView.indexing().LocalToGlobal(zoneIndex);
        const auto nmats = deviceMatsetView.numberOfMaterials(matZoneIndex);

        // Save some material information for later.
        matCountView[szIndex] = nmats;
        matZoneView[szIndex] = zoneIndex;

        // Sum total materials
        num_reduce += nmats;

        // The number of times we cut a zone is the number of materials in the zone minus one.
        reduce_maxcuts.max(nmats - 1);
      });
    const auto numFragments = num_reduce.get();
    const auto maxCuts = reduce_maxcuts.get();
    SLIC_ASSERT(numFragments > 0);
    SLIC_ASSERT(maxCuts > 0);
    SLIC_INFO(
      axom::fmt::format("ElviraAlgorithm: numFragments: {}, maxCuts: {}", numFragments, maxCuts));

#if defined(AXOM_ELVIRA_GATHER_INFO)
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      conduit::Node &n_group1 = n_result->operator[]("group1");
      n_group1["mixedZones"].set(mixedZonesView.data(), mixedZonesView.size());
      n_group1["matZone"].set(matZoneView.data(), matZoneView.size());
      n_group1["matCount"].set(matCountView.data(), matCountView.size());
      n_group1["numFragments"].set(numFragments);
    }
#endif
    AXOM_ANNOTATE_END("counting");

    //--------------------------------------------------------------------------
    // Sort the zones by the mat count. This should make adjacent zones in the
    // list more likely to have the same number of materials.
    AXOM_ANNOTATE_BEGIN("sorting");
    RAJA::stable_sort_pairs<loop_policy>(RAJA::make_span(matCountView.data(), nzones),
                                         RAJA::make_span(matZoneView.data(), nzones));
    AXOM_ANNOTATE_END("sorting");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("offsets");
    axom::Array<axom::IndexType> matOffset(nzones, nzones, allocatorID);
    auto matOffsetView = matOffset.view();
    axom::exclusive_scan<ExecSpace>(matCountView, matOffsetView);
#if defined(AXOM_ELVIRA_GATHER_INFO)
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      conduit::Node &n_group2 = n_result->operator[]("group2");
      n_group2["matZone"].set(matZoneView.data(), matZoneView.size());
      n_group2["matCount"].set(matCountView.data(), matCountView.size());
      n_group2["offsets"].set(matOffsetView.data(), matOffsetView.size());
    }
#endif
    AXOM_ANNOTATE_END("offsets");

    //--------------------------------------------------------------------------
    // NOTE: This makes zone centers for all zones. That's okay since we use
    //       these values in making stencils.
    AXOM_ANNOTATE_BEGIN("centroids");
    // _mir_utilities_makezonecenters_begin
    bputils::MakeZoneCenters<ExecSpace, TopologyView, CoordsetView> zc(m_topologyView,
                                                                       m_coordsetView);
    conduit::Node n_zcfield;
    zc.execute(n_topo, n_coordset, n_zcfield);
    // _mir_utilities_makezonecenters_end
    axom::ArrayView<CoordType> xview, yview, zview;
    xview = bputils::make_array_view<CoordType>(n_zcfield["values/x"]);
    yview = bputils::make_array_view<CoordType>(n_zcfield["values/y"]);
    if(n_zcfield.has_path("values/z"))
    {
      zview = bputils::make_array_view<CoordType>(n_zcfield["values/z"]);
    }
#if defined(AXOM_ELVIRA_GATHER_INFO)
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      conduit::Node &n_group3 = n_result->operator[]("group3");
      n_group3["xview"].set(xview.data(), xview.size());
      n_group3["yview"].set(yview.data(), yview.size());
      n_group3["zview"].set(zview.data(), zview.size());
    }
#endif
    AXOM_ANNOTATE_END("centroids");

    //--------------------------------------------------------------------------
    // Retrieve stencil data for each zone material.
    AXOM_ANNOTATE_BEGIN("stencil");
    const auto numFragmentsStencil = numFragments * StencilSize;

    axom::Array<double> fragmentVFStencil(numFragmentsStencil, numFragmentsStencil, allocatorID);
    auto fragmentVFStencilView = fragmentVFStencil.view();

    // Sorted material ids for each zone.
    axom::Array<typename MatsetView::IndexType> sortedMaterialIds(numFragments,
                                                                  numFragments,
                                                                  allocatorID);
    auto sortedMaterialIdsView = sortedMaterialIds.view();

    // Coordinate stencil data for each zone.
    const auto nzonesStencil = nzones * StencilSize;
    axom::Array<double> xcStencil(nzonesStencil, nzonesStencil, allocatorID);
    axom::Array<double> ycStencil(nzonesStencil, nzonesStencil, allocatorID);
    axom::Array<double> zcStencil(nzonesStencil, nzonesStencil, allocatorID);
    auto xcStencilView = xcStencil.view();
    auto ycStencilView = ycStencil.view();
    auto zcStencilView = zcStencil.view();

    // Traverse the selected zones based on how many materials there are in a zone.
    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        const auto matCount = matCountView[szIndex];
        // The index to use for the zone's material.
        const auto matZoneIndex = deviceTopologyView.indexing().LocalToGlobal(zoneIndex);
        // Where to begin writing this zone's fragment data.
        const auto offset = matOffsetView[szIndex];

        // Get materials for this zone from the matset.
        typename MatsetView::IDList ids;
        typename MatsetView::VFList vfs;
        deviceMatsetView.zoneMaterials(matZoneIndex, ids, vfs);

        // Reverse sort the materials by the volume fraction so the larger VFs are first.
        SLIC_ASSERT(ids.size() == matCount);
        axom::utilities::reverse_sort_multiple(vfs.data(), ids.data(), matCount);

        // Save sorted ids in sortedMaterialIdsView.
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          const auto fragmentIndex = offset + m;
          sortedMaterialIdsView[fragmentIndex] = ids[m];
        }

        // Retrieve the stencil data from neighbor zones.
        auto logical = deviceTopologyView.indexing().IndexToLogicalIndex(zoneIndex);
        for(int si = 0; si < StencilSize; si++)
        {
          // Stencil neighbor logical index.
          typename TopologyView::LogicalIndex neighbor(logical);

          // Neighbor offsets are in (-1, 0, 1) that get added to current zone's logical coordinate.
          const int neighborOffset[3] = {(si % 3) - 1, ((si % 9) / 3) - 1, (si / 9) - 1};
          for(int d = 0; d < NDIMS; d++)
          {
            neighbor[d] += neighborOffset[d];
          }

          // Clamp the neighbor to a zone that is inside the indexing space.
          neighbor = deviceTopologyView.indexing().clamp(neighbor);
          const auto neighborIndex = static_cast<typename MatsetView::ZoneIndex>(
            deviceTopologyView.indexing().LogicalIndexToIndex(neighbor));

          // Turn to a "global" logical index and transform it to an index to use in the material,
          // which for strided-structured can be larger than the mesh.
          const auto matNeighbor = deviceTopologyView.indexing().LocalToGlobal(neighbor);
          const auto matNeighborIndex = static_cast<typename MatsetView::ZoneIndex>(
            deviceTopologyView.indexing().GlobalToGlobal(matNeighbor));

          // Copy material vfs into the stencil.
          for(axom::IndexType m = 0; m < matCount; m++)
          {
            // Ask the neighbor zone for vf.
            const auto fragmentIndex = offset + m;
            typename MatsetView::FloatType vf = 0;

            deviceMatsetView.zoneContainsMaterial(matNeighborIndex,
                                                  sortedMaterialIdsView[fragmentIndex],
                                                  vf);

            // Store the vf into the stencil for the current material.
            const auto destIndex = fragmentIndex * StencilSize + si;
            fragmentVFStencilView[destIndex] = static_cast<double>(vf);
          }

          // The destination index for this coordinate stencil zone.
          const auto coordIndex = szIndex * StencilSize + si;

          // coord stencil
          xcStencilView[coordIndex] = xview[neighborIndex];
          ycStencilView[coordIndex] = yview[neighborIndex];
          zcStencilView[coordIndex] = zview.empty() ? 1. : zview[neighborIndex];
        }
      });
    // We're done with the zone centers.
    n_zcfield.reset();
    AXOM_ANNOTATE_END("stencil");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("vectors");
    const auto vecSize = numFragments * numVectorComponents;
    axom::Array<double> fragmentVectors(vecSize, vecSize, allocatorID);
    auto fragmentVectorsView = fragmentVectors.view();

#if defined(AXOM_ELVIRA_GATHER_INFO)
    conduit::Node *n_group4 = &(n_result->operator[]("group4"));
#endif

    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto matCount = matCountView[szIndex];
        // Where to begin writing this zone's fragment data.
        const auto offset = matOffsetView[szIndex];

        // The following lines are adapted from mira1.c:262

        // Compute Jacobian here since we have coordinate stencil data.
        double jac[3][3];
        const auto coordIndex = szIndex * StencilSize;
        const double *xcStencil = xcStencilView.data() + coordIndex;
        const double *ycStencil = ycStencilView.data() + coordIndex;
        const double *zcStencil = zcStencilView.data() + coordIndex;
        elvira::computeJacobian(xcStencil, ycStencil, zcStencil, NDIMS, jac);

        // The starting addresses for fragments in the current zone.
        const double *fragmentVFStencilStart = fragmentVFStencilView.data() + offset * StencilSize;
        double *fragmentVectorsStart = fragmentVectorsView.data() + offset * numVectorComponents;

        // Produce normal for each material in this zone.
        int iskip = matCount - 1;
        elvira::elvira<NDIMS>::execute(matCount, fragmentVFStencilStart, fragmentVectorsStart, iskip);

#if defined(AXOM_ELVIRA_GATHER_INFO) && !defined(AXOM_DEVICE_CODE)
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        conduit::Node &n_thisZone = n_group4->append();
        n_thisZone["szIndex"] = szIndex;
        n_thisZone["zone"] = zoneIndex;
        n_thisZone["matCount"] = matCount;
        n_thisZone["offset"] = offset;
        n_thisZone["xcStencil"].set(xcStencil, StencilSize);
        n_thisZone["ycStencil"].set(ycStencil, StencilSize);
        n_thisZone["zcStencil"].set(zcStencil, StencilSize);
        n_thisZone["jacobian"].set(&jac[0][0], 9);
        conduit::Node &n_mats = n_thisZone["mats"];
        const double *vf = fragmentVFStencilStart;
        double *n = fragmentVectorsStart;
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          conduit::Node &n_thismat = n_mats.append();
          n_thismat["mat"] = sortedMaterialIdsView[offset + m];
          n_thismat["stencil"].set(vf, StencilSize);
          n_thismat["normal"].set(n, 3);
          vf += StencilSize;
          n += numVectorComponents;
        }
#endif

        // Transform the normals.
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          double *normal = fragmentVectorsView.data() + ((offset + m) * numVectorComponents);
          elvira::transform(normal, jac);
#if defined(AXOM_ELVIRA_GATHER_INFO) && !defined(AXOM_DEVICE_CODE)
          n_thisZone["mats"][m]["transformed_normal"].set(normal, 3);
#endif
        }
      });
    AXOM_ANNOTATE_END("vectors");

    //--------------------------------------------------------------------------
#if defined(AXOM_ELVIRA_GATHER_INFO)
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      AXOM_ANNOTATE_SCOPE("saveInfo");
      conduit::relay::io::save(*n_result, "elvira.yaml", "yaml");
    }
    delete n_result;
#endif

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("allocate");

    // Free stencil data we no longer need.
    xcStencil.clear();
    ycStencil.clear();
    zcStencil.clear();

    // Make the builder that will set up the Blueprint output.
    Builder build;
    build.allocate(numFragments, maxCuts, n_newCoordset, n_newTopo, n_newFields, n_newMatset, n_options);
    if(n_matset.has_path("material_map"))
    {
      n_newMatset["material_map"].set(n_matset["material_map"]);
    }
    AXOM_ANNOTATE_END("allocate");

    //--------------------------------------------------------------------------
    makeFragments(build.view(),
                  matZoneView,
                  matCountView,
                  matOffsetView,
                  sortedMaterialIdsView,
                  fragmentVectorsView,
                  fragmentVFStencilView,
                  max_iterations,
                  tolerance);

    //--------------------------------------------------------------------------
    // Clean up the mesh in 3D so it has merged coordinates and merged faces.
    // This step does nothing in 2D at present.
    axom::Array<axom::IndexType> selectedIds;
    build.cleanMesh(n_newCoordset, point_tolerance, n_newTopo, selectedIds);

    //--------------------------------------------------------------------------
#if defined(AXOM_ELVIRA_DEBUG)
    {
      AXOM_ANNOTATE_SCOPE("verify");
      conduit::Node n_mesh;
      n_mesh[n_newCoordset.path()].set_external(n_newCoordset);
      n_mesh[n_newTopo.path()].set_external(n_newTopo);
      n_mesh[n_newFields.path()].set_external(n_newFields);
      n_mesh[n_newMatset.path()].set_external(n_newMatset);

      // Verify the MIR output.
      conduit::Node info;
      if(!conduit::blueprint::mesh::verify(n_mesh, info))
      {
        info.print();
      }
    }
#endif
  }

  /*!
   * \brief Use the normal vectors for each fragment to bulid the fragment shapes.
   *
   * \param buildView A view that lets us add a shape to the Blueprint output.
   * \param matZoneView A view containing zone ids sorted by material count in the zone.
   * \param matCountView A view containing material counts for the zones in \a matZoneView.
   * \param matOffsetView A view containing the offset where fragment data will be stored.
   * \param sortedMaterialIdsView A view containing the material ids for each fragment sorted by volume fraction (high to low).
   * \param fragmentVectorsView A view that contains the normals for each fragment.
   * \param fragmentVFStencilView A view that contains the VF stencil for each fragment.
   * \param max_iterations The max allowable iterations to find the fragment shape.
   * \param tolerance The volume tolerance for the fragment shape.
   */
  void makeFragments(BuilderView buildView,
                     axom::ArrayView<axom::IndexType> matZoneView,
                     axom::ArrayView<axom::IndexType> matCountView,
                     axom::ArrayView<axom::IndexType> matOffsetView,
                     axom::ArrayView<typename MatsetView::IndexType> sortedMaterialIdsView,
                     axom::ArrayView<double> fragmentVectorsView,
                     axom::ArrayView<double> fragmentVFStencilView,
                     int max_iterations,
                     double tolerance)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("fragments");

    const ShapeView deviceShapeView {m_topologyView, m_coordsetView};
    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        const auto matCount = matCountView[szIndex];
        // Where this zone's fragment data starts.
        const auto offset = matOffsetView[szIndex];

        // Get the starting shape.
        const auto inputShape = deviceShapeView.getShape(zoneIndex);

        // Get the zone's actual volume.
        const double zoneVol = bputils::ComputeShapeAmount<NDIMS>::execute(inputShape);

        ClipResultType remaining;

      // Make a fragment for each material. The biggest ones come first.
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
        SLIC_DEBUG("makeFragments: zoneIndex=" << zoneIndex << ", matCount=" << matCount);
#endif
        for(axom::IndexType m = 0; m < matCount - 1; m++)
        {
          const auto fragmentIndex = offset + m;
          // Get this material fragment's normal and material id.
          const auto matId = sortedMaterialIdsView[fragmentIndex];
          const double *normalPtr =
            fragmentVectorsView.data() + (fragmentIndex * numVectorComponents);

          // Compute the desired fragment volume.
          // Get current material vf from the stencil. (should be faster than material view)
          constexpr int StencilCenter = (NDIMS == 3) ? 13 : ((NDIMS == 2) ? 4 : 1);
          const auto si = fragmentIndex * StencilSize + StencilCenter;
          const auto matVolume = zoneVol * fragmentVFStencilView[si];

          // Make the normal
          VectorType normal;
          for(int d = 0; d < NDIMS; d++)
          {
            normal[d] = static_cast<CoordType>(normalPtr[d]);
          }

          ClipResultType clippedShape;
          PointType range[2];
          PointType pt {};

          if(m == 0)
          {
            // First time through, operate on the inputShape.

            // Compute start and end points along which to move the plane origin.
            detail::computeRange(inputShape, normal, range);
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
            SLIC_DEBUG("\tm=" << m << ", inputShape=" << inputShape << ", range={" << range[0]
                              << ", " << range[1] << "}");
#endif
            // Figure out the clipped shape that has the desired volume.
            clippedShape = detail::clipToVolume<ClipResultType>(inputShape,
                                                                normal,
                                                                range,
                                                                matVolume,
                                                                max_iterations,
                                                                tolerance,
                                                                pt);
          }
          else
          {
            // In subsequent iterations, the clippedShape is the input and it
            // can have a different type than inputShape.

            // Compute start and end points along which to move the plane origin.
            detail::computeRange(remaining, normal, range);
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
            SLIC_DEBUG("\tm=" << m << ", remaining=" << remaining << ", range={" << range[0] << ", "
                              << range[1] << "}");
#endif
            // Figure out the clipped shape that has the desired volume.
            clippedShape = detail::clipToVolume<ClipResultType>(remaining,
                                                                normal,
                                                                range,
                                                                matVolume,
                                                                max_iterations,
                                                                tolerance,
                                                                pt);
          }

          // Emit clippedShape as material matId
          buildView.addShape(zoneIndex, fragmentIndex, clippedShape, matId, normalPtr);

          // Clip in the other direction to get the remaining fragment for the next material.
          const auto P = PlaneType(normal, pt, false);
          if(m == 0)
          {
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
            SLIC_DEBUG("\tclip: P=" << P << ", before=" << inputShape);
#endif
            remaining = axom::primal::clip(inputShape, P);
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
            SLIC_DEBUG("\tclip: after=" << clippedShape);
#endif
          }
          else
          {
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
            SLIC_DEBUG("\tclip: P=" << P << ", before=" << remaining);
#endif
            remaining = axom::primal::clip(remaining, P);
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
            SLIC_DEBUG("\tclip: after=" << remaining);
#endif
          }
        }

        // Emit the last leftover fragment.
        const auto fragmentIndex = offset + matCount - 1;
        const auto matId = sortedMaterialIdsView[fragmentIndex];
        const double *normalPtr = fragmentVectorsView.data() + (fragmentIndex * numVectorComponents);
        buildView.addShape(zoneIndex, fragmentIndex, remaining, matId, normalPtr);
      });
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif
