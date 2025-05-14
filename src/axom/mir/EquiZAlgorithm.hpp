// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EQUIZ_ALGORITHM_HPP_
#define AXOM_MIR_EQUIZ_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/utilities/ExtractZones.hpp"
#include "axom/mir/utilities/MergeMeshes.hpp"
#include "axom/mir/utilities/NodeToZoneRelationBuilder.hpp"
#include "axom/mir/utilities/RecenterField.hpp"
#include "axom/mir/utilities/ZoneListBuilder.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/MaterialView.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <algorithm>
#include <string>

// This macro makes the EquiZ algorithm split processing into clean/mixed stages.
#define AXOM_EQUIZ_SPLIT_PROCESSING

// Uncomment to save inputs and outputs.
// #define AXOM_EQUIZ_DEBUG

#if defined(AXOM_EQUIZ_DEBUG)
  #include <conduit/conduit_relay_io_blueprint.hpp>
#endif

namespace axom
{
namespace mir
{
using MaterialID = int;
using MaterialIDArray = axom::Array<MaterialID>;
using MaterialIDView = axom::ArrayView<MaterialID>;
using MaterialVF = float;
using MaterialVFArray = axom::Array<MaterialVF>;
using MaterialVFView = axom::ArrayView<MaterialVF>;

constexpr static int NULL_MATERIAL = -1;
constexpr static MaterialVF NULL_MATERIAL_VF = -1.f;

namespace internal
{
/*!
 * \brief This class is an intersection policy compatible with ClipField. It
 *        helps determine clip cases and weights using material-aware logic.
 *
 * \tparam ConnectivityT The type of index we'd see in the associated mesh's
 *                       connectivity. We template on it so we can pass the
 *                       array views of connectivity (node lists) to methods here.
 * \tparam MAXMATERIALS The max number of materials to handle.
 */
template <typename ConnectivityT, int MAXMATERIALS = 10>
class MaterialIntersector
{
public:
  using ConnectivityType = ConnectivityT;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /*!
   * \brief This is a view class for MatsetIntersector that can be used in device code.
   */
  struct View
  {
    static constexpr int INVALID_INDEX = -1;

    /*!
     * \brief Determine the clipping case, taking into account the zone's material
     *        and the current material being added.
     *
     * \param zoneIndex The zone index in the zoneMaterialView field.
     * \param nodeIds A view containing node ids for the current zone.
     *
     * \return The clip case number for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineClipCase(axom::IndexType zoneIndex,
                                      const ConnectivityView &nodeIdsView) const
    {
      // Determine the matvf view index for the material that owns the zone.
      int backgroundIndex = INVALID_INDEX;
      int zoneMatID = m_zoneMatNumberView[zoneIndex];
      if(zoneMatID != NULL_MATERIAL) backgroundIndex = matNumberToIndex(zoneMatID);

      axom::IndexType clipcase = 0;
      const auto n = nodeIdsView.size();
      for(IndexType i = 0; i < n; i++)
      {
        const auto nid = nodeIdsView[i];
        SLIC_ASSERT_MSG(
          nid >= 0 && nid < m_matvfViews[0].size(),
          axom::fmt::format("Node id {} is not in range [0, {}).", nid, m_matvfViews[0].size()));

        // clang-format off
        MaterialVF vf1 = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][nid] : NULL_MATERIAL_VF;
        MaterialVF vf2 = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][nid] : 0;
        // clang-format on

        clipcase |= (vf2 > vf1) ? (1 << i) : 0;
      }
      return clipcase;
    }

    /*!
     * \brief Compute the weight of a clip value along an edge (id0, id1) using the clip field and value.
     *
     * \param id0 The mesh node at the start of the edge.
     * \param id1 The mesh node at the end of the edge.
     */
    AXOM_HOST_DEVICE
    float computeWeight(axom::IndexType zoneIndex, ConnectivityType id0, ConnectivityType id1) const
    {
      // Determine the matvf view index for the material that owns the zone.
      int backgroundIndex = INVALID_INDEX;
      int zoneMatID = m_zoneMatNumberView[zoneIndex];
      if(zoneMatID != NULL_MATERIAL) backgroundIndex = matNumberToIndex(zoneMatID);
      // Determine the matvf view index for the current material.
      SLIC_ASSERT_MSG(
        id0 >= 0 && id0 < m_matvfViews[0].size(),
        axom::fmt::format("Node id {} is not in range [0, {}).", id0, m_matvfViews[0].size()));
      SLIC_ASSERT_MSG(
        id1 >= 0 && id1 < m_matvfViews[0].size(),
        axom::fmt::format("Node id {} is not in range [0, {}).", id1, m_matvfViews[0].size()));

      // Get the volume fractions for mat1, mat2 at the edge endpoints id0, id1.
      MaterialVF vf1[2], vf2[2];
      // clang-format off
      vf1[0] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id0] : NULL_MATERIAL_VF;
      vf1[1] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id1] : NULL_MATERIAL_VF;
      vf2[0] = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][id0] : 0;
      vf2[1] = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][id1] : 0;
      // clang-format on

      float numerator = vf2[0] - vf1[0];
      float denominator = -vf1[0] + vf1[1] + vf2[0] - vf2[1];

      float t = 0.f;
      if(denominator != 0.f)
      {
        t = numerator / denominator;
      }
      t = axom::utilities::clampVal(t, 0.f, 1.f);

      return t;
    }

    /*!
     * \brief Return the volume fraction array index in m_matIndicesView for the
     *        given material number \a matNumber.
     *
     * \param matNumber A material number that occurs in the matset material ids.
     *
     * \return The m_matNumbersView index on success; INVALID_INDEX on failure.
     */
    AXOM_HOST_DEVICE
    inline int matNumberToIndex(int matNumber) const
    {
      auto index = axom::mir::utilities::bsearch(matNumber, m_matNumbersView);
      return (index != -1) ? m_matIndicesView[index] : INVALID_INDEX;
    }

    /// Helper initialization methods for the host.

    void addMaterial(const MaterialVFView &matvf) { m_matvfViews.push_back(matvf); }

    void setMaterialNumbers(const axom::ArrayView<int> &matNumbersView)
    {
      m_matNumbersView = matNumbersView;
    }

    void setMaterialIndices(const axom::ArrayView<int> &matIndicesView)
    {
      m_matIndicesView = matIndicesView;
    }

    void setZoneMaterialID(const axom::ArrayView<int> &zoneMatsView)
    {
      m_zoneMatNumberView = zoneMatsView;
    }

    void setCurrentMaterial(int matNumber, int matNumberIndex)
    {
      m_currentMaterial = matNumber;
      m_currentMaterialIndex = matNumberIndex;
    }

    axom::StaticArray<MaterialVFView, MAXMATERIALS> m_matvfViews {};  //!< Array of volume fraction views
    axom::ArrayView<int> m_matNumbersView {};  //!< Sorted array of material numbers.
    axom::ArrayView<int> m_matIndicesView {};  //!< Array of indices into m_matvfViews for the material numbers.
    axom::ArrayView<int> m_zoneMatNumberView {};  //!< Contains the current material number that owns each zone.
    int m_currentMaterial {};                     //!< The current material.
    int m_currentMaterialIndex {};  //!< The current material's index in the m_matvfViews.
  };

  /*!
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const conduit::Node &AXOM_UNUSED_PARAM(n_options),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_fields))
  { }

  /*!
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The clipping options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &AXOM_UNUSED_PARAM(n_input),
                              const conduit::Node &n_options) const
  {
    return n_options["topology"].as_string();
  }

  /// Set various attributes.

  void addMaterial(const MaterialVFView &matvf) { m_view.addMaterial(matvf); }

  void setMaterialNumbers(const axom::ArrayView<int> &matNumbers)
  {
    m_view.setMaterialNumbers(matNumbers);
  }

  void setMaterialIndices(const axom::ArrayView<int> &matIndices)
  {
    m_view.setMaterialIndices(matIndices);
  }

  void setZoneMaterialID(const axom::ArrayView<int> &zoneMatsView)
  {
    m_view.setZoneMaterialID(zoneMatsView);
  }

  void setCurrentMaterial(int matNumber, int matNumberIndex)
  {
    m_view.setCurrentMaterial(matNumber, matNumberIndex);
  }

  /*!
   * \brief Return a new instance of the view.
   * \return A new instance of the view.
   * \note Call this after all values are set.
   */
  View view() const { return m_view; }

private:
  View m_view {};
};

}  // end namespace internal

/*!
 * \accelerated
 * \brief Implements Meredith's Equi-Z algorithm on the GPU using Blueprint inputs/outputs.
 *
 * \tparam ExecSpace the execution space where the algorithm will run.
 * \tparam TopologyView A topology view to be used for accessing zones in the mesh.
 * \tparam CoordsetView A coordset view that accesses coordinates as primal::Point.
 * \tparam MatsetView A matset view that interfaces to the Blueprint material set.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
class EquiZAlgorithm : public axom::mir::MIRAlgorithm
{
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

public:
  using ConnectivityType = typename TopologyView::ConnectivityType;

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  EquiZAlgorithm(const TopologyView &topoView,
                 const CoordsetView &coordsetView,
                 const MatsetView &matsetView)
    : axom::mir::MIRAlgorithm()
    , m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~EquiZAlgorithm() = default;

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
   *                      These are documented in the Sphinx documentation.
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
    AXOM_ANNOTATE_SCOPE("EquizAlgorithm");

    // Copy the options.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

#if defined(AXOM_EQUIZ_DEBUG)
    // Save the MIR input.
    conduit::Node n_tmpInput;
    n_tmpInput[n_topo.path()].set_external(n_topo);
    n_tmpInput[n_coordset.path()].set_external(n_coordset);
    n_tmpInput[n_fields.path()].set_external(n_fields);
    n_tmpInput[n_matset.path()].set_external(n_matset);
    saveMesh(n_tmpInput, "debug_equiz_input");
#endif

#if defined(AXOM_EQUIZ_SPLIT_PROCESSING)
    // Come up with lists of clean/mixed zones.
    axom::Array<axom::IndexType> cleanZones, mixedZones;
    makeZoneLists(n_options_copy, cleanZones, mixedZones);
    SLIC_ASSERT((cleanZones.size() + mixedZones.size()) == m_topologyView.numberOfZones());
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
      conduit::Node &n_root_fields = n_root["fields"];
      for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
      {
        n_root_fields[n_fields[i].name()].set_external(n_fields[i]);
      }

      // Make the clean mesh.
      conduit::Node n_cleanOutput;
      makeCleanOutput(n_root, n_topo.name(), n_options_copy, cleanZones.view(), n_cleanOutput);

      // Add an original nodes field on the root mesh.
      addOriginal(n_root_fields[originalNodesFieldName()],
                  n_topo.name(),
                  "vertex",
                  m_coordsetView.numberOfNodes());
      // If there are fields in the options, make sure the new field is handled too.
      if(n_options_copy.has_child("fields"))
      {
        n_options_copy["fields/" + originalNodesFieldName()] = originalNodesFieldName();
      }

      // Process the mixed part of the mesh. We select just the mixed zones.
      n_options_copy["selectedZones"].set_external(mixedZones.data(), mixedZones.size());
      n_options_copy["newNodesField"] = newNodesFieldName();
      processMixedZones(n_root_topo,
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
  #if defined(AXOM_EQUIZ_DEBUG)
      saveMesh(n_mirOutput, "debug_equiz_mir");
      std::cout << "--- clean ---\n";
      printNode(n_cleanOutput);
      std::cout << "--- MIR ---\n";
      printNode(n_mirOutput);
  #endif

      // Merge the clean zones and MIR output
      conduit::Node n_merged;
      merge(n_newTopo.name(), n_cleanOutput, n_mirOutput, n_merged);
  #if defined(AXOM_EQUIZ_DEBUG)
      std::cout << "--- merged ---\n";
      printNode(n_merged);

      // Save merged output.
      saveMesh(n_merged, "debug_equiz_merged");
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
      processMixedZones(n_topo,
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
#else
    // Handle all zones via MIR.
    processMixedZones(n_topo,
                      n_coordset,
                      n_fields,
                      n_matset,
                      n_options_copy,
                      n_newTopo,
                      n_newCoordset,
                      n_newFields,
                      n_newMatset);
#endif
  }

#if defined(AXOM_EQUIZ_SPLIT_PROCESSING)
  /*!
   * \brief Examine the zones and make separate lists of clean and mixed zones
   *        so we can process them specially since mixed zones require much
   *        more work.
   */
  void makeZoneLists(const conduit::Node &n_options,
                     axom::Array<axom::IndexType> &cleanZones,
                     axom::Array<axom::IndexType> &mixedZones) const
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Call variants of the ZoneListBuilder methods that take into account adjacent
    // zones materials when determining if a zone should be mixed.

    // _mir_utilities_zlb_begin
    bputils::ZoneListBuilder<ExecSpace, TopologyView, MatsetView> zlb(m_topologyView, m_matsetView);
    if(n_options.has_child("selectedZones"))
    {
      auto selectedZonesView =
        bputils::make_array_view<axom::IndexType>(n_options.fetch_existing("selectedZones"));
      zlb.execute(m_coordsetView.numberOfNodes(), selectedZonesView, cleanZones, mixedZones);
    }
    else
    {
      zlb.execute(m_coordsetView.numberOfNodes(), cleanZones, mixedZones);
    }
    // _mir_utilities_zlb_end
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

    // Make node map and slice info for merging.
    axom::Array<axom::IndexType> nodeMap, nodeSlice;
    conduit::Node &n_mir_fields = n_mirOutput["fields"];
    createNodeMapAndSlice(n_mir_fields, nodeMap, nodeSlice);

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

    inputs[1].m_input = &n_mirOutput;
    inputs[1].m_nodeMapView = nodeMap.view();
    inputs[1].m_nodeSliceView = nodeSlice.view();

    conduit::Node mmOpts;
    mmOpts["topology"] = topoName;
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
   * \note This field is added to the mesh before feeding it through MIR so we will have an idea
   *       of which nodes are original nodes in the output. Blended nodes may not have good values
   *       but there is a mask field that can identify those nodes.
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
   * \param n_options Options to forward.
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

    // Make the clean mesh. Set compact=0 so it does not change the number of nodes.
    bputils::ExtractZonesAndMatset<ExecSpace, TopologyView, CoordsetView, MatsetView> ez(
      m_topologyView,
      m_coordsetView,
      m_matsetView);
    conduit::Node n_ezopts;
    n_ezopts["topology"] = topoName;
    n_ezopts["compact"] = 0;
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
  #if defined(AXOM_EQUIZ_DEBUG)
    AXOM_ANNOTATE_BEGIN("saveClean");
    saveMesh(n_cleanOutput, "debug_equiz_clean");
    AXOM_ANNOTATE_END("saveClean");
  #endif
  }

  /*!
   * \brief Create node map and node slice arrays for the MIR output that help
   *        merge it back with the clean output.
   *
   * \param n_newFields The fields from the MIR output.
   * \param[out] nodeMap An array used to map node ids from the MIR output to their node ids in the merged mesh.
   * \param[out] nodeSlice An array that identifies new blended node ids in the MIR output so they can be appended into coordsets and fields during merge.
   */
  void createNodeMapAndSlice(conduit::Node &n_newFields,
                             axom::Array<axom::IndexType> &nodeMap,
                             axom::Array<axom::IndexType> &nodeSlice) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("createNodeMapAndSlice");
    SLIC_ASSERT(n_newFields.has_child(originalNodesFieldName()));
    SLIC_ASSERT(n_newFields.has_child(newNodesFieldName()));

    const axom::IndexType numCleanNodes = m_coordsetView.numberOfNodes();

    // These are the original node ids.
    const conduit::Node &n_output_orig_nodes = n_newFields[originalNodesFieldName() + "/values"];
    auto numOutputNodes = n_output_orig_nodes.dtype().number_of_elements();
    auto outputOrigNodesView = bputils::make_array_view<ConnectivityType>(n_output_orig_nodes);

    // __equiz_new_nodes is the int mask field that identifies new nodes created from blending.
    const conduit::Node &n_new_nodes_values = n_newFields[newNodesFieldName() + "/values"];
    const auto maskView = bputils::make_array_view<int>(n_new_nodes_values);

    // Count new nodes created from blending.
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<int> maskOffset(numOutputNodes, numOutputNodes, allocatorID);
    auto maskOffsetsView = maskOffset.view();
    RAJA::ReduceSum<reduce_policy, int> mask_reduce(0);
    axom::for_all<ExecSpace>(
      numOutputNodes,
      AXOM_LAMBDA(axom::IndexType index) { mask_reduce += maskView[index]; });
    const auto numNewNodes = mask_reduce.get();

    // Make offsets.
    axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

    // Make a list of indices that we need to slice out of the node arrays.
    nodeSlice = axom::Array<axom::IndexType>(numNewNodes, numNewNodes, allocatorID);
    auto nodeSliceView = nodeSlice.view();
    axom::for_all<ExecSpace>(
      numOutputNodes,
      AXOM_LAMBDA(axom::IndexType index) {
        if(maskView[index] > 0)
        {
          nodeSliceView[maskOffsetsView[index]] = index;
        }
      });

    // Make a node map for mapping mixed connectivity into combined node numbering.
    nodeMap = axom::Array<axom::IndexType>(numOutputNodes, numOutputNodes, allocatorID);
    auto nodeMapView = nodeMap.view();
    axom::for_all<ExecSpace>(
      numOutputNodes,
      AXOM_LAMBDA(axom::IndexType index) {
        if(maskView[index] == 0)
        {
          nodeMapView[index] = static_cast<axom::IndexType>(outputOrigNodesView[index]);
        }
        else
        {
          nodeMapView[index] = numCleanNodes + maskOffsetsView[index];
        }
      });

    // Remove fields that are no longer needed.
    n_newFields.remove(originalNodesFieldName());
    n_newFields.remove(newNodesFieldName());
  }
#endif

  /*!
   * \brief Perform material interface reconstruction on mixed zones.
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
  void processMixedZones(const conduit::Node &n_topo,
                         const conduit::Node &n_coordset,
                         const conduit::Node &n_fields,
                         const conduit::Node &n_matset,
                         conduit::Node &n_options,
                         conduit::Node &n_newTopo,
                         conduit::Node &n_newCoordset,
                         conduit::Node &n_newFields,
                         conduit::Node &n_newMatset) const
  {
    AXOM_ANNOTATE_SCOPE("processMixedZones");

    // Make some nodes that will contain the inputs to subsequent iterations.
    // Store them under a single node so the nodes will have names.
    conduit::Node n_Input;
    conduit::Node &n_InputTopo = n_Input[n_topo.path()];
    conduit::Node &n_InputCoordset = n_Input[n_coordset.path()];
    conduit::Node &n_InputFields = n_Input[n_fields.path()];

    // Get the materials from the matset and determine which of them are clean/mixed.
    axom::mir::views::MaterialInformation allMats, cleanMats, mixedMats;
    classifyMaterials(n_matset, allMats, cleanMats, mixedMats);

    //--------------------------------------------------------------------------
    //
    // Make node-centered VF fields and add various working fields.
    //
    //--------------------------------------------------------------------------
    n_InputFields.reset();
    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      n_InputFields[n_fields[i].name()].set_external(n_fields[i]);
    }
    makeNodeCenteredVFs(n_topo, n_coordset, n_InputFields, mixedMats);
    makeWorkingFields(n_topo, n_InputFields, cleanMats, mixedMats);

    //--------------------------------------------------------------------------
    //
    // Iterate over mixed materials.
    //
    //--------------------------------------------------------------------------
    constexpr int first = 0;
    for(size_t i = first; i < mixedMats.size(); i++)
    {
      if(i == first)
      {
        // The first time through, we can use the supplied views.
        iteration(i,
                  m_topologyView,
                  m_coordsetView,

                  allMats,
                  mixedMats[i],

                  n_topo,
                  n_coordset,
                  n_InputFields,

                  n_options,

                  n_newTopo,
                  n_newCoordset,
                  n_newFields);

        // In later iterations, we do not want to pass selectedZones through
        // since they are only valid on the current input topology. Also, if they
        // were passed then the new topology only has those selected zones.
        if(n_options.has_child("selectedZones"))
        {
          n_options.remove("selectedZones");
        }
      }
      else
      {
        // Clear the inputs from the last iteration.
        n_InputTopo.reset();
        n_InputCoordset.reset();
        n_InputFields.reset();

        // Move the outputs of the last iteration to the inputs of this iteration.
        n_InputTopo.move(n_newTopo);
        n_InputCoordset.move(n_newCoordset);
        n_InputFields.move(n_newFields);

        // Create an appropriate coordset view.
        using CSDataType = typename CoordsetView::value_type;
        auto coordsetView =
          axom::mir::views::make_explicit_coordset<CSDataType, CoordsetView::dimension()>::view(
            n_InputCoordset);

        using ConnectivityType = typename TopologyView::ConnectivityType;
        // Dispatch to an appropriate topo view, taking into account the connectivity
        // type and the possible shapes that would be supported for the input topology.
        views::typed_dispatch_unstructured_topology<ConnectivityType,
                                                    views::view_traits<TopologyView>::selected_shapes()>(
          n_InputTopo,
          [&](const auto &AXOM_UNUSED_PARAM(shape), auto topologyView) {
            // Do the next iteration (uses new topologyView type).
            iteration(i,
                      topologyView,
                      coordsetView,

                      allMats,
                      mixedMats[i],

                      n_InputTopo,
                      n_InputCoordset,
                      n_InputFields,

                      n_options,

                      n_newTopo,
                      n_newCoordset,
                      n_newFields);
          });
      }
    }

    // Build the new matset.
    buildNewMatset(n_matset, n_newFields, n_newMatset);

    // Cleanup.
    {
      AXOM_ANNOTATE_SCOPE("cleanup");
      for(const auto &mat : allMats)
      {
        const std::string nodalMatName(nodalFieldName(mat.number));
        if(n_newFields.has_child(nodalMatName))
        {
          n_newFields.remove(nodalMatName);
        }
#if defined(AXOM_EQUIZ_DEBUG)
        const std::string zonalMatName(zonalFieldName(mat.number));
        if(n_newFields.has_child(zonalMatName))
        {
          n_newFields.remove(zonalMatName);
        }
#endif
      }
      n_newFields.remove(zonalMaterialIDName());
    }

#if defined(AXOM_EQUIZ_DEBUG)
    //--------------------------------------------------------------------------
    //
    // Save the MIR output.
    //
    //--------------------------------------------------------------------------
    conduit::Node n_output;
    n_output[n_newTopo.path()].set_external(n_newTopo);
    n_output[n_newCoordset.path()].set_external(n_newCoordset);
    n_output[n_newFields.path()].set_external(n_newFields);
    n_output[n_newMatset.path()].set_external(n_newMatset);
    saveMesh(n_output, "debug_equiz_output");
#endif
  }

  /*!
   * \brief Examine the materials and determine which are clean/mixed.
   *
   * \param n_matset A Conduit node containing the matset.
   * \param[out] allMats A vector of all of the materials.
   * \param[out] cleanMats A vector of the clean materials.
   * \param[out] mixedMats A vector of the mixed materials.
   */
  void classifyMaterials(const conduit::Node &n_matset,
                         axom::mir::views::MaterialInformation &allMats,
                         axom::mir::views::MaterialInformation &cleanMats,
                         axom::mir::views::MaterialInformation &mixedMats) const
  {
    AXOM_ANNOTATE_SCOPE("classifyMaterials");

    cleanMats.clear();
    mixedMats.clear();
    allMats = axom::mir::views::materials(n_matset);

    // TODO: actually determine which materials are clean/mixed. It's probably
    //       best to ask the matsetView since it takes some work to determine
    //       this.

    mixedMats = allMats;
  }

  /*!
   * \brief Return the name of the zonal material field for a given matId.
   * \return The name of the zonal material field.
   */
  std::string zonalFieldName(int matId) const
  {
    std::stringstream ss;
    ss << "__equiz_zonal_volume_fraction_" << matId;
    return ss.str();
  }

  /*!
   * \brief Return the name of the nodal material field for a given matId.
   * \return The name of the nodal material field.
   */
  std::string nodalFieldName(int matId) const
  {
    std::stringstream ss;
    ss << "__equiz_nodal_volume_fraction_" << matId;
    return ss.str();
  }

  /*!
   * \brief Return the name of the zonal material id field.
   * \return The name of the zonal material id field.
   */
  std::string zonalMaterialIDName() const { return "__equiz_zonalMaterialID"; }

  /*!
   * \brief Return the name of the original nodes field.
   * \return The name of the original nodes field.
   */
  std::string originalNodesFieldName() const { return "__equiz_original_node"; }

  /*!
   * \brief Return the name of the new nodes field that identifies blended nodes in the MIR output.
   * \return The name of the new nodes field.
   */
  std::string newNodesFieldName() const { return "__equiz_new_nodes"; }

  /*!
   * \brief Makes node-cenetered volume fractions for the materials in the matset
   *        and attaches them as fields.
   *
   * \param n_topo A Conduit node containing the input topology.
   * \param n_coordset A Conduit node containin the input coordset.
   * \param[inout] A Conduit node where the new fields will be added.
   * \param mixedMats A vector of mixed materials.
   */
  void makeNodeCenteredVFs(const conduit::Node &n_topo,
                           const conduit::Node &n_coordset,
                           conduit::Node &n_fields,
                           const axom::mir::views::MaterialInformation &mixedMats) const
  {
    AXOM_ANNOTATE_SCOPE("makeNodeCenteredVFs");

    namespace bputils = axom::mir::utilities::blueprint;
    // Make a node to zone relation so we know for each node, which zones it touches.
    conduit::Node relation;
    {
      AXOM_ANNOTATE_SCOPE("relation");
      bputils::NodeToZoneRelationBuilder<ExecSpace> rb;
      rb.execute(n_topo, n_coordset, relation);
      //printNode(relation);
      //std::cout.flush();
    }

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Make nodal VFs for each mixed material.
    const auto nzones = m_topologyView.numberOfZones();
    const auto nnodes = m_coordsetView.numberOfNodes();
    {
      AXOM_ANNOTATE_SCOPE("zonal");
      for(const auto &mat : mixedMats)
      {
        const int matNumber = mat.number;
        const std::string zonalName = zonalFieldName(matNumber);
        conduit::Node &n_zonalField = n_fields[zonalName];
        n_zonalField["topology"] = n_topo.name();
        n_zonalField["association"] = "element";
        n_zonalField["values"].set_allocator(c2a.getConduitAllocatorID());
        n_zonalField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nzones));
        auto zonalFieldView = bputils::make_array_view<MaterialVF>(n_zonalField["values"]);

        // Fill the zonal field from the matset.
        MatsetView deviceMatsetView(m_matsetView);
        axom::for_all<ExecSpace>(
          m_topologyView.numberOfZones(),
          AXOM_LAMBDA(axom::IndexType zoneIndex) {
            typename MatsetView::FloatType vf {};
            deviceMatsetView.zoneContainsMaterial(zoneIndex, matNumber, vf);
            zonalFieldView[zoneIndex] = static_cast<MaterialVF>(vf);
          });
      }
    }

    {
      AXOM_ANNOTATE_SCOPE("recenter");
      for(const auto &mat : mixedMats)
      {
        const int matNumber = mat.number;
        const std::string zonalName = zonalFieldName(matNumber);
        conduit::Node &n_zonalField = n_fields[zonalName];

        // Make a nodal field for the current material by recentering.
        const std::string nodalName = nodalFieldName(matNumber);
        conduit::Node &n_nodalField = n_fields[nodalName];
        n_nodalField["topology"] = n_topo.name();
        n_nodalField["association"] = "vertex";
        n_nodalField["values"].set_allocator(c2a.getConduitAllocatorID());
        n_nodalField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, nnodes));
        bputils::RecenterField<ExecSpace> z2n;
        z2n.execute(n_zonalField, relation, n_nodalField);

#if !defined(AXOM_EQUIZ_DEBUG)
        // Remove the zonal field that we don't normally need (unless we're debugging).
        n_fields.remove(zonalName);
#endif
      }
    }
  }

  /*!
   * \brief Set up the "working fields", mainly a zonalMaterialID that includes
   *        the contributions from the clean materials and the first mixed material.
   *
   * \param n_topo A Conduit node containing the input topology pre-MIR.
   * \param n_fields A Conduit node containing the fields pre-MIR.
   * \param cleanMats A vector of clean materials.
   * \param mixedMats A vector of mixed materials.
   */
  void makeWorkingFields(const conduit::Node &n_topo,
                         conduit::Node &n_fields,
                         const axom::mir::views::MaterialInformation &cleanMats,
                         const axom::mir::views::MaterialInformation &AXOM_UNUSED_PARAM(mixedMats)) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("makeWorkingFields");

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const auto nzones = m_topologyView.numberOfZones();

    // Make the zonal id field.
    conduit::Node &n_zonalIDField = n_fields[zonalMaterialIDName()];
    n_zonalIDField["topology"] = n_topo.name();
    n_zonalIDField["association"] = "element";
    n_zonalIDField["values"].set_allocator(c2a.getConduitAllocatorID());
    n_zonalIDField["values"].set(conduit::DataType(bputils::cpp2conduit<MaterialID>::id, nzones));
    auto zonalIDFieldView = bputils::make_array_view<MaterialID>(n_zonalIDField["values"]);

    // Fill all zones with NULL_MATERIAL.
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType nodeIndex) { zonalIDFieldView[nodeIndex] = NULL_MATERIAL; });

    // Fill in the clean zones.
    using FloatType = typename MatsetView::FloatType;
    MatsetView deviceMatsetView(m_matsetView);
    for(const auto &mat : cleanMats)
    {
      const int matNumber = mat.number;
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) {
          FloatType vf {};
          if(deviceMatsetView.zoneContainsMaterial(zoneIndex, matNumber, vf))
          {
            zonalIDFieldView[zoneIndex] = matNumber;
          }
        });
    }
  }

  /*!
   * \brief Perform one iteration of material clipping.
   *
   * \tparam ITopologyView The topology view type for the intermediate topology.
   * \tparam ICoordsetView The topology view type for the intermediate coordset.
   *
   * \param iter The iteration number.
   * \param topoView The topology view for the intermediate input topology.
   * \param coordsetView The coordset view for the intermediate input coordset.
   * \param allMats A vector of Material information (all materials).
   * \param currentMat A Material object for the current material.
   * \param n_topo A Conduit node containing the intermediate input topology.
   * \param n_fields A Conduit node containing the intermediate input fields.
   * \param n_options MIR options.
   * \param n_newTopo[out] A Conduit node to contain the new topology.
   * \param n_newCoordset[out] A Conduit node to contain the new coordset.
   * \param n_newFields[out] A Conduit node to contain the new fields.
   *
   * \note This algorithm uses a ClipField with a MaterialIntersector that gives
   *       it the ability to access nodal volume fraction fields and make intersection
   *       decisions with that data.
   */
  template <typename ITopologyView, typename ICoordsetView>
  void iteration(int iter,
                 const ITopologyView &topoView,
                 const ICoordsetView &coordsetView,

                 const axom::mir::views::MaterialInformation &allMats,
                 const axom::mir::views::Material &currentMat,

                 const conduit::Node &n_topo,
                 const conduit::Node &n_coordset,
                 conduit::Node &n_fields,

                 const conduit::Node &n_options,

                 conduit::Node &n_newTopo,
                 conduit::Node &n_newCoordset,
                 conduit::Node &n_newFields) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    namespace bpmeshutils = conduit::blueprint::mesh::utils;
    AXOM_ANNOTATE_SCOPE(axom::fmt::format("iteration {}", iter));

    const std::string colorField("__equiz__colors");

#if defined(AXOM_EQUIZ_DEBUG)
    //--------------------------------------------------------------------------
    //
    // Save the iteration inputs.
    //
    //--------------------------------------------------------------------------
    {
      AXOM_ANNOTATE_SCOPE("Saving input");
      conduit::Node n_mesh_input;
      n_mesh_input[n_topo.path()].set_external(n_topo);
      n_mesh_input[n_coordset.path()].set_external(n_coordset);
      n_mesh_input[n_fields.path()].set_external(n_fields);

      // save
      std::stringstream ss1;
      ss1 << "debug_equiz_input_iter." << iter;
      saveMesh(n_mesh_input, ss1.str());
    }
#else
    AXOM_UNUSED_VAR(iter);
#endif

    //--------------------------------------------------------------------------
    //
    // Make material intersector.
    //
    //--------------------------------------------------------------------------
    using ConnectivityType = typename ITopologyView::ConnectivityType;
    using IntersectorType = internal::MaterialIntersector<ConnectivityType, MatsetView::MaxMaterials>;

    IntersectorType intersector;
    int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const int nmats = static_cast<int>(allMats.size());
    axom::Array<int> matNumberDevice(nmats, nmats, allocatorID),
      matIndexDevice(nmats, nmats, allocatorID);
    {
      AXOM_ANNOTATE_SCOPE("Intersector setup");
      // Populate intersector, including making a number:index map
      axom::Array<int> matNumber, matIndex;
      for(int index = 0; index < nmats; index++)
      {
        // Add a matvf view to the intersector.
        const std::string matFieldName = nodalFieldName(allMats[index].number);
        auto matVFView =
          bputils::make_array_view<MaterialVF>(n_fields.fetch_existing(matFieldName + "/values"));
        intersector.addMaterial(matVFView);

        matNumber.push_back(allMats[index].number);
        matIndex.push_back(index);
      }
      // Sort indices by matNumber.
      std::sort(matIndex.begin(), matIndex.end(), [&](auto idx1, auto idx2) {
        return matNumber[idx1] < matNumber[idx2];
      });
      std::sort(matNumber.begin(), matNumber.end());
      // Get the current material's index in the number:index map.
      int currentMatIndex = 0;
      for(axom::IndexType i = 0; i < matNumber.size(); i++)
      {
        if(matNumber[i] == currentMat.number)
        {
          currentMatIndex = matIndex[i];
          break;
        }
      }

      // Store the number:index map into the intersector. The number:index map lets us
      // ask for the field index for a material number, allowing scattered material
      // numbers to be used in the matset.
      axom::copy(matNumberDevice.data(), matNumber.data(), sizeof(int) * nmats);
      axom::copy(matIndexDevice.data(), matIndex.data(), sizeof(int) * nmats);
      intersector.setMaterialNumbers(matNumberDevice.view());
      intersector.setMaterialIndices(matIndexDevice.view());
      intersector.setCurrentMaterial(currentMat.number, currentMatIndex);

      // Store the current zone material ids and current material number into the intersector.
      intersector.setZoneMaterialID(bputils::make_array_view<MaterialID>(
        n_fields.fetch_existing(zonalMaterialIDName() + "/values")));
    }

    //--------------------------------------------------------------------------
    //
    // Make clip options
    //
    //--------------------------------------------------------------------------
    conduit::Node options;
    options["inside"] = 1;
    options["outside"] = 1;
    options["colorField"] = colorField;
    if(n_options.has_child("selectedZones"))
    {
      // Pass selectedZones along in the clip options, if present.
      options["selectedZones"].set_external(n_options.fetch_existing("selectedZones"));
    }
    if(n_options.has_child("fields"))
    {
      // Pass along fields, if present.
      options["fields"].set_external(n_options.fetch_existing("fields"));
    }
    if(n_options.has_child("newNodesField"))
    {
      // Pass along newNodesField, if present.
      options["newNodesField"] = n_options.fetch_existing("newNodesField").as_string();
    }
    options["topology"] = n_options["topology"];

    //--------------------------------------------------------------------------
    //
    // Clip the topology using the material intersector.
    //
    //--------------------------------------------------------------------------
    {
      using ClipperType =
        axom::mir::clipping::ClipField<ExecSpace, ITopologyView, ICoordsetView, IntersectorType>;
      ClipperType clipper(topoView, coordsetView, intersector);
      clipper.execute(n_topo, n_coordset, n_fields, options, n_newTopo, n_newCoordset, n_newFields);
    }

    //--------------------------------------------------------------------------
    //
    // Update zoneMaterialID based on color field.
    //
    //--------------------------------------------------------------------------
    {
      AXOM_ANNOTATE_SCOPE("Update zonalMaterialID");

      const auto colorView =
        bputils::make_array_view<int>(n_newFields.fetch_existing(colorField + "/values"));
      const auto nzonesNew = colorView.size();

      // Get zonalMaterialID field so we can make adjustments.
      conduit::Node &n_zonalMaterialID =
        n_newFields.fetch_existing(zonalMaterialIDName() + "/values");
      auto zonalMaterialID = bputils::make_array_view<MaterialID>(n_zonalMaterialID);
      const int currentMatNumber = currentMat.number;
      axom::for_all<ExecSpace>(
        nzonesNew,
        AXOM_LAMBDA(axom::IndexType zoneIndex) {
          // Color the part we want with the current material.
          if(colorView[zoneIndex] == 1)
          {
            zonalMaterialID[zoneIndex] = currentMatNumber;
          }
        });
    }

#if defined(AXOM_EQUIZ_DEBUG)
    //--------------------------------------------------------------------------
    //
    // Save the clip results.
    //
    //--------------------------------------------------------------------------
    {
      AXOM_ANNOTATE_SCOPE("Saving output");
      conduit::Node mesh;
      mesh[n_newTopo.path()].set_external(n_newTopo);
      mesh[n_newCoordset.path()].set_external(n_newCoordset);
      mesh[n_newFields.path()].set_external(n_newFields);

      // save
      std::stringstream ss;
      ss << "debug_equiz_output_iter." << iter;
      saveMesh(mesh, ss.str());
    }
#endif

    // We do not want the color field to survive into the next iteration.
    n_newFields.remove(colorField);
  }

  /*!
   * \brief Build a new matset with only clean zones, representing the MIR output.
   *
   * \param n_matset n_matset The Conduit node that contains the input matset.
   * \param[inout] n_newFields The Conduit node that contains the fields for the MIR output.
   * \param[out] n_newMatset The node that contains the new matset.
   */
  void buildNewMatset(const conduit::Node &n_matset,
                      conduit::Node &n_newFields,
                      conduit::Node &n_newMatset) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("buildNewMatset");

    // Get the zonalMaterialID field that has our new material ids.
    conduit::Node &n_zonalMaterialID = n_newFields[zonalMaterialIDName() + "/values"];
    auto zonalMaterialID = bputils::make_array_view<MaterialID>(n_zonalMaterialID);
    const auto nzones = n_zonalMaterialID.dtype().number_of_elements();

    // Copy some information from the old matset to the new one.
    if(n_matset.has_child("topology"))
    {
      n_newMatset["topology"].set(n_matset.fetch_existing("topology"));
    }
    if(n_matset.has_child("material_map"))
    {
      n_newMatset["material_map"].set(n_matset.fetch_existing("material_map"));
    }

    // Make new nodes in the matset.
    conduit::Node &n_material_ids = n_newMatset["material_ids"];
    conduit::Node &n_volume_fractions = n_newMatset["volume_fractions"];
    conduit::Node &n_sizes = n_newMatset["sizes"];
    conduit::Node &n_offsets = n_newMatset["offsets"];
    conduit::Node &n_indices = n_newMatset["indices"];

    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set_allocator(c2a.getConduitAllocatorID());

    // We'll store the output matset in the same types as the input matset.
    using MIntType = typename MatsetView::IndexType;
    using MFloatType = typename MatsetView::FloatType;
    n_material_ids.set(conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));
    n_volume_fractions.set(conduit::DataType(bputils::cpp2conduit<MFloatType>::id, nzones));
    n_sizes.set(conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));
    n_offsets.set(conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));
    n_indices.set(conduit::DataType(bputils::cpp2conduit<MIntType>::id, nzones));

    auto material_ids_view = bputils::make_array_view<MIntType>(n_material_ids);
    auto volume_fractions_view = bputils::make_array_view<MFloatType>(n_volume_fractions);
    auto sizes_view = bputils::make_array_view<MIntType>(n_sizes);
    auto offsets_view = bputils::make_array_view<MIntType>(n_offsets);
    auto indices_view = bputils::make_array_view<MIntType>(n_indices);

    // Fill in the new matset data arrays.
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        material_ids_view[zoneIndex] = static_cast<MIntType>(zonalMaterialID[zoneIndex]);
        volume_fractions_view[zoneIndex] = 1;
        sizes_view[zoneIndex] = 1;
        offsets_view[zoneIndex] = static_cast<MIntType>(zoneIndex);
        indices_view[zoneIndex] = static_cast<MIntType>(zoneIndex);
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
