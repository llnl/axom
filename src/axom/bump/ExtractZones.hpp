// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_EXTRACT_ZONES_HPP
#define AXOM_BUMP_EXTRACT_ZONES_HPP

#include "axom/core.hpp"
#include "axom/bump/CoordsetBlender.hpp"
#include "axom/bump/CoordsetSlicer.hpp"
#include "axom/bump/FieldSlicer.hpp"
#include "axom/bump/MatsetSlicer.hpp"
#include "axom/bump/Options.hpp"
#include "axom/bump/views/view_traits.hpp"

namespace axom
{
namespace bump
{
/*!
 * \brief Make a new topology and coordset by extracting certain zones from the input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The topology view type on which the algorithm will run.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class ExtractZones
{
  using ConnectivityType = typename TopologyView::ConnectivityType;

public:
  using SelectedZonesView = axom::ArrayView<axom::IndexType>;
  using ZoneType = typename TopologyView::ShapeType;

  /*!
   * \brief Constructor
   *
   * \param topoView The input topology view.
   * \param coordsetView The input coordset view.
   * \param matsetView The input matset view.
   */
  ExtractZones(const TopologyView &topoView, const CoordsetView &coordsetView)
    : m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_zoneSlice()
  { }

  /*!
   * \brief Destructor
   *
   */
  virtual ~ExtractZones() = default;

  /*!
   * \brief Select zones from the input mesh by id and output them in the output mesh.
   *
   * \param selectedZonesView A view that contains the selected zone ids.
   * \param n_input The input mesh.
   * \param n_options The input options.
   * \param[out] n_output The output mesh.
   *
   * The \a n_options node controls how the algorithm works.
   *
   * \code{.yaml}
   * topology: mesh
   * compact: 1
   * extra:
   *  nodes: 0
   *  zones: 0
   *  connectivity: 0
   * \endcode
   *
   * The "topology" node contains a string that selects the name of the topology
   * to extract. The "compact" node causes the algorithm to restrict the output
   * coordset and vertex fields to only the nodes used by the selected zones. If
   * "compact" is set to 0 then the original coordset and vertex fields will retain
   * their size in the output mesh. The "extra" node is optional and it contains
   * 3 integer values for extra allocation to be made for nodes, zones, and connectivity.
   * This extra space can be filled in later by the application.
   */
  void execute(const SelectedZonesView &selectedZonesView,
               const conduit::Node &n_input,
               const conduit::Node &n_options,
               conduit::Node &n_output)
  {
    AXOM_ANNOTATE_SCOPE("ExtractZones");
    namespace utils = axom::bump::utilities;

    // Determine the dataSizes and map/slice information for nodes.
    axom::Array<ConnectivityType> old2new;
    axom::Array<axom::IndexType> nodeSlice;
    const auto extra = getExtra(n_options);
    Sizes dataSizes {};
    if(compact(n_options))
    {
      dataSizes = compactNodeMap(selectedZonesView, extra, old2new, nodeSlice);
    }
    else
    {
      dataSizes = nodeMap(selectedZonesView, extra, old2new, nodeSlice);
    }

    Options opts(n_options);

    // Make a new output topology.
    const conduit::Node &n_topologies = n_input.fetch_existing("topologies");
    const std::string topoName = topologyName(n_input, n_options);
    const conduit::Node &n_topo = n_topologies.fetch_existing(topoName);
    const std::string newTopoName = opts.topologyName(topoName);
    conduit::Node &n_newTopo = n_output["topologies/" + newTopoName];
    makeTopology(selectedZonesView, dataSizes, extra, old2new.view(), n_topo, n_options, n_newTopo);

    // Make a new coordset.
    SliceData nSlice;
    nSlice.m_indicesView = nodeSlice.view();
    const std::string coordsetName = n_topo.fetch_existing("coordset").as_string();
    const conduit::Node &n_coordset = n_input.fetch_existing("coordsets/" + coordsetName);
    const std::string newCoordsetName = opts.coordsetName(coordsetName);
    conduit::Node &n_newCoordset = n_output["coordsets/" + newCoordsetName];
    makeCoordset(nSlice, n_coordset, n_newCoordset);

    // Update the coordset name in the topo.
    n_newTopo["coordset"] = newCoordsetName;

    // Make new fields. If originalZones are present, they'll be sliced there.
    bool makeOriginalZones = true;
    if(n_input.has_child("fields"))
    {
      const conduit::Node &n_fields = n_input.fetch_existing("fields");
      conduit::Node &n_newFields = n_output["fields"];
      SliceData zSlice;
      zSlice.m_indicesView = zoneSliceView(selectedZonesView, extra);
      makeOriginalZones = !n_fields.has_child(opts.originalElementsField());
      makeFields(nSlice, zSlice, newTopoName, n_fields, n_newFields);
    }

    // Make originalElements.
    if(makeOriginalZones)
    {
      utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

      conduit::Node &n_outFields = n_output["fields"];
      conduit::Node &n_origElements = n_outFields[opts.originalElementsField()];
      n_origElements["topology"] = newTopoName;
      n_origElements["association"] = "element";
      n_origElements["values"].set_allocator(c2a.getConduitAllocatorID());
      n_origElements["values"].set(
        conduit::DataType(utils::cpp2conduit<axom::IndexType>::id, selectedZonesView.size()));
      axom::copy(n_origElements["values"].data_ptr(),
                 selectedZonesView.data(),
                 sizeof(axom::IndexType) * selectedZonesView.size());
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif
  /*!
   * \brief This struct contains extra amounts of storage that we might want to overallocate.
   */
  struct Sizes
  {
    axom::IndexType nodes {0};
    axom::IndexType zones {0};
    axom::IndexType connectivity {0};
  };

  /*!
   * \brief Create a zone slice view, building m_zoneSlice if needed.
   *
   * \param selectedZonesView A view that contains the selected zone ids.
   * \param extra A Sizes object containing any extra size that needs to be allocated.
   *
   * \return An array view containing the zone slice.
   */
  axom::ArrayView<axom::IndexType> zoneSliceView(const SelectedZonesView &selectedZonesView,
                                                 const Sizes &extra)
  {
    axom::ArrayView<axom::IndexType> view;
    if(extra.zones > 0)
    {
      // We need to make a zone slice array that contains selectedZonesView, plus extra indices.
      if(m_zoneSlice.size() == 0)
      {
        const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
        const auto n = selectedZonesView.size() + extra.zones;
        m_zoneSlice = axom::Array<axom::IndexType>(n, n, allocatorID);
        view = m_zoneSlice.view();
        axom::copy(view.data(),
                   selectedZonesView.data(),
                   sizeof(axom::IndexType) * selectedZonesView.size());
        axom::for_all<ExecSpace>(
          selectedZonesView.size(),
          n,
          AXOM_LAMBDA(axom::IndexType index) { view[index] = 0; });
      }
      view = m_zoneSlice.view();
    }
    else
    {
      view = selectedZonesView;
    }
    return view;
  }

  /*!
   * \brief Return a Sizes object initialized from the options.
   *
   * \param n_options The options node that contains extra sizes.
   *
   * \return A Sizes object that contains extra sizes. Values not present in the options will be 0.
   */
  Sizes getExtra(const conduit::Node &n_options) const
  {
    Sizes extra {};
    if(n_options.has_path("extra/nodes"))
    {
      extra.nodes = std::max(0, n_options["extra/nodes"].to_int());
    }
    if(n_options.has_path("extra/zones"))
    {
      extra.zones = std::max(0, n_options["extra/zones"].to_int());
    }
    if(n_options.has_path("extra/connectivity"))
    {
      extra.connectivity = std::max(0, n_options["extra/connectivity"].to_int());
    }
    return extra;
  }

  /*!
   * \brief Make node map and node slice information for the selected zones but
   *        do not limit the selection to only the used nodes.
   *
   * \param selectedZonesView The selected zones.
   * \param extra A Sizes object containing any extra sizes to use for allocation.
   * \param[out] old2new An array that contains the new node id for each node in the mesh.
   * \param[out] nodeSlice An array that contains the node ids that will be selected
   *                       from the input mesh when making coordsets and vertex fields.
   *
   * \return A Sizes object that contains the size of the nodes,zones,connectivity
   *         (excluding extra) for the output mesh.
   *
   * \note old2new is not used in this method.
   */
  Sizes nodeMap(const SelectedZonesView &selectedZonesView,
                const Sizes &extra,
                axom::Array<ConnectivityType> &AXOM_UNUSED_PARAM(old2new),
                axom::Array<axom::IndexType> &nodeSlice) const
  {
    AXOM_ANNOTATE_SCOPE("nodeMap");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    const auto nnodes = m_coordsetView.numberOfNodes();

    // Figure out the topology size based on selected zones.
    axom::ReduceSum<ExecSpace, int> connsize_reduce(0);
    const TopologyView deviceTopologyView(m_topologyView);
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto zone = deviceTopologyView.zone(zoneIndex);
        connsize_reduce += zone.numberOfNodes();
      });
    const auto newConnSize = connsize_reduce.get();

    Sizes sizes {};
    sizes.nodes = nnodes;
    sizes.zones = selectedZonesView.size();
    sizes.connectivity = newConnSize;

    nodeSlice =
      axom::Array<axom::IndexType>(sizes.nodes + extra.nodes, sizes.nodes + extra.nodes, allocatorID);
    auto nodeSliceView = nodeSlice.view();
    axom::for_all<ExecSpace>(
      sizes.nodes + extra.nodes,
      AXOM_LAMBDA(axom::IndexType index) {
        nodeSliceView[index] = (index < sizes.nodes) ? index : 0;
      });

    return sizes;
  }

  /*!
   * \brief Make node map and node slice information for the selected zones
   *        selecting only the used nodes.
   *
   * \param selectedZonesView The selected zones.
   * \param extra A Sizes object containing any extra sizes to use for allocation.
   * \param[out] old2new An array that contains the new node id for each node in the mesh.
   * \param[out] nodeSlice An array that contains the node ids that will be selected
   *                       from the input mesh when making coordsets and vertex fields.
   *
   * \return A Sizes object that contains the size of the nodes,zones,connectivity
   *         (excluding extra) for the output mesh.
   */
  Sizes compactNodeMap(const SelectedZonesView selectedZonesView,
                       const Sizes &extra,
                       axom::Array<ConnectivityType> &old2new,
                       axom::Array<axom::IndexType> &nodeSlice) const
  {
    AXOM_ANNOTATE_SCOPE("compactNodeMap");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // We need to figure out which nodes to keep.
    const auto nnodes = m_coordsetView.numberOfNodes();
    axom::Array<int> mask(nnodes, nnodes, allocatorID);
    auto maskView = mask.view();
    mask.fill(0);

    // Mark all the selected zones' nodes as 1. Multiple threads may write 1 to the same node.
    axom::ReduceSum<ExecSpace, int> connsize_reduce(0);
    TopologyView deviceTopologyView(m_topologyView);
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto zone = deviceTopologyView.zone(zoneIndex);
        const axom::IndexType nids = zone.numberOfNodes();
        for(axom::IndexType i = 0; i < nids; i++)
        {
          const auto nodeId = zone.getId(i);
          maskView[nodeId] = 1;
        }
        connsize_reduce += nids;
      });
    const auto newConnSize = connsize_reduce.get();

    // Count the used nodes.
    axom::ReduceSum<ExecSpace, int> mask_reduce(0);
    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType index) { mask_reduce += maskView[index]; });
    const int newNumNodes = mask_reduce.get();

    // Make a compact list of nodes.
    axom::Array<int> maskOffsets(nnodes, nnodes, allocatorID);
    auto maskOffsetsView = maskOffsets.view();
    axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

    // Make an array of original node ids that we can use to "slice" the nodal data.
    old2new = axom::Array<ConnectivityType>(nnodes, nnodes, allocatorID);
    nodeSlice =
      axom::Array<axom::IndexType>(newNumNodes + extra.nodes, newNumNodes + extra.nodes, allocatorID);
    auto old2newView = old2new.view();
    auto nodeSliceView = nodeSlice.view();
    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType index) {
        if(maskView[index] > 0)
        {
          nodeSliceView[maskOffsetsView[index]] = index;
          old2newView[index] = maskOffsetsView[index];
        }
      });
    if(extra.nodes > 0)
    {
      axom::for_all<ExecSpace>(
        nnodes,
        nnodes + extra.nodes,
        AXOM_LAMBDA(axom::IndexType index) { nodeSliceView[index] = 0; });
    }

    Sizes sizes {};
    sizes.nodes = newNumNodes;
    sizes.zones = selectedZonesView.size();
    sizes.connectivity = newConnSize;

    return sizes;
  }

  /*!
   * \brief Make the output topology for just the selected zones.
   *
   * \param selectedZonesView A view that contains the ids of the zones to extract.
   * \param dataSizes Array sizes for connectivity, sizes, etc.
   * \param extra Extra sizes for connectivity, sizes, etc.
   * \param old2newView A view that lets us map old node numbers to new node numbers.
   * \param n_topo The input topology.
   * \param n_options A node containing the options.
   * \param n_newTopo A node to contain the new topology.
   */
  virtual void makeTopology(const SelectedZonesView selectedZonesView,
                            const Sizes &dataSizes,
                            const Sizes &extra,
                            const axom::ArrayView<ConnectivityType> &old2newView,
                            const conduit::Node &n_topo,
                            const conduit::Node &n_options,
                            conduit::Node &n_newTopo) const
  {
    AXOM_ANNOTATE_SCOPE("makeTopology");
    namespace utils = axom::bump::utilities;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const std::string shape = outputShape(n_topo);
    if(shape == "polyhedron")
    {
      // TODO: Handle polyhedron shape.
      // NOTE: We could know whether we can have PH topos if the TopologyView
      //       handles PH zones. Maybe this method is separated out and partially
      //       specialized.
      SLIC_ERROR("Polyhedron is not handled yet.");
    }
    else
    {
      // The topology is not polyhedral. We'll output unstructured.
      n_newTopo["type"] = "unstructured";
      n_newTopo["coordset"] = n_topo["coordset"].as_string();
      n_newTopo["elements/shape"] = outputShape(n_topo);

      conduit::Node &n_conn = n_newTopo["elements/connectivity"];
      n_conn.set_allocator(c2a.getConduitAllocatorID());
      n_conn.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id,
                                   dataSizes.connectivity + extra.connectivity));
      auto connView = utils::make_array_view<ConnectivityType>(n_conn);

      conduit::Node &n_sizes = n_newTopo["elements/sizes"];
      n_sizes.set_allocator(c2a.getConduitAllocatorID());
      n_sizes.set(
        conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, dataSizes.zones + extra.zones));
      auto sizesView = utils::make_array_view<ConnectivityType>(n_sizes);

      conduit::Node &n_offsets = n_newTopo["elements/offsets"];
      n_offsets.set_allocator(c2a.getConduitAllocatorID());
      n_offsets.set(
        conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, dataSizes.zones + extra.zones));
      auto offsetsView = utils::make_array_view<ConnectivityType>(n_offsets);

      // Fill sizes, offsets
      const TopologyView deviceTopologyView(m_topologyView);
      axom::for_all<ExecSpace>(
        selectedZonesView.size(),
        AXOM_LAMBDA(axom::IndexType szIndex) {
          const auto zoneIndex = selectedZonesView[szIndex];
          const auto zone = deviceTopologyView.zone(zoneIndex);
          sizesView[szIndex] = zone.numberOfNodes();
        });

      if(extra.zones > 0)
      {
        axom::for_all<ExecSpace>(
          dataSizes.zones,
          dataSizes.zones + extra.zones,
          AXOM_LAMBDA(axom::IndexType index) { sizesView[index] = 0; });
      }
      axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);

      // Fill connectivity
      if(compact(n_options))
      {
        const axom::ArrayView<ConnectivityType> deviceOld2NewView(old2newView);
        axom::for_all<ExecSpace>(
          selectedZonesView.size(),
          AXOM_LAMBDA(axom::IndexType szIndex) {
            const auto zoneIndex = selectedZonesView[szIndex];
            const auto zone = deviceTopologyView.zone(zoneIndex);

            const int size = static_cast<int>(sizesView[szIndex]);
            const auto offset = offsetsView[szIndex];
            for(int i = 0; i < size; i++)
            {
              const auto oldNodeId = zone.getId(i);
              // When compact, we map node ids to the compact node ids.
              const auto newNodeId = deviceOld2NewView[oldNodeId];
              connView[offset + i] = newNodeId;
            }
          });
      }
      else
      {
        axom::for_all<ExecSpace>(
          selectedZonesView.size(),
          AXOM_LAMBDA(axom::IndexType szIndex) {
            const auto zoneIndex = selectedZonesView[szIndex];
            const auto zone = deviceTopologyView.zone(zoneIndex);

            const int size = static_cast<int>(sizesView[szIndex]);
            const auto offset = offsetsView[szIndex];
            for(int i = 0; i < size; i++)
            {
              connView[offset + i] = zone.getId(i);
            }
          });
      }
      if(extra.connectivity > 0)
      {
        axom::for_all<ExecSpace>(
          dataSizes.connectivity,
          dataSizes.connectivity + extra.connectivity,
          AXOM_LAMBDA(axom::IndexType index) { connView[index] = 0; });
      }

      // Handle shapes, if present.
      if(n_topo.has_path("elements/shapes"))
      {
        const conduit::Node &n_shapes = n_topo.fetch_existing("elements/shapes");
        auto shapesView = utils::make_array_view<ConnectivityType>(n_shapes);

        conduit::Node &n_newShapes = n_newTopo["elements/shapes"];
        n_newShapes.set_allocator(c2a.getConduitAllocatorID());
        n_newShapes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id,
                                          dataSizes.zones + extra.zones));
        auto newShapesView = utils::make_array_view<ConnectivityType>(n_newShapes);

        const SelectedZonesView deviceSelectedZonesView(selectedZonesView);
        axom::for_all<ExecSpace>(
          dataSizes.zones,
          AXOM_LAMBDA(axom::IndexType index) {
            newShapesView[index] = shapesView[deviceSelectedZonesView[index]];
          });
        if(extra.zones > 0)
        {
          axom::for_all<ExecSpace>(
            dataSizes.zones,
            dataSizes.zones + extra.zones,
            AXOM_LAMBDA(axom::IndexType index) { newShapesView[index] = 0; });
        }
      }
    }
  }

  /*!
   * \brief Make the new coordset using the blend data and the input coordset/coordsetview.
   *
   * \param nodeSlice Node slice information.
   * \param n_coordset The input coordset, which is passed for metadata.
   * \param[out] n_newCoordset The new coordset.
   */
  void makeCoordset(const SliceData &nodeSlice,
                    const conduit::Node &n_coordset,
                    conduit::Node &n_newCoordset) const
  {
    AXOM_ANNOTATE_SCOPE("makeCoordset");
    // _bump_utilities_coordsetslicer_begin
    axom::bump::CoordsetSlicer<ExecSpace, CoordsetView> cs(m_coordsetView);
    n_newCoordset.reset();
    cs.execute(nodeSlice, n_coordset, n_newCoordset);
    // _bump_utilities_coordsetslicer_end
  }

  /*!
   * \brief Make fields for the output mesh, as needed.
   *
   * \param nodeSlice Node slice information.
   * \param zoneSlice Zone slice information.
   * \param newTopoName The name of the new topology for the field.
   * \param n_fields The input fields.
   * \param n_newFields The output fields.
   */
  void makeFields(const SliceData &nodeSlice,
                  const SliceData &zoneSlice,
                  const std::string &newTopoName,
                  const conduit::Node &n_fields,
                  conduit::Node &n_newFields) const
  {
    AXOM_ANNOTATE_SCOPE("makeFields");

    // TODO: Change this to be more like ClipField::makeFields so it supports strided-structured.

    for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
    {
      const conduit::Node &n_field = n_fields[i];
      const std::string association = n_field["association"].as_string();
      conduit::Node &n_newField = n_newFields[n_field.name()];
      axom::bump::FieldSlicer<ExecSpace> fs;
      if(association == "element")
      {
        fs.execute(zoneSlice, n_field, n_newField);
      }
      else if(association == "vertex")
      {
        fs.execute(nodeSlice, n_field, n_newField);
      }
      n_newField["topology"] = newTopoName;
    }
  }

  /*!
   * \brief Get the topology name.
   *
   * \param n_input The input mesh.
   * \param n_options The options node that may contain a "topology" string.
   *
   * \return Returns the options topology name, if present. Otherwise, it returns the first topology name.
   */
  std::string topologyName(const conduit::Node &n_input, const conduit::Node &n_options) const
  {
    std::string name;
    if(n_options.has_path("topology"))
    {
      name = n_options["topology"].as_string();
    }
    else
    {
      const conduit::Node &n_topologies = n_input.fetch_existing("topologies");
      name = n_topologies[0].name();
    }
    return name;
  }

  /*!
   * \brief Return whether coordset/vertex compaction is desired.
   *
   * \param n_options The options node that may contain a "topology" string.
   *
   * \return True if compaction is on (the default), false otherwise.
   */
  bool compact(const conduit::Node &n_options) const
  {
    bool retval = true;
    if(n_options.has_path("compact"))
    {
      retval = n_options["compact"].to_int() != 0;
    }
    return retval;
  }

  /*!
   * \brief Return the name of the output shape type.
   *
   * \param n_topo The input topology node.
   *
   * \return The name of the output shape.
   */
  std::string outputShape(const conduit::Node &n_topo) const
  {
    std::string shape;
    if(n_topo["type"].as_string() == "unstructured")
    {
      shape = n_topo["elements/shape"].as_string();
    }
    else
    {
      if(TopologyView::dimension() == 3)
      {
        shape = "hex";
      }
      else if(TopologyView::dimension() == 2)
      {
        shape = "quad";
      }
      else
      {
        shape = "line";
      }
    }
    return shape;
  }

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  axom::Array<axom::IndexType> m_zoneSlice;
};

/*!
 * \brief Make a new topology and coordset by extracting certain zones from the input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The topology view type on which the algorithm will run.
 * \tparam MatsetView The matset view type on which the algorithm will run.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
class ExtractZonesAndMatset : public ExtractZones<ExecSpace, TopologyView, CoordsetView>
{
public:
  using SelectedZonesView = axom::ArrayView<axom::IndexType>;

  /*!
   * \brief Constructor
   *
   * \param topoView The input topology view.
   * \param coordsetView The input coordset view.
   * \param matsetView The input matset view.
   */
  ExtractZonesAndMatset(const TopologyView &topoView,
                        const CoordsetView &coordsetView,
                        const MatsetView &matsetView)
    : ExtractZones<ExecSpace, TopologyView, CoordsetView>(topoView, coordsetView)
    , m_matsetView(matsetView)
  { }

  /*!
   * \brief Destructor
   */
  virtual ~ExtractZonesAndMatset() = default;

  /*!
   * \brief Select zones from the input mesh by id and output them in the output mesh.
   *
   * \param selectedZonesView A view that contains the selected zone ids.
   * \param n_input The input mesh.
   * \param n_options The input options.
   * \param[out] n_output The output mesh.
   *
   * \note The \a n_options node contains a "topology" string that is selects the
   *       name of the topology to extract.
   */
  void execute(const SelectedZonesView &selectedZonesView,
               const conduit::Node &n_input,
               const conduit::Node &n_options,
               conduit::Node &n_output)
  {
    AXOM_ANNOTATE_SCOPE("ExtractZonesAndMatset");

    // Call base class to handle mesh/coordset/fields
    // _bump_utilities_extractzones_begin
    ExtractZones<ExecSpace, TopologyView, CoordsetView>::execute(selectedZonesView,
                                                                 n_input,
                                                                 n_options,
                                                                 n_output);
    // _bump_utilities_extractzones_end

    // Make new matset.
    const std::string topoName =
      ExtractZones<ExecSpace, TopologyView, CoordsetView>::topologyName(n_input, n_options);
    std::string mname = matsetName(n_input, topoName);
    if(!mname.empty())
    {
      const conduit::Node &n_matset = n_input.fetch_existing("matsets/" + mname);
      axom::bump::Options opts(n_options);

      const std::string newMatsetName = opts.matsetName(mname);
      conduit::Node &n_newMatset = n_output["matsets/" + newMatsetName];
      makeMatset(selectedZonesView, n_matset, n_newMatset);

      // Update the topology name in the matset.
      const std::string newTopoName = opts.topologyName(topoName);
      n_newMatset["topology"] = newTopoName;
    }
  }

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  /*!
   * \brief Return the matset for the named topology.
   *
   * \param n_input The input mesh.
   * \param topoName The name of the topology.
   *
   * \return The name of the matset for the topology or an empty string if no matset was found.
   */
  std::string matsetName(const conduit::Node &n_input, const std::string &topoName) const
  {
    std::string matset;
    if(n_input.has_child("matsets"))
    {
      const conduit::Node &n_matsets = n_input.fetch_existing("matsets");
      for(conduit::index_t i = 0; i < n_matsets.number_of_children(); i++)
      {
        const conduit::Node &n_matset = n_matsets[i];
        if(n_matset["topology"].as_string() == topoName)
        {
          matset = n_matset.name();
          break;
        }
      }
    }
    return matset;
  }

  /*!
   * \brief Make a new matset that covers the selected zones.
   *
   * \param selectedZonesView A view that contains the selected zones.
   * \param n_matset The input matset.
   * \param n_newMatset A node that will contain the new matset.
   */
  void makeMatset(const SelectedZonesView selectedZonesView,
                  const conduit::Node &n_matset,
                  conduit::Node &n_newMatset) const
  {
    AXOM_ANNOTATE_SCOPE("makeMatset");
    // _bump_utilities_matsetslicer_begin
    MatsetSlicer<ExecSpace, MatsetView> ms(m_matsetView);
    SliceData zSlice;
    zSlice.m_indicesView = selectedZonesView;
    ms.execute(zSlice, n_matset, n_newMatset);
    // _bump_utilities_matsetslicer_end
  }

  MatsetView m_matsetView;
};

}  // end namespace bump
}  // end namespace axom

#endif
