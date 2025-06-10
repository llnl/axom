// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_EXTRUDE_MESH_HPP_
#define AXOM_BUMP_EXTRUDE_MESH_HPP_

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/slic.hpp"

#include <conduit.hpp>

#include <iostream>

namespace axom
{
namespace bump
{

/*!
 * \brief Extrude 2D Blueprint topologies that consist of tri/quad zones.
 *
 * \tparam ExecSpace The execution space where the algorithm will execute.
 * \tparam TopologyView The topology view that wraps the Blueprint topology.
 * \tparam CoordsetView The coordset view that wraps the Blueprint coordset.
 *
 * \note Future work: fields, making PH zones from polygons.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class ExtrudeMesh
{
public:
  static_assert(CoordsetView::dimension() == 2, "Coordset must contain 2D coordinates.");

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view.
   * \param coordsetView The coordset view.
   */
  ExtrudeMesh(const TopologyView &topoView, const CoordsetView &coordsetView)
    : m_topologyView(topoView)
    , m_coordsetView(coordsetView)
  { }

  /*!
   * \brief Execute the extrusion algorithm.
   *
   * \param n_mesh A Conduit node containing the mesh (coordset, topologies, matset, etc.)
   * \param n_options A Conduit node containing options that control the algorithm.
   * \param n_output The Conduit node that will contain the extruded mesh.
   *
   * \verbatim
   * Supported options:
   *
   * nz: 10
   * z0: 0.
   * z1: 1.
   * topologyName: mesh
   * outputTopologyName: newmesh
   * outputCoordsetName: newcoordset
   * outputMatsetName: newmatset
   * \endverbatim
   */
  void execute(const conduit::Node &n_mesh, const conduit::Node &n_options, conduit::Node &n_output) const
  {
    namespace utils = axom::bump::utilities;
    namespace views = axom::bump::views;
    AXOM_ANNOTATE_SCOPE("Extrude::execute");
    int nz = n_options.has_child("nz") ? n_options["nz"].to_int() : 2;

    // Get some properties from the options.
    const std::string srcTopoName =
      n_options.has_child("topologyName") ? n_options["topologyName"].as_string() : "main";
    const conduit::Node &n_srcTopo = n_mesh.fetch_existing("topologies/" + srcTopoName);
    const std::string srcCoordsetName = n_srcTopo["coordset"].as_string();
    const std::string outputTopoName = n_options.has_child("outputTopologyName")
      ? n_options["outputTopologyName"].as_string()
      : srcTopoName;
    const std::string outputCoordsetName = n_options.has_child("outputCoordsetName")
      ? n_options["outputCoordsetName"].as_string()
      : srcCoordsetName;

    // Count the number of values needed to store the connectivity.
    AXOM_ANNOTATE_BEGIN("counts");
    const TopologyView topoView = m_topologyView;
    axom::ReduceSum<ExecSpace, int> connSizeReduce(0);
    axom::ReduceBitOr<ExecSpace, int> zoneTypeReduce(0);
    axom::for_all<ExecSpace>(
      m_topologyView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zi) {
        const auto zone = topoView.zone(zi);
        switch(zone.id())
        {
        case views::Tri_ShapeID:
          zoneTypeReduce |= static_cast<int>(1 << views::Wedge_ShapeID);
          connSizeReduce += 6;
          break;
        case views::Quad_ShapeID:
          zoneTypeReduce |= static_cast<int>(1 << views::Hex_ShapeID);
          connSizeReduce += 8;
          break;
        default:
          SLIC_ASSERT("Unsupported zone type");
        }
      });
    AXOM_ANNOTATE_END("counts");
    const auto shapes = zoneTypeReduce.get();
    const axom::IndexType totalConnSize = static_cast<axom::IndexType>(nz - 1) * connSizeReduce.get();
    const axom::IndexType totalZones =
      static_cast<axom::IndexType>(nz - 1) * m_topologyView.numberOfZones();
    const axom::IndexType totalNodes =
      static_cast<axom::IndexType>(nz) * m_coordsetView.numberOfNodes();

    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Create the new coordset.
    AXOM_ANNOTATE_BEGIN("coordset");
    const char *coordNames[] = {"values/x", "values/y", "values/z"};
    conduit::Node &n_outputCoordset = n_output["coordsets/" + outputCoordsetName];
    n_outputCoordset["type"] = "explicit";
    using value_type = typename CoordsetView::value_type;
    axom::StackArray<axom::ArrayView<value_type>, 3> values;
    for(int d = 0; d < 3; d++)
    {
      conduit::Node &n_value = n_outputCoordset[coordNames[d]];
      n_value.set_allocator(c2a.getConduitAllocatorID());
      n_value.set(conduit::DataType(utils::cpp2conduit<value_type>::id, totalNodes));
      values[d] = utils::make_array_view<value_type>(n_value);
    }
    const CoordsetView coordsetView = m_coordsetView;
    const auto nnodes = coordsetView.numberOfNodes();
    const value_type z0 = n_options.has_child("z0") ? n_options["z0"].to_double() : 0.;
    const value_type z1 = n_options.has_child("z1") ? n_options["z1"].to_double() : (nz - 1);
    for(int z = 0; z < nz; z++)
    {
      auto tz = static_cast<value_type>(z) / static_cast<value_type>(nz - 1);
      value_type zc = z0 + tz * (z1 - z0);
      axom::for_all<ExecSpace>(
        coordsetView.numberOfNodes(),
        AXOM_LAMBDA(int srcNodeIndex) {
          const auto destNodeIndex = z * nnodes + srcNodeIndex;
          const auto pt = coordsetView[srcNodeIndex];

          values[0][destNodeIndex] = pt[0];
          values[1][destNodeIndex] = pt[1];
          values[2][destNodeIndex] = zc;
        });
    }
    AXOM_ANNOTATE_END("coordset");

    // Create the new topology.
    AXOM_ANNOTATE_BEGIN("topology");
    conduit::Node &n_outputTopo = n_output["topologies/" + outputTopoName];
    n_outputTopo["type"] = "unstructured";
    n_outputTopo["coordset"] = outputCoordsetName;
    int count = axom::utilities::popcount(shapes);
    if(count > 1)
    {
      if(axom::utilities::bitIsSet(shapes, views::Wedge_ShapeID))
        n_outputTopo["elements/shape_map/wedge"] = views::Wedge_ShapeID;

      if(axom::utilities::bitIsSet(shapes, views::Hex_ShapeID))
        n_outputTopo["elements/shape_map/hex"] = views::Hex_ShapeID;

      n_outputTopo["elements/shape"] = "mixed";
    }
    else
    {
      if(axom::utilities::bitIsSet(shapes, views::Wedge_ShapeID))
        n_outputTopo["elements/shape"] = "wedge";

      if(axom::utilities::bitIsSet(shapes, views::Hex_ShapeID))
        n_outputTopo["elements/shape"] = "hex";
    }

    conduit::Node &n_connectivity = n_outputTopo["elements/connectivity"];
    conduit::Node &n_shapes = n_outputTopo["elements/shapes"];
    conduit::Node &n_sizes = n_outputTopo["elements/sizes"];
    conduit::Node &n_offsets = n_outputTopo["elements/offsets"];

    using ConnectivityType = typename TopologyView::ConnectivityType;
    n_connectivity.set_allocator(c2a.getConduitAllocatorID());
    n_shapes.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set_allocator(c2a.getConduitAllocatorID());

    n_connectivity.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, totalConnSize));
    n_shapes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, totalZones));
    n_sizes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, totalZones));
    n_offsets.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, totalZones));

    auto connView = utils::make_array_view<ConnectivityType>(n_connectivity);
    auto shapesView = utils::make_array_view<ConnectivityType>(n_shapes);
    auto sizesView = utils::make_array_view<ConnectivityType>(n_sizes);
    auto offsetsView = utils::make_array_view<ConnectivityType>(n_offsets);

    // Q: Would it be better to reverse loops like this?
    axom::IndexType zOffset = 0;
    for(int i = 0; i < nz - 1; i++)
    {
      axom::for_all<ExecSpace>(
        topoView.numberOfZones(),
        AXOM_LAMBDA(axom::IndexType zi) {
          const auto zone = topoView.zone(zi);
          const auto destIndex = zOffset + zi;
          switch(zone.id())
          {
          case views::Tri_ShapeID:
          {
            shapesView[destIndex] = views::Wedge_ShapeID;
            sizesView[destIndex] = 6;
          }
          break;
          case views::Quad_ShapeID:
            shapesView[destIndex] = views::Hex_ShapeID;
            sizesView[destIndex] = 8;
            break;
          default:
            SLIC_ASSERT("Unsupported zone type");
          }
        });
      zOffset += topoView.numberOfZones();
    }
    if(count <= 1)
    {
      n_outputTopo.remove("elements/shapes");
    }
    // Compute offsets.
    axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);
    // Add connectivity
    zOffset = 0;
    for(int i = 0; i < nz - 1; i++)
    {
      axom::for_all<ExecSpace>(
        topoView.numberOfZones(),
        AXOM_LAMBDA(axom::IndexType zi) {
          const auto zone = topoView.zone(zi);
          const auto offset = offsetsView[zOffset + zi];
          const auto lowNodeOffset = i * nnodes;
          const auto highNodeOffset = lowNodeOffset + nnodes;
          switch(zone.id())
          {
          case views::Tri_ShapeID:
          {
            connView[offset] = lowNodeOffset + zone.getId(0);
            connView[offset + 1] = lowNodeOffset + zone.getId(1);
            connView[offset + 2] = lowNodeOffset + zone.getId(2);
            connView[offset + 3] = highNodeOffset + zone.getId(0);
            connView[offset + 4] = highNodeOffset + zone.getId(1);
            connView[offset + 5] = highNodeOffset + zone.getId(2);
          }
          break;
          case views::Quad_ShapeID:
            connView[offset] = lowNodeOffset + zone.getId(0);
            connView[offset + 1] = lowNodeOffset + zone.getId(1);
            connView[offset + 2] = lowNodeOffset + zone.getId(2);
            connView[offset + 3] = lowNodeOffset + zone.getId(3);
            connView[offset + 4] = highNodeOffset + zone.getId(0);
            connView[offset + 5] = highNodeOffset + zone.getId(1);
            connView[offset + 6] = highNodeOffset + zone.getId(2);
            connView[offset + 7] = highNodeOffset + zone.getId(3);
            break;
          default:
            SLIC_ASSERT("Unsupported zone type");
          }
        });
      zOffset += topoView.numberOfZones();
    }
    AXOM_ANNOTATE_END("topology");

    // Extrude the matset if one exists.
    std::string matsetName = findMatset(n_mesh, srcTopoName);
    if(!matsetName.empty())
    {
      const conduit::Node &n_srcMatset = n_mesh.fetch_existing("matsets/" + matsetName);
      std::string outputMatsetName = n_options.has_child("outputMatsetName")
        ? n_options["outputMatsetName"].as_string()
        : matsetName;
      conduit::Node &n_outputMatset = n_output["matsets/" + outputMatsetName];
      extrudeMatset(n_srcMatset, n_outputMatset, outputTopoName, nz);
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Find a matset for the specified topology, if a matset exists.
   *
   * \param n_mesh The node that contains the mesh.
   * \param topoName The name of the topology whose matset we want to find.
   *
   * \return The name of the matset associated with the input topology, or an
   *         empty string if no matset exists.
   */
  std::string findMatset(const conduit::Node &n_mesh, const std::string &topoName) const
  {
    std::string matset;
    if(n_mesh.has_child("matsets"))
    {
      const conduit::Node &n_matsets = n_mesh["matsets"];
      for(conduit::index_t i = 0; i < n_matsets.number_of_children(); i++)
      {
        if(n_matsets[i]["topology"].as_string() == topoName)
        {
          matset = n_matsets[i].name();
          break;
        }
      }
    }
    return matset;
  }

  /*!
   * \brief Extrude the source matset (unibuffer only right now).
   *
   * \param n_srcMatset The node that contains the source matset.
   * \param n_outputMatset The node that will contain the extruded matset.
   * \param outputTopoName The name of the output topology associated with the output matset.
   * \param nz The number of coordinate plane repetitions in the extrusion. Zones are replicated (nz-1) times.
   *
   * \note In future work, we could use matset views to support more input matset types.
   */
  void extrudeMatset(const conduit::Node &n_srcMatset,
                     conduit::Node &n_outputMatset,
                     const std::string &outputTopoName,
                     int nz) const
  {
    namespace utils = axom::bump::utilities;
    namespace views = axom::bump::views;
    AXOM_ANNOTATE_SCOPE("matset");

    const conduit::Node &n_materialMap = n_srcMatset["material_map"];

    const conduit::Node &n_src_volume_fractions = n_srcMatset["volume_fractions"];
    const conduit::Node &n_src_material_ids = n_srcMatset["material_ids"];
    const conduit::Node &n_src_indices = n_srcMatset["indices"];
    const conduit::Node &n_src_sizes = n_srcMatset["sizes"];
    const conduit::Node &n_src_offsets = n_srcMatset["offsets"];

    // Make new matset nodes
    n_outputMatset["material_map"].set(n_materialMap);
    n_outputMatset["topology"].set(outputTopoName);

    conduit::Node &n_material_ids = n_outputMatset["material_ids"];
    conduit::Node &n_volume_fractions = n_outputMatset["volume_fractions"];
    conduit::Node &n_indices = n_outputMatset["indices"];
    conduit::Node &n_sizes = n_outputMatset["sizes"];
    conduit::Node &n_offsets = n_outputMatset["offsets"];

    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set_allocator(c2a.getConduitAllocatorID());

    n_material_ids.set(conduit::DataType(n_src_material_ids.dtype().id(),
                                         n_src_material_ids.dtype().number_of_elements() * (nz - 1)));
    n_volume_fractions.set(
      conduit::DataType(n_src_volume_fractions.dtype().id(),
                        n_src_volume_fractions.dtype().number_of_elements() * (nz - 1)));
    n_indices.set(conduit::DataType(n_src_indices.dtype().id(),
                                    n_src_indices.dtype().number_of_elements() * (nz - 1)));
    n_sizes.set(conduit::DataType(n_src_sizes.dtype().id(),
                                  n_src_sizes.dtype().number_of_elements() * (nz - 1)));
    n_offsets.set(conduit::DataType(n_src_offsets.dtype().id(),
                                    n_src_offsets.dtype().number_of_elements() * (nz - 1)));

    // Extrude the old arrays into the new arrays.
    views::FloatNode_to_ArrayView_same(n_src_volume_fractions,
                                       n_volume_fractions,
                                       [&](auto srcVolumeFractionsView, auto volumeFractionsView) {
                                         views::IndexNode_to_ArrayView_same(
                                           n_src_material_ids,
                                           n_src_indices,
                                           n_src_sizes,
                                           n_src_offsets,
                                           n_material_ids,
                                           n_indices,
                                           n_sizes,
                                           n_offsets,
                                           [&](auto srcMaterialIdsView,
                                               auto srcIndicesView,
                                               auto srcSizesView,
                                               auto srcOffsetsView,
                                               auto materialIdsView,
                                               auto indicesView,
                                               auto sizesView,
                                               auto offsetsView) {
                                             copyMatsetData(srcVolumeFractionsView,
                                                            volumeFractionsView,
                                                            srcMaterialIdsView,
                                                            srcIndicesView,
                                                            srcSizesView,
                                                            srcOffsetsView,
                                                            materialIdsView,
                                                            indicesView,
                                                            sizesView,
                                                            offsetsView,
                                                            nz);
                                           });
                                       });
  }

  /*!
   * \brief Copy material data of various types.
   *
   * \tparam FloatView The type of view wrapping float data.
   * \tparam IndexView The type of view wrapping index/int data.
   *
   * \note This is done as a template member method because we cannot invoke
   *       device lambda functions directly from inside other lambdas when
   *       using CUDA.
   */
  template <typename FloatView, typename IndexView>
  void copyMatsetData(FloatView srcVolumeFractionsView,
                      FloatView volumeFractionsView,
                      IndexView srcMaterialIdsView,
                      IndexView srcIndicesView,
                      IndexView srcSizesView,
                      IndexView srcOffsetsView,
                      IndexView materialIdsView,
                      IndexView indicesView,
                      IndexView sizesView,
                      IndexView offsetsView,
                      int nz) const
  {
    const auto nzones = srcSizesView.size();
    for(int z = 0; z < nz - 1; z++)
    {
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType zi) {
          const auto idxSize = srcIndicesView.size();
          const auto offset = srcOffsetsView[zi];
          const auto destZone = z * nzones + zi;
          sizesView[destZone] = srcSizesView[offset];
          offsetsView[destZone] = z * idxSize + srcOffsetsView[offset];
          using counter_type = typename decltype(srcSizesView)::value_type;
          for(counter_type i = 0; i < srcSizesView[offset]; i++)
          {
            const auto idx = srcIndicesView[offset + i];
            const auto outIdx = offsetsView[destZone] + i;
            volumeFractionsView[outIdx] = srcVolumeFractionsView[idx];
            materialIdsView[outIdx] = srcMaterialIdsView[idx];
            indicesView[outIdx] = outIdx;
          }
        });
    }
  }

  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
};

}  // namespace bump
}  // namespace axom

#endif
