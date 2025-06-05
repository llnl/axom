// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MAKE_POINT_MESH_
#define AXOM_BUMP_MAKE_POINT_MESH_

#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/MakeZoneCenters.hpp"
#include "axom/bump/Options.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{

/*!
 * \brief Create a point mesh representation of the input mesh using zone centers as
 *        the coordinates in a new coordset.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
struct MakePointMesh
{
  /*!
   * \brief Constructor.
   *
   * \param topologyView The topology view that describes the input topology.
   * \param coordsetView The coordset view that describes the input coordset.
   */
  MakePointMesh(const TopologyView &topologyView, const CoordsetView &coordsetView)
    : m_topologyView(topologyView)
    , m_coordsetView(coordsetView)
  { }

  /*!
   * \brief Create a new field from the input topology and place it in \a n_output.
   *
   * \param n_topology The node that contains the input topology.
   * \param n_coordset The input coordset that we're blending.
   * \param[out] n_outputField The output node that will contain the new zone centers field.
   *
   * \note The coordset view must agree with the coordset in n_coordset. We pass both
   *       a view and the coordset node since the view may not be able to contain
   *       some coordset metadata and remain trivially copyable.
   */
  void execute(const conduit::Node &n_topology,
               const conduit::Node &n_coordset,
               const conduit::Node &n_options,
               conduit::Node &n_output) const
  {
    const auto numZones = m_topologyView.numberOfZones();
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    // Select all zones.
    axom::Array<axom::IndexType> selectedZones(numZones, numZones, allocatorID);
    auto selectedZonesView = selectedZones.view();
    axom::for_all<ExecSpace>(
      numZones,
      AXOM_LAMBDA(axom::IndexType index) { selectedZonesView[index] = index; });
    // Make the point mesh.
    execute(selectedZonesView, n_topology, n_coordset, n_options, n_output);
  }

  /*!
   * \brief Create a new point mesh using a subset of the input mesh zones.
   *
   * \param selectedZonesView A view that contains the selected zone ids.
   * \param n_topology A node that contains the input topology.
   * \param n_coordset A node that contains the input coordset.
   * \param n_options A node that contains options.
   * \param[out] n_output A node that will contain the new point mesh.
   */
  void execute(axom::ArrayView<axom::IndexType> selectedZonesView,
               const conduit::Node &n_topology,
               const conduit::Node &n_coordset,
               const conduit::Node &n_options,
               conduit::Node &n_output) const
  {
    AXOM_ANNOTATE_SCOPE("ConvertToPointMesh");
    namespace utils = axom::bump::utilities;
    using ConnectivityType = typename TopologyView::ConnectivityType;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    Options opts(n_options);

    // Make zone centers to use for the new coordset.
    MakeZoneCenters<ExecSpace, TopologyView, CoordsetView> zc(m_topologyView, m_coordsetView);
    conduit::Node zcfield;
    zc.execute(selectedZonesView, n_topology, n_coordset, zcfield);

    // Make the zone centers be the new coordset values in the output coordset.
    AXOM_ANNOTATE_BEGIN("allocate");
    conduit::Node &n_output_coordset = n_output["coordsets/" + opts.coordsetName(n_coordset.name())];
    n_output_coordset.reset();
    n_output_coordset["type"] = "explicit";
    n_output_coordset["values"].move(zcfield["values"]);

    // Allocate point mesh data.
    conduit::Node &n_output_topo = n_output["topologies/" + opts.topologyName(n_topology.name())];
    const auto numPoints = selectedZonesView.size();
    n_output_topo.reset();
    n_output_topo["type"] = "unstructured";
    n_output_topo["coordset"] = opts.coordsetName(n_coordset.name());
    n_output_topo["elements/shape"] = "point";
    conduit::Node &n_conn = n_output_topo["elements/connectivity"];
    n_conn.set_allocator(c2a.getConduitAllocatorID());
    n_conn.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numPoints));
    auto connectivity = utils::make_array_view<ConnectivityType>(n_conn);

    conduit::Node &n_sizes = n_output_topo["elements/sizes"];
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numPoints));
    auto sizes = utils::make_array_view<ConnectivityType>(n_sizes);

    conduit::Node &n_offsets = n_output_topo["elements/offsets"];
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numPoints));
    auto offsets = utils::make_array_view<ConnectivityType>(n_offsets);
    AXOM_ANNOTATE_END("allocate");

    // Build the point mesh
    AXOM_ANNOTATE_BEGIN("build");
    axom::for_all<ExecSpace>(
      numPoints,
      AXOM_LAMBDA(axom::IndexType index) {
        connectivity[index] = index;
        sizes[index] = 1;
      });
    axom::exclusive_scan<ExecSpace>(sizes, offsets);
    AXOM_ANNOTATE_END("build");
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
};

}  // end namespace bump
}  // end namespace axom

#endif
