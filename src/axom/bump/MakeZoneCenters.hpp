// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MAKE_ZONE_CENTERS_HPP_
#define AXOM_BUMP_MAKE_ZONE_CENTERS_HPP_

#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{

/*!
 * \brief Makes a centroids field using the input topology and coordset views.
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam TopologyView The topology view type.
 * \tparam CoordsetView The coordset view type.
 *
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class MakeZoneCenters
{
public:
  /*!
   * \brief Constructor
   *
   * \param topologyView The view for the input topology.
   * \param coordsetView The view for the input coordset.
   */
  MakeZoneCenters(const TopologyView &topologyView, const CoordsetView &coordsetView)
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
               conduit::Node &n_outputField) const
  {
    const auto numZones = m_topologyView.numberOfZones();
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    // Select all zones.
    axom::Array<axom::IndexType> selectedZones(numZones, numZones, allocatorID);
    auto selectedZonesView = selectedZones.view();
    axom::for_all<ExecSpace>(
      numZones,
      AXOM_LAMBDA(axom::IndexType index) { selectedZonesView[index] = index; });
    // Make the zone centers.
    execute(selectedZonesView, n_topology, n_coordset, n_outputField);
  }

  /*!
   * \brief Create a new field from the input topology and place it in \a n_output.
   *
   * \param selectedZonesView A view that contains a list of selected zones.
   * \param n_topology The node that contains the input topology.
   * \param n_coordset The input coordset that we're blending.
   * \param[out] n_outputField The output node that will contain the new zone centers field.
   *
   * \note The coordset view must agree with the coordset in n_coordset. We pass both
   *       a view and the coordset node since the view may not be able to contain
   *       some coordset metadata and remain trivially copyable.
   *
   * \note When passing selectedZones that is a subset of the zones in the mesh,
   *       be aware that the generated field may not be the right length for the
   *       input topology. This is okay because we may be using this routine to
   *       generate a field that is repurposed some other way.
   */
  void execute(axom::ArrayView<axom::IndexType> selectedZonesView,
               const conduit::Node &n_topology,
               const conduit::Node &n_coordset,
               conduit::Node &n_outputField) const
  {
    using value_type = typename CoordsetView::value_type;
    using PointType = typename CoordsetView::PointType;
    using VectorType = axom::primal::Vector<value_type, PointType::DIMENSION>;
    namespace utils = axom::bump::utilities;

    // Get the axis names for the output components.
    std::vector<std::string> axes(utils::coordsetAxes(n_coordset));

    const auto nComponents = axes.size();
    SLIC_ASSERT(PointType::DIMENSION == nComponents);

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    n_outputField.reset();
    n_outputField["association"] = "element";
    n_outputField["topology"] = n_topology.name();
    conduit::Node &n_values = n_outputField["values"];

    // Determine output size.
    const auto outputSize = selectedZonesView.size();

    // Make output nodes using axis names from the input coordset. Make array views too.
    axom::StackArray<axom::ArrayView<value_type>, PointType::DIMENSION> compViews;
    for(size_t i = 0; i < nComponents; i++)
    {
      // Allocate data in the Conduit node and make a view.
      conduit::Node &comp = n_values[axes[i]];
      comp.set_allocator(c2a.getConduitAllocatorID());
      comp.set(conduit::DataType(utils::cpp2conduit<value_type>::id, outputSize));
      compViews[i] = utils::make_array_view<value_type>(comp);
    }

    const TopologyView deviceTopoView(m_topologyView);
    const CoordsetView deviceCoordsetView(m_coordsetView);

    // Blend the nodes in each zone to make a center point.
    axom::for_all<ExecSpace>(
      outputSize,
      AXOM_LAMBDA(axom::IndexType zi) {
        const auto zoneIndex = selectedZonesView[zi];
        const auto zone = deviceTopoView.zone(zoneIndex);
        const axom::IndexType nnodes = zone.numberOfNodes();

        VectorType blended {};

        // Blend points for this zone.
        for(IndexType i = 0; i < nnodes; i++)
        {
          const auto index = zone.getId(i);
          blended += VectorType(deviceCoordsetView[index]);
        }
        blended = blended / static_cast<value_type>(nnodes);

        // Store the point into the Conduit component arrays.
        for(int comp = 0; comp < PointType::DIMENSION; comp++)
        {
          compViews[comp][zi] = blended[comp];
        }
      });
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
};

}  // end namespace bump
}  // end namespace axom

#endif
