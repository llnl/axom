// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MAKE_ZONE_VOLUMES_HPP_
#define AXOM_BUMP_MAKE_ZONE_VOLUMES_HPP_

#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/PrimalAdaptor.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{

/*!
 * \brief Makes a new element field with zone areas or volumes (depending on dimension)
 *        using the input topology and coordset views.
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam TopologyView The topology view type.
 * \tparam CoordsetView The coordset view type.
 *
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class MakeZoneVolumes
{
public:
  using value_type = double;

  /*!
   * \brief Constructor
   *
   * \param topologyView The view for the input topology.
   * \param coordsetView The view for the input coordset.
   */
  MakeZoneVolumes(const TopologyView &topologyView, const CoordsetView &coordsetView)
    : m_topologyView(topologyView)
    , m_coordsetView(coordsetView)
  { }

  /*!
   * \brief Create a new field from the input topology and place it in \a n_output.
   *
   * \param n_topology The node that contains the input topology.
   * \param n_coordset The input coordset that we're blending.
   * \param[out] n_outputField The output node that will contain the new field.
   *
   */
  void execute(const conduit::Node &n_topology,
               const conduit::Node &AXOM_UNUSED_PARAM(n_coordset),
               conduit::Node &n_outputField) const
  {
    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    namespace utils = axom::bump::utilities;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Determine output size.
    const auto outputSize = m_topologyView.numberOfZones();

    // Make output field.
    n_outputField.reset();
    n_outputField["association"] = "element";
    n_outputField["topology"] = n_topology.name();
    conduit::Node &n_values = n_outputField["values"];
    n_values.set_allocator(c2a.getConduitAllocatorID());
    n_values.set(conduit::DataType(utils::cpp2conduit<value_type>::id, outputSize));
    auto valuesView = utils::make_array_view<value_type>(n_values);

    // _bump_utilities_makezonevolumes_begin
    // Get the zone as a primal shape and compute area or volume, as needed.
    using ShapeView = PrimalAdaptor<TopologyView, CoordsetView>;
    const ShapeView deviceShapeView {m_topologyView, m_coordsetView};
    axom::for_all<ExecSpace>(
      m_topologyView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto shape = deviceShapeView.getShape(zoneIndex);

        // Get the area or volume of the target shape (depends on the dimension).
        double amount = utils::ComputeShapeAmount<CoordsetView::dimension()>::execute(shape);

        valuesView[zoneIndex] = amount;
      });
    // _bump_utilities_makezonevolumes_end
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
};

}  // end namespace bump
}  // end namespace axom

#endif
