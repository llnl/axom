// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MAKE_ZONE_VOLUMES_HPP_
#define AXOM_BUMP_MAKE_ZONE_VOLUMES_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/PrimalAdaptor.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

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
    , m_allocator_id(axom::execution_space<ExecSpace>::allocatorID())
  { }

  /*!
   * \brief Set the allocator id to use when allocating memory.
   *
   * \param allocator_id The allocator id to use when allocating memory.
   */
  void setAllocatorID(int allocator_id)
  {
    SLIC_ERROR_IF(!axom::isValidAllocatorID(allocator_id), "Invalid allocator id.");
    SLIC_ERROR_IF(!axom::execution_space<ExecSpace>::usesAllocId(allocator_id),
                  "Allocator id is not compatible with execution space.");
    m_allocator_id = allocator_id;
  }

  /*!
   * \brief Get the allocator id to use when allocating memory.
   *
   * \return The allocator id to use when allocating memory.
   */
  int getAllocatorID() const { return m_allocator_id; }

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
    namespace utils = axom::bump::utilities;
    const auto conduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(getAllocatorID());

    // Determine output size.
    const auto outputSize = m_topologyView.numberOfZones();

    // Make output field.
    n_outputField.reset();
    n_outputField["association"] = "element";
    n_outputField["topology"] = n_topology.name();
    conduit::Node &n_values = n_outputField["values"];
    n_values.set_allocator(conduitAllocatorId);
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
  int m_allocator_id;
};

}  // end namespace bump
}  // end namespace axom

#endif
