// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#pragma once

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/sidre/core/ConduitMemory.hpp"

#include <conduit/conduit.hpp>

#include <string>

namespace axom
{
namespace bump
{

/*!
 * \brief This class computes volume for 3D and area for 2D and stores the values in a field.
 *
 * \tparam ExecSpace The execution space where the compute happens.
 * \tparam Adaptor A PrimalAdaptor.
 */
template <typename ExecSpace, typename Adaptor>
class ComputeMeasure
{
public:
  /*!
   * \brief Constructor.
   *
   * \param adaptor The adaptor object to use to compute the measure.
   */
  ComputeMeasure(Adaptor &adaptor)
    : m_adaptor(adaptor)
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
   * \brief Compute the area or volume (depending on shape dimension) and store
   *        it in the field.
   *
   * \param topoName The topology name for the field.
   * \param n_field The node that will contain the new field.
   */
  void execute(const std::string &topoName, conduit::Node &n_field)
  {
    const auto conduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(getAllocatorID());

    n_field["topology"] = topoName;
    n_field["association"] = "element";
    conduit::Node &n_values = n_field["values"];
    n_values.set_allocator(conduitAllocatorId);
    n_values.set(conduit::DataType::float64(m_adaptor.numberOfZones()));
    auto valuesView = bump::utilities::make_array_view<double>(n_values);

    // Use the Adaptor on device to compute area or volume.
    const Adaptor deviceAdaptor(m_adaptor);
    axom::for_all<ExecSpace>(
      deviceAdaptor.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto shape = deviceAdaptor.getShape(zoneIndex);

        double value = 0.;
        if constexpr(Adaptor::dimension() == 3)
        {
          value = shape.volume();
        }
        else if constexpr(Adaptor::dimension() == 2)
        {
          value = shape.area();
        }

        valuesView[zoneIndex] = value;
      });
  }

private:
  Adaptor m_adaptor;
  int m_allocator_id;
};

}  // end namespace bump
}  // end namespace axom
