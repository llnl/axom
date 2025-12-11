// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internal.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_FIELD_INTERSECTOR_HPP_
#define AXOM_BUMP_FIELD_INTERSECTOR_HPP_

#include "axom/core.hpp"
#include "axom/bump/extraction/FieldOptions.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/utilities/utilities.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
namespace extraction
{

/*!
 * \brief This class helps TableBasedExtractor determine intersection cases and
 *        weights using a field designated by the options.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class FieldIntersector
{
public:
  using FieldType = float;
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /*!
   * \brief This is a view class for FieldIntersector that can be used in device code.
   */
  struct View
  {
    /*!
     * \brief Given a zone index and the node ids that comprise the zone, return
     *        the appropriate table case, taking into account the field and
     *        value.
     *
     * \param zoneIndex The zone index.
     * \param nodeIds A view containing node ids for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineTableCase(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                                       const ConnectivityView &nodeIds) const
    {
      axom::IndexType caseNumber = 0, numIds = nodeIds.size();
      for(IndexType i = 0; i < numIds; i++)
      {
        const auto id = nodeIds[i];
        const auto distance = m_fieldView[id] - m_fieldValue;
        caseNumber |= (distance > 0) ? (1 << i) : 0;
      }
      return caseNumber;
    }

    /*!
     * \brief Compute the weight of a clip value along an edge (id0, id1) using the clip field and value.
     *
     * \param id0 The mesh node at the start of the edge.
     * \param id1 The mesh node at the end of the edge.
     *
     * \return A parametric position t [0,1] where we locate \a clipValues in [d0,d1].
     */
    AXOM_HOST_DEVICE
    FieldType computeWeight(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                            ConnectivityType id0,
                            ConnectivityType id1) const
    {
      const FieldType d0 = m_fieldView[id0];
      const FieldType d1 = m_fieldView[id1];
      constexpr FieldType tiny = 1.e-09;
      return axom::utilities::clampVal(
        axom::utilities::abs(m_fieldValue - d0) / (axom::utilities::abs(d1 - d0) + tiny),
        FieldType(0),
        FieldType(1));
    }

    axom::ArrayView<FieldType> m_fieldView {};
    FieldType m_fieldValue {};
  };

  /*!
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const TopologyView &AXOM_UNUSED_PARAM(topologyView),
                  const CoordsetView &AXOM_UNUSED_PARAM(coordsetView),
                  const conduit::Node &n_options,
                  const conduit::Node &AXOM_UNUSED_PARAM(n_topology),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_coordset),
                  const conduit::Node &n_fields)
  {
    namespace utils = axom::bump::utilities;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Get the field name and value.
    FieldOptions opts(n_options);
    m_view.m_fieldValue = opts.value();

    // Make sure the clipField is the right data type and store access to it in the view.
    const conduit::Node &n_field = n_fields.fetch_existing(opts.field());
    const conduit::Node &n_field_values = n_field["values"];
    SLIC_ASSERT(n_field["association"].as_string() == "vertex");
    SLIC_ASSERT(!n_field_values.dtype().is_object());
    if(n_field_values.dtype().id() == utils::cpp2conduit<FieldType>::id)
    {
      // Make a view.
      m_view.m_fieldView = utils::make_array_view<FieldType>(n_field_values);
    }
    else
    {
      // Convert to FieldType.
      const IndexType n = static_cast<IndexType>(n_field_values.dtype().number_of_elements());
      m_fieldData = axom::Array<FieldType>(n, n, allocatorID);
      m_view.m_fieldView = m_fieldData.view();
      views::Node_to_ArrayView(n_field_values,
                               [&](auto clipFieldViewSrc) { copyValues(clipFieldViewSrc); });
    }
  }

  /*!
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &n_input, const conduit::Node &n_options) const
  {
    // Get the topo name.
    FieldOptions opts(n_options);
    const conduit::Node &n_fields = n_input.fetch_existing("fields");
    const conduit::Node &n_field = n_fields.fetch_existing(opts.field());
    return n_field["topology"].as_string();
  }

  /*!
   * \brief Return a new instance of the view.
   * \return A new instance of the view.
   */
  View view() const { return m_view; }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Copy values from srcView into m_fieldData.
   *
   * \param srcView The source data view.
   */
  template <typename DataView>
  void copyValues(DataView srcView)
  {
    auto clipFieldView = m_fieldData.view();
    axom::for_all<ExecSpace>(
      srcView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        clipFieldView[index] = static_cast<FieldType>(srcView[index]);
      });
  }

  axom::Array<FieldType> m_fieldData {};
  View m_view {};
};

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
