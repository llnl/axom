// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_CLIP_FIELD_FILTER_HPP_
#define AXOM_MIR_CLIP_FIELD_FILTER_HPP_

#include "axom/core.hpp"
#include "axom/mir.hpp"

namespace axom
{
namespace mir
{
namespace clipping
{
/**
 * \brief This class runs a ClipField algorithm.
 *
 */
class ClipFieldFilter
{
public:
  /// Constructor
  ClipFieldFilter() : m_runtime(axom::runtime_policy::Policy::seq) { }

  /**
   * \brief Set the runtime policy.
   * \param value The new runtime policy.
   */
  void setRuntime(axom::runtime_policy::Policy value) { m_runtime = value; }

  /**
   * \brief Get the runtime policy.
   * 
   */
  axom::runtime_policy::Policy runtime() const { return m_runtime; }

  /**
   * \brief Execute the clipping operation using the specified options.
   *
   * \param[in] n_input The Conduit node that contains the topology, coordsets, and fields.
   * \param[in] n_options A Conduit node that contains clipping options.
   * \param[out] n_output A Conduit node that will hold the clipped output mesh. This should be a different node from \a n_input.
   *
   * \note The clipField field must currently be vertex-associated.
   */
  void execute(const conduit::Node &n_input, const conduit::Node &n_options, conduit::Node &n_output);

  /**
   * \brief Execute the clipping operation using the specified options.
   *
   * \param[in] n_topo The node that contains the input mesh topology.
   * \param[in] n_coordset The node that contains the input mesh coordset.
   * \param[in] n_fields The node that contains the input fields.
   * \param[in] n_options A Conduit node that contains clipping options.
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   *
   * \note The clipField field must currently be vertex-associated. Also, the output topology will be an unstructured topology with mixed shape types.
   */
  void execute(const conduit::Node &n_topo,
               const conduit::Node &n_coordset,
               const conduit::Node &n_fields,
               const conduit::Node &n_options,
               conduit::Node &n_newTopo,
               conduit::Node &n_newCoordset,
               conduit::Node &n_newFields);

private:
  axom::runtime_policy::Policy m_runtime;
};

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom

#endif
