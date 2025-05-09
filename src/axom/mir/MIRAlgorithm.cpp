// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/MIROptions.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"
#include "axom/slic.hpp"

#include <conduit_blueprint_mesh.hpp>
#include <conduit_blueprint_mesh_utils.hpp>
#include <conduit_relay_io.hpp>
#include <conduit_relay_io_blueprint.hpp>

namespace axom
{
namespace mir
{
void MIRAlgorithm::execute(const conduit::Node &n_input,
                           const conduit::Node &n_options,
                           conduit::Node &n_output)
{
  const auto domains = conduit::blueprint::mesh::domains(n_input);
  if(domains.size() > 1)
  {
    // Handle multiple domains
    for(const auto &dom_ptr : domains)
    {
      const conduit::Node &n_domain = *dom_ptr;
      conduit::Node &n_newDomain = n_output.append();

      executeSetup(n_domain, n_options, n_newDomain);
    }
  }
  else if(domains.size() > 0)
  {
    // Handle single domain
    const conduit::Node &n_domain = *domains[0];
    executeSetup(n_domain, n_options, n_output);
  }
}

void MIRAlgorithm::executeSetup(const conduit::Node &n_domain,
                                const conduit::Node &n_options,
                                conduit::Node &n_newDomain)
{
  MIROptions options(n_options);

  // Get the matset that we'll operate on.
  const std::string matset = options.matset();

  // Which topology is that matset defined on?
  const conduit::Node &n_matsets = n_domain.fetch_existing("matsets");
  const conduit::Node &n_matset = n_matsets.fetch_existing(matset);
  const conduit::Node *n_topo =
    conduit::blueprint::mesh::utils::find_reference_node(n_matset, "topology");
  SLIC_ASSERT(n_topo != nullptr);

  // Which coordset is used by that topology?
  const conduit::Node *n_coordset =
    conduit::blueprint::mesh::utils::find_reference_node(*n_topo, "coordset");
  SLIC_ASSERT(n_coordset != nullptr);

  // Get the names of the output items.
  const std::string newTopoName = options.topologyName(n_topo->name());
  const std::string newCoordsetName = options.coordsetName(n_coordset->name());
  const std::string newMatsetName = options.matsetName(matset);

  // Make some new nodes in the output.
  conduit::Node &newCoordset = n_newDomain["coordsets/" + newCoordsetName];
  conduit::Node &newTopo = n_newDomain["topologies/" + newTopoName];
  newTopo["coordset"] = newCoordsetName;
  conduit::Node &newMatset = n_newDomain["matsets/" + newMatsetName];
  newMatset["topology"] = newTopoName;

  // Execute the algorithm on the domain.
  if(n_domain.has_path("state")) copyState(n_domain["state"], n_newDomain["state"]);
  if(n_domain.has_path("fields"))
  {
    conduit::Node &newFields = n_newDomain["fields"];
    executeDomain(*n_topo,
                  *n_coordset,
                  n_domain["fields"],
                  n_matset,
                  n_options,
                  newTopo,
                  newCoordset,
                  newFields,
                  newMatset);
  }
  else
  {
    // There are no input fields, but make sure n_fields has a name.
    conduit::Node tmp;
    conduit::Node &n_fields = tmp["fields"];
    // MIR is likely to output some created fields.
    conduit::Node &newFields = n_newDomain["fields"];
    executeDomain(*n_topo,
                  *n_coordset,
                  n_fields,
                  n_matset,
                  n_options,
                  newTopo,
                  newCoordset,
                  newFields,
                  newMatset);
  }
}

void MIRAlgorithm::copyState(const conduit::Node &srcState, conduit::Node &destState) const
{
  for(conduit::index_t i = 0; i < srcState.number_of_children(); i++)
    destState[srcState[i].name()].set(srcState[i]);
}

void MIRAlgorithm::printNode(const conduit::Node &n) const
{
  conduit::Node options;
  options["num_children_threshold"] = 10000;
  options["num_elements_threshold"] = 10000;

  // Make sure data are on host.
  conduit::Node n_host;
  axom::mir::utilities::blueprint::copy<axom::SEQ_EXEC>(n_host, n);
  n_host.to_summary_string_stream(std::cout, options);
}

void MIRAlgorithm::saveMesh(const conduit::Node &n_mesh, const std::string &filebase) const
{
  // Make sure data are on host.
  conduit::Node n_mesh_host;
  axom::mir::utilities::blueprint::copy<axom::SEQ_EXEC>(n_mesh_host, n_mesh);

  conduit::relay::io::save(n_mesh_host, filebase + ".yaml", "yaml");
#if defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(n_mesh_host, filebase, "hdf5");
#endif
}

}  // namespace mir
}  // namespace axom
