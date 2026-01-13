// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/Options.hpp"
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
    SLIC_ERROR("The input node contains multiple domains. Pass a single domain at a time instead.");
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
  axom::bump::Options options(n_options);

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
  conduit::Node &newFields = n_newDomain["fields"];

  // Execute the algorithm on the domain.
  if(n_domain.has_path("state"))
  {
    copyState(n_domain["state"], n_newDomain["state"]);
  }
  if(n_domain.has_path("fields"))
  {
    executeDomain(*n_topo,
                  *n_coordset,
                  n_domain["fields"],
                  n_matset,
                  n_options,
                  newTopo,
                  newCoordset,
                  newFields,
                  newMatset);
    updateNames(n_topo->name(),
                newTopoName,
                n_coordset->name(),
                newCoordsetName,
                matset,
                newMatsetName,
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
    executeDomain(*n_topo,
                  *n_coordset,
                  n_fields,
                  n_matset,
                  n_options,
                  newTopo,
                  newCoordset,
                  newFields,
                  newMatset);
    updateNames(n_topo->name(),
                newTopoName,
                n_coordset->name(),
                newCoordsetName,
                matset,
                newMatsetName,
                newTopo,
                newCoordset,
                newFields,
                newMatset);
  }
}

void MIRAlgorithm::updateNames(const std::string &origTopoName,
                               const std::string &newTopoName,
                               const std::string &origCoordsetName,
                               const std::string &newCoordsetName,
                               const std::string &AXOM_UNUSED_PARAM(origMatsetName),
                               const std::string &AXOM_UNUSED_PARAM(newMatsetName),
                               conduit::Node &n_newTopo,
                               conduit::Node &AXOM_UNUSED_PARAM(n_newCoordset),
                               conduit::Node &n_newFields,
                               conduit::Node &n_newMatset)
{
  // If the coordset was renamed in the output, make sure it the new topology references that new name.
  if(origCoordsetName != newCoordsetName)
  {
    n_newTopo["coordset"] = newCoordsetName;
  }
  // If the topology was renamed in the output, make sure the matset and any fields reference that new name.
  if(origTopoName != newTopoName)
  {
    n_newMatset["topology"] = newTopoName;

    for(conduit::index_t i = 0; i < n_newFields.number_of_children(); i++)
    {
      conduit::Node &n_field = n_newFields[i];
      if(n_field["topology"].as_string() == origTopoName)
      {
        n_field["topology"] = newTopoName;
      }
    }
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
  axom::bump::utilities::copy<axom::SEQ_EXEC>(n_host, n);
  n_host.to_summary_string_stream(std::cout, options);
}

void MIRAlgorithm::saveMesh(const conduit::Node &n_mesh, const std::string &filebase) const
{
  // Make sure data are on host.
  conduit::Node n_mesh_host;
  axom::bump::utilities::copy<axom::SEQ_EXEC>(n_mesh_host, n_mesh);

  // Check the mesh we're saving.
  conduit::Node info;
  if(!conduit::blueprint::mesh::verify(n_mesh_host, info))
  {
    printNode(n_mesh_host);
  }

  conduit::relay::io::save(n_mesh_host, filebase + ".yaml", "yaml");
#if defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(n_mesh_host, filebase, "hdf5");
#endif
}

std::string MIRAlgorithm::localPath(const conduit::Node &obj) const
{
  std::string path(obj.path());
  const auto dpos = path.find("domain");
  const auto spos = path.find("/");
  if(dpos == 0 && spos != std::string::npos && spos > dpos && obj.parent() != nullptr)
  {
    path = path.substr(spos + 1, path.size() - spos - 1);
  }
  return path;
}

}  // namespace mir
}  // namespace axom
