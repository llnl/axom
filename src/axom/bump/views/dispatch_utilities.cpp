// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/bump/views/dispatch_utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"

#include <conduit/conduit_blueprint_mesh.hpp>

namespace axom
{
namespace bump
{
namespace views
{

static bool needToPerformCheck(const conduit::Node &obj, const std::string &protocol)
{
  bool performCheck = true;

  // NOTE: Some verification functions in conduit::blueprint::mesh::verify access
  //       array data. We need to skip verification for such things for now when
  //       data are on the GPU.

  if(protocol == "topology")
  {
    const conduit::Node &n_topo = obj;
    if(n_topo["type"].as_string() == "unstructured" && n_topo.has_path("elements/shape"))
    {
      if(n_topo["elements/shape"].as_string() == "mixed" && n_topo.has_path("elements/shapes"))
      {
        performCheck = !axom::bump::utilities::isDeviceAllocated(n_topo["elements/shapes"]);
      }
    }
  }

  return performCheck;
}

void verify(const conduit::Node &obj, const std::string &protocol)
{
  conduit::Node info;

  const bool performCheck = needToPerformCheck(obj, protocol);

  if(protocol.empty())
  {
    // Check the mesh
    if(performCheck && !conduit::blueprint::mesh::verify(obj, info))
    {
      SLIC_ERROR(info.to_summary_string());
    }
  }
  else
  {
    // Protocol is not empty so check a specific thing in the mesh.
    if(performCheck && !conduit::blueprint::mesh::verify(protocol, obj, info))
    {
      SLIC_ERROR(info.to_summary_string());
    }
  }
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom
