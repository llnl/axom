// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/bump/views/dispatch_utilities.hpp"
#include <conduit/conduit_blueprint_mesh.hpp>

namespace axom
{
namespace bump
{
namespace views
{

void verify(const conduit::Node &obj, const std::string &protocol)
{
  conduit::Node info;
  if(protocol.empty())
  {
    // Check the mesh
    if(!conduit::blueprint::mesh::verify(obj, info))
    {
      SLIC_ERROR(info.to_summary_string());
    }
  }
  else
  {
    // Protocol is not empty so check a specific thing in the mesh.
    if(!conduit::blueprint::mesh::verify(protocol, obj, info))
    {
      SLIC_ERROR(info.to_summary_string());
    }
  }
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom
