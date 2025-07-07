// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"

#include <conduit/conduit_blueprint.hpp>

#include <string>
#include <vector>

namespace axom
{
namespace bump
{
namespace utilities
{

std::vector<std::string> coordsetAxes(const conduit::Node &n_input)
{
  std::vector<std::string> axes;
  // Get the axis names for the output coordset. For uniform, prefer x,y,z
  // instead of i,j,k since we're making an explicit coordset.
  if(n_input.fetch_existing("type").as_string() == "uniform")
  {
    if(n_input.has_path("dims/i")) axes.push_back("x");
    if(n_input.has_path("dims/j")) axes.push_back("y");
    if(n_input.has_path("dims/k")) axes.push_back("z");
  }
  else
  {
    axes = conduit::blueprint::mesh::utils::coordset::axes(n_input);
  }
  return axes;
}

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom
