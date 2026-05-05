// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/bump/views/MaterialView.hpp"

namespace axom
{
namespace bump
{
namespace views
{
MaterialInformation materials(const conduit::Node &matset)
{
  MaterialInformation info;
  if(matset.has_child("material_map"))
  {
    const conduit::Node &mm = matset["material_map"];
    for(conduit::index_t i = 0; i < mm.number_of_children(); i++)
    {
      info.push_back(Material {mm[i].to_int(), mm[i].name()});
    }
  }
  return info;
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom
