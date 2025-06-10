// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/bump/views/UnstructuredTopologyMixedShapeView.hpp"

namespace axom
{
namespace bump
{
namespace views
{
ShapeMap buildShapeMap(const conduit::Node &n_topo,
                       axom::Array<IndexType> &values,
                       axom::Array<IndexType> &ids,
                       int allocatorID)
{
  // Make the map from the Conduit shape_map. Use std::map to sort the key values.
  // The shape_map nodes should be in host memory since the int values can fit
  // in a Conduit::Node.
  std::map<IndexType, IndexType> sm;
  const conduit::Node &n_shape_map = n_topo.fetch_existing("elements/shape_map");
  for(conduit::index_t i = 0; i < n_shape_map.number_of_children(); i++)
  {
    const auto value = static_cast<IndexType>(n_shape_map[i].to_int());
    sm[value] = axom::bump::views::shapeNameToID(n_shape_map[i].name());
  }

  // Store the map in 2 vectors so data are contiguous.
  std::vector<IndexType> valuesvec, idsvec;
  valuesvec.reserve(sm.size());
  idsvec.reserve(sm.size());
  for(auto it = sm.begin(); it != sm.end(); it++)
  {
    valuesvec.push_back(it->first);
    idsvec.push_back(it->second);
  }

  // Copy the map values to the device memory.
  const axom::IndexType n = static_cast<axom::IndexType>(sm.size());
  values = axom::Array<IndexType>(axom::ArrayOptions::Uninitialized(), n, n, allocatorID);
  ids = axom::Array<IndexType>(axom::ArrayOptions::Uninitialized(), n, n, allocatorID);
  axom::copy(values.data(), valuesvec.data(), n * sizeof(IndexType));
  axom::copy(ids.data(), idsvec.data(), n * sizeof(IndexType));

  return ShapeMap(values.view(), ids.view());
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom
