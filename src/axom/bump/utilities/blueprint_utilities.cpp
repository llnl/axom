// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"

#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include <conduit/conduit_blueprint.hpp>

#include <string>
#include <vector>

namespace axom
{
namespace bump
{
namespace utilities
{

std::vector<std::string> coordsetAxes(const conduit::Node& n_input)
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

namespace
{

/*!
 * \brief Copy an integer metadata array from \a path into host memory.
 *
 * \param[in] n Node to query.
 * \param[in] path Path to the metadata array within \a n.
 * \param[out] values Host vector that receives the values.
 *
 * \return true when \a path exists. In that case, the function fills \a values.
 */
bool readIndexMetadata(const conduit::Node& n, const std::string& path, std::vector<int>& values)
{
  if(!n.has_path(path))
  {
    return false;
  }
  // The Conduit node may refer to device memory, but the accessor below runs on the host.
  conduit::Node hostNode;
  axom::bump::utilities::copy<axom::SEQ_EXEC>(hostNode, n.fetch_existing(path));

  const auto acc = hostNode.as_int_accessor();
  values.resize(static_cast<std::size_t>(acc.number_of_elements()));
  for(std::size_t i = 0; i < values.size(); i++)
  {
    values[i] = acc[static_cast<conduit::index_t>(i)];
  }
  return true;
}

std::string formatIndices(const std::vector<int>& values)
{
  std::string s("[");
  for(std::size_t i = 0; i < values.size(); i++)
  {
    s += (i ? ", " : "") + std::to_string(values[i]);
  }
  return s + "]";
}

}  // end anonymous namespace

void validateVertexFieldIndexing(const conduit::Node& n_topology,
                                 const conduit::Node& n_field,
                                 const std::string& fieldName)
{
  std::vector<int> topoOffsets, topoStrides;
  const bool hasTopoOffsets = readIndexMetadata(n_topology, "elements/dims/offsets", topoOffsets);
  const bool hasTopoStrides = readIndexMetadata(n_topology, "elements/dims/strides", topoStrides);

  // StructuredTopologyView::zone() constructs adjacent corner ids by adding 1.
  // Any other i-stride produces incorrect zones.
  SLIC_ERROR_IF(
    hasTopoStrides && !topoStrides.empty() && topoStrides[0] != 1,
    axom::fmt::format(
      "Bump requires a structured topology whose i-stride is 1, but the topology has strides {}.",
      formatIndices(topoStrides)));

  std::vector<int> fieldOffsets, fieldStrides;
  const bool hasFieldOffsets = readIndexMetadata(n_field, "offsets", fieldOffsets);
  const bool hasFieldStrides = readIndexMetadata(n_field, "strides", fieldStrides);

  if(!hasFieldOffsets && !hasFieldStrides)
  {
    // Without field-specific layout metadata, topology node ids index the values directly.
    return;
  }

  SLIC_ERROR_IF(
    !hasTopoOffsets && !hasTopoStrides,
    axom::fmt::format("Field '{}' has Blueprint offsets or strides, but its topology does not. "
                      "Bump kernels use topology node ids and do not apply the field's layout.",
                      fieldName));

  SLIC_ERROR_IF(
    hasFieldOffsets && fieldOffsets != topoOffsets,
    axom::fmt::format("Field '{}' has offsets {}, but its topology has offsets {}. Bump kernels "
                      "use topology node ids and do not apply separate field offsets.",
                      fieldName,
                      formatIndices(fieldOffsets),
                      formatIndices(topoOffsets)));

  SLIC_ERROR_IF(
    hasFieldStrides && fieldStrides != topoStrides,
    axom::fmt::format("Field '{}' has strides {}, but its topology has strides {}. Bump kernels "
                      "use topology node ids and do not apply separate field strides.",
                      fieldName,
                      formatIndices(fieldStrides),
                      formatIndices(topoStrides)));
}

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom
