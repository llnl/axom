// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/IndexingPolicies.hpp"

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

// Format a strided layout for diagnostics.
template <typename Indexing>
std::string formatLayout(const Indexing& indexing)
{
  constexpr int NDIMS = Indexing::dimension();
  return axom::fmt::format(
    "offsets [{}] and strides [{}]",
    axom::fmt::join(&indexing.m_offsets[0], &indexing.m_offsets[0] + NDIMS, ", "),
    axom::fmt::join(&indexing.m_strides[0], &indexing.m_strides[0] + NDIMS, ", "));
}

// Compare a vertex field's resolved layout with the structured topology's node layout.
template <int NDIMS>
void validateStructuredVertexField(const conduit::Node& n_dims,
                                   const conduit::Node& n_hostDims,
                                   const conduit::Node& n_field,
                                   const std::string& fieldName)
{
  using Indexing = axom::bump::views::StridedStructuredIndexing<axom::IndexType, NDIMS>;
  using LogicalIndex = typename Indexing::LogicalIndex;

  // elements/dims holds zone counts. A vertex field has one more node along each axis.
  const char* axes[] = {"i", "j", "k"};
  LogicalIndex nodeDims;
  for(int d = 0; d < NDIMS; d++)
  {
    nodeDims[d] = n_hostDims.fetch_existing(axes[d]).to_index_t() + 1;
  }
  Indexing topoLayout(nodeDims);
  axom::bump::detail::setLayoutFromNode(n_dims, topoLayout);

  // StructuredTopologyView::zone() constructs adjacent corner ids by adding 1.
  SLIC_ERROR_IF(topoLayout.m_strides[0] != 1,
                axom::fmt::format("Bump requires a structured topology whose i-stride is 1, "
                                  "but the topology has {}.",
                                  formatLayout(topoLayout)));

  if(!n_field.has_path("offsets") && !n_field.has_path("strides"))
  {
    // Without field layout metadata, topology node ids index the values directly.
    return;
  }

  Indexing fieldLayout(nodeDims);
  axom::bump::detail::setLayoutFromNode(n_field, fieldLayout);
  SLIC_ERROR_IF(
    fieldLayout.m_offsets != topoLayout.m_offsets || fieldLayout.m_strides != topoLayout.m_strides,
    axom::fmt::format(
      "Field '{}' has {}, but the topology node layout has {}. "
      "Bump kernels index this field by topology node id, so the layouts must match.",
      fieldName,
      formatLayout(fieldLayout),
      formatLayout(topoLayout)));
}

}  // end anonymous namespace

void validateVertexFieldIndexing(const conduit::Node& n_topology,
                                 const conduit::Node& n_field,
                                 const std::string& fieldName)
{
  const bool fieldHasLayout = n_field.has_path("offsets") || n_field.has_path("strides");
  if(!n_topology.has_path("elements/dims"))
  {
    SLIC_ERROR_IF(fieldHasLayout,
                  axom::fmt::format("Field '{}' has Blueprint offsets or strides, but its topology "
                                    "is not structured. Bump kernels use topology node ids and do "
                                    "not apply the field layout.",
                                    fieldName));
    return;
  }

  // Read the dimensions on the host, since the topology can be in device memory.
  const conduit::Node& n_dims = n_topology.fetch_existing("elements/dims");
  conduit::Node n_hostDims;
  axom::bump::utilities::copy<axom::SEQ_EXEC>(n_hostDims, n_dims);
  if(n_hostDims.has_child("k"))
  {
    validateStructuredVertexField<3>(n_dims, n_hostDims, n_field, fieldName);
  }
  else if(n_hostDims.has_child("j"))
  {
    validateStructuredVertexField<2>(n_dims, n_hostDims, n_field, fieldName);
  }
  else
  {
    validateStructuredVertexField<1>(n_dims, n_hostDims, n_field, fieldName);
  }
}

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom
