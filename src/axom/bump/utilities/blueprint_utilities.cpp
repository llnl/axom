// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/IndexingPolicies.hpp"

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include <conduit/conduit_blueprint.hpp>

#include <string>
#include <vector>

namespace axom::bump::utilities
{

const conduit::Node* findReferenceNode(const conduit::Node& node,
                                       const std::string& refKey)
{
  SLIC_ERROR_IF(refKey.empty(), "Reference key must not be empty.");

  const conduit::Node* result = nullptr;

  const std::string refSection =
    (refKey.back() != 'y') ? refKey + "s" : refKey.substr(0, refKey.length() - 1) + "ies";

  if(node.has_child(refKey))
  {
    const std::string refValue = hostReadString(node.fetch_existing(refKey));

    const conduit::Node* traverseNode = node.parent();
    while(traverseNode != nullptr)
    {
      if(traverseNode->has_child(refSection))
      {
        const conduit::Node& refParent = traverseNode->fetch_existing(refSection);
        if(refParent.has_child(refValue))
        {
          result = &refParent.fetch_existing(refValue);
        }
        break;
      }

      traverseNode = traverseNode->parent();
    }
  }

  return result;
}

conduit::Node* findReferenceNode(conduit::Node& node, const std::string& refKey)
{
  return const_cast<conduit::Node*>(findReferenceNode(static_cast<const conduit::Node&>(node),
                                                      refKey));
}

std::vector<std::string> coordsetAxes(const conduit::Node& n_coordset)
{
  std::vector<std::string> axes;

  // Get the axis names for the output coordset. For uniform, prefer x,y,z
  // instead of i,j,k since we're making an explicit coordset.
  if(coordsetType(n_coordset) == "uniform")
  {
    if(n_coordset.has_path("dims/i")) axes.push_back("x");
    if(n_coordset.has_path("dims/j")) axes.push_back("y");
    if(n_coordset.has_path("dims/k")) axes.push_back("z");
  }
  else
  {
    if(isDeviceAllocated(n_coordset.fetch_existing("type")))
    {
      const conduit::Node &n_values = n_coordset.fetch_existing("values");
      for(conduit::index_t i = 0; i < n_values.number_of_children(); i++)
      {
        axes.push_back(n_values[i].name());
      }
    }
    else
    {
      axes = conduit::blueprint::mesh::utils::coordset::axes(n_coordset);
    }
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
void validateStructuredVertexField(const conduit::Node& n_hostDims,
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
  axom::bump::detail::setLayoutFromNode(n_hostDims, topoLayout);

  // StructuredTopologyView::zone() constructs adjacent corner ids by adding 1.
  SLIC_ERROR_IF(topoLayout.m_strides[0] != 1,
                axom::fmt::format("Bump requires a structured topology whose i-stride is 1, "
                                  "but the topology has {}.",
                                  formatLayout(topoLayout)));

  if(!n_field.has_path("offsets") && !n_field.has_path("strides"))
  {
    // Topology node ids directly index fields without layout metadata.
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

  // Copy dimensions to the host because the topology may reside in device memory.
  const conduit::Node& n_dims = n_topology.fetch_existing("elements/dims");
  conduit::Node n_host;
  conduit::Node& n_hostDims = n_host[n_dims.path()];
  axom::bump::utilities::copy<axom::SEQ_EXEC>(n_hostDims, n_dims);
  if(n_hostDims.has_child("k"))
  {
    validateStructuredVertexField<3>(n_hostDims, n_field, fieldName);
  }
  else if(n_hostDims.has_child("j"))
  {
    validateStructuredVertexField<2>(n_hostDims, n_field, fieldName);
  }
  else
  {
    validateStructuredVertexField<1>(n_hostDims, n_field, fieldName);
  }
}

std::string coordsetType(const conduit::Node &n_coordset)
{
  SLIC_ERROR_IF(!n_coordset.has_path("type"), "The coordset must contain \"type\".");
  return hostReadString(n_coordset["type"]);
}

int coordsetDimension(const conduit::Node &n_coordset)
{
  const auto axes = coordsetAxes(n_coordset);
  return static_cast<int>(axes.size());
}

std::string topologyType(const conduit::Node &n_topo)
{
  SLIC_ERROR_IF(!n_topo.has_path("type"), "The topology must contain \"type\".");
  return hostReadString(n_topo["type"]);
}

std::string topologyCoordset(const conduit::Node &n_topo)
{
  SLIC_ERROR_IF(!n_topo.has_path("coordset"), "The topology must contain \"coordset\".");
  return hostReadString(n_topo["coordset"]);
}

int topologyDimension(const conduit::Node &n_topo)
{
  int dims = 0;
  const auto type = topologyType(n_topo);
  if(type != "unstructured")
  {
    const conduit::Node *n_coordset = findReferenceNode(n_topo, "coordset");
    SLIC_ERROR_IF(n_coordset == nullptr, "Coordset not found.");
    dims = coordsetDimension(*n_coordset);
  }
  else
  {
    const auto shape = unstructuredTopologyShape(n_topo);
    if(shape == "mixed")
    {
      SLIC_ERROR_IF(!n_topo.has_path("elements/shape_map"), "Mixed topologies need shape_map");
      // Iterate over the shape names in the shape map and use ShapeType to
      // determine the dimension.
      const conduit::Node& n_shape_map = n_topo["elements/shape_map"];
      for(conduit::index_t i = 0; i < n_shape_map.number_of_children(); i++)
      {
        conduit::blueprint::mesh::utils::ShapeType S(n_shape_map[i].name());
        dims = axom::utilities::max(dims, static_cast<int>(S.dim));
      }
    }
    else
    {
      conduit::blueprint::mesh::utils::ShapeType S(shape);
      dims = static_cast<int>(S.dim);
    }
  }

  return dims;
}

conduit::blueprint::mesh::utils::ShapeType topologyShapeType(const conduit::Node &n_topo)
{
  // NOTE: This function is based on ShapeType::ShapeType(const std::string &).
  conduit::blueprint::mesh::utils::ShapeType S;

  const auto type = topologyType(n_topo);
  const int dim = topologyDimension(n_topo);
  if(type == "unstructured")
  {
    const auto shape = unstructuredTopologyShape(n_topo);
    S = conduit::blueprint::mesh::utils::ShapeType(shape);
    if(shape == "mixed")
    {
      S.dim = dim;
    }
  }
  else if(type == "points")
  {
    S = conduit::blueprint::mesh::utils::ShapeType("point");
  }
  else
  {
    const char *dimShapes[] = {"point", "line", "quad", "hex"};
    S = conduit::blueprint::mesh::utils::ShapeType(dimShapes[dim]);
  }

  return S;
}

std::string unstructuredTopologyShape(const conduit::Node &n_topo)
{
  SLIC_ERROR_IF(!n_topo.has_path("elements/shape"), "The topology must contain \"elements/shape\".");
  return hostReadString(n_topo["elements/shape"]);
}

std::string matsetTopology(const conduit::Node &n_matset)
{
  SLIC_ERROR_IF(!n_matset.has_path("topology"), "The matset must contain \"topology\".");
  return hostReadString(n_matset["topology"]);
}

std::string fieldAssociation(const conduit::Node &n_field)
{
  SLIC_ERROR_IF(!n_field.has_path("association"), "The field must contain \"association\".");
  return hostReadString(n_field["association"]);
}

std::string fieldTopology(const conduit::Node &n_field)
{
  SLIC_ERROR_IF(!n_field.has_path("topology"), "The field must contain \"topology\".");
  return hostReadString(n_field["topology"]);
}

}  // end namespace axom::bump::utilities
