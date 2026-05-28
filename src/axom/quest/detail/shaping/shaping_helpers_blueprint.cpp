// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers_blueprint.hpp"

#if defined(AXOM_USE_CONDUIT)

#include "conduit_blueprint_mesh.hpp"

#if defined(AXOM_USE_BUMP)
  #include "GenerateQuadratureMesh.hpp"
  #include "axom/bump/views/dispatch_topology.hpp"
  #include "axom/bump/views/dispatch_unstructured_topology.hpp"
#endif

#include <vector>

namespace axom
{
namespace quest
{
namespace shaping
{

namespace
{

constexpr const char* QUADRATURE_COORDSET_NAME = "quadrature_points";
constexpr const char* QUADRATURE_TOPOLOGY_NAME = "quadrature_points";
constexpr const char* ORIGINAL_ELEMENTS_FIELD_NAME = "originalElements";
constexpr const char* QUADRATURE_WEIGHTS_FIELD_NAME = "quadratureWeights";
constexpr const char* QUADRATURE_PHYSICAL_WEIGHTS_FIELD_NAME = "quadraturePhysicalWeights";

numerics::QuadratureRule getBlueprintQuadratureRule(axom::numerics::QuadratureType quadratureType,
                                                    int npts,
                                                    int allocatorID)
{
  SLIC_ERROR_IF(npts < 1, axom::fmt::format("Invalid sample resolution {}.", npts));
  SLIC_ERROR_IF(
    !axom::numerics::is_supported_quadrature_type(quadratureType),
    axom::fmt::format("Quadrature type {} is not yet supported for Blueprint quadrature meshes.",
                      static_cast<int>(quadratureType)));

  return numerics::get_quadrature_rule(quadratureType, npts, allocatorID);
}
std::string getBlueprintCellShapeImpl(const conduit::Node& topoNode)
{
  const std::string topoType = topoNode.fetch_existing("type").as_string();
  if(topoNode.has_path("elements/shape"))
  {
    return topoNode.fetch_existing("elements/shape").as_string();
  }

  if(topoType == "structured")
  {
    const conduit::Node& dimsNode = topoNode.fetch_existing("elements/dims");
    if(dimsNode.has_child("k"))
    {
      return "hex";
    }
    if(dimsNode.has_child("j"))
    {
      return "quad";
    }
    if(dimsNode.has_child("i"))
    {
      return "line";
    }

    SLIC_ERROR("Structured Blueprint topology is missing recognizable element dims.");
  }

  SLIC_ERROR(axom::fmt::format("Blueprint topology type '{}' is missing 'elements/shape'.", topoType));
  return "";
}

#if defined(AXOM_USE_BUMP)
template <typename ExecSpace, typename CoordsetView>
void buildBlueprintQuadratureMesh(const conduit::Node& topoNode,
                                  const conduit::Node& coordsetNode,
                                  const CoordsetView& coordsetView,
                                  int allocatorID,
                                  const numerics::QuadratureRule& ruleX,
                                  const numerics::QuadratureRule& ruleY,
                                  const numerics::QuadratureRule& ruleZ,
                                  conduit::Node& meshNode)
{
  namespace views = axom::bump::views;
  constexpr int SupportedShapes = views::select_shapes(views::Quad_ShapeID, views::Hex_ShapeID);

  views::dispatch_topology<views::select_dimensions(2, 3), SupportedShapes>(
    topoNode,
    [&](const auto&, auto topoView) {
      GenerateQuadratureMesh<ExecSpace, decltype(topoView), CoordsetView> generator(topoView,
                                                                                    coordsetView);
      generator.setAllocatorID(allocatorID);
      generator.execute(topoNode,
                        coordsetNode,
                        QUADRATURE_TOPOLOGY_NAME,
                        QUADRATURE_COORDSET_NAME,
                        ORIGINAL_ELEMENTS_FIELD_NAME,
                        QUADRATURE_WEIGHTS_FIELD_NAME,
                        QUADRATURE_PHYSICAL_WEIGHTS_FIELD_NAME,
                        ruleX,
                        ruleY,
                        ruleZ,
                        meshNode);
    });
}
#endif

}  // namespace

std::string getBlueprintCellShape(const conduit::Node& topoNode)
{
  return getBlueprintCellShapeImpl(topoNode);
}

#if defined(AXOM_USE_BUMP)
void printRegisteredFieldNames(const BlueprintState& bpState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling AXOM_UNUSED_PARAM(vfSampling),
                               const std::string& initialMessage)
{
  auto extractChildren = [](const conduit::Node& node) {
    std::vector<std::string> names;
    if(node.dtype().is_object())
    {
      names.reserve(node.number_of_children());
      for(conduit::index_t i = 0; i < node.number_of_children(); ++i)
      {
        names.push_back(node.child(i).name());
      }
    }
    return names;
  };

  auto extractMatchingFields = [&](const std::string& prefix) {
    std::vector<std::string> names;
    if(bpState.m_internal_node.has_path("fields"))
    {
      const conduit::Node& fieldsNode = bpState.m_internal_node.fetch_existing("fields");
      for(conduit::index_t i = 0; i < fieldsNode.number_of_children(); ++i)
      {
        const std::string name = fieldsNode.child(i).name();
        if(axom::utilities::string::startsWith(name, prefix))
        {
          names.push_back(name);
        }
      }
    }
    return names;
  };

  auto extractOtherFields = [&]() {
    std::vector<std::string> names;
    if(bpState.m_internal_node.has_path("fields"))
    {
      const conduit::Node& fieldsNode = bpState.m_internal_node.fetch_existing("fields");
      for(conduit::index_t i = 0; i < fieldsNode.number_of_children(); ++i)
      {
        const std::string name = fieldsNode.child(i).name();
        if(!axom::utilities::string::startsWith(name, "inout_") &&
           !axom::utilities::string::startsWith(name, "mat_inout_") &&
           !axom::utilities::string::startsWith(name, "vol_frac_"))
        {
          names.push_back(name);
        }
      }
    }
    return names;
  };

  const std::vector<std::string> topologyNames = bpState.m_internal_node.has_path("topologies")
    ? extractChildren(bpState.m_internal_node.fetch_existing("topologies"))
    : std::vector<std::string> {};
  const std::vector<std::string> coordsetNames = bpState.m_internal_node.has_path("coordsets")
    ? extractChildren(bpState.m_internal_node.fetch_existing("coordsets"))
    : std::vector<std::string> {};
  const std::vector<std::string> fieldNames = bpState.m_internal_node.has_path("fields")
    ? extractChildren(bpState.m_internal_node.fetch_existing("fields"))
    : std::vector<std::string> {};

  axom::fmt::memory_buffer out;
  axom::fmt::format_to(std::back_inserter(out),
                       "List of registered fields in the SamplingShaper {}"
                       "\n\t* Blueprint topologies: {}"
                       "\n\t* Blueprint coordsets: {}"
                       "\n\t* Blueprint fields: {}"
                       "\n\t* Known materials: {}"
                       "\n\t* Shape inout fields: {}"
                       "\n\t* Mat inout fields: {}"
                       "\n\t* Volume fraction fields: {}"
                       "\n\t* Other Blueprint fields: {}",
                       initialMessage,
                       axom::fmt::join(topologyNames, ", "),
                       axom::fmt::join(coordsetNames, ", "),
                       axom::fmt::join(fieldNames, ", "),
                       axom::fmt::join(knownMaterials, ", "),
                       axom::fmt::join(extractMatchingFields("inout_"), ", "),
                       axom::fmt::join(extractMatchingFields("mat_inout_"), ", "),
                       axom::fmt::join(extractMatchingFields("vol_frac_"), ", "),
                       axom::fmt::join(extractOtherFields(), ", "));

  SLIC_INFO_ROOT(axom::fmt::to_string(out));
}

void generateQuadraturePointMesh(conduit::Node& bpMeshNode,
                                 const std::string& topologyName,
                                 int allocatorID,
                                 axom::ArrayView<int> sampleResolution,
                                 axom::numerics::QuadratureType quadratureType)
{
  if(bpMeshNode.has_path(axom::fmt::format("topologies/{}", QUADRATURE_TOPOLOGY_NAME)))
  {
    return;
  }

  const conduit::Node& topoNode =
    bpMeshNode.fetch_existing("topologies").fetch_existing(topologyName);
  const std::string topoType = topoNode.fetch_existing("type").as_string();
  SLIC_ERROR_IF(
    topoType != "unstructured" && topoType != "structured",
    axom::fmt::format("Unsupported Blueprint topology type '{}' for quadrature mesh generation.",
                      topoType));

  const std::string shape = shaping::getBlueprintCellShape(topoNode);
  SLIC_ERROR_IF(
    shape != "quad" && shape != "hex",
    axom::fmt::format("Unsupported Blueprint element shape '{}' for quadrature mesh generation.",
                      shape));

  const std::string coordsetName = topoNode.fetch_existing("coordset").as_string();
  const conduit::Node& coordsetNode =
    bpMeshNode.fetch_existing("coordsets").fetch_existing(coordsetName);
  const std::string coordsetType = coordsetNode.fetch_existing("type").as_string();
  SLIC_ERROR_IF(
    coordsetType != "explicit",
    axom::fmt::format("Unsupported Blueprint coordset type '{}' for quadrature mesh generation.",
                      coordsetType));

  int selectedAllocatorID = allocatorID;
  if(!axom::execution_space<seq_exec>::usesAllocId(selectedAllocatorID) &&
     !axom::execution_space<omp_exec>::usesAllocId(selectedAllocatorID)
  #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
     && !axom::execution_space<cuda_exec>::usesAllocId(selectedAllocatorID)
  #endif
  #if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
     && !axom::execution_space<hip_exec>::usesAllocId(selectedAllocatorID)
  #endif
  )
  {
    selectedAllocatorID = axom::execution_space<seq_exec>::allocatorID();
  }

  auto ruleX = getBlueprintQuadratureRule(quadratureType, sampleResolution[0], selectedAllocatorID);
  auto ruleY = getBlueprintQuadratureRule(quadratureType, sampleResolution[1], selectedAllocatorID);
  const int nz = (sampleResolution.size() > 2) ? sampleResolution[2] : 1;
  auto ruleZ = getBlueprintQuadratureRule(quadratureType, nz, selectedAllocatorID);

  axom::bump::views::dispatch_explicit_coordset(coordsetNode, [&](auto coordsetView) {
  #if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    if(axom::execution_space<hip_exec>::usesAllocId(selectedAllocatorID))
    {
      buildBlueprintQuadratureMesh<hip_exec>(topoNode,
                                             coordsetNode,
                                             coordsetView,
                                             selectedAllocatorID,
                                             ruleX,
                                             ruleY,
                                             ruleZ,
                                             bpMeshNode);
      return;
    }
  #endif
  #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    if(axom::execution_space<cuda_exec>::usesAllocId(selectedAllocatorID))
    {
      buildBlueprintQuadratureMesh<cuda_exec>(topoNode,
                                              coordsetNode,
                                              coordsetView,
                                              selectedAllocatorID,
                                              ruleX,
                                              ruleY,
                                              ruleZ,
                                              bpMeshNode);
      return;
    }
  #endif
    if(axom::execution_space<omp_exec>::usesAllocId(selectedAllocatorID))
    {
      buildBlueprintQuadratureMesh<omp_exec>(topoNode,
                                             coordsetNode,
                                             coordsetView,
                                             selectedAllocatorID,
                                             ruleX,
                                             ruleY,
                                             ruleZ,
                                             bpMeshNode);
      return;
    }

    buildBlueprintQuadratureMesh<seq_exec>(topoNode,
                                           coordsetNode,
                                           coordsetView,
                                           selectedAllocatorID,
                                           ruleX,
                                           ruleY,
                                           ruleZ,
                                           bpMeshNode);
  });
}

void generateSamplingPositions(BlueprintState& bpState,
                               axom::ArrayView<int> sampleResolution,
                               axom::numerics::QuadratureType quadratureType)
{
  AXOM_ANNOTATE_SCOPE("generateSamplingPositions");
  checkSampleResolution(bpState, sampleResolution, quadratureType);

  if(bpState.m_internal_node.has_path(axom::fmt::format("topologies/{}", QUADRATURE_TOPOLOGY_NAME)))
  {
    return;
  }

  generateQuadraturePointMesh(bpState.m_internal_node,
                              bpState.m_topology_name,
                              bpState.m_allocator_id,
                              sampleResolution,
                              quadratureType);
}

void importInitialVolumeFractions(BlueprintState& bpState,
                                  const std::map<std::string, conduit::Node*>& initialVolumeFractions)
{
  conduit::Node& n_mesh = bpState.getBlueprintMeshNode();
  const std::string quadName("quadrature_points");
  const conduit::Node& n_quad_points = n_mesh.fetch_existing("coordsets/" + quadName);
  const auto totalQuadPoints = conduit::blueprint::mesh::coordset::length(n_quad_points);

  // Get the topology we want to sample.
  const conduit::Node& n_topo = bpState.getBlueprintTopologyNode();
  const auto totalZones = conduit::blueprint::mesh::topology::length(n_topo);

  const auto samplesPerZone = totalQuadPoints / totalZones;

  for(auto& entry : initialVolumeFractions)
  {
    const auto& name = entry.first;
    auto* field_ptr = entry.second;

    SLIC_INFO_ROOT(axom::fmt::format("Importing volume fraction field for '{}' material", name));

    if(field_ptr == nullptr)
    {
      SLIC_WARNING(
        axom::fmt::format("Skipping missing volume fraction field for material '{}'", name));
      continue;
    }

    // Get the source field.
    const auto srcPath = axom::fmt::format("fields/vol_frac_{}", name);
    conduit::Node& n_src_field = n_mesh.fetch_existing(srcPath);
    SLIC_ERROR_IF(n_src_field.fetch_existing("association").as_string() != "element",
                  "The imported field must have element association.");
    const auto src_values = n_src_field["values"].as_double_accessor();

    // Make the new quadrature field.
    const auto destPath = axom::fmt::format("fields/mat_inout_{}", name);
    conduit::Node& n_dest_field = n_mesh.fetch(destPath);
    n_dest_field["topology"] = quadName;
    n_dest_field["association"] = "element";
    conduit::Node& n_dest_values = n_dest_field["values"];
    n_dest_values.set(conduit::DataType::float64(totalQuadPoints));
    double* dptr = n_dest_values.as_double_ptr();

    // Copy the source field into the dest field. We just copy samplesPerZone values
    // from the source into the dest since each block of samplesPerZone points in
    // the quadrature mesh corresponds to a zone in the source mesh.
    for(conduit::index_t i = 0; i < totalZones; i++)
    {
      const auto src_value = src_values[i];
      for(conduit::index_t c = 0; c < samplesPerZone; c++)
      {
        *dptr++ = src_value;
      }
    }
  }
}

void computeVolumeFractionsForMaterial(BlueprintState& bpState, const std::string& matField)
{
  AXOM_ANNOTATE_SCOPE("computeVolumeFractionsForMaterial");

  SLIC_ASSERT(axom::utilities::string::startsWith(matField, "mat_inout_"));

  conduit::Node* inout = bpState.getMaterialFunction(matField);
  SLIC_ERROR_IF(
    inout == nullptr,
    axom::fmt::format("Missing Blueprint material field '{}' for volume fraction projection.",
                      matField));

  conduit::Node& bpMeshNode = bpState.m_internal_node;
  SLIC_ERROR_IF(!bpMeshNode.has_path("fields/originalElements/values"),
                "Missing Blueprint originalElements field for volume fraction projection.");
  SLIC_ERROR_IF(!bpMeshNode.has_path("fields/quadraturePhysicalWeights/values") &&
                  !bpMeshNode.has_path("fields/quadratureWeights/values"),
                "Missing Blueprint quadrature weight field for volume fraction projection.");

  const conduit::Node& topoNode =
    bpMeshNode.fetch_existing("topologies").fetch_existing(bpState.m_topology_name);

  const axom::IndexType numZones = conduit::blueprint::mesh::topology::length(topoNode);

  namespace utils = axom::bump::utilities;
  const auto originalElements =
    utils::make_array_view<conduit::index_t>(bpMeshNode["fields/originalElements/values"]);
  const conduit::Node& quadratureWeightsNode =
    bpMeshNode.has_path("fields/quadraturePhysicalWeights/values")
    ? bpMeshNode["fields/quadraturePhysicalWeights/values"]
    : bpMeshNode["fields/quadratureWeights/values"];
  const auto quadratureWeights = utils::make_array_view<double>(quadratureWeightsNode);
  const auto inoutValues = utils::make_array_view<double>(inout->fetch_existing("values"));

  SLIC_ASSERT(originalElements.size() == quadratureWeights.size());
  SLIC_ASSERT(originalElements.size() == inoutValues.size());

  const std::string vfName = axom::fmt::format("vol_frac_{}", matField.substr(10));
  conduit::Node& vfNode = bpMeshNode["fields/" + vfName];
  vfNode.reset();
  vfNode["association"] = "element";
  vfNode["topology"] = bpState.m_topology_name;

  const auto conduitAllocatorId =
    axom::sidre::ConduitMemory::axomAllocIdToConduit(bpState.m_allocator_id);
  conduit::Node& valuesNode = vfNode["values"];
  valuesNode.set_allocator(conduitAllocatorId);
  valuesNode.set(conduit::DataType::float64(numZones));
  auto vfValues = utils::make_array_view<double>(valuesNode);
  axom::Array<double> totalWeights(numZones, numZones, bpState.m_allocator_id);
  auto totalWeightsView = totalWeights.view();

  for(axom::IndexType zoneIdx = 0; zoneIdx < vfValues.size(); ++zoneIdx)
  {
    vfValues[zoneIdx] = 0.;
    totalWeightsView[zoneIdx] = 0.;
  }

  for(axom::IndexType pointIdx = 0; pointIdx < inoutValues.size(); ++pointIdx)
  {
    const conduit::index_t zoneIdx = originalElements[pointIdx];
    SLIC_ASSERT(zoneIdx >= 0);
    SLIC_ASSERT(zoneIdx < vfValues.size());
    vfValues[zoneIdx] += inoutValues[pointIdx] * quadratureWeights[pointIdx];
    totalWeightsView[zoneIdx] += quadratureWeights[pointIdx];
  }

  for(axom::IndexType zoneIdx = 0; zoneIdx < vfValues.size(); ++zoneIdx)
  {
    SLIC_ERROR_IF(
      axom::utilities::isNearlyEqual(totalWeightsView[zoneIdx], 0.0),
      axom::fmt::format(
        "Blueprint quadrature weights sum to zero in zone {} during volume fraction projection.",
        zoneIdx));
    vfValues[zoneIdx] /= totalWeightsView[zoneIdx];
  }
}

void replaceMaterial(conduit::Node* shapeNode, conduit::Node* materialNode, bool shapeReplacesMaterial)
{
  SLIC_ASSERT(shapeNode != nullptr);
  SLIC_ASSERT(materialNode != nullptr);

  namespace utils = axom::bump::utilities;
  auto shapeValues = utils::make_array_view<double>(shapeNode->fetch_existing("values"));
  auto materialValues = utils::make_array_view<double>(materialNode->fetch_existing("values"));

  SLIC_ASSERT(shapeValues.size() == materialValues.size());

  for(axom::IndexType i = 0; i < materialValues.size(); ++i)
  {
    if(shapeReplacesMaterial)
    {
      materialValues[i] = shapeValues[i] > 0. ? 0. : materialValues[i];
    }
    else
    {
      shapeValues[i] = materialValues[i] > 0. ? 0. : shapeValues[i];
    }
  }
}

void copyShapeIntoMaterial(const conduit::Node* shapeNode,
                           conduit::Node* materialNode,
                           bool reuseExisting)
{
  SLIC_ASSERT(shapeNode != nullptr);
  SLIC_ASSERT(materialNode != nullptr);

  namespace utils = axom::bump::utilities;
  const auto shapeValues = utils::make_array_view<double>(shapeNode->fetch_existing("values"));
  auto materialValues = utils::make_array_view<double>(materialNode->fetch_existing("values"));

  SLIC_ASSERT(shapeValues.size() == materialValues.size());

  if(reuseExisting)
  {
    for(axom::IndexType i = 0; i < materialValues.size(); ++i)
    {
      materialValues[i] = shapeValues[i] > 0. ? 1. : materialValues[i];
    }
  }
  else
  {
    for(axom::IndexType i = 0; i < materialValues.size(); ++i)
    {
      materialValues[i] = shapeValues[i];
    }
  }
}

conduit::Node* cloneInOutFunction(const conduit::Node* node)
{
  SLIC_ASSERT(node != nullptr);
  return new conduit::Node(*node);
}

#endif  // defined(AXOM_USE_BUMP)

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // defined(AXOM_USE_CONDUIT)
