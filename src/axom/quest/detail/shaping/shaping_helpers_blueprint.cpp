// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers_blueprint.hpp"

#if defined(AXOM_USE_CONDUIT)

  #include "axom/quest/util/mesh_helpers.hpp"
  #include "conduit_blueprint_mesh.hpp"

  #if defined(AXOM_USE_BUMP)
    #include "axom/bump/GenerateQuadratureMesh.hpp"
    #include "axom/bump/utilities/conduit_memory.hpp"
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
  constexpr int SupportedShapes = (CoordsetView::dimension() == 2)
    ? views::select_shapes(views::Quad_ShapeID)
    : views::select_shapes(views::Hex_ShapeID);

  views::dispatch_topology<views::select_dimensions(CoordsetView::dimension()), SupportedShapes>(
    topoNode,
    [&](const auto&, auto topoView) {
      axom::bump::GenerateQuadratureMesh<ExecSpace, decltype(topoView), CoordsetView> generator(
        topoView,
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

void BlueprintState::refreshBlueprintMeshNode()
{
  if(isSidreBacked())
  {
    m_internal_node.reset();
    m_group_ptr->createNativeLayout(m_internal_node);
  }
}

int BlueprintState::meshDimension() const
{
  const std::string shapeType = cellShape();

  if(shapeType == "quad")
  {
    return 2;
  }
  if(shapeType == "hex")
  {
    return 3;
  }

  SLIC_ERROR(axom::fmt::format("Unsupported Blueprint cell shape '{}'.", shapeType));
  return -1;
}

std::vector<std::string> BlueprintState::childNames(const std::string& path) const
{
  std::vector<std::string> names;
  const conduit::Node& bpMeshNode = getBlueprintMeshNode();
  if(!bpMeshNode.has_path(path))
  {
    return names;
  }

  const conduit::Node& node = bpMeshNode.fetch_existing(path);
  if(!node.dtype().is_object())
  {
    return names;
  }

  names.reserve(node.number_of_children());
  for(conduit::index_t i = 0; i < node.number_of_children(); ++i)
  {
    names.push_back(node.child(i).name());
  }
  return names;
}

void BlueprintState::ensureUnstructured(axom::runtime_policy::Policy execPolicy)
{
  const std::string topoType = getBlueprintTopologyNode().fetch_existing("type").as_string();
  if(topoType == "unstructured")
  {
    return;
  }

  if(isSidreBacked())
  {
    // Sidre
    const std::string shapeType = shaping::getBlueprintCellShape(getBlueprintTopologyNode());
    if(shapeType == "hex")
    {
      axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_3d(m_group_ptr,
                                                                                  topologyName(),
                                                                                  execPolicy);
    }
    else if(shapeType == "quad")
    {
      axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_2d(m_group_ptr,
                                                                                  topologyName(),
                                                                                  execPolicy);
    }
    else
    {
      SLIC_ERROR("Axom Internal error: Unhandled shape type.");
    }

    refreshBlueprintMeshNode();
  }
  else
  {
    // Conduit
    axom::quest::util::convert_blueprint_structured_explicit_to_unstructured(getBlueprintMeshNode(),
                                                                             topologyName(),
                                                                             execPolicy);
  }
}

void BlueprintState::deleteShapeFunction(const std::string& name)
{
  if(isSidreBacked())
  {
    const std::string fieldPath = axom::fmt::format("fields/{}", name);
    if(m_group_ptr->hasGroup(fieldPath))
    {
      m_group_ptr->destroyGroupAndData(fieldPath);
      refreshBlueprintMeshNode();
    }
    return;
  }

  conduit::Node& bpMeshNode = getBlueprintMeshNode();
  if(bpMeshNode.has_path("fields"))
  {
    conduit::Node& n_fields = bpMeshNode["fields"];
    if(n_fields.has_path(name))
    {
      n_fields.remove(name);
    }
  }
}

axom::ArrayView<double> BlueprintState::createField(const std::string& name,
                                                    const std::string& topologyName,
                                                    axom::IndexType size,
                                                    bool addVolumeDependent,
                                                    bool volumeDependent)
{
  if(isSidreBacked())
  {
    const std::string fieldPath = axom::fmt::format("fields/{}", name);
    if(m_group_ptr->hasGroup(fieldPath))
    {
      m_group_ptr->destroyGroupAndData(fieldPath);
    }

    auto* fieldGrp = m_group_ptr->createGroup(fieldPath);
    SLIC_ASSERT(fieldGrp != nullptr);
    fieldGrp->createViewString("association", "element");
    fieldGrp->createViewString("topology", topologyName);
    if(addVolumeDependent)
    {
      fieldGrp->createViewString("volume_dependent", volumeDependent ? "true" : "false");
    }

    auto* valuesView =
      fieldGrp->createViewAndAllocate("values", axom::sidre::DataTypeId::FLOAT64_ID, size);
    SLIC_ASSERT(valuesView != nullptr);
    refreshBlueprintMeshNode();
    return axom::ArrayView<double>(static_cast<double*>(valuesView->getVoidPtr()), size);
  }

  conduit::Node& fieldNode = getBlueprintMeshNode()["fields/" + name];
  fieldNode.reset();
  fieldNode["association"] = "element";
  fieldNode["topology"] = topologyName;
  if(addVolumeDependent)
  {
    fieldNode["volume_dependent"] = volumeDependent ? "true" : "false";
  }

  const auto conduitAllocatorId = axom::sidre::ConduitMemory::axomAllocIdToConduit(m_allocator_id);
  conduit::Node& valuesNode = fieldNode["values"];
  valuesNode.set_allocator(conduitAllocatorId);
  valuesNode.set(conduit::DataType::float64(size));
  return axom::ArrayView<double>(valuesNode.as_double_ptr(), valuesNode.dtype().number_of_elements());
}

axom::ArrayView<double> BlueprintState::getScalarFieldView(const std::string& name,
                                                           axom::IndexType size)
{
  conduit::Node& fieldNode = getField(name);
  SLIC_ASSERT(fieldNode.fetch_existing("association").as_string() == std::string("element"));
  SLIC_ASSERT(fieldNode.fetch_existing("topology").as_string() == topologyName());

  conduit::Node& valuesNode = fieldNode.fetch_existing("values");
  SLIC_ASSERT(valuesNode.dtype().id() == conduit::DataType::float64(size).id());
  SLIC_ASSERT(valuesNode.dtype().number_of_elements() == size);

  return axom::ArrayView<double>(valuesNode.as_double_ptr(), size);
}

  #if defined(AXOM_USE_BUMP)
void BlueprintState::importQuadraturePointMesh(const conduit::Node& quadratureMesh)
{
  auto replaceSubtree = [&](const std::string& path, const conduit::Node& node) {
    if(isSidreBacked())
    {
      if(m_group_ptr->hasGroup(path))
      {
        m_group_ptr->destroyGroupAndData(path);
      }

      auto* group = m_group_ptr->createGroup(path);
      SLIC_ERROR_IF(group == nullptr,
                    axom::fmt::format("Failed to create Sidre group for Blueprint path '{}'.", path));
      const bool importSuccess = group->importConduitTree(node);
      SLIC_ERROR_IF(!importSuccess,
                    axom::fmt::format("Failed to import Blueprint subtree '{}'.", path));
      return;
    }

    conduit::Node& outputNode = getBlueprintMeshNode()[path];
    outputNode.reset();
    outputNode.update(node);
  };

  const std::string quadratureCoordsetPath =
    axom::fmt::format("coordsets/{}", QUADRATURE_COORDSET_NAME);
  const std::string quadratureTopologyPath =
    axom::fmt::format("topologies/{}", QUADRATURE_TOPOLOGY_NAME);
  const std::string originalElementsPath =
    axom::fmt::format("fields/{}", ORIGINAL_ELEMENTS_FIELD_NAME);
  const std::string quadratureWeightsPath =
    axom::fmt::format("fields/{}", QUADRATURE_WEIGHTS_FIELD_NAME);
  const std::string quadraturePhysicalWeightsPath =
    axom::fmt::format("fields/{}", QUADRATURE_PHYSICAL_WEIGHTS_FIELD_NAME);

  SLIC_ERROR_IF(!quadratureMesh.has_path(quadratureCoordsetPath),
                "Quadrature mesh is missing its Blueprint quadrature coordset.");
  SLIC_ERROR_IF(!quadratureMesh.has_path(quadratureTopologyPath),
                "Quadrature mesh is missing its Blueprint quadrature topology.");

  replaceSubtree(quadratureCoordsetPath, quadratureMesh.fetch_existing(quadratureCoordsetPath));
  replaceSubtree(quadratureTopologyPath, quadratureMesh.fetch_existing(quadratureTopologyPath));
  replaceSubtree(originalElementsPath, quadratureMesh.fetch_existing(originalElementsPath));
  replaceSubtree(quadratureWeightsPath, quadratureMesh.fetch_existing(quadratureWeightsPath));

  if(quadratureMesh.has_path(quadraturePhysicalWeightsPath))
  {
    replaceSubtree(quadraturePhysicalWeightsPath,
                   quadratureMesh.fetch_existing(quadraturePhysicalWeightsPath));
  }

  if(isSidreBacked())
  {
    refreshBlueprintMeshNode();
  }
}

conduit::Node* BlueprintState::createMaterialFunction(const std::string& name)
{
  SLIC_ERROR_IF(
    !getBlueprintMeshNode().has_path(
      axom::fmt::format("coordsets/{}/values", QUADRATURE_COORDSET_NAME)),
    std::string("Cannot create material function '") + name + "' without quadrature points.");

  const conduit::Node& values =
    getBlueprintCoordsetNode(QUADRATURE_COORDSET_NAME).fetch_existing("values");
  const auto numValues = values.child(0).dtype().number_of_elements();

  auto fieldValues = createField(name, QUADRATURE_TOPOLOGY_NAME, numValues);
  for(axom::IndexType i = 0; i < fieldValues.size(); ++i)
  {
    fieldValues[i] = 0.;
  }

  return &getField(name);
}

void printRegisteredFieldNames(const BlueprintState& bpState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling AXOM_UNUSED_PARAM(vfSampling),
                               const std::string& initialMessage)
{
  auto extractMatchingFields = [&](const std::string& prefix) {
    std::vector<std::string> names;
    for(const auto& name : bpState.fieldNames())
    {
      if(axom::utilities::string::startsWith(name, prefix))
      {
        names.push_back(name);
      }
    }
    return names;
  };

  auto extractOtherFields = [&]() {
    std::vector<std::string> names;
    for(const auto& name : bpState.fieldNames())
    {
      if(!shaping::isVolumeFractionFieldName(name) && !shaping::isMaterialInOutFieldName(name) &&
         !shaping::isShapeInOutFieldName(name))
      {
        names.push_back(name);
      }
    }
    return names;
  };

  const std::vector<std::string> topologyNames = bpState.topologyNames();
  const std::vector<std::string> coordsetNames = bpState.coordsetNames();
  const std::vector<std::string> fieldNames = bpState.fieldNames();

  axom::fmt::memory_buffer out;
  axom::fmt::format_to(
    std::back_inserter(out),
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
    axom::fmt::join(extractMatchingFields(shaping::shapeInOutFieldName("")), ", "),
    axom::fmt::join(extractMatchingFields(shaping::materialInOutFieldName("")), ", "),
    axom::fmt::join(extractMatchingFields(shaping::volumeFractionFieldName("")), ", "),
    axom::fmt::join(extractOtherFields(), ", "));

  SLIC_INFO_ROOT(axom::fmt::to_string(out));
}

void generateQuadraturePointMesh(const conduit::Node& bpMeshNode,
                                 conduit::Node& outputMeshNode,
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
                                             outputMeshNode);
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
                                              outputMeshNode);
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
                                             outputMeshNode);
      return;
    }

    buildBlueprintQuadratureMesh<seq_exec>(topoNode,
                                           coordsetNode,
                                           coordsetView,
                                           selectedAllocatorID,
                                           ruleX,
                                           ruleY,
                                           ruleZ,
                                           outputMeshNode);
  });
}

void generateSamplingPositions(BlueprintState& bpState,
                               axom::ArrayView<int> sampleResolution,
                               axom::numerics::QuadratureType quadratureType)
{
  AXOM_ANNOTATE_SCOPE("generateSamplingPositions");
  checkSampleResolution(bpState, sampleResolution, quadratureType);

  conduit::Node& bpMeshNode = bpState.getBlueprintMeshNode();
  if(bpMeshNode.has_path(axom::fmt::format("topologies/{}", QUADRATURE_TOPOLOGY_NAME)))
  {
    return;
  }

  if(bpState.isSidreBacked())
  {
    conduit::Node quadratureMesh;
    generateQuadraturePointMesh(bpMeshNode,
                                quadratureMesh,
                                bpState.topologyName(),
                                bpState.allocatorId(),
                                sampleResolution,
                                quadratureType);
    bpState.importQuadraturePointMesh(quadratureMesh);
  }
  else
  {
    generateQuadraturePointMesh(bpMeshNode,
                                bpMeshNode,
                                bpState.topologyName(),
                                bpState.allocatorId(),
                                sampleResolution,
                                quadratureType);
  }
}

void importInitialVolumeFractions(BlueprintState& bpState,
                                  const std::map<std::string, conduit::Node*>& initialVolumeFractions)
{
  conduit::Node& n_mesh = bpState.getBlueprintMeshNode();
  const conduit::Node& n_quad_points =
    n_mesh.fetch_existing(axom::fmt::format("coordsets/{}", QUADRATURE_COORDSET_NAME));
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
    const auto srcPath = axom::fmt::format("fields/{}", shaping::volumeFractionFieldName(name));
    conduit::Node& n_src_field = n_mesh.fetch_existing(srcPath);
    SLIC_ERROR_IF(n_src_field.fetch_existing("association").as_string() != "element",
                  "The imported field must have element association.");
    const auto src_values = n_src_field["values"].as_double_accessor();

    // Make the new quadrature field.
    auto destValues = bpState.createField(shaping::materialInOutFieldName(name),
                                          QUADRATURE_TOPOLOGY_NAME,
                                          totalQuadPoints);
    double* dptr = destValues.data();

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

  const std::string materialName = shaping::materialNameFromMaterialInOutFieldName(matField);
  SLIC_ASSERT(!materialName.empty());

  conduit::Node* inout = bpState.getMaterialFunction(matField);
  SLIC_ERROR_IF(
    inout == nullptr,
    axom::fmt::format("Missing Blueprint material field '{}' for volume fraction projection.",
                      matField));

  conduit::Node& bpMeshNode = bpState.getBlueprintMeshNode();
  const std::string originalElementsPath =
    axom::fmt::format("fields/{}/values", ORIGINAL_ELEMENTS_FIELD_NAME);
  const std::string quadraturePhysicalWeightsPath =
    axom::fmt::format("fields/{}/values", QUADRATURE_PHYSICAL_WEIGHTS_FIELD_NAME);
  const std::string quadratureWeightsPath =
    axom::fmt::format("fields/{}/values", QUADRATURE_WEIGHTS_FIELD_NAME);

  SLIC_ERROR_IF(!bpMeshNode.has_path(originalElementsPath),
                "Missing Blueprint originalElements field for volume fraction projection.");
  SLIC_ERROR_IF(!bpMeshNode.has_path(quadraturePhysicalWeightsPath) &&
                  !bpMeshNode.has_path(quadratureWeightsPath),
                "Missing Blueprint quadrature weight field for volume fraction projection.");

  const conduit::Node& topoNode = bpState.getBlueprintTopologyNode();

  const axom::IndexType numZones = conduit::blueprint::mesh::topology::length(topoNode);

  namespace utils = axom::bump::utilities;
  const auto originalElements =
    utils::make_array_view<conduit::index_t>(bpMeshNode.fetch_existing(originalElementsPath));
  const conduit::Node& quadratureWeightsNode = bpMeshNode.has_path(quadraturePhysicalWeightsPath)
    ? bpMeshNode.fetch_existing(quadraturePhysicalWeightsPath)
    : bpMeshNode.fetch_existing(quadratureWeightsPath);
  const auto quadratureWeights = utils::make_array_view<double>(quadratureWeightsNode);
  const auto inoutValues = utils::make_array_view<double>(inout->fetch_existing("values"));

  SLIC_ASSERT(originalElements.size() == quadratureWeights.size());
  SLIC_ASSERT(originalElements.size() == inoutValues.size());

  const std::string vfName = shaping::volumeFractionFieldName(materialName);
  auto vfValues = bpState.createField(vfName, bpState.topologyName(), numZones);
  axom::Array<double> totalWeights(numZones, numZones, bpState.allocatorId());
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
