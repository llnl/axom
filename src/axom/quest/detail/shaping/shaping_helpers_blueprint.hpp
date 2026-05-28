// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SHAPING_HELPERS_BLUEPRINT__HPP_
#define AXOM_QUEST_SHAPING_HELPERS_BLUEPRINT__HPP_

#include "shaping_helpers.hpp"

#if defined(AXOM_USE_CONDUIT)

#include "axom/fmt.hpp"

#if defined(AXOM_USE_BUMP)
  #include "axom/bump/utilities/conduit_memory.hpp"
  #include "axom/bump/views/dispatch_coordset.hpp"
#endif

#include "conduit_node.hpp"

#include <set>
#include <string>
#include <vector>

namespace axom
{
namespace quest
{
namespace shaping
{

std::string getBlueprintCellShape(const conduit::Node& topoNode);

struct BlueprintState
{
  virtual ~BlueprintState() = default;

  axom::sidre::Group* m_group_ptr {nullptr};
  int m_allocator_id {axom::getDefaultAllocatorID()};
  std::string m_topology_name;
  conduit::Node* m_external_node_ptr {nullptr};
  conduit::Node m_internal_node;

  int meshDimension() const
  {
    const std::string shapeType = shaping::getBlueprintCellShape(getBlueprintTopologyNode());

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

  const conduit::Node& getBlueprintTopologyNode() const
  {
    return m_internal_node.fetch_existing("topologies").fetch_existing(m_topology_name);
  }

  conduit::Node* getShapeFunction(const std::string& name)
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

  const conduit::Node* getShapeFunction(const std::string& name) const
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

  void deleteShapeFunction(const std::string& name)
  {
    if(m_internal_node.has_path("fields"))
    {
      conduit::Node& n_fields = m_internal_node["fields"];
      if(n_fields.has_path(name))
      {
        n_fields.remove(name);
      }
    }
  }

  conduit::Node* getMaterialFunction(const std::string& name)
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

  const conduit::Node* getMaterialFunction(const std::string& name) const
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

#if defined(AXOM_USE_BUMP)
  conduit::Node* createMaterialFunction(const std::string& name)
  {
    constexpr const char* quadratureTopologyName = "quadrature_points";
    SLIC_ERROR_IF(!m_internal_node.has_path("coordsets/quadrature_points/values"),
                  std::string("Cannot create material function '") + name +
                    "' without quadrature points.");

    conduit::Node& fieldNode = m_internal_node["fields/" + name];
    fieldNode.reset();
    fieldNode["association"] = "element";
    fieldNode["topology"] = quadratureTopologyName;

    const auto conduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(m_allocator_id);
    conduit::Node& valuesNode = fieldNode["values"];
    valuesNode.set_allocator(conduitAllocatorId);

    const conduit::Node& values =
      m_internal_node["coordsets/quadrature_points"].fetch_existing("values");
    const auto numValues = values.child(0).dtype().number_of_elements();
    valuesNode.set(conduit::DataType::float64(numValues));

    auto fieldValues = axom::bump::utilities::make_array_view<double>(valuesNode);
    for(axom::IndexType i = 0; i < fieldValues.size(); ++i)
    {
      fieldValues[i] = 0.;
    }

    return &fieldNode;
  }
#endif
};

#if defined(AXOM_USE_BUMP)
void printRegisteredFieldNames(const BlueprintState& bpState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling vfSampling,
                               const std::string& initialMessage);

void replaceMaterial(conduit::Node* shapeNode,
                     conduit::Node* materialNode,
                     bool shouldReplace);

void copyShapeIntoMaterial(const conduit::Node* shapeNode,
                           conduit::Node* materialNode,
                           bool reuseExisting = true);

conduit::Node* cloneInOutFunction(const conduit::Node* node);

void generateQuadraturePointMesh(conduit::Node& bpMeshNode,
                                 const std::string& topologyName,
                                 int allocatorID,
                                 axom::ArrayView<int> sampleResolution,
                                 axom::numerics::QuadratureType quadratureType);

void generateSamplingPositions(BlueprintState& bpState,
                               axom::ArrayView<int> sampleResolution,
                               axom::numerics::QuadratureType quadratureType);

void computeVolumeFractionsForMaterial(BlueprintState& bpState, const std::string& matField);

template <int FromDim, int ToDim, typename InsideFunc>
void sampleInOutField(const std::string& shapeName,
                      shaping::BlueprintState& bpState,
                      InsideFunc&& checkInside,
                      PointProjector<FromDim, ToDim> projector = {})
{
  using FromPoint = primal::Point<double, FromDim>;
  using ToPoint = primal::Point<double, ToDim>;
  AXOM_ANNOTATE_SCOPE("sampleInOutField");

  SLIC_ERROR_IF(FromDim != ToDim && !projector,
                "A projector callback function is required when FromDim != ToDim");

  constexpr const char* quadratureCoordsetName = "quadrature_points";
  constexpr const char* quadratureTopologyName = "quadrature_points";
  const std::string inoutName = axom::fmt::format("inout_{}", shapeName);

  conduit::Node& bpMeshNode = bpState.m_internal_node;
  SLIC_ERROR_IF(!bpMeshNode.has_path("coordsets/quadrature_points"),
                "Missing Blueprint quadrature coordset. Generate sampling positions first.");
  SLIC_ERROR_IF(!bpMeshNode.has_path("topologies/quadrature_points"),
                "Missing Blueprint quadrature topology. Generate sampling positions first.");

  conduit::Node& inoutNode = bpMeshNode["fields/" + inoutName];
  inoutNode.reset();
  inoutNode["association"] = "element";
  inoutNode["topology"] = quadratureTopologyName;

  namespace utils = axom::bump::utilities;
  const auto conduitAllocatorId =
    axom::sidre::ConduitMemory::axomAllocIdToConduit(bpState.m_allocator_id);
  conduit::Node& valuesNode = inoutNode["values"];
  valuesNode.set_allocator(conduitAllocatorId);

  axom::utilities::Timer timer(true);
  axom::IndexType numQueryPoints = 0;
  axom::bump::views::dispatch_explicit_coordset(
    bpMeshNode["coordsets/" + std::string(quadratureCoordsetName)], [&](auto coordsetView) {
      using CoordsetView = typename std::decay<decltype(coordsetView)>::type;

      SLIC_ERROR_IF(CoordsetView::dimension() != FromDim,
                    axom::fmt::format("Expected {}D quadrature point coordset, got {}D.",
                                      FromDim,
                                      CoordsetView::dimension()));

      numQueryPoints = coordsetView.size();
      valuesNode.set(conduit::DataType::float64(numQueryPoints));
      auto inoutValues = utils::make_array_view<double>(valuesNode);

      for(axom::IndexType i = 0; i < numQueryPoints; ++i)
      {
        FromPoint fromPt;
        const auto coordsetPoint = coordsetView[i];
        for(int d = 0; d < FromDim; ++d)
        {
          fromPt[d] = coordsetPoint[d];
        }

        const ToPoint queryPt = projector ? projector(fromPt) : ToPoint(fromPt.data());
        inoutValues[i] = checkInside(queryPt) ? 1. : 0.;
      }
    });
  timer.stop();

  SLIC_INFO_ROOT(axom::fmt::format(
    axom::utilities::locale(),
    "\t Sampling inout field '{}' took {:.3Lf} seconds (@ {:L} queries per second)",
    inoutName,
    timer.elapsed(),
    static_cast<int>(numQueryPoints / timer.elapsed())));
}

#endif  // defined(AXOM_USE_BUMP)

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // defined(AXOM_USE_CONDUIT)

#endif  // AXOM_QUEST_SHAPING_HELPERS_BLUEPRINT__HPP_
