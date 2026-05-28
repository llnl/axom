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

  #include <map>
  #include <set>
  #include <string>
  #include <vector>

namespace axom
{
namespace quest
{
namespace shaping
{
/*!
 * \brief Return the cell shape for a Blueprint topology.
 *
 * \param topoNode The Blueprint topology being queried.
 *
 * \return A string containing the cell shape for the topology.
 */
std::string getBlueprintCellShape(const conduit::Node& topoNode);

/// A class that contains Blueprint mesh and field state for SamplingShaper class.
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

  conduit::Node& getBlueprintMeshNode() { return m_internal_node; }

  const conduit::Node& getBlueprintTopologyNode() const
  {
    return m_internal_node.fetch_existing("topologies").fetch_existing(m_topology_name);
  }

  conduit::Node* getShapeFunction(const std::string& name)
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name] : nullptr;
  }

  const conduit::Node* getShapeFunction(const std::string& name) const
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name] : nullptr;
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
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name] : nullptr;
  }

  const conduit::Node* getMaterialFunction(const std::string& name) const
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name] : nullptr;
  }

#if defined(AXOM_USE_BUMP)
  conduit::Node* createMaterialFunction(const std::string& name)
  {
    constexpr const char* quadratureTopologyName = "quadrature_points";
    SLIC_ERROR_IF(
      !m_internal_node.has_path("coordsets/quadrature_points/values"),
      std::string("Cannot create material function '") + name + "' without quadrature points.");

    conduit::Node& fieldNode = m_internal_node["fields/" + name];
    fieldNode.reset();
    fieldNode["association"] = "element";
    fieldNode["topology"] = quadratureTopologyName;

    const auto conduitAllocatorId = axom::sidre::ConduitMemory::axomAllocIdToConduit(m_allocator_id);
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
/*!
 * \brief Print the registered field names in the \a bpState.
 *
 * \param bpState The Blueprint state.
 * \param knownMaterials A set of known material names.
 * \param vfSampling The type of volume fraction sampling being performed.
 * \param initialMessage A string to prepend to the printed message.
 */
void printRegisteredFieldNames(const BlueprintState& bpState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling vfSampling,
                               const std::string& initialMessage);

/*!
 * Utility function to zero out inout quadrature points for a material replaced by a shape
 *
 * Each location in space can only be covered by one material.
 * When \a shouldReplace is true, we clear all values in \a materialQFunc 
 * that are set in \a shapeQFunc. When it is false, we do the opposite.
 *
 * \param shapeNode The node that contains the shape function.
 * \param materialNode The node that contains the material function.
 * \param shouldReplace Flag for whether the shape replaces the material 
 *   or whether the material remains and we should zero out the shape sample (when false)
 */
void replaceMaterial(conduit::Node* shapeNode, conduit::Node* materialNode, bool shouldReplace);

/*!
 * \brief Utility function to copy inout quadrature point values from \a shapeNode to \a materialNode
 *
 * \param shapeNode The inout samples field for the current shape
 * \param materialNode The inout samples field for the material we're writing into
 * \param reuseExisting When a value is not set in \a shapeNode, should we retain existing values 
 * from \a materialNode or overwrite them based on \a shapeNode. The default is to retain values
 */
void copyShapeIntoMaterial(const conduit::Node* shapeNode,
                           conduit::Node* materialNode,
                           bool reuseExisting = true);

/*!
 * \brief Create a copy of the supplied field.
 *
 * \param node A pointer to the field to clone.
 *
 * \return A pointer to a new copy of the supplied field.
 */
conduit::Node* cloneInOutFunction(const conduit::Node* node);

/*!
 * \brief Generate sampling positions within each zone based on element quadrature, creating a new topology.
 *
 * \param bpMeshNode The node that will contain the new quadrature point mesh topology.
 * \param topologyName The name of the new topology to create.
 * \param allocatorID The allocator Id to use for allocating memory.
 * \param sampleResolution The number of samples in each dimension.
 * \param quadratureType The quadrature type that determines the sample locations.
 */
void generateQuadraturePointMesh(conduit::Node& bpMeshNode,
                                 const std::string& topologyName,
                                 int allocatorID,
                                 axom::ArrayView<int> sampleResolution,
                                 axom::numerics::QuadratureType quadratureType);

/*!
 * \brief Generates sampling positions within each zone based on element quadrature.
 *
 * \param bpState The Blueprint state.
 * \param sampleResolution The number of samples in each dimension.
 * \param quadratureType The quadrature type that determines the sample locations.
 *
 * \note The sample points are stored as a new quadrature_points topology.
 */
void generateSamplingPositions(BlueprintState& bpState,
                               axom::ArrayView<int> sampleResolution,
                               axom::numerics::QuadratureType quadratureType);

/*!
 * \brief Import initial volume fractions from the map into the quadrature
 *        "mat_inout_" fields in \a bpState.
 *
 * \param bpState The Blueprint state.
 * \param initialVolumeFractions A map of initial volume fraction fields used to
 *                               initialize mat_inout fields over the quadrature
 *                               points.
 */
void importInitialVolumeFractions(BlueprintState& bpState,
                                  const std::map<std::string, conduit::Node*>& initialVolumeFractions);

/*!
 * \brief Create volume fractions for a material using the existing material field
 *        (mat_inout_{matField}) to make the new field (vol_fract_{matField}).
 *
 * \param bpState The Blueprint state that contains the mesh and functions.
 * \param matField The name of the material field.
 */
void computeVolumeFractionsForMaterial(BlueprintState& bpState, const std::string& matField);

/*!
  * \brief Samples the inout field over the indexed geometry, possibly using a
  * callback function to project the input points (from the computational mesh)
  * to query points on the spatial index
  *
  * \tparam FromDim The dimension of points from the input mesh
  * \tparam ToDim The dimension of points on the indexed shape
  * \tparam InsideFunc A function that takes a point and returns a bool indicating whether the
  *                    point is inside or outside of relevant shapes.
  *
  * \param [in] shapeName The name of the shape used in making data array names.
  * \param [in] mfemState The data collection containing the mesh, associated query points
  *                       and a collection of quadrature functions for the shape and material
  *                       inout samples.
  * \param [in] checkInside The function that determines whether a point is inside.
  * \param [in] projector A callback function to apply to points from the input mesh
  * before querying them on the spatial index
  *
  * \note A projector callback must be supplied when \a FromDim is not equal
  *       to \a ToDim.
  */
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
    bpMeshNode["coordsets/" + std::string(quadratureCoordsetName)],
    [&](auto coordsetView) {
      using CoordsetView = typename std::decay<decltype(coordsetView)>::type;

      // Limit to handling coordsets whose dimensions match FromDim.
      if constexpr (CoordsetView::dimension() == FromDim)
      {
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
