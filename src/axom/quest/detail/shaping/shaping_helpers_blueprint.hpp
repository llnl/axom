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

#if defined(AXOM_USE_BUMP)
constexpr const char* QUADRATURE_COORDSET_NAME = "quadrature_points";
constexpr const char* QUADRATURE_TOPOLOGY_NAME = "quadrature_points";
constexpr const char* ORIGINAL_ELEMENTS_FIELD_NAME = "originalElements";
constexpr const char* QUADRATURE_WEIGHTS_FIELD_NAME = "quadratureWeights";
constexpr const char* QUADRATURE_PHYSICAL_WEIGHTS_FIELD_NAME = "quadraturePhysicalWeights";
#endif

/*!
 * \brief Stores Blueprint mesh state and backend-specific access for Quest shapers.
 *
 * A `BlueprintState` may be backed either by a caller-owned `sidre::Group`
 * or by a caller-owned `conduit::Node`. The helper methods on this struct
 * provide a single access layer for reading and updating Blueprint mesh data
 * without forcing the shaper code to know which storage backend is active.
 */
struct BlueprintState
{
  /// Destructor.
  virtual ~BlueprintState() = default;

  /// The caller-owned Sidre group backing the mesh, when Sidre-backed.
  axom::sidre::Group* m_group_ptr {nullptr};
  /// Allocator id used for newly-created Blueprint array data.
  int m_allocator_id {axom::getDefaultAllocatorID()};
  /// Name of the active Blueprint topology.
  std::string m_topology_name;
  /// The caller-owned Blueprint mesh node, when Conduit-backed.
  conduit::Node* m_external_node_ptr {nullptr};
  /// Cached native Conduit layout for Sidre-backed meshes.
  conduit::Node m_internal_node;

  /// Return whether the active Blueprint mesh is backed by Sidre storage.
  bool isSidreBacked() const { return m_group_ptr != nullptr; }
  /// Return whether the active Blueprint mesh is backed by a Conduit node.
  bool isConduitBacked() const { return m_external_node_ptr != nullptr; }

  /*!
   * \brief Refresh the cached native Conduit layout for Sidre-backed meshes.
   *
   * This is a no-op for Conduit-backed meshes.
   */
  void refreshBlueprintMeshNode();

  /*!
   * \brief Return the dimension implied by the active Blueprint cell shape.
   *
   * \return `2` for quadrilateral meshes or `3` for hexahedral meshes.
   */
  int meshDimension() const;

  /*!
   * \brief Convert a structured Blueprint mesh to unstructured topology in place.
   *
   * \param execPolicy Runtime policy used by the conversion helper.
   */
  void ensureUnstructured(axom::runtime_policy::Policy execPolicy);

  /// Return the active Blueprint mesh node for the current backing store.
  conduit::Node& getBlueprintMeshNode()
  {
    return isConduitBacked() ? *m_external_node_ptr : m_internal_node;
  }

  /// Return the active Blueprint mesh node for the current backing store.
  const conduit::Node& getBlueprintMeshNode() const
  {
    return isConduitBacked() ? *m_external_node_ptr : m_internal_node;
  }

  /// Return the active Blueprint topology node.
  const conduit::Node& getBlueprintTopologyNode() const
  {
    return getBlueprintMeshNode().fetch_existing("topologies").fetch_existing(m_topology_name);
  }

  /// Return a named Blueprint coordset node.
  conduit::Node& getBlueprintCoordsetNode(const std::string& name)
  {
    return getBlueprintMeshNode().fetch_existing("coordsets").fetch_existing(name);
  }

  /// Return a named Blueprint coordset node.
  const conduit::Node& getBlueprintCoordsetNode(const std::string& name) const
  {
    return getBlueprintMeshNode().fetch_existing("coordsets").fetch_existing(name);
  }

  /// Return whether a named Blueprint field is present.
  bool hasField(const std::string& name) const
  {
    return getBlueprintMeshNode().has_path(axom::fmt::format("fields/{}", name));
  }

  /// Return a named Blueprint field node.
  conduit::Node& getField(const std::string& name)
  {
    return getBlueprintMeshNode().fetch_existing("fields").fetch_existing(name);
  }

  /// Return a named Blueprint field node.
  const conduit::Node& getField(const std::string& name) const
  {
    return getBlueprintMeshNode().fetch_existing("fields").fetch_existing(name);
  }

  /// Return a shape in/out field, if present.
  conduit::Node* getShapeFunction(const std::string& name)
  {
    return hasField(name) ? &getField(name) : nullptr;
  }

  /// Return a shape in/out field, if present.
  const conduit::Node* getShapeFunction(const std::string& name) const
  {
    return hasField(name) ? &getField(name) : nullptr;
  }

  /*!
   * \brief Remove a shape in/out field from the active Blueprint mesh.
   *
   * \param name The field name to remove.
   */
  void deleteShapeFunction(const std::string& name);

  /// Return a material in/out field, if present.
  conduit::Node* getMaterialFunction(const std::string& name)
  {
    return hasField(name) ? &getField(name) : nullptr;
  }

  /// Return a material in/out field, if present.
  const conduit::Node* getMaterialFunction(const std::string& name) const
  {
    return hasField(name) ? &getField(name) : nullptr;
  }

  /*!
   * \brief Create or replace an element-associated Blueprint field.
   *
   * \param name The field name.
   * \param topologyName The associated topology name.
   * \param size The number of scalar values to allocate.
   * \param addVolumeDependent Whether to add the `volume_dependent` metadata.
   * \param volumeDependent Value to store in `volume_dependent` when requested.
   *
   * \return A writable view over the allocated field values.
   */
  axom::ArrayView<double> createField(const std::string& name,
                                      const std::string& topologyName,
                                      axom::IndexType size,
                                      bool addVolumeDependent = false,
                                      bool volumeDependent = false);

  #if defined(AXOM_USE_BUMP)
  /*!
   * \brief Import generated quadrature Blueprint objects into the active mesh.
   *
   * \param quadratureMesh A mesh node containing the generated quadrature
   *        coordset, topology, and support fields.
   */
  void importQuadraturePointMesh(const conduit::Node& quadratureMesh);

  /*!
   * \brief Create a zero-initialized material in/out field on the quadrature topology.
   *
   * \param name The field name to create.
   *
   * \return The created field node.
   */
  conduit::Node* createMaterialFunction(const std::string& name);
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
void generateQuadraturePointMesh(const conduit::Node& bpMeshNode,
                                 conduit::Node& outputMeshNode,
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

  conduit::Node& bpMeshNode = bpState.getBlueprintMeshNode();
  SLIC_ERROR_IF(!bpMeshNode.has_path(axom::fmt::format("coordsets/{}", QUADRATURE_COORDSET_NAME)),
                "Missing Blueprint quadrature coordset. Generate sampling positions first.");
  SLIC_ERROR_IF(!bpMeshNode.has_path(axom::fmt::format("topologies/{}", QUADRATURE_TOPOLOGY_NAME)),
                "Missing Blueprint quadrature topology. Generate sampling positions first.");

  const std::string inoutName = shaping::shapeInOutFieldName(shapeName);

  axom::utilities::Timer timer(true);
  axom::IndexType numQueryPoints = 0;
  axom::bump::views::dispatch_explicit_coordset(
    bpMeshNode["coordsets/" + std::string(QUADRATURE_COORDSET_NAME)],
    [&](auto coordsetView) {
      using CoordsetView = typename std::decay<decltype(coordsetView)>::type;

      // Limit to handling coordsets whose dimensions match FromDim.
      if constexpr(CoordsetView::dimension() == FromDim)
      {
        numQueryPoints = coordsetView.size();
        auto inoutValues =
          bpState.createField(inoutName, QUADRATURE_TOPOLOGY_NAME, numQueryPoints);

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
