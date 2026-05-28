// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_GENERATE_QUADRATURE_MESH_HPP_
#define AXOM_QUEST_GENERATE_QUADRATURE_MESH_HPP_

#include "axom/config.hpp"

#if defined(AXOM_USE_CONDUIT)

  #include "MappedZoneUtilities.hpp"
  #include "axom/bump/utilities/blueprint_utilities.hpp"
  #include "axom/bump/utilities/conduit_memory.hpp"
  #include "axom/core.hpp"
  #include "axom/core/numerics/quadrature.hpp"
  #include "axom/primal.hpp"
  #include "axom/sidre/core/ConduitMemory.hpp"
  #include "axom/slic.hpp"

  #include <conduit/conduit.hpp>
  #include <string>
  #include <vector>

/*!
 * \file GenerateQuadratureMesh.hpp
 *
 * \brief Builds a derived Blueprint point mesh from tensor-product quadrature
 *        samples over low-order quad/hex zones.
 */

namespace axom
{
namespace quest
{
namespace shaping
{

/*!
 * \brief Generates a Blueprint point mesh of quadrature samples over an input
 *        topology view.
 *
 * The generated mesh stores one point element per sampled quadrature point and
 * publishes fields that map those points back to their source zones.
 *
 * \tparam ExecSpace The execution space used to populate the generated data.
 * \tparam TopologyView The bump topology view type.
 * \tparam CoordsetView The bump coordset view type.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class GenerateQuadratureMesh
{
public:
  using CoordsetType = typename CoordsetView::value_type;
  using PointType = primal::Point<CoordsetType, CoordsetView::dimension()>;

  /// Struct for capturing views.
  struct ViewPackage
  {
    TopologyView topologyView;
    CoordsetView coordsetView;
  };

  /*!
   * \brief Constructs the generator from a topology and coordset view.
   *
   * \param [in] topologyView The source topology view.
   * \param [in] coordsetView The source coordset view.
   */
  GenerateQuadratureMesh(const TopologyView& topologyView, const CoordsetView& coordsetView)
    : m_topologyView(topologyView)
    , m_coordsetView(coordsetView)
    , m_allocator_id(axom::execution_space<ExecSpace>::allocatorID())
  { }

  /*!
   * \brief Sets the allocator used for generated Conduit-backed storage.
   *
   * \param [in] allocator_id The allocator to use for generated arrays.
   */
  void setAllocatorID(int allocator_id)
  {
    SLIC_ERROR_IF(!axom::isValidAllocatorID(allocator_id), "Invalid allocator id.");
    SLIC_ERROR_IF(!axom::execution_space<ExecSpace>::usesAllocId(allocator_id),
                  "Allocator id is not compatible with execution space.");
    m_allocator_id = allocator_id;
  }

  int getAllocatorID() const { return m_allocator_id; }

  /*!
   * \brief Executes the quadrature-point mesh generation.
   *
   * \param [in] n_topology The source topology node.
   * \param [in] n_coordset The source coordset node.
   * \param [in] outputTopologyName The generated point-topology name.
   * \param [in] outputCoordsetName The generated coordset name.
   * \param [in] originalElementsFieldName The generated provenance field name.
   * \param [in] quadratureWeightsFieldName The generated reference-weight field
   *            name.
   * \param [in] quadraturePhysicalWeightsFieldName The generated physical-weight
   *            field name.
   * \param [in] ruleX The quadrature rule in the first logical direction.
   * \param [in] ruleY The quadrature rule in the second logical direction.
   * \param [in] ruleZ The quadrature rule in the third logical direction.
   * \param [in,out] n_output The Blueprint mesh tree to augment.
   */
  void execute(const conduit::Node& AXOM_UNUSED_PARAM(n_topology),
               const conduit::Node& n_coordset,
               const std::string& outputTopologyName,
               const std::string& outputCoordsetName,
               const std::string& originalElementsFieldName,
               const std::string& quadratureWeightsFieldName,
               const std::string& quadraturePhysicalWeightsFieldName,
               const numerics::QuadratureRule& ruleX,
               const numerics::QuadratureRule& ruleY,
               const numerics::QuadratureRule& ruleZ,
               conduit::Node& n_output) const
  {
    namespace utils = axom::bump::utilities;

    const auto numZones = m_topologyView.numberOfZones();
    const int dim = CoordsetView::dimension();
    const int npts = detail::quadraturePointCount(ruleX, ruleY, ruleZ, dim);
    const IndexType numPoints = numZones * static_cast<IndexType>(npts);
    const auto conduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(getAllocatorID());

    const std::vector<std::string> axes = utils::coordsetAxes(n_coordset);
    SLIC_ASSERT(static_cast<int>(axes.size()) == dim);

    conduit::Node& n_outputCoordset = n_output["coordsets/" + outputCoordsetName];
    n_outputCoordset.reset();
    n_outputCoordset["type"] = "explicit";

    // Store the sampled coordinates as plain explicit coordset components.
    axom::StackArray<axom::ArrayView<CoordsetType>, CoordsetView::dimension()> coordViews;
    for(int d = 0; d < dim; ++d)
    {
      conduit::Node& comp = n_outputCoordset["values/" + axes[d]];
      comp.set_allocator(conduitAllocatorId);
      comp.set(conduit::DataType(utils::cpp2conduit<CoordsetType>::id, numPoints));
      coordViews[d] = utils::make_array_view<CoordsetType>(comp);
    }

    conduit::Node& n_outputTopo = n_output["topologies/" + outputTopologyName];
    n_outputTopo.reset();
    n_outputTopo["type"] = "unstructured";
    n_outputTopo["coordset"] = outputCoordsetName;
    n_outputTopo["elements/shape"] = "point";

    // The derived topology is a point mesh, so connectivity is the identity.
    conduit::Node& n_connectivity = n_outputTopo["elements/connectivity"];
    n_connectivity.set_allocator(conduitAllocatorId);
    n_connectivity.set(conduit::DataType::index_t(numPoints));
    auto connectivity = utils::make_array_view<conduit::index_t>(n_connectivity);

    conduit::Node& n_sizes = n_outputTopo["elements/sizes"];
    n_sizes.set_allocator(conduitAllocatorId);
    n_sizes.set(conduit::DataType::index_t(numPoints));
    auto sizes = utils::make_array_view<conduit::index_t>(n_sizes);

    conduit::Node& n_offsets = n_outputTopo["elements/offsets"];
    n_offsets.set_allocator(conduitAllocatorId);
    n_offsets.set(conduit::DataType::index_t(numPoints));
    auto offsets = utils::make_array_view<conduit::index_t>(n_offsets);

    conduit::Node& n_originalElements = n_output["fields/" + originalElementsFieldName];
    n_originalElements.reset();
    n_originalElements["association"] = "element";
    n_originalElements["topology"] = outputTopologyName;
    conduit::Node& n_originalValues = n_originalElements["values"];
    n_originalValues.set_allocator(conduitAllocatorId);
    n_originalValues.set(conduit::DataType::index_t(numPoints));
    auto originalElements = utils::make_array_view<conduit::index_t>(n_originalValues);

    conduit::Node& n_quadratureWeights = n_output["fields/" + quadratureWeightsFieldName];
    n_quadratureWeights.reset();
    n_quadratureWeights["association"] = "element";
    n_quadratureWeights["topology"] = outputTopologyName;
    conduit::Node& n_weightValues = n_quadratureWeights["values"];
    n_weightValues.set_allocator(conduitAllocatorId);
    n_weightValues.set(conduit::DataType::float64(numPoints));
    auto quadratureWeights = utils::make_array_view<double>(n_weightValues);

    conduit::Node& n_physicalQuadratureWeights =
      n_output["fields/" + quadraturePhysicalWeightsFieldName];
    n_physicalQuadratureWeights.reset();
    n_physicalQuadratureWeights["association"] = "element";
    n_physicalQuadratureWeights["topology"] = outputTopologyName;
    conduit::Node& n_physicalWeightValues = n_physicalQuadratureWeights["values"];
    n_physicalWeightValues.set_allocator(conduitAllocatorId);
    n_physicalWeightValues.set(conduit::DataType::float64(numPoints));
    auto physicalQuadratureWeights = utils::make_array_view<double>(n_physicalWeightValues);

    // Package these views into a struct to help with device access.
    const ViewPackage deviceViews {m_topologyView, m_coordsetView};

    axom::for_all<ExecSpace>(
      numZones,
      AXOM_LAMBDA(IndexType zoneIndex) {
        const auto zone = deviceViews.topologyView.zone(zoneIndex);
        IndexType pointIndex = zoneIndex * static_cast<IndexType>(npts);

        for(int kz = 0; kz < (dim == 3 ? ruleZ.getNumPoints() : 1); ++kz)
        {
          const double zeta = dim == 3 ? ruleZ.node(kz) : 0.0;
          const double wz = dim == 3 ? ruleZ.weight(kz) : 1.0;
          for(int jy = 0; jy < ruleY.getNumPoints(); ++jy)
          {
            const double eta = ruleY.node(jy);
            const double wy = ruleY.weight(jy);
            for(int ix = 0; ix < ruleX.getNumPoints(); ++ix)
            {
              const double xi = ruleX.node(ix);
              const double wx = ruleX.weight(ix);

              PointType pt;
              double physicalMeasure = 0.;
              if constexpr(CoordsetView::dimension() == 2)
              {
                pt = detail::mapToPhysicalPoint(zone, deviceViews.coordsetView, xi, eta);
                physicalMeasure =
                  detail::computePhysicalMeasureFactor(zone, deviceViews.coordsetView, xi, eta);
              }
              else
              {
                pt = detail::mapToPhysicalPoint(zone, deviceViews.coordsetView, xi, eta, zeta);
                physicalMeasure =
                  detail::computePhysicalMeasureFactor(zone, deviceViews.coordsetView, xi, eta, zeta);
              }

              // Retain both the reference-space tensor-product weights and the
              // Jacobian-weighted physical weights for downstream consumers.
              const double referenceWeight = wx * wy * wz;
              for(int d = 0; d < dim; ++d)
              {
                coordViews[d][pointIndex] = pt[d];
              }
              connectivity[pointIndex] = pointIndex;
              sizes[pointIndex] = 1;
              offsets[pointIndex] = pointIndex;
              originalElements[pointIndex] = zoneIndex;
              quadratureWeights[pointIndex] = referenceWeight;
              physicalQuadratureWeights[pointIndex] = referenceWeight * physicalMeasure;
              ++pointIndex;
            }
          }
        }
      });
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  int m_allocator_id;
};

}  // namespace shaping
}  // namespace quest
}  // namespace axom

#endif

#endif
