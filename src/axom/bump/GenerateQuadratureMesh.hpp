// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"

#include "axom/bump/MappedZoneUtilities.hpp"
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
namespace bump
{
namespace detail
{

/*!
 * \brief Returns the number of tensor-product quadrature points per zone.
 *
 * \param ruleX Quadrature rule in the reference x direction.
 * \param ruleY Quadrature rule in the reference y direction.
 * \param ruleZ Quadrature rule in the reference z direction.
 * \param dim Mesh dimension, expected to be 2 or 3.
 *
 * \return The number of quadrature points generated for one zone.
 */
inline int quadraturePointCount(const numerics::QuadratureRule& ruleX,
                                const numerics::QuadratureRule& ruleY,
                                const numerics::QuadratureRule& ruleZ,
                                int dim)
{
  return dim == 2 ? ruleX.getNumPoints() * ruleY.getNumPoints()
                  : ruleX.getNumPoints() * ruleY.getNumPoints() * ruleZ.getNumPoints();
}

}  // namespace detail

/*!
 * \brief Builds a Blueprint point mesh of quadrature samples over a low-order
 *        quad or hex mesh.
 *
 * The generated mesh stores one point element for each tensor-product
 * quadrature point in each zone. It also creates fields that record the
 * originating zone and both reference-space and physical-space quadrature
 * weights for each generated point.
 *
 * \tparam ExecSpace Execution space used to populate the output arrays.
 * \tparam TopologyView View type for the source topology.
 * \tparam CoordsetView View type for the source coordset.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class GenerateQuadratureMesh
{
public:
  using CoordsetType = typename CoordsetView::value_type;
  using PointType = primal::Point<CoordsetType, CoordsetView::dimension()>;

  /*!
   * \brief Bundles the topology and coordset views for device access.
   */
  struct ViewPackage
  {
    TopologyView topologyView;
    CoordsetView coordsetView;
  };

  /*!
   * \brief Constructs a quadrature-mesh generator over the supplied mesh views.
   *
   * \param topologyView View of the source Blueprint topology.
   * \param coordsetView View of the source Blueprint coordset.
   */
  GenerateQuadratureMesh(const TopologyView& topologyView, const CoordsetView& coordsetView)
    : m_topologyView(topologyView)
    , m_coordsetView(coordsetView)
    , m_allocator_id(axom::execution_space<ExecSpace>::allocatorID())
  { }

  /*!
   * \brief Sets the allocator used for output Conduit arrays.
   *
   * The allocator must be valid and accessible from \a ExecSpace.
   *
   * \param allocator_id Axom allocator identifier for output storage.
   */
  void setAllocatorID(int allocator_id)
  {
    SLIC_ERROR_IF(!axom::isValidAllocatorID(allocator_id), "Invalid allocator id.");
    SLIC_ERROR_IF(!axom::execution_space<ExecSpace>::usesAllocId(allocator_id),
                  "Allocator id is not compatible with execution space.");
    m_allocator_id = allocator_id;
  }

  /*!
   * \brief Returns the allocator used for output Conduit arrays.
   *
   * \return The Axom allocator identifier for output storage.
   */
  int getAllocatorID() const { return m_allocator_id; }

  /*!
   * \brief Generates a Blueprint point mesh containing quadrature samples.
   *
   * The output mesh contains explicit point coordinates, a point-element
   * topology, the source zone index for each point, and both reference and
   * physical quadrature weights.
   *
   * \param n_topology Source topology node. Currently unused because the
   *        topology is accessed through \a m_topologyView.
   * \param n_coordset Source coordset node, used to preserve axis naming.
   * \param outputTopologyName Name of the generated point topology.
   * \param outputCoordsetName Name of the generated explicit coordset.
   * \param originalElementsFieldName Name of the field storing source zone ids.
   * \param quadratureWeightsFieldName Name of the field storing reference
   *        quadrature weights.
   * \param quadraturePhysicalWeightsFieldName Name of the field storing
   *        physical quadrature weights.
   * \param ruleX Quadrature rule in the reference x direction.
   * \param ruleY Quadrature rule in the reference y direction.
   * \param ruleZ Quadrature rule in the reference z direction.
   * \param n_output Output Blueprint node populated with the generated mesh.
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

#if !defined(__CUDACC__)
private:
#endif
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  int m_allocator_id;
};

}  // namespace bump
}  // namespace axom
