// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_GENERATE_QUADRATURE_MESH_HPP_
#define AXOM_QUEST_GENERATE_QUADRATURE_MESH_HPP_

#include "axom/config.hpp"

#if defined(AXOM_USE_CONDUIT)

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

namespace axom
{
namespace quest
{
namespace shaping
{
namespace detail
{

template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE primal::Point<typename CoordsetView::value_type, 2> mapToPhysicalPoint(
  const ShapeType& zone,
  const CoordsetView& coordsetView,
  double u,
  double v)
{
  using PointType = primal::Point<typename CoordsetView::value_type, 2>;
  const auto p0 = coordsetView[zone.getId(0)];
  const auto p1 = coordsetView[zone.getId(1)];
  const auto p2 = coordsetView[zone.getId(2)];
  const auto p3 = coordsetView[zone.getId(3)];

  const double n0 = (1.0 - u) * (1.0 - v);
  const double n1 = u * (1.0 - v);
  const double n2 = u * v;
  const double n3 = (1.0 - u) * v;

  PointType pt;
  for(int d = 0; d < 2; ++d)
  {
    pt[d] = n0 * p0[d] + n1 * p1[d] + n2 * p2[d] + n3 * p3[d];
  }
  return pt;
}

template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE primal::Point<typename CoordsetView::value_type, 3> mapToPhysicalPoint(
  const ShapeType& zone,
  const CoordsetView& coordsetView,
  double u,
  double v,
  double w)
{
  using PointType = primal::Point<typename CoordsetView::value_type, 3>;
  const auto p0 = coordsetView[zone.getId(0)];
  const auto p1 = coordsetView[zone.getId(1)];
  const auto p2 = coordsetView[zone.getId(2)];
  const auto p3 = coordsetView[zone.getId(3)];
  const auto p4 = coordsetView[zone.getId(4)];
  const auto p5 = coordsetView[zone.getId(5)];
  const auto p6 = coordsetView[zone.getId(6)];
  const auto p7 = coordsetView[zone.getId(7)];

  const double a = 1.0 - u;
  const double b = 1.0 - v;
  const double c = 1.0 - w;

  const double n0 = a * b * c;
  const double n1 = u * b * c;
  const double n2 = u * v * c;
  const double n3 = a * v * c;
  const double n4 = a * b * w;
  const double n5 = u * b * w;
  const double n6 = u * v * w;
  const double n7 = a * v * w;

  PointType pt;
  for(int d = 0; d < 3; ++d)
  {
    pt[d] = n0 * p0[d] + n1 * p1[d] + n2 * p2[d] + n3 * p3[d] + n4 * p4[d] + n5 * p5[d] +
      n6 * p6[d] + n7 * p7[d];
  }
  return pt;
}

inline int quadraturePointCount(const numerics::QuadratureRule& ruleX,
                                const numerics::QuadratureRule& ruleY,
                                const numerics::QuadratureRule& ruleZ,
                                int dim)
{
  return dim == 2 ? ruleX.getNumPoints() * ruleY.getNumPoints()
                  : ruleX.getNumPoints() * ruleY.getNumPoints() * ruleZ.getNumPoints();
}

}  // namespace detail

template <typename ExecSpace, typename TopologyView, typename CoordsetView>
class GenerateQuadratureMesh
{
public:
  using CoordsetType = typename CoordsetView::value_type;
  using PointType = primal::Point<CoordsetType, CoordsetView::dimension()>;

  GenerateQuadratureMesh(const TopologyView& topologyView, const CoordsetView& coordsetView)
    : m_topologyView(topologyView)
    , m_coordsetView(coordsetView)
    , m_allocator_id(axom::execution_space<ExecSpace>::allocatorID())
  { }

  void setAllocatorID(int allocator_id)
  {
    SLIC_ERROR_IF(!axom::isValidAllocatorID(allocator_id), "Invalid allocator id.");
    SLIC_ERROR_IF(!axom::execution_space<ExecSpace>::usesAllocId(allocator_id),
                  "Allocator id is not compatible with execution space.");
    m_allocator_id = allocator_id;
  }

  int getAllocatorID() const { return m_allocator_id; }

  void execute(const conduit::Node& n_topology,
               const conduit::Node& n_coordset,
               const std::string& outputTopologyName,
               const std::string& outputCoordsetName,
               const std::string& originalElementsFieldName,
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

    const TopologyView deviceTopoView(m_topologyView);
    const CoordsetView deviceCoordsetView(m_coordsetView);

    axom::for_all<ExecSpace>(
      numZones,
      AXOM_LAMBDA(IndexType zoneIndex) {
        const auto zone = deviceTopoView.zone(zoneIndex);
        IndexType pointIndex = zoneIndex * static_cast<IndexType>(npts);

        for(int kz = 0; kz < (dim == 3 ? ruleZ.getNumPoints() : 1); ++kz)
        {
          const double zeta = dim == 3 ? ruleZ.node(kz) : 0.0;
          for(int jy = 0; jy < ruleY.getNumPoints(); ++jy)
          {
            const double eta = ruleY.node(jy);
            for(int ix = 0; ix < ruleX.getNumPoints(); ++ix)
            {
              const double xi = ruleX.node(ix);

              PointType pt;
              if constexpr(CoordsetView::dimension() == 2)
              {
                pt = detail::mapToPhysicalPoint(zone, deviceCoordsetView, xi, eta);
              }
              else
              {
                pt = detail::mapToPhysicalPoint(zone, deviceCoordsetView, xi, eta, zeta);
              }

              for(int d = 0; d < dim; ++d)
              {
                coordViews[d][pointIndex] = pt[d];
              }
              connectivity[pointIndex] = pointIndex;
              sizes[pointIndex] = 1;
              offsets[pointIndex] = pointIndex;
              originalElements[pointIndex] = zoneIndex;
              ++pointIndex;
            }
          }
        }
      });

    AXOM_UNUSED_VAR(n_topology);
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  int m_allocator_id;
};

}  // namespace shaping
}  // namespace quest
}  // namespace axom

#endif

#endif
