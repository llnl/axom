// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file HandleMesh.cpp
 *
 * \brief Demonstrates SLAM sets, relations, and maps on a mesh with strongly typed handles.
 *
 * This example separates the positions used to address SLAM containers
 * from the application-level elements stored in them. 
 * ZoneSet and NodeSet use different PositionTypes 
 * and return distinct ZoneHandle and NodeHandle ElementTypes. 
 * A zone-to-node relation stores positions in the NodeSet, 
 * while maps attach temperature data to entities
 * and interpolation weights to individual zone/node incidences.
 */

#include "axom/core/Array.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <cmath>
#include <cstdint>
#include <ostream>
#include <type_traits>
#include <utility>

namespace slam = axom::slam;

namespace
{
/**
 * \class Handle
 * \brief A strongly typed identifier for a mesh entity.
 *
 * EntityTag makes handles for different entity kinds incompatible even when
 * their underlying integer representations are the same.
 */
template <typename EntityTag, typename IndexType>
struct Handle
{
  using EntityType = EntityTag;
  using IdType = IndexType;

  Handle() = default;
  explicit Handle(IndexType id) : mID(id) { }
  Handle(const Handle&) = default;
  Handle& operator=(const Handle& h) = default;

  bool operator==(const Handle& h) const { return mID == h.mID; }

  static Handle make_handle(IndexType id) { return Handle(id); }

  void print(std::ostream& os) const { os << "Handle(" << mID << ")"; }

  IndexType id() const { return mID; }

  IndexType mID {};
};

template <typename EntityTag, typename IndexType>
std::ostream& operator<<(std::ostream& os, const Handle<EntityTag, IndexType>& h)
{
  h.print(os);
  return os;
}

/**
 * \class HandleMesh
 * \brief Owns typed zone/node handles and their zone-to-node connectivity.
 *
 * The arrays must be declared before the sets and relation because those SLAM
 * objects retain views or pointers into the arrays rather than owning them.
 */
class HandleMesh
{
public:
  // Positions address entries in the sets; handles are the elements returned
  // when those positions are dereferenced.
  struct ZoneTag;
  using ZonePosition = std::int32_t;
  using ZoneHandle = Handle<ZoneTag, ZonePosition>;
  using ZoneSet = slam::ArrayViewIndirectionSet<ZonePosition, ZoneHandle>;

  struct NodeTag;
  using NodePosition = std::int64_t;
  using NodeHandle = Handle<NodeTag, NodePosition>;
  using NodeSet = slam::ArrayViewIndirectionSet<NodePosition, NodeHandle>;

  // A relation row is selected by ZonePosition and contains NodePosition
  // entries. Its flattened CSR storage uses a type wide enough for both.
  using FlatPosition = std::common_type_t<ZonePosition, NodePosition>;
  using ZoneToNodeRelation = slam::VariableRelation<ZoneSet, NodeSet, FlatPosition, NodePosition>;
  using ConnectivitySet = typename slam::RelationSet<ZoneToNodeRelation>::ConcreteSet;

  // Univariate maps attach values to entities; the bivariate map attaches a
  // value to each zone/node incidence in the connectivity relation.
  using NodeTemperatureMap = slam::Map<double, NodeSet>;
  using ZoneTemperatureMap = slam::Map<double, ZoneSet>;
  using ConnectivityWeightMap = slam::BivariateMap<double, ConnectivitySet>;

  // The following asserts validate that our types model the concepts
  static_assert(slam::SetLike<ZoneSet>);
  static_assert(slam::OrderedSetLike<ZoneSet>);
  static_assert(slam::SetLike<NodeSet>);
  static_assert(slam::OrderedSetLike<NodeSet>);
  static_assert(std::is_trivially_copyable_v<ZoneHandle>);
  static_assert(std::is_trivially_copyable_v<NodeHandle>);
  static_assert(!std::is_same_v<ZoneHandle, NodeHandle>);

  static_assert(slam::RelationLike<ZoneToNodeRelation>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::FromPositionType, ZonePosition>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::ToPositionType, NodePosition>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::FlatPositionType, FlatPosition>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::FromSetType::ElementType, ZoneHandle>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::ToSetType::ElementType, NodeHandle>);

  static_assert(slam::BivariateSetLike<ConnectivitySet>);
  static_assert(
    std::is_same_v<typename ConnectivitySet::ElementType, std::pair<ZonePosition, NodePosition>>);
  static_assert(std::is_same_v<typename ConnectivitySet::FirstSetType::ElementType, ZoneHandle>);
  static_assert(std::is_same_v<typename ConnectivitySet::SecondSetType::ElementType, NodeHandle>);

  static_assert(slam::MapOver<NodeTemperatureMap, NodeSet>);
  static_assert(slam::MapOver<ZoneTemperatureMap, ZoneSet>);
  static_assert(slam::BivariateMapLike<ConnectivityWeightMap>);
  static_assert(slam::MapOver<ConnectivityWeightMap, ConnectivitySet>);
  static_assert(std::is_same_v<typename ConnectivityWeightMap::BivariateSetType::FirstSetType::ElementType,
                               ZoneHandle>);
  static_assert(std::is_same_v<typename ConnectivityWeightMap::BivariateSetType::SecondSetType::ElementType,
                               NodeHandle>);

private:
  void invalidateConnectivity()
  {
    m_zoneToNode = ZoneToNodeRelation {};
    m_relationBegins.clear();
    m_nodePositions.clear();
  }

public:
  HandleMesh() = default;

  HandleMesh(const HandleMesh&) = delete;
  HandleMesh& operator=(const HandleMesh&) = delete;
  HandleMesh(HandleMesh&&) = delete;
  HandleMesh& operator=(HandleMesh&&) = delete;

  /// Replace the zone handles and rebuild the view-based set.
  /// Existing connectivity is invalidated because its from-set has changed.
  void setZones(axom::Array<ZoneHandle> zoneHandles)
  {
    invalidateConnectivity();
    m_zoneHandles = std::move(zoneHandles);
    m_zones = slam::make_indirection_set<ZonePosition>(m_zoneHandles);
  }

  /// Replace the node handles and rebuild the view-based set.
  /// Existing connectivity is invalidated because its to-set has changed.
  void setNodes(axom::Array<NodeHandle> nodeHandles)
  {
    invalidateConnectivity();
    m_nodeHandles = std::move(nodeHandles);
    m_nodes = slam::make_indirection_set<NodePosition>(m_nodeHandles);
  }

  /**
   * \brief Set the zone-to-node connectivity from CSR begins and node positions.
   *
   * The relation stores NodePosition values rather than NodeHandles; projecting
   * a relation entry through nodes() produces the corresponding NodeHandle.
   *
   * \pre The zone and node sets have already been populated.
   * \pre relationBegins.size() == zones.size() + 1.
   */
  void setConnectivity(axom::Array<FlatPosition> relationBegins,
                       axom::Array<NodePosition> nodePositions)
  {
    m_relationBegins = std::move(relationBegins);
    m_nodePositions = std::move(nodePositions);
    m_zoneToNode = slam::make_variable_relation(m_zones, m_nodes, m_relationBegins, m_nodePositions);
  }

  const ZoneSet& zones() const { return m_zones; }
  const NodeSet& nodes() const { return m_nodes; }
  ZoneToNodeRelation& zoneToNode() { return m_zoneToNode; }
  const ZoneToNodeRelation& zoneToNode() const { return m_zoneToNode; }

private:
  // These arrays own the storage referenced by the sets and relation below.
  // Their declaration order ensures that they outlive those non-owning views.
  axom::Array<ZoneHandle> m_zoneHandles;
  axom::Array<NodeHandle> m_nodeHandles;
  axom::Array<FlatPosition> m_relationBegins;
  axom::Array<NodePosition> m_nodePositions;

  ZoneSet m_zones;
  NodeSet m_nodes;
  ZoneToNodeRelation m_zoneToNode;
};

}  // end anonymous namespace

int main(int, char**)
{
  axom::slic::SimpleLogger logger;

  using ZonePosition = HandleMesh::ZonePosition;
  using ZoneHandle = HandleMesh::ZoneHandle;
  using NodePosition = HandleMesh::NodePosition;
  using NodeHandle = HandleMesh::NodeHandle;

  // Build a small mixed-element mesh: a triangle and a quadrilateral share
  // nodes. Handle IDs deliberately differ from their zero-based set positions.
  HandleMesh mesh;
  mesh.setZones({ZoneHandle::make_handle(100), ZoneHandle::make_handle(200)});
  mesh.setNodes({NodeHandle::make_handle(10),
                 NodeHandle::make_handle(20),
                 NodeHandle::make_handle(30),
                 NodeHandle::make_handle(40),
                 NodeHandle::make_handle(50)});
  // Zone position 0 uses node positions {0,1,2}; zone position 1 uses
  // {1,3,4,2}. The CSR begins array therefore contains {0,3,7}.
  mesh.setConnectivity({0, 3, 7}, {0, 1, 2, 1, 3, 4, 2});

  // RelationSet presents the relation as a flat bivariate domain, allowing a
  // BivariateMap to store one value for every connected zone/node pair.
  HandleMesh::ConnectivitySet connectivity(&mesh.zoneToNode());

  const auto& zones = mesh.zones();
  const auto& nodes = mesh.nodes();
  const auto& zoneToNode = mesh.zoneToNode();

  // A map over NodeSet owns one temperature value per node position. Range-for
  // is sufficient here because only the values, not their positions, are needed.
  HandleMesh::NodeTemperatureMap nodeTemperature(&nodes);
  double temperature = 10.0;
  for(double& value : nodeTemperature)
  {
    value = temperature;
    temperature += 10.0;
  }

  HandleMesh::ZoneTemperatureMap zoneTemperature(&zones, 0.0);
  HandleMesh::ConnectivityWeightMap interpolationWeight(connectivity, 0.0);

  // Gather nodal temperatures to zone-centered values. The explicit outer
  // iterator provides both ZonePosition (for indexing) and ZoneHandle (for
  // application-facing output). A relation row can use range-for because its
  // elements are already NodePosition values.
  for(auto zoneIt = zones.begin(); zoneIt != zones.end(); ++zoneIt)
  {
    const ZonePosition zonePos = zoneIt.index();
    const ZoneHandle zone = *zoneIt;

    const auto connectedNodes = zoneToNode[zonePos];
    const double weight = 1.0 / static_cast<double>(connectedNodes.size());

    SLIC_INFO("Gathering temperature for zone " << zone << ":");
    for(NodePosition nodePos : connectedNodes)
    {
      // Project the relation's position through NodeSet to obtain its element.
      const NodeHandle node = nodes[nodePos];
      interpolationWeight(zonePos, nodePos) = weight;
      zoneTemperature[zonePos] += weight * nodeTemperature[nodePos];

      SLIC_INFO("  " << node << " has temperature " << nodeTemperature[nodePos]
                     << " and interpolation weight " << weight);
    }

    SLIC_INFO("  resulting zone temperature: " << zoneTemperature[zonePos]);
  }

  // Bivariate-set iteration yields endpoint positions. Projecting each endpoint
  // through its set recovers the strongly typed handles for presentation.
  for(const auto coordinate : connectivity)
  {
    const ZoneHandle zone = zones[coordinate.first];
    const NodeHandle node = nodes[coordinate.second];
    SLIC_INFO("Connectivity " << zone << " -> " << node << " has interpolation weight "
                              << interpolationWeight(coordinate.first, coordinate.second));
  }

  // Range-for over a bivariate map visits its per-incidence values directly.
  // Each zone's interpolation weights sum to one, so the global sum is two.
  double totalInterpolationWeight = 0.0;
  for(const double weight : interpolationWeight)
  {
    totalInterpolationWeight += weight;
  }

  SLIC_ASSERT(zoneToNode.isValid(true));
  SLIC_ASSERT(connectivity.isValid(true));
  SLIC_ASSERT(std::abs(zoneTemperature[0] - 20.0) < 1.e-12);
  SLIC_ASSERT(std::abs(zoneTemperature[1] - 35.0) < 1.e-12);
  SLIC_ASSERT(std::abs(totalInterpolationWeight - 2.0) < 1.e-12);

  return 0;
}
