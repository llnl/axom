// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file HandleMesh.cpp
 *
 * \brief Example where the ElementType of a set is a handle that wraps a mesh index.
 */

#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <cstdint>
#include <ostream>
#include <type_traits>
#include <vector>

namespace slam = axom::slam;

namespace
{
/**
 * \class Handle
 * \brief A handle is a wrapper around a variable.
 *
 * This class is a stub for some mesh data structures that use typed index handles
 * rather than (untyped) indices to refer to mesh elements.
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

}  // end anonymous namespace

int main(int, char**)
{
  axom::slic::SimpleLogger logger;

  struct ZoneTag;
  using ZonePosition = std::int32_t;
  using ZoneHandle = Handle<ZoneTag, ZonePosition>;
  using ZoneSet = slam::VectorIndirectionSet<ZonePosition, ZoneHandle>;
  static_assert(slam::SetLike<ZoneSet>);
  static_assert(std::is_trivially_copyable_v<ZoneHandle>);
  static_assert(slam::OrderedSetLike<ZoneSet>);
  std::vector<ZoneHandle> zoneHandles {ZoneHandle::make_handle(100), ZoneHandle::make_handle(200)};
  ZoneSet zones = slam::make_indirection_set<ZonePosition>(zoneHandles);

  struct NodeTag;
  using NodePosition = std::int64_t;
  using NodeHandle = Handle<NodeTag, NodePosition>;
  using NodeSet = slam::VectorIndirectionSet<NodePosition, NodeHandle>;
  static_assert(slam::SetLike<NodeSet>);
  static_assert(std::is_trivially_copyable_v<NodeHandle>);
  static_assert(slam::OrderedSetLike<NodeSet>);

  std::vector<NodeHandle> nodeHandles {NodeHandle::make_handle(10),
                                       NodeHandle::make_handle(20),
                                       NodeHandle::make_handle(30),
                                       NodeHandle::make_handle(40)};
  NodeSet nodes = slam::make_indirection_set<NodePosition>(nodeHandles);

  static_assert(!std::is_same_v<ZoneHandle, NodeHandle>);

  // Set up connectivity. 
  // Begins and flat positions use the common position type.
  // Relation entries retain NodeSet::PositionType.
  using FlatPosition = std::common_type_t<ZonePosition, NodePosition>;
  std::vector<FlatPosition> begins {0, 3, 6};
  std::vector<NodePosition> nodePositions {0, 1, 2, 1, 2, 3};

  auto zoneToNode = slam::make_variable_relation(zones, nodes, begins, nodePositions);
  using ZoneToNodeRelation = decltype(zoneToNode);

  static_assert(slam::RelationLike<ZoneToNodeRelation>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::FromPositionType, ZonePosition>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::ToPositionType, NodePosition>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::FlatPositionType, FlatPosition>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::FromSetType::ElementType, ZoneHandle>);
  static_assert(std::is_same_v<typename ZoneToNodeRelation::ToSetType::ElementType, NodeHandle>);

  using ConnectivitySet = typename slam::RelationSet<ZoneToNodeRelation>::ConcreteSet;
  ConnectivitySet connectivity(&zoneToNode);

  static_assert(slam::BivariateSetLike<ConnectivitySet>);
  static_assert(
    std::is_same_v<typename ConnectivitySet::ElementType, std::pair<ZonePosition, NodePosition>>);
  static_assert(std::is_same_v<typename ConnectivitySet::FirstSetType::ElementType, ZoneHandle>);
  static_assert(std::is_same_v<typename ConnectivitySet::SecondSetType::ElementType, NodeHandle>);

  using ConnectivityMap = slam::BivariateMap<double, ConnectivitySet>;
  static_assert(slam::BivariateMapLike<ConnectivityMap>);
  static_assert(slam::MapOver<ConnectivityMap, ConnectivitySet>);
  static_assert(
    std::is_same_v<typename ConnectivityMap::BivariateSetType::FirstSetType::ElementType, ZoneHandle>);
  static_assert(
    std::is_same_v<typename ConnectivityMap::BivariateSetType::SecondSetType::ElementType, NodeHandle>);

  ConnectivityMap weights(connectivity, 0.0);

  // Process the mesh by traversing raw connectivity positions and projecting
  // them through the endpoint sets to obtain distinct, strongly typed handles.
  NodePosition nodeIdSum = 0;
  for(ZonePosition zonePos = 0; zonePos < zones.size(); ++zonePos)
  {
    const ZoneHandle zone = zones[zonePos];
    const auto relatedNodes = zoneToNode[zonePos];

    SLIC_INFO("Zone " << zone << " contains:");
    for(FlatPosition local = 0; local < relatedNodes.size(); ++local)
    {
      const NodePosition nodePos = relatedNodes[local];
      const NodeHandle node = nodes[nodePos];
      const FlatPosition flat = connectivity.findElementFlatIndex(zonePos, nodePos);

      weights.flatValue(flat) = static_cast<double>(zone.id() + node.id());
      nodeIdSum += node.id();
      SLIC_INFO("  node position " << nodePos << " -> " << node << ", weight "
                                   << weights.flatValue(flat));
    }
  }

  SLIC_ASSERT(zoneToNode.isValid(true));
  SLIC_ASSERT(connectivity.isValid(true));
  SLIC_ASSERT(nodeIdSum == 150);

  // Flat bivariate iteration yields the correctly typed endpoint coordinates.
  for(const auto coordinate : connectivity)
  {
    const ZoneHandle zone = zones[coordinate.first];
    const NodeHandle node = nodes[coordinate.second];
    SLIC_INFO("Connectivity pair: " << zone << " -> " << node);
  }

  return 0;
}
