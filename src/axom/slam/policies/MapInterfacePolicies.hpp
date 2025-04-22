// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_MapIfacePolicies_HPP
#define SLAM_MapIfacePolicies_HPP

#include <type_traits>

#include "axom/slam/Set.hpp"
#include "axom/slam/MapBase.hpp"
#include "axom/slam/policies/InterfacePolicies.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
namespace detail
{
template <typename InterfacePolicy, typename PosType>
struct MapInterfaceSelector
{
  static_assert(std::is_same<InterfacePolicy, ConcreteInterface>::value,
                "InterfacePolicy must be one of policies::ConcreteInterface or "
                "policies::VirtualInterface.");
  using Type = ConcreteInterface;
};

template <typename PosType>
struct MapInterfaceSelector<VirtualInterface, PosType>
{
  using Type = MapBase<PosType>;
};
}  // namespace detail

template <typename InterfacePolicy, typename SetPositionType>
using MapInterface = typename detail::MapInterfaceSelector<InterfacePolicy, SetPositionType>::Type;

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_MapIfacePolicies_HPP
