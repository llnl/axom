// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_SetIfacePolicies_HPP
#define SLAM_SetIfacePolicies_HPP

#include <type_traits>

#include "axom/slam/Set.hpp"
#include "axom/slam/policies/InterfacePolicies.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
namespace detail
{
/*!
 * \brief Templated helper class to select the correct base class to inherit
 *  for instances of OrderedSet.
 *
 *  This class helps select a virtual or non-virtual interface depending on
 *  the interface type:
 *  * In the ConcreteInterface case, we use ConcreteInterface as an empty
 *    base class to avoid virtual function calls in the derived instance.
 *  * In the VirtualInterface case, the base class is Set, which contains a
 *    common virtual interface between instances of OrderedSet.
 *
 * \tparam InterfacePolicy the interface type to select (virtual or concrete)
 *
 * \see InterfacePolicies.hpp
 */
template <typename InterfacePolicy, typename PosType, typename ElemType>
struct SetInterfaceSelector;

template <typename InterfacePolicy, typename PosType, typename ElemType>
struct SetInterfaceSelector
{
  static_assert(std::is_same<InterfacePolicy, ConcreteInterface>::value,
                "InterfacePolicy must be one of policies::ConcreteInterface or "
                "policies::VirtualInterface.");
  using Type = ConcreteInterface;
};

template <typename PosType, typename ElemType>
struct SetInterfaceSelector<VirtualInterface, PosType, ElemType>
{
  using Type = Set<PosType, ElemType>;
};
}  // namespace detail

template <typename InterfacePolicy, typename PosType, typename ElemType>
using SetInterface = typename detail::SetInterfaceSelector<InterfacePolicy, PosType, ElemType>::Type;

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_SetIfacePolicies_HPP
