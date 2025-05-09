// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_InterfacePolicies_HPP
#define SLAM_InterfacePolicies_HPP

/**
 * \file InterfacePolicies.hpp
 *
 * \brief Interface policies for SLAM
 *
 * Interface policies, when used to instantiate supported SLAM types, allow for
 * selecting between a virtual interface or a concrete interface.
 *
 * Using a virtual interface has the advantage of being able to switch at
 * runtime between different types that correspond to the same base interface.
 * However, there may be some disadvantages:
 *  * Virtual function calls may be slower, if performance is a concern.
 *  * For CUDA/HIP code, virtual SLAM objects cannot cross a host/device
 *    boundary, since the virtual function table is only correct for either the
 *    host or the device, depending on where the object was constructed.
 *
 * Using a concrete interface trades the convenience of runtime polymorphism
 * for improved performance by avoiding virtual function calls. For CUDA/HIP
 * code, SLAM types instantiated with a concrete interface are usually able to
 * cross a host/device boundary, i.e. when captured in a RAJA lambda-based loop.
 *
 * \see SetInterfacePolicies.hpp
 * \set MapInterfacePolicies.hpp
 * \see BivariateSetInterfacePolicies.hpp
 */

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \class VirtualInterface
 *
 * \brief Policy to use a virtual interface with a given Slam type.
 */
struct VirtualInterface
{ };

/**
 * \class VirtualInterface
 *
 * \brief Policy to use a concrete, CRTP-based interface with a given Slam type.
 */
struct ConcreteInterface
{ };

}  // namespace policies
}  // namespace slam
}  // namespace axom

#endif  // SLAM_InterfacePolicies_HPP
