// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MapBase.hpp
 *
 * \brief Contains an Abstract class MapBase
 *
 */

#ifndef SLAM_MAPBASE_HPP_
#define SLAM_MAPBASE_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

#include "axom/slam/Set.hpp"

namespace axom
{
namespace slam
{
/**
 * \class   MapBase
 *
 * \brief   A base class specifying methods that support operations of a Map,
 *          associating value(s) to each element in a Set. MapBase can be used
 *          as a base class pointer to a templated Map object.
 * \see     Map
 *
 */

template <typename SetPositionType = slam::DefaultPositionType>
class MapBase
{
public:
  using SetPosition = SetPositionType;

public:
  AXOM_HOST_DEVICE
  virtual ~MapBase() { };

  /**
   * \brief Get the number of entities in the set used by this map
   * \return The number of entities in the set used in the map.
   */
  AXOM_HOST_DEVICE virtual SetPosition size() const = 0;

  /**
   * \brief Checks whether the Map is valid.
   * \return   True if valid, false otherwise.
   */
  virtual bool isValid(bool verboseOutput) const = 0;

private:
  /**
   * \brief Verifies that the provided SetPosition is in a valid range.
   */
  virtual void verifyPosition(SetPosition) const = 0;
};

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_MAPBASE_HPP_
