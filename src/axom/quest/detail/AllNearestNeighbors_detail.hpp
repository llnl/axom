// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_ALL_NEAREST_NEIGHBORS_DETAIL_HPP_
#define AXOM_QUEST_ALL_NEAREST_NEIGHBORS_DETAIL_HPP_

namespace axom
{
namespace quest
{
namespace detail
{
//------------------------------------------------------------------------------
inline double squared_distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;

  return dx * dx + dy * dy + dz * dz;
}

}  // end namespace detail
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_ALL_NEAREST_NEIGHBORS_DETAIL_HPP_
