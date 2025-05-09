// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file AllNearestNeighbors.hpp
 * \brief Defines all-nearest-neighbor queries
 */

#ifndef AXOM_QUEST_ALL_NEAREST_NEIGHBORS_HPP_
#define AXOM_QUEST_ALL_NEAREST_NEIGHBORS_HPP_

namespace axom
{
namespace quest
{
enum SearchStatus
{
  NEIGHBOR_NOT_FOUND = -1
};

/// \name Nearest Neighbor query
/// @{

/*!
 * \brief Given a list of point locations and regions, for each point, find
 *   the closest point in a different region within a given search radius.
 * \param [in] x X-coordinates of input points
 * \param [in] y Y-coordinates of input points
 * \param [in] z Z-coordinates of input points
 * \param [in] region Region of each point
 * \param [in] n Number of points
 * \param [in] limit Max distance for all-nearest-neighbors query
 * \param [out] neighbor Index of nearest neighbor not in the same class
 *    (or NEIGHBOR_NOT_FOUND)
 * \param [out] sqdistance Squared distance to nearest neighbor
 * \pre x, y, z, and region have n entries
 * \pre neighbor is allocated with room for n entries
 *
 * This method inserts all points p at (x[i], y[i], z[i]) into a UniformGrid
 * index. Then for each point p, it gets the UniformGrid bins that overlap
 * the box (p - (limit, limit, limit), p + (limit, limit, limit).  The method
 * compares p to each point in this list of bins and returns the index of the
 * closest point.
 *
 * We expect the use of the UniformGrid  will result in a substantial time
 * savings over a brute-force all-to-all algorithm, but the query's run time
 * is dependent on the point distribution.
 */
void all_nearest_neighbors(const double* x,
                           const double* y,
                           const double* z,
                           const int* region,
                           int n,
                           double limit,
                           int* neighbor,
                           double* sqdistance);

/// @}

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_ALL_NEAREST_NEIGHBORS_HPP_
