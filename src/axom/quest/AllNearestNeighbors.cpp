// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/AllNearestNeighbors.hpp"
#include "axom/quest/detail/AllNearestNeighbors_detail.hpp"

#include <cfloat>  // for DBL_MAX

#include "axom/config.hpp"
#include "axom/spin/UniformGrid.hpp"

namespace axom
{
namespace quest
{
/* Given a list of point locations and regions, for each point, find
 * the closest point in a different region within a given search radius.
 */
void all_nearest_neighbors(const double* x,
                           const double* y,
                           const double* z,
                           const int* region,
                           int n,
                           double limit,
                           int* neighbor,
                           double* sqdistance)
{
  // Indexed approach.  For each point i, test distance to all other
  // points in this and neighboring UniformGrid bins (out to distance limit)
  // and report result.

  using GridType = spin::UniformGrid<int, 3>;
  using BoxType = GridType::BoxType;
  using PointType = GridType::PointType;
  using VectorType = BoxType::VectorType;

  double sqlimit = limit * limit;

  PointType pmin, pmax;
  for(int i = 0; i < n; ++i)
  {
    sqdistance[i] = DBL_MAX;
    neighbor[i] = NEIGHBOR_NOT_FOUND;
    pmin[0] = std::min(pmin[0], x[i]);
    pmax[0] = std::max(pmax[0], x[i]);
    pmin[1] = std::min(pmin[1], y[i]);
    pmax[1] = std::max(pmax[1], y[i]);
    pmin[2] = std::min(pmin[2], z[i]);
    pmax[2] = std::max(pmax[2], z[i]);
  }
  BoxType allpointsbox(pmin, pmax);

  int res[3];
  VectorType boxrange = allpointsbox.range();
  for(int i = 0; i < 3; ++i)
  {
    res[i] = std::max(1, (int)(boxrange[i] / limit + 0.5));
  }

  // 1. Build an index, inserting each point individually
  GridType ugrid(allpointsbox, res);
  for(int i = 0; i < n; ++i)
  {
    ugrid.insert(BoxType(PointType::make_point(x[i], y[i], z[i])), i);
  }

  // 2. For point a,
  for(int i = 0; i < n; ++i)
  {
    // 3. For each other bin B less than limit distance away from a, for each
    // point b in B,
    PointType qmin = PointType::make_point(x[i] - limit, y[i] - limit, z[i] - limit);
    PointType qmax = PointType::make_point(x[i] + limit, y[i] + limit, z[i] + limit);
    BoxType qbox(qmin, qmax);
    const std::vector<int> qbins = ugrid.getBinsForBbox(qbox);
    const size_t querybincount = qbins.size();
    for(size_t binidx = 0; binidx < querybincount; ++binidx)
    {
      const axom::ArrayView<int> bs = ugrid.getBinContents(qbins[binidx]);
      const size_t binsize = bs.size();
      for(size_t bj = 0; bj < binsize; ++bj)
      {
        // 4. Compare distances to find the closest distance d = |ab|
        int j = bs[bj];
        if(region[i] != region[j])
        {
          double sqdist = detail::squared_distance(x[i], y[i], z[i], x[j], y[j], z[j]);
          if(sqdist < sqdistance[i] && sqdist < sqlimit)
          {
            sqdistance[i] = sqdist;
            neighbor[i] = j;
          }
        }
      }
    }
  }
}

}  // end namespace quest
}  // end namespace axom
