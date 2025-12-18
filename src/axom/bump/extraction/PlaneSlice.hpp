// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internal.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_PLANE_SLICE_HPP_
#define AXOM_BUMP_PLANE_SLICE_HPP_

#include "axom/bump/extraction/TableBasedExtractor.hpp"
#include "axom/bump/extraction/CutTableManager.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{
/*!
 * \brief This class slices a topology using a plane and puts the new topology into a new Conduit node.
 *
 * \tparam ExecSpace    The execution space where the compute-heavy kernels run.
 * \tparam TopologyView The topology view that can operate on the Blueprint topology.
 * \tparam CoordsetView The coordset view that can operate on the Blueprint coordset.
 * \tparam IntersectPolicy The intersector policy that can helps with cases and weights.
 * \tparam NamingPolicy The policy for making names from arrays of ids.
 */
template <typename ExecSpace,
          typename TopologyView,
          typename CoordsetView,
          typename IntersectPolicy = axom::bump::extraction::PlaneIntersector<TopologyView, CoordsetView>,
          typename NamingPolicy = axom::bump::HashNaming<axom::IndexType>>
using PlaneSlice =
  TableBasedExtractor<ExecSpace, CutTableManager, TopologyView, CoordsetView, IntersectPolicy, NamingPolicy, false>;

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
