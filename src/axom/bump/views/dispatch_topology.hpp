// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_DISPATCH_TOPOLOGY_HPP_
#define AXOM_BUMP_DISPATCH_TOPOLOGY_HPP_

#include "axom/bump/views/StructuredTopologyView.hpp"
#include "axom/bump/views/dispatch_uniform_topology.hpp"
#include "axom/bump/views/dispatch_rectilinear_topology.hpp"
#include "axom/bump/views/dispatch_structured_topology.hpp"
#include "axom/bump/views/dispatch_unstructured_topology.hpp"
#include "axom/bump/views/dispatch_utilities.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
namespace views
{
/*!
 * \brief Creates a topology view and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 * \tparam ShapeTypes A bitset containing the shape types that will be supported for unstructured.
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), int ShapeTypes = AnyShape, typename FuncType>
void dispatch_topology(const conduit::Node &topo, FuncType &&func)
{
  verify(topo, "topology");
  const auto type = topo.fetch_existing("type").as_string();

  if(type == "unstructured")
  {
    dispatch_unstructured_topology<ShapeTypes>(topo, func);
  }
  else if(type == "uniform" || type == "rectilinear" || type == "structured")
  {
    // Dispatch the structured topologies together because we can be smarter about the views.
    dispatch_structured_topologies<SelectedDimensions>(topo, func);
  }
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
