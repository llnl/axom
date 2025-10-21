// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_DISPATCH_UNIFORM_TOPOLOGY_HPP_
#define AXOM_BUMP_DISPATCH_UNIFORM_TOPOLOGY_HPP_

#include "axom/bump/views/StructuredTopologyView.hpp"
#include "axom/bump/views/dispatch_utilities.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

namespace axom
{
namespace bump
{
namespace views
{
/*!
 * \brief Base template for uniform topology creation
 */
template <int NDIMS>
struct make_uniform_topology
{ };

/*!
 * \brief Create a 3D structured topology view with normal structured indexing.
 */
template <>
struct make_uniform_topology<3>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 3>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /*!
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    verify(topo, "topology");
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const conduit::Node &n_dims = coordset->fetch_existing("dims");
    LogicalIndex zoneDims;
    zoneDims[0] = n_dims[0].to_index_t() - 1;
    zoneDims[1] = n_dims[1].to_index_t() - 1;
    zoneDims[2] = n_dims[2].to_index_t() - 1;
    return Indexing(zoneDims);
  }

  /*!
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo) { return TopoView(indexing(topo)); }
};

/*!
 * \brief Create a 2D structured topology view with normal structured indexing.
 */
template <>
struct make_uniform_topology<2>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /*!
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    verify(topo, "topology");
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const conduit::Node &n_dims = coordset->fetch_existing("dims");
    LogicalIndex zoneDims;
    zoneDims[0] = n_dims[0].to_index_t() - 1;
    zoneDims[1] = n_dims[1].to_index_t() - 1;
    return Indexing(zoneDims);
  }

  /*!
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo) { return TopoView(indexing(topo)); }
};

/*!
 * \brief Create a 1D structured topology view with normal structured indexing.
 */
template <>
struct make_uniform_topology<1>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 1>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /*!
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    verify(topo, "topology");
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const conduit::Node &n_dims = coordset->fetch_existing("dims");
    LogicalIndex zoneDims;
    zoneDims[0] = n_dims[0].to_index_t() - 1;
    return Indexing(zoneDims);
  }

  /*!
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo) { return TopoView(indexing(topo)); }
};

namespace internal
{
/*!
 * \brief Base template for dispatching uniform topology.
 */
template <bool enabled, int NDIMS, typename FuncType>
struct dispatch_one_uniform_topology
{
  static void execute(const conduit::Node &AXOM_UNUSED_PARAM(topo),
                      FuncType &&AXOM_UNUSED_PARAM(func))
  { }
};

/*!
 * \brief Partial specialization to dispatch 3D uniform topology.
 */
template <typename FuncType>
struct dispatch_one_uniform_topology<true, 3, FuncType>
{
  /*!
   * \brief Make a proper view type for the uniform topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    auto topoView = make_uniform_topology<3>::view(topo);
    const std::string shape("hex");
    func(shape, topoView);
  }
};

/*!
 * \brief Partial specialization to dispatch 2D uniform topology.
 */
template <typename FuncType>
struct dispatch_one_uniform_topology<true, 2, FuncType>
{
  /*!
   * \brief Make a proper view type for the uniform topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    auto topoView = make_uniform_topology<2>::view(topo);
    const std::string shape("quad");
    func(shape, topoView);
  }
};

/*!
 * \brief Partial specialization to dispatch 1D uniform topology.
 */
template <typename FuncType>
struct dispatch_one_uniform_topology<true, 1, FuncType>
{
  /*!
   * \brief Make a proper view type for the uniform topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    auto topoView = make_uniform_topology<1>::view(topo);
    const std::string shape("line");
    func(shape, topoView);
  }
};

}  // end namespace internal

/*!
 * \brief Creates a topology view compatible with uniform topologies and passes that view to the supplied function.
 *
 * \tparam FuncType            The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set.
 *
 * \param topo The node that contains the topology.
 * \param func The function to invoke using the view.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_uniform_topology(const conduit::Node &topo, FuncType &&func)
{
  verify(topo, "topology");
  const conduit::Node *coordset =
    conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
  SLIC_ASSERT(coordset != nullptr);
  const conduit::Node &n_dims = coordset->fetch_existing("dims");
  switch(n_dims.dtype().number_of_elements())
  {
  case 3:
    internal::dispatch_one_uniform_topology<dimension_selected(SelectedDimensions, 3), 3, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  case 2:
    internal::dispatch_one_uniform_topology<dimension_selected(SelectedDimensions, 2), 2, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  case 1:
    internal::dispatch_one_uniform_topology<dimension_selected(SelectedDimensions, 1), 1, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  default:
    break;
  }
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
