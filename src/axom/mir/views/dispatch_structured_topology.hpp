// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/StructuredIndexing.hpp"
#include "axom/mir/views/StridedStructuredIndexing.hpp"
#include "axom/mir/views/dispatch_utilities.hpp"
#include "axom/mir/views/dispatch_uniform_topology.hpp"
#include "axom/mir/views/dispatch_rectilinear_topology.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"

#include <conduit/conduit.hpp>

#include <utility>

namespace axom
{
namespace mir
{
namespace views
{
/*!
 * \brief Base template for strided structured topology creation
 */
template <int NDIMS>
struct make_strided_structured
{ };

/*!
 * \brief Create a 3D structured topology view with strided structured indexing.
 */
template <>
struct make_strided_structured<3>
{
  using Indexing = views::StridedStructuredIndexing<axom::IndexType, 3>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /*!
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const std::string offsetsKey("elements/dims/offsets");
    const std::string stridesKey("elements/dims/strides");

    LogicalIndex zoneDims;
    zoneDims[0] = topo.fetch_existing("elements/dims/i").to_int();
    zoneDims[1] = topo.fetch_existing("elements/dims/j").to_int();
    zoneDims[2] = topo.fetch_existing("elements/dims/k").to_int();

    LogicalIndex offsets {{0, 0, 0}}, strides {{1, 1, 1}};
    bputils::fillFromNode(topo, offsetsKey, offsets, true);
    if(bputils::fillFromNode(topo, stridesKey, strides, true))
    {
      // Make zone striding.
      strides[1]--;
      strides[2] = strides[1] * (zoneDims[1] - 1);
    }
    else
    {
      strides[1] = zoneDims[0];
      strides[2] = zoneDims[0] * zoneDims[1];
    }

    return Indexing(zoneDims, offsets, strides);
  }

  /*!
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo) { return TopoView(indexing(topo)); }
};

/*!
 * \brief Create a 2D structured topology view with strided structured indexing.
 */
template <>
struct make_strided_structured<2>
{
  using Indexing = views::StridedStructuredIndexing<axom::IndexType, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /*!
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const std::string offsetsKey("elements/dims/offsets");
    const std::string stridesKey("elements/dims/strides");
    LogicalIndex zoneDims;
    zoneDims[0] = topo.fetch_existing("elements/dims/i").to_int();
    zoneDims[1] = topo.fetch_existing("elements/dims/j").to_int();

    LogicalIndex offsets {{0, 0}}, strides {{1, 1}};
    bputils::fillFromNode(topo, offsetsKey, offsets, true);
    if(bputils::fillFromNode(topo, stridesKey, strides, true))
    {
      // Make zone striding.
      strides[1]--;
    }
    else
    {
      strides[1] = zoneDims[0];
    }

    return Indexing(zoneDims, offsets, strides);
  }

  /*!
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo) { return TopoView(indexing(topo)); }
};

/*!
 * \brief Create a 1D structured topology view with strided structured indexing.
 */
template <>
struct make_strided_structured<1>
{
  using Indexing = views::StridedStructuredIndexing<axom::IndexType, 1>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /*!
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const std::string offsetsKey("elements/dims/offsets");
    const std::string stridesKey("elements/dims/strides");

    LogicalIndex zoneDims;
    zoneDims[0] = topo.fetch_existing("elements/dims/i").to_int();

    LogicalIndex offsets {0}, strides {1};
    bputils::fillFromNode(topo, offsetsKey, offsets, true);
    bputils::fillFromNode(topo, stridesKey, strides, true);

    return Indexing(zoneDims, offsets, strides);
  }

  /*!
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo) { return TopoView(indexing(topo)); }
};

/*!
 * \brief Base template for structured topology creation
 */
template <int NDIMS>
struct make_structured
{ };

/*!
 * \brief Create a 3D structured topology view with normal structured indexing.
 */
template <>
struct make_structured<3>
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
    LogicalIndex zoneDims;
    zoneDims[0] = topo.fetch_existing("elements/dims/i").to_int();
    zoneDims[1] = topo.fetch_existing("elements/dims/j").to_int();
    zoneDims[2] = topo.fetch_existing("elements/dims/k").to_int();

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
struct make_structured<2>
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
    LogicalIndex zoneDims;
    zoneDims[0] = topo.fetch_existing("elements/dims/i").to_int();
    zoneDims[1] = topo.fetch_existing("elements/dims/j").to_int();
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
struct make_structured<1>
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
    LogicalIndex zoneDims;
    zoneDims[0] = topo.fetch_existing("elements/dims/i").to_int();

    return Indexing(zoneDims);
  }

  /*!
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo) { return TopoView(indexing(topo)); }
};

//------------------------------------------------------------------------------
namespace internal
{
/*!
 * \brief Base template for dispatching structured topology.
 */
template <bool enabled, int NDIMS, typename FuncType>
struct dispatch_only_structured_topology
{
  static void execute(const conduit::Node &AXOM_UNUSED_PARAM(topo),
                      FuncType &&AXOM_UNUSED_PARAM(func))
  { }
};

/*!
 * \brief Partial specialization to dispatch 3D structured topology.
 */
template <typename FuncType>
struct dispatch_only_structured_topology<true, 3, FuncType>
{
  /*!
   * \brief Make a proper view type for the structured topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    const std::string offsetsKey("elements/dims/offsets");
    const std::string stridesKey("elements/dims/strides");
    const std::string shape("hex");
    if(topo.has_path(offsetsKey) || topo.has_path(stridesKey))
    {
      auto topoView = make_strided_structured<3>::view(topo);
      func(shape, topoView);
    }
    else
    {
      auto topoView = make_structured<3>::view(topo);
      func(shape, topoView);
    }
  }
};

/*!
 * \brief Partial specialization to dispatch 2D structured topology.
 */
template <typename FuncType>
struct dispatch_only_structured_topology<true, 2, FuncType>
{
  /*!
   * \brief Make a proper view type for the structured topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    const std::string offsetsKey("elements/dims/offsets");
    const std::string stridesKey("elements/dims/strides");
    const std::string shape("quad");
    if(topo.has_path(offsetsKey) || topo.has_path(stridesKey))
    {
      auto topoView = make_strided_structured<2>::view(topo);
      func(shape, topoView);
    }
    else
    {
      auto topoView = make_structured<2>::view(topo);
      func(shape, topoView);
    }
  }
};

/*!
 * \brief Partial specialization to dispatch 1D structured topology.
 */
template <typename FuncType>
struct dispatch_only_structured_topology<true, 1, FuncType>
{
  /*!
   * \brief Make a proper view type for the structured topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    const std::string offsetsKey("elements/dims/offsets");
    const std::string stridesKey("elements/dims/strides");
    const std::string shape("line");
    if(topo.has_path(offsetsKey) || topo.has_path(stridesKey))
    {
      auto topoView = make_strided_structured<1>::view(topo);
      func(shape, topoView);
    }
    else
    {
      auto topoView = make_structured<1>::view(topo);
      func(shape, topoView);
    }
  }
};

//------------------------------------------------------------------------------
/*!
 * \brief Base template for dispatching structured topology.
 */
template <bool enabled, int NDIMS, typename FuncType>
struct dispatch_any_structured_topology
{
  static void execute(const conduit::Node &AXOM_UNUSED_PARAM(topo),
                      FuncType &&AXOM_UNUSED_PARAM(func))
  { }
};

/*!
 * \brief Partial specialization to dispatch 3D structured topologies.
 */
template <typename FuncType>
struct dispatch_any_structured_topology<true, 3, FuncType>
{
  /*!
   * \brief Make a proper view type for the structured topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    const std::string offsetsKey("offsets"), stridesKey("strides");
    const std::string type = topo.fetch_existing("type").as_string();
    const std::string shape("hex");

    if(type == "structured" && topo.has_path(offsetsKey) && topo.has_path(stridesKey))
    {
      auto topoView = make_strided_structured<3>::view(topo);
      func(shape, topoView);
    }
    else
    {
      // Make these topology types share the same func dispatch.
      StructuredTopologyView<StructuredIndexing<axom::IndexType, 3>> topoView;
      if(type == "uniform")
        topoView = make_uniform<3>::view(topo);
      else if(type == "rectilinear")
        topoView = make_rectilinear<3>::view(topo);
      else if(type == "structured")
        topoView = make_structured<3>::view(topo);
      func(shape, topoView);
    }
  }
};

/*!
 * \brief Partial specialization to dispatch 2D structured topologies.
 */
template <typename FuncType>
struct dispatch_any_structured_topology<true, 2, FuncType>
{
  /*!
   * \brief Make a proper view type for the structured topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    const std::string offsetsKey("offsets"), stridesKey("strides");
    const std::string type = topo.fetch_existing("type").as_string();
    const std::string shape("quad");
    if(type == "structured" && topo.has_path(offsetsKey) && topo.has_path(stridesKey))
    {
      auto topoView = make_strided_structured<2>::view(topo);
      func(shape, topoView);
    }
    else
    {
      // Make these topology types share the same func dispatch.
      StructuredTopologyView<StructuredIndexing<axom::IndexType, 2>> topoView;
      if(type == "uniform")
        topoView = make_uniform<2>::view(topo);
      else if(type == "rectilinear")
        topoView = make_rectilinear<2>::view(topo);
      else if(type == "structured")
        topoView = make_structured<2>::view(topo);
      func(shape, topoView);
    }
  }
};

/*!
 * \brief Partial specialization to dispatch 1D structured topologies.
 */
template <typename FuncType>
struct dispatch_any_structured_topology<true, 1, FuncType>
{
  /*!
   * \brief Make a proper view type for the structured topology and pass the
   *        view to the supplied kernel.
   *
   * \param topo The node that contains the topology.
   * \param func The kernel to be invoked.
   */
  static void execute(const conduit::Node &topo, FuncType &&func)
  {
    const std::string offsetsKey("offsets"), stridesKey("strides");
    const std::string type = topo.fetch_existing("type").as_string();
    const std::string shape("line");

    // Make these topology types share the same func dispatch.
    StructuredTopologyView<StructuredIndexing<axom::IndexType, 1>> topoView;
    if(type == "uniform")
      topoView = make_uniform<1>::view(topo);
    else if(type == "rectilinear")
      topoView = make_rectilinear<1>::view(topo);
    else if(type == "structured")
      topoView = make_structured<1>::view(topo);
    func(shape, topoView);
  }
};

}  // end namespace internal

/*!
 * \brief Creates a topology view compatible with structured topologies and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_structured_topology(const conduit::Node &topo, FuncType &&func)
{
  int ndims = 1;
  ndims += topo.has_path("elements/dims/j") ? 1 : 0;
  ndims += topo.has_path("elements/dims/k") ? 1 : 0;

  switch(ndims)
  {
  case 3:
    internal::dispatch_only_structured_topology<dimension_selected(SelectedDimensions, 3), 3, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  case 2:
    internal::dispatch_only_structured_topology<dimension_selected(SelectedDimensions, 2), 2, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  case 1:
    internal::dispatch_only_structured_topology<dimension_selected(SelectedDimensions, 1), 1, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  default:
    break;
  }
}

/*!
 * \brief Creates a topology view compatible with various logically "structured" topologies (uniform, rectilinear, structured) and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 *
 * \param topo     The node that contains the topology.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.
 *
 * \note We try to initialize the topoView for each dimension and share the dispatch.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_structured_topologies(const conduit::Node &topo, FuncType &&func)
{
  const auto ndims = conduit::blueprint::mesh::utils::topology::dims(topo);
  switch(ndims)
  {
  case 3:
    internal::dispatch_any_structured_topology<dimension_selected(SelectedDimensions, 3), 3, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  case 2:
    internal::dispatch_any_structured_topology<dimension_selected(SelectedDimensions, 2), 2, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  case 1:
    internal::dispatch_any_structured_topology<dimension_selected(SelectedDimensions, 1), 1, FuncType>::execute(
      topo,
      std::forward<FuncType>(func));
    break;
  default:
    break;
  }
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
