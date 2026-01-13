// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_DISPATCH_UNSTRUCTURED_TOPOLOGY_HPP_
#define AXOM_BUMP_DISPATCH_UNSTRUCTURED_TOPOLOGY_HPP_

#include "axom/core.hpp"
#include "axom/bump/views/UnstructuredTopologySingleShapeView.hpp"
#include "axom/bump/views/UnstructuredTopologyPolyhedralView.hpp"
#include "axom/bump/views/UnstructuredTopologyMixedShapeView.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/views/Shapes.hpp"
#include "axom/bump/views/dispatch_utilities.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"

#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace bump
{
namespace views
{
// Turn on all bits so all shapes will be enabled.
constexpr int AnyShape = -1;

/*!
 * \brief This struct instantiates a topology view for an unstructured mesh that
 *        contains a single shape.
 *
 * \tparam ShapeType A shape class such as axom::bump::views::TriShape<int>
 */
template <typename ShapeType>
struct make_unstructured_single_shape_topology
{
  using TopologyView = UnstructuredTopologySingleShapeView<ShapeType>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

  /*!
   * \brief Instantiate a topology view from a Conduit node.
   *
   * \param n_topo The node that contains the topology.
   *
   * \return The topology view that wraps the Conduit data.
   */
  static TopologyView view(const conduit::Node &n_topo)
  {
    namespace utils = axom::bump::utilities;
    verify(n_topo, "topology");

    const std::string shape = n_topo["elements/shape"].as_string();
    SLIC_ERROR_IF(n_topo["type"].as_string() != "unstructured", "Type must be unstructured");
    SLIC_ERROR_IF(shape != ShapeType::name(), "Incompatible shape type");

    // Connectivity must exist.
    auto connView = utils::make_array_view<ConnectivityType>(n_topo["elements/connectivity"]);

    // Use sizes and offsets if they are present.
    axom::ArrayView<ConnectivityType> sizesView, offsetsView;
    if(n_topo.has_path("elements/sizes") && n_topo.has_path("elements/offsets"))
    {
      sizesView = utils::make_array_view<ConnectivityType>(n_topo["elements/sizes"]);
      offsetsView = utils::make_array_view<ConnectivityType>(n_topo["elements/offsets"]);
      SLIC_ASSERT(sizesView.size() == offsetsView.size());
    }
    else
    {
      // Variable-size shapes must have specified sizes, offsets.
      SLIC_ERROR_IF(ShapeType::is_variable_size(),
                    "A variable-size shape was provided without sizes of offsets.");
    }
    return TopologyView(connView, sizesView, offsetsView);
  }
};

/*!
 * \brief This struct instantiates a topology view for an unstructured mesh that
 *        contains polyhedra.
 *
 * \tparam ConnType An integer type.
 */
template <typename ConnType>
struct make_unstructured_polyhedral_topology
{
  using TopologyView = UnstructuredTopologyPolyhedralView<ConnType>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

  /*!
   * \brief Instantiate a topology view from a Conduit node.
   *
   * \param n_topo The node that contains the topology.
   *
   * \return The topology view that wraps the Conduit data.
   */
  static TopologyView view(const conduit::Node &n_topo)
  {
    namespace utils = axom::bump::utilities;
    verify(n_topo, "topology");

    const std::string shape = n_topo["elements/shape"].as_string();
    SLIC_ERROR_IF(n_topo["type"].as_string() != "unstructured", "Type must be unstructured");
    SLIC_ERROR_IF(shape != TopologyView::ShapeType::name(), "Incompatible shape type");

    // _bump_views_ph_topoview_begin
    auto topoView =
      TopologyView(utils::make_array_view<ConnectivityType>(n_topo["subelements/connectivity"]),
                   utils::make_array_view<ConnectivityType>(n_topo["subelements/sizes"]),
                   utils::make_array_view<ConnectivityType>(n_topo["subelements/offsets"]),
                   utils::make_array_view<ConnectivityType>(n_topo["elements/connectivity"]),
                   utils::make_array_view<ConnectivityType>(n_topo["elements/sizes"]),
                   utils::make_array_view<ConnectivityType>(n_topo["elements/offsets"]));
    // _bump_views_ph_topoview_end
    return topoView;
  }
};

/*!
 * \brief This function dispatches a Conduit polyhedral unstructured topology.
 *
 * \tparam FuncType The function/lambda type that will be invoked on the view.
 *
 * \param topo The node that contains the topology.
 * \param func The function/lambda to call with the topology view.
 */
///@{
template <typename FuncType>
void dispatch_unstructured_polyhedral_topology(const conduit::Node &topo, FuncType &&func)
{
  verify(topo, "topology");
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "polyhedral")
  {
    IndexNode_to_ArrayView_same(  //
      topo["subelements/connectivity"],
      topo["subelements/sizes"],
      topo["subelements/offsets"],
      topo["elements/connectivity"],
      topo["elements/sizes"],
      topo["elements/offsets"],
      [&](auto seConnView,
          auto seSizesView,
          auto seOffsetsView,
          auto connView,
          auto sizesView,
          auto offsetsView) {
        using ConnType = typename decltype(seConnView)::value_type;
        UnstructuredTopologyPolyhedralView<ConnType> ugView(seConnView,
                                                            seSizesView,
                                                            seOffsetsView,
                                                            connView,
                                                            sizesView,
                                                            offsetsView);
        func(shape, ugView);
      });
  }
}

template <typename ConnType, typename FuncType>
void typed_dispatch_unstructured_polyhedral_topology(const conduit::Node &topo, FuncType &&func)
{
  namespace utils = axom::bump::utilities;
  verify(topo, "topology");
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "polyhedral")
  {
    auto ugView = make_unstructured_polyhedral_topology<ConnType>::view(topo);
    func(shape, ugView);
  }
}
///@}

/*!
 * \brief This function dispatches a Conduit mixed unstructured topology.
 *
 * \tparam FuncType The function/lambda type that will be invoked on the view.
 *
 * \param topo The node that contains the topology.
 * \param func The function/lambda to call with the topology view.
 *
 * \note When this function makes the view, the view keeps a reference to
 *       the shape_map within the topology so we can build our own shape map
 *       later in the for_all_zones method.
 */
///@{
template <typename FuncType>
void dispatch_unstructured_mixed_topology(const conduit::Node &topo, FuncType &&func)
{
  verify(topo, "topology");
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "mixed")
  {
    IndexNode_to_ArrayView_same(
      topo["elements/connectivity"],
      topo["elements/shapes"],
      topo["elements/sizes"],
      topo["elements/offsets"],
      [&](auto connView, auto shapesView, auto sizesView, auto offsetsView) {
        using ConnType = typename decltype(connView)::value_type;

        // Get the allocator that allocated the connectivity. The shape map data
        // need to go into the same memory space.
        const int allocatorID =
          axom::getAllocatorIDFromPointer(topo["elements/connectivity"].data_ptr());

        // Make the shape map.
        axom::Array<IndexType> values, ids;
        auto shapeMap = buildShapeMap(topo, values, ids, allocatorID);

        UnstructuredTopologyMixedShapeView<ConnType> ugView(connView,
                                                            shapesView,
                                                            sizesView,
                                                            offsetsView,
                                                            shapeMap);
        func(shape, ugView);
      });
  }
}

template <typename ConnType, typename FuncType>
void typed_dispatch_unstructured_mixed_topology(const conduit::Node &topo, FuncType &&func)
{
  namespace utils = axom::bump::utilities;
  verify(topo, "topology");
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "mixed")
  {
    auto connView = utils::make_array_view<ConnType>(topo["elements/connectivity"]);
    auto shapesView = utils::make_array_view<ConnType>(topo["elements/shapes"]);
    auto sizesView = utils::make_array_view<ConnType>(topo["elements/sizes"]);
    auto offsetsView = utils::make_array_view<ConnType>(topo["elements/offsets"]);

    // Get the allocator that allocated the connectivity. The shape map data
    // need to go into the same memory space.
    volatile int allocatorID =
      axom::getAllocatorIDFromPointer(topo["elements/connectivity"].data_ptr());

    // Make the shape map.
    axom::Array<IndexType> values, ids;
    auto shapeMap = buildShapeMap(topo, values, ids, allocatorID);

    UnstructuredTopologyMixedShapeView<ConnType> ugView(connView,
                                                        shapesView,
                                                        sizesView,
                                                        offsetsView,
                                                        shapeMap);
    func(shape, ugView);
  }
}
///@}

template <typename... Args>
constexpr int encode_shapes(Args... args)
{
  return (... | args);
}

/*!
 * \brief This function turns a list of shapeID values into a bitfield that
 *        encodes the shapes. We use this in templating to limit which
 *        shapes get supported in dispatch instantiation.
 *
 * \param args A template parameter pack that contains ShapeID values.
 *
 * \return An integer that encodes the shape ids.
 */
template <typename... Args>
constexpr int select_shapes(Args... args)
{
  return encode_shapes((1 << args)...);
}

//------------------------------------------------------------------------------
namespace internal
{
/*!
 * \brief Base template for dispatching various shapes conditionally.
 */
template <bool enabled, typename ConnType, typename ShapeType, typename FuncType>
struct dispatch_shape
{
  /*!
   * \brief Execute method that gets generated when a shape is not enabled or supported. Do nothing.
   */
  static void execute(bool &AXOM_UNUSED_PARAM(eligible),
                      const std::string &AXOM_UNUSED_PARAM(shape),
                      const axom::ArrayView<ConnType> &AXOM_UNUSED_PARAM(connView),
                      const axom::ArrayView<ConnType> &AXOM_UNUSED_PARAM(sizesView),
                      const axom::ArrayView<ConnType> &AXOM_UNUSED_PARAM(offsetsView),
                      FuncType &&AXOM_UNUSED_PARAM(func))
  { }

  /*!
   * \brief Execute method that gets generated when a shape is not enabled or supported. Do nothing.
   */
  static void execute(bool &AXOM_UNUSED_PARAM(eligible),
                      const std::string &AXOM_UNUSED_PARAM(shape),
                      const conduit::Node &AXOM_UNUSED_PARAM(topo),
                      FuncType &&AXOM_UNUSED_PARAM(func))
  { }
};

// Partial specializations that make views for various shape types.

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, TriShape<ConnType>, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const axom::ArrayView<ConnType> &connView,
                      const axom::ArrayView<ConnType> &sizesView,
                      const axom::ArrayView<ConnType> &offsetsView,
                      FuncType &&func)
  {
    if(eligible && shape == "tri")
    {
      UnstructuredTopologySingleShapeView<TriShape<ConnType>> ugView(connView, sizesView, offsetsView);
      func(shape, ugView);
      eligible = false;
    }
  }
};

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, QuadShape<ConnType>, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const axom::ArrayView<ConnType> &connView,
                      const axom::ArrayView<ConnType> &sizesView,
                      const axom::ArrayView<ConnType> &offsetsView,
                      FuncType &&func)
  {
    if(eligible && shape == "quad")
    {
      UnstructuredTopologySingleShapeView<QuadShape<ConnType>> ugView(connView, sizesView, offsetsView);
      func(shape, ugView);
      eligible = false;
    }
  }
};

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, PolygonShape<ConnType>, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const axom::ArrayView<ConnType> &connView,
                      const axom::ArrayView<ConnType> &sizesView,
                      const axom::ArrayView<ConnType> &offsetsView,
                      FuncType &&func)
  {
    if(eligible && shape == "polygonal")
    {
      UnstructuredTopologySingleShapeView<PolygonShape<ConnType>> ugView(connView,
                                                                         sizesView,
                                                                         offsetsView);
      func(shape, ugView);
      eligible = false;
    }
  }
};

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, TetShape<ConnType>, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const axom::ArrayView<ConnType> &connView,
                      const axom::ArrayView<ConnType> &sizesView,
                      const axom::ArrayView<ConnType> &offsetsView,
                      FuncType &&func)
  {
    if(eligible && shape == "tet")
    {
      UnstructuredTopologySingleShapeView<TetShape<ConnType>> ugView(connView, sizesView, offsetsView);
      func(shape, ugView);
      eligible = false;
    }
  }
};

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, PyramidShape<ConnType>, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const axom::ArrayView<ConnType> &connView,
                      const axom::ArrayView<ConnType> &sizesView,
                      const axom::ArrayView<ConnType> &offsetsView,
                      FuncType &&func)
  {
    if(eligible && shape == "pyramid")
    {
      UnstructuredTopologySingleShapeView<PyramidShape<ConnType>> ugView(connView,
                                                                         sizesView,
                                                                         offsetsView);
      func(shape, ugView);
      eligible = false;
    }
  }
};

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, WedgeShape<ConnType>, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const axom::ArrayView<ConnType> &connView,
                      const axom::ArrayView<ConnType> &sizesView,
                      const axom::ArrayView<ConnType> &offsetsView,
                      FuncType &&func)
  {
    if(eligible && shape == "wedge")
    {
      UnstructuredTopologySingleShapeView<WedgeShape<ConnType>> ugView(connView,
                                                                       sizesView,
                                                                       offsetsView);
      func(shape, ugView);
      eligible = false;
    }
  }
};

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, HexShape<ConnType>, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const axom::ArrayView<ConnType> &connView,
                      const axom::ArrayView<ConnType> &sizesView,
                      const axom::ArrayView<ConnType> &offsetsView,
                      FuncType &&func)
  {
    if(eligible && shape == "hex")
    {
      UnstructuredTopologySingleShapeView<HexShape<ConnType>> ugView(connView, sizesView, offsetsView);
      func(shape, ugView);
      eligible = false;
    }
  }
};

struct SelectMixedShape
{ };

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, SelectMixedShape, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const conduit::Node &topo,
                      FuncType &&func)
  {
    if(eligible && shape == "mixed")
    {
      typed_dispatch_unstructured_mixed_topology<ConnType>(topo, std::forward<FuncType>(func));
      eligible = false;
    }
  }
};

struct SelectPHShape
{ };

template <typename ConnType, typename FuncType>
struct dispatch_shape<true, ConnType, SelectPHShape, FuncType>
{
  static void execute(bool &eligible,
                      const std::string &shape,
                      const conduit::Node &topo,
                      FuncType &&func)
  {
    if(eligible && shape == "polyhedral")
    {
      typed_dispatch_unstructured_polyhedral_topology<ConnType>(topo, std::forward<FuncType>(func));
      eligible = false;
    }
  }
};

}  // end namespace internal
//------------------------------------------------------------------------------

/*!
 * \brief This function dispatches a Conduit topology to the right view type
 *        and passes that view to the supplied function/lambda.
 *
 * \tparam ShapeTypes Allows us to limit which shape types get compiled in.
 * \tparam FuncType The function/lambda type that will be invoked on the view.
 *
 * \param topo The node that contains the topology.
 * \param func The function/lambda to call with the topology view.
 */
template <typename ConnType, int ShapeTypes = AnyShape, typename FuncType>
void typed_dispatch_unstructured_topology(const conduit::Node &topo, FuncType &&func)
{
  namespace utils = axom::bump::utilities;
  verify(topo, "topology");
  const std::string type = topo["type"].as_string();
  if(type == "unstructured")
  {
    const std::string shape = topo["elements/shape"].as_string();
    const auto connView = utils::make_array_view<ConnType>(topo["elements/connectivity"]);
    bool eligible = true;

    // Conditionally add polyhedron support.
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Polyhedron_ShapeID),
                             ConnType,
                             internal::SelectPHShape,
                             FuncType>::execute(eligible, shape, topo, std::forward<FuncType>(func));

    // Conditionally add mixed shape support.
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Mixed_ShapeID),
                             ConnType,
                             internal::SelectMixedShape,
                             FuncType>::execute(eligible, shape, topo, std::forward<FuncType>(func));

    // Make sizes / offsets views if the values are present.
    axom::ArrayView<ConnType> sizesView, offsetsView;
    if(topo.has_path("elements/sizes"))
      sizesView = utils::make_array_view<ConnType>(topo.fetch_existing("elements/sizes"));
    if(topo.has_path("elements/offsets"))
      offsetsView = utils::make_array_view<ConnType>(topo.fetch_existing("elements/offsets"));

    // Conditionally add support for other shapes.
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Tri_ShapeID),
                             ConnType,
                             TriShape<ConnType>,
                             FuncType>::execute(eligible,
                                                shape,
                                                connView,
                                                sizesView,
                                                offsetsView,
                                                std::forward<FuncType>(func));
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Quad_ShapeID),
                             ConnType,
                             QuadShape<ConnType>,
                             FuncType>::execute(eligible,
                                                shape,
                                                connView,
                                                sizesView,
                                                offsetsView,
                                                std::forward<FuncType>(func));
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Polygon_ShapeID),
                             ConnType,
                             PolygonShape<ConnType>,
                             FuncType>::execute(eligible,
                                                shape,
                                                connView,
                                                sizesView,
                                                offsetsView,
                                                std::forward<FuncType>(func));
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Tet_ShapeID),
                             ConnType,
                             TetShape<ConnType>,
                             FuncType>::execute(eligible,
                                                shape,
                                                connView,
                                                sizesView,
                                                offsetsView,
                                                std::forward<FuncType>(func));
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Pyramid_ShapeID),
                             ConnType,
                             PyramidShape<ConnType>,
                             FuncType>::execute(eligible,
                                                shape,
                                                connView,
                                                sizesView,
                                                offsetsView,
                                                std::forward<FuncType>(func));
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Wedge_ShapeID),
                             ConnType,
                             WedgeShape<ConnType>,
                             FuncType>::execute(eligible,
                                                shape,
                                                connView,
                                                sizesView,
                                                offsetsView,
                                                std::forward<FuncType>(func));
    internal::dispatch_shape<axom::utilities::bitIsSet(ShapeTypes, Hex_ShapeID),
                             ConnType,
                             HexShape<ConnType>,
                             FuncType>::execute(eligible,
                                                shape,
                                                connView,
                                                sizesView,
                                                offsetsView,
                                                std::forward<FuncType>(func));

    // TODO: points, lines, polygon
  }
}

/// Dispatch in a way that does not care about the connectivity type.
template <int ShapeTypes = AnyShape, typename FuncType>
void dispatch_unstructured_topology(const conduit::Node &topo, FuncType &&func)
{
  verify(topo, "topology");
  IndexNode_to_ArrayView(topo["elements/connectivity"], [&](auto connView) {
    using ConnType = typename decltype(connView)::value_type;
    typed_dispatch_unstructured_topology<ConnType, ShapeTypes>(topo, func);
  });
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
