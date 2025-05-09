// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_COORDSET_HPP_
#define AXOM_MIR_DISPATCH_COORDSET_HPP_

#include "axom/mir/views/ExplicitCoordsetView.hpp"
#include "axom/mir/views/UniformCoordsetView.hpp"
#include "axom/mir/views/RectilinearCoordsetView.hpp"
#include "axom/mir/views/NodeArrayView.hpp"

namespace axom
{
namespace mir
{
namespace views
{
/*!
 * \brief Base template for creating a rectilinear coordset view.
 */
template <typename DataType, int NDIMS>
struct make_rectilinear_coordset
{ };

/*!
 * \brief Partial specialization for creating 3D rectilinear coordset view.
 */
template <typename DataType>
struct make_rectilinear_coordset<DataType, 3>
{
  using CoordsetView = axom::mir::views::RectilinearCoordsetView3<DataType>;

  /*!
   * \brief Create the coordset view and initialize it from the coordset.
   * \param topo The node containing the coordset.
   * \return The coordset view.
   */
  static CoordsetView view(const conduit::Node &coordset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const conduit::Node &values = coordset.fetch_existing("values");
    auto xView = bputils::make_array_view<DataType>(values[0]);
    auto yView = bputils::make_array_view<DataType>(values[1]);
    auto zView = bputils::make_array_view<DataType>(values[2]);
    return CoordsetView(xView, yView, zView);
  }
};

/*!
 * \brief Partial specialization for creating 2D rectilinear coordset view.
 */
template <typename DataType>
struct make_rectilinear_coordset<DataType, 2>
{
  using CoordsetView = axom::mir::views::RectilinearCoordsetView2<DataType>;

  /*!
   * \brief Create the coordset view and initialize it from the coordset.
   * \param topo The node containing the coordset.
   * \return The coordset view.
   */
  static CoordsetView view(const conduit::Node &coordset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const conduit::Node &values = coordset.fetch_existing("values");
    auto xView = bputils::make_array_view<DataType>(values[0]);
    auto yView = bputils::make_array_view<DataType>(values[1]);
    return CoordsetView(xView, yView);
  }
};

/*!
 * \brief Base template for creating a rectilinear coordset view.
 */
template <int NDIMS>
struct make_uniform_coordset
{ };

/*!
 * \brief Partial specialization for creating 3D uniform coordset view.
 */
template <>
struct make_uniform_coordset<3>
{
  using CoordsetView = axom::mir::views::UniformCoordsetView<double, 3>;

  /*!
   * \brief Create the coordset view and initialize it from the coordset.
   * \param topo The node containing the coordset.
   * \return The coordset view.
   */
  static CoordsetView view(const conduit::Node &coordset)
  {
    const std::string keys[] = {"i", "j", "k"};
    const conduit::Node &n_dims = coordset["dims"];
    axom::StackArray<axom::IndexType, 3> dims;
    axom::StackArray<double, 3> origin {0., 0., 0.}, spacing {1., 1., 1.};
    for(int i = 0; i < 3; i++)
    {
      dims[i] = n_dims.fetch_existing(keys[i]).to_int();
      if(coordset.has_child("origin")) origin[i] = coordset["origin"][i].to_double();
      if(coordset.has_child("spacing")) spacing[i] = coordset["spacing"][i].to_double();
    }
    return CoordsetView(dims, origin, spacing);
  }
};

/*!
 * \brief Partial specialization for creating 2D rectilinear coordset view.
 */
template <>
struct make_uniform_coordset<2>
{
  using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;

  /*!
   * \brief Create the coordset view and initialize it from the coordset.
   * \param topo The node containing the coordset.
   * \return The coordset view.
   */
  static CoordsetView view(const conduit::Node &coordset)
  {
    const std::string keys[] = {"i", "j"};
    const conduit::Node &n_dims = coordset["dims"];
    axom::StackArray<axom::IndexType, 2> dims;
    axom::StackArray<double, 2> origin {0., 0.}, spacing {1., 1.};
    for(int i = 0; i < 2; i++)
    {
      dims[i] = n_dims.fetch_existing(keys[i]).to_int();
      if(coordset.has_child("origin")) origin[i] = coordset["origin"][i].to_double();
      if(coordset.has_child("spacing")) spacing[i] = coordset["spacing"][i].to_double();
    }
    return CoordsetView(dims, origin, spacing);
  }
};

/*!
 * \brief Dispatch an uniform coordset to a function.
 *
 * \tparam FuncType The type of the function / lambda to invoke. It is expected
 *                  that the callable accepts an auto argument for a coordset view.
 *
 * \param coordset The Conduit node that contains the coordset.
 * \param func     The function/lambda to invoke using the coordset view.
 */
template <typename FuncType>
void dispatch_uniform_coordset(const conduit::Node &coordset, FuncType &&func)
{
  const conduit::Node &n_dims = coordset["dims"];
  const conduit::index_t ndims = n_dims.number_of_children();
  if(ndims == 2)
  {
    auto coordsetView = make_uniform_coordset<2>::view(coordset);
    func(coordsetView);
  }
  else if(ndims == 3)
  {
    auto coordsetView = make_uniform_coordset<3>::view(coordset);
    func(coordsetView);
  }
}

/*!
 * \brief Dispatch a rectilinear coordset to a function.
 *
 * \tparam FuncType The type of the function / lambda to invoke. It is expected
 *                  that the callable accepts an auto argument for a coordset view.
 *
 * \param coordset The Conduit node that contains the coordset.
 * \param func     The function/lambda to invoke using the coordset view.
 */
template <typename FuncType>
void dispatch_rectilinear_coordset(const conduit::Node &coordset, FuncType &&func)
{
  const conduit::Node &values = coordset["values"];
  if(values.number_of_children() == 2)
  {
    axom::mir::views::FloatNode_to_ArrayView_same(values[0], values[1], [&](auto xView, auto yView) {
      RectilinearCoordsetView2<typename decltype(xView)::value_type> coordView(xView, yView);
      func(coordView);
    });
  }
  else if(values.number_of_children() == 3)
  {
    axom::mir::views::FloatNode_to_ArrayView_same(
      values[0],
      values[1],
      values[2],
      [&](auto xView, auto yView, auto zView) {
        RectilinearCoordsetView3<typename decltype(xView)::value_type> coordView(xView, yView, zView);
        func(coordView);
      });
  }
}

/*!
 * \brief Base template for creating a explicit coordset view.
 */
template <typename DataType, int NDIMS>
struct make_explicit_coordset
{ };

/*!
 * \brief Partial specialization for creating 3D explicit coordset view.
 */
template <typename DataType>
struct make_explicit_coordset<DataType, 3>
{
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<DataType, 3>;

  /*!
   * \brief Create the coordset view and initialize it from the coordset.
   * \param topo The node containing the coordset.
   * \return The coordset view.
   */
  static CoordsetView view(const conduit::Node &coordset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const conduit::Node &values = coordset.fetch_existing("values");
    auto x = bputils::make_array_view<DataType>(values[0]);
    auto y = bputils::make_array_view<DataType>(values[1]);
    auto z = bputils::make_array_view<DataType>(values[2]);
    return CoordsetView(x, y, z);
  }
};

/*!
 * \brief Partial specialization for creating 2D explicit coordset view.
 */
template <typename DataType>
struct make_explicit_coordset<DataType, 2>
{
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<DataType, 2>;

  /*!
   * \brief Create the coordset view and initialize it from the coordset.
   * \param topo The node containing the coordset.
   * \return The coordset view.
   */
  static CoordsetView view(const conduit::Node &coordset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const conduit::Node &values = coordset.fetch_existing("values");
    auto x = bputils::make_array_view<DataType>(values[0]);
    auto y = bputils::make_array_view<DataType>(values[1]);
    return CoordsetView(x, y);
  }
};

/*!
 * \brief Dispatch an explicit coordset to a function.
 *
 * \tparam FuncType The type of the function / lambda to invoke. It is expected
 *                  that the callable accepts an auto argument for a coordset view.
 *
 * \param coordset The Conduit node that contains the coordset.
 * \param func     The function/lambda to invoke using the coordset view.
 */
template <typename FuncType>
void dispatch_explicit_coordset(const conduit::Node &coordset, FuncType &&func)
{
  const conduit::Node &values = coordset["values"];
  if(values.number_of_children() == 2)
  {
    axom::mir::views::FloatNode_to_ArrayView_same(values[0], values[1], [&](auto xView, auto yView) {
      ExplicitCoordsetView<typename decltype(xView)::value_type, 2> coordView(xView, yView);
      func(coordView);
    });
  }
  else if(values.number_of_children() == 3)
  {
    axom::mir::views::FloatNode_to_ArrayView_same(
      values[0],
      values[1],
      values[2],
      [&](auto xView, auto yView, auto zView) {
        ExplicitCoordsetView<typename decltype(xView)::value_type, 3> coordView(xView, yView, zView);
        func(coordView);
      });
  }
}

/*!
 * \brief Given a Conduit/Blueprint coordset, create an appropriate view and
 *        call the supplied function, passing the coordset view to it.
 *
 * \tparam FuncType The type of the function / lambda to invoke. It is expected
 *                  that the callable accepts an auto argument for a coordset view.
 *
 * \param coordset The Conduit node that contains the coordset.
 * \param func     The function/lambda to invoke using the coordset view.
 */
template <typename FuncType>
void dispatch_coordset(const conduit::Node &coordset, FuncType &&func)
{
  const std::string cstype = coordset["type"].as_string();
  if(cstype == "uniform")
  {
    dispatch_uniform_coordset(coordset, func);
  }
  else if(cstype == "rectilinear")
  {
    dispatch_rectilinear_coordset(coordset, func);
  }
  else if(cstype == "explicit")
  {
    // TODO: get the axis names.
    dispatch_explicit_coordset(coordset, func);
  }
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
