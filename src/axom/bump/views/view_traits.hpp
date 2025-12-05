// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_VIEW_TRAITS_HPP_
#define AXOM_BUMP_VIEW_TRAITS_HPP_

#include "axom/core/utilities/BitUtilities.hpp"
#include "axom/bump/views/StructuredTopologyView.hpp"
#include "axom/bump/views/StructuredIndexing.hpp"
#include "axom/bump/views/StridedStructuredIndexing.hpp"
#include "axom/bump/views/Shapes.hpp"

namespace axom
{
namespace bump
{
namespace views
{
/// Return which shapes to expect for a given dimension.
static constexpr int shapes_for_dimension(int dimension)
{
  int shapes = 0;
  switch(dimension)
  {
  case 0:
    axom::utilities::setBitOn(shapes, Point_ShapeID);
    break;
  case 1:
    axom::utilities::setBitOn(shapes, Line_ShapeID);
    break;
  case 2:
    axom::utilities::setBitOn(shapes, Tri_ShapeID);
    axom::utilities::setBitOn(shapes, Quad_ShapeID);
    axom::utilities::setBitOn(shapes, Polygon_ShapeID);
    axom::utilities::setBitOn(shapes, Mixed_ShapeID);
    break;
  case 3:
    axom::utilities::setBitOn(shapes, Tet_ShapeID);
    axom::utilities::setBitOn(shapes, Pyramid_ShapeID);
    axom::utilities::setBitOn(shapes, Wedge_ShapeID);
    axom::utilities::setBitOn(shapes, Hex_ShapeID);
    axom::utilities::setBitOn(shapes, Mixed_ShapeID);
    break;
  }
  return shapes;
}

/// General traits for topology views.
template <typename TopologyView>
struct view_traits
{
  static constexpr bool supports_strided_structured() { return false; }
  static constexpr int selected_shapes() { return shapes_for_dimension(TopologyView::dimension()); }
};

/// If StructuredTopologyView was instantiated with StridedStructuredIndexing
/// (of varying dimensions) then say that strided structured is supported.
template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 3>>>
{
  static constexpr bool supports_strided_structured() { return true; }
  static constexpr int selected_shapes() { return shapes_for_dimension(3); }
};

template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 2>>>
{
  static constexpr bool supports_strided_structured() { return true; }
  static constexpr int selected_shapes() { return shapes_for_dimension(2); }
};

template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 1>>>
{
  static constexpr bool supports_strided_structured() { return true; }
  static constexpr int selected_shapes() { return shapes_for_dimension(1); }
};

/*!
 * \brief Base template for some ArrayView traits.
 */
template <typename T>
struct array_view_traits
{ };

/*!
 * \brief This macro defines some template specializations that help us access
 *        ArrayView<> names. This can be helpful when the ArrayView comes into
 *        a lambda as an auto argument.
 */
#define AXOM_MAKE_TRAIT(TYPE)                                                          \
  template <>                                                                          \
  struct array_view_traits<axom::ArrayView<TYPE>>                                      \
  {                                                                                    \
    using value_type = TYPE;                                                           \
    static constexpr const char* name() { return "axom::ArrayView<" #TYPE ">"; }       \
  };                                                                                   \
  template <>                                                                          \
  struct array_view_traits<axom::ArrayView<const TYPE>>                                \
  {                                                                                    \
    using value_type = const TYPE;                                                     \
    static constexpr const char* name() { return "axom::ArrayView<const " #TYPE ">"; } \
  };

AXOM_MAKE_TRAIT(std::uint8_t)
AXOM_MAKE_TRAIT(std::uint16_t)
AXOM_MAKE_TRAIT(std::uint32_t)
AXOM_MAKE_TRAIT(std::uint64_t)
AXOM_MAKE_TRAIT(std::int8_t)
AXOM_MAKE_TRAIT(std::int16_t)
AXOM_MAKE_TRAIT(std::int32_t)
AXOM_MAKE_TRAIT(std::int64_t)
AXOM_MAKE_TRAIT(float)
AXOM_MAKE_TRAIT(double)
#undef AXOM_MAKE_TRAIT

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
