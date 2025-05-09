// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_UTILITIES_PRIMAL_ADAPTOR_HPP_
#define AXOM_MIR_UTILITIES_PRIMAL_ADAPTOR_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir/utilities/VariableShape.hpp"
#include "axom/primal.hpp"

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{

/*!
 * \brief This class is a view that exposes shapes in MIR topology and
 *        coordset views as primal shapes. SFINAE is used for the getShape()
 *        method.
 */
template <typename TopologyView, typename CoordsetView>
struct PrimalAdaptor
{
  using value_type = typename CoordsetView::value_type;
  using Polygon =
    axom::primal::Polygon<value_type, CoordsetView::dimension(), axom::primal::PolygonArray::Static>;
  using Tetrahedron = axom::primal::Tetrahedron<value_type, CoordsetView::dimension()>;
  using Hexahedron = axom::primal::Hexahedron<value_type, CoordsetView::dimension()>;
  using BoundingBox = axom::primal::BoundingBox<value_type, CoordsetView::dimension()>;

  /*!
   * \brief Return the number of zones in the associated topology view.
   * \return The number of zones in the associated topology view.
   */
  AXOM_HOST_DEVICE axom::IndexType numberOfZones() const { return m_topologyView.numberOfZones(); }

  /*!
   * \brief Return the bounding box for the zi'th zone.
   *
   * \param zi The index of the zone whose bounding box we want.
   *
   * \return The bounding box of the zi'th zone.
   */
  AXOM_HOST_DEVICE BoundingBox getBoundingBox(axom::IndexType zi) const
  {
    const auto zone = m_topologyView.zone(zi);
    const auto nnodes = zone.numberOfNodes();
    BoundingBox b;
    for(axom::IndexType i = 0; i < nnodes; i++)
    {
      b.addPoint(m_coordsetView[zone.getId(i)]);
    }
    return b;
  }

  // SFINAE methods for returning a zone as a primal shape (or VariableShape
  // if primal lacks the shape type).

  /*!
   * \brief Get the zone \a zi as a polygon. This is enabled when the input topology is 2D.
   *
   * \param zi The index of the zone we want to return.
   *
   * \return A Polygon that represents the zone from the input topology.
   */
  template <int TDIM = CoordsetView::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, Polygon>::type getShape(axom::IndexType zi) const
  {
    const auto zone = m_topologyView.zone(zi);
    Polygon p;
    for(axom::IndexType i = 0; i < zone.numberOfNodes(); i++)
    {
      p.addVertex(m_coordsetView[zone.getId(i)]);
    }
    return p;
  }

  /*!
   * \brief Get the zone \a zi as a Tetrahedron. This is enabled when the input topology contains tets.
   *
   * \param zi The index of the zone we want to return.
   *
   * \return A Tetrahedron that represents the zone from the input topology.
   */
  template <int TDIM = CoordsetView::dimension(), typename ShapeType = typename TopologyView::ShapeType>
  AXOM_HOST_DEVICE typename std::enable_if<
    TDIM == 3 &&
      std::is_same<ShapeType, axom::mir::views::TetShape<typename ShapeType::ConnectivityStorage>>::value,
    Tetrahedron>::type
  getShape(axom::IndexType zi) const
  {
    const auto zone = m_topologyView.zone(zi);
    return Tetrahedron(m_coordsetView[zone.getId(0)],
                       m_coordsetView[zone.getId(1)],
                       m_coordsetView[zone.getId(2)],
                       m_coordsetView[zone.getId(3)]);
  }

  /*!
   * \brief Get the zone \a zi as a Hexahedron. This is enabled when the input topology contains hexs.
   *
   * \param zi The index of the zone we want to return.
   *
   * \return A Hexahedron that represents the zone from the input topology.
   */
  template <int TDIM = CoordsetView::dimension(), typename ShapeType = typename TopologyView::ShapeType>
  AXOM_HOST_DEVICE typename std::enable_if<
    TDIM == 3 &&
      std::is_same<ShapeType, axom::mir::views::HexShape<typename ShapeType::ConnectivityStorage>>::value,
    Hexahedron>::type
  getShape(axom::IndexType zi) const
  {
    const auto zone = m_topologyView.zone(zi);
    return Hexahedron(m_coordsetView[zone.getId(0)],
                      m_coordsetView[zone.getId(1)],
                      m_coordsetView[zone.getId(2)],
                      m_coordsetView[zone.getId(3)],
                      m_coordsetView[zone.getId(4)],
                      m_coordsetView[zone.getId(5)],
                      m_coordsetView[zone.getId(6)],
                      m_coordsetView[zone.getId(7)]);
  }

  /*!
   * \brief Get the zone \a zi as a VariableShape. This is enabled when the input
   *        topology contains wedges, pyramids, or mixed shapes.
   *
   * \param zi The index of the zone we want to return.
   *
   * \return A VariableShape that represents the zone from the input topology.
   *         VariableShape is used because primal lacks Pyramid, Wedge classes and
   *         because if we're using this, we might have a mesh with mixed zone types.
   */
  template <int TDIM = CoordsetView::dimension(), typename ShapeType = typename TopologyView::ShapeType>
  AXOM_HOST_DEVICE typename std::enable_if<
    (TDIM == 3) &&
      (std::is_same<ShapeType, axom::mir::views::PyramidShape<typename ShapeType::ConnectivityType>>::value ||
       std::is_same<ShapeType, axom::mir::views::WedgeShape<typename ShapeType::ConnectivityType>>::value ||
       std::is_same<ShapeType, axom::mir::views::VariableShape<typename ShapeType::ConnectivityType>>::value),
    VariableShape<value_type, 3>>::type
  getShape(axom::IndexType zi) const
  {
    const auto zone = m_topologyView.zone(zi);
    VariableShape<value_type, 3> shape;
    shape.m_shapeId = zone.id();
    for(int i = 0; i < zone.numberOfNodes(); i++)
    {
      shape.push_back(m_coordsetView[zone.getId(i)]);
    }
    return shape;
  }

  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
};

}  // namespace blueprint
}  // namespace utilities
}  // namespace mir
}  // namespace axom

#endif
