// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_UTILITIES_PRIMAL_ADAPTOR_HPP_
#define AXOM_BUMP_UTILITIES_PRIMAL_ADAPTOR_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/bump/VariableShape.hpp"
#include "axom/primal.hpp"

namespace axom
{
namespace bump
{

/*!
 * \brief This struct is a simpler way of representing a polyhedron if we only need to
 *        know about it as set of planes.
 */
template <typename T>
struct PolyhedralFaces
{
  using PlaneType = axom::primal::Plane<T, 3>;
  static constexpr int MAX_PLANES = 64;

  AXOM_HOST_DEVICE inline int size() const { return m_planes.size(); }
  AXOM_HOST_DEVICE inline const PlaneType &operator[](size_t i) const { return m_planes[i]; }
  AXOM_HOST_DEVICE inline PlaneType &operator[](size_t i) { return m_planes[i]; }
  AXOM_HOST_DEVICE inline void push_back(const PlaneType &plane) { m_planes.push_back(plane); }
  AXOM_HOST_DEVICE axom::ArrayView<PlaneType> getFaces() const
  {
    return axom::ArrayView<PlaneType>(const_cast<PlaneType *>(m_planes.data()), m_planes.size());
  }

  axom::StaticArray<PlaneType, MAX_PLANES> m_planes;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const PolyhedralFaces<T> &obj)
{
  os << "PolyhedralFaces\n";
  for(int i = 0; i < obj.size(); i++)
  {
    os << "\t" << i << ": " << obj[i] << std::endl;
  }
  return os;
}

/*!
 * \brief Average a series of points.
 *
 * \tparam T The point precision.
 * \tparam NDIMS The number of dimensions in the point.
 */
template <typename T, int NDIMS>
struct AveragePoints
{
  using PointType = axom::primal::Point<T, NDIMS>;
  using VectorType = axom::primal::Vector<T, NDIMS>;

  /*!
   * \brief Add a point to the running sum.
   *
   * \param pt The point to add.
   */
  AXOM_HOST_DEVICE inline void add(const PointType &pt)
  {
    sum += VectorType(pt);
    numPoints++;
  }

  /*!
   * \brief Clear the sum.
   */
  AXOM_HOST_DEVICE inline void clear()
  {
    numPoints = 0;
    sum = VectorType {};
  }

  /*!
   * \brief Get the average point value.
   *
   * \return The average point value.
   */
  AXOM_HOST_DEVICE inline PointType get() const
  {
    return (numPoints > 0)
      ? PointType(sum.array() * (static_cast<T>(1.) / static_cast<T>(numPoints)))
      : PointType {};
  }

  VectorType sum {};
  int numPoints {0};
};

/*!
 * \brief Make a primal::Polyhedron from a polyhedron shape.
 */
template <typename TopologyView, typename CoordsetView, bool makeFaces>
struct AdaptPolyhedron
{
  using value_type = typename CoordsetView::value_type;
  using Polyhedron = axom::primal::Polyhedron<value_type, 3>;
  using PolyhedralRepresentation = Polyhedron;

  /*!
   * \brief Return a zone as primal::Polyhedron.
   *
   * \param topologyView The topology view.
   * \param coordsetView The coordset view.
   * \param zoneIndex The index of the zone to convert.
   *
   * \note This is a lot of work to convert Blueprint to primal::Polyhedron. It
   *       probably needs checking to ensure proper ordering of the vertex neighbors,
   *       which if you look at the point, need to be added in a counter-clockwise order.
   *
   * \return A representation of the polyhedral zone.
   */
  AXOM_HOST_DEVICE static PolyhedralRepresentation convert(const TopologyView &topologyView,
                                                           const CoordsetView &coordsetView,
                                                           size_t zoneIndex)
  {
    const auto zone = topologyView.zone(zoneIndex);
    const auto uniqueNodeIds = zone.getUniqueIds();
    Polyhedron poly;

    // Sort the ids.
    const axom::IndexType nnodes = uniqueNodeIds.size();
    StaticArray<int, Polyhedron::MAX_VERTS> ids;
    SLIC_ASSERT(nnodes < Polyhedron::MAX_VERTS);
    for(axom::IndexType i = 0; i < nnodes; i++)
    {
      ids.push_back(uniqueNodeIds[i]);
    }
    axom::utilities::Sorting<int, Polyhedron::MAX_VERTS>::sort(ids.data(), nnodes);

    // Add vertices in ids sorted order.
    for(axom::IndexType i = 0; i < nnodes; i++)
    {
      poly.addVertex(coordsetView[ids[i]]);
    }

    // Get the center point of the vertices.
    auto center = poly.vertexMean();

    // Add neighbors for each vertex.
    for(axom::IndexType i = 0; i < nnodes; i++)
    {
      const auto currentNodeId = ids[i];

      StaticArray<int, Polyhedron::MAX_VERTS> seenNeighbors;

      // Emit neigbor nodes for the current node. We have to scan through the
      // faces though.
      for(axom::IndexType f = 0; f < zone.numberOfFaces(); f++)
      {
        const auto faceIds = zone.getFace(f);
        const auto lastIndex = faceIds.size() - 1;

        // Check which way the plane normal points relative to the zone center.
        // We want to point towards it.
        const auto p0 = coordsetView[faceIds[0]];
        const auto p1 = coordsetView[faceIds[1]];
        const auto p2 = coordsetView[faceIds[2]];
        auto plane = axom::primal::make_plane(p0, p1, p2);
        const bool reverseOrder = (plane.signedDistance(center) > 0.);

        for(int fi = 0; fi < faceIds.size(); fi++)
        {
          if(faceIds[fi] == currentNodeId)
          {
            // Neighbors for currentNodeId
            int candidates[2];
            candidates[0] = static_cast<int>(faceIds[(fi == 0) ? lastIndex : (fi - 1)]);
            candidates[1] = static_cast<int>(faceIds[(fi == lastIndex) ? 0 : (fi + 1)]);
            if(reverseOrder)
            {
              axom::utilities::swap(candidates[0], candidates[1]);
            }

            for(int ci = 0; ci < 2; ci++)
            {
              // Check whether this neighbor has been seen before.
              bool found = false;
              for(axom::IndexType ni = 0; ni < seenNeighbors.size() && !found; ni++)
              {
                found |= (seenNeighbors[ni] == candidates[ci]);
              }

              // If the neighbor has not been seen before, record it, and add
              // it to the polyhedron.
              if(!found)
              {
                seenNeighbors.push_back(candidates[ci]);

                // Look up the index of the candidate point in the sorted indices.
                auto neighborIndex = axom::utilities::binary_search(ids, candidates[ci]);
                SLIC_ASSERT(neighborIndex != -1);

                poly.addNeighbors(i, neighborIndex);
              }
            }

            break;
          }
        }
      }
    }
    return poly;
  }
};

/*!
 * \brief Make PolyhedralFaces from a polyhedron shape.
 *
 * \note This is a partial template specialization for makeFaces=true that makes
 *       a PolyhedralFaces object instead of a primal::Polyhedron.
 */
template <typename TopologyView, typename CoordsetView>
struct AdaptPolyhedron<TopologyView, CoordsetView, true>
{
  using value_type = typename CoordsetView::value_type;
  using PolyhedralRepresentation = PolyhedralFaces<value_type>;

  /*!
   * \brief Return a zone as PolyhedralFaces.
   *
   * \param topologyView The topology view.
   * \param coordsetView The coordset view.
   * \param zoneIndex The index of the zone to convert.
   *
   * \return A representation of the polyhedral zone.
   */
  AXOM_HOST_DEVICE static PolyhedralRepresentation convert(const TopologyView &topologyView,
                                                           const CoordsetView &coordsetView,
                                                           size_t zoneIndex)
  {
    PolyhedralRepresentation faces;
    const auto zone = topologyView.zone(zoneIndex);
    const int numFaces = static_cast<int>(zone.numberOfFaces());
    const value_type eps = std::is_same<value_type, float>::value ? 1.e-6 : 1.e-10;

    // Make a zone center.
    AveragePoints<value_type, 3> avg;
    const auto nodeIds = zone.getIds();
    for(const auto id : nodeIds)
    {
      avg.add(coordsetView[id]);
    }
    const auto center = avg.get();

    // Go through the faces and add them to the representation.
    for(int f = 0; f < numFaces; f++)
    {
      const auto faceIds = zone.getFace(f);

      // Make the first plane in the face.
      const auto p0 = coordsetView[faceIds[0]];
      const auto p1 = coordsetView[faceIds[1]];
      const auto p2 = coordsetView[faceIds[2]];

      // Check which way the plane normal points relative to the zone center.
      // We want to point towards it.
      auto plane = axom::primal::make_plane(p0, p1, p2);
      if(plane.signedDistance(center) < 0.)
      {
        plane.flip();
      }

      faces.push_back(plane);

      // We've already made 1 plane for the face. Treat the face as a triangle
      // fan and add the rest of the planes, if they are different enough from
      // the previous plane.
      const int extraPlanes = faceIds.size() - 3;
      for(int ep = 0; ep < extraPlanes; ep++)
      {
        // Get the last point in the triangle.
        const auto ep2 = coordsetView[faceIds[3 + ep]];

        // Check its signed distance vs the previous plane. If the point is far
        // enough away, we make the plane.
        const auto dist = plane.signedDistance(ep2);
        if(axom::utilities::abs(dist) > eps)
        {
          const auto ep1 = coordsetView[faceIds[2 + ep]];
          auto plane2 = axom::primal::make_plane(p0, ep1, ep2);
          if(plane2.signedDistance(center) < 0.)
          {
            plane2.flip();
          }
          faces.push_back(plane2);
          plane = plane2;
        }
      }
    }

    return faces;
  }
};

/*!
 * \brief This class is a view that exposes shapes in MIR topology and
 *        coordset views as primal shapes. SFINAE is used for the getShape()
 *        method.
 *
 * \tparam TopologyView The topology view type.
 * \tparam CoordsetView The coordset view type.
 * \tparam MAX_VERTS_2D The maximum number of vertices allowed in a polygon.
 * \tparam makeFaces Whether to make faces for polyhedral shapes or to make primal::Polyhedron.
 */
template <typename TopologyView, typename CoordsetView, int MAX_VERTS_2D = 12, bool makeFaces = false>
struct PrimalAdaptor
{
  using value_type = typename CoordsetView::value_type;
  using Polygon =
    axom::primal::Polygon<value_type, CoordsetView::dimension(), axom::primal::PolygonArray::Static, MAX_VERTS_2D>;
  using Tetrahedron = axom::primal::Tetrahedron<value_type, CoordsetView::dimension()>;
  using Hexahedron = axom::primal::Hexahedron<value_type, CoordsetView::dimension()>;
  using Polyhedron =
    typename AdaptPolyhedron<TopologyView, CoordsetView, makeFaces>::PolyhedralRepresentation;
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
    const auto nodeIds = zone.getIds();
    const auto nnodes = nodeIds.size();
    BoundingBox b;
    for(axom::IndexType i = 0; i < nnodes; i++)
    {
      b.addPoint(m_coordsetView[nodeIds[i]]);
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
      std::is_same<ShapeType, axom::bump::views::TetShape<typename ShapeType::ConnectivityStorage>>::value,
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
      std::is_same<ShapeType, axom::bump::views::HexShape<typename ShapeType::ConnectivityStorage>>::value,
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
      (std::is_same<ShapeType, axom::bump::views::PyramidShape<typename ShapeType::ConnectivityType>>::value ||
       std::is_same<ShapeType, axom::bump::views::WedgeShape<typename ShapeType::ConnectivityType>>::value ||
       std::is_same<ShapeType, axom::bump::views::VariableShape<typename ShapeType::ConnectivityType>>::value),
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

  /*!
   * \brief Get the zone \a zi as a Polyhedron. This is enabled when the input topology contains polyhedra.
   *
   * \param zi The index of the zone we want to return.
   *
   * \return A Polyhedron that represents the zone from the input topology.
   */
  template <int TDIM = CoordsetView::dimension(), typename ShapeType = typename TopologyView::ShapeType>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3 && ShapeType::is_polyhedral(), Polyhedron>::type
  getShape(axom::IndexType zi) const
  {
    // Delegate out to the AdaptPolyhedron classes to make the polyhedron.
    return AdaptPolyhedron<TopologyView, CoordsetView, makeFaces>::convert(m_topologyView,
                                                                           m_coordsetView,
                                                                           zi);
  }

  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
};

}  // namespace bump
}  // namespace axom

#endif
