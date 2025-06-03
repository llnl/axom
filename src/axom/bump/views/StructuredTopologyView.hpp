// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_
#define AXOM_BUMP_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/views/Shapes.hpp"

#include <type_traits>

namespace axom
{
namespace bump
{
namespace views
{
/*!
 * \brief This class provides a view for Conduit/Blueprint structured grid types.
 *
 * \tparam IndexPolicy The policy for making/using indices.
 */
template <typename IndexPolicy>
class StructuredTopologyView
{
public:
  using IndexingPolicy = IndexPolicy;
  using IndexType = typename IndexingPolicy::IndexType;
  using LogicalIndex = typename IndexingPolicy::LogicalIndex;
  using ConnectivityType = IndexType;
  using Shape1D = LineShape<axom::StackArray<ConnectivityType, 2>>;
  using Shape2D = QuadShape<axom::StackArray<ConnectivityType, 4>>;
  using Shape3D = HexShape<axom::StackArray<ConnectivityType, 8>>;
  using ShapeType = typename std::conditional<
    IndexingPolicy::dimension() == 3,
    Shape3D,
    typename std::conditional<IndexingPolicy::dimension() == 2, Shape2D, Shape1D>::type>::type;

  /*!
   * \brief Return the number of dimensions.
   *
   * \return The number of dimensions.
   */
  AXOM_HOST_DEVICE constexpr static int dimension() { return IndexingPolicy::dimension(); }

  /*!
   * \brief Constructor
   */
  AXOM_HOST_DEVICE StructuredTopologyView() : m_zoneIndexing(), m_nodeIndexing() { }

  /*!
   * \brief Constructor
   *
   * \param indexing The indexing policy for the topology (num zones in each dimension).
   */
  AXOM_HOST_DEVICE StructuredTopologyView(const IndexingPolicy &indexing)
    : m_zoneIndexing(indexing)
    , m_nodeIndexing(indexing.expand())
  { }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE inline IndexType size() const { return m_zoneIndexing.size(); }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE inline IndexType numberOfZones() const { return size(); }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE IndexType connectivitySize() const
  {
    IndexType nodesPerElem = 1;
    for(int d = 0; d < dimension(); d++)
    {
      nodesPerElem *= 2;
    }
    return numberOfZones() * nodesPerElem;
  }

  /*!
   * \brief Return the mesh logical dimensions.
   *
   * \return The mesh logical dimensions.
   */
  AXOM_HOST_DEVICE const LogicalIndex &logicalDimensions() const
  {
    return m_zoneIndexing.logicalDimensions();
  }

  /*!
   * \brief Return indexing object.
   *
   * \return The indexing object.
   */
  AXOM_HOST_DEVICE inline IndexingPolicy &indexing() { return m_zoneIndexing; }

  /*!
   * \brief Return indexing object.
   *
   * \return The indexing object.
   */
  AXOM_HOST_DEVICE inline const IndexingPolicy &indexing() const { return m_zoneIndexing; }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The index of the zone to return.
   *
   * \return The requested zone.
   *
   * \note 3D implementation.
   */
  template <int _ndims = IndexingPolicy::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, Shape3D>::type zone(axom::IndexType zoneIndex) const
  {
    SLIC_ASSERT(zoneIndex < numberOfZones());

    const auto localLogical = m_zoneIndexing.IndexToLogicalIndex(zoneIndex);
    const auto jp = m_nodeIndexing.jStride();
    const auto kp = m_nodeIndexing.kStride();

    Shape3D shape;
    auto &data = shape.getIdsStorage();
    data[0] = m_nodeIndexing.GlobalToGlobal(m_nodeIndexing.LocalToGlobal(localLogical));
    data[1] = data[0] + 1;
    data[2] = data[1] + jp;
    data[3] = data[2] - 1;
    data[4] = data[0] + kp;
    data[5] = data[1] + kp;
    data[6] = data[2] + kp;
    data[7] = data[3] + kp;

    return shape;
  }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The index of the zone to return.
   *
   * \return The requested zone.
   *
   * \note 2D implementation.
   */
  template <int _ndims = IndexingPolicy::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 2, Shape2D>::type zone(axom::IndexType zoneIndex) const
  {
    SLIC_ASSERT(zoneIndex < numberOfZones());

    const auto localLogical = m_zoneIndexing.IndexToLogicalIndex(zoneIndex);
    const auto jp = m_nodeIndexing.jStride();

    Shape2D shape;
    auto &data = shape.getIdsStorage();
    data[0] = m_nodeIndexing.GlobalToGlobal(m_nodeIndexing.LocalToGlobal(localLogical));
    data[1] = data[0] + 1;
    data[2] = data[1] + jp;
    data[3] = data[2] - 1;

    return shape;
  }

  /*!
   * \brief Return a zone.
   *
   * \param index The index of the zone to return.
   *
   * \return The requested zone.
   *
   * \note 1D implementation.
   */
  template <int _ndims = IndexingPolicy::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 1, Shape1D>::type zone(axom::IndexType zoneIndex) const
  {
    SLIC_ASSERT(zoneIndex < numberOfZones());

    const auto localLogical = m_zoneIndexing.IndexToLogicalIndex(zoneIndex);

    Shape1D shape;
    auto &data = shape.getIdsStorage();
    data[0] = m_nodeIndexing.GlobalToGlobal(m_nodeIndexing.LocalToGlobal(localLogical));
    data[1] = data[0] + 1;

    return shape;
  }

private:
  IndexingPolicy m_zoneIndexing;
  IndexingPolicy m_nodeIndexing;
};

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
