// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_
#define AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mir/views/BasicIndexing.hpp"
#include "axom/mir/views/Shapes.hpp"

namespace axom
{
namespace mir
{
namespace views
{
namespace detail
{

// Base template
template <typename ShapeType, bool VariableSize>
struct NumberOfZones
{ };

// Specialization for fixed-size shapes (user might not provide sizes, offsets)
template <typename ShapeType>
struct NumberOfZones<ShapeType, false>
{
  /*!
   * \brief Compute number of zones using either sizes or connSize and knowledge
   *        of a shape with a static number of nodes per zone.
   */
  AXOM_HOST_DEVICE static axom::IndexType execute(axom::IndexType sizesSize, axom::IndexType connSize)
  {
    return (sizesSize != 0) ? sizesSize : (connSize / ShapeType::numberOfNodes());
  }
};

// Specialization for variable-size shapes.
template <typename ShapeType>
struct NumberOfZones<ShapeType, true>
{
  /*!
   * \brief Compute number of zones using sizes. A variable-sized shape must have provided sizes.
   */
  AXOM_HOST_DEVICE static axom::IndexType execute(axom::IndexType sizesSize,
                                                  axom::IndexType AXOM_UNUSED_PARAM(connSize))
  {
    return sizesSize;
  }
};

}  // end namespace detail

/*!
 * \brief This class provides a view for Conduit/Blueprint single shape unstructured grids.
 *
 * \tparam IndexT The index type that will be used for connectivity, etc.
 * \tparam ShapeT The shape type.
 */
template <typename ShapeT>
class UnstructuredTopologySingleShapeView
{
public:
  using ShapeType = ShapeT;
  using ConnectivityType = typename ShapeType::ConnectivityType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  using IndexingPolicy = BasicIndexing;

  /*!
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   */
  AXOM_HOST_DEVICE
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn)
    : m_connectivityView(conn)
    , m_sizesView()
    , m_offsetsView()
    , m_indexing()
  {
    SLIC_ASSERT(!ShapeType::is_variable_size());
    m_indexing.m_size = numberOfZones();
  }

  /*!
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   * \param sizes The number of nodes in each zone.
   * \param offsets The offset to each zone in the connectivity.
   */
  AXOM_HOST_DEVICE
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn,
                                      const ConnectivityView &sizes,
                                      const ConnectivityView &offsets)
    : m_connectivityView(conn)
    , m_sizesView(sizes)
    , m_offsetsView(offsets)
    , m_indexing(sizes.size())
  {
    SLIC_ASSERT(m_offsetsView.size() == m_sizesView.size());
  }

  /*!
   * \brief Return the dimension of the shape.
   *
   * \return The dimension of the shape.
   */
  AXOM_HOST_DEVICE static constexpr int dimension() { return ShapeType::dimension(); }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE inline IndexType numberOfZones() const
  {
    return detail::NumberOfZones<ShapeType, ShapeType::is_variable_size()>::execute(
      m_sizesView.size(),
      m_connectivityView.size());
  }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE inline IndexType connectivitySize() const { return m_connectivityView.size(); }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE inline const IndexingPolicy &indexing() const { return m_indexing; }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The requested zone.
   *
   * \return The requested zone.
   */
  /// @{
  template <bool _variable_size = ShapeType::is_variable_size()>
  AXOM_HOST_DEVICE typename std::enable_if<_variable_size, ShapeType>::type zone(
    axom::IndexType zoneIndex) const
  {
    SLIC_ASSERT(zoneIndex < numberOfZones());

    return ShapeType(ConnectivityView(m_connectivityView.data() + m_offsetsView[zoneIndex],
                                      m_sizesView[zoneIndex]));
  }

  template <bool _variable_size = ShapeType::is_variable_size()>
  AXOM_HOST_DEVICE typename std::enable_if<!_variable_size, ShapeType>::type zone(
    axom::IndexType zoneIndex) const
  {
    SLIC_ASSERT(zoneIndex < numberOfZones());

    ConnectivityView shapeIdsView {};
    if(m_sizesView.empty())
    {
      shapeIdsView = ConnectivityView(m_connectivityView.data() + ShapeType::zoneOffset(zoneIndex),
                                      ShapeType::numberOfNodes());
    }
    else
    {
      shapeIdsView = ConnectivityView(m_connectivityView.data() + m_offsetsView[zoneIndex],
                                      m_sizesView[zoneIndex]);
    }
    return ShapeType(shapeIdsView);
  }
  /// @}

private:
  ConnectivityView m_connectivityView;
  ConnectivityView m_sizesView;
  ConnectivityView m_offsetsView;
  IndexingPolicy m_indexing;
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
