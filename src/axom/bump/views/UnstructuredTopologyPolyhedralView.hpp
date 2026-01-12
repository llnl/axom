// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_VIEWS_UNSTRUCTURED_TOPOLOGY_POLYHEDRAL_VIEW_HPP_
#define AXOM_BUMP_VIEWS_UNSTRUCTURED_TOPOLOGY_POLYHEDRAL_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/views/BasicIndexing.hpp"
#include "axom/bump/views/Shapes.hpp"

namespace axom
{
namespace bump
{
namespace views
{

/*!
 * \brief This class implements a view for Blueprint polyhedral topologies.
 */
template <typename ConnType>
class UnstructuredTopologyPolyhedralView
{
public:
  using ConnectivityType = ConnType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  using IndexingPolicy = BasicIndexing;

  /*!
   * \brief This struct contains views that hold polyhedral connectivity.
   */
  struct PolyhedronData
  {
    /// Constructor
    AXOM_HOST_DEVICE
    PolyhedronData(const ConnectivityView &subelement_conn,
                   const ConnectivityView &subelement_sizes,
                   const ConnectivityView &subelement_offsets,
                   const ConnectivityView &element_conn,
                   const ConnectivityView &element_sizes,
                   const ConnectivityView &element_offsets)
      : m_subelement_conn(subelement_conn)
      , m_subelement_sizes(subelement_sizes)
      , m_subelement_offsets(subelement_offsets)
      , m_element_conn(element_conn)
      , m_element_sizes(element_sizes)
      , m_element_offsets(element_offsets)
      , m_indexing(element_sizes.size())
    {
#if !defined(AXOM_DEVICE_CODE)
      SLIC_ERROR_IF(subelement_conn.empty() && !subelement_sizes.empty(),
                    "Array views for subelements are inconsistent.");
      SLIC_ERROR_IF(subelement_sizes.size() != subelement_offsets.size(),
                    "Array views for subelement sizes,offsets have different sizes.");
      SLIC_ERROR_IF(element_conn.empty() && !element_sizes.empty(),
                    "Array views for elements are inconsistent.");
      SLIC_ERROR_IF(element_sizes.size() != element_offsets.size(),
                    "Array views for element sizes,offsets have different sizes.");
#endif
    }

    /// Copy Constructor
    AXOM_HOST_DEVICE
    PolyhedronData(const PolyhedronData &obj)
      : m_subelement_conn(obj.m_subelement_conn)
      , m_subelement_sizes(obj.m_subelement_sizes)
      , m_subelement_offsets(obj.m_subelement_offsets)
      , m_element_conn(obj.m_element_conn)
      , m_element_sizes(obj.m_element_sizes)
      , m_element_offsets(obj.m_element_offsets)
      , m_indexing(obj.m_indexing)
    { }

    ConnectivityView m_subelement_conn;
    ConnectivityView m_subelement_sizes;
    ConnectivityView m_subelement_offsets;
    ConnectivityView m_element_conn;
    ConnectivityView m_element_sizes;
    ConnectivityView m_element_offsets;
    IndexingPolicy m_indexing;
  };

  /*!
   * \brief This struct provides data about Zone i's shape.
   */
  struct PolyhedronShape : public PolyhedronTraits
  {
    constexpr static IndexType MaximumNumberOfIds = 20 * 3;

    /// Constructor.
    AXOM_HOST_DEVICE PolyhedronShape(const PolyhedronData &obj, axom::IndexType zi)
      : m_data(obj)
      , m_zoneIndex(zi)
      , m_ids()
    { }

    /// This implementation does not return unique number of nodes.
    AXOM_HOST_DEVICE IndexType numberOfNodes() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
      {
        nnodes += getFace(f).size();
      }
      return nnodes;
    }

    AXOM_HOST_DEVICE IndexType numberOfFaces() const { return m_data.m_element_sizes[m_zoneIndex]; }

    AXOM_HOST_DEVICE IndexType numberOfNodesInFace(int faceIndex) const
    {
      return getFace(faceIndex).size();
    }

    AXOM_HOST_DEVICE ConnectivityView getIds() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
      {
        const auto faceIds = getFace(f);
        for(axom::IndexType i = 0; i < faceIds.size(); i++)
        {
          if(nnodes + 1 <= MaximumNumberOfIds)
            m_ids[nnodes++] = faceIds[i];
          else
          {
#if !defined(AXOM_DEVICE_CODE)
            SLIC_ERROR("m_ids is not large enough to hold all node ids.");
#endif
            break;
          }
        }
      }
      return ConnectivityView(m_ids.m_data, nnodes);
    }

    AXOM_HOST_DEVICE ConnectivityType getId(size_t index) const
    {
      const auto allIds = getIds();
      return allIds[index];
    }

    AXOM_HOST_DEVICE ConnectivityView getUniqueIds() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
      {
        const auto faceIds = getFace(f);
        for(axom::IndexType i = 0; i < faceIds.size(); i++)
        {
          if(!find(m_ids.m_data, nnodes, faceIds[i]))
          {
            if(nnodes + 1 <= MaximumNumberOfIds)
              m_ids[nnodes++] = faceIds[i];
            else
            {
#if !defined(AXOM_DEVICE_CODE)
              SLIC_ERROR("m_ids is not large enough to hold all node ids.");
#endif
              break;
            }
          }
        }
      }
      return ConnectivityView(m_ids.m_data, nnodes);
    }

    AXOM_HOST_DEVICE ConnectivityView getFace(int faceIndex) const
    {
      const ConnectivityView element_face_ids(
        m_data.m_element_conn.data() + m_data.m_element_offsets[m_zoneIndex],
        m_data.m_element_sizes[m_zoneIndex]);
      const auto faceId = element_face_ids[faceIndex];

      return ConnectivityView(m_data.m_subelement_conn.data() + m_data.m_subelement_offsets[faceId],
                              m_data.m_subelement_sizes[faceId]);
    }

    AXOM_HOST_DEVICE void getFace(int faceIndex, ConnectivityType *ids, axom::IndexType &numIds) const
    {
      const auto faceIds = getFace(faceIndex);
      numIds = faceIds.size();
      for(axom::IndexType i = 0; i < numIds; i++)
      {
        ids[i] = faceIds[i];
      }
    }

  private:
    AXOM_HOST_DEVICE bool find(const ConnectivityType *arr,
                               axom::IndexType n,
                               ConnectivityType value) const
    {
      bool found = false;
      for(axom::IndexType i = 0; i < n && !found; i++)
      {
        found = arr[i] == value;
      }
      return found;
    }

    PolyhedronData m_data;
    IndexType m_zoneIndex {0};
    mutable axom::StackArray<ConnectivityType, MaximumNumberOfIds> m_ids;
  };
  //----------------------------------------------------------------------------

  using ShapeType = PolyhedronShape;

  /*!
   * \brief Constructor.
   */
  AXOM_HOST_DEVICE
  UnstructuredTopologyPolyhedralView(const ConnectivityView &subelement_conn,
                                     const ConnectivityView &subelement_sizes,
                                     const ConnectivityView &subelement_offsets,
                                     const ConnectivityView &element_conn,
                                     const ConnectivityView &element_sizes,
                                     const ConnectivityView &element_offsets)
    : m_data(subelement_conn,
             subelement_sizes,
             subelement_offsets,
             element_conn,
             element_sizes,
             element_offsets)
  { }

  /*!
   * \brief Return the number of zones in the mesh.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE
  IndexType numberOfZones() const { return m_data.m_element_sizes.size(); }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE
  IndexType connectivitySize() const { return m_data.element_conn.size(); }

  /*!
   * \brief Return the dimension of the shape.
   *
   * \return The dimension of the shape.
   */
  AXOM_HOST_DEVICE static constexpr int dimension() { return 3; }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE inline const IndexingPolicy &indexing() const { return m_data.m_indexing; }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The requested zone.
   *
   * \return The requested zone.
   */
  AXOM_HOST_DEVICE ShapeType zone(axom::IndexType zoneIndex) const
  {
    SLIC_ASSERT(zoneIndex < numberOfZones());

    return ShapeType(m_data, zoneIndex);
  }

private:
  PolyhedronData m_data;
};

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
