// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_VIEWS_SHAPES_HPP_
#define AXOM_BUMP_VIEWS_SHAPES_HPP_

#include "axom/core/ArrayView.hpp"
#include "axom/slic.hpp"

#include <cassert>
#include <iostream>
#include <type_traits>

namespace axom
{
namespace bump
{
namespace views
{
/*!
 * \brief Shape ids. These are used to identify shapes. These are used as
 *        indices in bit fields in some algorithms.
 */
enum
{
  Point_ShapeID = 0,
  Line_ShapeID = 1,
  Tri_ShapeID = 2,
  Quad_ShapeID = 3,
  Polygon_ShapeID = 4,
  Tet_ShapeID = 5,
  Pyramid_ShapeID = 6,
  Wedge_ShapeID = 7,
  Hex_ShapeID = 8,
  Polyhedron_ShapeID = 9,
  Mixed_ShapeID = 10,

  Invalid_ShapeID = 20
};

/// prototype
AXOM_HOST_DEVICE constexpr IndexType shapeDimension(int shapeId);

/*!
 \brief Point type traits.

\verbatim

  0*

\endverbatim
 */
struct PointTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Point_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 0; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int AXOM_UNUSED_PARAM(faceIndex))
  {
    return 1;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 1; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 0; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 0; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return zoneIndex; }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 1> getFace(IndexType faceIndex)
  {
#if !defined(AXOM_DEBUG)
    AXOM_UNUSED_VAR(faceIndex);
#endif
    assert(faceIndex == 0);
    return StackArray<IndexType, 1> {0};
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(
    int AXOM_UNUSED_PARAM(edgeIndex))
  {
    return axom::StackArray<IndexType, 2> {0, 0};
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "point"; }
};

/*!
 \brief Line type traits.

\verbatim

  0*-----------* 1

\endverbatim
 */
struct LineTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Line_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 1; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 2; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int AXOM_UNUSED_PARAM(faceIndex))
  {
    return 0;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 0; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex)
  {
    return numberOfNodes() * zoneIndex;
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getFace(IndexType faceIndex)
  {
#if !defined(AXOM_DEBUG)
    AXOM_UNUSED_VAR(faceIndex);
#endif
    assert(faceIndex == 0);
    return StackArray<IndexType, 2> {0, 1};
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int /*edgeIndex*/)
  {
    return axom::StackArray<IndexType, 2> {0, 1};
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "line"; }
};

/*!
 \brief Triangle type traits.

\verbatim

  2*
   |\
   | \
   |  \
   |   \
   |    \
  0*-----* 1

\endverbatim
 */
struct TriTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Tri_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int AXOM_UNUSED_PARAM(faceIndex))
  {
    return 3;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex)
  {
    return numberOfNodes() * zoneIndex;
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 3> getFace(IndexType faceIndex)
  {
#if !defined(AXOM_DEBUG)
    AXOM_UNUSED_VAR(faceIndex);
#endif
    assert(faceIndex == 0);
    return StackArray<IndexType, 3> {0, 1, 2};
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0, 1}, {1, 2}, {2, 0}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "tri"; }
};

/*!
 \brief Quad type traits.

\verbatim

  3*-----------* 2
   |           |
   |           |
   |           |
   |           |
   |           |
  0*-----------* 1

\endverbatim
 */
struct QuadTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Quad_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int AXOM_UNUSED_PARAM(faceIndex))
  {
    return 4;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex)
  {
    return numberOfNodes() * zoneIndex;
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 4> getFace(IndexType faceIndex)
  {
#if !defined(AXOM_DEBUG)
    AXOM_UNUSED_VAR(faceIndex);
#endif
    assert(faceIndex == 0);
    return StackArray<IndexType, 4> {0, 1, 2, 3};
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "quad"; }
};

/*!
 \brief Tet type traits.

\verbatim

      3
      *
     /|\         face 0: 0,2,1
    / | \        face 1: 0,1,3
   /  |  \       face 2: 1,2,3
  /   |   \      face 3: 2,0,3
0*----|----* 2
  \   |   /      edge 0: 0,1
   \  |  /       edge 1: 1,2
    \ | /        edge 2: 2,0
     \|/         edge 3: 0,3
      *          edge 4: 1,3
       1         edge 5: 2,3

\endverbatim
 */
struct TetTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Tet_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int AXOM_UNUSED_PARAM(faceIndex))
  {
    return 3;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex)
  {
    return numberOfNodes() * zoneIndex;
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 3> getFace(IndexType faceIndex)
  {
    const axom::StackArray<IndexType, 3> faces[] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};
    assert(faceIndex >= 0 && faceIndex < numberOfFaces());
    return faces[faceIndex];
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "tet"; }
};

/*!
 \brief Pyramid type traits.

\verbatim

3*-----------* 2  face 0: 3,2,1,0
 |\         /|    face 1: 0,1,4
 | \       / |    face 2: 1,2,4
 |  \     /  |    face 3: 2,3,4
 |   \   /   |    face 4: 3,0,4
 |    \ /    |
 |     * 4   |    edge 0: 0,1
 |    / \    |    edge 1: 1,2
 |   /   \   |    edge 2: 2,3
 |  /     \  |    edge 3: 3,0
 | /       \ |    edge 4: 0,4
 |/         \|    edge 5: 1,4
0*-----------* 1  edge 6: 2,4
                  edge 7: 3,4

\endverbatim
 */
struct PyramidTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Pyramid_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int faceIndex)
  {
    return faceIndex == 0 ? 4 : 3;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex)
  {
    return numberOfNodes() * zoneIndex;
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 4> getFace(IndexType faceIndex)
  {
    const axom::StackArray<IndexType, 4> faces[] = {{3, 2, 1, 0},
                                                    {0, 1, 4, -1},
                                                    {1, 2, 4, -1},
                                                    {2, 3, 4, -1},
                                                    {3, 0, 4, -1}};
    assert(faceIndex >= 0 && faceIndex < numberOfFaces());
    return faces[faceIndex];
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] =
      {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "pyramid"; }
};

/*!
 \brief Wedge type traits.

\verbatim

3*---------* 5  face 0: 0,2,1
 |\       /|    face 1: 3,4,5
 | \     / |    face 2: 0,1,4,3
 |  \   /  |    face 3: 1,2,5,4
 |   \ /   |    face 4: 2,0,3,5
 |    *4   |
 |    |    |    edge 0: 0,1
0*----|----* 2  edge 1: 1,2
  \   |   /     edge 2: 2,0
   \  |  /      edge 3: 3,4
    \ | /       edge 4: 4,5
     \|/        edge 5: 5,3
      * 1       edge 6: 0,3
                edge 7: 1,4
                edge 8: 2,3

\endverbatim
 */
struct WedgeTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Wedge_ShapeID; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int faceIndex)
  {
    return (faceIndex == 0 || faceIndex == 1) ? 3 : 4;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 9; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex)
  {
    return numberOfNodes() * zoneIndex;
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 4> getFace(IndexType faceIndex)
  {
    const axom::StackArray<IndexType, 4> faces[] = {{0, 2, 1, -1},
                                                    {3, 4, 5, -1},
                                                    {0, 1, 4, 3},
                                                    {1, 2, 5, 4},
                                                    {2, 0, 3, 5}};
    assert(faceIndex >= 0 && faceIndex < numberOfFaces());
    return faces[faceIndex];
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] =
      {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "wedge"; }
};

/*!
 \brief Hex type traits.

\verbatim

      4*------------* 7
      /|           /|     
     / |          / |
    /  |         /  |
  5*------------*6  |
   |   |        |   |
   |   |        |   |
   |  0*--------|---* 3
   |  /         |  /
   | /          | /
   |/           |/
   *------------*
   1            2

\endverbatim
 */
struct HexTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Hex_ShapeID; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int AXOM_UNUSED_PARAM(faceIndex))
  {
    return 4;
  }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 12; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex)
  {
    return numberOfNodes() * zoneIndex;
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 4> getFace(IndexType faceIndex)
  {
    const axom::StackArray<IndexType, 4> faces[] =
      {{3, 0, 4, 7}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};
    assert(faceIndex >= 0 && faceIndex < numberOfFaces());
    return faces[faceIndex];
  }

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] =
      {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "hex"; }
};

/*!

\verbatim

n-1 *-... * 2
    |     |
    |     |
  0 *-----* 1

\endverbatim
 */
struct PolygonTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Polygon_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return true; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 20; }
  AXOM_HOST_DEVICE constexpr static const char *name() { return "polygonal"; }
};

/*!
 * \brief Polyhedron shape traits.
 *
 * \note This struct contains a subset of the interface used in some other shape
 *       traits.
 */
struct PolyhedronTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Polyhedron_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return true; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return true; }
  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }
  AXOM_HOST_DEVICE constexpr static const char *name() { return "polyhedral"; }
};

/*!
 * \brief This struct represents a polygon zone.
 */
template <typename ConnType>
struct PolygonShape : public PolygonTraits
{
  using ConnectivityType = ConnType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  using ConnectivityStorage = ConnectivityType;
  using ConnectivityStorageRef = ConnectivityType &;
  using ConnectivityStorageConstRef = const ConnectivityType &;

  /*!
   * \brief Construct a shape.
   */
  AXOM_HOST_DEVICE PolygonShape(const ConnectivityView &ids) : m_ids(ids) { }

  /*!
   * \brief Return the number of nodes in the polygon.
   *
   * \return The number of nodes in the polygon.
   */
  AXOM_HOST_DEVICE IndexType numberOfNodes() const { return m_ids.size(); }

  /*!
   * \brief Get a specific id that makes up this shape.
   *
   * \return The i'th id that makes up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityType getId(size_t index) const
  {
    SLIC_ASSERT(index < static_cast<size_t>(m_ids.size()));
    return m_ids[index];
  }

  /*!
   * \brief Get the storage for the ids that make up this shape.
   *
   * \return The container for the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityStorageRef getIdsStorage() { return m_ids; }

  /*!
   * \brief Get the storage for the ids that make up this shape.
   *
   * \return The container for the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityStorageConstRef getIdsStorage() const { return m_ids; }

  /*!
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityView getIds() const { return m_ids; }

  /*!
   * \brief Get the ids for the requested face.
   *
   * \param faceIndex The index of the desired face.
   * \param[out] ids A buffer that will contain the ids.
   * \param[out] numIds The number of ids returned for the face.
   *
   */
  AXOM_HOST_DEVICE void getFace(int AXOM_UNUSED_PARAM(faceIndex),
                                ConnectivityType *ids,
                                axom::IndexType &numIds) const
  {
    numIds = m_ids.size();
    for(axom::IndexType i = 0; i < numIds; i++)
    {
      ids[i] = m_ids[i];
    }
  }

  /*!
   * \brief Return the number of edges in the polygon.
   * \return The number of edges.
   */
  AXOM_HOST_DEVICE IndexType numberOfEdges() const { return m_ids.size(); }

  /*!
   * \brief Return the number of nodes in the face.
   * \return The number of nodes in the 1 face.
   */
  AXOM_HOST_DEVICE IndexType numberOfNodesInFace(int) const { return m_ids.size(); }

  /*!
   * \brief Return the specified edge.
   * \param edgeIndex The index of the edge to return.
   * \return The 2 point ids in the edge.
   */
  AXOM_HOST_DEVICE axom::StackArray<IndexType, 2> getEdge(int edgeIndex) const
  {
    const auto p0 = edgeIndex % m_ids.size();
    const auto p1 = (edgeIndex + 1) % m_ids.size();
    return axom::StackArray<IndexType, 2> {p0, p1};
  }

private:
  ConnectivityView m_ids;
};

/*!
 * \brief This class extends the ShapeTraits with object state so it can represent a zone.
 *
 * \tparam ShapeTraits A shape traits class from which to inherit.
 * \tparam ConnStorage A view or container that contains connectivity.
 *
 */
template <typename ShapeTraits, typename ConnStorage>
struct Shape : public ShapeTraits
{
  using ConnectivityStorage = ConnStorage;
  using ConnectivityStorageRef = ConnStorage &;
  using ConnectivityStorageConstRef = const ConnStorage &;
  using ConnectivityType = typename ConnStorage::value_type;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /*!
   * \brief Construct a shape.
   */
  AXOM_HOST_DEVICE Shape() : m_ids() { }

  /*!
   * \brief Construct a shape.
   *
   * \param ids A reference to connectivity storage for this shape.
   */
  AXOM_HOST_DEVICE Shape(ConnectivityStorageConstRef ids) : m_ids(ids)
  {
    SLIC_ASSERT(m_ids.size() == ShapeTraits::numberOfNodes());
  }

  /*!
   * \brief Get a specific id that makes up this shape.
   *
   * \return The i'th id that makes up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityType getId(size_t index) const
  {
    SLIC_ASSERT(index < static_cast<size_t>(m_ids.size()));
    return m_ids[index];
  }

  /*!
   * \brief Get the storage for the ids that make up this shape.
   *
   * \return The container for the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityStorageRef getIdsStorage() { return m_ids; }

  /*!
   * \brief Get the storage for the ids that make up this shape.
   *
   * \return The container for the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityStorageConstRef getIdsStorage() const { return m_ids; }

  /*!
   * \brief Get the ids that make up this shape as a view.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityView getIds() const
  {
    return ConnectivityView(const_cast<ConnectivityType *>(m_ids.data()), m_ids.size());
  }

  /*!
   * \brief Get the unique ids that make up this shape. For basic shapes, assume they are unique.
   *
   * \return The unique ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityStorageConstRef getUniqueIds() const { return m_ids; }

  /*!
   * \brief Get the ids for the requested face.
   *
   * \param faceIndex The index of the desired face.
   * \param[out] ids A buffer that will contain the ids. The buffer must be large enough
   *                 to accommodate the largest number of ids that could be returned, which
   *                 can be obtained in most cases from the maxNodesInFace() method.
   * \param[out] numIds The number of ids returned for the face.
   */
  /// @{
  template <int _ndims = ShapeTraits::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<(_ndims == 2), void>::type getFace(
    axom::IndexType AXOM_UNUSED_PARAM(faceIndex),
    ConnectivityType *ids,
    axom::IndexType &numIds) const
  {
    numIds = static_cast<axom::IndexType>(m_ids.size());
    for(axom::IndexType i = 0; i < numIds; i++)
    {
      ids[i] = m_ids[i];
    }
  }

  template <int _ndims = ShapeTraits::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims != 2, void>::type getFace(axom::IndexType faceIndex,
                                                                            ConnectivityType *ids,
                                                                            axom::IndexType &numIds) const
  {
    numIds = ShapeTraits::numberOfNodesInFace(faceIndex);
    const auto faceIds = ShapeTraits::getFace(faceIndex);
    for(IndexType i = 0; i < numIds; i++)
    {
      ids[i] = m_ids[faceIds[i]];
    }
  }
  /// @}

private:
  ConnectivityStorage m_ids;
};

/*!
 * \brief Some concrete shape classes based on the shape traits.
 *
 * \tparam ConnType A type of the connectivity values or a type that "stores"
 *                  the data either in actuality, or as a view. If an integral
 *                  type is passed, an axom::ArrayView will be used. 
 */
/// @{
template <typename ConnType>
using PointShape =
  Shape<PointTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;

template <typename ConnType>
using LineShape =
  Shape<LineTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;

template <typename ConnType>
using TriShape =
  Shape<TriTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;

template <typename ConnType>
using QuadShape =
  Shape<QuadTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;

template <typename ConnType>
using TetShape =
  Shape<TetTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;

template <typename ConnType>
using PyramidShape =
  Shape<PyramidTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;

template <typename ConnType>
using WedgeShape =
  Shape<WedgeTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;

template <typename ConnType>
using HexShape =
  Shape<HexTraits,
        typename std::conditional<std::is_integral<ConnType>::value, axom::ArrayView<ConnType>, ConnType>::type>;
/// @}

/*!
 * \brief This is a shape that can act as any of the other shapes.
 *
 * \tparam ConnType type of the connectivity values.
 *
 * \note This is a substitute for polymorphism so we can run on device.
 */
template <typename ConnType>
struct VariableShape
{
  using ConnectivityStorage = axom::ArrayView<ConnType>;
  using ConnectivityStorageRef = ConnectivityStorage &;
  using ConnectivityStorageConstRef = const ConnectivityStorage &;

  using ConnectivityType = ConnType;
  using ConnectivityView = ConnectivityStorage;

  /*!
   * \brief Constructor
   *
   * \param shapeId The shape id that describes the points.
   * \param ids The ids that describe the shape.
   */
  AXOM_HOST_DEVICE
  VariableShape(int shapeId, ConnectivityStorageConstRef ids) : m_shapeId(shapeId), m_ids(ids)
  {
    SLIC_ASSERT(shapeId >= Point_ShapeID && shapeId <= Hex_ShapeID);
  }

  /*!
   * \brief Returns the shape id of the actual shape represented by the variable shape.
   * \return The actual shape represented.
   */
  AXOM_HOST_DEVICE int id() const { return m_shapeId; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return true; }

  AXOM_HOST_DEVICE IndexType dimension() const { return shapeDimension(m_shapeId); }

  AXOM_HOST_DEVICE IndexType numberOfNodes() const { return m_ids.size(); }

  AXOM_HOST_DEVICE IndexType numberOfNodesInFace(int faceIndex) const
  {
    IndexType nnodes = 0;
    switch(m_shapeId)
    {
    case Point_ShapeID:
      nnodes = PointTraits::numberOfNodesInFace(faceIndex);
      break;
    case Line_ShapeID:
      nnodes = LineTraits::numberOfNodesInFace(faceIndex);
      break;
    case Tri_ShapeID:
      nnodes = TriTraits::numberOfNodesInFace(faceIndex);
      break;
    case Quad_ShapeID:
      nnodes = QuadTraits::numberOfNodesInFace(faceIndex);
      break;
    case Polygon_ShapeID:
      nnodes = (faceIndex == 0) ? m_ids.size() : 0;
      break;
    case Tet_ShapeID:
      nnodes = TetTraits::numberOfNodesInFace(faceIndex);
      break;
    case Pyramid_ShapeID:
      nnodes = PyramidTraits::numberOfNodesInFace(faceIndex);
      break;
    case Wedge_ShapeID:
      nnodes = WedgeTraits::numberOfNodesInFace(faceIndex);
      break;
    case Hex_ShapeID:
      nnodes = HexTraits::numberOfNodesInFace(faceIndex);
      break;
    }
    return nnodes;
  }

  AXOM_HOST_DEVICE static constexpr IndexType maxNodesInFace()
  {
    IndexType nnodes = PointTraits::maxNodesInFace();
    nnodes = std::max(nnodes, LineTraits::maxNodesInFace());
    nnodes = std::max(nnodes, TriTraits::maxNodesInFace());
    nnodes = std::max(nnodes, QuadTraits::maxNodesInFace());
    nnodes = std::max(nnodes, PolygonTraits::maxNodesInFace());
    nnodes = std::max(nnodes, TetTraits::maxNodesInFace());
    nnodes = std::max(nnodes, PyramidTraits::maxNodesInFace());
    nnodes = std::max(nnodes, WedgeTraits::maxNodesInFace());
    nnodes = std::max(nnodes, HexTraits::maxNodesInFace());
    return nnodes;
  }

  AXOM_HOST_DEVICE IndexType numberOfFaces() const
  {
    IndexType nfaces = 0;
    switch(m_shapeId)
    {
    case Point_ShapeID:
      nfaces = PointTraits::numberOfFaces();
      break;
    case Line_ShapeID:
      nfaces = LineTraits::numberOfFaces();
      break;
    case Tri_ShapeID:
      nfaces = TriTraits::numberOfFaces();
      break;
    case Quad_ShapeID:
      nfaces = QuadTraits::numberOfFaces();
      break;
    case Polygon_ShapeID:
      nfaces = 1;
      break;
    case Tet_ShapeID:
      nfaces = TetTraits::numberOfFaces();
      break;
    case Pyramid_ShapeID:
      nfaces = PyramidTraits::numberOfFaces();
      break;
    case Wedge_ShapeID:
      nfaces = WedgeTraits::numberOfFaces();
      break;
    case Hex_ShapeID:
      nfaces = HexTraits::numberOfFaces();
      break;
    }
    return nfaces;
  }

  AXOM_HOST_DEVICE void getFace(axom::IndexType faceIndex,
                                ConnectivityType *ids,
                                axom::IndexType &numIds) const
  {
    switch(m_shapeId)
    {
    case Point_ShapeID:
      // Falls through
    case Line_ShapeID:
      // Falls through
    case Tri_ShapeID:
      // Falls through
    case Quad_ShapeID:
      // Falls through
    case Polygon_ShapeID:
    {
      numIds = m_ids.size();
      for(axom::IndexType i = 0; i < numIds; i++)
      {
        ids[i] = m_ids[i];
      }
    }
    break;
    case Tet_ShapeID:
    {
      numIds = TetTraits::numberOfNodesInFace(faceIndex);
      const auto faceIds = TetTraits::getFace(faceIndex);
      for(IndexType i = 0; i < numIds; i++)
      {
        ids[i] = m_ids[faceIds[i]];
      }
    }
    break;
    case Pyramid_ShapeID:
    {
      numIds = PyramidTraits::numberOfNodesInFace(faceIndex);
      const auto faceIds = PyramidTraits::getFace(faceIndex);
      for(IndexType i = 0; i < numIds; i++)
      {
        ids[i] = m_ids[faceIds[i]];
      }
    }
    break;
    case Wedge_ShapeID:
    {
      numIds = WedgeTraits::numberOfNodesInFace(faceIndex);
      const auto faceIds = WedgeTraits::getFace(faceIndex);
      for(IndexType i = 0; i < numIds; i++)
      {
        ids[i] = m_ids[faceIds[i]];
      }
    }
    break;
    case Hex_ShapeID:
    {
      numIds = HexTraits::numberOfNodesInFace(faceIndex);
      const auto faceIds = HexTraits::getFace(faceIndex);
      for(IndexType i = 0; i < numIds; i++)
      {
        ids[i] = m_ids[faceIds[i]];
      }
    }
    break;
    }
  }

  AXOM_HOST_DEVICE IndexType numberOfEdges() const
  {
    IndexType nedges = 0;
    switch(m_shapeId)
    {
    case Point_ShapeID:
      nedges = PointTraits::numberOfEdges();
      break;
    case Line_ShapeID:
      nedges = LineTraits::numberOfEdges();
      break;
    case Tri_ShapeID:
      nedges = TriTraits::numberOfEdges();
      break;
    case Quad_ShapeID:
      nedges = QuadTraits::numberOfEdges();
      break;
    case Polygon_ShapeID:
      nedges = m_ids.size();
      break;
    case Tet_ShapeID:
      nedges = TetTraits::numberOfEdges();
      break;
    case Pyramid_ShapeID:
      nedges = PyramidTraits::numberOfEdges();
      break;
    case Wedge_ShapeID:
      nedges = WedgeTraits::numberOfEdges();
      break;
    case Hex_ShapeID:
      nedges = HexTraits::numberOfEdges();
      break;
    }
    return nedges;
  }

  AXOM_HOST_DEVICE axom::StackArray<IndexType, 2> getEdge(int edgeIndex) const
  {
    axom::StackArray<IndexType, 2> edge;
    switch(m_shapeId)
    {
    case Point_ShapeID:
      edge = PointTraits::getEdge(edgeIndex);
      break;
    case Line_ShapeID:
      edge = LineTraits::getEdge(edgeIndex);
      break;
    case Tri_ShapeID:
      edge = TriTraits::getEdge(edgeIndex);
      break;
    case Quad_ShapeID:
      edge = QuadTraits::getEdge(edgeIndex);
      break;
    case Polygon_ShapeID:
    {
      const auto n = m_ids.size();
      edge[0] = edgeIndex % n;
      edge[1] = (edgeIndex + 1) % n;
      break;
    }
    case Tet_ShapeID:
      edge = TetTraits::getEdge(edgeIndex);
      break;
    case Pyramid_ShapeID:
      edge = PyramidTraits::getEdge(edgeIndex);
      break;
    case Wedge_ShapeID:
      edge = WedgeTraits::getEdge(edgeIndex);
      break;
    case Hex_ShapeID:
      edge = HexTraits::getEdge(edgeIndex);
      break;
    }
    return edge;
  }

  /*!
   * \brief Get a specific id that makes up this shape.
   *
   * \return The i'th id that makes up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityType getId(IndexType index) const { return m_ids[index]; }

  /*!
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityView getIds() const { return m_ids; }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "mixed"; }

private:
  int m_shapeId;
  ConnectivityStorage m_ids;
};

/*!
 * \brief Given a shape name (matches Blueprint shape name), return the Shape id() value.
 *
 * \param name The shape name.
 *
 * \return The shape id that matches the name, or 0 if there is no match.
 */
inline int shapeNameToID(const std::string &name)
{
  int id = Invalid_ShapeID;
  if(name == PointTraits::name())
    id = Point_ShapeID;
  else if(name == LineTraits::name())
    id = Line_ShapeID;
  else if(name == TriTraits::name())
    id = Tri_ShapeID;
  else if(name == QuadTraits::name())
    id = Quad_ShapeID;
  else if(name == PolygonTraits::name())
    id = Polygon_ShapeID;
  else if(name == TetTraits::name())
    id = Tet_ShapeID;
  else if(name == PyramidTraits::name())
    id = Pyramid_ShapeID;
  else if(name == WedgeTraits::name())
    id = Wedge_ShapeID;
  else if(name == HexTraits::name())
    id = Hex_ShapeID;
  else if(name == "polyhedral")
    id = Polyhedron_ShapeID;
  else if(name == "mixed")
    id = Mixed_ShapeID;
  return id;
}

/*!
 * \brief Determine whether a value is a valid ShapeID.
 *
 * \param shapeID The value to test for validity.
 *
 * \return True if the value is a valid ShapeID; False otherwise.
 */
template <typename T>
AXOM_HOST_DEVICE constexpr bool isValidShapeID(T shapeID)
{
  return shapeID >= static_cast<T>(Point_ShapeID) && shapeID <= static_cast<T>(Mixed_ShapeID);
}

/*!
 * \brief Return the dimension for the specified shape id.
 *
 * \param shapeId The shape id.
 *
 * \return The dimension for the shape id.
 *
 * \note We can't tell the dimension for mixed shapes from the id alone. Return -1.
 */
AXOM_HOST_DEVICE constexpr IndexType shapeDimension(int shapeId)
{
  IndexType dim = 2;
  switch(shapeId)
  {
  case Point_ShapeID:
    dim = PointTraits::dimension();
    break;
  case Line_ShapeID:
    dim = LineTraits::dimension();
    break;
  case Tri_ShapeID:
    dim = TriTraits::dimension();
    break;
  case Quad_ShapeID:
    dim = QuadTraits::dimension();
    break;
  case Polygon_ShapeID:
    dim = PolygonTraits::dimension();
    break;
  case Tet_ShapeID:
    dim = TetTraits::dimension();
    break;
  case Pyramid_ShapeID:
    dim = PyramidTraits::dimension();
    break;
  case Wedge_ShapeID:
    dim = WedgeTraits::dimension();
    break;
  case Hex_ShapeID:
    dim = HexTraits::dimension();
    break;
  case Polyhedron_ShapeID:
    dim = 3;
    break;
  case Mixed_ShapeID:
    dim = -1;
    break;
  }
  return dim;
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
