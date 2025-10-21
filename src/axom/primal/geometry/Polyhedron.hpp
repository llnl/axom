// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Polyhedron.hpp
 *
 * \brief A Polyhedron primitive for primal
 */

#ifndef AXOM_PRIMAL_POLYHEDRON_HPP_
#define AXOM_PRIMAL_POLYHEDRON_HPP_

#include "axom/core/StackArray.hpp"

#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class Polyhedron;

/*! \brief Overloaded output operator for polyhedrons */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Polyhedron<T, NDIMS>& poly);

/*!
 * \class NeighborCollection
 *
 * \brief Represents a collection of neighbor relations between vertices.
 */
class NeighborCollection
{
public:
  static constexpr int MAX_VERTS = 32;
  static constexpr int MAX_NBRS_PER_VERT = 8;

  using VertexNbrs = axom::StackArray<std::int8_t, MAX_NBRS_PER_VERT>;

public:
  /*!
   * \brief Constructs an empty NeighborCollection.
   */
  NeighborCollection() = default;

  /*!
   * \brief Clears the set of neighbors.
   */
  AXOM_HOST_DEVICE void clear()
  {
    for(int i = 0; i < MAX_VERTS; i++)
    {
      num_nbrs[i] = 0;
    }
  }

  /*!
   * \brief Returns the number of neighbors of a given vertex.
   *
   * \param [in] vtx index of the vertex to check.
   * \return num_nbrs the number of neighbors the vertex has.
   */
  AXOM_HOST_DEVICE int getNumNeighbors(int vtx) const
  {
#if !defined(__HIP_DEVICE_COMPILE__)
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
#endif
    return num_nbrs[vtx];
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE const VertexNbrs& operator[](int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return getNeighbors(vtx);
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE VertexNbrs& operator[](int vtx)
  {
#if !defined(__HIP_DEVICE_COMPILE__)
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
#endif
    return getNeighbors(vtx);
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE const VertexNbrs& getNeighbors(int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return nbrs[vtx];
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE VertexNbrs& getNeighbors(int vtx)
  {
#if !defined(__HIP_DEVICE_COMPILE__)
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
#endif
    return nbrs[vtx];
  }

  /*!
   * \brief Adds a set of neighbors to a given vertex.
   *
   * \param [in] vtx the index of the vertex to add new neighbors to
   * \param [in] nbr a set of indices of neighboring vertices
   *
   * \pre vtx < MAX_VERTS
   */
  AXOM_HOST_DEVICE void addNeighbors(std::int8_t vtx, std::initializer_list<std::int8_t> nbrIds)
  {
#if !defined(__HIP_DEVICE_COMPILE__)
    SLIC_ASSERT(num_nbrs[vtx] + nbrIds.size() <= MAX_NBRS_PER_VERT);
#endif
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    for(std::int8_t nbr : nbrIds)
    {
      std::int8_t idx_insert = num_nbrs[vtx];
      nbrs[vtx][idx_insert] = nbr;
      num_nbrs[vtx]++;
    }
  };

  /*!
   * \brief Adds a single neighbor to a given vertex.
   *
   * \param [in] vtx the index of the vertex to add the new neighbor to
   * \param [in] nbr an index of a neighboring vertex
   *
   * \pre vtx < MAX_VERTS
   */
  AXOM_HOST_DEVICE void addNeighbors(std::int8_t vtx, std::int8_t nbrId)
  {
#if !defined(__HIP_DEVICE_COMPILE__)
    SLIC_ASSERT(num_nbrs[vtx] + 1 <= MAX_NBRS_PER_VERT);
#endif
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    std::int8_t idx_insert = num_nbrs[vtx];
    nbrs[vtx][idx_insert] = nbrId;
    num_nbrs[vtx]++;
  };

  /*!
   * \brief Inserts a new neighbor for a vertex at a given index in the
   *  neighbor list.
   *
   * \param [in] vtx the index of the vertex to add a new neighbor to
   * \param [in] nbr the index of the neighboring vertex
   * \param [in] pos the position in the neighbor list to insert into
   *
   * \pre vtx < MAX_VERTS
   * \pre pos <= num_nbrs[vtx]
   */
  AXOM_HOST_DEVICE void insertNeighborAtPos(std::int8_t vtx, std::int8_t nbr, std::int8_t pos)
  {
#if !defined(__HIP_DEVICE_COMPILE__)
    SLIC_ASSERT(num_nbrs[vtx] + 1 <= MAX_NBRS_PER_VERT);
    SLIC_ASSERT(pos <= num_nbrs[vtx]);
#endif

    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);

    std::uint8_t old_nbrs[MAX_NBRS_PER_VERT];
    // copy elements from [pos, nnbrs)
    for(int ip = pos; ip < num_nbrs[vtx]; ip++)
    {
      old_nbrs[ip - pos] = nbrs[vtx][ip];
    }
    nbrs[vtx][pos] = nbr;
    for(int ip = pos; ip < num_nbrs[vtx]; ip++)
    {
      nbrs[vtx][ip + 1] = old_nbrs[ip - pos];
    }
    num_nbrs[vtx]++;
  }

  /*!
   * \brief Compacts all neighbor sets by removing neighbors with index -1.
   */
  AXOM_HOST_DEVICE void pruneNeighbors()
  {
    for(int iv = 0; iv < MAX_VERTS; iv++)
    {
      int curr_idx = 0;
      for(int inbr = 0; inbr < num_nbrs[iv]; inbr++)
      {
        if(nbrs[iv][inbr] == -1)
        {
          continue;
        }
        nbrs[iv][curr_idx] = nbrs[iv][inbr];
        curr_idx++;
      }
      num_nbrs[iv] = curr_idx;
    }
  }

private:
  std::int8_t num_nbrs[MAX_VERTS];
  VertexNbrs nbrs[MAX_VERTS];
};

/*!
 * \class Polyhedron
 *
 * \brief Represents a convex polyhedron defined by an array of points and
 *        optionally their neighbors (a point is a neighbor of another point if
 *        there is an edge between them)
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * \note The Polyhedron functions do not check that points defining a face are
 *       coplanar. It is the responsibility of the caller to pass a
 *       valid set of points and neighbors representing a convex polyhedron.
 *
 * \note The polyhedron neighbors should be ordered counter clockwise
 *       as when viewing the polyhedron externally.
 *
 *       <pre>
 *
 *          3--------2          +y
 *         /|       /|
 *        / |      / |           ^
 *       7--------6  |           |
 *       |  0-----|--1           |
 *       | /      | /            -----> +x
 *       |/       |/            /
 *       4--------5            /
 *                            <
 *                           +z
 *
 *       </pre>
 *
 *       For example, vertex 5 above should have neighbors (4,1,6).
 *
 * \note The Polyhedron functions do not check that neighbors are ordered
 *       counter clockwise. It is the responsibility of the caller to pass a
 *       valid neighbors ordering.
 */
template <typename T, int NDIMS = 3>
class Polyhedron
{
public:
  constexpr static int MAX_VERTS = NeighborCollection::MAX_VERTS;
  constexpr static int MAX_PLANES = MAX_VERTS;

  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using NumArrayType = axom::NumericArray<T, NDIMS>;
  using PlaneType = Plane<T, NDIMS>;
  using PlaneArrayType = StackArray<PlaneType, MAX_PLANES>;

private:
  using Coords = StackArray<PointType, MAX_VERTS>;
  using Neighbors = NeighborCollection;

public:
  /*! Default constructor for an empty polyhedron   */
  Polyhedron() = default;

  /*! Return the number of vertices in the polyhedron */
  AXOM_HOST_DEVICE int numVertices() const { return m_num_vertices; }

  /*!
   * \brief Appends a vertex to the list of vertices.
   *
   * \return The index where the vertex was inserted into.
   */
  AXOM_HOST_DEVICE int addVertex(const PointType& pt)
  {
    SLIC_ASSERT(m_num_vertices + 1 < MAX_VERTS);
    m_vertices[m_num_vertices] = pt;
    m_num_vertices++;
    return m_num_vertices - 1;
  }

  /*!
   * \brief Stores nbrs as the neighbors of vertex pt.
   *        If vertex pt is not in the polyhedron, neighbors are not added.
   *
   * \note pt should be a vertex of this polyhedron. Otherwise, pt may not
   *       be found due to floating point precision differences.
   *
   * \param [in] pt The vertex to add neighbors
   * \param [in] nbrs The neighbors to add to the list of neighbors
   */
  AXOM_HOST_DEVICE
  void addNeighbors(const PointType& pt, std::initializer_list<std::int8_t> nbrs)
  {
    for(int i = 0; i < m_num_vertices; i++)
    {
      if(m_vertices[i] == pt)
      {
        m_neighbors.addNeighbors(i, nbrs);
      }
    }
  }

  /*!
   * \brief Stores nbrs as the neighbors of vertex at a given index.
   *
   * \note Caller is responsible for ensuring a given vertex has enough space
   *       remaining to fit the new neighbors.
   *
   * \param [in] vtxId The vertex id to add neighbors
   * \param [in] nbrs The neighbors to add to the list of neighbors
   *
   * \pre vtxId < getVertices()
   */
  AXOM_HOST_DEVICE
  void addNeighbors(int vtxId, std::initializer_list<std::int8_t> nbrs)
  {
    m_neighbors.addNeighbors(vtxId, nbrs);
  }

  /*!
   * \brief Adds a single nbr as the neighbor of vertex at a given index.
   *
   * \note Caller is responsible for ensuring a given vertex has enough space
   *       remaining to fit the new neighbor.
   *
   * \param [in] vtxId The vertex id to add the neighbor
   * \param [in] nbr The neighbor to add to the list of neighbors
   *
   * \pre vtxId < getVertices()
   */
  AXOM_HOST_DEVICE
  void addNeighbors(int vtxId, int nbr) { m_neighbors.addNeighbors(vtxId, nbr); }

  /*! Clears the list of vertices and neighbors */
  AXOM_HOST_DEVICE void clear()
  {
    m_num_vertices = 0;
    m_neighbors.clear();
  }

  /*! Retrieves the vertices */
  Coords& getVertices() { return m_vertices; }

  /*! Retrieves the vertex at index idx */
  AXOM_HOST_DEVICE PointType& operator[](int idx) { return m_vertices[idx]; }
  /*! Retrieves the vertex at index idx */
  AXOM_HOST_DEVICE const PointType& operator[](int idx) const { return m_vertices[idx]; }

  /*! Retrieves the neighbors */
  AXOM_HOST_DEVICE Neighbors& getNeighbors() { return m_neighbors; }

  AXOM_HOST_DEVICE int getNumNeighbors(int i) const { return m_neighbors.getNumNeighbors(i); }

  /*! Retrieves the neighbors for vertex i */
  AXOM_HOST_DEVICE
  typename Neighbors::VertexNbrs& getNeighbors(int i) { return m_neighbors[i]; }
  /*! Retrieves the neighbors for vertex i */
  AXOM_HOST_DEVICE
  const typename Neighbors::VertexNbrs& getNeighbors(int i) const { return m_neighbors[i]; }

  /*!
   * \brief Computes the vertex mean as the average of the polyhedron's vertex
   *  positions
   *
   * \return The vertex mean of the polyhedron's vertices
   *
   * \pre  polyhedron.isValid() is true
   */
  AXOM_HOST_DEVICE
  PointType vertexMean() const
  {
    SLIC_ASSERT(isValid());

    NumArrayType sum;

    const int nVerts = numVertices();
    if(nVerts > 0)
    {
      for(int i = 0; i < nVerts; ++i)
      {
        sum += m_vertices[i].array();
      }
      sum /= nVerts;
    }
    return PointType(sum);
  }

  /*!
   * \brief Helper function to find the faces of the Polyhedron, assuming the
   *        vertex neighbors are in counter-clockwise ordering.
   *
   * \tparam ConnectivityType The type to generate for the polyhedral face data.
   *
   * \param [out] faces is the vertex indices for faces
   * \param [out] face_offset is the offset for each face
   * \param [out] face_size is the number of vertices for each face
   * \param [out] face_count is the number of faces
   *
   * \warning Function is experimental, input parameters and/or output may
   *          change in the future.
   *
   * \note This function does not check that the provided buffers are
   *       large enough to hold all values. It is the responsibility of the
   *       caller to know ahead of time and pass buffers with appropriate size.
   *
   * \note Function is based off extractFaces() in Mike Owen's PolyClipper.
   *
   * \pre polyhedron vertex neighbors are defined
   */
  template <typename ConnectivityType = int>
  AXOM_HOST_DEVICE void getFaces(ConnectivityType* faces,
                                 ConnectivityType* face_size,
                                 ConnectivityType* face_offset,
                                 axom::IndexType& face_count) const
  {
    std::int8_t curFaceIndex = 0;
    std::int8_t checkedSize = 0;
    std::int8_t facesAdded = 0;
    // # edges * (# vertices per edge) * (# orientation per edge)
    std::int8_t checkedEdges[MAX_VERTS * 2 * 2] = {0};

    // Check each vertex
    for(int i = 0; i < numVertices(); ++i)
    {
      // Check each neighbor index (edge) of the vertex
      for(int j = 0; j < m_neighbors.getNumNeighbors(i); j++)
      {
        // Check if edge has not been visited
        int ni = m_neighbors[i][j];
        int edgeToCheck[2] = {i, ni};

        bool alreadyChecked = false;

        for(int edge = 0; edge < checkedSize; edge++)
        {
          if(edgeToCheck[0] == checkedEdges[edge * 2] &&
             edgeToCheck[1] == checkedEdges[(edge * 2) + 1])
          {
            alreadyChecked = true;
            break;
          }
        }

        if(!alreadyChecked)
        {
          face_offset[facesAdded] = curFaceIndex;
          faces[curFaceIndex++] = i;
          std::int8_t curFaceSize = 1;
          std::int8_t vstart = i;
          std::int8_t vnext = ni;
          std::int8_t vprev = i;

          // Add neighboring vertices until we reach the starting vertex.
          while(vnext != vstart)
          {
            faces[curFaceIndex++] = vnext;
            curFaceSize++;
            checkedEdges[checkedSize * 2] = vprev;
            checkedEdges[(checkedSize * 2) + 1] = vnext;
            checkedSize++;

            int numNeighbors = m_neighbors.getNumNeighbors(vnext);
            int itr = -1;
            for(int k = 0; k < m_neighbors.getNumNeighbors(vnext); k++)
            {
              if(m_neighbors[vnext][k] == vprev)
              {
                itr = k;
                break;
              }
            }
            vprev = vnext;
            if(itr == 0)
            {
              vnext = m_neighbors[vnext][numNeighbors - 1];
            }
            else
            {
              vnext = m_neighbors[vnext][itr - 1];
            }
          }  // end of while loop

          // Add last edge connecting the face
          checkedEdges[checkedSize * 2] = vprev;
          checkedEdges[(checkedSize * 2) + 1] = vnext;
          checkedSize++;

          face_size[facesAdded++] = curFaceSize++;
        }
      }
    }
    face_count = facesAdded;
  }

  /*!
   * \brief Get the polyhedron's faces as planes.
   *
   * \param[out] The number of faces.
   * \param divideNonPlanarFaces Whether to divide non-planar faces into separate planes.
   * \param eps The tolerance used to decide whether points are close enough to a plane.
   *
   * \return An array of planes that describe the faces.
   */
  AXOM_HOST_DEVICE PlaneArrayType getFaces(axom::IndexType& numFaces,
                                           bool divideNonPlanarFaces = true,
                                           double eps = 1.e-8) const
  {
    PlaneArrayType facePlanes;
    numFaces = 0;

    // Get faces
    int faces[MAX_VERTS * MAX_VERTS];
    int face_size[MAX_VERTS * 2];
    int face_offset[MAX_VERTS * 2];
    axom::IndexType face_count;
    getFaces<int>(faces, face_size, face_offset, face_count);

    // Turn the faces to planes.
    numFaces = 0;
    for(axom::IndexType i = 0; i < face_count; ++i)
    {
      const int i_offset = face_offset[i];
      const auto p0 = m_vertices[faces[i_offset]];
      const auto p1 = m_vertices[faces[i_offset + 1]];
      const auto p2 = m_vertices[faces[i_offset + 2]];

      // Make a plane.
      auto plane = axom::primal::make_plane(p0, p1, p2);

      // We've already made 1 plane for the face. Treat the face as a triangle
      // fan and add the rest of the planes, if they are different enough from
      // the previous plane.
      if(divideNonPlanarFaces)
      {
        facePlanes[numFaces] = plane;
        numFaces = (numFaces + 1 < MAX_PLANES) ? (numFaces + 1) : numFaces;

        const int extraPlanes = face_size[i] - 3;
        for(int ep = 0; ep < extraPlanes; ep++)
        {
          // Get the last point in the triangle.
          const auto ep2 = m_vertices[faces[i_offset + 3 + ep]];

          // Check its signed distance vs the previous plane. If the point is far
          // enough away, we make the plane.
          const auto dist = plane.signedDistance(ep2);
          if(axom::utilities::abs(dist) > eps)
          {
            const auto ep1 = m_vertices[faces[i_offset + 2 + ep]];
            const auto plane2 = axom::primal::make_plane(p0, ep1, ep2);

            facePlanes[numFaces] = plane2;
            numFaces = (numFaces + 1 < MAX_PLANES) ? (numFaces + 1) : numFaces;

            plane = plane2;
          }
        }
      }
      else
      {
        facePlanes[numFaces] = plane;
        numFaces = (numFaces + 1 < MAX_PLANES) ? (numFaces + 1) : numFaces;
      }
    }
    return facePlanes;
  }

  /*!
   * \brief Computes the moments of the polyhedron. The 0th moment is the
   *        volume of the polyhedron, the 1st moment is the centroid.
   *
   * \param [out] volume The volume of the polyhedron (0th moment)
   * \param [out] centroid The centroid of the polyhedron (1st moment)
   * \param [in]  should_compute_centroid If true, computes the centroid
   *              of the polyhedron. Default is true.
   *
   * \note Function is based off moments() in Mike Owen's PolyClipper.
   *
   * \pre polyhedron vertex neighbors are defined, and polyhedron is 3D
   *
   * \sa volume()
   * \sa centroid()
   */
  AXOM_HOST_DEVICE
  void moments(double& volume, PointType& centroid, bool should_compute_centroid = true) const
  {
    volume = 0.0;

    VectorType centroid_vector;

    if(!isValid() || hasDuplicateVertices())
    {
      return;
    }

    // Computes the signed volume of tetrahedra formed from vertices of the
    // Polyhedron faces and an arbitrary origin (the first vertex), as well
    // as the centroid of the Polyhedron.
    else
    {
      SLIC_CHECK_MSG(hasNeighbors(), "Polyhedron::moments() is only valid with vertex neighbors.");

      // faces is an overestimation
      int faces[MAX_VERTS * MAX_VERTS];
      int face_size[MAX_VERTS * 2];
      int face_offset[MAX_VERTS * 2];
      axom::IndexType face_count;
      getFaces(faces, face_size, face_offset, face_count);

      const PointType& origin = m_vertices[0];

      for(axom::IndexType i = 0; i < face_count; ++i)
      {
        const int N = face_size[i];
        const int i_offset = face_offset[i];
        const VectorType v0 = m_vertices[faces[i_offset]] - origin;

        for(int j = 1, k = 2; j < N - 1; ++j, ++k)
        {
          VectorType v1 = m_vertices[faces[i_offset + j]] - origin;
          VectorType v2 = m_vertices[faces[i_offset + k]] - origin;
          const T curVol = VectorType::scalar_triple_product(v0, v1, v2);

          volume += curVol;
          if(should_compute_centroid)
          {
            centroid_vector += (v0 + v1 + v2) * curVol;
          }
        }
      }

      volume /= 6.;

      if(should_compute_centroid)
      {
        centroid_vector /= (volume != 0.0) ? (24.0 * volume) : axom::primal::PRIMAL_TINY;
        centroid = centroid_vector + origin;
      }
    }
  }

  /*!
   * \brief Computes the signed volume of the polyhedron.
   *
   * \return The signed volume of the polyhedron
   *
   * \pre polyhedron vertex neighbors are defined, and polyhedron is 3D
   *
   * \sa volume()
   */
  AXOM_HOST_DEVICE
  double signedVolume() const
  {
    double volume;
    PointType centroid;
    moments(volume, centroid, false);
    return volume;
  }

  /*!
   * \brief Returns the absolute (unsigned) volume of the polyhedron
   * \sa signedVolume()
   */
  AXOM_HOST_DEVICE
  double volume() const { return axom::utilities::abs(signedVolume()); }

  /*!
   * \brief Computes the centroid as the center of mass of the polyhedron
   *
   * \return The centroid of the polyhedron
   *
   * \pre  polyhedron.isValid() is true
   */
  AXOM_HOST_DEVICE
  PointType centroid() const
  {
    double volume;
    PointType centroid;
    moments(volume, centroid);
    return centroid;
  }

  /*!
   * \brief Simple formatted print of a polyhedron instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int sz = numVertices();

    os << "{" << sz << "-vertex polyhedron:\n";
    for(int i = 0; i < sz - 1; ++i)
    {
      os << "Vertex " << i << ": " << m_vertices[i] << "\n";
      if(hasNeighbors() && m_neighbors.getNumNeighbors(i) > 0)
      {
        os << "Neighbors of vertex " << i << ": ";
        int nz = m_neighbors.getNumNeighbors(i);
        for(int j = 0; j < nz; j++)
        {
          int curN = getNeighbors(i)[j];
          os << curN << " ";
        }
        os << "\n";
      }
    }
    if(sz >= 2)
    {
      os << "Vertex " << sz - 1 << ": " << m_vertices[sz - 1] << "\n";
      if(hasNeighbors() && m_neighbors.getNumNeighbors(sz - 1) > 0)
      {
        os << "Neighbors of vertex " << sz - 1 << ": ";
        int nz = m_neighbors.getNumNeighbors(sz - 1);
        for(int j = 0; j < nz; j++)
        {
          int curN = getNeighbors(sz - 1)[j];
          os << curN << " ";
        }
        os << "\n";
      }
    }
    os << "}";

    return os;
  }

  /*!
   * \brief Simple check for validity of a polyhedron
   *
   * Initial check is that the polyhedron has four or more vertices
   * \return True, if the polyhedron is valid, False otherwise
   */
  AXOM_HOST_DEVICE
  bool isValid() const { return m_num_vertices >= 4; }

  /*!
   * \brief Check if any vertices are duplicate, within a tolerance
   *
   * \param [in] eps The tolerance
   *
   * \return True, if the polyhedron has duplicate vertices, False otherwise
   */
  AXOM_HOST_DEVICE
  bool hasDuplicateVertices(double eps = 1.e-10) const
  {
    const auto typedEPS = static_cast<T>(eps);
    for(int i = 0; i < m_num_vertices; i++)
    {
      for(int j = i + 1; j < m_num_vertices; j++)
      {
        if(axom::utilities::isNearlyEqual(m_vertices[i][0], m_vertices[j][0], typedEPS) &&
           axom::utilities::isNearlyEqual(m_vertices[i][1], m_vertices[j][1], typedEPS) &&
           axom::utilities::isNearlyEqual(m_vertices[i][2], m_vertices[j][2], typedEPS))
        {
          return true;
        }
      }
    }
    return false;
  }

  /*!
   * \brief Check if vertex neighbors information is available
   *
   * Checks if neighbors is not empty
   * \return True, if the polyhedron has neighbors info for each vertex,
   *         False otherwise
   */
  AXOM_HOST_DEVICE
  bool hasNeighbors() const
  {
    if(m_num_vertices <= 1)
    {
      return false;
    }
    bool has_nbrs = true;
    for(int i = 0; i < m_num_vertices; i++)
    {
      has_nbrs = has_nbrs && (m_neighbors.getNumNeighbors(i) > 0);
    }
    return has_nbrs;
  }

  /*!
   * \brief Creates a Polyhedron from a given Hexahedron's vertices.
   *
   * \param [in] hex The hexahedron
   * \param [in] tryFixOrientation If true, checks if the signed volume of the
   *             Polyhedron is negative and swaps the order of some vertices
   *             in that shape to try to obtain a nonnegative signed volume.
   *             Defaults to false.
   *
   * \return A Polyhedron with the Hexahedron's vertices and added
   *         vertex neighbors
   *
   * \note The Hexahedron is assumed to have a specific vertex order:
   * \verbatim
   *
   *          4--------7          +z
   *         /|       /|               +y
   *        / |      / |           ^  >
   *       5--------6  |           | /
   *       |  0-----|--3           |/
   *       | /      | /            -----> +x
   *       |/       |/
   *       1--------2
   *
   * \endverbatim
   *
   *       The Polyhedron's vertex neighbors are created assuming this vertex
   *       ordering.
   *
   * \warning tryFixOrientation flag does not guarantee the Polyhedron's vertex order
   *          will be valid. It is the responsiblity of the caller to pass
   *          a Hexahedron with a valid vertex order. Otherwise, if the
   *          Hexahedron has an invalid vertex order, the returned Polyhedron
   *          will have a non-positive and/or unexpected volume.
   *
   * \warning If tryFixOrientation flag is false and some of the shapes have
   *          a negative signed volume, the returned Polyhedron
   *          will have a non-positive and/or unexpected volume.
   *
   */
  AXOM_HOST_DEVICE
  static Polyhedron from_primitive(const Hexahedron<T, NDIMS>& hex, bool tryFixOrientation = false)
  {
    // Initialize our polyhedron to return
    Polyhedron<T, NDIMS> poly;

    poly.addVertex(hex[0]);
    poly.addVertex(hex[1]);
    poly.addVertex(hex[2]);
    poly.addVertex(hex[3]);
    poly.addVertex(hex[4]);
    poly.addVertex(hex[5]);
    poly.addVertex(hex[6]);
    poly.addVertex(hex[7]);

    poly.addNeighbors(0, {1, 4, 3});
    poly.addNeighbors(1, {0, 2, 5});
    poly.addNeighbors(2, {1, 3, 6});
    poly.addNeighbors(3, {2, 0, 7});
    poly.addNeighbors(4, {0, 5, 7});
    poly.addNeighbors(5, {1, 6, 4});
    poly.addNeighbors(6, {2, 7, 5});
    poly.addNeighbors(7, {3, 4, 6});

    // Reverses order of vertices 1,3 and 5,7 if signed volume is negative
    if(tryFixOrientation)
    {
      if(poly.signedVolume() < 0)
      {
        axom::utilities::swap<Point<T, NDIMS>>(poly[1], poly[3]);
        axom::utilities::swap<Point<T, NDIMS>>(poly[5], poly[7]);
      }
    }

    return poly;
  }

  /*!
   * \brief Creates a Polyhedron from a given Octahedron's vertices.
   *
   * \param [in] oct The octahedron
   * \param [in] tryFixOrientation If true, checks if the signed volume of the
   *             Polyhedron is negative and swaps the order of some vertices
   *             in that shape to try to obtain a nonnegative signed volume.
   *             Defaults to false.
   *
   * \return A Polyhedron with the Octahedron's vertices and added
   *         vertex neighbors
   *
   * \note The Octahedron is assumed to have a specific vertex order:
   *       (view looking down from +z axis):
   *
   * \verbatim
   *
   *            4                +z
   *            /\                    +y
   *       0 --/  \-- 2           ^  >
   *         \/    \ /            | /
   *         /      \             |/
   *       5 -------- 3           -----> +x
   *            \/
   *            1
   *
   * \endverbatim
   *
   *       The Polyhedron's vertex neighbors are created assuming this vertex
   *       ordering.
   *
   * \warning tryFixOrientation flag does not guarantee the Polyhedron's vertex order
   *          will be valid. It is the responsiblity of the caller to pass
   *          an Octahedron with a valid vertex order. Otherwise, if the
   *          Octahedron has an invalid vertex order, the returned Polyhedron
   *          will have a non-positive and/or unexpected volume.
   *
   * \warning If tryFixOrientation flag is false and some of the shapes have
   *          a negative signed volume, the returned Polyhedron
   *          will have a non-positive and/or unexpected volume.
   *
   */
  AXOM_HOST_DEVICE
  static Polyhedron from_primitive(const Octahedron<T, NDIMS>& oct, bool tryFixOrientation = false)
  {
    // Initialize our polyhedron to return
    Polyhedron<T, NDIMS> poly;

    poly.addVertex(oct[0]);
    poly.addVertex(oct[1]);
    poly.addVertex(oct[2]);
    poly.addVertex(oct[3]);
    poly.addVertex(oct[4]);
    poly.addVertex(oct[5]);

    poly.addNeighbors(0, {1, 5, 4, 2});
    poly.addNeighbors(1, {0, 2, 3, 5});
    poly.addNeighbors(2, {0, 4, 3, 1});
    poly.addNeighbors(3, {1, 2, 4, 5});
    poly.addNeighbors(4, {0, 5, 3, 2});
    poly.addNeighbors(5, {0, 1, 3, 4});

    // Reverses order of vertices 1,2 and 4,5 if volume is negative.
    if(tryFixOrientation)
    {
      if(poly.signedVolume() < 0)
      {
        axom::utilities::swap<PointType>(poly[1], poly[2]);
        axom::utilities::swap<PointType>(poly[4], poly[5]);
      }
    }

    return poly;
  }

  /*!
   * \brief Creates a Polyhedron from a given Tetrahedron's vertices.
   *
   * \param [in] tet The tetrahedron
   * \param [in] tryFixOrientation If true, checks if the signed volume of the
   *             Polyhedron is negative and swaps the order of some vertices
   *             in that shape to try to obtain a nonnegative signed volume.
   *             Defaults to false.
   *
   * \return A Polyhedron with the Tetrahedron's vertices and added
   *         vertex neighbors
   *
   * \note The Tetrahedron is assumed to have a specific vertex order:
   * \verbatim
   *
   *              3                    +z
   *             / \\                       +y
   *            /   \ \                 ^  >
   *           /     \  \               | /
   *          /       \   \             |/
   *         /         \    2           -----> +x
   *        /           \  /
   *       /_____________\/
   *      0               1
   *
   * \endverbatim
   *
   *       The Polyhedron's vertex neighbors are created assuming this vertex
   *       ordering.
   *
   * \warning tryFixOrientation flag does not guarantee the Polyhedron's vertex
   *          order will be valid. It is the responsiblity of the caller to
   *          pass a Tetrahedron with a valid vertex order. Otherwise, if the
   *          Tetrahedron has an invalid vertex order, the returned Polyhedron
   *          will have a non-positive and/or unexpected volume.
   *
   * \warning If tryFixOrientation flag is false and some of the shapes have
   *          a negative signed volume, the returned Polyhedron
   *          will have a non-positive and/or unexpected volume.
   *
   */
  AXOM_HOST_DEVICE
  static Polyhedron from_primitive(const Tetrahedron<T, NDIMS>& tet, bool tryFixOrientation = false)
  {
    // Initialize our polyhedron to return
    Polyhedron<T, NDIMS> poly;

    poly.addVertex(tet[0]);
    poly.addVertex(tet[1]);
    poly.addVertex(tet[2]);
    poly.addVertex(tet[3]);

    poly.addNeighbors(0, {1, 3, 2});
    poly.addNeighbors(1, {0, 2, 3});
    poly.addNeighbors(2, {0, 3, 1});
    poly.addNeighbors(3, {0, 1, 2});

    // Reverses order of vertices 1 and 2 if signed volume is negative
    if(tryFixOrientation)
    {
      if(tet.signedVolume() < 0)
      {
        axom::utilities::swap<PointType>(poly[1], poly[2]);
      }
    }

    return poly;
  }

private:
  int m_num_vertices {0};
  Coords m_vertices;
  Neighbors m_neighbors {};
};

//------------------------------------------------------------------------------
/// Free functions implementing Polyhedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Polyhedron<T, NDIMS>& poly)
{
  poly.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Polyhedron using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Polyhedron<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_POLYHEDRON_HPP_
