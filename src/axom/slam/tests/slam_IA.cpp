// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_IA.cpp
 *
 * \brief This file contains tests for slam's IA mesh
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/mesh_struct/IA.hpp"
#include "axom/slam/Utilities.hpp"

namespace slam = axom::slam;

namespace
{
using PointType = slam::util::Point3<double>;
using IndexType = axom::IndexType;

struct BasicTriMeshData
{
  //This is the surface of a regular cube composed of 12 triangles
  std::vector<double> points;
  std::vector<IndexType> elem;

  static const double point_arr[];
  static const IndexType tri_arr[];
  static const IndexType vert_to_el_num[];  // for verification purpose
  static const IndexType el_nbr_rel[];

  BasicTriMeshData() : points(point_arr, point_arr + 24), elem(tri_arr, tri_arr + 36) { }

  static int numVertices() { return 8; }
  static int numTriangles() { return 12; }
};

// clang-format off
const double BasicTriMeshData::point_arr [] = {
  0.0, 0.0, 0.0,
  0.0, 0.0, 1.0,
  0.0, 1.0, 0.0,
  0.0, 1.0, 1.0,
  1.0, 0.0, 0.0,
  1.0, 0.0, 1.0,
  1.0, 1.0, 0.0,
  1.0, 1.0, 1.0
};
const IndexType BasicTriMeshData::tri_arr [] = {
  0, 6, 4,
  0, 2, 6,
  0, 3, 2,
  2, 7, 6,
  2, 3, 7,
  4, 6, 7,
  4, 7, 5,
  0, 4, 5,
  0, 1, 3,
  0, 5, 1,
  1, 5, 7,
  1, 7, 3
};
const IndexType BasicTriMeshData::vert_to_el_num [] = {
  6, 4, 4, 4, 4, 4, 4, 6
};
const IndexType BasicTriMeshData::el_nbr_rel [] = {
  1,5,7,
  2,3,0,
  8,4,1,
  4,5,1,
  2,11,3,
  0,3,6,
  5,10,7,
  0,6,9,
  9,11,2,
  7,10,8,
  9,6,11,
  10,4,8
};
// clang-format on

struct BasicTetMeshData
{
  //cube divided into six tets

  std::vector<double> points;
  std::vector<IndexType> elem;

  static const double point_arr[];
  static const IndexType tet_arr[];
  static const IndexType vert_to_el_num[];
  static const IndexType el_nbr_rel[];

  BasicTetMeshData() : points(point_arr, point_arr + 24), elem(tet_arr, tet_arr + 24) { }

  static int numVertices() { return 8; }
  static int numTetrahedra() { return 6; }
};

// clang-format off
const double BasicTetMeshData::point_arr[] = {
  -1.0, -1.0, -1.0,
  -1.0, -1.0,  1.0,
  -1.0,  1.0, -1.0,
  -1.0,  1.0,  1.0,
  1.0, -1.0, -1.0,
  1.0, -1.0,  1.0,
  1.0,  1.0, -1.0,
  1.0,  1.0,  1.0
};
const IndexType BasicTetMeshData::tet_arr [] = {
  3,2,4,0,
  3,1,4,0,
  3,6,2,4,
  3,6,7,4,
  3,5,1,4,
  3,5,7,4
};
const IndexType BasicTetMeshData::vert_to_el_num [] = {
  2, 2, 2, 6, 6, 2, 2, 2
};
const IndexType BasicTetMeshData::el_nbr_rel [] = {
  2,-1,1,-1,
  4,-1,0,-1,
  -1,-1,0,3,
  -1,-1,5,2,
  -1,-1,1,5,
  -1,-1,3,4
};
// clang-format on

/* Check if vertex with ID \a vert_id is in boundary of
 * element with ID \a elem_id in mesh \a ia_mesh
 */
template <typename IAMeshType>
bool isInBoundary(const IAMeshType& ia_mesh, IndexType vert_id, IndexType elem_id)
{
  for(auto bdry : ia_mesh.boundaryVertices(elem_id))
  {
    if(bdry == vert_id)
    {
      return true;
    }
  }
  return false;
}

/* Check if element with ID \a el_1 is adjacent to
 * element with ID \a el_2 in mesh \a ia_mesh
 */
template <typename IAMeshType>
bool isAdjacent(const IAMeshType& ia_mesh, IndexType el_1, IndexType el_2)
{
  for(auto adj : ia_mesh.adjacentElements(el_1))
  {
    if(adj == el_2)
    {
      return true;
    }
  }
  return false;
}

}  // end anonymous namespace

TEST(slam_IA, empty_mesh)
{
  SLIC_INFO("Testing creating an empty IA mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  IAMeshType ia_mesh;  //empty mesh

  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_TRUE(ia_mesh.isEmpty());
  EXPECT_EQ(0, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(0, ia_mesh.getNumberOfValidVertices());

  EXPECT_FALSE(ia_mesh.isValidElement(0));
  EXPECT_FALSE(ia_mesh.isValidVertex(0));

  // Check coboundary relation on non-existent vertex 0 (w/ warning)
  {
    int sz = ia_mesh.vertexStar(0).size();
    EXPECT_EQ(0, sz);
  }

  // Removing invalid vertices/elements is a no-op (w/ warning)
  ia_mesh.removeElement(0);
  ia_mesh.removeVertex(0);
}

TEST(slam_IA, basic_tri_mesh)
{
  SLIC_INFO("Testing constructing basic triangle mesh...");

  constexpr int TDIM = 2;
  constexpr int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;
  using IndexArray = IAMeshType::IndexArray;

  constexpr int vert_per_elem = IAMeshType::VERTS_PER_ELEM;
  constexpr int coord_per_vert = IAMeshType::COORDS_PER_VERT;
  constexpr int numAdjacentElems = vert_per_elem;

  EXPECT_EQ(TDIM + 1, vert_per_elem);
  EXPECT_EQ(SDIM, coord_per_vert);

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_FALSE(ia_mesh.isEmpty());

  const int numVerts = ia_mesh.getNumberOfValidVertices();
  EXPECT_EQ(basic_mesh_data.numVertices(), numVerts);

  const int numElems = ia_mesh.getNumberOfValidElements();
  EXPECT_EQ(basic_mesh_data.numTriangles(), numElems);

  // Test Element boundary relation
  for(auto e : ia_mesh.elements())
  {
    auto bdry = ia_mesh.boundaryVertices(e);
    EXPECT_EQ(vert_per_elem, bdry.size());
    for(auto it = bdry.begin(); it != bdry.end(); ++it)
    {
      EXPECT_EQ(basic_mesh_data.elem[e * vert_per_elem + it.index()], *it);
    }
  }

  // Test Vertex co-boundary relation
  for(auto v : ia_mesh.vertices())
  {
    IndexArray star_elems = ia_mesh.vertexStar(v);
    EXPECT_EQ(basic_mesh_data.vert_to_el_num[v], star_elems.size());

    for(auto el : star_elems)
    {
      EXPECT_TRUE(isInBoundary(ia_mesh, v, el))
        << "Vertex co-boundary relation indicates that " << v
        << " should be in boundary of element " << el << " but it is not";
    }
  }

  // Test Element adjacency relation
  for(auto e : ia_mesh.elements())
  {
    auto neighbors = ia_mesh.adjacentElements(e);
    EXPECT_EQ(numAdjacentElems, neighbors.size());
    for(auto nbr : neighbors)
    {
      EXPECT_TRUE(isAdjacent(ia_mesh, nbr, e))
        << "Element adjacency relation indicates that " << nbr << " should be adjacent to element "
        << e << " but it is not";
    }
  }

  ia_mesh.print_all();

  SLIC_INFO("Done");
}

TEST(slam_IA, dynamically_build_tri_mesh)
{
  SLIC_INFO("Testing dynamically modifying a triangle mesh...");

  BasicTriMeshData basic_mesh_data;
  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  const int vert_per_elem = IAMeshType::VERTS_PER_ELEM;
  EXPECT_EQ(3, vert_per_elem);

  const int coord_per_vert = IAMeshType::COORDS_PER_VERT;
  EXPECT_EQ(3, coord_per_vert);

  const int numAdjacentElems = vert_per_elem;

  // Build the mesh from nothing
  IAMeshType ia_mesh;  //empty mesh
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_TRUE(ia_mesh.isEmpty());

  // Add the vertices
  for(int v_i = 0; v_i < basic_mesh_data.numVertices(); ++v_i)
  {
    EXPECT_EQ(v_i, ia_mesh.getNumberOfValidVertices());

    PointType pt(&basic_mesh_data.points[v_i * coord_per_vert]);
    ia_mesh.addVertex(pt);
    EXPECT_TRUE(ia_mesh.isValid());
  }

  const int numVerts = ia_mesh.getNumberOfValidVertices();
  EXPECT_EQ(basic_mesh_data.numVertices(), numVerts);
  EXPECT_FALSE(ia_mesh.isEmpty());

  // Add the elements
  for(int e_i = 0; e_i < basic_mesh_data.numTriangles(); ++e_i)
  {
    EXPECT_EQ(e_i, ia_mesh.getNumberOfValidElements());

    int startIdx = e_i * vert_per_elem;
    ia_mesh.addElement(basic_mesh_data.elem[startIdx + 0],
                       basic_mesh_data.elem[startIdx + 1],
                       basic_mesh_data.elem[startIdx + 2]);

    EXPECT_TRUE(ia_mesh.isValid());
  }
  const int numTris = ia_mesh.getNumberOfValidElements();
  EXPECT_EQ(basic_mesh_data.numTriangles(), numTris);
  EXPECT_FALSE(ia_mesh.isEmpty());

  // Check that the element boundary relation is correct
  for(auto e : ia_mesh.elements())
  {
    auto bdry = ia_mesh.boundaryVertices(e);
    EXPECT_EQ(vert_per_elem, bdry.size());

    for(auto idx : bdry.positions())
    {
      const int expVert = basic_mesh_data.elem[e * vert_per_elem + idx];
      EXPECT_EQ(expVert, bdry[idx]);
    }
  }

  // Check that the vertex co-boundary relation is correct
  for(auto v : ia_mesh.vertices())
  {
    IndexType el = ia_mesh.coboundaryElement(v);
    EXPECT_TRUE(isInBoundary(ia_mesh, v, el))
      << "Vertex co-boundary relation indicates that " << v << " should be in boundary of element "
      << el << " but it is not";
  }

  // Check the element adjacencies are correct
  for(auto e : ia_mesh.elements())
  {
    auto neighbors = ia_mesh.adjacentElements(e);
    EXPECT_EQ(numAdjacentElems, neighbors.size());
    for(auto nbr : neighbors)
    {
      EXPECT_TRUE(isAdjacent(ia_mesh, nbr, e))
        << "Element adjacency relation indicates that " << nbr << " should be adjacent to element "
        << e << " but it is not";
    }
  }
}

TEST(slam_IA, tri_mesh_remove_verts_and_elems)
{
  SLIC_INFO("Testing removing elements and vertices from a triangle mesh...");

  constexpr int TDIM = 2;
  constexpr int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;
  constexpr int vert_per_elem = IAMeshType::VERTS_PER_ELEM;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  // Check that we begin with the correct number of verts and tris
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles(), ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  //removing elements bigger than 7
  ia_mesh.removeElement(8);
  ia_mesh.removeElement(9);
  ia_mesh.removeElement(10);
  ia_mesh.removeElement(11);

  EXPECT_FALSE(ia_mesh.isEmpty());

  EXPECT_TRUE(ia_mesh.isValid(true));
  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfValidElements());

  //removing vertex 1, which the removed elements contain
  ia_mesh.removeVertex(1);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_FALSE(ia_mesh.isEmpty());

  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfValidVertices());

  // Check the validity of the set entries
  for(auto v_i : ia_mesh.vertices().positions())
  {
    EXPECT_EQ(ia_mesh.isValidVertex(v_i), v_i != 1);
  }

  for(auto e_i : ia_mesh.elements().positions())
  {
    EXPECT_EQ(ia_mesh.isValidElement(e_i), e_i < 8);
  }

  // Check incidence and adjacency relations
  for(auto e_idx : ia_mesh.elements().positions())
  {
    auto neighbors = ia_mesh.adjacentElements(e_idx);
    EXPECT_EQ(vert_per_elem, neighbors.size());
    for(int n_idx : neighbors.positions())
    {
      int orig_nbr = basic_mesh_data.el_nbr_rel[e_idx * vert_per_elem + n_idx];
      if(orig_nbr > 7 || e_idx > 7)
      {
        EXPECT_EQ(neighbors[n_idx], IAMeshType::ElementAdjacencyRelation::INVALID_INDEX);
      }
      else
      {
        EXPECT_EQ(neighbors[n_idx], orig_nbr);
      }
    }
  }

  //removing vertex 0
  ia_mesh.removeVertex(0);
  EXPECT_TRUE(ia_mesh.isValid());

  //check that the elements with the removed vertices 0 are also removed.
  for(auto e_idx : ia_mesh.elements().positions())
  {
    bool bDeleted = false;
    for(int j = 0; j < vert_per_elem; j++)
    {
      bDeleted |= basic_mesh_data.elem[e_idx * vert_per_elem + j] == 0;
    }
    EXPECT_EQ(ia_mesh.isValidElement(e_idx), !bDeleted && e_idx <= 7);
  }
}

TEST(slam_IA, tri_mesh_remove_elem_and_compact)
{
  SLIC_INFO("Testing removing an element and compacting a triangle mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles(), ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  ia_mesh.removeElement(4);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 1, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 1, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());
}

TEST(slam_IA, tri_mesh_remove_vert_and_compact)
{
  SLIC_INFO("Testing removing a vertex and compacting a triangle mesh...");

  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTriMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles(), ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  ia_mesh.removeVertex(3);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfValidVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTriangles() - 4, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfValidVertices());
}

TEST(slam_IA, basic_tet_mesh)
{
  SLIC_INFO("Testing constructing basic tetrahedral mesh...");

  constexpr int vert_per_elem = 4;
  constexpr int coord_per_vert = 3;
  using IAMeshType = slam::IAMesh<vert_per_elem - 1, coord_per_vert, PointType>;

  BasicTetMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_EQ(ia_mesh.VERTS_PER_ELEM, vert_per_elem);
  EXPECT_EQ(ia_mesh.COORDS_PER_VERT, coord_per_vert);

  EXPECT_TRUE(ia_mesh.isValid());

  EXPECT_EQ(basic_mesh_data.numTetrahedra(), ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  // check EV relation
  for(auto el_i : ia_mesh.elements().positions())
  {
    auto bdry = ia_mesh.boundaryVertices(el_i);
    EXPECT_EQ(vert_per_elem, bdry.size());
    for(auto idx : bdry.positions())
    {
      EXPECT_EQ(bdry[idx], basic_mesh_data.elem[el_i * vert_per_elem + idx]);
    }
  }

  // check VE relation -- expanding to full star of each vertex
  for(auto vert_i : ia_mesh.vertices().positions())
  {
    auto star_elems = ia_mesh.vertexStar(vert_i);
    EXPECT_EQ(basic_mesh_data.vert_to_el_num[vert_i], star_elems.size());
    for(auto el_i : star_elems)
    {
      bool bContains = false;
      for(int j = 0; j < vert_per_elem; j++)
      {
        if(basic_mesh_data.elem[el_i * vert_per_elem + j] == vert_i)
        {
          bContains = true;
          break;
        }
      }
      EXPECT_TRUE(bContains);
    }
  }

  ia_mesh.print_all();

  SLIC_INFO("Done");
}

TEST(slam_IA, dynamically_build_tet_mesh)
{
  SLIC_INFO("Testing dynamically modifying a tetrahedral mesh...");

  BasicTetMeshData basic_mesh_data;
  const int vert_per_elem = 4;
  const int coord_per_vert = 3;

  using IAMeshType = slam::IAMesh<vert_per_elem - 1, coord_per_vert, PointType>;
  IAMeshType ia_mesh;  //empty mesh

  EXPECT_TRUE(ia_mesh.isValid());

  //build the mesh from nothing

  //adding the vertices
  for(int vert_i = 0; vert_i < (int)basic_mesh_data.points.size(); vert_i += coord_per_vert)
  {
    PointType pt(&basic_mesh_data.points[vert_i]);
    ia_mesh.addVertex(pt);
    EXPECT_TRUE(ia_mesh.isValid());
  }

  //adding the elements
  for(int elem_i = 0; elem_i < (int)basic_mesh_data.elem.size(); elem_i += vert_per_elem)
  {
    ia_mesh.addElement(basic_mesh_data.elem[elem_i],
                       basic_mesh_data.elem[elem_i + 1],
                       basic_mesh_data.elem[elem_i + 2],
                       basic_mesh_data.elem[elem_i + 3]);

    EXPECT_TRUE(ia_mesh.isValid());
  }

  //check that the Element-Vertex boundary entries are correct
  for(auto e : ia_mesh.elements())
  {
    auto e_verts = ia_mesh.boundaryVertices(e);
    EXPECT_EQ(e_verts.size(), vert_per_elem);

    for(auto idx : e_verts.positions())
    {
      EXPECT_EQ(e_verts[idx], basic_mesh_data.elem[e * vert_per_elem + idx]);
    }
  }

  // Check that the Vertex-Element coboundary entries are correct
  for(auto v : ia_mesh.vertices())
  {
    auto e_i = ia_mesh.coboundaryElement(v);
    bool bContains = false;
    for(int j = 0; j < vert_per_elem; ++j)
    {
      bContains |= basic_mesh_data.elem[e_i * vert_per_elem + j];
    }
    EXPECT_TRUE(bContains);
  }

  //check that the Element-Element adjacency entries are correct
  for(auto e : ia_mesh.elements())
  {
    auto neighbors = ia_mesh.adjacentElements(e);
    EXPECT_EQ(vert_per_elem, neighbors.size());
    for(auto idx : neighbors.positions())
    {
      EXPECT_EQ(neighbors[idx], basic_mesh_data.el_nbr_rel[e * vert_per_elem + idx]);
    }
  }

  //removing elements bigger than 3
  ia_mesh.removeElement(4);
  ia_mesh.removeElement(5);

  //removing vertex 1, which only the removed elements contain
  ia_mesh.removeVertex(5);

  ia_mesh.isValid();

  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 2, ia_mesh.getNumberOfValidElements());

  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfValidVertices());

  //check the validity of the set entries
  for(auto idx : ia_mesh.elements().positions())
  {
    EXPECT_EQ(ia_mesh.isValidElement(idx), idx < 4);
  }

  //check the ee_rel entries are correct after removing the elements
  for(auto idx : ia_mesh.elements().positions())
  {
    auto neighbors = ia_mesh.adjacentElements(idx);
    EXPECT_EQ(vert_per_elem, neighbors.size());
    for(auto n_idx : neighbors.positions())
    {
      int orig_nbr = basic_mesh_data.el_nbr_rel[idx * vert_per_elem + n_idx];
      if(orig_nbr < 4 && idx < 4)
      {
        EXPECT_EQ(neighbors[n_idx], orig_nbr);
      }
      else
      {
        EXPECT_EQ(neighbors[n_idx], IAMeshType::ElementAdjacencyRelation::INVALID_INDEX);
      }
    }
  }

  //removing vertex 6
  ia_mesh.removeVertex(6);

  //check that the elements with the removed vertices are also removed.
  for(int i = 0; i < basic_mesh_data.numTetrahedra(); ++i)
  {
    bool bDeleted = false;
    for(int j = 0; j < vert_per_elem; ++j)
    {
      bDeleted |= basic_mesh_data.elem[i * vert_per_elem + j] == 6;
    }
    EXPECT_EQ(ia_mesh.isValidElement(i), !bDeleted && i < 4);
  }

  SLIC_INFO("Done");
}

TEST(slam_IA, compact_mesh)
{
  const int TDIM = 2;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;
  using IndexType = IAMeshType::IndexType;

  IAMeshType mesh;

  mesh.addVertex(PointType(0, 0, 1));
  mesh.addVertex(PointType(0, 1, 0));
  mesh.addVertex(PointType(1, 0, 0));
  mesh.addVertex(PointType(1, 1, 1));
  mesh.addElement(0, 1, 2);
  mesh.addElement(1, 2, 3);

  mesh.removeVertex(0);

  EXPECT_TRUE(mesh.isValid());

  auto v_before = mesh.getNumberOfValidVertices();
  auto e_before = mesh.getNumberOfValidElements();
  EXPECT_EQ(IndexType(3), v_before);
  EXPECT_EQ(IndexType(1), e_before);

  mesh.compact();

  auto v_after = mesh.getNumberOfValidVertices();
  auto e_after = mesh.getNumberOfValidElements();

  EXPECT_TRUE(mesh.isValid(true));

  EXPECT_EQ(v_before, v_after);
  EXPECT_EQ(e_before, e_after);
}

TEST(slam_IA, tet_mesh_remove_elem_and_compact)
{
  SLIC_INFO("Testing removing an element and compacting a tetrahedral mesh...");

  const int TDIM = 3;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTetMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra(), ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  ia_mesh.removeElement(4);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 1, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 1, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());
}

TEST(slam_IA, tet_mesh_remove_vert_and_compact)
{
  SLIC_INFO("Testing removing a vertex and compacting a tetrahedral mesh...");

  const int TDIM = 3;
  const int SDIM = 3;
  using IAMeshType = slam::IAMesh<TDIM, SDIM, PointType>;

  BasicTetMeshData basic_mesh_data;
  IAMeshType ia_mesh(basic_mesh_data.points, basic_mesh_data.elem);

  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra(), ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices(), ia_mesh.getNumberOfValidVertices());

  ia_mesh.removeVertex(2);
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 2, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfValidVertices());

  ia_mesh.compact();
  EXPECT_TRUE(ia_mesh.isValid());
  EXPECT_EQ(basic_mesh_data.numTetrahedra() - 2, ia_mesh.getNumberOfValidElements());
  EXPECT_EQ(basic_mesh_data.numVertices() - 1, ia_mesh.getNumberOfValidVertices());
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
