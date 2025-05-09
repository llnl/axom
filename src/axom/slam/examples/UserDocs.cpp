// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file UserDocs.cpp
 *
 * \brief Example code to support slam user doc examples
 */

/* This example code contains snippets used in the Slam Sphinx documentation.
 * They begin and end with comments of the form
 *
 * <example-tag>_start
 * <example-tag>_end
 *
 * each prepended with an underscore.
 */

#include "axom/core.hpp"
#include "axom/slic.hpp"

// _quadmesh_example_import_header_start
#include "axom/slam.hpp"
// _quadmesh_example_import_header_end

#include <sstream>
#include <cmath>
#include <fstream>

// _quadmesh_example_slam_namespace_start
namespace slam = axom::slam;
// _quadmesh_example_slam_namespace_end

/**
 * \brief Simple 2D Point class for example
 */
struct Point2
{
  Point2(double x = 0., double y = 0.) : m_x(x), m_y(y) { }

  Point2(const Point2& other) : m_x(other.m_x), m_y(other.m_y) { }

  Point2& operator=(const Point2& other)
  {
    m_x = other.m_x;
    m_y = other.m_y;
    return *this;
  }

  Point2& operator+=(const Point2& other)
  {
    m_x += other.m_x;
    m_y += other.m_y;
    return *this;
  }

  Point2& operator/=(double val)
  {
    m_x /= val;
    m_y += val;
    return *this;
  }

  double& operator[](int i) { return (i == 0) ? m_x : m_y; }
  const double& operator[](int i) const { return (i == 0) ? m_x : m_y; }

  friend std::ostream& operator<<(std::ostream& os, const Point2& pt)
  {
    return os << "{x:" << pt.m_x << ", y:" << pt.m_y << "}";
  }

  double m_x, m_y;
};

/**
 * \brief Simple example class to define a star-shaped quad mesh with eleven
 * vertices and five elements.
 *
 * \note This code is in support of the slam user docs.
 */

struct SimpleQuadMesh
{
  /// Type aliases for sets
  // _quadmesh_example_set_typedefs_start
  using PosType = slam::DefaultPositionType;
  using ElemType = slam::DefaultElementType;
  using VertSet = slam::PositionSet<PosType, ElemType>;
  using ElemSet = slam::PositionSet<PosType, ElemType>;
  // _quadmesh_example_set_typedefs_end

  // _quadmesh_example_common_typedefs_start
  using ArrayIndir = slam::policies::CArrayIndirection<PosType, ElemType>;
  // _quadmesh_example_common_typedefs_end

  /// Type aliases for relations
  // _quadmesh_example_bdry_relation_typedefs_start
  // Type aliases for element-to-vertex boundary relation
  enum
  {
    VertsPerElem = 4
  };
  using CTStride = slam::policies::CompileTimeStride<PosType, VertsPerElem>;
  using ConstCard = slam::policies::ConstantCardinality<PosType, CTStride>;
  using ElemToVertRelation =
    slam::StaticRelation<PosType, ElemType, ConstCard, ArrayIndir, ElemSet, VertSet>;
  // _quadmesh_example_bdry_relation_typedefs_end

  // _quadmesh_example_cobdry_relation_typedefs_start
  // Type aliases for vertex-to-element coboundary relation
  using VarCard = slam::policies::VariableCardinality<PosType, ArrayIndir>;
  using VertToElemRelation =
    slam::StaticRelation<PosType, ElemType, VarCard, ArrayIndir, VertSet, ElemSet>;
  // _quadmesh_example_cobdry_relation_typedefs_end

  /// Type alias for position map
  // _quadmesh_example_maps_typedefs_start
  using BaseSet = slam::Set<PosType, ElemType>;
  using ScalarMap = slam::Map<Point2, BaseSet>;
  using PointMap = slam::Map<Point2, BaseSet>;
  using VertPositions = PointMap;
  // _quadmesh_example_maps_typedefs_end

  SimpleQuadMesh()
  {
    // data for element-vertex boundary relation
    evInds = {
      0, 3, 2,  1,  // elem 0
      0, 5, 4,  3,  // elem 1
      0, 7, 6,  5,  // elem 2
      0, 9, 8,  7,  // elem 3
      0, 1, 10, 9   // elem 4
    };

    // data for vertex-element coboundary relation
    veInds = {
      0, 1, 2, 3, 4,  // vert  0, card 5, start 0
      0, 4,           // vert  1, card 2, start 5
      0,              // vert  2, card 1, start 7
      0, 1,           // vert  3, card 2, start 8
      1,              // vert  4, card 1, start 10
      1, 2,           // vert  5, card 2, start 11
      2,              // vert  6, card 1, start 13
      2, 3,           // vert  7, card 2, start 14
      3,              // vert  8, card 1, start 16
      3, 4,           // vert  9, card 2, start 17
      4,              // vert 10, card 1, start 19, end 20
    };
    veBegins = {0, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20};
  }

  /// Constructs the mesh vertex and element sets
  void constructMeshSets()
  {
    // construct Sets
    // _quadmesh_example_construct_sets_start
    verts = VertSet(11);  // Construct vertex set with 11 vertices
    elems = ElemSet(5);   // Construct the element set with 5 elements
    // _quadmesh_example_construct_sets_end

    // _quadmesh_example_set_isvalid_start
    SLIC_ASSERT_MSG(verts.isValid(), "Vertex set is not valid.");
    SLIC_ASSERT_MSG(elems.isValid(), "Elment set is not valid.");
    // _quadmesh_example_set_isvalid_end

    SLIC_INFO("Mesh has " << verts.size() << " vertices");
    SLIC_INFO("Mesh has " << elems.size() << " elems");
  }

  /// Constructs the mesh boundary and coboundary relations
  void constructMeshRelations()
  {
    {
      // _quadmesh_example_construct_bdry_relation_start
      // construct boundary relation from elements to vertices
      using RelationBuilder = ElemToVertRelation::RelationBuilder;
      bdry = RelationBuilder().fromSet(&elems).toSet(&verts).indices(
        RelationBuilder::IndicesSetBuilder().size(static_cast<int>(evInds.size())).data(evInds.data()));
      // _quadmesh_example_construct_bdry_relation_end
    }

    {
      // _quadmesh_example_construct_cobdry_relation_start
      // construct coboundary relation from vertices to elements
      using RelationBuilder = VertToElemRelation::RelationBuilder;
      cobdry = RelationBuilder()
                 .fromSet(&verts)
                 .toSet(&elems)
                 .begins(RelationBuilder::BeginsSetBuilder().size(verts.size()).data(veBegins.data()))
                 .indices(RelationBuilder::IndicesSetBuilder()
                            .size(static_cast<int>(veInds.size()))
                            .data(veInds.data()));
      // _quadmesh_example_construct_cobdry_relation_end
    }

    SLIC_ASSERT_MSG(bdry.isValid(), "Boundary relation is not valid.");
    SLIC_ASSERT_MSG(cobdry.isValid(), "Coboundary relation is not valid.");

    SLIC_INFO("Elem-Vert relation has size " << bdry.totalSize());
    SLIC_INFO("Vert-Elem relation has size " << cobdry.totalSize());
  }

  /// Constucts the positions map on the vertices
  void constructMeshMaps()
  {
    // _quadmesh_example_vert_positions_start
    // construct the position map on the vertices
    position = VertPositions(&verts);

    // first vertex is at origin
    position[0] = Point2(0., 0.);

    // remaining vertices lie within annulus around unit disk
    // in cw order, starting at angleOffset
    constexpr double rInner = 0.8;
    constexpr double rOuter = 1.2;
    constexpr double angleOffset = 0.75;
    const double N = verts.size() - 1;

    for(int i = 1; i < verts.size(); ++i)
    {
      const double angle = -(i - 1) / N * 2 * M_PI + angleOffset;
      const double mag = axom::utilities::random_real(rInner, rOuter);

      position[i] = Point2(mag * std::cos(angle), mag * std::sin(angle));
    }
    // _quadmesh_example_vert_positions_end

    SLIC_ASSERT_MSG(position.isValid(), "Position map is not valid.");

    SLIC_INFO("-- Vertex positions:");
    for(int vID = 0; vID < verts.size(); ++vID)
    {
      SLIC_INFO("Position of vert " << vID << " is " << position[vID]);
    }
  }

  /// Computes the distance of each vertex from the origin
  void computeDistances()
  {
    // _quadmesh_example_vert_distances_start
    // Create a Map of scalars over the vertices
    ScalarMap distances(&verts);

    for(int i = 0; i < distances.size(); ++i)  // <-- Map::size()
    {
      auto vID = verts[i];               // <-- Set::operator[]
      const Point2& pt = position[vID];  // <-- Map::operator[]

      distances[i] = std::sqrt(pt[0] * pt[0]  // <-- Map::operator[]
                               + pt[1] * pt[1]);
    }
    // _quadmesh_example_vert_distances_end

    SLIC_INFO("-- Vertex distances to origin:");
    for(int i = 0; i < verts.size(); ++i)
    {
      auto vID = verts[i];
      SLIC_INFO("Distance of vert " << vID << ": " << distances[i]);
    }
  }

  /// Computes the centroid of each element using boundary relation
  void computeCentroids()
  {
    // _quadmesh_example_elem_centroids_start
    // Create a Map of Point2 over the mesh elements
    using ElemCentroidMap = PointMap;
    ElemCentroidMap centroid = ElemCentroidMap(&elems);

    // for each element...
    for(int eID = 0; eID < elems.size(); ++eID)  // <-- Set::size()
    {
      Point2 ctr;

      auto elVerts = bdry[eID];  // <-- Relation::operator[]

      // find average position of incident vertices
      for(int i = 0; i < elVerts.size(); ++i)  // <-- Set::size()
      {
        auto vID = elVerts[i];  // <-- Set::operator[]
        ctr += position[vID];   // <-- Map::operator[]
      }
      ctr /= elVerts.size();  // <-- Set::size())
      centroid[eID] = ctr;    // <-- Map::operator[]
    }
    // _quadmesh_example_elem_centroids_end

    SLIC_INFO("-- Element centroids:");
    for(int eID = 0; eID < elems.size(); ++eID)
    {
      SLIC_INFO("Centroid of elem " << eID << " is " << centroid[eID]);
    }
  }

  void outputVTKMesh()
  {
    // _quadmesh_example_output_vtk_start
    std::ofstream meshfile;
    meshfile.open("quadMesh.vtk");
    std::ostream_iterator<PosType> out_it(meshfile, " ");

    // write header
    meshfile << "# vtk DataFile Version 3.0\n"
             << "vtk output\n"
             << "ASCII\n"
             << "DATASET UNSTRUCTURED_GRID\n\n"
             << "POINTS " << verts.size() << " double\n";

    // write positions
    for(auto pos : position)  // <-- Uses range-based for on position map
    {
      meshfile << pos[0] << " " << pos[1] << " 0\n";
    }

    // write elem-to-vert boundary relation
    meshfile << "\nCELLS " << elems.size() << " " << 5 * elems.size();
    for(auto e : elems)  // <-- uses range-based for on element set
    {
      meshfile << "\n4 ";
      std::copy(bdry.begin(e),  // <-- uses relation's iterators
                bdry.end(e),
                out_it);
    }

    // write element types ( 9 == VKT_QUAD )
    meshfile << "\n\nCELL_TYPES " << elems.size() << "\n";
    for(int i = 0; i < elems.size(); ++i)
    {
      meshfile << "9 ";
    }

    // write element ids
    meshfile << "\n\nCELL_DATA " << elems.size() << "\nSCALARS cellIds int 1"
             << "\nLOOKUP_TABLE default \n";
    for(int i = 0; i < elems.size(); ++i)
    {
      meshfile << elems[i] << " ";  // <-- uses size() and operator[] on set
    }

    // write vertex ids
    meshfile << "\n\nPOINT_DATA " << verts.size() << "\nSCALARS vertIds int 1"
             << "\nLOOKUP_TABLE default \n";
    for(int i = 0; i < verts.size(); ++i)
    {
      meshfile << verts[i] << " ";
    }
    meshfile << "\n";
    // _quadmesh_example_output_vtk_end
  }

private:
  // _quadmesh_example_set_variables_start
  VertSet verts;  // The set of vertices in the mesh
  ElemSet elems;  // The set of elements in the mesh
  //  _quadmesh_example_set_variables_end

  // _quadmesh_example_relation_variables_start
  ElemToVertRelation bdry;    // Boundary relation from elements to vertices
  VertToElemRelation cobdry;  // Coboundary relation from vertices to elements
  // _quadmesh_example_relation_variables_end

  // _quadmesh_example_map_variables_start
  VertPositions position;  // vertex position
  // _quadmesh_example_map_variables_end

  // support data for mesh connectivity
  std::vector<PosType> evInds;
  std::vector<PosType> veInds;
  std::vector<PosType> veBegins;
};

void quadMeshExample()
{
  SimpleQuadMesh quadMesh;
  quadMesh.constructMeshSets();
  quadMesh.constructMeshRelations();
  quadMesh.constructMeshMaps();

  quadMesh.computeDistances();
  quadMesh.computeCentroids();

  quadMesh.outputVTKMesh();
}

int main()
{
  axom::slic::SimpleLogger logger;
  quadMeshExample();
}
