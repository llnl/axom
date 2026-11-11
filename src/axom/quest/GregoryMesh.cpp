// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file GregoryMesh.cpp
 * \brief Build Gregory surfaces from a Blueprint mesh and write sampled VTK output.
 */

#include "axom/quest/GregoryMesh.hpp"

#include "axom/core/ArrayView.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/slic.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "GregoryMesh usage requires Conduit"
#endif

#include "conduit_blueprint.hpp"
#include "conduit_node.hpp"

#include <array>
#include <cstddef>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <utility>
#include <vector>

namespace axom::quest
{
namespace
{
///@{
//! @name Methods in an anonymous namespace to facilitate GregoryMesh initialization

using VectorType = GregoryMesh::VectorType;
using PointType = GregoryMesh::PointType;

/*!
 * \brief Append one facet's Blueprint connectivity to the construction records.
 *
 * The original coordinate indices are retained for later normal construction
 * and adjacency checks.
 *
 * \param[in,out] facet_node_ids The list of facet connectivities for this topology.
 * \param[in] connectivity Input Blueprint element connectivity.
 * \param[in] offset First connectivity entry for this facet.
 * \param[in] num_nodes Number of vertices in the facet.
 * \param[in] num_points Number of points in the Blueprint coordinate set.
 *
 * \return True after appending the facet.
 */
bool appendFacet(axom::Array<axom::Array<int>>& facet_node_ids,
                 const conduit::index_t_accessor& connectivity,
                 conduit::index_t offset,
                 int num_nodes)
{
  // Keep original coordinate indices for normal sharing and later adjacency checks.
  axom::Array<int> this_facet_ids(0, num_nodes);
  for(int i = 0; i < num_nodes; ++i)
  {
    const conduit::index_t node_id = connectivity[offset + i];
    this_facet_ids.push_back(static_cast<int>(node_id));
  }

  facet_node_ids.push_back(std::move(this_facet_ids));
  return true;
}

/*!
 * \brief Decode one Blueprint topology into temporary triangle and quad index arrays.
 *
 * Handles either fixed-shape or mixed connectivity during GregoryMesh construction.
 * An unsupported shape or mixed element causes the topology/surface to be skipped.
 *
 * \param[in] topology Blueprint topology to decode.
 * \param[in] num_points Number of points in the shared Blueprint coordinate set.
 *
 * \return An array of facet node IDs for a single surface/topology
 */
axom::Array<axom::Array<int>> readFacets(const conduit::Node& topology)
{
  const std::string name = topology.name();
  const conduit::Node& elements = topology.fetch_existing("elements");
  const std::string shape = elements.fetch_existing("shape").as_string();
  const auto connectivity = elements.fetch_existing("connectivity").as_index_t_accessor();
  axom::Array<axom::Array<int>> facet_node_ids;

  if(shape == "mixed")
  {
    // Only needed for mixed topology
    const auto shapes = elements.fetch_existing("shapes").as_index_t_accessor();
    const auto sizes = elements.fetch_existing("sizes").as_index_t_accessor();
    const auto offsets = elements.fetch_existing("offsets").as_index_t_accessor();

    // Mixed elements supply a shape and connectivity offset for each facet.
    const conduit::Node& shape_map = elements.fetch_existing("shape_map");
    const conduit::index_t tri_shape =
      shape_map.has_child("tri") ? shape_map.fetch_existing("tri").to_index_t() : -1;
    const conduit::index_t quad_shape =
      shape_map.has_child("quad") ? shape_map.fetch_existing("quad").to_index_t() : -1;

    facet_node_ids.reserve(shapes.number_of_elements());
    for(conduit::index_t i = 0; i < shapes.number_of_elements(); ++i)
    {
      const int num_nodes = shapes[i] == tri_shape ? 3 : shapes[i] == quad_shape ? 4 : 0;
      if(num_nodes == 0 || sizes[i] != num_nodes)
      {
        SLIC_WARNING("Topology '" << name << "' has an unsupported mixed element. Skipping.");
        return {};
      }

      appendFacet(facet_node_ids, connectivity, offsets[i], num_nodes);
    }
  }
  else
  {
    // Fixed-shape connectivity stores facet node IDs consecutively in one array.
    const int num_nodes = shape == "tri" ? 3 : shape == "quad" ? 4 : 0;
    if(num_nodes == 0 || connectivity.number_of_elements() % num_nodes != 0)
    {
      SLIC_WARNING("Topology '" << name << "' has an unsupported shape or connectivity. Skipping.");
      return {};
    }

    facet_node_ids.reserve(connectivity.number_of_elements() / num_nodes);
    for(conduit::index_t offset = 0; offset < connectivity.number_of_elements(); offset += num_nodes)
    {
      appendFacet(facet_node_ids, connectivity, offset, num_nodes);
    }
  }

  return facet_node_ids;
}

/*!
 * \brief Compute one facet corner's unit normal for GregoryMesh construction.
 *
 * Triangle facets use the same face normal at each corner.
 * Quadrilateral facets use adjacent edges (i.e., the normals from a bilinear interpolant).
 * All non-zero normals are unitized before returning.
 *
 * The result contributes to the average at the corresponding surface vertex.
 *
 * \param[in] node_ids Coordinate indices for the triangle or quadrilateral.
 * \param[in] local_node_id Corner whose normal is requested.
 * \param[in] x X coordinates of the shared coordinate set.
 * \param[in] y Y coordinates of the shared coordinate set.
 * \param[in] z Z coordinates of the shared coordinate set.
 *
 * \return Unit corner normal, or a zero vector when the local geometry is degenerate.
 */
VectorType localNormal(const axom::Array<int>& node_ids,
                       int local_node_id,
                       const conduit::double_accessor& x,
                       const conduit::double_accessor& y,
                       const conduit::double_accessor& z)
{
  const int num_nodes = static_cast<int>(node_ids.size());

  // Fix triangles to use node 0 as corner so each gets a numerically identical face normal.
  //  Note: it's a little redundant to compute the same thing multiple times for the same
  //  triangular facet, but it isn't too expensive and makes for an easier execution
  const int corner = num_nodes == 3 ? 0 : local_node_id;

  // Find the IDs of the corner and adjacent vertices
  const int p_id = node_ids[corner];
  const int next_id = node_ids[(corner + 1) % num_nodes];
  const int prev_id = node_ids[(corner + num_nodes - 1) % num_nodes];

  // Find the coordinates of the three coordinates
  const PointType p {x[p_id], y[p_id], z[p_id]};
  const PointType next {x[next_id], y[next_id], z[next_id]};
  const PointType prev {x[prev_id], y[prev_id], z[prev_id]};

  // Compute the normal, but ignore it if it's zero
  const VectorType normal = VectorType::cross_product(VectorType(p, next), VectorType(p, prev));

  return (normal.squared_norm() < primal::PRIMAL_TINY) ? VectorType {0.0, 0.0, 0.0}
                                                       : normal.unitVector();
}

/*!
 * \brief Build one named surface from its temporary facet node IDs.
 *
 * "Averages" corner normals at shared coordinate indices within this surface,
 * then uses them to construct the owned Gregory triangles and patches.
 *
 * Nonzero corner normals are normalized after averaging. A zero sum remains zero
 * and is passed to the Gregory primitive for its degenerate geometry handling.
 *
 * \param[in] name Blueprint topology name for the surface.
 * \param[in] facet_node_ids Original coordinate indices for each facet.
 * \param[in] x X coordinates of the shared coordinate set.
 * \param[in] y Y coordinates of the shared coordinate set.
 * \param[in] z Z coordinates of the shared coordinate set.
 * \return Named surface containing the constructed Gregory facets.
 */
GregorySurface buildSurface(const std::string& name,
                            const axom::Array<axom::Array<int>>& facet_node_ids,
                            const conduit::double_accessor& x,
                            const conduit::double_accessor& y,
                            const conduit::double_accessor& z)
{
  GregorySurface surface;

  surface.name = name;
  surface.facets.reserve(facet_node_ids.size());

  // Find the sum of all incident unit normal vectors at each vertex
  axom::FlatMap<int, VectorType> normal_sums;
  for(const axom::Array<int>& node_ids : facet_node_ids)
  {
    // For each facet, iterate through the IDs of its nodes
    for(int i = 0; i < node_ids.size(); ++i)
    {
      // Insert a zero vector into this node's position, if not already there
      auto inserted = normal_sums.emplace(node_ids[i], VectorType {0.0, 0.0, 0.0});

      const VectorType normal = localNormal(node_ids, i, x, y, z);
      inserted.first->second += normal;
    }
  }

  // Build each facet using the averaged normals from the first pass.
  for(const axom::Array<int>& node_ids : facet_node_ids)
  {
    const int num_nodes = static_cast<int>(node_ids.size());
    GregoryFacet facet;
    facet.shape = num_nodes == 3 ? GregoryFacetShape::Tri : GregoryFacetShape::Quad;

    axom::Array<PointType> corners(0, num_nodes);
    axom::Array<VectorType> normals(0, num_nodes);
    for(int i = 0; i < num_nodes; ++i)
    {
      const int node_id = node_ids[i];

      // Keep track of source IDs for later adjacency checks.
      facet.node_ids[i] = node_id;
      corners.push_back(PointType {x[node_id], y[node_id], z[node_id]});

      // Get normal from the map
      const auto normal_it = normal_sums.find(node_id);
      VectorType normal = normal_it->second;
      if(normal.squared_norm() < primal::PRIMAL_TINY)
      {
        // Add a zero vector directly instead of normalizing it
        normals.push_back(VectorType {0.0, 0.0, 0.0});
      }
      else
      {
        // Set nonzero facet geometry (normal)
        normals.push_back(normal.unitVector());
      }
    }

    // Build geometry from the temporary corner arrays.
    if(num_nodes == 3)
    {
      facet.geometry = primal::GregoryTriangle<double>(corners, normals);
    }
    else
    {
      facet.geometry = primal::GregoryPatch<double>(corners, normals);
    }
    surface.facets.push_back(std::move(facet));
  }

  return surface;
}
///@}

}  // namespace

/*!
 * \brief Build named Gregory surfaces from a Blueprint mesh.
 *
 * Uses the first Blueprint domain and an explicit 3D coordinate set shared by
 * supported unstructured tri, quad, or mixed tri/quad topologies. Each topology
 * becomes one surface, with corner normals averaged over its shared node IDs.
 * Unsupported topologies and topologies using another coordinate set are skipped;
 * degenerate facets retain their original node IDs and geometry.
 *
 * Invalid Blueprint input and a missing usable coordinate set are reported through
 * SLIC_ERROR. If abort on error is disabled, construction returns an empty mesh.
 *
 * \param[in] mesh Input Conduit Blueprint mesh.
 */
GregoryMesh::GregoryMesh(const conduit::Node& mesh)
{
  conduit::Node info;

  // Use the first domain; all named surface topologies must share one explicit coordset.
  const auto domains = conduit::blueprint::mesh::domains(mesh);
  if(domains.empty())
  {
    SLIC_ERROR("Blueprint mesh contains no domains");
    return;
  }
  if(domains.size() > 1)
  {
    SLIC_WARNING(
      "Quest: More than one domain observed. The GregoryMesh will be constructed "
      "only from domain 0, and will ignore the rest.");
  }

  const conduit::Node& domain = *domains[0];
  if(!conduit::blueprint::mesh::verify(domain, info))
  {
    SLIC_ERROR("Blueprint verification failed for domain 0: " << info.to_yaml());
    return;
  }

  const conduit::Node& topologies = domain.fetch_existing("topologies");
  if(topologies.number_of_children() == 0)
  {
    SLIC_ERROR("Mesh contains no topologies");
    return;
  }
  // Select the shared coordset from the first unstructured topology that can produce Gregory facets.
  const conduit::Node* first_surface_topology = nullptr;
  for(const conduit::Node& topology : topologies.children())
  {
    // Skip any not-unstructured topologies.
    if(topology.fetch_existing("type").as_string() != "unstructured")
    {
      continue;
    }

    // Skip any unstructured topologies with non tri, quad, or mixed elements
    const std::string shape = topology.fetch_existing("elements/shape").as_string();
    if(shape == "tri" || shape == "quad" || shape == "mixed")
    {
      first_surface_topology = &topology;
      break;
    }
  }

  if(first_surface_topology == nullptr)
  {
    SLIC_ERROR("Mesh contains no supported surface topologies.");
    return;
  }

  const std::string coordset_name = first_surface_topology->fetch_existing("coordset").as_string();
  const conduit::Node& coordset = domain.fetch_existing("coordsets/" + coordset_name);
  if(coordset.fetch_existing("type").as_string() != "explicit")
  {
    SLIC_ERROR("Coordinate set must be explicit");
    return;
  }

  // Access coordinates in place.
  const conduit::Node& values = coordset.fetch_existing("values");
  const auto x = values.fetch_existing("x").as_double_accessor();
  const auto y = values.fetch_existing("y").as_double_accessor();
  const auto z = values.fetch_existing("z").as_double_accessor();

  // Each topology becomes one named surface, even if its facets are disconnected.
  m_surfaces.reserve(topologies.number_of_children());
  for(const conduit::Node& topology : topologies.children())
  {
    // Skip any topology which doesn't use the first valid one seen
    if(topology.fetch_existing("coordset").as_string() != coordset_name)
    {
      SLIC_WARNING("Topology '" << topology.name() << "' uses a different coordset. Skipping.");
      continue;
    }

    // Derive facet IDs from Blueprint connectivity
    const auto facet_node_ids = readFacets(topology);
    if(facet_node_ids.empty())
    {
      // Skip topologies without facets
      continue;
    }

    // Build the geometry of each facet by averaging shared-vertex normals.
    m_surfaces.push_back(buildSurface(topology.name(), facet_node_ids, x, y, z));
  }
}

/*!
 * \brief Write the Gregory surfaces as triangles in a legacy ASCII VTK file.
 *
 * Each facet is sampled independently, so adjacent facets have separate output
 * point IDs. A triangle produces `subdivisions * subdivisions` VTK triangles;
 * a quadrilateral produces twice that number. The cell centered `surface_id`
 * scalar stores the index of the originating surface in getSurfaces(), while
 * `facet_id` stores a zero based facet index across surfaces in that order.
 *
 * Invalid subdivisions, an empty facet collection, and file I/O failures are
 * reported through SLIC_ERROR.
 *
 * \param[in] file_path Destination VTK file path.
 * \param[in] subdivisions Positive number of divisions along each parameter edge.
 */
void GregoryMesh::writeVtkFile(const std::string& file_path, int subdivisions) const
{
  if(subdivisions <= 0)
  {
    SLIC_ERROR("GregoryMesh VTK subdivisions must be positive");
    return;
  }

  std::vector<PointType> points;
  std::vector<std::array<axom::IndexType, 3>> triangles;
  std::vector<int> surface_ids;
  std::vector<int> facet_ids;
  const double step = 1.0 / subdivisions;

  // Evaluate a facet at a given point, append it to the list
  const auto add_point = [&](const GregoryFacet& facet, double u, double v) {
    points.push_back(facet.evaluate(u, v));
    return static_cast<axom::IndexType>(points.size() - 1);
  };

  // Push a triangle to the array by node IDs and original surface and facet IDs
  const auto add_triangle =
    [&](axom::IndexType a, axom::IndexType b, axom::IndexType c, int surface_id, int facet_id) {
      triangles.push_back({a, b, c});
      surface_ids.push_back(surface_id);
      facet_ids.push_back(facet_id);
    };

  // Iterate through each surface
  int next_facet_id = 0;
  for(axom::IndexType surface_index = 0; surface_index < m_surfaces.size(); ++surface_index)
  {
    // Iterate through each facet of the surface
    const int surface_id = static_cast<int>(surface_index);
    for(const GregoryFacet& facet : m_surfaces[surface_index].facets)
    {
      const int facet_id = next_facet_id++;
      if(facet.shape == GregoryFacetShape::Tri)
      {
        // Add the coordinates of each subdivision...
        std::vector<axom::IndexType> nodes;
        std::vector<std::size_t> row_starts;
        for(int i = 0; i <= subdivisions; ++i)
        {
          row_starts.push_back(nodes.size());
          for(int j = 0; j <= subdivisions - i; ++j)
          {
            nodes.push_back(add_point(facet, i * step, j * step));
          }
        }

        // ...and the indices which correspond to each coordinate
        const auto node = [&](int i, int j) { return nodes[row_starts[i] + j]; };
        for(int i = 0; i < subdivisions; ++i)
        {
          for(int j = 0; j < subdivisions - i; ++j)
          {
            // GregoryTriangle parameters put corner 1 on the v axis and corner 2 on u.
            add_triangle(node(i, j), node(i, j + 1), node(i + 1, j), surface_id, facet_id);
            if(i + j + 1 < subdivisions)
            {
              add_triangle(node(i + 1, j), node(i, j + 1), node(i + 1, j + 1), surface_id, facet_id);
            }
          }
        }
      }
      else
      {
        // Add the coordinates of each subdivision...
        std::vector<axom::IndexType> nodes;
        for(int j = 0; j <= subdivisions; ++j)
        {
          for(int i = 0; i <= subdivisions; ++i)
          {
            nodes.push_back(add_point(facet, i * step, j * step));
          }
        }

        // ...and the indices which correspond to each coordinate
        const auto node = [&](int i, int j) { return nodes[j * (subdivisions + 1) + i]; };
        for(int j = 0; j < subdivisions; ++j)
        {
          for(int i = 0; i < subdivisions; ++i)
          {
            add_triangle(node(i, j), node(i + 1, j), node(i + 1, j + 1), surface_id, facet_id);
            add_triangle(node(i, j), node(i + 1, j + 1), node(i, j + 1), surface_id, facet_id);
          }
        }
      }
    }
  }

  if(triangles.empty())
  {
    SLIC_ERROR("Cannot write a GregoryMesh VTK file with no facets");
    return;
  }

  std::ofstream output(file_path);
  if(!output.is_open())
  {
    SLIC_ERROR("Could not open GregoryMesh VTK file at path '" << file_path << "'");
    return;
  }

  // Dump the VTK data

  output << "# vtk DataFile Version 3.0\n"
         << "Axom Quest Gregory mesh visualization\n"
         << "ASCII\n"
         << "DATASET UNSTRUCTURED_GRID\n"
         << "POINTS " << points.size() << " double\n";
  output << std::setprecision(17);
  for(const PointType& point : points)
  {
    output << point[0] << ' ' << point[1] << ' ' << point[2] << '\n';
  }

  output << "CELLS " << triangles.size() << ' ' << 4 * triangles.size() << '\n';
  for(const auto& triangle : triangles)
  {
    output << "3 " << triangle[0] << ' ' << triangle[1] << ' ' << triangle[2] << '\n';
  }

  output << "CELL_TYPES " << triangles.size() << '\n';
  for(std::size_t cell_id = 0; cell_id < triangles.size(); ++cell_id)
  {
    output << "5\n";  // VTK_TRIANGLE
  }

  output << "CELL_DATA " << surface_ids.size() << '\n'
         << "SCALARS surface_id int 1\n"
         << "LOOKUP_TABLE default\n";
  for(const int surface_id : surface_ids)
  {
    output << surface_id << '\n';
  }
  output << "SCALARS facet_id int 1\n"
         << "LOOKUP_TABLE default\n";
  for(const int facet_id : facet_ids)
  {
    output << facet_id << '\n';
  }

  output.close();
  if(!output)
  {
    SLIC_ERROR("Could not write GregoryMesh VTK file at path '" << file_path << "'");
  }
}

}  // namespace axom::quest
