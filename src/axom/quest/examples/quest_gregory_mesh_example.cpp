// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_gregory_mesh_example.cpp
 * \brief Load a Blueprint surface mesh, summarize its Gregory facets, and optionally write VTK.
 */

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "quest_gregory_mesh_example.cpp requires Conduit"
#endif

#include "axom/CLI11.hpp"
#include "axom/quest/GregoryMesh.hpp"
#include "axom/slic.hpp"

#include "conduit_node.hpp"
#include "conduit_relay_io.hpp"
#include "conduit_relay_io_blueprint.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <string>

//------------------------------------------------------------------------------
// CLI input
//------------------------------------------------------------------------------

/// Command line options for the GregoryMesh example.
class Input
{
public:
  std::string inputFile;
  std::string outputVtk;
  int subdivisions {8};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-i,--input", inputFile, "Path to a Blueprint mesh YAML or root file")
      ->required()
      ->check(axom::CLI::ExistingFile);
    app.add_option("-o,--output-vtk",
                   outputVtk,
                   "Write a triangulated Gregory mesh to this VTK file");
    app.add_option("--subdivisions", subdivisions, "Subdivisions per facet edge for VTK output")
      ->check(axom::CLI::PositiveNumber);

    app.parse(argc, argv);
  }
};

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;

  Input input;
  axom::CLI::App app {"Construct Gregory surfaces from a Conduit Blueprint mesh"};
  try
  {
    input.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& error)
  {
    return app.exit(error);
  }

  // Read in the mesh
  conduit::Node blueprint_mesh;
  conduit::relay::io::blueprint::load_mesh(input.inputFile, blueprint_mesh);

  // Construct the Gregory mesh
  const axom::quest::GregoryMesh gregory_mesh(blueprint_mesh);
  const auto& surfaces = gregory_mesh.getSurfaces();
  std::size_t total_triangles = 0;
  std::size_t total_quadrilaterals = 0;
  axom::quest::GregoryFacet::BoundingBoxType mesh_bounds;

  // For each surface, record the number and type of each facet
  std::cout << "Gregory surfaces: " << surfaces.size() << '\n';
  for(const auto& surface : surfaces)
  {
    std::size_t triangles = 0;
    std::size_t quadrilaterals = 0;
    axom::quest::GregoryFacet::BoundingBoxType surface_bounds;
    for(const auto& facet : surface.facets)
    {
      surface_bounds.addBox(facet.boundingBox());
      if(facet.shape == axom::quest::GregoryFacetShape::Tri)
      {
        ++triangles;
      }
      else
      {
        ++quadrilaterals;
      }
    }

    total_triangles += triangles;
    total_quadrilaterals += quadrilaterals;
    mesh_bounds.addBox(surface_bounds);

    // Report per-surface statistics
    std::cout << "  " << surface.name << ": " << triangles << " triangles, " << quadrilaterals
              << " quadrilaterals, bounds " << surface_bounds.getMin() << " to "
              << surface_bounds.getMax() << '\n';
  }

  // Report the total statistics
  std::cout << "Total facets: " << total_triangles + total_quadrilaterals << " (" << total_triangles
            << " triangles, " << total_quadrilaterals << " quadrilaterals)\n";
  if(mesh_bounds.isValid())
  {
    std::cout << "Mesh bounds: " << mesh_bounds.getMin() << " to " << mesh_bounds.getMax() << '\n';
  }

  // Optionally write to a file
  if(!input.outputVtk.empty())
  {
    gregory_mesh.writeVtkFile(input.outputVtk, input.subdivisions);
    std::cout << "Wrote triangulated Gregory mesh to " << input.outputVtk << '\n';
  }
  return 0;
}
