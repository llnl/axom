// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
/*! 
 * \file validated_inlet_metadata.cpp
 * \brief Example of how to use Inlet to parse and validate inlet metadata from YAML input
 * using a user-defined MeshMetadata struct with a templated FromInlet specialization.
 *
 * Example run:
 * ./inlet_metadata input.yaml
 *
 * Where input.yaml contains:
 * mesh:
 *   bounding_box:
 *     min:
 *       x: 0.0
 *       y: 0.0
 *     max:
 *       x: 1.0
 *       y: 1.5
 *   resolution:
 *     x: 15
 *     y: 25
 */
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/inlet.hpp"
#include "axom/sidre.hpp"
#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using axom::inlet::Inlet;
using axom::inlet::YAMLReader;

namespace inlet = axom::inlet;

// Definition of the MeshMetadata struct
struct MeshMetadata
{
  struct BoundingBox
  {
    double min_x, min_y;
    double max_x, max_y;
  };

  struct Resolution
  {
    int x, y;
  };

  BoundingBox bounding_box;
  Resolution resolution;

  // Define schema for MeshMetadata with validation
  static void defineSchema(inlet::Container& mesh_schema)
  {
    // setup bounding box info. min values must be less than max values.
    auto& bb = mesh_schema.addStruct("bounding_box", "Mesh bounding box").required();

    auto& min = bb.addStruct("min", "Minimum coordinates").required();
    min.addDouble("x", "Minimum x coordinate").required();
    min.addDouble("y", "Minimum y coordinate").required();

    auto& max = bb.addStruct("max", "Maximum coordinates").required();
    max.addDouble("x", "Maximum x coordinate").required();
    max.addDouble("y", "Maximum y coordinate").required();

    // each resolution value must be positive
    auto& res = mesh_schema.addStruct("resolution", "Mesh resolution").required();
    res.addInt("x", "Resolution in x direction").required().range(1, std::numeric_limits<int>::max());
    res.addInt("y", "Resolution in y direction").required().range(1, std::numeric_limits<int>::max());

    // Add constraints to ensure min < max for each coordinate
    bb.registerVerifier([](const inlet::Container& input) {
      const double min_x = input["min/x"];
      const double max_x = input["max/x"];
      const double min_y = input["min/y"];
      const double max_y = input["max/y"];

      SLIC_WARNING_IF(
        min_x >= max_x,
        axom::fmt::format("Invalid bounding box range for x-coordinate: {} >= {}", min_x, max_x));
      SLIC_WARNING_IF(
        min_y >= max_y,
        axom::fmt::format("Invalid bounding box range for y-coordinate: {} >= {}", min_y, max_y));
      return (min_x < max_x) && (min_y < max_y);
    });
  }
};

// Initialize a MeshMetadata from inlet
template <>
struct FromInlet<MeshMetadata>
{
  MeshMetadata operator()(const inlet::Container& input_data)
  {
    MeshMetadata result;

    auto bb = input_data["bounding_box"];
    result.bounding_box.min_x = bb["min/x"];
    result.bounding_box.min_y = bb["min/y"];
    result.bounding_box.max_x = bb["max/x"];
    result.bounding_box.max_y = bb["max/y"];

    auto res = input_data["resolution"];
    result.resolution.x = res["x"];
    result.resolution.y = res["y"];

    return result;
  }
};

void print_metadata(const MeshMetadata& metadata)
{
  axom::fmt::print("Bounding Box Min: ({}, {})\n",
                   metadata.bounding_box.min_x,
                   metadata.bounding_box.min_y);
  axom::fmt::print("Bounding Box Max: ({}, {})\n",
                   metadata.bounding_box.max_x,
                   metadata.bounding_box.max_y);
  axom::fmt::print("Resolution: ({}, {})\n", metadata.resolution.x, metadata.resolution.y);
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;
  axom::CLI::App app {"Inlet Metadata Setup"};

  std::string inputFilename;
  app.add_option("input_file", inputFilename, "YAML input file with inlet metadata")->required();
  CLI11_PARSE(app, argc, argv);

  // Parse YAML directly to MeshMetadata
  auto reader = std::make_unique<YAMLReader>();
  reader->parseFile(inputFilename);
  Inlet inlet(std::move(reader));

  // Define schema at top level
  auto& mesh_schema = inlet.addStruct("mesh", "Mesh metadata").required();
  MeshMetadata::defineSchema(mesh_schema);

  // Validate the input
  if(!inlet.verify())
  {
    SLIC_WARNING("Error: Input validation failed.");
    SLIC_WARNING("Missing required fields or invalid data.");
    return 1;
  }

  // Initialize a MeshMetadata from inlet and print its values
  MeshMetadata metadata = inlet["mesh"].get<MeshMetadata>();
  print_metadata(metadata);

  return 0;
}
