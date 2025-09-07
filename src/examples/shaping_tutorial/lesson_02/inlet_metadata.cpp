#include "axom/config.hpp"
#include "axom/inlet.hpp"
#include "axom/sidre.hpp"
#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <iostream>
#include <fstream>
#include <vector>

using axom::inlet::Inlet;
using axom::inlet::YAMLReader;

/*! 
 * \file inlet_metadata.cpp
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

  // Define schema for MeshMetadata
  static void defineSchema(inlet::Container& mesh_schema)
  {
    auto& bb = mesh_schema.addStruct("bounding_box", "Mesh bounding box");
    auto& min = bb.addStruct("min", "Minimum coordinates");
    min.addDouble("x", "Minimum x coordinate").defaultValue(0.0);
    min.addDouble("y", "Minimum y coordinate").defaultValue(0.0);

    auto& max = bb.addStruct("max", "Maximum coordinates");
    max.addDouble("x", "Maximum x coordinate").defaultValue(1.0);
    max.addDouble("y", "Maximum y coordinate").defaultValue(1.0);

    auto& res = mesh_schema.addStruct("resolution", "Mesh resolution");
    res.addInt("x", "Resolution in x direction").defaultValue(10);
    res.addInt("y", "Resolution in y direction").defaultValue(20);
  }
};

// Specialization of FromInlet for MeshMetadata
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
  fmt::print("Bounding Box Min: ({}, {})\n", metadata.bounding_box.min_x, metadata.bounding_box.min_y);
  fmt::print("Bounding Box Max: ({}, {})\n", metadata.bounding_box.max_x, metadata.bounding_box.max_y);
  fmt::print("Resolution: ({}, {})\n", metadata.resolution.x, metadata.resolution.y);
}

int main(int argc, char** argv)
{
  axom::CLI::App app {"Inlet Metadata Setup"};

  std::string inputFilename;
  app.add_option("input_file", inputFilename, "YAML input file with inlet metadata")->required();
  CLI11_PARSE(app, argc, argv);

  // Parse YAML directly to MeshMetadata
  auto reader = std::make_unique<YAMLReader>();
  reader->parseFile(inputFilename);
  Inlet inlet(std::move(reader));

  auto& mesh_schema = inlet.addStruct("mesh", "Mesh metadata");
  MeshMetadata::defineSchema(mesh_schema);

  if(!inlet.verify())
  {
    SLIC_ERROR("badness");
  }

  MeshMetadata metadata = inlet["mesh"].get<MeshMetadata>();
  print_metadata(metadata);

  return 0;
}
