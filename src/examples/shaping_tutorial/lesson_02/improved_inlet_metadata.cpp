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
 *   dim: 2
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
 *
 * For 3D:
 * mesh:
 *   dim: 3
 *   bounding_box:
 *     min:
 *       x: 0.0
 *       y: 0.0
 *       z: 0.0
 *     max:
 *       x: 1.0
 *       y: 1.5
 *       z: 2.0
 *   resolution:
 *     x: 15
 *     y: 25
 *     z: 30
 */

namespace inlet = axom::inlet;

// Definition of the MeshMetadata struct
struct MeshMetadata
{
  struct BoundingBox
  {
    std::vector<double> min;
    std::vector<double> max;
  };

  int dim;
  BoundingBox bounding_box;
  std::vector<int> resolution;

  // Define schema for MeshMetadata with validation
  static void defineSchema(inlet::Container& mesh_schema)
  {
    // Add dimension (either 2 or 3)
    mesh_schema.addInt("dim", "Dimension (2 or 3)").required().range(2, 3);

    // setup bounding box info. min values must be less than max values.
    auto& bb = mesh_schema.addStruct("bounding_box", "Mesh bounding box").required();

    auto& min = bb.addStruct("min", "Minimum coordinates").required();
    min.addDouble("x", "Minimum x coordinate").required();
    min.addDouble("y", "Minimum y coordinate").required();
    min.addDouble("z", "Minimum z coordinate (only specify when dim is 3)");

    auto& max = bb.addStruct("max", "Maximum coordinates").required();
    max.addDouble("x", "Maximum x coordinate").required();
    max.addDouble("y", "Maximum y coordinate").required();
    max.addDouble("z", "Maximum z coordinate (only specify when dim is 3)");

    // each resolution value must be positive
    auto& res = mesh_schema.addStruct("resolution", "Mesh resolution").required();
    res.addInt("x", "Resolution in x direction").required().range(1, std::numeric_limits<int>::max());
    res.addInt("y", "Resolution in y direction").required().range(1, std::numeric_limits<int>::max());
    res.addInt("z", "Resolution in z direction (only specify when dim is 3)")
      .range(1, std::numeric_limits<int>::max());

    // Add constraints to ensure min < max for each coordinate
    bb.registerVerifier([](const inlet::Container& input) {
      bool valid = true;
      for(const std::string& axis : {"x", "y", "z"})
      {
        const std::string min_str = axom::fmt::format("min/{}", axis);
        const std::string max_str = axom::fmt::format("max/{}", axis);
        if(axis == "z" && (!input.contains(min_str) && !input.contains(max_str)))  // skip for 2d inputs
        {
          continue;
        }

        if(const double min_val = input[min_str], max_val = input[max_str]; min_val >= max_val)
        {
          SLIC_WARNING(axom::fmt::format("Invalid bounding box range for {}-coordinate: {} >= {}",
                                         axis,
                                         min_val,
                                         max_val));
          valid = false;
        }
      }
      return valid;
    });

    // Add constraint to ensure z values are only provided when dim is 3
    mesh_schema.registerVerifier([](const inlet::Container& input) {
      const int dim = input["dim"];
      bool valid = true;

      for(const auto& field : {"bounding_box/min/z", "bounding_box/max/z", "resolution/z"})
      {
        if(dim == 3)
        {
          if(!input.contains(field))
          {
            SLIC_WARNING(
              axom::fmt::format("Z-coordinate for '{}' is required when dimension is 3", field));
            valid = false;
          }
        }
        else if(dim == 2)
        {
          if(input.contains(field))
          {
            SLIC_WARNING(
              axom::fmt::format("Z-coordinate for '{}' should not be provided when dimension is 2",
                                field));
            valid = false;
          }
        }
      }

      return valid;
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

    result.dim = input_data["dim"];

    // Initialize vectors with appropriate size based on dimension
    result.bounding_box.min.resize(result.dim);
    result.bounding_box.max.resize(result.dim);
    result.resolution.resize(result.dim);

    auto bb = input_data["bounding_box"];
    result.bounding_box.min[0] = bb["min/x"];
    result.bounding_box.min[1] = bb["min/y"];

    result.bounding_box.max[0] = bb["max/x"];
    result.bounding_box.max[1] = bb["max/y"];

    auto res = input_data["resolution"];
    result.resolution[0] = res["x"];
    result.resolution[1] = res["y"];

    // Only grab z values when dimension is 3
    if(result.dim == 3)
    {
      result.bounding_box.min[2] = bb["min/z"];
      result.bounding_box.max[2] = bb["max/z"];
      result.resolution[2] = res["z"];
    }

    return result;
  }
};

void print_metadata(const MeshMetadata& metadata)
{
  fmt::print("Dimension: {}\n", metadata.dim);

  if(metadata.dim == 2)
  {
    fmt::print("Bounding Box Min: ({}, {})\n",
               metadata.bounding_box.min[0],
               metadata.bounding_box.min[1]);
    fmt::print("Bounding Box Max: ({}, {})\n",
               metadata.bounding_box.max[0],
               metadata.bounding_box.max[1]);
    fmt::print("Resolution: ({}, {})\n", metadata.resolution[0], metadata.resolution[1]);
  }
  else
  {
    fmt::print("Bounding Box Min: ({}, {}, {})\n",
               metadata.bounding_box.min[0],
               metadata.bounding_box.min[1],
               metadata.bounding_box.min[2]);
    fmt::print("Bounding Box Max: ({}, {}, {})\n",
               metadata.bounding_box.max[0],
               metadata.bounding_box.max[1],
               metadata.bounding_box.max[2]);
    fmt::print("Resolution: ({}, {}, {})\n",
               metadata.resolution[0],
               metadata.resolution[1],
               metadata.resolution[2]);
  }
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
