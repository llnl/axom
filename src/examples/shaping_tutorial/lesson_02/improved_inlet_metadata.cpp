// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
/*! 
 * \file improved_inlet_metadata.cpp
 * \brief Example of how to use Inlet to parse and validate inlet metadata from YAML or Lua input
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
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
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
  int dim;
  axom::Array<double> bb_min;
  axom::Array<double> bb_max;
  axom::Array<int> resolution;

  // Convert MeshMetadata to axom::primal::BoundingBox
  template <int DIM>
  axom::primal::BoundingBox<double, DIM> getBoundingBox() const
  {
    static_assert(DIM == 2 || DIM == 3, "Invalid dimension");
    SLIC_ASSERT_MSG(DIM == dim, "Template dimension must match MeshMetadata dimension");

    if constexpr(DIM == 2)
    {
      return axom::primal::BoundingBox<double, DIM>({bb_min[0], bb_min[1]}, {bb_max[0], bb_max[1]});
    }
    else  // DIM == 3
    {
      return axom::primal::BoundingBox<double, DIM>({bb_min[0], bb_min[1], bb_min[2]},
                                                    {bb_max[0], bb_max[1], bb_max[2]});
    }
  }

  // Define schema for MeshMetadata with validation
  static void defineSchema(inlet::Container& mesh_schema)
  {
    // Add dimension (either 2 or 3)
    mesh_schema.addInt("dim", "Dimension (2 or 3)").required().range(2, 3);

    // set up bounding box info. min values must be less than max values.
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
      for(const std::string axis : {"x", "y", "z"})
      {
        const std::string min_str = axom::fmt::format("min/{}", axis);
        const std::string max_str = axom::fmt::format("max/{}", axis);
        if(axis == "z" && (!input.contains(min_str) && !input.contains(max_str)))  // skip z for 2d inputs
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
    result.bb_min.resize(result.dim);
    result.bb_max.resize(result.dim);
    result.resolution.resize(result.dim);

    auto bb = input_data["bounding_box"];
    result.bb_min[0] = bb["min/x"];
    result.bb_min[1] = bb["min/y"];

    result.bb_max[0] = bb["max/x"];
    result.bb_max[1] = bb["max/y"];

    auto res = input_data["resolution"];
    result.resolution[0] = res["x"];
    result.resolution[1] = res["y"];

    // Only grab z values when dimension is 3
    if(result.dim == 3)
    {
      result.bb_min[2] = bb["min/z"];
      result.bb_max[2] = bb["max/z"];
      result.resolution[2] = res["z"];
    }

    return result;
  }
};

void print_metadata(const MeshMetadata& metadata)
{
  SLIC_INFO(
    axom::fmt::format("Parsed mesh metadata:"
                      "\n  - dimension: {}"
                      "\n  - bounding Box: {}"
                      "\n  - resolution: {}",
                      metadata.dim,
                      metadata.dim == 2 ? axom::fmt::format("{}", metadata.getBoundingBox<2>())
                                        : axom::fmt::format("{}", metadata.getBoundingBox<3>()),
                      metadata.resolution));
}

int main(int argc, char** argv)
{
  // define the list of supported extensions
  const std::vector<std::string> supported_extensions = []() {
    std::vector<std::string> vec = {".yaml", ".yml"};
#ifdef AXOM_USE_LUA
    vec.push_back(".lua");
#endif
    return vec;
  }();

  axom::slic::SimpleLogger logger;

  // parse and validate input
  axom::CLI::App app {"Inlet Metadata Setup"};
  std::string inputFilename;
  app.add_option("input_file", inputFilename)
    ->description(axom::fmt::format("input file containing inlet metadata, supported formats: {}",
                                    supported_extensions))
    ->required()
    ->check(axom::CLI::ExistingFile)
    ->check([&supported_extensions](const std::string& filename) {
      for(auto& ext : supported_extensions)
      {
        if(axom::utilities::string::endsWith(filename, ext))
        {
          return std::string();
        }
      }
      return axom::fmt::format("Invalid extension for file '{}'; supported extensions: {}",
                               filename,
                               supported_extensions);
    });

  CLI11_PARSE(app, argc, argv);

  // Create appropriate reader based on file extension
  std::unique_ptr<axom::inlet::Reader> reader;

  if(axom::utilities::string::endsWith(inputFilename, ".yaml") ||
     axom::utilities::string::endsWith(inputFilename, ".yml"))
  {
    reader = std::make_unique<axom::inlet::YAMLReader>();
  }
#ifdef AXOM_USE_LUA
  else if(axom::utilities::string::endsWith(inputFilename, ".lua"))
  {
    reader = std::make_unique<axom::inlet::LuaReader>();
  }
#endif

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
