// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
/*! 
 * \file klee_operators_and_validation.cpp
 * \brief Loads and validates a Klee input file
 */
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/klee.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <limits>

/**
 * Get the list of materials that the given shape will replace out of a given list of materials.
 *
 * @param shape The shape to check against.
 * @param materials The list of materials to check.
 * @return A vector containing the materials from the input list that the shape will replace.
 */
inline std::vector<std::string> getMaterialsReplacedBy(const axom::klee::Shape& shape,
                                                       const std::vector<std::string>& materials)
{
  std::vector<std::string> replacedMaterials;
  for(const auto& material : materials)
  {
    if(shape.replaces(material))
    {
      replacedMaterials.push_back(material);
    }
  }
  return replacedMaterials;
}

// Function to print information about klee ShapeSet
void printShapeSetInfo(const axom::klee::ShapeSet& shapeSet)
{
  // lambda to help format a klee::Dimensions
  auto dimensionsToString = [](axom::klee::Dimensions dim) -> std::string {
    return dim == axom::klee::Dimensions::Two ? "2D"
      : dim == axom::klee::Dimensions::Three  ? "3D"
                                              : "<unknown>";
  };

  // lambda to help format a klee::LengthUnit
  auto lengthUnitToString = [](axom::klee::LengthUnit unit) -> std::string {
    switch(unit)
    {
    case axom::klee::LengthUnit::km:
      return "km";
    case axom::klee::LengthUnit::m:
      return "m";
    case axom::klee::LengthUnit::dm:
      return "dm";
    case axom::klee::LengthUnit::cm:
      return "cm";
    case axom::klee::LengthUnit::mm:
      return "mm";
    case axom::klee::LengthUnit::um:
      return "um";
    case axom::klee::LengthUnit::nm:
      return "nm";
    case axom::klee::LengthUnit::angstrom:
      return "A";
    case axom::klee::LengthUnit::miles:
      return "miles";
    case axom::klee::LengthUnit::feet:
      return "feet";
    case axom::klee::LengthUnit::inches:
      return "inches";
    case axom::klee::LengthUnit::mils:
      return "mils";
    case axom::klee::LengthUnit::unspecified:
    default:
      return "<unspecified>";
    }
  };

  // lambda to help format a parir of klee::LengthUnits
  auto formattedUnits = [&lengthUnitToString](axom::klee::LengthUnit startUnits,
                                              axom::klee::LengthUnit endUnits) -> std::string {
    return (startUnits == endUnits) ? axom::fmt::format("'{}'", lengthUnitToString(startUnits))
                                    : axom::fmt::format("'{}' -> '{}'",
                                                        lengthUnitToString(startUnits),
                                                        lengthUnitToString(endUnits));
  };

  // lambda to help format a parir of klee::Dimensions
  auto formattedDimensions = [&dimensionsToString](axom::klee::Dimensions startDims,
                                                   axom::klee::Dimensions endDims) -> std::string {
    return (startDims == endDims)
      ? dimensionsToString(startDims)
      : axom::fmt::format("{} -> {}", dimensionsToString(startDims), dimensionsToString(endDims));
  };

  axom::fmt::memory_buffer buffer;
  axom::fmt::format_to(std::back_inserter(buffer), "Klee ShapeSet Information:\n");

  axom::fmt::format_to(std::back_inserter(buffer),
                       "  Overall dimensions: {}\n",
                       dimensionsToString(shapeSet.getDimensions()));

  // Collect and print sorted list of material names
  std::vector<std::string> materials = [&shapeSet]() {
    std::set<std::string> materialSet;
    for(const auto& shape : shapeSet.getShapes())
    {
      materialSet.insert(shape.getMaterial());
    }
    return std::vector<std::string>(materialSet.begin(), materialSet.end());
  }();
  axom::fmt::format_to(std::back_inserter(buffer),
                       "\n  Unique materials ({}): {}\n",
                       materials.size(),
                       materials);

  // Print information about each shape
  axom::fmt::format_to(std::back_inserter(buffer),
                       "\n  Details for the {} shapes:\n",
                       shapeSet.getShapes().size());
  for(const auto& shape : shapeSet.getShapes())
  {
    const auto& geom = shape.getGeometry();

    axom::fmt::format_to(
      std::back_inserter(buffer),
      "  - name: '{}'\n"
      "    material: '{}'\n"
      "    format: '{}'\n"
      "    units: {}\n"
      "    dimensions: {}\n"
      "    replaces materials: {}\n\n",
      shape.getName(),
      shape.getMaterial(),
      geom.getFormat(),
      formattedUnits(geom.getStartProperties().units, geom.getEndProperties().units),
      formattedDimensions(geom.getInputDimensions(), geom.getOutputDimensions()),
      getMaterialsReplacedBy(shape, materials));
  }

  SLIC_INFO(axom::fmt::to_string(buffer));
  axom::slic::flushStreams();
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;

  // CLI
  axom::CLI::App app {"Klee Input Validator and Summary"};
  std::string inputFilename;
  app.add_option("input", inputFilename)
    ->description("Klee input file")
    ->required()
    ->check(axom::CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  // Load the klee shape file and extract some information
  try
  {
    auto shapeSet = axom::klee::readShapeSet(inputFilename);
  }
  catch(axom::klee::KleeError& error)
  {
    std::vector<std::string> errs;
    for(auto verificationError : error.getErrors())
    {
      errs.push_back(axom::fmt::format(" - '{}': {}",
                                       static_cast<std::string>(verificationError.path),
                                       verificationError.message));
    }

    SLIC_WARNING(
      axom::fmt::format("Error during parsing klee input. Found the following errors:\n{}",
                        axom::fmt::join(errs, "\n")));
    exit(1);
  }

  auto shapeSet = axom::klee::readShapeSet(inputFilename);
  printShapeSetInfo(shapeSet);

  return 0;
}