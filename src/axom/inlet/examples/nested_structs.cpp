// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "axom/slic/core/SimpleLogger.hpp"
#include "axom/primal.hpp"

namespace inlet = axom::inlet;
using Vector = inlet::FunctionType::Vector;

// Used to convert a std::vector (the type Inlet uses to represent arrays)
// to a geometric vector type
Vector toVector(const std::vector<double>& vec)
{
  return Vector {vec.data(), static_cast<int>(vec.size())};
}

// A union of the members required for each of the operations is stored for simplicity
struct Operator
{
  double rotate;
  Vector axis;
  Vector center;
  Vector translate;
  double x;
  double y;
  double z;
  Vector origin;
  enum class Type
  {
    Translate,
    Rotate,
    Slice
  } type;

  // Again, the union of the necessary members are defined as part of the schema
  static void defineSchema(inlet::Container& container)
  {
    // The rotation is in degrees
    container.addDouble("rotate", "Degrees of rotation").range(-180, 180);
    // Vectors are defined as arrays of doubles
    container.addDoubleArray("axis", "Axis on which to rotate");
    container.addDoubleArray("center", "Center of rotation");
    container.addDoubleArray("translate", "Translation vector");
    // The slice operation can have sub-entries, so we represent it as a struct
    // Note that Inlet does not require a 1-1 correspondence (or any correspondence)
    // between structures defined in the schema via addStruct and structures extracted
    // from the input file with FromInlet specializations (defined below)
    auto& slice = container.addStruct("slice", "Options for a slice operation");
    slice.addDouble("x", "x-axis point to slice on");
    slice.addDouble("y", "y-axis point to slice on");
    slice.addDouble("z", "z-axis point to slice on");
    slice.addDoubleArray("origin", "Origin for the slice operation");

    // Verify that exactly one type of operator is defined
    container.registerVerifier([](const inlet::Container& container) {
      const bool is_translate = container.contains("translate");
      const bool is_rotate = container.contains("rotate");
      const bool is_slice = container.contains("slice");

      // There can be only one
      if((is_translate && is_rotate) || (is_translate && is_slice) || (is_rotate && is_slice))
      {
        return false;
      }

      return is_translate || is_rotate || is_slice;
    });
  }
};

template <>
struct FromInlet<Operator>
{
  Operator operator()(const inlet::Container& base)
  {
    Operator result;
    // Even though all the possible members are part of the schema, the
    // "contains" operation can be used to check what is actually present
    // in the input file
    if(base.contains("translate"))
    {
      result.type = Operator::Type::Translate;
      result.translate = toVector(base["translate"]);
    }
    else if(base.contains("rotate"))
    {
      result.type = Operator::Type::Rotate;
      result.rotate = base["rotate"];
      result.axis = toVector(base["axis"]);
      result.center = toVector(base["center"]);
    }
    else if(base.contains("slice"))
    {
      result.type = Operator::Type::Slice;
      // Grab the substructure corresponding to the slice operation
      // and use it to populate the Operator instance
      auto slice = base["slice"];
      if(slice.contains("x"))
      {
        result.x = slice["x"];
      }
      if(slice.contains("y"))
      {
        result.y = slice["y"];
      }
      if(slice.contains("z"))
      {
        result.z = slice["z"];
      }
      if(slice.contains("origin"))
      {
        result.origin = toVector(slice["origin"]);
      }
    }
    return result;
  }
};

std::ostream& operator<<(std::ostream& os, const Operator& op)
{
  switch(op.type)
  {
  case Operator::Type::Translate:
    os << "   Translation operator:\n";
    os << axom::fmt::format("      with vector: {0}\n", op.translate);
    break;
  case Operator::Type::Rotate:
    os << "   Rotation operator:\n";
    os << axom::fmt::format("      with axis: {0}\n", op.axis);
    os << axom::fmt::format("      with center: {0}\n", op.center);
    break;
  case Operator::Type::Slice:
    os << "   Slice operator:\n";
    os << axom::fmt::format("      with x-coord: {0}\n", op.x);
    os << axom::fmt::format("      with y-coord: {0}\n", op.y);
    os << axom::fmt::format("      with z-coord: {0}\n", op.z);
    os << axom::fmt::format("       with origin: {0}\n", op.origin);
    break;
  default:
    SLIC_ERROR("Operator had unknown type");
  }
  return os;
}

struct Geometry
{
  std::vector<Operator> operators;
  std::string format;
  std::string path;
  enum class Units
  {
    Centimeters,
    Meters
  } units;
  int start_dim;
  static void defineSchema(inlet::Container& container)
  {
    container.addString("format", "File format for the shape").required();
    container.addString("path", "Path to the shape file").required();
    // A string is used to represent the enumeration
    container.addString("units", "Units for length").defaultValue("cm").validValues({"cm", "m"});
    container.addInt("start_dimensions", "Dimension in which to begin applying operations")
      .defaultValue(3);
    // addStructArray is used to represent std::vector<T> where
    // T is any non-primitive type
    auto& ops_schema =
      container.addStructArray("operators", "List of shape operations to apply").required();
    Operator::defineSchema(ops_schema);
  }
};

template <>
struct FromInlet<Geometry>
{
  Geometry operator()(const inlet::Container& base)
  {
    Geometry result;
    result.format = base["format"];
    result.path = base["path"];
    result.start_dim = base["start_dimensions"];

    std::string units = base["units"];
    if(units == "cm")
    {
      result.units = Geometry::Units::Centimeters;
    }
    else if(units == "m")
    {
      result.units = Geometry::Units::Meters;
    }

    result.operators = base["operators"].get<std::vector<Operator>>();
    return result;
  }
};

std::ostream& operator<<(std::ostream& os, const Geometry& geom)
{
  os << axom::fmt::format("Geometry in format: '{0}'\n", geom.format);
  os << axom::fmt::format("  with path: '{0}'\n", geom.path);
  for(const auto& op : geom.operators)
  {
    os << op;
  }
  return os;
}

struct Shape
{
  std::string name;
  enum class Material
  {
    Steel,
    Wood,
    Plastic
  } material;
  Geometry geom;
  static void defineSchema(inlet::Container& container)
  {
    container.addString("name", "Name of the shape").required();
    container.addString("material", "Material of the shape")
      .validValues({"steel", "wood", "plastic"})
      .required();
    // addStruct is used for a single instance of a user-defined type
    auto& geom_schema =
      container.addStruct("geometry", "Geometric information on the shape").required();
    Geometry::defineSchema(geom_schema);
  }
};

std::ostream& operator<<(std::ostream& os, const Shape& shape)
{
  os << axom::fmt::format("Shape: '{0}'\n", shape.name);
  os << shape.geom;
  return os;
}

template <>
struct FromInlet<Shape>
{
  Shape operator()(const inlet::Container& base)
  {
    Shape result;
    result.name = base["name"];
    std::string material = base["material"];
    if(material == "steel")
    {
      result.material = Shape::Material::Steel;
    }
    else if(material == "wood")
    {
      result.material = Shape::Material::Wood;
    }
    else if(material == "plastic")
    {
      result.material = Shape::Material::Plastic;
    }
    result.geom = base["geometry"].get<Geometry>();
    return result;
  }
};

// "rubber" is not one of the allowed shape materials, so this will fail verification
const std::string input = R"(
dimensions: 2
shapes:
  - name: wheel
    material: steel
    geometry:
      format: test_format
      path: path/to/file.format
      units: cm
      start_dimensions: 3
      operators:
        - translate: [10, 20, 30]
        - rotate: 45
          axis: [1, 2, 3]
          center: [4, 5, 6]
        - slice:
            x: 10
            origin: [40, 50, 60]
        - translate: [5, 6]
  - name: tire
    material: rubber
    geometry:
      format: test_format
      path: path/to/tire.format
      units: cm
      operators:
        - translate: [10, 20]
)";

int main(int argc, char** argv)
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  axom::CLI::App app {"Example of Axom's Inlet component for nested structures"};
  bool docsEnabled {false};
  app.add_flag("--docs", docsEnabled, "Enables documentation generation");
  bool strictVerification {false};
  app.add_flag("--strict", strictVerification, "Warns if any unexpected fields are provided");
  CLI11_PARSE(app, argc, argv);

  auto reader = std::unique_ptr<inlet::YAMLReader>(new inlet::YAMLReader());
  reader->parseString(input);
  inlet::Inlet inlet(std::move(reader));

  // _inlet_nested_struct_array_start
  auto& shapes_container = inlet.addStructArray("shapes");
  Shape::defineSchema(shapes_container);
  // _inlet_nested_struct_array_end

  if(strictVerification)
  {
    // Mark the entire input as "strict"
    inlet.getGlobalContainer().strict();
  }

  if(inlet.verify())
  {
    SLIC_INFO("Verification was successful");
  }
  else
  {
    SLIC_INFO("Verification was unsuccessful");
  }

  auto shapes = inlet["shapes"].get<std::unordered_map<int, Shape>>();
  for(const auto& entry : shapes)
  {
    std::cout << entry.second << "\n";
  }

  if(docsEnabled)
  {
    const std::string docFileName = "nested_structs";
    inlet.write(inlet::SphinxWriter(docFileName + ".rst"));
    inlet.write(inlet::JSONSchemaWriter(docFileName + ".json"));
    SLIC_INFO("Documentation was written to " << docFileName << " (rst and json)");
  }
}
