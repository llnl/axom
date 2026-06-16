// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"
#include "axom/fmt.hpp"

#include <string>
#include <type_traits>
#include <variant>

namespace inlet = axom::inlet;

// _inlet_user_defined_variant_start
struct Circle
{
  double radius;
};

struct Box
{
  double width;
  double height;
};

using Shape = std::variant<Circle, Box>;

template <>
struct FromInlet<Circle>
{
  Circle operator()(const inlet::Container& input_data) { return {input_data["radius"]}; }
};

template <>
struct FromInlet<Box>
{
  Box operator()(const inlet::Container& input_data)
  {
    return {input_data["width"], input_data["height"]};
  }
};

template <>
struct FromInlet<Shape>
{
  // Shape is a std::variant, so it cannot be read through Container::get<T>().
  // Use the public Proxy returned by inlet["shape"] to access its fields instead.
  Shape operator()(const inlet::Proxy& input_data)
  {
    const std::string kind = input_data["kind"];
    if(kind == "circle")
    {
      return Circle {input_data["radius"]};
    }

    return Box {input_data["width"], input_data["height"]};
  }
};

void defineShapeSchema(inlet::Container& shape)
{
  shape.addString("kind", "Shape variant discriminator").required().validValues({"circle", "box"});
  shape.addDouble("radius", "Circle radius").required(false);
  shape.addDouble("width", "Box width").required(false);
  shape.addDouble("height", "Box height").required(false);

  shape.registerVerifier([](const inlet::Container& input_data) {
    if(!input_data.isUserProvided("kind"))
    {
      return false;
    }

    const std::string kind = input_data["kind"];
    if(kind == "circle")
    {
      return input_data.isUserProvided("radius") && !input_data.isUserProvided("width") &&
        !input_data.isUserProvided("height");
    }

    return input_data.isUserProvided("width") && input_data.isUserProvided("height") &&
      !input_data.isUserProvided("radius");
  });
}
// _inlet_user_defined_variant_end

const std::string input = R"(
  -- _inlet_user_defined_variant_input_start
  shape = {
    kind = "box",
    width = 3.0,
    height = 4.0
  }
  -- _inlet_user_defined_variant_input_end
)";

void printShape(const Shape& shape)
{
  std::visit(
    [](const auto& concrete_shape) {
      using ShapeType = std::decay_t<decltype(concrete_shape)>;
      if constexpr(std::is_same_v<ShapeType, Circle>)
      {
        SLIC_INFO(axom::fmt::format("circle radius = {}", concrete_shape.radius));
      }
      else
      {
        SLIC_INFO(axom::fmt::format("box width = {}, height = {}",
                                    concrete_shape.width,
                                    concrete_shape.height));
      }
    },
    shape);
}

int main()
{
  // Initialize Axom's logger
  axom::slic::SimpleLogger logger;

  // Create Inlet object with the Lua Reader and parse the input file snippet
  auto lr = std::make_unique<inlet::LuaReader>();
  lr->parseString(input);
  inlet::Inlet inlet(std::move(lr));

  // Define the input file schema
  // _inlet_user_defined_variant_schema_usage_start
  auto& shape_schema = inlet.addStruct("shape", "A single user-defined variant");
  defineShapeSchema(shape_schema);
  // _inlet_user_defined_variant_schema_usage_end

  // Verify input file validates against the schema
  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Create a shape object from inlet container
  // _inlet_user_defined_variant_access_start
  const Shape shape = FromInlet<Shape> {}(inlet["shape"]);
  // _inlet_user_defined_variant_access_end

  printShape(shape);

  return 0;
}
