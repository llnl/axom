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
#include <vector>

namespace inlet = axom::inlet;

// _inlet_variant_struct_collections_start
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
  Circle operator()(const inlet::Container& input_data)
  {
    return {input_data["radius"]};
  }
};

template <>
struct FromInlet<Box>
{
  Box operator()(const inlet::Container& input_data)
  {
    return {input_data["width"], input_data["height"]};
  }
};

void defineShapeSchema(inlet::VariantStructCollection<Shape>& shapes)
{
  shapes.addAlternative<Circle>("circle", [](inlet::Container& circle) {
    circle.addDouble("radius").required();
  });
  shapes.addAlternative<Box>("box", [](inlet::Container& box) {
    box.addDouble("width").required();
    box.addDouble("height").required();
  });
}
// _inlet_variant_struct_collections_end

const std::string input = R"(
  shapes = {
    { kind = "circle", radius = 2.5 },
    { kind = "box", width = 3.0, height = 4.0 }
  }
)";

int main()
{
  axom::slic::SimpleLogger logger;

  auto lr = std::make_unique<inlet::LuaReader>();
  lr->parseString(input);
  inlet::Inlet inlet(std::move(lr));

  // _inlet_variant_struct_collections_usage_start
  auto shapes_schema = inlet.addVariantStructArray<Shape>("shapes", "kind");
  defineShapeSchema(shapes_schema);

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  const std::vector<Shape> shapes = inlet["shapes"].get<std::vector<Shape>>();
  // _inlet_variant_struct_collections_usage_end

  for(const Shape& shape : shapes)
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

  return 0;
}
