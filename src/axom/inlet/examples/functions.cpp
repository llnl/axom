// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/fmt.hpp"
#include "axom/inlet.hpp"
#include "axom/slic.hpp"

#include <memory>
#include <string>
#include <utility>

namespace inlet = axom::inlet;

double readScale(const std::string& lua_input, double time)
{
  auto reader = std::make_unique<inlet::LuaReader>();
  reader->parseString(lua_input);
  inlet::Inlet input(std::move(reader));

  // _inlet_or_function_schema_start
  // Accept a number or a time-dependent function that returns a number
  input.addDouble("scale").orFunction({inlet::FunctionTag::Double}).required();
  // _inlet_or_function_schema_end

  if(!input.verify())
  {
    return 0.0;
  }

  // _inlet_or_function_access_start
  // Evaluate the function at the provided time, or return the number
  const double scale = input["scale"].evaluate<double>(time);
  // _inlet_or_function_access_end

  // _inlet_or_function_type_start
  const bool is_function = input["scale"].type() == inlet::InletType::Function;
  // _inlet_or_function_type_end
  SLIC_INFO(
    axom::fmt::format("scale = {0} (supplied as a {1})", scale, is_function ? "function" : "value"));

  return scale;
}

int main()
{
  axom::slic::SimpleLogger logger;

  const double concrete = readScale("scale = 2.0", 1.5);
  const double callback = readScale("scale = function(t) return 2.0 * t end", 1.5);

  SLIC_ERROR_IF(concrete != 2.0,
                axom::fmt::format("Expected a concrete scale of 2.0, got {0}", concrete));
  SLIC_ERROR_IF(callback != 3.0,
                axom::fmt::format("Expected a callback scale of 3.0, got {0}", callback));
  return 0;
}
