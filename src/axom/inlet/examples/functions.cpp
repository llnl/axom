// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/fmt.hpp"
#include "axom/inlet.hpp"
#include "axom/slic.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

#include <memory>
#include <string>
#include <utility>

namespace inlet = axom::inlet;

double readScale(const std::string& luaInput)
{
  auto reader = std::make_unique<inlet::LuaReader>();
  reader->parseString(luaInput);
  inlet::Inlet input(std::move(reader));

  // _inlet_function_value_alternative_schema_start
  // The alternative is declared before the concrete entry it applies to
  input.addFunctionAsValueAlternative("scale", inlet::FunctionTag::Double, {});
  input.addDouble("scale");
  // _inlet_function_value_alternative_schema_end

  if(!input.verify())
  {
    return 0.0;
  }

  // _inlet_function_value_alternative_access_start
  const auto& root = input.getGlobalContainer();
  return root.containsFunctionValueAlternative("scale")
    ? root.getFunctionValueAlternative("scale").call<double>()
    : input["scale"].get<double>();
  // _inlet_function_value_alternative_access_end
}

int main()
{
  axom::slic::SimpleLogger logger;

  const double concrete = readScale("scale = 2.0");
  const double callback = readScale("scale = function() return 3.0 end");

  SLIC_ERROR_IF(concrete != 2.0,
                axom::fmt::format("Expected a concrete scale of 2.0, got {0}", concrete));
  SLIC_ERROR_IF(callback != 3.0,
                axom::fmt::format("Expected a callback scale of 3.0, got {0}", callback));
  return 0;
}
