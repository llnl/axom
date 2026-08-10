// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"
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
  input.addDouble("scale");
  input.addFunctionAsValueAlternative("scale", inlet::FunctionTag::Double, {});
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

  const bool concreteWorks = readScale("scale = 2.0") == 2.0;
  const bool callbackWorks = readScale("scale = function() return 3.0 end") == 3.0;
  return concreteWorks && callbackWorks ? 0 : 1;
}
