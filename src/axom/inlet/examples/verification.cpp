// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <iostream>
#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

// _inlet_verification_input_start
const std::string input = R"(
  dimensions = 2
  dim = 2 -- An example typo, should be "dimensions" and not "dim"
  vector = {
    x = 1.0,
    y = 2.0,
    z = 3.0 -- Only 2 component vectors are supported
  }
)";
// _inlet_verification_input_end

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::SimpleLogger logger;

  // Initialize Inlet
  auto lr = std::make_unique<axom::inlet::LuaReader>();
  lr->parseString(input);
  axom::inlet::Inlet myInlet(std::move(lr));

  // _inlet_workflow_defining_schema_start
  // defines a required global field named "dimensions" with a default value of 2
  myInlet.addInt("dimensions").required(true).defaultValue(2);

  // _inlet_verification_container_start
  // defines a required container named vector with an internal field named 'x'
  auto& v = myInlet.addStruct("vector").required(true);
  // _inlet_verification_container_end
  v.addDouble("x");
  // _inlet_workflow_defining_schema_end

  // _inlet_verification_strict_start
  v.strict();
  // _inlet_verification_strict_end

  // _inlet_workflow_verification_start
  v.registerVerifier([&myInlet](const axom::inlet::Container& container) -> bool {
    int dim = myInlet["dimensions"];
    bool x_present =
      container.contains("x") && (container["x"].type() == axom::inlet::InletType::Double);
    bool y_present =
      container.contains("y") && (container["y"].type() == axom::inlet::InletType::Double);
    bool z_present =
      container.contains("z") && (container["z"].type() == axom::inlet::InletType::Double);
    if(dim == 1 && x_present)
    {
      return true;
    }
    else if(dim == 2 && x_present && y_present)
    {
      return true;
    }
    else if(dim == 3 && x_present && y_present && z_present)
    {
      return true;
    }
    return false;
  });

  std::string msg;
  // We expect verification to be unsuccessful since the only Field
  // in vector is x but 2 dimensions are expected
  SLIC_INFO("This should fail due to a missing dimension:");
  myInlet.verify() ? msg = "Verification was successful\n" : msg = "Verification was unsuccessful\n";
  SLIC_INFO(msg);

  // Add required dimension to schema
  v.addDouble("y");

  // We expect the verification to succeed because vector now contains
  // both x and y to match the 2 dimensions
  SLIC_INFO("After adding the required dimension:");
  myInlet.verify() ? msg = "Verification was successful\n" : msg = "Verification was unsuccessful\n";
  SLIC_INFO(msg);
  // _inlet_workflow_verification_end

  // _inlet_workflow_accessing_data_start

  // Get dimensions if it was present in input file
  auto proxy = myInlet["dimensions"];
  if(proxy.type() == axom::inlet::InletType::Integer)
  {
    msg = "Dimensions = " + std::to_string(proxy.get<int>()) + "\n";
    SLIC_INFO(msg);
  }

  // Get vector information if it was present in input file
  bool x_found = myInlet["vector/x"].type() == axom::inlet::InletType::Double;
  bool y_found = myInlet["vector/y"].type() == axom::inlet::InletType::Double;
  if(x_found && y_found)
  {
    msg = "Vector = " + std::to_string(myInlet["vector/x"].get<double>()) + "," +
      std::to_string(myInlet["vector/y"].get<double>()) + "\n";
    SLIC_INFO(msg);
  }
  // _inlet_workflow_accessing_data_end

  // _inlet_verification_toplevel_unexpected_start
  const std::vector<std::string> all_unexpected_names =
    myInlet.unexpectedNames();  // {"dim", "vector/z"}
  // _inlet_verification_toplevel_unexpected_end

  for(const auto& name : all_unexpected_names)
  {
    SLIC_INFO("Entry '" << name << "' was not expected");
  }

  // _inlet_verification_container_unexpected_start
  const std::vector<std::string> vector_unexpected_names = v.unexpectedNames();  // {"vector/z"}
  // _inlet_verification_container_unexpected_end

  for(const auto& name : vector_unexpected_names)
  {
    SLIC_INFO("Within the vector, entry '" << name << "' was not expected");
  }

  return 0;
}
