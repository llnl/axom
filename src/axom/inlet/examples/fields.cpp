// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <iostream>
#include <string>

#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

/* Input file snippet used for documentation
//_inlet_simple_types_fields_input_start

a_simple_bool = true
a_simple_int = 5
a_simple_double = 7.5
a_simple_string = 'such simplicity'

//_inlet_simple_types_fields_input_end
*/

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::SimpleLogger logger;

  const std::string input = R"(
    a_simple_bool = true
    a_simple_int = 5
    a_simple_double = 7.5
    a_simple_string = 'such simplicity'
  )";

  // Create Inlet Reader that supports Lua input files
  auto lr = std::make_unique<axom::inlet::LuaReader>();

  // Parse example input file string
  // Note: the Reader class also supports parseFile(std::string filepath)
  lr->parseString(input);

  // Inlet stores all input file in the Sidre DataStore
  axom::sidre::DataStore ds;

  // Create Inlet with LuaReader and the Sidre Group which Inlet will use
  axom::inlet::Inlet inlet(std::move(lr), ds.getRoot());

  // _inlet_simple_types_fields_add_start

  // Define and store the values in the input file

  // Add an optional top-level boolean
  inlet.addBool("a_simple_bool", "A description of a_simple_bool");

  // Add an optional top-level integer
  inlet.addInt("a_simple_int", "A description of a_simple_int");

  // Add a required top-level double
  inlet.addDouble("a_simple_double", "A description of a_simple_double").required();

  // Add an optional top-level string
  inlet.addString("a_simple_string", "A description of a_simple_string");

  // Add an optional top-level integer with a default value of 17 if not defined by the user
  inlet.addInt("a_defaulted_int", "An int that has a default value").defaultValue(17);

  // Add an optional top-level string not defined in the input file for example purposes
  inlet.addString("does_not_exist", "Shows that not all fields need to be present in input file");
  // _inlet_simple_types_fields_add_end

  // _inlet_simple_types_fields_access_start
  // Access values stored in the Datastore via Inlet

  // Check if input file contained info before accessing optional fields
  if(inlet.contains("a_simple_bool"))
  {
    // Access field via "[]" operator, save value first to avoid type ambiquity
    bool a_simple_bool = inlet["a_simple_bool"];
    std::cout << "a_simple_bool = " << a_simple_bool << std::endl;
  }

  if(inlet.contains("a_simple_int"))
  {
    // Access field via `get<T>` directly, no ambiquity
    std::cout << "a_simple_int = " << inlet.get<int>("a_simple_int") << std::endl;
  }

  // Because this field was marked required, we do not have to call contains before accessing
  std::cout << "a_simple_double = " << inlet.get<double>("a_simple_double") << std::endl;

  if(inlet.contains("a_simple_string"))
  {
    std::string a_simple_string = inlet["a_simple_string"];
    std::cout << "a_simple_string = " << a_simple_string << std::endl;
  }

  // If the user did not provide a value, the default value will be used.  Safe to use
  // without checking contains
  int a_defaulted_int = inlet["a_defaulted_int"];
  std::cout << "a_defaulted_int = " << a_defaulted_int << std::endl;
  // We can also verify that the user did not provided a value
  std::cout << "a_defaulted_int provided by user: " << inlet.isUserProvided("a_defaulted_int")
            << std::endl;

  // _inlet_simple_types_fields_access_end

  return 0;
}
