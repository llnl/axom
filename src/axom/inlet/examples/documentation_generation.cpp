// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// usage : ./inlet_documentation_generation_example --enableDocs --fil lua_file.lua

#include "axom/inlet.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

#include "axom/CLI11.hpp"
#include <iostream>

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::inlet::SphinxWriter;
using axom::sidre::DataStore;

void findStr(std::string path, const Inlet& inlet)
{
  auto proxy = inlet[path];
  if(proxy.type() == axom::inlet::InletType::String)
  {
    std::cout << "found " << proxy.get<std::string>();
  }
  else
  {
    std::cout << "not found ";
  }
  std::cout << std::endl;
}

void findInt(std::string path, const Inlet& inlet)
{
  auto proxy = inlet[path];
  if(proxy.type() == axom::inlet::InletType::Integer)
  {
    std::cout << "found " << proxy.get<int>();
  }
  else
  {
    std::cout << "not found ";
  }
  std::cout << std::endl;
}

void findDouble(std::string path, const Inlet& inlet)
{
  auto proxy = inlet[path];
  if(proxy.type() == axom::inlet::InletType::Double)
  {
    std::cout << "found " << proxy.get<double>();
  }
  else
  {
    std::cout << "not found ";
  }
  std::cout << std::endl;
}

void defineSchema(Inlet& inlet)
{
  // Add the description to the thermal_solver/mesh/filename Field
  auto& filename_field = inlet.addString("thermal_solver/mesh/filename", "mesh filename");
  // Set the field's required property to true
  filename_field.required();

  inlet.addInt("thermal_solver/mesh/serial", "number of serial refinements")
    .range(0, axom::numeric_limits<int>::max())
    .defaultValue(1);

  // The description for thermal_solver/mesh/parallel is left unspecified
  inlet.addInt("thermal_solver/mesh/parallel").range(1, axom::numeric_limits<int>::max()).defaultValue(1);

  inlet.addInt("thermal_solver/order", "polynomial order")
    .required()
    .range(1, axom::numeric_limits<int>::max());

  auto& timestep_field = inlet.addString("thermal_solver/timestepper", "thermal solver timestepper");
  timestep_field.defaultValue("quasistatic")
    .validValues({"quasistatic", "forwardeuler", "backwardeuler"});

  auto& coef_type_field = inlet.addString("thermal_solver/u0/type", "description for u0 type");
  coef_type_field.defaultValue("constant").validValues({"constant", "function"});

  inlet.addString("thermal_solver/u0/func", "description for u0 func").required();

  inlet.addString("thermal_solver/kappa/type", "description for kappa type")
    .required()
    .validValues({"constant", "function"});

  inlet.addDouble("thermal_solver/kappa/constant", "thermal conductivity constant").required();

  // Add description to solver container by using the addStruct function
  auto& solver_schema = inlet.addStruct("thermal_solver/solver", "linear equation solver options");

  // You can also add fields through a container

  auto& rel_tol_field = solver_schema.addDouble("rel_tol", "solver relative tolerance");
  rel_tol_field.required(false);
  rel_tol_field.defaultValue(1.e-6);
  rel_tol_field.range(0.0, axom::numeric_limits<double>::max());

  auto& abs_tol_field = solver_schema.addDouble("abs_tol", "solver absolute tolerance");
  abs_tol_field.required(true);
  abs_tol_field.defaultValue(1.e-12);
  abs_tol_field.range(0.0, axom::numeric_limits<double>::max());

  auto& print_level_field = solver_schema.addInt("print_level", "solver print/debug level");
  print_level_field.required(true);
  print_level_field.defaultValue(0);
  print_level_field.range(0, 3);

  auto& max_iter_field = solver_schema.addInt("max_iter", "maximum iteration limit");
  max_iter_field.required(false);
  max_iter_field.defaultValue(100);
  max_iter_field.range(1, axom::numeric_limits<int>::max());

  auto& dt_field = solver_schema.addDouble("dt", "time step");
  dt_field.required(true);
  dt_field.defaultValue(1);
  dt_field.range(0.0, axom::numeric_limits<double>::max());

  auto& steps_field = solver_schema.addInt("steps", "number of steps/cycles to take");
  steps_field.required(true);
  steps_field.defaultValue(1);
  steps_field.range(1, axom::numeric_limits<int>::max());
}

// Checking the contents of the passed inlet
void checkValues(const Inlet& inlet)
{
  findStr("thermal_solver/mesh/filename", inlet);
  findStr("thermal_solver/timestepper", inlet);
  findStr("thermal_solver/u0/type", inlet);
  findStr("thermal_solver/u0/func", inlet);
  findStr("thermal_solver/kappa/type", inlet);

  findInt("thermal_solver/mesh/serial", inlet);
  findInt("thermal_solver/mesh/parallel", inlet);
  findInt("thermal_solver/order", inlet);
  findInt("thermal_solver/solver/print_level", inlet);
  findInt("thermal_solver/solver/max_iter", inlet);
  findInt("thermal_solver/solver/steps", inlet);

  findDouble("thermal_solver/solver/dt", inlet);
  findDouble("thermal_solver/solver/abs_tol", inlet);
  findDouble("thermal_solver/solver/rel_tol", inlet);
  findDouble("thermal_solver/kappa/constant", inlet);

  // Verify that contents of Inlet meet the requirements of the specified schema
  if(inlet.verify())
  {
    SLIC_INFO("Inlet verify successful.");
  }
  else
  {
    SLIC_INFO("Inlet verify failed.");
  }
}

int main(int argc, char** argv)
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::SimpleLogger logger;

  // Handle command line arguments
  axom::CLI::App app {"Basic example of Axom's Inlet component"};
  bool docsEnabled {false};
  app.add_flag("--enableDocs", docsEnabled, "Enables documentation generation");

  std::string inputFileName;
  auto opt = app.add_option("--file", inputFileName, "Path to input file");
  opt->check(axom::CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  // Create inlet and parse input file data into the inlet

  DataStore ds;
  auto lr = std::make_unique<LuaReader>();
  lr->parseFile(inputFileName);
  Inlet inlet(std::move(lr), ds.getRoot(), docsEnabled);

  defineSchema(inlet);
  checkValues(inlet);

  // Generate the documentation
  // _inlet_documentation_generation_start
  inlet.write(SphinxWriter("example_doc.rst"));
  // _inlet_documentation_generation_end

  inlet.write(axom::inlet::JSONSchemaWriter("example_doc.json"));

  if(docsEnabled)
  {
    SLIC_INFO("Sphinx documentation was written to example_doc.rst\n");
    SLIC_INFO("A JSON schema was written to example_doc.json\n");
  }

  return 0;
}
