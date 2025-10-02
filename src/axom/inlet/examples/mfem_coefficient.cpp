// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/slic/core/SimpleLogger.hpp"
#include "mfem.hpp"

#include "axom/CLI11.hpp"

#include <unordered_map>
#include <iostream>

using axom::inlet::FunctionType;
using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

namespace inlet = axom::inlet;

// This example relies on mfem's support for std::function-based coefficients
#if MFEM_VERSION_MAJOR >= 4 && MFEM_VERSION_MINOR >= 2
  #define MFEM_STDFUNCTION_COEF
#endif

#ifdef MFEM_STDFUNCTION_COEF

struct BoundaryCondition
{
  std::unordered_map<int, int> attrs;
  // These act as a sort of variant without C++17
  std::unique_ptr<mfem::Coefficient> coef;
  std::unique_ptr<mfem::VectorCoefficient> vec_coef;

  // The data from the input file needed to construct a BC
  // Needed as the mesh dimension won't be in the input file, which is needed
  // to construct a coefficient
  struct InputInfo
  {
    std::function<double(const mfem::Vector&, double)> scalar_func;
    std::function<void(const mfem::Vector&, double, mfem::Vector&)> vec_func;
    std::unordered_map<int, int> attrs;
  };

  // The object gets constructed with the input info and any user-provided info
  // in this case, the dimension
  BoundaryCondition(InputInfo&& info, const int dim) : attrs(std::move(info.attrs))
  {
    if(info.scalar_func)
    {
      coef = std::make_unique<mfem::FunctionCoefficient>(std::move(info.scalar_func));
      SLIC_INFO("Created an mfem::Coefficient");
    }
    else if(info.vec_func)
    {
      vec_coef = std::make_unique<mfem::VectorFunctionCoefficient>(dim, std::move(info.vec_func));
      SLIC_INFO("Created an mfem::VectorCoefficient with dimension " << vec_coef->GetVDim());
    }
    else
    {
      SLIC_WARNING("No function present in BC input info, no coef was created");
    }
  }

  static void defineSchema(inlet::Container& schema)
  {
    schema.addIntArray("attrs", "List of boundary attributes");
    // Inlet does not support sum types, so options are added to the schema
    // for vector/scalar coefficients
    // Supported function parameter/return types are Double and Vector

    schema.addFunction(
      "vec_coef",
      inlet::FunctionTag::Vector,
      {inlet::FunctionTag::Vector, inlet::FunctionTag::Double},  // Multiple argument types
      "The function representing the BC coefficient");

    // _inlet_mfem_func_coef_start
    schema
      .addFunction("coef",
                   inlet::FunctionTag::Double,
                   {inlet::FunctionTag::Vector, inlet::FunctionTag::Double},
                   "The function representing the BC coefficient")
      // _inlet_mfem_func_coef_end
      .registerVerifier([](const inlet::Function& func) {
        // An arbitrary restriction, but this calls the function and checks its result
        return func.call<double>(inlet::FunctionType::Vector {1, 1, 1}, 1.0) < 15;
      });
  }
};

inlet::InletVector toInletVector(const mfem::Vector& vec)
{
  return inlet::InletVector {vec.GetData(), vec.Size()};
}

// Uses out-params to match MFEM semantics
void toMFEMVector(const inlet::InletVector& input, mfem::Vector& result)
{
  std::copy(input.vec.data(), input.vec.data() + input.dim, result.GetData());
}

/**
 * Example Lua definition:
 * \code{.lua}
 * bc = {
 *   attrs = {
 *      3, 4, 6, 9
 *   }
 *   coef = function (v, t)
 *     return v.x * 0.12
 *   end
 * }
 * -- or, for vector coefficients:
 * [8] = {
 *   attrs = { [4] = 14, [8] = 62, [6] = 11},
 *   vec_coef = function (v, t)
 *     scale = 0.12
 *     return v * scale
 *   end
 * }
 * \endcode
 */
template <>
struct FromInlet<BoundaryCondition::InputInfo>
{
  BoundaryCondition::InputInfo operator()(const inlet::Container& base)
  {
    BoundaryCondition::InputInfo result;
    result.attrs = base["attrs"];
    if(base.contains("coef"))
    {
      // _inlet_mfem_coef_simple_retrieve_start
      // Retrieve the function (makes a copy) to be moved into the lambda
      auto func = base["coef"].get<std::function<double(FunctionType::Vector, double)>>();
      // _inlet_mfem_coef_simple_retrieve_end
      result.scalar_func = [func(std::move(func))](const mfem::Vector& vec, double t) {
        return func(toInletVector(vec), t);
      };
    }
    else if(base.contains("vec_coef"))
    {
      auto func =
        base["vec_coef"].get<std::function<FunctionType::Vector(FunctionType::Vector, double)>>();
      result.vec_func =
        [func(std::move(func))](const mfem::Vector& input, double t, mfem::Vector& output) {
          auto ret = func(toInletVector(input), t);
          toMFEMVector(ret, output);
        };
    }
    else
    {
      SLIC_ERROR("Container did not contain a coefficient function: " << base.name());
    }
    return result;
  }
};

#endif  // MFEM_STDFUNCTION_COEF

int main(int argc, char** argv)
{
#ifdef MFEM_STDFUNCTION_COEF
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  axom::CLI::App app {"Example of Axom's Inlet component with user-defined types"};
  // Intended to be used with mfem_coef.lua
  std::string inputFileName;
  auto opt = app.add_option("--file", inputFileName, "Path to input file");
  opt->check(axom::CLI::ExistingFile);

  bool docsEnabled {false};
  app.add_flag("--docs", docsEnabled, "Enables documentation generation");

  CLI11_PARSE(app, argc, argv);

  DataStore ds;
  auto lr = std::make_unique<LuaReader>();
  lr->parseFile(inputFileName);
  Inlet inlet(std::move(lr), ds.getRoot());

  // We only need the boundary condition sub-table
  auto& bc_schema = inlet.addStructArray("bcs", "List of boundary conditions");
  BoundaryCondition::defineSchema(bc_schema);

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Read all the data into a set of boundary conditions
  auto bc_infos = inlet["bcs"].get<std::unordered_map<int, BoundaryCondition::InputInfo>>();

  // Then construct the actual boundary conditions once the mesh dimension is known
  const int dim = 3;
  std::unordered_map<int, BoundaryCondition> bcs;
  mfem::Vector input_vec(dim);
  input_vec[0] = 1;
  input_vec[1] = 2;
  input_vec[2] = 3;
  const double t = 0.0;
  mfem::Vector output_vec(dim);
  for(auto&& info : bc_infos)
  {
    if(info.second.scalar_func)
    {
      const double result = info.second.scalar_func(input_vec, t);
      SLIC_INFO(axom::fmt::format("Calling scalar function with {0} returned: {1}",
                                  toInletVector(input_vec),
                                  result));
    }
    else if(info.second.vec_func)
    {
      info.second.vec_func(input_vec, t, output_vec);
      SLIC_INFO(axom::fmt::format("Calling vector function with {0} returned: {1}",
                                  toInletVector(input_vec),
                                  toInletVector(output_vec)));
    }
    bcs.emplace(info.first, BoundaryCondition {std::move(info.second), dim});
  }

  if(docsEnabled)
  {
    const std::string docFileName = "mfem_coefficient.rst";
    inlet.write(inlet::SphinxWriter(docFileName));
    SLIC_INFO("Documentation was written to " << docFileName);
  }

  return 0;
#else   // MFEM_STDFUNCTION_COEF

  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  return 0;
#endif  // MFEM_STDFUNCTION_COEF
}
