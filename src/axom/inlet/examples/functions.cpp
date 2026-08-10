// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

#include <functional>
#include <unordered_map>

namespace inlet = axom::inlet;

// _inlet_function_value_alternative_access_start
struct Foo
{
  double bar;
};

template <>
struct FromInlet<Foo>
{
  Foo operator()(const inlet::Container& input)
  {
    if(input.containsFunctionValueAlternative("bar"))
    {
      return {input.getFunctionValueAlternative("bar").call<double>()};
    }

    return {input["bar"].get<double>()};
  }
};
// _inlet_function_value_alternative_access_end

bool runNestedCallbackExample()
{
  auto reader = std::make_unique<inlet::LuaReader>();
  reader->parseString(
    "foo = { [7] = { bar = 2 }, "
    "        [12] = { bar = function () return 3 end } }");
  inlet::Inlet inlet(std::move(reader));

  // _inlet_nested_callback_alternative_start
  auto& foo = inlet.addStructArray("foo");
  foo.addDouble("bar");
  foo.addFunctionAsValueAlternative(
    "bar",
    inlet::FunctionTag::Double,
    {});
  // _inlet_nested_callback_alternative_end

  if(!inlet.verify())
  {
    return false;
  }

  const auto values = inlet["foo"].get<std::unordered_map<int, Foo>>();
  return values.at(7).bar == 2.0 && values.at(12).bar == 3.0;
}

bool runCallbackLifetimeExample()
{
  // _inlet_callback_after_destruction_start
  std::function<double(double)> callback;
  {
    auto reader = std::make_unique<inlet::LuaReader>();
    reader->parseString("offset = 3.0; function foo (value) return value + offset end");
    inlet::Inlet inlet(std::move(reader));
    inlet.addFunction("foo",
                      inlet::FunctionTag::Double,
                      {inlet::FunctionTag::Double});
    callback = inlet["foo"].get<std::function<double(double)>>();
  }

  const double result = callback(4.0);
  // _inlet_callback_after_destruction_end
  return result == 7.0;
}

int main()
{
  axom::slic::SimpleLogger logger;

  return runNestedCallbackExample() && runCallbackLifetimeExample() ? 0 : 1;
}
