// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"
#include "axom/fmt.hpp"

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace inlet = axom::inlet;

/* Input file snippets used for documentation
// _inlet_simple_types_variant_arrays_input_start

contiguous = { 42, "hello", true, 3.14 }

indexed = {
  [10] = 42,
  [20] = "hello",
  [30] = true,
  [40] = 3.14
}

// _inlet_simple_types_variant_arrays_input_end

// _inlet_simple_types_variant_dictionary_input_start

keyed = {
  [1] = 42,
  message = "hello",
  enabled = true,
  pi = 3.14
}

// _inlet_simple_types_variant_dictionary_input_end
*/

const std::string input = R"(
  contiguous = { 42, "hello", true, 3.14 }

  indexed = {
    [10] = 42,
    [20] = "hello",
    [30] = true,
    [40] = 3.14
  }

  keyed = {
    [1] = 42,
    message = "hello",
    enabled = true,
    pi = 3.14
  }
)";

int main()
{
  axom::slic::SimpleLogger logger;

  auto lr = std::make_unique<inlet::LuaReader>();
  lr->parseString(input);
  inlet::Inlet inlet(std::move(lr));

  // _inlet_simple_types_variant_collections_add_start
  inlet.addVariantArray("contiguous", "Contiguous mixed POD values");
  inlet.addVariantArray("indexed", "Integer-keyed mixed POD values");
  inlet.addVariantDictionary("keyed", "Mixed-key mixed POD values");
  // _inlet_simple_types_variant_collections_add_end

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // _inlet_simple_types_variant_arrays_access_start
  SLIC_INFO("Contiguous VariantArray as std::vector:");
  const std::vector<inlet::VariantValue> contiguous_values =
    inlet["contiguous"].get<std::vector<inlet::VariantValue>>();
  for(const inlet::VariantValue& value : contiguous_values)
  {
    std::visit([](const auto& concrete_value) {
      SLIC_INFO(axom::fmt::format("{}", concrete_value));
    }, value);
  }

  SLIC_INFO("Integer-keyed VariantArray as std::unordered_map<int, VariantValue>:");
  const std::unordered_map<int, inlet::VariantValue> indexed_values =
    inlet["indexed"].get<std::unordered_map<int, inlet::VariantValue>>();
  std::map<int, inlet::VariantValue> sorted_indexed_values(indexed_values.begin(),
                                                           indexed_values.end());
  for(const auto& entry : sorted_indexed_values)
  {
    std::visit([&entry](const auto& concrete_value) {
      SLIC_INFO(axom::fmt::format("{} = {}", entry.first, concrete_value));
    }, entry.second);
  }
  // _inlet_simple_types_variant_arrays_access_end

  // _inlet_simple_types_variant_dictionary_access_start
  SLIC_INFO("Mixed-key VariantDictionary as std::unordered_map<VariantKey, VariantValue>:");
  const std::unordered_map<inlet::VariantKey, inlet::VariantValue> keyed_values =
    inlet["keyed"].get<std::unordered_map<inlet::VariantKey, inlet::VariantValue>>();
  for(const auto& entry : keyed_values)
  {
    std::visit([&entry](const auto& concrete_value) {
      SLIC_INFO(axom::fmt::format("{} = {}", entry.first, concrete_value));
    }, entry.second);
  }
  // _inlet_simple_types_variant_dictionary_access_end

  return 0;
}
