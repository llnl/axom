// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"
#include "axom/fmt.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace inlet = axom::inlet;

/* Input file snippets used for documentation
// _inlet_simple_types_homogeneous_arrays_input_start

contiguous = { 10, 20, 30 }

indexed = {
  [10] = 100,
  [20] = 200,
  [30] = 300
}

// _inlet_simple_types_homogeneous_arrays_input_end

// _inlet_simple_types_homogeneous_dictionary_input_start

keyed = {
  [1] = 100,
  low = 10,
  medium = 20,
  high = 30
}

// _inlet_simple_types_homogeneous_dictionary_input_end
*/

const std::string input = R"(
  contiguous = { 10, 20, 30 }

  indexed = {
    [10] = 100,
    [20] = 200,
    [30] = 300
  }

  keyed = {
    [1] = 100,
    low = 10,
    medium = 20,
    high = 30
  }
)";

int main()
{
  axom::slic::SimpleLogger logger;

  auto lr = std::make_unique<inlet::LuaReader>();
  lr->parseString(input);
  inlet::Inlet inlet(std::move(lr));

  // _inlet_simple_types_homogeneous_collections_add_start
  inlet.addIntArray("contiguous", "Contiguous integer values");
  inlet.addIntArray("indexed", "Integer-keyed integer values");
  inlet.addIntDictionary("keyed", "Mixed-key integer values");
  // _inlet_simple_types_homogeneous_collections_add_end

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // _inlet_simple_types_homogeneous_arrays_access_start
  SLIC_INFO("Contiguous int array as std::vector<int>:");
  const std::vector<int> contiguous_values = inlet["contiguous"].get<std::vector<int>>();
  for(const int value : contiguous_values)
  {
    SLIC_INFO(axom::fmt::format("{}", value));
  }

  SLIC_INFO("Integer-keyed int array as std::unordered_map<int, int>:");
  const std::unordered_map<int, int> indexed_values =
    inlet["indexed"].get<std::unordered_map<int, int>>();
  for(const auto& entry : indexed_values)
  {
    SLIC_INFO(axom::fmt::format("{} = {}", entry.first, entry.second));
  }
  // _inlet_simple_types_homogeneous_arrays_access_end

  // _inlet_simple_types_homogeneous_dictionary_access_start
  SLIC_INFO("Mixed-key int dictionary as std::unordered_map<VariantKey, int>:");
  const std::unordered_map<inlet::VariantKey, int> keyed_values =
    inlet["keyed"].get<std::unordered_map<inlet::VariantKey, int>>();
  for(const auto& entry : keyed_values)
  {
    SLIC_INFO(axom::fmt::format("{} = {}", entry.first, entry.second));
  }
  // _inlet_simple_types_homogeneous_dictionary_access_end

  return 0;
}
