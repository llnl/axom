// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"
#include "axom/slic.hpp"

int main(void)
{
  // Initialize slic
  axom::slic::initialize();

  // Create IDs for both Task 22 and Run 1024
  axom::sina::ID task22 {"Task_22", axom::sina::IDType::Global};
  axom::sina::ID run1024 {"Run_1024", axom::sina::IDType::Global};

  // Create the relationship
  axom::sina::Relationship myRelationship {task22, "contains", run1024};

  // Compare the actual output to the expected JSON output, then print it to console
  std::string actualJsonString = myRelationship.toNode().to_json();
  std::string expectedJsonString = R"(
{
  "predicate": "contains",
  "subject": "Task_22",
  "object": "Run_1024"
})";
  SLIC_ASSERT_MSG(actualJsonString.compare(expectedJsonString) == 0,
                  "JSON output does not match expected structure.");
  std::cout << actualJsonString << std::endl;

  // Finalize slic
  axom::slic::finalize();
}
