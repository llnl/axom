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

  // Create the record
  axom::sina::ID myID {"my_record", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> myRecord {new axom::sina::Record {myID, "my_type"}};

  // Create the datum with an array of strings
  std::vector<std::string> myTags {"input"};
  axom::sina::Datum myDatum {myTags};

  // Add the datum to the record
  myRecord->add("my_scalar", std::move(myDatum));

  // Compare the actual output to the expected JSON output, then print it to console
  std::string actualJsonString = myRecord->toNode().to_json();
  std::string expectedJsonString = R"(
{
  "data": 
  {
    "my_scalar": 
    {
      "value": 
      [
        "input"
      ]
    }
  },
  "type": "my_type",
  "local_id": "my_record"
})";
  SLIC_ASSERT_MSG(actualJsonString.compare(expectedJsonString) == 0,
                  "JSON output does not match expected structure.");
  std::cout << actualJsonString << std::endl;

  // Finalize slic
  axom::slic::finalize();
}