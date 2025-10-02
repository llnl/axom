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

  // Create 2 different files
  axom::sina::File myFile {"/path/to/file.png"};
  myFile.setMimeType("image/png");
  axom::sina::File myOtherFile {"/path/to/other/file.txt"};
  myOtherFile.setTags({"these", "are", "tags"});

  // Create a record to store the files
  axom::sina::ID myID {"my_record", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> myRecord {new axom::sina::Record {myID, "my_type"}};

  // Add the files to the record
  myRecord->add(myFile);
  myRecord->add(myOtherFile);

  // Compare the actual output to the expected JSON output, then print it to console
  std::string actualJsonString = myRecord->toNode().to_json();
  std::string expectedJsonString = R"(
{
  "type": "my_type",
  "local_id": "my_record",
  "files": 
  {
    "/path/to/other/file.txt": 
    {
      "tags": 
      [
        "these",
        "are",
        "tags"
      ]
    },
    "/path/to/file.png": 
    {
      "mimetype": "image/png"
    }
  }
})";
  SLIC_ASSERT_MSG(actualJsonString.compare(expectedJsonString) == 0,
                  "JSON output does not match expected structure.");
  std::cout << actualJsonString << std::endl;

  // Finalize slic
  axom::slic::finalize();
}