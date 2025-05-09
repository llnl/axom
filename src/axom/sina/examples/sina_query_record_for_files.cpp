// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
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

  // Query the record for files
  auto& files = myRecord->getFiles();
  for(const auto& file : files)
  {
    std::cout << "File with URI '" << file.getUri() << "' has mimetype '" << file.getMimeType()
              << "' and tags '";
    for(const auto& tag : file.getTags())
    {
      std::cout << tag << " ";
    }
    std::cout << "'" << std::endl;
  }
}