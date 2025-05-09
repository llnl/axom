// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
  // Define 3 different datums
  axom::sina::Datum myDatum {12.34};
  std::string value = "foobar";
  axom::sina::Datum myOtherDatum {value};
  std::vector<double> scalars = {1, 2, 20.0};
  axom::sina::Datum myArrayDatum {scalars};

  // Create a record to store the datum
  axom::sina::ID myID {"my_record", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> myRecord {new axom::sina::Record {myID, "my_type"}};

  // Add the datum instances to the record
  myRecord->add("datum1", std::move(myDatum));
  myRecord->add("datum2", std::move(myOtherDatum));
  myRecord->add("datum3", std::move(myArrayDatum));

  // Query the datum from the record
  auto& data = myRecord->getData();

  // Print the keys and type of datum
  for(const auto& pair : data)
  {
    std::cout << pair.first << " is type: "
              << static_cast<std::underlying_type<axom::sina::ValueType>::type>(pair.second.getType())
              << std::endl;
  }
}