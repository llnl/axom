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

  // Define 3 different datums
  double scalarValue = 12.34;
  axom::sina::Datum myDatum {scalarValue};
  std::string stringValue = "foobar";
  axom::sina::Datum myOtherDatum {stringValue};
  std::vector<double> scalarArrayValue = {1, 2, 20.0};
  axom::sina::Datum myArrayDatum {scalarArrayValue};

  // Create a record to store the datum
  axom::sina::ID myID {"my_record", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> myRecord {new axom::sina::Record {myID, "my_type"}};

  // Add the datum instances to the record
  myRecord->add("datum1", std::move(myDatum));
  myRecord->add("datum2", std::move(myOtherDatum));
  myRecord->add("datum3", std::move(myArrayDatum));

  // Query the datum
  auto& data = myRecord->getData();

  // Print the datum values
  double datum1Val = data.at("datum1").getScalar();
  SLIC_ASSERT_MSG(datum1Val == scalarValue, "Data stored in record at 'datum1' is unexpected.");
  std::cout << "datum1: " << datum1Val << std::endl;

  std::string datum2Val = data.at("datum2").getValue();
  SLIC_ASSERT_MSG(datum2Val == stringValue, "Data stored in record at 'datum2' is unexpected.");
  std::cout << "datum2: " << datum2Val << std::endl;

  std::vector<double> datum3Val = data.at("datum3").getScalarArray();
  SLIC_ASSERT_MSG(datum3Val == scalarArrayValue, "Data stored in record at 'datum3' is unexpected.");
  std::cout << "datum3: ";
  for(const auto& value : datum3Val)
  {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  // Finalize slic
  axom::slic::finalize();
}