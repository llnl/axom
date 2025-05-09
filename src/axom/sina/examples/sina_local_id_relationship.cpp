// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
  // Create IDs for Task 22 and Run 1024
  axom::sina::ID task22 {"Task_22", axom::sina::IDType::Global};
  axom::sina::ID run1024 {"Run_1024", axom::sina::IDType::Global};

  // Create the relationship and print it out
  axom::sina::Relationship myRelationship {task22, "contains", run1024};
  std::cout << myRelationship.toNode().to_json() << std::endl;

  // Create a new ID with local scope and use it to create a Record and Relationship
  axom::sina::ID myLocalID {"my_local_run", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> myRun {
    new axom::sina::Run {myLocalID, "My Sim Code", "1.2.3", "jdoe"}};
  axom::sina::Relationship myLocalRelationship {task22, "containts", myLocalID};
}