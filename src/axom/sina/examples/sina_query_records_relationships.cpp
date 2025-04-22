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

  // Create a new document
  axom::sina::Document document;

  // Create a record of this specific study
  // This study has an ID of "study1", which has to be unique to this file
  axom::sina::ID studyID {"study1", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> study {new axom::sina::Record {studyID, "UQ study"}};

  // Create a run of "My Sim Code" version "1.2.3", which was run by "jdoe".
  // The run has an ID of "run1", which has to be unique to this file.
  axom::sina::ID runID {"run1", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> run {
    new axom::sina::Run {runID, "My Sim Code", "1.2.3", "jdoe"}};

  // Create a relationship between the study and the run
  // Here we're saying that the study contains the run
  axom::sina::Relationship relationship {studyID, "contains", runID};

  // Add the run, study record, and relationship to the document
  document.add(std::move(run));
  document.add(std::move(study));
  document.add(relationship);

  // Query for a list of records and relationships
  auto &records = document.getRecords();
  auto &relationships = document.getRelationships();

  SLIC_ASSERT_MSG(records.size() == 2, "Unexpected number of records found.");
  std::cout << "Number of Records: " << records.size() << std::endl;
  SLIC_ASSERT_MSG(relationships.size() == 1, "Unexpected number of relationships found.");
  std::cout << "Number of Relationships: " << relationships.size() << std::endl;

  // Finalize slic
  axom::slic::finalize();
}