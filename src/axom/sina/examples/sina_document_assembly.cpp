// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/sina.hpp"

int main(void)
{
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

  // Save the document directly to a file.
  // since we gave saveDocument no optional protocol parameter, it will default to JSON
  axom::sina::saveDocument(document, "MySinaData.json");

  // You can also add a document's contents to an already-dumped one, useful for capturing
  // snapshots from a code. Appending has quite a bit of functionality, see the documentation
  // for documents for more info.
  axom::sina::Document document2;
  axom::sina::ID studyID2 {"study2", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> study2 {new axom::sina::Record {studyID2, "UQ study"}};
  document2.add(std::move(study2));
  axom::sina::appendDocumentToJson("MySinaData.json", document2, 1, false);

#ifdef AXOM_USE_HDF5
  // We will also save a copy of the document as an HDF5 file
  // which can be done by passing the protocol as HDF5
  axom::sina::saveDocument(document, "MySinaData.hdf5", axom::sina::Protocol::HDF5);
#endif
}
