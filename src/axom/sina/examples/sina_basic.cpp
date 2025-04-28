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
  // Create a run of "My Sim Code" version "1.2.3", which was run by "jdoe".
  // The run has an ID of "run1", which has to be unique to this file.
  axom::sina::ID runID {"run1", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> run {
    new axom::sina::Run {runID, "My Sim Code", "1.2.3", "jdoe"}};
  // Add the run to the document
  document.add(std::move(run));
  // Save the document directly to a file.
  // by specifying Protocol::JSON, we set the file to save as a JSON.
  // The third argument is optional here, it defaults to JSON.
  axom::sina::saveDocument(document, "MySinaData.json", axom::sina::Protocol::JSON);

#ifdef AXOM_USE_HDF5
  // by specifying Protocol::HDF5, we also save a copy as an HDF5 file.
  axom::sina::saveDocument(document, "MySinaData.hdf5", axom::sina::Protocol::HDF5);
#endif
}
