// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "axom/core/utilities/About.hpp"

#include <sstream>

namespace nb = nanobind;

namespace {

std::string about_string()
{
  std::ostringstream oss;
  axom::about(oss);
  return oss.str();
}

void about_print()
{
  axom::about();
}

} // namespace

NB_MODULE(pycore, m)
{
  m.doc() = "Python bindings for a small tutorial-facing slice of Axom core.";

  m.def("getVersion", &axom::getVersion, "Return Axom version string.");
  m.def("gitSHA", &axom::gitSHA, "Return Axom git SHA (may be empty).");
  m.def("about", &about_print, "Print Axom about() to stdout.");
  m.def("about_str", &about_string, "Return Axom about information as a string.");
}
