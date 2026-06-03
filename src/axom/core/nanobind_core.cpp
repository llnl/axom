// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "axom/core/utilities/About.hpp"

namespace nb = nanobind;

namespace axom
{

NB_MODULE(pycore, m_core)
{
  m_core.doc() = "Python bindings for a small tutorial-facing slice of Axom core.";

  m_core.def("getVersion", &axom::getVersion, "Return Axom version string.");
  m_core.def("gitSHA", &axom::gitSHA, "Return Axom git SHA (may be empty).");
  m_core.def("about",
             nb::overload_cast<>(&axom::about),
             " Prints info about how Axom was configured and built.");
}

}  // namespace axom