// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
  // Define 2 different datums
  axom::sina::Datum myDatum {12.34};
  std::vector<double> scalars = {1, 2, 20.0};
  axom::sina::Datum myArrayDatum {scalars};

  // Set the units for one datum and the tags for the other
  myDatum.setUnits("km/s");
  std::vector<std::string> tags = {"input", "core"};
  myArrayDatum.setTags(tags);
}