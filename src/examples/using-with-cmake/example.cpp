// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: example.cpp
///
//-----------------------------------------------------------------------------

#include "axom/core.hpp"
#include "axom/fmt.hpp"
#include "axom/slic.hpp"

#include <iostream>

int main()
{
#ifdef AXOM_EXAMPLE_TEST_SLIC_DEBUG
  axom::slic::SimpleLogger logger;
  int debug_evaluations = 0;

  const auto verify_debug_param = [](int AXOM_DEBUG_PARAM(value)) {
    SLIC_ASSERT_MSG(value == 4, "Unexpected value: " << value);
  };

  verify_debug_param(4);
  SLIC_DEBUG(++debug_evaluations);

  if(debug_evaluations != AXOM_EXAMPLE_EXPECT_SLIC_DEBUG)
  {
    std::cerr << "SLIC_DEBUG evaluation mismatch: expected " << AXOM_EXAMPLE_EXPECT_SLIC_DEBUG
              << ", got " << debug_evaluations << std::endl;
    return 1;
  }
#endif

  // Using fmt library exported by axom
  std::cout << axom::fmt::format("Example of using an installed version of Axom {}",
                                 axom::getVersion())
            << std::endl
            << std::endl;

  // Uses installed axom library
  axom::about();

  return 0;
}
