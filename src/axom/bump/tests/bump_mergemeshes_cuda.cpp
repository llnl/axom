// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/config.hpp"
#include "axom/core.hpp"

#include "axom/bump/tests/bump_mergemeshes_impl.hpp"

TEST(bump_mergemeshes, mergemeshes_cuda) { test_mergemeshes<axom::CUDA_EXEC<256>>::test(); }

int main(int argc, char* argv[])
{
  std::string annotationMode("none");
  bool handlerEnabled = false;

  axom::CLI::App app;
#if defined(AXOM_USE_CALIPER)
  app.add_option("--caliper", annotationMode)
    ->description("caliper annotation mode, e.g. 'report'. Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  app.add_flag("--handler", handlerEnabled, "Enable Conduit handler.");
  CLI11_PARSE(app, argc, argv);

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(annotationMode);
  axom::slic::SimpleLogger logger;
  if(handlerEnabled)
  {
    conduit::utils::set_error_handler(conduit_debug_err_handler);
  }

  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  return result;
}
