// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/STEPReader.hpp"
#include "axom/primal/operators/evaluate_integral.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <string>

namespace primal = axom::primal;
namespace quest = axom::quest;

TEST(quest_step_quadrature, tet_surface_and_volume)
{
  using Point3D = primal::Point<double, 3>;

  const std::string file = std::string(AXOM_DATA_DIR) + "/quest/step/tet.step";

  quest::STEPReader reader;
  reader.setFileName(file);
  EXPECT_EQ(reader.read(false), 0);

  const auto& patches = reader.getPatchArray();
  EXPECT_EQ(patches.size(), 4);

  auto const_integrand = [](Point3D /*x*/) -> double { return 1.0; };

  constexpr int npts = 6;
  constexpr double abs_tol = 1e-10;

  EXPECT_NEAR(primal::evaluate_surface_integral(patches, const_integrand, npts),
              8.0 * std::sqrt(3.0),
              abs_tol);
  EXPECT_NEAR(primal::evaluate_volume_integral(patches, const_integrand, npts), 8.0 / 3.0, abs_tol);
}
