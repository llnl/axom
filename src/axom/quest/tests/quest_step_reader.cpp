// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <string>

namespace fs = axom::utilities::filesystem;
namespace primal = axom::primal;
namespace quest = axom::quest;

//------------------------------------------------------------------------------
std::string pjoin(const std::string& str) { return str; }

std::string pjoin(const char* str) { return std::string(str); }

template <typename... Args>
std::string pjoin(const std::string& str, Args... args)
{
  return fs::joinPath(str, pjoin(args...));
}

template <typename... Args>
std::string pjoin(const char* str, Args... args)
{
  return fs::joinPath(std::string(str), pjoin(args...));
}

//------------------------------------------------------------------------------
bool isNearInteger(double value, double eps)
{
  const double nearest = std::round(value);
  return std::abs(value - nearest) <= eps;
}

//------------------------------------------------------------------------------
void runStepFileTest(const std::string& stepFile)
{
  const std::string fileName = pjoin(AXOM_DATA_DIR, "quest", "step", stepFile);
  SLIC_INFO(axom::fmt::format("Testing STEP file '{}'", fileName));

  quest::STEPReader stepReader;
  stepReader.setFileName(fileName);

  constexpr bool validate = false;
  const int ret = stepReader.read(validate);
  if(ret != 0)
  {
    SLIC_ERROR(axom::fmt::format("Failed to read STEP file '{}'", fileName));
  }

  const auto& patches = stepReader.getPatchArray();
  EXPECT_GT(patches.size(), 0) << "No NURBS patches were extracted from '" << fileName << "'";
  if(patches.empty())
  {
    return;
  }

  const auto shapeBbox = stepReader.getBRepBoundingBox();
  auto bboxMin = shapeBbox.getMin();
  auto bboxMax = shapeBbox.getMax();
  const auto bboxDiag = bboxMax.array() - bboxMin.array();

  const primal::WindingTolerances tol;
  constexpr double integer_eps = 1e-3;

  axom::Array<primal::Point<double, 3>> query_arr(0, 27);
  for(const double fx : {0.25, 0.5, 0.75})
  {
    for(const double fy : {0.25, 0.5, 0.75})
    {
      for(const double fz : {0.25, 0.5, 0.75})
      {
        query_arr.emplace_back(primal::Point<double, 3>({bboxMin[0] + fx * bboxDiag[0],
                                                         bboxMin[1] + fy * bboxDiag[1],
                                                         bboxMin[2] + fz * bboxDiag[2]}));
      }
    }
  }

  const auto gwn_arr = primal::winding_number(query_arr,
                                              patches,
                                              tol.edge_tol,
                                              tol.ls_tol,
                                              tol.quad_tol,
                                              tol.disk_size,
                                              tol.EPS);

  EXPECT_EQ(gwn_arr.size(), query_arr.size());
  for(int i = 0; i < gwn_arr.size() && i < query_arr.size(); ++i)
  {
    const double wn = gwn_arr[i];
    EXPECT_NEAR(wn, std::round(wn), integer_eps);
  }
}

//------------------------------------------------------------------------------
TEST(quest_step_reader, gwn_is_near_integer_on_bbox_grid)
{
  runStepFileTest("sliced_cylinder.step");
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}