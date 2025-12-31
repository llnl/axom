// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 * \brief Unit tests for quest's IntersectionShaper class replacement rules.
 *
 */

#include "gtest/gtest.h"

// Uncomment this macro to regenerate baseline YAML files.
//#define GENERATE_BASELINES

// Uncomment this macro to save MFEM datasets for use in VisIt.
//#define VISUALIZE_DATASETS

#include "axom/core.hpp"
#include "quest_intersection_shaper_utils.hpp"

#ifndef AXOM_USE_MFEM
  #error "Quest's IntersectionShaper tests on mfem meshes require mfem library."
#endif

ShapingTestApplication testApp;

std::vector<std::string> case1 {"shaping/case1/case1_012.yaml",
                                "shaping/case1/case1_021.yaml",
                                "shaping/case1/case1_102.yaml",
                                "shaping/case1/case1_120.yaml",
                                "shaping/case1/case1_201.yaml",
                                "shaping/case1/case1_210.yaml"};

std::vector<std::string> case2 {"shaping/case2/case2_012.yaml",
                                "shaping/case2/case2_021.yaml",
                                "shaping/case2/case2_102.yaml",
                                "shaping/case2/case2_120.yaml",
                                "shaping/case2/case2_201.yaml",
                                "shaping/case2/case2_210.yaml"};

std::vector<std::string> case3 {"shaping/case3/case3_012.yaml",
                                "shaping/case3/case3_021.yaml",
                                "shaping/case3/case3_102.yaml",
                                "shaping/case3/case3_120.yaml",
                                "shaping/case3/case3_201.yaml",
                                "shaping/case3/case3_210.yaml"};

std::vector<std::string> case4 {"shaping/case4/case4.yaml", "shaping/case4/case4_overwrite.yaml"};

std::vector<std::string> proeCase {"shaping/proeCase/proeCase1.yaml",
                                   "shaping/proeCase/proeCase2.yaml"};

constexpr double tolerance = 1.e-10;

//---------------------------------------------------------------------------
TEST(IntersectionShaperFileReadTest, loadShape_missing_stl_file_aborts)
{
  // Tests Klee shape file referencing non-existant STL mesh; should fail
  const std::string shape_yaml = R"(
dimensions: 3

shapes:
- name: missing_stl
  material: mat
  geometry:
    format: stl
    path: missing.stl
)";

  axom::utilities::filesystem::TempFile shape_file("missing_stl", ".yaml");
  shape_file.write(shape_yaml);

  const auto shapeSet = klee::readShapeSet(shape_file.getPath());
  EXPECT_FALSE(shapeSet.getShapes().empty());

  sidre::MFEMSidreDataCollection dc("missing_stl_dc", nullptr, true);
  makeTestMesh(dc, false);

  const auto policy = RuntimePolicy::seq;
  const int alloc = axom::policyToDefaultAllocatorID(policy);
  quest::IntersectionShaper shaper(policy, alloc, shapeSet, &dc);

  const auto& shape = shapeSet.getShapes().front();
  slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(shaper.loadShape(shape), slic::SlicAbortException);
}

TEST(IntersectionShaperFileReadTest, loadShape_missing_c2c_file_aborts)
{
  // Tests Klee shape file referencing non-existant c2c file; should fail
  const std::string shape_yaml = R"(
dimensions: 3

shapes:
- name: missing_c2c
  material: mat
  geometry:
    format: c2c
    path: missing.contour
)";

  axom::utilities::filesystem::TempFile shape_file("missing_c2c", ".yaml");
  shape_file.write(shape_yaml);

  const auto shapeSet = klee::readShapeSet(shape_file.getPath());
  EXPECT_FALSE(shapeSet.getShapes().empty());

  sidre::MFEMSidreDataCollection dc("missing_c2c_dc", nullptr, true);
  makeTestMesh(dc, false);

  const auto policy = RuntimePolicy::seq;
  const int alloc = axom::policyToDefaultAllocatorID(policy);
  quest::IntersectionShaper shaper(policy, alloc, shapeSet, &dc);

  const auto& shape = shapeSet.getShapes().front();
  slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(shaper.loadShape(shape), slic::SlicAbortException);
}

//---------------------------------------------------------------------------
// Define testing functions for different modes.
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case1_seq)
{
  if(testApp.selected("seq", 1))
  {
    replacementRuleTestSet(case1, "seq", RuntimePolicy::seq, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, case1_omp)
{
  if(testApp.selected("omp", 1))
  {
    replacementRuleTestSet(case1, "omp", RuntimePolicy::omp, tolerance);

    // Include a version that has some initial materials.
    replacementRuleTestSet(case1, "omp", RuntimePolicy::omp, tolerance, true);
  }
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, case1_cuda)
{
  if(testApp.selected("cuda", 1))
  {
    replacementRuleTestSet(case1, "cuda", RuntimePolicy::cuda, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, case1_hip)
{
  if(testApp.selected("hip", 1))
  {
    replacementRuleTestSet(case1, "hip", RuntimePolicy::hip, tolerance);
  }
}
#endif

// case2
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case2_seq)
{
  if(testApp.selected("seq", 2))
  {
    replacementRuleTestSet(case2, "seq", RuntimePolicy::seq, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, case2_omp)
{
  if(testApp.selected("omp", 2))
  {
    replacementRuleTestSet(case2, "omp", RuntimePolicy::omp, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, case2_cuda)
{
  if(testApp.selected("cuda", 2))
  {
    replacementRuleTestSet(case2, "cuda", RuntimePolicy::cuda, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, case2_hip)
{
  if(testApp.selected("hip", 2))
  {
    replacementRuleTestSet(case2, "hip", RuntimePolicy::hip, tolerance);
  }
}
#endif

// case3
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case3_seq)
{
  if(testApp.selected("seq", 3))
  {
    replacementRuleTestSet(case3, "seq", RuntimePolicy::seq, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, case3_omp)
{
  if(testApp.selected("omp", 3))
  {
    replacementRuleTestSet(case3, "omp", RuntimePolicy::omp, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, case3_cuda)
{
  if(testApp.selected("cuda", 3))
  {
    replacementRuleTestSet(case3, "cuda", RuntimePolicy::cuda, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, case3_hip)
{
  if(testApp.selected("hip", 3))
  {
    replacementRuleTestSet(case3, "hip", RuntimePolicy::hip, tolerance);
  }
}
#endif

// case4
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case4_seq)
{
  if(testApp.selected("seq", 4))
  {
    replacementRuleTestSet(case4, "seq", RuntimePolicy::seq, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, case4_omp)
{
  if(testApp.selected("omp", 4))
  {
    replacementRuleTestSet(case4, "omp", RuntimePolicy::omp, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, case4_cuda)
{
  if(testApp.selected("cuda", 4))
  {
    replacementRuleTestSet(case4, "cuda", RuntimePolicy::cuda, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, case4_hip)
{
  if(testApp.selected("hip", 4))
  {
    replacementRuleTestSet(case4, "hip", RuntimePolicy::hip, tolerance);
  }
}
#endif

// proeCase
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, proeCase_seq)
{
  if(testApp.selected("seq", 5))
  {
    replacementRuleTestSet(proeCase, "seq", RuntimePolicy::seq, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, proeCase_omp)
{
  if(testApp.selected("omp", 5))
  {
    replacementRuleTestSet(proeCase, "omp", RuntimePolicy::omp, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, proeCase_cuda)
{
  if(testApp.selected("cuda", 5))
  {
    replacementRuleTestSet(proeCase, "cuda", RuntimePolicy::cuda, tolerance);
  }
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, proeCase_hip)
{
  if(testApp.selected("hip", 5))
  {
    replacementRuleTestSet(proeCase, "hip", RuntimePolicy::hip, tolerance);
  }
}
#endif

//---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  int result = testApp.execute(argc, argv);

  return result;
}
