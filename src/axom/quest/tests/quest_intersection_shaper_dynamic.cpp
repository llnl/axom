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

#include "quest_intersection_shaper_utils.hpp"

#ifndef AXOM_USE_MFEM
  #error "Quest's IntersectionShaper tests on mfem meshes require mfem library."
#endif

//---------------------------------------------------------------------------
void dynamicRefinementTest_Line(const std::string &policyName, RuntimePolicy policy)
{
  const std::string contour = R"(piece = line(start=(2cm,0cm), end=(2cm,2cm))
)";

  const std::string yaml = R"(# Order 0, 1, 2
dimensions: 3

shapes:
- name: line
  material: line
  geometry:
    format: c2c
    path: line.contour
)";
  const std::string filebase = "line";
  const double expectedRevolvedVolume = 25.132741228718345;

  const std::vector<double> percentError {1., 0.1, 0.01};
  const std::vector<int> refinementLevel {7, 7, 7};
  for(size_t i = 0; i < percentError.size(); i++)
  {
    IntersectionWithErrorTolerances(filebase,
                                    contour,
                                    yaml,
                                    expectedRevolvedVolume,
                                    refinementLevel[i],
                                    percentError[i],
                                    policyName,
                                    policy);
  }
}

//---------------------------------------------------------------------------
void dynamicRefinementTest_Cone(const std::string &policyName, RuntimePolicy policy)
{
  const std::string contour = R"(piece = line(start=(2cm,0cm), end=(3cm,2cm))
)";

  const std::string yaml = R"(# Order 0, 1, 2
dimensions: 3

shapes:
- name: cone
  material: cone
  geometry:
    format: c2c
    path: cone.contour
)";
  const std::string filebase = "cone";
  const double expectedRevolvedVolume = 39.79350694547071;

  const std::vector<double> percentError {1., 0.1, 0.01};
  const std::vector<int> refinementLevel {7, 7, 7};
  for(size_t i = 0; i < percentError.size(); i++)
  {
    IntersectionWithErrorTolerances(filebase,
                                    contour,
                                    yaml,
                                    expectedRevolvedVolume,
                                    refinementLevel[i],
                                    percentError[i],
                                    policyName,
                                    policy);
  }
}

//---------------------------------------------------------------------------
void dynamicRefinementTest_Spline(const std::string &policyName, RuntimePolicy policy)
{
  const std::string contour = R"(piece = rz(units=cm,
  rz=2 0
     3 2
     3 3
)
)";

  const std::string yaml = R"(# Order 0, 1, 2
dimensions: 3

shapes:
- name: spline
  material: spline
  geometry:
    format: c2c
    path: spline.contour
)";
  const std::string filebase = "spline";
  const double expectedRevolvedVolume = 71.53270589320874;

  const std::vector<double> percentError {1., 0.1, 0.01};
  const std::vector<int> refinementLevel {7, 7, 7};
  for(size_t i = 0; i < percentError.size(); i++)
  {
    IntersectionWithErrorTolerances(filebase,
                                    contour,
                                    yaml,
                                    expectedRevolvedVolume,
                                    refinementLevel[i],
                                    percentError[i],
                                    policyName,
                                    policy,
                                    0.04);
  }
}

//---------------------------------------------------------------------------
void dynamicRefinementTest_Circle(const std::string &policyName, RuntimePolicy policy)
{
  const std::string contour =
    R"(piece = circle(origin=(0cm,0cm), radius=8cm, start=0deg, end=180deg)
)";

  const std::string yaml = R"(# Order 0, 1, 2
dimensions: 3

shapes:
- name: circle
  material: circle
  geometry:
    format: c2c
    path: circle.contour
)";
  const std::string filebase = "circle";
  const double expectedRevolvedVolume = 2144.660584850632;

  const std::vector<double> percentError {1., 0.1, 0.01};
  const std::vector<int> refinementLevel {7, 7, 7};
  for(size_t i = 0; i < percentError.size(); i++)
  {
    IntersectionWithErrorTolerances(filebase,
                                    contour,
                                    yaml,
                                    expectedRevolvedVolume,
                                    refinementLevel[i],
                                    percentError[i],
                                    policyName,
                                    policy,
                                    0.1);
  }
}

//---------------------------------------------------------------------------
void dynamicRefinementTest_LineTranslate(const std::string &policyName, RuntimePolicy policy)
{
  const std::string contour = R"(piece = line(start=(2cm,0cm), end=(2cm,2cm))
)";

  const std::string yaml = R"(# Order 0, 1, 2
dimensions: 3

shapes:
- name: line
  material: line
  geometry:
    format: c2c
    path: line.contour
    start_units: cm
    end_units: cm
    operators:
      - translate: [1., 1., 0.]
)";
  const std::string filebase = "line";
  const double expectedRevolvedVolume = 56.548667764616276;

  const std::vector<double> percentError {1., 0.1, 0.01};
  const std::vector<int> refinementLevel {7, 7, 7};
  for(size_t i = 0; i < percentError.size(); i++)
  {
    IntersectionWithErrorTolerances(filebase,
                                    contour,
                                    yaml,
                                    expectedRevolvedVolume,
                                    refinementLevel[i],
                                    percentError[i],
                                    policyName,
                                    policy);
  }
}

//---------------------------------------------------------------------------
void dynamicRefinementTest_LineScale(const std::string &policyName, RuntimePolicy policy)
{
  const std::string contour = R"(piece = line(start=(2cm,0cm), end=(2cm,2cm))
)";

  const std::string yaml = R"(# Order 0, 1, 2
dimensions: 3

shapes:
- name: line
  material: line
  geometry:
    format: c2c
    path: line.contour
    start_units: cm
    end_units: cm
    operators:
      - scale: 2.
)";
  const std::string filebase = "line";
  const double expectedRevolvedVolume = 201.06192982974676;

  const std::vector<double> percentError {1., 0.1, 0.01};
  const std::vector<int> refinementLevel {7, 7, 7};
  for(size_t i = 0; i < percentError.size(); i++)
  {
    IntersectionWithErrorTolerances(filebase,
                                    contour,
                                    yaml,
                                    expectedRevolvedVolume,
                                    refinementLevel[i],
                                    percentError[i],
                                    policyName,
                                    policy);
  }
}

//---------------------------------------------------------------------------
void dynamicRefinementTest_LineRotate(const std::string &policyName, RuntimePolicy policy)
{
  const std::string contour = R"(piece = line(start=(2cm,0cm), end=(2cm,2cm))
)";

  const std::string yaml = R"(# Order 0, 1, 2
dimensions: 3

shapes:
- name: line
  material: line
  geometry:
    format: c2c
    path: line.contour
    start_units: cm
    end_units: cm
    operators:
      - rotate: 45
        center: [0., 2., 0.]
        axis: [0., 0., 1.]
)";
  const std::string filebase = "line";
  const double expectedRevolvedVolume = 33.299824325764874;

  const std::vector<double> percentError {1., 0.1, 0.01};
  const std::vector<int> refinementLevel {7, 7, 7};
  for(size_t i = 0; i < percentError.size(); i++)
  {
    IntersectionWithErrorTolerances(filebase,
                                    contour,
                                    yaml,
                                    expectedRevolvedVolume,
                                    refinementLevel[i],
                                    percentError[i],
                                    policyName,
                                    policy);
  }
}

//---------------------------------------------------------------------------
// Line
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_seq) { dynamicRefinementTest_Line("seq", RuntimePolicy::seq); }
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, line_omp) { dynamicRefinementTest_Line("omp", RuntimePolicy::omp); }
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, line_cuda) { dynamicRefinementTest_Line("cuda", RuntimePolicy::cuda); }
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, line_hip) { dynamicRefinementTest_Line("hip", RuntimePolicy::hip); }
#endif

//---------------------------------------------------------------------------
// Cone
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, cone_seq) { dynamicRefinementTest_Cone("seq", RuntimePolicy::seq); }
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, cone_omp) { dynamicRefinementTest_Cone("omp", RuntimePolicy::omp); }
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, cone_cuda) { dynamicRefinementTest_Cone("cuda", RuntimePolicy::cuda); }
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, cone_hip) { dynamicRefinementTest_Cone("hip", RuntimePolicy::hip); }
#endif

//---------------------------------------------------------------------------
// Spline
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, spline_seq)
{
  dynamicRefinementTest_Spline("seq", RuntimePolicy::seq);
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, spline_omp)
{
  dynamicRefinementTest_Spline("omp", RuntimePolicy::omp);
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, spline_cuda)
{
  dynamicRefinementTest_Spline("cuda", RuntimePolicy::cuda);
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, spline_hip)
{
  dynamicRefinementTest_Spline("hip", RuntimePolicy::hip);
}
#endif

//---------------------------------------------------------------------------
// Circle
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, circle_seq)
{
  dynamicRefinementTest_Circle("seq", RuntimePolicy::seq);
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, circle_omp)
{
  dynamicRefinementTest_Circle("omp", RuntimePolicy::omp);
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, circle_cuda)
{
  dynamicRefinementTest_Circle("cuda", RuntimePolicy::cuda);
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, circle_hip)
{
  dynamicRefinementTest_Circle("hip", RuntimePolicy::hip);
}
#endif

//---------------------------------------------------------------------------
// LineTranslate
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_translate_seq)
{
  dynamicRefinementTest_LineTranslate("seq", RuntimePolicy::seq);
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, line_translate_omp)
{
  dynamicRefinementTest_LineTranslate("omp", RuntimePolicy::omp);
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, line_translate_cuda)
{
  dynamicRefinementTest_LineTranslate("cuda", RuntimePolicy::cuda);
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, line_translate_hip)
{
  dynamicRefinementTest_LineTranslate("hip", RuntimePolicy::hip);
}
#endif

//---------------------------------------------------------------------------
// LineScale
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_scale_seq)
{
  dynamicRefinementTest_LineScale("seq", RuntimePolicy::seq);
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, line_scale_omp)
{
  dynamicRefinementTest_LineScale("omp", RuntimePolicy::omp);
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, line_scale_cuda)
{
  dynamicRefinementTest_LineScale("cuda", RuntimePolicy::cuda);
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, line_scale_hip)
{
  dynamicRefinementTest_LineScale("hip", RuntimePolicy::hip);
}
#endif

//---------------------------------------------------------------------------
// LineRotate
#if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_rotate_seq)
{
  dynamicRefinementTest_LineRotate("seq", RuntimePolicy::seq);
}
#endif
#if defined(RUN_AXOM_OMP_TESTS)
TEST(IntersectionShaperTest, line_rotate_omp)
{
  dynamicRefinementTest_LineRotate("omp", RuntimePolicy::omp);
}
#endif
#if defined(RUN_AXOM_CUDA_TESTS)
TEST(IntersectionShaperTest, line_rotate_cuda)
{
  dynamicRefinementTest_LineRotate("cuda", RuntimePolicy::cuda);
}
#endif
#if defined(RUN_AXOM_HIP_TESTS)
TEST(IntersectionShaperTest, line_rotate_hip)
{
  dynamicRefinementTest_LineRotate("hip", RuntimePolicy::hip);
}
#endif

//---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
#ifdef AXOM_USE_MPI
  // This is needed because of Axom's c2c reader.
  MPI_Init(&argc, &argv);
  // See if this aborts right away.
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#endif

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  result = RUN_ALL_TESTS();
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif
  return result;
}
