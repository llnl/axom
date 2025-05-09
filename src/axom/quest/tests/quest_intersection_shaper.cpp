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

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/klee.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"
#include "axom/quest/IntersectionShaper.hpp"
#include "axom/quest/util/mesh_helpers.hpp"
#include "conduit_relay_io.hpp"

#ifndef AXOM_USE_MFEM
  #error "Quest's IntersectionShaper tests on mfem meshes require mfem library."
#endif
#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif
#include <cmath>
#include <string>
#include <vector>

// Uncomment this macro to run sequential tests (they take a long time).
//#define RUN_AXOM_SEQ_TESTS

// Uncomment this macro to regenerate baseline YAML files.
//#define GENERATE_BASELINES

// Uncomment this macro to save MFEM datasets for use in VisIt.
//#define VISUALIZE_DATASETS

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

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;
namespace slic = axom::slic;

using RuntimePolicy = axom::runtime_policy::Policy;

std::string pjoin(const std::string &path, const std::string &filename)
{
  return axom::utilities::filesystem::joinPath(path, filename);
}

void psplit(const std::string &filepath, std::string &path, std::string &filename)
{
  axom::Path p(filepath);
  path = p.dirName();
  filename = p.baseName();
}

std::string dataDirectory() { return AXOM_DATA_DIR; }

std::string testData(const std::string &filename) { return pjoin(dataDirectory(), filename); }

std::string baselineDirectory()
{
  return pjoin(pjoin(pjoin(dataDirectory(), "quest"), "regression"), "quest_intersection_shaper");
}

std::string yamlRoot(const std::string &filepath)
{
  std::string retval, path, filename;
  psplit(filepath, path, filename);
  auto idx = filename.rfind(".");
  if(idx != std::string::npos)
  {
    retval = filename.substr(0, idx);
  }
  else
  {
    retval = filename;
  }
  return retval;
}

// The caller is responsible for freeing the returned grid function.
mfem::GridFunction *newGridFunction(mfem::Mesh *mesh)
{
  const int vfOrder = 0;
  const int dim = mesh->Dimension();
  mfem::L2_FECollection *coll = new mfem::L2_FECollection(vfOrder, dim, mfem::BasisType::Positive);
  mfem::FiniteElementSpace *fes = new mfem::FiniteElementSpace(mesh, coll);
  mfem::GridFunction *gf = new mfem::GridFunction(fes);
  gf->MakeOwner(coll);
  // Initialize the values to 0.
  *gf = 0;
  return gf;
}

void makeTestMesh(sidre::MFEMSidreDataCollection &dc, bool initialMats)
{
  const int polynomialOrder = 1;
  const auto celldims = axom::NumericArray<int, 3> {20, 20, 1};
  const auto bbox = primal::BoundingBox<double, 3> {{0., 0., 0.}, {1., 1., .25}};

  auto mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, celldims, polynomialOrder, false);

  dc.SetMeshNodesName("positions");
  dc.SetMesh(mesh);

  // This mode will make 2 clean materials in the mesh.
  if(initialMats)
  {
    mfem::GridFunction *mata = newGridFunction(mesh);
    mfem::GridFunction *matb = newGridFunction(mesh);
    for(int k = 0; k < celldims[2]; k++)
    {
      for(int j = 0; j < celldims[1]; j++)
      {
        for(int i = 0; i < celldims[0]; i++)
        {
          int id = k * celldims[1] * celldims[0] + j * celldims[0] + i;
          (*mata)(id) = (i < celldims[0] / 2);
          (*matb)(id) = (i >= celldims[0] / 2);
        }
      }
    }
    // Register the fields. The dc will own them now.
    dc.RegisterField("vol_frac_a", mata);
    dc.RegisterField("vol_frac_b", matb);
  }
}

// Save Sidre as VisIt
void saveVisIt(const std::string &path, const std::string &filename, sidre::MFEMSidreDataCollection &dc)
{
  // Wrap mesh and grid functions in a VisItDataCollection and save it.
  mfem::VisItDataCollection vdc(filename, dc.GetMesh());
  if(!path.empty())
  {
    vdc.SetPrefixPath(path);
  }
  vdc.SetOwnData(false);
  vdc.SetFormat(mfem::DataCollection::SERIAL_FORMAT);
  for(auto it : dc.GetFieldMap())
  {
    if(it.first.find("vol_frac_") != std::string::npos)
    {
      vdc.RegisterField(it.first, it.second);
    }
  }
  vdc.Save();
}

// Load VisIt as Sidre
void loadVisIt(mfem::VisItDataCollection &vdc, sidre::MFEMSidreDataCollection &dc)
{
  // Wrap mesh and grid functions in a VisItDataCollection and save it.
  vdc.SetFormat(mfem::DataCollection::SERIAL_FORMAT);
  vdc.Load();
  dc.SetOwnData(false);
  dc.SetMesh(vdc.GetMesh());
  for(auto it : vdc.GetFieldMap())
  {
    if(it.first.find("vol_frac_") != std::string::npos)
    {
      dc.RegisterField(it.first, it.second);
    }
  }
}

// Turn a MFEMSidreDataCollection's fields into a simple Conduit node so
// I/O is not so problematic.
void dcToConduit(sidre::MFEMSidreDataCollection &dc, conduit::Node &n)
{
  for(auto it : dc.GetFieldMap())
  {
    // Just compare vol_frac_ grid functions.
    if(it.first.find("vol_frac_") != std::string::npos)
    {
      n[it.first].set(it.second->GetData(), it.second->Size());
    }
  }
}

bool compareConduit(const conduit::Node &n1,
                    const conduit::Node &n2,
                    double tolerance,
                    conduit::Node &info)
{
  bool same = true;
  if(n1.dtype().id() == n2.dtype().id() && n1.dtype().is_floating_point())
  {
    const auto a1 = n1.as_double_accessor();
    const auto a2 = n2.as_double_accessor();
    double maxdiff = 0.;
    for(int i = 0; i < a1.number_of_elements() && same; i++)
    {
      double diff = fabs(a1[i] - a2[i]);
      maxdiff = std::max(diff, maxdiff);
      same &= diff <= tolerance;
      if(!same)
      {
        info.append().set(axom::fmt::format("\"{}\" fields differ at index {}.", n1.name(), i));
      }
    }
    info["maxdiff"][n1.name()] = maxdiff;
  }
  else
  {
    for(int i = 0; i < n1.number_of_children() && same; i++)
    {
      const auto &n1c = n1.child(i);
      const auto &n2c = n2.fetch_existing(n1c.name());
      same &= compareConduit(n1c, n2c, tolerance, info);
    }
  }
  return same;
}

// NOTE: The baselines are read/written using Conduit directly because the
//       various data collections in Sidre, MFEM, VisIt all exhibited problems
//       either saving or loading the data.
void saveBaseline(const std::string &filename, const conduit::Node &n)
{
  std::string file_with_ext(filename + ".yaml");
  SLIC_INFO(axom::fmt::format("Save baseline ", file_with_ext));
  conduit::relay::io::save(n, file_with_ext, "yaml");
}

bool loadBaseline(const std::string &filename, conduit::Node &n)
{
  bool loaded = false;
  std::string file_with_ext(filename + ".yaml");
  SLIC_INFO(axom::fmt::format("Load baseline {}", file_with_ext));
  // Check before we read because Sidre installs a conduit error handler
  // that terminates.
  if(axom::utilities::filesystem::pathExists(file_with_ext))
  {
    conduit::relay::io::load(file_with_ext, "yaml", n);
    loaded = true;
  }
  return loaded;
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
void replacementRuleTest(const std::string &shapeFile,
                         const std::string &policyName,
                         RuntimePolicy policy,
                         double tolerance,
                         bool initialMats = false)
{
  // Make potential baseline filenames for this test. Make a policy-specific
  // baseline that we can check first. If it is not present, the next baseline
  // is tried.
  std::string baselineName(yamlRoot(shapeFile));
  if(initialMats)
  {
    baselineName += "_initial_mats";
  }
  std::vector<std::string> baselinePaths;
  // Example /path/to/axom/src/quest/tests/baseline/quest_intersection_shaper/cuda
  baselinePaths.push_back(pjoin(baselineDirectory(), policyName));
  // Example: /path/to/axom/src/quest/tests/baseline/quest_intersection_shaper/
  baselinePaths.push_back(baselineDirectory());

  // Need to make a target mesh
  SLIC_INFO(axom::fmt::format("Creating dc {}", baselineName));
  sidre::MFEMSidreDataCollection dc(baselineName, nullptr, true);
  makeTestMesh(dc, initialMats);

  // Set up shapes.
  SLIC_INFO(axom::fmt::format("Reading shape set from {}", shapeFile));
  klee::ShapeSet shapeSet(klee::readShapeSet(shapeFile));

  // Need to do the pipeline of the shaping driver.
  SLIC_INFO(axom::fmt::format("Shaping materials..."));
  const int refinementLevel = 7;
  #ifdef AXOM_USE_MPI
  // This has to happen here because the shaper gets its communicator from it.
  // If we do it before the mfem mesh is added to the data collection then the
  // data collection communicator gets set to MPI_COMM_NULL, which is bad for
  // the C2C reader.
  dc.SetComm(MPI_COMM_WORLD);
  #endif
  quest::IntersectionShaper shaper(policy, axom::INVALID_ALLOCATOR_ID, shapeSet, &dc);
  shaper.setLevel(refinementLevel);

  // Borrowed from shaping_driver.
  const klee::Dimensions shapeDim = shapeSet.getDimensions();
  for(const auto &shape : shapeSet.getShapes())
  {
    SLIC_INFO(axom::fmt::format("\tshape {} -> material {}", shape.getName(), shape.getMaterial()));

    // Load the shape from file
    shaper.loadShape(shape);
    slic::flushStreams();

    // Generate a spatial index over the shape
    shaper.prepareShapeQuery(shapeDim, shape);
    slic::flushStreams();

    // Query the mesh against this shape
    shaper.runShapeQuery(shape);
    slic::flushStreams();

    // Apply the replacement rules for this shape against the existing materials
    shaper.applyReplacementRules(shape);
    slic::flushStreams();

    // Finalize data structures associated with this shape and spatial index
    shaper.finalizeShapeQuery();
    slic::flushStreams();
  }

  // Wrap the parts of the dc data we want in the baseline as a conduit node.
  conduit::Node current;
  dcToConduit(dc, current);

  #ifdef VISUALIZE_DATASETS
  saveVisIt("", baselineName, dc);
  #endif
  #ifdef GENERATE_BASELINES
  for(const auto &path : baselinePaths)
  {
    SLIC_INFO(axom::fmt::format("Saving baseline to {}", path));
    axom::utilities::filesystem::makeDirsForPath(path);
    std::string filename(pjoin(path, baselineName));
    saveBaseline(filename, current);
  }
  #endif

  // TODO: I might want an auto compare for generating baselines so I know if I need a policy-specific baseline.

  // Need to get the MFEM mesh out and compare to expected results
  bool success = false;
  for(const auto &path : baselinePaths)
  {
    try
    {
      // Load the baseline file.
      conduit::Node info, baselineNode;
      std::string filename(pjoin(path, baselineName));
      if(loadBaseline(filename, baselineNode))
      {
        // Compare the baseline to the current DC.
        SLIC_INFO(axom::fmt::format("Comparing to baseline ", filename));
        success = compareConduit(baselineNode, current, tolerance, info);
        info.print();
        break;
      }
    }
    catch(...)
    {
      SLIC_INFO(axom::fmt::format("Could not load {} from {}!", baselineName, path));
    }
  }
  EXPECT_EQ(success, true);
}

void replacementRuleTestSet(const std::vector<std::string> &cases,
                            const std::string &policyName,
                            RuntimePolicy policy,
                            double tolerance,
                            bool initialMats = false)
{
  for(const auto &c : cases)
  {
    replacementRuleTest(testData(c), policyName, policy, tolerance, initialMats);
  }
}

void IntersectionWithErrorTolerances(const std::string &filebase,
                                     const std::string &contour,
                                     const std::string &shapeYAML,
                                     double expectedRevolvedVolume,
                                     int refinementLevel,
                                     double targetPercentError,
                                     const std::string &policyName,
                                     RuntimePolicy policy,
                                     double revolvedVolumeEPS = 1.e-4)
{
  SLIC_INFO(axom::fmt::format("Testing {} with {}", filebase, policyName));

  // Save the contour and YAML data to files so klee can read them.
  std::vector<std::string> filenames;
  filenames.emplace_back(filebase + ".contour");
  filenames.emplace_back(filebase + ".yaml");

  std::ofstream ofs;
  ofs.open(filenames[0].c_str(), std::ofstream::out);
  ofs << contour;
  ofs.close();

  ofs.open(filenames[1].c_str(), std::ofstream::out);
  ofs << shapeYAML;
  ofs.close();

  // Need to make a target mesh
  SLIC_INFO(axom::fmt::format("Creating dc {}", filebase));
  sidre::MFEMSidreDataCollection dc(filebase, nullptr, true);
  bool initialMats = false;
  makeTestMesh(dc, initialMats);

  // Set up shapes.
  SLIC_INFO(axom::fmt::format("Reading shape set from {}", filenames[1]));
  klee::ShapeSet shapeSet(klee::readShapeSet(filenames[1]));

  // Need to do the pipeline of the shaping driver.
  SLIC_INFO(axom::fmt::format("Shaping materials..."));
  #ifdef AXOM_USE_MPI
  // This has to happen here because the shaper gets its communicator from it.
  // If we do it before the mfem mesh is added to the data collection then the
  // data collection communicator gets set to MPI_COMM_NULL, which is bad for
  // the C2C reader.
  dc.SetComm(MPI_COMM_WORLD);
  #endif
  quest::IntersectionShaper shaper(policy, axom::INVALID_ALLOCATOR_ID, shapeSet, &dc);
  shaper.setLevel(refinementLevel);
  shaper.setPercentError(targetPercentError);
  shaper.setRefinementType(quest::DiscreteShape::RefinementDynamic);

  // Borrowed from shaping_driver (there should just be one shape)
  const klee::Dimensions shapeDim = shapeSet.getDimensions();
  for(const auto &shape : shapeSet.getShapes())
  {
    SLIC_INFO(axom::fmt::format("\tshape {} -> material {}", shape.getName(), shape.getMaterial()));

    // Load the shape from file
    shaper.loadShape(shape);
    slic::flushStreams();

    // Refine the shape to tolerance
    shaper.prepareShapeQuery(shapeDim, shape);
    slic::flushStreams();

    // NOTE: We do not actually run the query. We're mainly interested
    //       in how the shape was refined and whether we hit the percent
    //       error.

    // Now check the analytical revolved volume vs the value we expect. This makes
    // sure the quadrature-computed value is "close enough".
    double revolvedVolume = shaper.getRevolvedVolume();
    EXPECT_TRUE(
      axom::utilities::isNearlyEqual(revolvedVolume, expectedRevolvedVolume, revolvedVolumeEPS));

    // Now check the precent error derived from the revolved volume and the
    // linearized revolved volume
    double actualPercentError = 100. * (1. - shaper.getApproximateRevolvedVolume() / revolvedVolume);
    EXPECT_LT(actualPercentError, targetPercentError);

    // Finalize data structures associated with this shape and spatial index
    shaper.finalizeShapeQuery();
    slic::flushStreams();
  }

  // Clean up files.
  for(const auto &filename : filenames)
  {
    EXPECT_EQ(axom::utilities::filesystem::removeFile(filename), 0);
  }
}

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
#endif

//---------------------------------------------------------------------------
// Define testing functions for different modes.
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case1_seq)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case1, "seq", RuntimePolicy::seq, tolerance);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, case1_omp)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case1, "omp", RuntimePolicy::omp, tolerance);

  // Include a version that has some initial materials.
  replacementRuleTestSet(case1, "omp", RuntimePolicy::omp, tolerance, true);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, case1_cuda)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case1, "cuda", RuntimePolicy::cuda, tolerance);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, case1_hip)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case1, "hip", RuntimePolicy::hip, tolerance);
}
  #endif
#endif

// case2
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case2_seq)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case2, "seq", RuntimePolicy::seq, tolerance);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, case2_omp)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case2, "omp", RuntimePolicy::omp, tolerance);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, case2_cuda)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case2, "cuda", RuntimePolicy::cuda, tolerance);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, case2_hip)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case2, "hip", RuntimePolicy::hip, tolerance);
}
  #endif
#endif

// case3
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case3_seq)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case3, "seq", RuntimePolicy::seq, tolerance);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, case3_omp)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case3, "omp", RuntimePolicy::omp, tolerance);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, case3_cuda)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case3, "cuda", RuntimePolicy::cuda, tolerance);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, case3_hip)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case3, "hip", RuntimePolicy::hip, tolerance);
}
  #endif
#endif

// case4
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, case4_seq)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case4, "seq", RuntimePolicy::seq, tolerance);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, case4_omp)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case4, "omp", RuntimePolicy::omp, tolerance);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, case4_cuda)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case4, "cuda", RuntimePolicy::cuda, tolerance);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, case4_hip)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(case4, "hip", RuntimePolicy::hip, tolerance);
}
  #endif
#endif

// proeCase
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, proeCase_seq)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(proeCase, "seq", RuntimePolicy::seq, tolerance);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, proeCase_omp)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(proeCase, "omp", RuntimePolicy::omp, tolerance);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, proeCase_cuda)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(proeCase, "cuda", RuntimePolicy::cuda, tolerance);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, proeCase_hip)
{
  constexpr double tolerance = 1.e-10;
  replacementRuleTestSet(proeCase, "hip", RuntimePolicy::hip, tolerance);
}
  #endif
#endif

//---------------------------------------------------------------------------
// Line
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_seq) { dynamicRefinementTest_Line("seq", RuntimePolicy::seq); }
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, line_omp) { dynamicRefinementTest_Line("omp", RuntimePolicy::omp); }
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, line_cuda) { dynamicRefinementTest_Line("cuda", RuntimePolicy::cuda); }
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, line_hip) { dynamicRefinementTest_Line("hip", RuntimePolicy::hip); }
  #endif
#endif

//---------------------------------------------------------------------------
// Cone
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, cone_seq) { dynamicRefinementTest_Cone("seq", RuntimePolicy::seq); }
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, cone_omp) { dynamicRefinementTest_Cone("omp", RuntimePolicy::omp); }
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, cone_cuda) { dynamicRefinementTest_Cone("cuda", RuntimePolicy::cuda); }
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, cone_hip) { dynamicRefinementTest_Cone("hip", RuntimePolicy::hip); }
  #endif
#endif

//---------------------------------------------------------------------------
// Spline
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, spline_seq)
{
  dynamicRefinementTest_Spline("seq", RuntimePolicy::seq);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, spline_omp)
{
  dynamicRefinementTest_Spline("omp", RuntimePolicy::omp);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, spline_cuda)
{
  dynamicRefinementTest_Spline("cuda", RuntimePolicy::cuda);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, spline_hip)
{
  dynamicRefinementTest_Spline("hip", RuntimePolicy::hip);
}
  #endif
#endif

//---------------------------------------------------------------------------
// Circle
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, circle_seq)
{
  dynamicRefinementTest_Circle("seq", RuntimePolicy::seq);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, circle_omp)
{
  dynamicRefinementTest_Circle("omp", RuntimePolicy::omp);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, circle_cuda)
{
  dynamicRefinementTest_Circle("cuda", RuntimePolicy::cuda);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, circle_hip)
{
  dynamicRefinementTest_Circle("hip", RuntimePolicy::hip);
}
  #endif
#endif

//---------------------------------------------------------------------------
// LineTranslate
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_translate_seq)
{
  dynamicRefinementTest_LineTranslate("seq", RuntimePolicy::seq);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, line_translate_omp)
{
  dynamicRefinementTest_LineTranslate("omp", RuntimePolicy::omp);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, line_translate_cuda)
{
  dynamicRefinementTest_LineTranslate("cuda", RuntimePolicy::cuda);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, line_translate_hip)
{
  dynamicRefinementTest_LineTranslate("hip", RuntimePolicy::hip);
}
  #endif
#endif

//---------------------------------------------------------------------------
// LineScale
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_scale_seq)
{
  dynamicRefinementTest_LineScale("seq", RuntimePolicy::seq);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, line_scale_omp)
{
  dynamicRefinementTest_LineScale("omp", RuntimePolicy::omp);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, line_scale_cuda)
{
  dynamicRefinementTest_LineScale("cuda", RuntimePolicy::cuda);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, line_scale_hip)
{
  dynamicRefinementTest_LineScale("hip", RuntimePolicy::hip);
}
  #endif
#endif

//---------------------------------------------------------------------------
// LineRotate
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(RUN_AXOM_SEQ_TESTS)
TEST(IntersectionShaperTest, line_rotate_seq)
{
  dynamicRefinementTest_LineRotate("seq", RuntimePolicy::seq);
}
  #endif
  #if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, line_rotate_omp)
{
  dynamicRefinementTest_LineRotate("omp", RuntimePolicy::omp);
}
  #endif
  #if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, line_rotate_cuda)
{
  dynamicRefinementTest_LineRotate("cuda", RuntimePolicy::cuda);
}
  #endif
  #if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, line_rotate_hip)
{
  dynamicRefinementTest_LineRotate("hip", RuntimePolicy::hip);
}
  #endif
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
