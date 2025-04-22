// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/quest.hpp"

#include "quest_test_utilities.hpp"

// Google test include
#include "gtest/gtest.h"

// C++ includes
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>  // for std::stoi

namespace mint = axom::mint;
namespace quest = axom::quest;

using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

std::string vecToString(const std::vector<int>& v)
{
  std::stringstream retval;
  for(unsigned int i = 0; i < v.size(); ++i)
  {
    retval << v[i] << "  ";
  }
  return retval.str();
}

std::string vecToString(const std::vector<std::pair<int, int>>& v)
{
  std::stringstream retval;
  for(unsigned int i = 0; i < v.size(); ++i)
  {
    retval << "(" << v[i].first << " " << v[i].second << ")  ";
  }
  return retval.str();
}

template <typename T>
void reportVectorMismatch(const std::vector<T>& standard,
                          const std::vector<T>& result,
                          const std::string& label)
{
  std::vector<T> missing, unexpected;

  std::set_difference(standard.begin(),
                      standard.end(),
                      result.begin(),
                      result.end(),
                      std::inserter(missing, missing.begin()));
  std::set_difference(result.begin(),
                      result.end(),
                      standard.begin(),
                      standard.end(),
                      std::inserter(unexpected, unexpected.begin()));

  EXPECT_TRUE(missing.size() == 0) << "Missing " << missing.size() << " " << label << ":" << std::endl
                                   << vecToString(missing);
  EXPECT_TRUE(unexpected.size() == 0)
    << "Unexpectedly, " << unexpected.size() << " extra " << label << ":" << std::endl
    << vecToString(unexpected);
}

void runIntersectTest(const std::string& tname,
                      UMesh* surface_mesh,
                      const std::vector<std::pair<int, int>>& expisect,
                      const std::vector<int>& expdegen)
{
  SCOPED_TRACE(tname);

  SLIC_INFO("Intersection test " << tname);

  std::vector<int> degenerate;
  std::vector<std::pair<int, int>> collisions;
  // Later, perhaps capture the return value as a status and report it.
  (void)quest::findTriMeshIntersections(surface_mesh, collisions, degenerate);

  // report discrepancies
  std::sort(collisions.begin(), collisions.end());
  std::sort(degenerate.begin(), degenerate.end());

  reportVectorMismatch(expisect, collisions, "triangle collisions");
  reportVectorMismatch(expdegen, degenerate, "degenerate triangles");
}

void splitStringToIntPairs(std::string& pairs, std::vector<std::pair<int, int>>& dat)
{
  if(!pairs.empty())
  {
    std::istringstream iss(pairs);
    while(iss.good())
    {
      std::pair<int, int> p;
      iss >> p.first >> p.second;
      dat.push_back(p);
    }
  }
}

void splitStringToInts(std::string& ints, std::vector<int>& dat)
{
  if(!ints.empty())
  {
    std::istringstream iss(ints);
    while(iss.good())
    {
      int i;
      iss >> i;
      dat.push_back(i);
    }
  }
}

std::string readIntersectTest(std::string& test,
                              std::string& tfname,
                              std::vector<std::pair<int, int>>& expisect,
                              std::vector<int>& expdegen,
                              quest::WatertightStatus& expwatertight,
                              int& expgenus)
{
  // given a test file path in argument test,
  // return the display name for the test (from the first line of the file).
  //
  // Output arg tfname supplies the mesh file to read in
  // (second line, path relative to test)
  //
  // Output arg expwatertight (third line) indicates if the mesh
  // is expected to be watertight after welding
  //
  // Output arg expgenus (fourth line) indicates the expected genus of the
  // mesh after welding
  //
  // Output arg expisect (fifth line) supplies the expected intersecting
  // triangles
  // Output arg expdegen (sixth line) supplies the expected degenerate
  // triangles

  std::string testdir;
  axom::utilities::filesystem::getDirName(testdir, test);

  // Read test description and file name
  std::ifstream testfile(test.c_str());
  std::string retval;
  std::getline(testfile, retval);
  std::getline(testfile, tfname);
  tfname = axom::utilities::filesystem::joinPath(testdir, tfname);

  // Read watertightness value
  {
    std::string watertight;
    std::getline(testfile, watertight);
    if(watertight == "0")
    {
      expwatertight = quest::WatertightStatus::WATERTIGHT;
    }
    else if(watertight == "1")
    {
      expwatertight = quest::WatertightStatus::NOT_WATERTIGHT;
    }
    else
    {
      expwatertight = quest::WatertightStatus::CHECK_FAILED;
    }
  }

  // Read genus value
  {
    std::string genus;
    std::getline(testfile, genus);

    expgenus = std::stoi(genus);
  }

  std::string splitline;

  // Read list of expected intersecting triangle pairs
  {
    std::getline(testfile, splitline);
    splitStringToIntPairs(splitline, expisect);
    std::sort(expisect.begin(), expisect.end());
  }

  // Read list of expected degenerate triangles
  {
    std::getline(testfile, splitline);
    splitStringToInts(splitline, expdegen);
    std::sort(expdegen.begin(), expdegen.end());
  }

  return retval;
}

std::vector<std::string> findIntersectTests()
{
  std::vector<std::string> tests;

#ifdef AXOM_DATA_DIR
  namespace fs = axom::utilities::filesystem;
  std::string catalogue = fs::joinPath(AXOM_DATA_DIR, "quest/meshtester/catalogue.txt");
  std::string testdir;
  axom::utilities::filesystem::getDirName(testdir, catalogue);

  // open file, and put each of its lines into return value tests.
  std::ifstream catfile(catalogue.c_str());
  std::string line;
  while(std::getline(catfile, line))
  {
    tests.push_back(axom::utilities::filesystem::joinPath(testdir, line));
  }
#endif

  return tests;
}

TEST(quest_mesh_tester, surfacemesh_self_intersection_intrinsic)
{
  std::vector<std::pair<int, int>> intersections;
  std::vector<int> degenerate;
  UMesh* surface_mesh = nullptr;
  std::string testname;
  std::string testdescription;

  {
    testname = "tetrahedron";
    testdescription = "Tetrahedron with no errors";

    // Construct and fill the mesh.
    // There are many ways to do this, some nicer than others.  Whether the
    // mesh has nice de-duplicated nodes or is a tiresome STL-style triangle
    // soup, it should not matter.  We will test the deduplicated triangles.
    // Nice (non-duplicated) vertices
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_tetrahedron_mesh());

    // No self-intersections or degenerate triangles
    intersections.clear();
    degenerate.clear();
    runIntersectTest(testdescription, surface_mesh, intersections, degenerate);
    delete surface_mesh;
  }

  {
    testname = "cracked tetrahedron";
    testdescription =
      "Tetrahedron with a crack but no self-intersections or degenerate "
      "triangles";

    // Construct and fill the mesh.
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_crackedtet_mesh());

    // No self-intersections or degenerate triangles
    intersections.clear();
    degenerate.clear();
    runIntersectTest(testdescription, surface_mesh, intersections, degenerate);
    delete surface_mesh;
  }

  {
    testname = "caved-in tetrahedron";
    testdescription =
      "Tetrahedron with one side intersecting two others, no degenerate "
      "triangles";

    // Construct and fill the mesh.
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_cavedtet_mesh());

    intersections.clear();
    intersections.push_back(std::make_pair(0, 1));
    intersections.push_back(std::make_pair(0, 2));
    // No degenerate triangles
    degenerate.clear();
    runIntersectTest(testdescription, surface_mesh, intersections, degenerate);
    delete surface_mesh;
  }

  {
    testname = "caved-in tet with added degenerate tris";
    testdescription =
      "Tetrahedron with one side intersecting two others, some degenerate "
      "triangles";

    // Construct and fill the mesh.
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_degen_cavedtet_mesh());

    intersections.clear();
    intersections.push_back(std::make_pair(0, 1));
    intersections.push_back(std::make_pair(0, 2));
    degenerate.clear();
    degenerate.push_back(4);
    degenerate.push_back(5);
    runIntersectTest(testdescription, surface_mesh, intersections, degenerate);
    delete surface_mesh;
  }
}

TEST(quest_mesh_tester, surfacemesh_self_intersection_ondisk)
{
  std::vector<std::string> tests = findIntersectTests();

  if(tests.size() < 1)
  {
    SLIC_INFO("*** No surface mesh self intersection tests found.");

    SUCCEED();
  }

  std::vector<std::string>::iterator it = tests.begin();
  for(; it != tests.end(); ++it)
  {
    std::string& test = *it;
    if(!axom::utilities::filesystem::pathExists(test))
    {
      SLIC_INFO("Test file does not exist; skipping: " << test);
    }
    else
    {
      std::vector<std::pair<int, int>> expisect;
      std::vector<int> expdegen;
      std::string tfname;
      quest::WatertightStatus expwatertight;
      int expgenus;

      std::string tname =
        readIntersectTest(test, tfname, expisect, expdegen, expwatertight, expgenus);

      // read in the test file into a Mesh
      quest::STLReader reader;
      reader.setFileName(tfname);
      reader.read();

      // Get surface mesh
      UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
      reader.getMesh(surface_mesh);

      runIntersectTest(tname, surface_mesh, expisect, expdegen);
      delete surface_mesh;
    }
  }
}

TEST(quest_mesh_tester, surfacemesh_watertight_intrinsic)
{
  constexpr int ON_BOUNDARY = 1;
  constexpr int INTERNAL = 0;

  UMesh* surface_mesh = nullptr;

  {
    SCOPED_TRACE("Closed tetrahedron");
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_tetrahedron_mesh());
    EXPECT_EQ(quest::WatertightStatus::WATERTIGHT, quest::isSurfaceMeshWatertight(surface_mesh));
    EXPECT_TRUE(surface_mesh->hasField("boundary", mint::CELL_CENTERED));

    // check boundary flag
    int* boundary = surface_mesh->getFieldPtr<int>("boundary", mint::CELL_CENTERED);
    const axom::IndexType numCells = surface_mesh->getNumberOfCells();
    for(axom::IndexType icell = 0; icell < numCells; ++icell)
    {
      EXPECT_TRUE(boundary[icell] == INTERNAL);
    }  // END for all cells

    delete surface_mesh;
  }

  {
    SCOPED_TRACE("Cracked tetrahedron");
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_crackedtet_mesh());
    EXPECT_EQ(quest::WatertightStatus::NOT_WATERTIGHT, quest::isSurfaceMeshWatertight(surface_mesh));
    EXPECT_TRUE(surface_mesh->hasField("boundary", mint::CELL_CENTERED));

    // check boundary flag
    int* boundary = surface_mesh->getFieldPtr<int>("boundary", mint::CELL_CENTERED);

    const axom::IndexType numCells = surface_mesh->getNumberOfCells();
    EXPECT_EQ(numCells, 4);
    EXPECT_EQ(boundary[0], ON_BOUNDARY);
    EXPECT_EQ(boundary[1], ON_BOUNDARY);
    EXPECT_EQ(boundary[2], ON_BOUNDARY);
    EXPECT_EQ(boundary[3], INTERNAL);

    delete surface_mesh;
  }

  {
    SCOPED_TRACE("Caved-in tetrahedron");
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_cavedtet_mesh());
    EXPECT_EQ(quest::WatertightStatus::NOT_WATERTIGHT, quest::isSurfaceMeshWatertight(surface_mesh));
    EXPECT_TRUE(surface_mesh->hasField("boundary", mint::CELL_CENTERED));

    // check boundary flag
    int* boundary = surface_mesh->getFieldPtr<int>("boundary", mint::CELL_CENTERED);

    const axom::IndexType numCells = surface_mesh->getNumberOfCells();
    EXPECT_EQ(numCells, 4);
    EXPECT_EQ(boundary[0], ON_BOUNDARY);
    EXPECT_EQ(boundary[1], ON_BOUNDARY);
    EXPECT_EQ(boundary[2], ON_BOUNDARY);
    EXPECT_EQ(boundary[3], INTERNAL);

    delete surface_mesh;
  }

  {
    SCOPED_TRACE("Caved-in tetrahedron with degenerate triangles");
    surface_mesh = static_cast<UMesh*>(quest::utilities::make_degen_cavedtet_mesh());
    EXPECT_EQ(quest::WatertightStatus::CHECK_FAILED, quest::isSurfaceMeshWatertight(surface_mesh));
    EXPECT_FALSE(surface_mesh->hasField("boundary", mint::CELL_CENTERED));

    delete surface_mesh;
  }
}

TEST(quest_mesh_tester, surfacemesh_watertight_ondisk)
{
  // Get the list of test cases
  std::vector<std::string> tests = findIntersectTests();

  // Test against several welding threshold value
  std::vector<double> epsilons = {1e-4, 1e-8, 1e-16};

  if(tests.size() < 1)
  {
    SLIC_INFO("*** No surface mesh watertightness tests found.");

    SUCCEED();
  }

  for(auto& test : tests)  // for each intersection test
  {
    for(double EPS : epsilons)  // for each value of epsilon
    {
      if(!axom::utilities::filesystem::pathExists(test))
      {
        SLIC_INFO("Test file does not exist; skipping: " << test);
      }
      else
      {
        std::vector<std::pair<int, int>> expisect;
        std::vector<int> expdegen;
        std::string tfname;
        quest::WatertightStatus expwatertight;
        int expgenus;

        std::string tname =
          readIntersectTest(test, tfname, expisect, expdegen, expwatertight, expgenus);

        SLIC_INFO("Running watertightness check on '" << tname << "'"
                                                      << " with EPS = " << EPS);

        // Read in the test file into a Mesh
        quest::STLReader reader;
        reader.setFileName(tfname);
        reader.read();

        // Get surface mesh
        UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
        reader.getMesh(surface_mesh);

        {
          SCOPED_TRACE(tname);

          int numOrigVerts = surface_mesh->getNumberOfNodes();
          int numOrigTris = surface_mesh->getNumberOfCells();
          EXPECT_EQ(3 * numOrigTris, numOrigVerts);

          // First weld vertices
          quest::weldTriMeshVertices(&surface_mesh, EPS);

          // Then check for holes (for STL, only meaningful after welding)
          EXPECT_EQ(expwatertight, quest::isSurfaceMeshWatertight(surface_mesh));

          /// Perform some additional checks on the welded mesh
          int numWeldedVerts = surface_mesh->getNumberOfNodes();
          int numWeldedEdges = surface_mesh->getNumberOfFaces();
          int numWeldedTris = surface_mesh->getNumberOfCells();

          // Triangle count should equal original count minus degenerate count
          EXPECT_EQ(numWeldedTris, numOrigTris - expdegen.size());

          // Check Euler characteristic for watertight meshes
          // These meshes have no boundaries
          if(expwatertight == quest::WatertightStatus::WATERTIGHT)
          {
            // Computed from genus, g, as: 2- 2g
            int expEulerCharacteristic = 2 - 2 * expgenus;

            // Computed from mesh element counts as: |V| - |E| + |F|
            int actualEulerCharacteristic = numWeldedVerts - numWeldedEdges + numWeldedTris;

            EXPECT_EQ(expEulerCharacteristic, actualEulerCharacteristic);
          }
        }

        delete surface_mesh;
      }
    }
  }
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}
