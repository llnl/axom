// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: naive_self_intersections.cpp
///
/// This example loads an STL mesh and checks pairs of triangles for self-intersections
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "../patch/hip_patch.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <memory>
#include <array>

//-----------------------------------------------------------------------------
/// Basic RAII utility class for initializing and finalizing slic logger
//-----------------------------------------------------------------------------
struct BasicLogger
{
  BasicLogger()
  {
    namespace slic = axom::slic;

    // Initialize the SLIC logger
    slic::initialize();
    slic::setLoggingMsgLevel(slic::message::Debug);

    // Customize logging levels and formatting
    const std::string slicFormatStr = "[lesson_02: <LEVEL>] <MESSAGE> \n";

    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr), slic::message::Error);
    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr, slicFormatStr),
                              slic::message::Warning);

    auto* compactStream = new slic::GenericOutputStream(&std::cout, slicFormatStr);
    slic::addStreamToMsgLevel(compactStream, slic::message::Info);
    slic::addStreamToMsgLevel(compactStream, slic::message::Debug);
  }

  ~BasicLogger() { axom::slic::finalize(); }
};

//-----------------------------------------------------------------------------
/// Struct to help with parsing and storing command line args
//-----------------------------------------------------------------------------
struct Input
{
  std::string mesh_file {""};
  bool verboseOutput {false};
  double weldThreshold {1e-6};
  double intersectionThreshold {1e-08};
  bool useBoundingBoxes {false};

  void parse(int argc, char** argv, axom::CLI::App& app);
  bool isVerbose() const { return verboseOutput; }
};

void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  app.add_option("-i, --infile", mesh_file)
    ->description("The input STL mesh file")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_flag("-v,--verbose", verboseOutput)
    ->description("Increase logging verbosity?")
    ->capture_default_str();

  app.add_option("--weld-threshold", weldThreshold)
    ->description(
      "Threshold to use when welding vertices.\n"
      "Will skip if not strictly positive.")
    ->capture_default_str();

  app.add_option("--intersection-threshold", intersectionThreshold)
    ->description("Threshold to use when testing for intersecting triangles")
    ->capture_default_str();

  app.add_flag("--use-bounding-boxes", useBoundingBoxes)
    ->description("Use bounding boxes to accelerate intersection query?")
    ->capture_default_str();

  app.get_formatter()->column_width(40);

  app.parse(argc, argv);  // Could throw an exception

  // Output parsed information
  SLIC_INFO(axom::fmt::format(R"(
     Parsed parameters:
      * STL mesh: '{}'
      * Threshold for welding: {}
      * Skip welding: {}
      * Threshold for intersections: {}
      * verbose logging: {}
      * use bounding boxes to accelerate query: {}
      )",
                              mesh_file,
                              weldThreshold,
                              (weldThreshold <= 0.),
                              intersectionThreshold,
                              verboseOutput,
                              useBoundingBoxes));
}

//-----------------------------------------------------------------------------
/// Basic triangle mesh to be used in our application
//-----------------------------------------------------------------------------
struct TriangleMesh
{
  using Point = axom::primal::Point<double, 3>;
  using Triangle = axom::primal::Triangle<double, 3>;
  using BoundingBox = axom::primal::BoundingBox<double, 3>;

  axom::IndexType numTriangles() const { return m_triangles.size(); }
  axom::Array<Triangle>& triangles() { return m_triangles; }
  const axom::Array<Triangle>& triangles() const { return m_triangles; }

  axom::IndexType numDegenerateTriangles() const
  {
    return static_cast<axom::IndexType>(
      std::count_if(m_degeneracies.begin(), m_degeneracies.end(), [](bool b) { return b == true; }));
  }

  bool isTriangleDegenerate(axom::IndexType idx) const { return m_degeneracies[idx]; }

  BoundingBox& meshBoundingBox() { return m_meshBoundingBox; }
  const BoundingBox& meshBoundingBox() const { return m_meshBoundingBox; }

  axom::Array<BoundingBox>& triangleBoundingBoxes() { return m_triangleBoundingBoxes; }
  const axom::Array<BoundingBox>& triangleBoundingBoxes() const { return m_triangleBoundingBoxes; }

  axom::Array<Triangle> m_triangles;
  axom::Array<bool> m_degeneracies;
  axom::Array<BoundingBox> m_triangleBoundingBoxes;
  BoundingBox m_meshBoundingBox;
};

TriangleMesh makeTriangleMesh(const std::string& stl_mesh_path, double weldThreshold)
{
  TriangleMesh triMesh;

  // load STL mesh into a mint unstructured mesh
  auto* surface_mesh =
    new axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>(3, axom::mint::TRIANGLE);
  {
    axom::utilities::Timer timer(true);

    auto reader = std::make_unique<axom::quest::STLReader>();
    reader->setFileName(stl_mesh_path);
    reader->read();
    reader->getMesh(surface_mesh);

    timer.stop();
    SLIC_INFO(axom::fmt::format("Loading the mesh took {:4.3} seconds.", timer.elapsedTimeInSec()));
  }

  // optionally weld triangle mesh
  if(weldThreshold > 0.)
  {
    axom::utilities::Timer timer(true);
    axom::quest::weldTriMeshVertices(&surface_mesh, weldThreshold);
    timer.stop();

    SLIC_INFO(axom::fmt::format("Vertex welding took {:4.3} seconds.", timer.elapsedTimeInSec()));
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "After welding, mesh has {:L} vertices and {:L} triangles.",
                                surface_mesh->getNumberOfNodes(),
                                surface_mesh->getNumberOfCells()));
  }

  // extract triangles into an axom::Array
  const int numCells = surface_mesh->getNumberOfCells();
  triMesh.m_triangles.reserve(numCells);
  TriangleMesh::Triangle tri;
  std::array<axom::IndexType, 3> triCell;
  for(int i = 0; i < numCells; ++i)
  {
    surface_mesh->getCellNodeIDs(i, triCell.data());
    surface_mesh->getNode(triCell[0], tri[0].data());
    surface_mesh->getNode(triCell[1], tri[1].data());
    surface_mesh->getNode(triCell[2], tri[2].data());

    triMesh.m_triangles.emplace_back(tri);
  }

  delete surface_mesh;
  surface_mesh = nullptr;

  // keep track of degenerate triangles
  triMesh.m_degeneracies.reserve(numCells);
  for(int i = 0; i < numCells; ++i)
  {
    const bool is_degenerate = triMesh.m_triangles[i].degenerate();
    triMesh.m_degeneracies.push_back(is_degenerate);
  }

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Mesh has {:L} degenerate triangles.",
                              triMesh.numDegenerateTriangles()));

  // compute and store triangle bounding boxes and mesh bounding box
  triMesh.m_triangleBoundingBoxes.reserve(numCells);
  for(const auto& tri : triMesh.triangles())
  {
    triMesh.m_triangleBoundingBoxes.emplace_back(axom::primal::compute_bounding_box(tri));
    triMesh.m_meshBoundingBox.addBox(triMesh.m_triangleBoundingBoxes.back());
  }

  SLIC_INFO(axom::fmt::format("Mesh bounding box is {}.", triMesh.meshBoundingBox()));

  return triMesh;
}

using IndexPair = std::pair<axom::IndexType, axom::IndexType>;

template <typename IntersectionLambda>
axom::Array<IndexPair> naiveFindIntersections(const TriangleMesh& triMesh,
                                              IntersectionLambda&& trianglesIntersect,
                                              bool verboseOutput = false)
{
  axom::Array<IndexPair> intersectionPairs;

  const int numTriangles = triMesh.numTriangles();

  for(axom::IndexType idx1 = 0; idx1 < numTriangles - 1; ++idx1)
  {
    if(triMesh.isTriangleDegenerate(idx1))
    {
      continue;
    }

    SLIC_INFO_IF(verboseOutput && idx1 % 100 == 0,
                 axom::fmt::format(axom::utilities::locale(), "Outer index {:L}", idx1));

    for(axom::IndexType idx2 = idx1 + 1; idx2 < numTriangles; ++idx2)
    {
      if(triMesh.isTriangleDegenerate(idx2))
      {
        continue;
      }

      if(trianglesIntersect(idx1, idx2))
      {
        intersectionPairs.emplace_back(std::make_pair(idx1, idx2));
      }
    }
  }

  return intersectionPairs;
}

int main(int argc, char** argv)
{
  // Initialize logger; use RAII so it will finalize at the end of the application
  BasicLogger logger;

  // Parse the command line arguments
  Input params;
  {
    axom::CLI::App app {"Naive triangle mesh intersection tester"};
    try
    {
      params.parse(argc, argv, app);
    }
    catch(const axom::CLI::ParseError& e)
    {
      return app.exit(e);
    }
  }

  // Update the logging level based on verbosity flag
  axom::slic::setLoggingMsgLevel(params.isVerbose() ? axom::slic::message::Debug
                                                    : axom::slic::message::Info);

  // Load STL mesh into local TriangleMesh struct
  SLIC_INFO(axom::fmt::format("Reading file: '{}'...\n", params.mesh_file));
  TriangleMesh mesh = makeTriangleMesh(params.mesh_file, params.weldThreshold);

  // Define a lambda to perform the intersection test on a pair of triangles in the mesh
  auto checkIntersect = [=, tol = params.intersectionThreshold, &mesh](axom::IndexType idx1,
                                                                       axom::IndexType idx2) {
    constexpr bool includeBoundaries = false;  // only use triangle interiors
    const auto& tris = mesh.triangles();

    return axom::primal::intersect(tris[idx1], tris[idx2], includeBoundaries, tol);
  };

  auto checkIntersectWithBoundingBoxes =
    [=, tol = params.intersectionThreshold, &mesh](axom::IndexType idx1, axom::IndexType idx2) {
      constexpr bool includeBoundaries = false;  // only use triangle interiors
      const auto& tris = mesh.triangles();
      const auto& bboxes = mesh.triangleBoundingBoxes();

      return axom::primal::intersect(bboxes[idx1], bboxes[idx2]) &&
        axom::primal::intersect(tris[idx1], tris[idx2], includeBoundaries, tol);
    };

  // Check for self-intersections; results are returned as an array of index pairs
  axom::utilities::Timer timer(true);
  auto intersectionPairs = params.useBoundingBoxes
    ? naiveFindIntersections(mesh, checkIntersectWithBoundingBoxes, params.isVerbose())
    : naiveFindIntersections(mesh, checkIntersect, params.isVerbose());
  timer.stop();

  SLIC_INFO(
    axom::fmt::format("Computing intersections {} took {:4.3} seconds.",
                      params.useBoundingBoxes ? "with bounding boxes" : "without bounding boxes",
                      timer.elapsedTimeInSec()));
  SLIC_INFO(axom::fmt::format("Mesh had {} intersection pairs", intersectionPairs.size()));

  SLIC_INFO_IF(
    intersectionPairs.size() > 0 && params.isVerbose(),
    axom::fmt::format("Intersecting pairs: {}\n", axom::fmt::join(intersectionPairs, ", ")));

  return 0;
}
