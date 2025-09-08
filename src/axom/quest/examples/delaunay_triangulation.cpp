// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file delaunay_triangulation.cpp
 * \brief This file contains an example that creates Delaunay triangulation
 * incrementally from a set of random points contained within an axis-aligned bounding box.
 *
 * This example can generate a triangle mesh in 3D or a tetrahedral mesh in 3D.
 * It uses SLAM IA, a topological mesh data structure.
 */

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

/// Struct to parse and contain command line arguments
struct Input
{
  std::string outputVTKFile {"delaunay_output"};
  int numRandPoints {20};
  int numOutputSteps {0};
  int dimension {2};
  std::vector<double> boundsMin;
  std::vector<double> boundsMax;

public:
  bool shouldOutputSteps() const { return numOutputSteps != 0; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-n,--nrandpt", numRandPoints)
      ->description("The number of random points")
      ->capture_default_str();

    app.add_option("-d,--dim", dimension)
      ->description(
        "The dimension of the mesh. 2 for triangle mesh in 2D; "
        "3 for tetrahedral mesh in 3D")
      ->capture_default_str();

    app.add_option("-s,--outputsteps", numOutputSteps)
      ->description(
        "Number of intermediate steps to write out as VTK files. "
        "None by default; Use -1 to write one file per insterted point")
      ->capture_default_str();

    app.add_option("-o,--outfile", outputVTKFile)
      ->description("The VTK output file")
      ->capture_default_str();

    // Optional bounding box for query region
    auto* minbb = app.add_option("--min", boundsMin)
                    ->description("Min bounds for query box (x,y[,z])")
                    ->expected(2, 3);
    auto* maxbb = app.add_option("--max", boundsMax)
                    ->description(
                      "Max bounds for query box (x,y[,z]). "
                      "Defaults to unit square/cube when not provided")
                    ->expected(2, 3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    app.get_formatter()->column_width(45);

    // could throw an exception
    app.parse(argc, argv);

    // Update number of output steps
    if(numOutputSteps == -1 || numOutputSteps > numRandPoints)
    {
      numOutputSteps = numRandPoints;
    }

    // If user doesn't provide bounds, default to unit cube
    if(boundsMin.empty())
    {
      boundsMin.clear();
      boundsMax.clear();
      for(int i = 0; i < dimension; ++i)
      {
        boundsMin.push_back(0.);
        boundsMax.push_back(1.);
      }
    }

    SLIC_INFO(axom::fmt::format(R"(Using parameter values:
    {{
      dimension: {}
      nrandpt: {}
      bounding box min: {{{}}}
      bounding box max: {{{}}}
      outfile = '{}'
      intermediate output steps: {}
    }})",
                                dimension,
                                numRandPoints,
                                axom::fmt::join(boundsMin, ", "),
                                axom::fmt::join(boundsMax, ", "),
                                outputVTKFile,
                                numOutputSteps));
  }
};

/// Utility function to create Delaunay Triangulation over random points
template <int DIM>
void run_delaunay(const Input& params)
{
  using axom::utilities::random_real;
  using axom::utilities::string::removeSuffix;

  using Delaunay = axom::quest::Delaunay<DIM>;
  using PointType = typename Delaunay::PointType;
  using BoundingBox = typename Delaunay::BoundingBox;

  const int numPoints = params.numRandPoints;
  const int numOutputVTKsteps = params.numOutputSteps;
  const std::string& outputVTKFile = removeSuffix(params.outputVTKFile, ".vtk");

  // Use a slam::ModularInt to help with bookkeeping
  const int dumpFreq = params.shouldOutputSteps() ? numPoints / numOutputVTKsteps : numPoints;
  axom::slam::ModularInt<> dumperMod(0, dumpFreq);

  BoundingBox bbox {PointType(params.boundsMin.data()), PointType(params.boundsMax.data())};

  axom::utilities::Timer timer(true);

  // Create initial Delaunay triangulation over bounding box
  Delaunay dt;
  dt.initializeBoundary(bbox);

  // Incrementally insert random points within bounding box
  for(int i = 0; i < numPoints; ++i)
  {
    PointType new_pt;
    for(int d = 0; d < DIM; ++d)
    {
      new_pt[d] = random_real(bbox.getMin()[d], bbox.getMax()[d]);
    }

    // Insert the point into the triangulation
    dt.insertPoint(new_pt);

    // Optionally, dump an intermediate mesh file
    if(params.shouldOutputSteps() && (dumperMod++ == 0 || i == numPoints - 1))
    {
      std::string fname = axom::fmt::format("{}_{:06}.vtk", outputVTKFile, i);
      dt.writeToVTKFile(fname);
    }
  }

  //Remove the starting rectangular box
  dt.removeBoundary();

  timer.stop();

  SLIC_INFO(
    axom::fmt::format("It took {} seconds to create a Delaunay complex with {} "
                      "points. Mesh has {} {}. Insertion rate of {:.1f} points per second.",
                      timer.elapsedTimeInSec(),
                      dt.getMeshData()->getNumberOfValidVertices(),
                      dt.getMeshData()->getNumberOfValidElements(),
                      DIM == 2 ? "triangles" : "tetrahedra",
                      numPoints / timer.elapsedTimeInSec()));

  // Check that the Delaunay complex is valid
  SLIC_INFO("Checking validity of Delaunay complex and underlying mesh...");
  {
    timer.reset();
    timer.start();
    dt.getMeshData()->isValid(true);
    dt.isValid(true);
    timer.stop();
    SLIC_INFO(axom::fmt::format("Validation took {} seconds", timer.elapsedTimeInSec()));
  }

  // Write the final mesh to a vtk file
  {
    std::string fname = axom::fmt::format("{}.vtk", outputVTKFile);
    SLIC_INFO(axom::fmt::format("Writing out final Delaunay complex to file '{}'", fname));
    dt.writeToVTKFile(fname);
  }

  SLIC_INFO("Done!");
}

int main(int argc, char** argv)
{
  // Initialize the SLIC logger
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  // Initialize default parameters and update with command line arguments:
  Input params;
  axom::CLI::App app {"Delaunay triangulation of a set of points in 2D or 3D"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  // Run the delaunay algorithm
  switch(params.dimension)
  {
  case 2:
    run_delaunay<2>(params);
    break;
  case 3:
    run_delaunay<3>(params);
    break;
  }

  return 0;
}
