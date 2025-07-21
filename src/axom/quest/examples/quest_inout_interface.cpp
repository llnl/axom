// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_inout_interface.cpp
 *
 * Exercises the C-style quest_inout interface which determines if a point
 * is contained within an enclosed volume defined by a surface mesh.
 * Supports 3D queries against STL meshes and 2D queries against C2C contour files.
 * Note: 2D queries are only supported when Axom is configured with C2C.
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

// _quest_inout_interface_include_start
#include "axom/quest/interface/inout.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif
// _quest_inout_interface_include_end

// C/C++ includes
#include <vector>

namespace slic = axom::slic;
namespace primal = axom::primal;
namespace quest = axom::quest;

using Point2D = primal::Point<double, 2>;
using Box2D = primal::BoundingBox<double, 2>;
using Point3D = primal::Point<double, 3>;
using Box3D = primal::BoundingBox<double, 3>;
using CoordsVec = std::vector<Point3D>;

//------------------------------------------------------------------------------

/*!
 * \brief Utility function to generate a random set of query points
 * within a given bounding box
 *
 * \param [out] queryPoints The generated set of query points
 * \param [in] bbox The bounding box for the query points
 * \param [in] numPoints The number of points to generate
 * \param [in] dim The number of coordinates to generate per point (2 or 3)
 */
void generateQueryPoints(CoordsVec& queryPoints, Box3D const& bbox, int numPoints, int dim)
{
  using axom::utilities::random_real;

  queryPoints.clear();
  queryPoints.reserve(numPoints);
  for(int i = 0; i < numPoints; ++i)
  {
    double x = random_real(bbox.getMin()[0], bbox.getMax()[0]);
    double y = random_real(bbox.getMin()[1], bbox.getMax()[1]);
    double z = dim == 3 ? random_real(bbox.getMin()[2], bbox.getMax()[2]) : 0.;
    queryPoints.emplace_back(Point3D {x, y, z});
  }
}

/*!
 * \brief Utility function to initialize the logger
 */
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(axom::slic::message::Info);

  slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
  const int RLIMIT = 8;
  logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
  logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
#else
  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif  // AXOM_USE_MPI

  slic::addStreamToAllMsgLevels(logStream);
}

/*!
 * \brief Utility function to finalize the logger
 */
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

/*!
 * \brief Utility function for a clean abort that finalizes
 * quest_inout, slic and mpi
 */
void cleanAbort()
{
  if(quest::inout_initialized())
  {
    quest::inout_finalize();
  }

  finalizeLogger();

  axom::utilities::processAbort();
}

/// \brief A simple example illustrating the use of the Quest C-Style interface.
int main(int argc, char** argv)
{
  // -- Initialize MPI
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);

  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#else
  int my_rank = 0;
  int num_ranks = 1;
  AXOM_UNUSED_VAR(num_ranks);
#endif

  // -- Initialize logger
  initializeLogger();

  // -- Set up and parse command line options
  std::string fileName;
  bool isVerbose = false;
  int nQueryPoints = 100000;
  int segmentsPerKnotSpan = 25;
  double weldThresh = 1E-9;

  axom::CLI::App app {"Driver for containment query using inout API"};
  app.add_option("-i,--input", fileName)
    ->description("The input file describing a closed surface in 2D or 3D")
    ->check(axom::CLI::ExistingFile)
    ->required();
  app.add_flag("-v,--verbose", isVerbose)
    ->description("Enable/disable verbose output")
    ->capture_default_str();
  app.add_option("-t,--weld-threshold", weldThresh)
    ->description("Threshold for welding")
    ->check(axom::CLI::NonNegativeNumber)
    ->capture_default_str();
  app.add_option("-q,--num-query-points", nQueryPoints)
    ->description("Number of query points")
    ->check(axom::CLI::NonNegativeNumber)
    ->capture_default_str();
  app.add_option("-n,--segments-per-knot-span", segmentsPerKnotSpan)
    ->description("(2D only) Number of linear segments to generate per NURBS knot span")
    ->capture_default_str()
    ->check(axom::CLI::PositiveNumber);

  app.get_formatter()->column_width(50);

  try
  {
    app.parse(argc, argv);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(my_rank == 0)
    {
      retval = app.exit(e);
    }
    finalizeLogger();

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(retval);
  }

  int rc = quest::QUEST_INOUT_SUCCESS;

  // _quest_inout_interface_parameters_start
  // -- Set quest_inout parameters
  rc = quest::inout_set_verbose(isVerbose);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  // Note: contour files are only supported when Axom is configured with C2C
  using axom::utilities::string::endsWith;
  const int dim = endsWith(fileName, ".contour") ? 2 : 3;
  rc = quest::inout_set_dimension(dim);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  rc = quest::inout_set_vertex_weld_threshold(weldThresh);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  rc = quest::inout_set_segments_per_knot_span(segmentsPerKnotSpan);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }
  // _quest_inout_interface_parameters_end

  // -- Initialize quest_inout
  axom::utilities::Timer timer;
  SLIC_INFO("Initializing quest_inout...");
  timer.start();
  {
// _quest_inout_interface_init_start
#ifdef AXOM_USE_MPI
    rc = quest::inout_init(fileName, MPI_COMM_WORLD);
#else
    rc = quest::inout_init(fileName);
#endif
    // _quest_inout_interface_init_end
  }
  timer.stop();
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }
  else
  {
    SLIC_INFO("  initialization took " << timer.elapsed() << " seconds.");
  }

  // -- Query mesh bounding box and center of mass
  double bbMin[3], bbMax[3], cMass[3];
  rc = quest::inout_mesh_min_bounds(bbMin);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  rc = quest::inout_mesh_max_bounds(bbMax);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  rc = quest::inout_mesh_center_of_mass(cMass);
  if(rc != quest::QUEST_INOUT_SUCCESS)
  {
    cleanAbort();
  }

  switch(quest::inout_get_dimension())
  {
  case 2:
    SLIC_INFO("Mesh bounding box: " << Box2D(Point2D(bbMin), Point2D(bbMax)));
    SLIC_INFO("Mesh center of mass: " << Point2D(cMass));
    break;
  case 3:
  default:
    SLIC_INFO("Mesh bounding box: " << Box3D(Point3D(bbMin), Point3D(bbMax)));
    SLIC_INFO("Mesh center of mass: " << Point3D(cMass));
    break;
  }
  slic::flushStreams();

  // -- Generate query points
  Box3D bbox {Point3D(bbMin), Point3D(bbMax)};
  CoordsVec queryPoints;
  generateQueryPoints(queryPoints, bbox, nQueryPoints, dim);

  // -- Run the queries (the z-coordinate is ignored for 2D queries)
  int numInside = 0;
  SLIC_INFO(axom::fmt::format("Querying mesh with {} query points...", nQueryPoints));
  timer.reset();
  timer.start();
  for(auto& pt : queryPoints)
  {
    // _quest_inout_interface_test_start
    const bool ins = quest::inout_evaluate(pt[0], pt[1], pt[2]);
    numInside += ins ? 1 : 0;
    // _quest_inout_interface_test_end
  }
  timer.stop();

  // -- Output some query statistics
  {
    SLIC_INFO(axom::fmt::format("  queries took {} seconds.", timer.elapsed()));
    SLIC_INFO(axom::fmt::format("  query rate: {} queries per second.",
                                queryPoints.size() / timer.elapsed()));
    SLIC_INFO(
      axom::fmt::format("  {} of {} ({}%) of the query points were contained in the surface.",
                        numInside,
                        queryPoints.size(),
                        (100 * numInside) / static_cast<double>(queryPoints.size())));
  }

  // -- Finalize quest_inout
  // _quest_inout_interface_finalize_start
  quest::inout_finalize();
  // _quest_inout_interface_finalize_end

  // -- Finalize logger
  finalizeLogger();

  // -- Finalize MPI
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
