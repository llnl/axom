// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file containment_driver.cpp
 * \brief Demo of InOutOctree point containment acceleration structure over surfaces.
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/spin.hpp"
#include "axom/mint.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;

using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

//------------------------------------------------------------------------------

template <int DIM>
class ContainmentDriver
{
public:
  using CellVertIndices = primal::Point<axom::IndexType, DIM>;

  using InOutOctreeType = quest::InOutOctree<DIM>;

  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;
  using SpacePt = typename InOutOctreeType::SpacePt;
  using SpaceVector = typename InOutOctreeType::SpaceVector;
  using GridPt = typename InOutOctreeType::GridPt;
  using BlockIndex = typename InOutOctreeType::BlockIndex;
  using SpaceCell = typename InOutOctreeType::SpaceCell;

// Determine an appropriate execution policy
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)
  using ExecPolicy = axom::OMP_EXEC;
#else
  using ExecPolicy = axom::SEQ_EXEC;
#endif

  ContainmentDriver() = default;

  ~ContainmentDriver()
  {
    m_surfaceMesh.reset();
    delete m_octree;
  }

#ifdef AXOM_USE_C2C
  void loadContourMesh(const std::string& inputFile, int segmentsPerKnotSpan)
  {
    AXOM_ANNOTATE_SCOPE("load c2c");
    quest::C2CReader reader;
    reader.setFileName(inputFile);
    reader.read();

    // Create surface mesh
    m_surfaceMesh.reset(new UMesh(2, mint::SEGMENT));
    quest::LinearizeCurves lin;
    lin.getLinearMeshUniform(reader.getCurvesView(),
                             static_cast<UMesh*>(m_surfaceMesh.get()),
                             segmentsPerKnotSpan);
  }
#else
  void loadContourMesh(const std::string& inputFile, int segmentsPerKnotSpan)
  {
    AXOM_UNUSED_VAR(inputFile);
    AXOM_UNUSED_VAR(segmentsPerKnotSpan);
    SLIC_ERROR(
      "Configuration error: Loading contour files is only supported when Axom "
      "is configured with C2C support.");
  }
#endif  // AXOM_USE_C2C

  void loadSTLMesh(const std::string& inputFile)
  {
    AXOM_ANNOTATE_SCOPE("load stl");
    quest::STLReader reader;
    reader.setFileName(inputFile);
    reader.read();

    // Create surface mesh
    m_surfaceMesh.reset(new UMesh(3, mint::TRIANGLE));
    reader.getMesh(static_cast<UMesh*>(m_surfaceMesh.get()));
  }

  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh.get(); }

  int dimension() const { return DIM; }

  /// Computes the bounding box of the surface mesh
  void computeBounds()
  {
    AXOM_ANNOTATE_SCOPE("compute bounds");
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    GeometricBoundingBox meshBB;
    SpacePt pt;

    for(int i = 0; i < m_surfaceMesh->getNumberOfNodes(); ++i)
    {
      m_surfaceMesh->getNode(i, pt.data());
      m_meshBB.addPoint(pt);
    }

    SLIC_ASSERT(m_meshBB.isValid());
  }

  void initializeQueryBox(const std::vector<double>& mins, const std::vector<double>& maxs)
  {
    SLIC_ERROR_IF(m_octree == nullptr,
                  "Need to initialize InOutOctree before setting the bounding box");

    m_queryBB.clear();
    if(mins.size() == DIM && maxs.size() == DIM)
    {
      m_queryBB.addPoint(SpacePt(mins.data(), DIM));
      m_queryBB.addPoint(SpacePt(maxs.data(), DIM));
    }
    else
    {
      m_queryBB = m_octree->boundingBox();
    }

    SLIC_INFO("Bounding box for query points: " << m_queryBB);
  }

  void initializeInOutOctree()
  {
    AXOM_ANNOTATE_SCOPE("generate octree");
    m_octree = new InOutOctreeType(m_meshBB, m_surfaceMesh);
    m_octree->generateIndex();
  }

  /**
  * Query the inOutOctree using uniform grid of resolution \a gridRes
  * in region defined by the query bounding box (initialized in \a initializeQueryBox() )
  */
  void testContainmentOnRegularGrid(int gridRes, bool isBatched, bool shouldOutputMeshes)
  {
    AXOM_ANNOTATE_SCOPE(axom::fmt::format("test containment regular grid resolution {}", gridRes));

    const double* low = m_queryBB.getMin().data();
    const double* high = m_queryBB.getMax().data();
    mint::UniformMesh* umesh = (this->dimension() == 2)
      ? new mint::UniformMesh(low, high, gridRes, gridRes)
      : new mint::UniformMesh(low, high, gridRes, gridRes, gridRes);

    const int nnodes = umesh->getNumberOfNodes();
    int* containment = umesh->createField<int>("containment", mint::NODE_CENTERED);
    SLIC_ASSERT(containment != nullptr);

    axom::utilities::Timer timer;
    if(!isBatched)
    {
      timer.start();
      for(int inode = 0; inode < nnodes; ++inode)
      {
        SpacePt pt;
        umesh->getNode(inode, pt.data());

        containment[inode] = m_octree->within(pt) ? 1 : 0;
      }
      timer.stop();
    }
    else
    {
      (dimension() == 2) ? batchedPointContainment2D(umesh, timer)
                         : batchedPointContainment3D(umesh, timer);
    }

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "\tQuerying {}^{} containment field took {:.3Lf} "
                                "seconds (@ {:.0Lf} queries per second)",
                                gridRes,
                                DIM,
                                timer.elapsed(),
                                nnodes / timer.elapsed()));

    if(shouldOutputMeshes)
    {
      AXOM_ANNOTATE_SCOPE("write vtk");
      std::string fname = axom::fmt::format("gridContainment_{}.vtk", gridRes);
      mint::write_vtk(umesh, fname);
    }

    delete umesh;
  }

  void batchedPointContainment2D(mint::UniformMesh* umesh, axom::utilities::Timer& timer)
  {
    SLIC_ASSERT(umesh->getDimension() == 2);

    timer.start();

    // Allocate space for the coordinate arrays
    const int nnodes = umesh->getNumberOfNodes();
    double* x = axom::allocate<double>(nnodes);
    double* y = axom::allocate<double>(nnodes);

    // Fill the coordinate arrays
    mint::for_all_nodes<ExecPolicy, mint::xargs::xy>(
      umesh,
      AXOM_LAMBDA(axom::IndexType idx, double xx, double yy) {
        x[idx] = xx;
        y[idx] = yy;
      });

    // Loop through the points using ExecPolicy
    int* containment = umesh->getFieldPtr<int>("containment", mint::NODE_CENTERED);
    SLIC_ASSERT(containment != nullptr);

    axom::for_all<ExecPolicy>(0, nnodes, [&](axom::IndexType idx) {
      const bool inside = m_octree->within(SpacePt {x[idx], y[idx]});
      containment[idx] = inside ? 1 : 0;
    });

    // Deallocate the coordinate arrays
    axom::deallocate(x);
    axom::deallocate(y);

    timer.stop();
  }

  void batchedPointContainment3D(mint::UniformMesh* umesh, axom::utilities::Timer& timer)
  {
    SLIC_ASSERT(umesh->getDimension() == 3);

    timer.start();

    // Allocate space for the coordinate arrays
    const int nnodes = umesh->getNumberOfNodes();
    double* x = axom::allocate<double>(nnodes);
    double* y = axom::allocate<double>(nnodes);
    double* z = axom::allocate<double>(nnodes);

    // Fill the coordinate arrays
    mint::for_all_nodes<ExecPolicy, mint::xargs::xyz>(
      umesh,
      AXOM_LAMBDA(axom::IndexType idx, double xx, double yy, double zz) {
        x[idx] = xx;
        y[idx] = yy;
        z[idx] = zz;
      });

    // Loop through the points using ExecPolicy
    int* containment = umesh->getFieldPtr<int>("containment", mint::NODE_CENTERED);
    SLIC_ASSERT(containment != nullptr);

    axom::for_all<ExecPolicy>(0, nnodes, [&](axom::IndexType idx) {
      const bool inside = m_octree->within(SpacePt {x[idx], y[idx], z[idx]});
      containment[idx] = inside ? 1 : 0;
    });

    // Deallocate the coordinate arrays
    axom::deallocate(x);
    axom::deallocate(y);
    axom::deallocate(z);

    timer.stop();
  }

private:
  /**
  * \brief Extracts the vertex indices of cell \a cellIndex from the mesh
  */
  CellVertIndices getCellVertIndices(axom::IndexType cellIndex) const
  {
    SLIC_ASSERT(m_surfaceMesh != nullptr);
    SLIC_ASSERT(cellIndex >= 0 && cellIndex < m_surfaceMesh->getNumberOfCells());

    CellVertIndices tvInd;
    m_surfaceMesh->getCellNodeIDs(cellIndex, tvInd.data());
    return tvInd;
  }

  /**
  * \brief Extracts the positions of a cell's vertices from the mesh
  * \return The cell vertex positions in a SpaceCell instance
  */
  SpaceCell getMeshCell(const CellVertIndices& vertIndices) const
  {
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    SpaceCell cell;
    for(int i = 0; i < DIM; ++i)
    {
      m_surfaceMesh->getNode(vertIndices[i], cell[i].data());
    }
    return cell;
  }

  double cellMeasure(const primal::Segment<double, 2>& cell) const { return cell.length(); }
  double cellMeasure(const primal::Triangle<double, 3>& cell) const { return cell.area(); }

public:
  /**
  * \brief Computes some statistics about the surface mesh.
  *
  * Specifically, computes histograms (and ranges) of the edge lengths and
  * cell areas on a logarithmic scale and logs the results
  */
  void printSurfaceStats() const
  {
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Mesh has {:L} nodes and {:L} cells.",
                                m_surfaceMesh->getNumberOfNodes(),
                                m_surfaceMesh->getNumberOfCells()));
    SLIC_INFO("Mesh bounding box: " << m_meshBB);

    SpacePt pt;

    using MinMaxRange = primal::BoundingBox<double, 1>;
    using LengthType = MinMaxRange::PointType;

    MinMaxRange meshEdgeLenRange;
    MinMaxRange meshCellAreaRange;
    const int nCells = m_surfaceMesh->getNumberOfCells();
    using CellIdxSet = std::set<axom::IndexType>;
    CellIdxSet badCells;

    // simple binning based on the exponent
    using LogHistogram = std::map<int, int>;
    LogHistogram edgeLenHist;  // Create histogram of edge lengths (log scale)
    LogHistogram areaHist;     // Create histogram of triangle areas (log scale)

    using LogRangeMap = std::map<int, MinMaxRange>;
    LogRangeMap edgeLenRangeMap;  // Tracks range of edge lengths at each scale
    LogRangeMap areaRangeMap;     // Tracks range of triangle areas at each scale

    int expBase2;

    // Traverse mesh cells and bin the edge lengths and areas
    for(int i = 0; i < nCells; ++i)
    {
      // Get the indices and positions of the cell's vertices
      CellVertIndices vertIndices = getCellVertIndices(i);
      SpaceCell cell = getMeshCell(vertIndices);

      if(dimension() == 3)
      {
        // Compute edge stats -- note edges are double counted
        for(int j = 0; j < DIM; ++j)
        {
          double len = SpaceVector(cell[j], cell[(j + 1) % DIM]).norm();
          if(axom::utilities::isNearlyEqual(len, 0.))
          {
            badCells.insert(i);
          }
          else
          {
            LengthType edgeLen(len);
            meshEdgeLenRange.addPoint(edgeLen);
            std::frexp(len, &expBase2);
            edgeLenHist[expBase2]++;
            edgeLenRangeMap[expBase2].addPoint(edgeLen);
          }
        }
      }

      // Compute cell area stats
      double area = cellMeasure(cell);
      if(axom::utilities::isNearlyEqual(area, 0.))
      {
        badCells.insert(i);
      }
      else
      {
        LengthType cellArea(area);
        meshCellAreaRange.addPoint(cellArea);
        std::frexp(area, &expBase2);
        areaHist[expBase2]++;
        areaRangeMap[expBase2].addPoint(cellArea);
      }
    }

    // Log the results
    const int nVerts = m_surfaceMesh->getNumberOfNodes();
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Mesh has {:L} vertices  and {:L} cells.",
                                nVerts,
                                nCells));

    SLIC_INFO_IF(dimension() == 3, "Edge length range: " << meshEdgeLenRange);
    SLIC_INFO("Cell area range is: " << meshCellAreaRange);

    if(dimension() == 3)
    {
      axom::fmt::memory_buffer edgeHistStr;
      axom::fmt::format_to(std::back_inserter(edgeHistStr), "Edge length histogram (lg-arithmic): ");
      for(auto it = edgeLenHist.begin(); it != edgeLenHist.end(); ++it)
      {
        axom::fmt::format_to(std::back_inserter(edgeHistStr),
                             "\n\texp: {}\tcount: {}\tRange: {}",
                             it->first,
                             it->second / 2,
                             edgeLenRangeMap[it->first]);
      }
      SLIC_DEBUG(axom::fmt::to_string(edgeHistStr));
    }

    axom::fmt::memory_buffer cellHistStr;
    axom::fmt::format_to(std::back_inserter(cellHistStr), "Cell areas histogram (lg-arithmic): ");
    for(auto it = areaHist.begin(); it != areaHist.end(); ++it)
    {
      axom::fmt::format_to(std::back_inserter(cellHistStr),
                           "\n\texp: {}\tcount: {}\tRange: {}",
                           it->first,
                           it->second,
                           areaRangeMap[it->first]);
    }
    SLIC_DEBUG(axom::fmt::to_string(cellHistStr));

    if(!badCells.empty())
    {
      axom::fmt::memory_buffer badCellStr;
      axom::fmt::format_to(std::back_inserter(badCellStr),
                           "The following cell(s) have zero area/edge lengths:");
      for(auto it = badCells.begin(); it != badCells.end(); ++it)
      {
        axom::fmt::format_to(std::back_inserter(badCellStr), "\n\tCell {}", *it);
        CellVertIndices vertIndices;
        m_surfaceMesh->getCellNodeIDs(*it, vertIndices.data());

        SpacePt vertPos;
        for(int j = 0; j < DIM; ++j)
        {
          m_surfaceMesh->getNode(vertIndices[j], vertPos.data());
          axom::fmt::format_to(std::back_inserter(badCellStr),
                               "\n\t\t vId: {} @ position: {}",
                               vertIndices[j],
                               vertPos);
        }
      }
      SLIC_DEBUG(axom::fmt::to_string(badCellStr));
    }
  }

private:
  std::shared_ptr<mint::Mesh> m_surfaceMesh {nullptr};
  InOutOctreeType* m_octree {nullptr};
  GeometricBoundingBox m_meshBB;
  GeometricBoundingBox m_queryBB;
};

/** Struct to parse and store the input parameters */
struct Input
{
public:
  std::string inputFile;
  int maxQueryLevel {7};
  int samplesPerKnotSpan {25};
  std::vector<double> queryBoxMins;
  std::vector<double> queryBoxMaxs;
  std::string annotationMode {"none"};

private:
  bool m_verboseOutput {false};
  bool m_use_batched_query {false};

public:
  Input()
  {
// set default stl file, when the 'data' submodule is available
#ifdef AXOM_DATA_DIR
    {
      namespace fs = axom::utilities::filesystem;
      const auto dir = fs::joinPath(AXOM_DATA_DIR, "quest");
      inputFile = fs::joinPath(dir, "plane_simp.stl");
    }
#endif

// increase default for max query resolution in release builds
#ifndef AXOM_DEBUG
    maxQueryLevel += 2;
#endif
  }

  bool isInput2D() const
  {
    using axom::utilities::string::endsWith;
    return endsWith(inputFile, ".contour");
  }

  bool isVerbose() const { return m_verboseOutput; }

  bool useBatchedQuery() const { return m_use_batched_query; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-i,--input", inputFile, "Path to input file")->check(axom::CLI::ExistingFile);

    app
      .add_flag("-v,--verbose",
                m_verboseOutput,
                "Enable/disable verbose output, "
                "including outputting generated containment grids.")
      ->capture_default_str();

    app
      .add_option("-l,--levels",
                  maxQueryLevel,
                  "Max query resolution. \n"
                  "Will query uniform grids at levels 1 through the provided level")
      ->capture_default_str()
      ->check(axom::CLI::PositiveNumber);

    // Optional bounding box for query region
    auto* minbb =
      app.add_option("--min", queryBoxMins, "Min bounds for query box (x,y[,z])")->expected(2, 3);
    auto* maxbb =
      app.add_option("--max", queryBoxMaxs, "Max bounds for query box (x,y[,z])")->expected(2, 3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    app
      .add_flag("--batched",
                m_use_batched_query,
                "uses a single batched query on all points instead of many "
                "individual queries")
      ->capture_default_str();

    app
      .add_option("-n,--segments-per-knot-span",
                  samplesPerKnotSpan,
                  "(2D only) Number of linear segments to generate per NURBS knot span")
      ->capture_default_str()
      ->check(axom::CLI::PositiveNumber);

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    app.get_formatter()->column_width(48);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug : slic::message::Info);
  }
};

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);

  axom::slic::SimpleLogger logger;
  // slic::debug::checksAreErrors = true;

  // Set up and parse command line arguments
  Input params;
  axom::CLI::App app {"Driver for In/Out surface containment query"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(params.annotationMode);

  AXOM_ANNOTATE_BEGIN("quest containment example");

  const bool is2D = params.isInput2D();

  ContainmentDriver<2> driver2D;
  ContainmentDriver<3> driver3D;

  /// Load mesh file
  SLIC_INFO(axom::fmt::format("{:-^80}", " Loading the mesh "));
  SLIC_INFO(axom::fmt::format("Reading file: '{}'...", params.inputFile));

  AXOM_ANNOTATE_METADATA("dimension", is2D ? 2 : 3, "");
  AXOM_ANNOTATE_BEGIN("init");

  if(is2D)
  {
    driver2D.loadContourMesh(params.inputFile, params.samplesPerKnotSpan);
  }
  else
  {
    driver3D.loadSTLMesh(params.inputFile);
  }

  /// Compute mesh bounding box and log some stats about the surface
  if(is2D)
  {
    driver2D.computeBounds();
    driver2D.printSurfaceStats();
  }
  else
  {
    driver3D.computeBounds();
    driver3D.printSurfaceStats();
  }

  /// Create octree over mesh's bounding box
  SLIC_INFO(axom::fmt::format("{:-^80}", " Generating the octree "));
  if(is2D)
  {
    driver2D.initializeInOutOctree();
    driver2D.printSurfaceStats();

    mint::write_vtk(driver2D.getSurfaceMesh(), "meldedSegmentMesh.vtk");
  }
  else
  {
    driver3D.initializeInOutOctree();
    driver3D.printSurfaceStats();

    mint::write_vtk(driver3D.getSurfaceMesh(), "meldedTriMesh.vtk");
  }

  AXOM_ANNOTATE_END("init");
  AXOM_ANNOTATE_BEGIN("query");

  /// Query the octree over mesh's bounding box
  SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the octree "));
  if(is2D)
  {
    driver2D.initializeQueryBox(params.queryBoxMins, params.queryBoxMaxs);

    // Query the mesh
    for(int i = 1; i < params.maxQueryLevel; ++i)
    {
      const int res = 1 << i;
      driver2D.testContainmentOnRegularGrid(res, params.useBatchedQuery(), params.isVerbose());
    }
  }
  else
  {
    driver3D.initializeQueryBox(params.queryBoxMins, params.queryBoxMaxs);

    // Query the mesh
    for(int i = 1; i < params.maxQueryLevel; ++i)
    {
      const int res = 1 << i;
      driver3D.testContainmentOnRegularGrid(res, params.useBatchedQuery(), params.isVerbose());
    }
  }

  AXOM_ANNOTATE_END("query");
  SLIC_INFO(axom::fmt::format("{:-^80}", ""));
  axom::slic::flushStreams();

  AXOM_ANNOTATE_END("quest containment example");

  return 0;
}
