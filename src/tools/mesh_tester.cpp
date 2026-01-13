// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/mint/mesh/FieldVariable.hpp"
// _read_stl_include2_start
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
// _read_stl_include2_end
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/utils/vtk_utils.hpp"

#include "axom/primal.hpp"
#include "axom/spin/UniformGrid.hpp"

// _read_stl_include1_start
#include "axom/quest/io/STLReader.hpp"
#include "axom/quest/io/STLWriter.hpp"
// _read_stl_include1_end
// _check_repair_include_start
#include "axom/quest/MeshTester.hpp"
// _check_repair_include_end

#include "axom/fmt.hpp"

#include "axom/core/execution/nested_for_exec.hpp"

using seq_exec = axom::SEQ_EXEC;

// clang-format off
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  using omp_exec = axom::OMP_EXEC;
#else
  using omp_exec = seq_exec;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  constexpr int CUDA_BLOCK_SIZE = 256;
  using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
#else
  using cuda_exec = seq_exec;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  constexpr int HIP_BLOCK_SIZE = 256;
  using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
#else
  using hip_exec = seq_exec;
#endif
// clang-format on

#include "axom/CLI11.hpp"

// C/C++ includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include <iomanip>
#include <memory>

// _read_stl_typedefs_start
using namespace axom;

using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
// _read_stl_typedefs_end
using Triangle3 = primal::Triangle<double, 3>;

using Point3 = primal::Point<double, 3>;
using SpatialBoundingBox = primal::BoundingBox<double, 3>;
using UniformGrid3 = spin::UniformGrid<int, 3>;
using Vector3 = primal::Vector<double, 3>;
using Segment3 = primal::Segment<double, 3>;

enum RuntimePolicy
{
  seq = 0,
  raja_seq = 1,
  raja_omp = 2,
  raja_cuda = 3,
  raja_hip = 4
};

template <>
struct axom::fmt::formatter<RuntimePolicy> : axom::fmt::formatter<std::string>
{
  template <typename FormatContext>
  auto format(const RuntimePolicy& policy, FormatContext& ctx) const
  {
    std::string name = "unknown";
    switch(policy)
    {
    case seq:
      name = "seq";
      break;
    case raja_seq:
      name = "raja_seq";
      break;
    case raja_omp:
      name = "raja_omp";
      break;
    case raja_cuda:
      name = "raja_cuda";
      break;
    case raja_hip:
      name = "raja_hip";
      break;
    default:
      name = "unknown";
      break;
    }
    return axom::fmt::formatter<std::string>::format(name, ctx);
  }
};

struct Input
{
  static const std::set<std::string> s_validMethods;
  static const std::set<std::string> s_validFormats;
  static const std::map<std::string, RuntimePolicy> s_validPolicies;

  std::string stlInput {""};
  std::string fileOutput {""};
  std::string fileFormat {"vtk"};
  std::string method {"uniform"};
  RuntimePolicy policy {seq};
  std::string annotationMode {"none"};

  int resolution {0};
  double weldThreshold {1e-6};
  double intersectionThreshold {1e-08};
  bool skipWeld {false};
  bool verboseOutput {false};
  bool binary {false};

  Input() = default;

  void parse(int argc, char** argv, axom::CLI::App& app);

  std::string collisionsMeshName() { return fileOutput + ".collisions." + fileFormat; }
  std::string collisionsTextName() { return fileOutput + ".collisions.txt"; }
  std::string weldMeshName() { return fileOutput + ".weld." + fileFormat; }

private:
  void fixOutfilePath();
};

// clang-format off
const std::set<std::string> Input::s_validMethods({
  "bvh",
  "uniform",
  "implicit",
  "naive"
});

const std::set<std::string> Input::s_validFormats({
  "stl",
  "vtk"
});

const std::map<std::string, RuntimePolicy> Input::s_validPolicies({
  {"seq", seq}
  #ifdef AXOM_USE_RAJA
  , {"raja_seq", raja_seq}
    #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  , {"raja_omp", raja_omp}
    #endif
    #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  , {"raja_cuda", raja_cuda}
    #endif
    #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  , {"raja_hip", raja_hip}
    #endif
  #endif
});
// clang-format on

void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  app
    .add_option("-m, --method",
                method,
                "Method to use. \n"
                "Set to 'bvh' to use the bounding volume hierarchy spatial index.\n"
                "Set to 'naive' to use the naive algorithm (without a spatial index).\n"
                "Set to 'uniform' to use the uniform grid spatial index.\n"
                "Set to 'implicit' to use the implicit grid spatial index.")
    ->capture_default_str()
    ->check(axom::CLI::IsMember {Input::s_validMethods});

  app
    .add_option("-r,--resolution",
                resolution,
                "With '-m uniform', set resolution of uniform grid. \n"
                "Set to less than 1 to use the uniform grid spatial index\n"
                "with a resolution of the cube root of the number of\n"
                "triangles.")
    ->capture_default_str();

  std::stringstream pol_sstr;
  pol_sstr << "With '-m bvh' or '-m naive', set runtime policy. \n"
           << "Set to 'seq' or 0 to use the sequential algorithm "
           << "(w/o RAJA).";
#ifdef AXOM_USE_RAJA
  pol_sstr << "\nSet to 'raja_seq' or 1 to use the RAJA sequential policy.";
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  pol_sstr << "\nSet to 'raja_omp' or 2 to use the RAJA OpenMP policy.";
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  pol_sstr << "\nSet to 'raja_cuda' or 3 to use the RAJA CUDA policy.";
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  pol_sstr << "\nSet to 'raja_hip' or 4 to use the RAJA HIP policy.";
  #endif
#endif

  app.add_option("-p, --policy", policy, pol_sstr.str())
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(Input::s_validPolicies));

  app.add_option("-f, --format", fileFormat, "Output file format, one of: 'stl', 'vtk'")
    ->capture_default_str()
    ->check(axom::CLI::IsMember {Input::s_validFormats});

  app.add_option("-i,--infile", stlInput, "The STL input file")->required()->check(axom::CLI::ExistingFile);

  app.add_option("-o,--outfile",
                 fileOutput,
                 "Output file name for collisions and welded mesh. \n"
                 "Defaults to a file in the CWD with the input file name. \n"
                 "Collisions mesh will end with '.collisions.{format}' and \n"
                 "welded mesh will end with '.welded.{format}', where {format} is\n"
                 "the selected file format.");

  app.add_option("--weldThresh", weldThreshold, "Distance threshold for welding vertices.")
    ->capture_default_str();

  app
    .add_option("--intersectionThresh",
                intersectionThreshold,
                "Tolerance threshold to use when testing for intersecting triangles")
    ->capture_default_str();

  app.add_flag("--skipWeld",
               skipWeld,
               "Don't weld vertices (useful for testing, not helpful otherwise).");

  app.add_flag("--binary", binary, "Write binary output files in supported formats.");

  app.add_flag("-v,--verbose", verboseOutput, "Increase logging verbosity.")->capture_default_str();

#ifdef AXOM_USE_CALIPER
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

  app.get_formatter()->column_width(76);

  // Could throw an exception
  app.parse(argc, argv);

  // Fix any options that need fixing
  fixOutfilePath();

  // clang-format off
  const std::string method_str =
    axom::fmt::format("\n  method = {} {}", method, [this]() -> std::string {
      if(this->method == "uniform")    { return axom::fmt::format(" (use uniform grid)\n  resolution = {}", this->resolution); }
      else if(this->method == "naive") { return " (use naive algorithm)"; }
      else if(this->method == "bvh")   { return " (use bounding volume hierarchy)"; }
      else                             { return "";}
    }());

  const std::string policy_str = (method == "naive" || method == "bvh")
    ? axom::fmt::format("\n  policy = {} {}", static_cast<int>(policy),
      [this]() -> std::string {
        switch(this->policy)
        {
        case seq:       return " (use sequential policy)";
        case raja_omp:  return " (use RAJA OpenMP policy)";
        case raja_cuda: return " (use RAJA CUDA policy)";
        case raja_hip:  return " (use RAJA HIP policy)";
        case raja_seq:  return " (use RAJA sequential policy)";
        default:        return "";
        }
      }())
    : "";
  // clang-format on

#ifdef AXOM_USE_CALIPER
  const std::string annotation_str = axom::fmt::format("\n  annotation mode = '{}'", annotationMode);
#else
  const std::string annotation_str {};
#endif

  // Output parsed information
  SLIC_INFO("Using parameter values: "  //
            << method_str               //
            << policy_str               //
            << annotation_str
            << (skipWeld ? "\n  skipping weld"
                         : axom::fmt::format("\n  weld threshold = {}", weldThreshold))
            << axom::fmt::format("\n  intersection tolerance = {}", intersectionThreshold)
            << axom::fmt::format("\n  infile = '{}'", stlInput)
            << axom::fmt::format("\n  collisions outfile = '{}'", collisionsMeshName())
            << axom::fmt::format("\n  weld outfile = '{}'", weldMeshName()));
}

void Input::fixOutfilePath()
{
  namespace futil = axom::utilities::filesystem;

  // Extract the stem of the input file, to output files in the CWD
  std::string inFileDir;
  futil::getDirName(inFileDir, stlInput);
  int sepSkip = (int)(inFileDir.size() > 0);
  std::string inFileStem = stlInput.substr(inFileDir.size() + sepSkip);
  std::string outFileBase = futil::joinPath(futil::getCWD(), inFileStem);

  // set output file name when not provided
  int sz = static_cast<int>(fileOutput.size());
  if(sz < 1)
  {
    fileOutput = outFileBase;
    sz = static_cast<int>(fileOutput.size());
  }

  // ensure that output file does not end with '.vtk'
  if(sz > 4)
  {
    std::string ext = fileOutput.substr(sz - 4, 4);
    if(ext == ".vtk" || ext == ".stl")
    {
      fileOutput = fileOutput.substr(0, sz - ext.size());
    }
  }
}

inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS);

AXOM_HOST_DEVICE bool checkTT(Triangle3& t1, Triangle3& t2, double EPS);

std::vector<std::pair<int, int>> naiveIntersectionAlgorithm(mint::Mesh* surface_mesh,
                                                            std::vector<int>& degenerate,
                                                            double EPS);

void announceMeshProblems(int triangleCount, int intersectPairCount, int degenerateCount);

void saveProblemFlagsToMesh(mint::Mesh* surface_mesh,
                            const std::vector<std::pair<int, int>>& c,
                            const std::vector<int>& d);

bool writeAnnotatedMesh(mint::Mesh* surface_mesh,
                        const std::string& outfile,
                        const std::string& fileFormat,
                        bool binary);

bool writeCollisions(const std::vector<std::pair<int, int>>& c,
                     const std::vector<int>& d,
                     std::string basename);

inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS = 1.0e-9)
{
  return axom::utilities::isNearlyEqual(p1[0], p2[0], EPS) &&
    axom::utilities::isNearlyEqual(p1[1], p2[1], EPS) &&
    axom::utilities::isNearlyEqual(p1[2], p2[2], EPS);
}

AXOM_HOST_DEVICE
bool checkTT(Triangle3& t1, Triangle3& t2, double EPS)
{
  if(t2.degenerate())
  {
    return false;
  }

  constexpr bool includeBoundaries = false;  // only check internal intersections
  if(primal::intersect(t1, t2, includeBoundaries, EPS))
  {
    return true;
  }
  return false;
}

inline Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getCellType(i) == mint::TRIANGLE);
  primal::Point<axom::IndexType, 3> triCell;
  Triangle3 tri;
  surface_mesh->getCellNodeIDs(i, triCell.data());

  surface_mesh->getNode(triCell[0], tri[0].data());
  surface_mesh->getNode(triCell[1], tri[1].data());
  surface_mesh->getNode(triCell[2], tri[2].data());

  return tri;
}

std::vector<std::pair<int, int>> naiveIntersectionAlgorithm(mint::Mesh* surface_mesh,
                                                            std::vector<int>& degenerate,
                                                            double EPS)
{
  SLIC_INFO("Running naive intersection algorithm.");
  AXOM_ANNOTATE_SCOPE("weldTriMeshVertices");

  // For each triangle, check for intersection against
  // every other triangle with a greater index in the mesh, excluding
  // degenerate triangles.
  std::vector<std::pair<int, int>> retval;

  const int ncells = surface_mesh->getNumberOfCells();
  SLIC_INFO(
    axom::fmt::format(axom::utilities::locale(), "Checking mesh with a total of {:L} cells.", ncells));

  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();

  // For each triangle in the mesh
  for(int i = 0; i < ncells; i++)
  {
    t1 = getMeshTriangle(i, surface_mesh);

    // Skip if degenerate
    if(t1.degenerate())
    {
      degenerate.push_back(i);
      continue;
    }

    // If the triangle is not degenerate, test against all other
    // triangles that this triangle has not been checked against
    for(int j = i + 1; j < ncells; j++)
    {
      t2 = getMeshTriangle(j, surface_mesh);
      if(checkTT(t1, t2, EPS))
      {
        retval.push_back(std::make_pair(i, j));
      }
    }
  }

  return retval;
}

#if defined(AXOM_USE_RAJA)
template <typename ExecSpace>
std::vector<std::pair<int, int>> naiveIntersectionAlgorithm(mint::Mesh* surface_mesh,
                                                            std::vector<int>& degenerate,
                                                            double EPS)
{
  SLIC_INFO("Running naive intersection algorithm in execution Space: "
            << axom::execution_space<ExecSpace>::name());
  AXOM_ANNOTATE_SCOPE("naiveIntersectionAlgorithm");

  // Get allocators
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecSpace>::allocatorID();

  std::vector<std::pair<int, int>> retval;

  const int ncells = surface_mesh->getNumberOfCells();
  SLIC_INFO(
    axom::fmt::format(axom::utilities::locale(), "Checking mesh with a total of {:L} cells.", ncells));

  axom::Array<Triangle3> tris_h(ncells, ncells, host_allocator);

  // Get each triangle in the mesh and check for degeneracies
  for(int i = 0; i < ncells; i++)
  {
    tris_h[i] = getMeshTriangle(i, surface_mesh);
    if(tris_h[i].degenerate())
    {
      degenerate.push_back(i);
    }
  }

  axom::StackArray<axom::IndexType, 2> row_range {{0, static_cast<axom::IndexType>(ncells)}};
  axom::StackArray<axom::IndexType, 2> col_range {{0, static_cast<axom::IndexType>(ncells)}};

  axom::ReduceSum<ExecSpace, int> numIntersect(0);

  // Copy triangles to device
  axom::Array<Triangle3> tris_d =
    on_device ? axom::Array<Triangle3>(tris_h, device_allocator) : axom::Array<Triangle3>();

  auto tris_v = on_device ? tris_d.view() : tris_h.view();

  // Compute the number of intersections
  axom::for_all<ExecSpace>(
    col_range,
    row_range,
    AXOM_LAMBDA(int col, int row) {
      if(row > col)
      {
        if(checkTT(tris_v[row], tris_v[col], EPS))
        {
          numIntersect += 1;
        }
      }
    });

  // Allocation to hold intersection pairs and counter to know where to store
  axom::Array<int> intersections_d(numIntersect.get() * 2, numIntersect.get() * 2, device_allocator);
  axom::Array<int> counter_h(1, 1, host_allocator);
  counter_h[0] = 0;
  axom::Array<int> counter_d =
    on_device ? axom::Array<int>(counter_h, device_allocator) : axom::Array<int>();

  auto intersections_v = intersections_d.view();
  auto counter_v = on_device ? counter_d.view() : counter_h.view();

  // RAJA loop to populate with intersections
  axom::for_all<ExecSpace>(
    col_range,
    row_range,
    AXOM_LAMBDA(int col, int row) {
      if(row > col)
      {
        if(checkTT(tris_v[row], tris_v[col], EPS))
        {
          auto idx = axom::atomicAdd<ExecSpace>(counter_v.data(), 2);
          intersections_v[idx] = row;
          intersections_v[idx + 1] = col;
        }
      }
    });

  // Copy intersections to host
  axom::Array<int> intersections_h =
    on_device ? axom::Array<int>(intersections_d, host_allocator) : std::move(intersections_d);

  // Initialize pairs of clashes
  for(auto i = 0; i < numIntersect.get() * 2; i += 2)
  {
    retval.push_back(std::make_pair(intersections_h[i], intersections_h[i + 1]));
  }

  return retval;
}
#endif

void announceMeshProblems(int triangleCount, int intersectPairCount, int degenerateCount)
{
  SLIC_INFO(
    axom::fmt::format(axom::utilities::locale(),
                      "{:L} triangles, with {:L} intersecting tri pairs {:L} degenerate tris.",
                      triangleCount,
                      intersectPairCount,
                      degenerateCount));
}

void saveProblemFlagsToMesh(mint::Mesh* mesh,
                            const std::vector<std::pair<int, int>>& c,
                            const std::vector<int>& d)
{
  // Create new Field variables to hold degenerate and intersecting info
  const int num_cells = mesh->getNumberOfCells();
  int* intersectptr = mesh->createField<int>("nbr_intersection", mint::CELL_CENTERED);
  int* dgnptr = mesh->createField<int>("degenerate_triangles", mint::CELL_CENTERED);

  // Initialize everything to 0
  for(int i = 0; i < num_cells; ++i)
  {
    intersectptr[i] = 0;
    dgnptr[i] = 0;
  }

  // Fill in intersect flag
  size_t csize = c.size();
  for(size_t i = 0; i < csize; ++i)
  {
    std::pair<int, int> theC = c[i];
    intersectptr[theC.first] += 1;
    intersectptr[theC.second] += 1;
  }

  // Fill in degenerate flag
  size_t dsize = d.size();
  for(size_t i = 0; i < dsize; ++i)
  {
    dgnptr[d[i]] = 1;
  }
}

bool writeAnnotatedMesh(mint::Mesh* surface_mesh,
                        const std::string& outfile,
                        const std::string& fileFormat,
                        bool binary)
{
  bool retval = false;
  if(fileFormat == "vtk")
  {
    AXOM_ANNOTATE_SCOPE("save vtk mesh");
    retval = mint::write_vtk(surface_mesh, outfile) == 0;
  }
  else if(fileFormat == "stl")
  {
    AXOM_ANNOTATE_SCOPE("save stl mesh");
    // NOTE: The STL writer is co-located with the STL reader in quest.
    retval = quest::write_stl(surface_mesh, outfile, binary) == 0;
  }
  return retval;
}

bool writeCollisions(const std::vector<std::pair<int, int>>& c,
                     const std::vector<int>& d,
                     std::string filename)
{
  std::ofstream outf(filename.c_str());
  if(!outf)
  {
    return false;
  }

  outf << c.size() << " intersecting triangle pairs:" << std::endl;
  for(size_t i = 0; i < c.size(); ++i)
  {
    outf << c[i].first << " " << c[i].second << std::endl;
  }

  outf << d.size() << " degenerate triangles:" << std::endl;
  for(size_t i = 0; i < d.size(); ++i)
  {
    outf << d[i] << std::endl;
  }

  return true;
}

/// Simple RAII struct to set up and tear down slic-based logger
struct MeshTesterRAIILogger
{
  MeshTesterRAIILogger()
  {
    // Initialize the SLIC logger
    slic::initialize();
    slic::setLoggingMsgLevel(axom::slic::message::Debug);

    // Customize logging levels and formatting
    std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";
    slic::GenericOutputStream* defaultStream = new slic::GenericOutputStream(&std::cout);
    slic::GenericOutputStream* compactStream =
      new slic::GenericOutputStream(&std::cout, slicFormatStr);

    slic::addStreamToMsgLevel(defaultStream, axom::slic::message::Error);
    slic::addStreamToMsgLevel(compactStream, axom::slic::message::Warning);
    slic::addStreamToMsgLevel(compactStream, axom::slic::message::Info);
    slic::addStreamToMsgLevel(compactStream, axom::slic::message::Debug);
  }

  ~MeshTesterRAIILogger() { slic::finalize(); }
};

/*!
 * The mesh tester checks a triangulated surface mesh for several problems.
 *
 * Currently the mesh tester checks for intersecting triangles.  This is
 * implemented in three ways.  First, a naive algorithm tests each triangle
 * against each other triangle.  This is easy to understand and verify,
 * but slow. Second, the same naive algorithm is run using raja. A third
 * algorithm uses a UniformGrid to index the bounding box of the mesh,
 * and checks each triangle for intersection with the triangles in all the
 * bins the triangle's bounding box falls into.
 *
 * Currently, the mesh tester works only with Triangle meshes.
 */

int main(int argc, char** argv)
{
  MeshTesterRAIILogger raii_logger;

  SLIC_INFO(axom::fmt::format("Axom Version: [{}]", axom::getVersion()));

  // Parse the command line arguments
  Input params;
  axom::CLI::App app {"MeshTester example"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(params.annotationMode);

  AXOM_ANNOTATE_SCOPE("mesh_tester");

  // _read_stl_file_start
  // Read file
  UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
  {
    SLIC_INFO("Reading file: '" << params.stlInput << "'...\n");
    AXOM_ANNOTATE_SCOPE("load stl mesh");
    auto reader = std::make_unique<quest::STLReader>();
    reader->setFileName(params.stlInput);
    reader->read();

    reader->getMesh(surface_mesh);
  }

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Mesh has {:L} vertices and {:L} triangles",
                              surface_mesh->getNumberOfNodes(),
                              surface_mesh->getNumberOfCells()));
  // _read_stl_file_end

  // Vertex welding
  if(!params.skipWeld)
  {
    AXOM_ANNOTATE_SCOPE("weld triangles");
    axom::utilities::Timer timer(true);

    // _check_repair_weld_start
    quest::weldTriMeshVertices(&surface_mesh, params.weldThreshold);
    // _check_repair_weld_end

    timer.stop();
    SLIC_INFO("Vertex welding took " << timer.elapsedTimeInSec() << " seconds.");
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "After welding, mesh has {:L} vertices and {:L} triangles",
                                surface_mesh->getNumberOfNodes(),
                                surface_mesh->getNumberOfCells()));
  }

  // Detect collisions
  {
    AXOM_ANNOTATE_SCOPE("detect collisions");

    // _check_repair_intersections_containers_start
    std::vector<std::pair<int, int>> collisions;
    std::vector<int> degenerate;
    // _check_repair_intersections_containers_end

    axom::utilities::Timer timer(true);
    if(params.method == "naive")
    {
      switch(params.policy)
      {
      case seq:
        collisions =
          naiveIntersectionAlgorithm(surface_mesh, degenerate, params.intersectionThreshold);
        break;
#if defined(AXOM_USE_RAJA)
      case raja_seq:
        collisions =
          naiveIntersectionAlgorithm<seq_exec>(surface_mesh, degenerate, params.intersectionThreshold);
        break;
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case raja_omp:
        collisions =
          naiveIntersectionAlgorithm<omp_exec>(surface_mesh, degenerate, params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case raja_cuda:
        collisions = naiveIntersectionAlgorithm<cuda_exec>(surface_mesh,
                                                           degenerate,
                                                           params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case raja_hip:
        collisions =
          naiveIntersectionAlgorithm<hip_exec>(surface_mesh, degenerate, params.intersectionThreshold);
        break;
  #endif
#endif  // AXOM_USE_RAJA && AXOM_USE_UMPIRE

      default:
        SLIC_ERROR("Unhandled runtime policy case " << params.policy);
        break;
      }
    }  // end of if method == 'naive'
    else if(params.method == "bvh")
    {
      switch(params.policy)
      {
      case seq:
#ifdef AXOM_USE_RAJA
        SLIC_INFO("BVH was compiled with RAJA - seq and raja_seq execution will be equivalent");
#endif
        quest::findTriMeshIntersectionsBVH<seq_exec, double>(surface_mesh,
                                                             collisions,
                                                             degenerate,
                                                             params.intersectionThreshold);
        break;
#if defined(AXOM_USE_RAJA)
      case raja_seq:
        quest::findTriMeshIntersectionsBVH<seq_exec, double>(surface_mesh,
                                                             collisions,
                                                             degenerate,
                                                             params.intersectionThreshold);
        break;
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case raja_omp:
        quest::findTriMeshIntersectionsBVH<omp_exec, double>(surface_mesh,
                                                             collisions,
                                                             degenerate,
                                                             params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case raja_cuda:
        quest::findTriMeshIntersectionsBVH<cuda_exec, double>(surface_mesh,
                                                              collisions,
                                                              degenerate,
                                                              params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case raja_hip:
        quest::findTriMeshIntersectionsBVH<hip_exec, double>(surface_mesh,
                                                             collisions,
                                                             degenerate,
                                                             params.intersectionThreshold);
        break;
  #endif
#endif  // AXOM_USE_RAJA && AXOM_USE_UMPIRE
      default:
        SLIC_ERROR("Unhandled runtime policy case " << params.policy);
        break;
      }
    }  // end of if method == 'bvh'
    else if(params.method == "implicit")
    {
      switch(params.policy)
      {
      case seq:
#ifdef AXOM_USE_RAJA
        SLIC_INFO(
          "ImplicitGrid was compiled with RAJA - seq and raja_seq execution will be equivalent");
#endif
        quest::findTriMeshIntersectionsImplicitGrid<seq_exec, double>(surface_mesh,
                                                                      collisions,
                                                                      degenerate,
                                                                      params.resolution,
                                                                      params.intersectionThreshold);
        break;
#ifdef AXOM_USE_RAJA
      case raja_seq:
        quest::findTriMeshIntersectionsImplicitGrid<seq_exec, double>(surface_mesh,
                                                                      collisions,
                                                                      degenerate,
                                                                      params.resolution,
                                                                      params.intersectionThreshold);
        break;
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case raja_omp:
        quest::findTriMeshIntersectionsImplicitGrid<omp_exec, double>(surface_mesh,
                                                                      collisions,
                                                                      degenerate,
                                                                      params.resolution,
                                                                      params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case raja_cuda:
        quest::findTriMeshIntersectionsImplicitGrid<cuda_exec, double>(surface_mesh,
                                                                       collisions,
                                                                       degenerate,
                                                                       params.resolution,
                                                                       params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case raja_hip:
        quest::findTriMeshIntersectionsImplicitGrid<hip_exec, double>(surface_mesh,
                                                                      collisions,
                                                                      degenerate,
                                                                      params.resolution,
                                                                      params.intersectionThreshold);
        break;
  #endif
#endif  // AXOM_USE_RAJA
      default:
        SLIC_ERROR("Unhandled runtime policy case " << params.policy);
        break;
      }
    }
    else
    {
      switch(params.policy)
      {
      case seq:
        // _check_repair_intersections_start
        // Use a uniform grid spatial index
        quest::findTriMeshIntersections(surface_mesh,
                                        collisions,
                                        degenerate,
                                        params.resolution,
                                        params.intersectionThreshold);
        // _check_repair_intersections_end
        break;
#ifdef AXOM_USE_RAJA
      case raja_seq:
        quest::findTriMeshIntersectionsUniformGrid<seq_exec, double>(surface_mesh,
                                                                     collisions,
                                                                     degenerate,
                                                                     params.resolution,
                                                                     params.intersectionThreshold);
        break;
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case raja_omp:
        quest::findTriMeshIntersectionsUniformGrid<omp_exec, double>(surface_mesh,
                                                                     collisions,
                                                                     degenerate,
                                                                     params.resolution,
                                                                     params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case raja_cuda:
        quest::findTriMeshIntersectionsUniformGrid<cuda_exec, double>(surface_mesh,
                                                                      collisions,
                                                                      degenerate,
                                                                      params.resolution,
                                                                      params.intersectionThreshold);
        break;
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case raja_hip:
        quest::findTriMeshIntersectionsUniformGrid<hip_exec, double>(surface_mesh,
                                                                     collisions,
                                                                     degenerate,
                                                                     params.resolution,
                                                                     params.intersectionThreshold);
        break;
  #endif
#endif  // AXOM_USE_RAJA
      default:
        SLIC_ERROR("Unhandled runtime policy case " << params.policy);
        break;
      }
    }
    timer.stop();
    SLIC_INFO("Detecting intersecting triangles took " << timer.elapsedTimeInSec() << " seconds.");

    {
      AXOM_ANNOTATE_SCOPE("print collisions");
      announceMeshProblems(surface_mesh->getNumberOfCells(),
                           static_cast<int>(collisions.size()),
                           static_cast<int>(degenerate.size()));

      saveProblemFlagsToMesh(surface_mesh, collisions, degenerate);

      if(!writeAnnotatedMesh(surface_mesh, params.collisionsMeshName(), params.fileFormat, params.binary))
      {
        SLIC_ERROR("Couldn't write results to " << params.collisionsMeshName());
      }

      if(!writeCollisions(collisions, degenerate, params.collisionsTextName()))
      {
        SLIC_ERROR("Couldn't write results to " << params.collisionsTextName());
      }

      if(params.verboseOutput && !collisions.empty())
      {
        SLIC_INFO("Intersecting triangle pairs:");
        // Initialize pairs of clashes
        for(auto i : collisions)
        {
          auto t1 = getMeshTriangle(i.first, surface_mesh);
          auto t2 = getMeshTriangle(i.second, surface_mesh);

          SLIC_INFO("  Triangle " << i.first << " -- " << std::setprecision(17) << t1);
          SLIC_INFO("  Triangle " << i.second << " -- " << std::setprecision(17) << t2 << "\n");
        }
      }
    }
  }

  // Look for holes---must be welded
  if(params.skipWeld)
  {
    SLIC_INFO("Watertight check depends on vertex welding, which was skipped.");
  }
  else
  {
    AXOM_ANNOTATE_SCOPE("check watertight");

    SLIC_INFO("Checking for watertight mesh.");
    axom::utilities::Timer timer2(true);
    // _check_watertight_start
    quest::WatertightStatus wtstat = quest::isSurfaceMeshWatertight(surface_mesh);
    // _check_watertight_end
    timer2.stop();
    // _report_watertight_start
    switch(wtstat)
    {
    case quest::WatertightStatus::WATERTIGHT:
      SLIC_INFO("  The mesh is watertight.");
      break;
    case quest::WatertightStatus::NOT_WATERTIGHT:
      SLIC_INFO("  The mesh is not watertight: at least one boundary edge was detected.");
      break;
    default:
      SLIC_INFO(
        "  An error was encountered while checking. This may be due to a non-manifold mesh.");
      break;
    }
    // _report_watertight_end
    SLIC_INFO("Testing for watertightness took " << timer2.elapsedTimeInSec() << " seconds.");

    writeAnnotatedMesh(surface_mesh, params.weldMeshName(), params.fileFormat, params.binary);
  }

  // Delete the mesh
  delete surface_mesh;
  surface_mesh = nullptr;

  return 0;
}
