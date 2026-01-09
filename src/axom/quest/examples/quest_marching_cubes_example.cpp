// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 \file marching_cubes_example.cpp
 \brief Driver and test for a marching cubes isocontour generation

  The test can generate planar and round contours.  Planar contours
  can be checked to machine-zero accuracy, but it doesn't test a great
  variety of contour-mesh intersection types.  Round contours can
  check more intersection types but requires a tolerance to allow
  for the function not varying linearly along mesh lines.
*/

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesFullParallel.hpp requires conduit"
#endif

// Axom includes
#include "axom/core.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/core/MDMapping.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/MeshViewUtil.hpp"
#if defined(AXOM_USE_SIDRE)
  #include "axom/sidre.hpp"
#endif
#include "axom/core/Types.hpp"
#include "axom/core/numerics/floating_point_limits.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_relay_io_blueprint.hpp"
#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
  #include "conduit_relay_mpi_io_blueprint.hpp"
#endif

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// C/C++ includes
#include <string>
#include <map>
#include <vector>
#include <cmath>

namespace quest = axom::quest;
namespace slic = axom::slic;
#if defined(AXOM_USE_SIDRE)
namespace sidre = axom::sidre;
#endif
namespace primal = axom::primal;
namespace mint = axom::mint;
namespace numerics = axom::numerics;

using RuntimePolicy = axom::runtime_policy::Policy;

///////////////////////////////////////////////////////////////
// converts the input string into an 80 character string
// padded on both sides with '=' symbols
std::string banner(const std::string& str) { return axom::fmt::format("{:=^80}", str); }

///////////////////////////////////////////////////////////////
/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;
  std::string fieldsFile {"fields"};

  // Center of round contour function
  std::vector<double> fcnCenter;

  // Scaling factor for gyroid function
  std::vector<double> gyroidScale;

  // Parameters for planar contour function
  std::vector<double> inPlane;
  std::vector<double> perpDir;

  std::size_t ndim {0};

  double contourVal {1.0};

  bool checkResults {false};

  RuntimePolicy policy {RuntimePolicy::seq};

  quest::MarchingCubesDataParallelism dataParallelism = quest::MarchingCubesDataParallelism::byPolicy;

  // Distinct MarchingCubes objects count.
  int objectRepCount = 1;
  // Contour generation count for each MarchingCubes objects.
  int contourGenCount = 1;
  // Number of masking cycles.
  int maskCount = 1;

  std::string annotationMode {"none"};

private:
  bool _verboseOutput {false};

  // clang-format off
  const std::map<std::string, quest::MarchingCubesDataParallelism> s_validImplChoices
  {
    {"byPolicy", quest::MarchingCubesDataParallelism::byPolicy}
    , {"hybridParallel", quest::MarchingCubesDataParallelism::hybridParallel}
    , {"fullParallel", quest::MarchingCubesDataParallelism::fullParallel}
  };
  // clang-format on

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-p, --policy", policy)
      ->description("Set runtime policy for point query method")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.add_option("--dataParallelism", dataParallelism)
      ->description("Set full or partial data-parallelism, or by-policy")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(s_validImplChoices));

    app.add_option("-m,--mesh-file", meshFile)
      ->description(
        "Path to multidomain computational mesh following conduit blueprint "
        "convention.")
      ->check(axom::CLI::ExistingFile);

    app.add_option("-s,--fields-file", fieldsFile)
      ->description("Name of output mesh file with all its fields.")
      ->capture_default_str();

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    auto* distanceFunctionOption =
      app.add_option_group("distanceFunctionOption", "Options for specifying a distance function.");

    auto* distFromPtOption = distanceFunctionOption->add_option_group(
      "distFromPtOption",
      "Options for setting up distance-from-point function");
    distFromPtOption->add_option("--center", fcnCenter)
      ->description("Center for distance-from-point function (x,y[,z])")
      ->expected(2, 3);

    auto* gyroidOption =
      distanceFunctionOption->add_option_group("gyroidOption",
                                               "Options for setting up gyroid function");
    gyroidOption->add_option("--scale", gyroidScale)
      ->description("Scaling factor for gyroid function (x,y[,z])")
      ->expected(2, 3);

    auto* distFromPlaneOption = distanceFunctionOption->add_option_group(
      "distFromPlaneOption",
      "Options for setting up distance-from-plane function");
    auto* perpDirOption =
      distFromPlaneOption->add_option("--dir", perpDir)
        ->description("Positive direction for distance-from-plane function (x,y[,z])")
        ->expected(2, 3);
    distFromPlaneOption->add_option("--inPlane", inPlane)
      ->description("In-plane point for distance-from-plane function (x,y[,z])")
      ->expected(2, 3)
      ->needs(perpDirOption);

    // Require at least one distance function, and allow all three.
    distanceFunctionOption->require_option(1, 3);

    app.add_option("--contourVal", contourVal)->description("Contour value")->capture_default_str();

    app.add_flag("-c,--check-results,!--no-check-results", checkResults)
      ->description("Enable/disable checking results against analytical solution")
      ->capture_default_str();

    //
    // Number of repetitions to run
    //

    app.add_option("--objectReps", objectRepCount)
      ->description("Number of MarchingCube object repetitions to run")
      ->capture_default_str();

    app.add_option("--contourGenReps", contourGenCount)
      ->description("Number of contour repetitions to run for each MarchingCubes object")
      ->capture_default_str();

    app.add_option("--maskCount", maskCount)
      ->description(
        "Group the cells using this many masking groups to test masking "
        "(default to 1).")
      ->capture_default_str()
      ->check(axom::CLI::Range(1, std::numeric_limits<int>::max()));

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(_verboseOutput ? slic::message::Debug : slic::message::Info);

    ndim = std::max({ndim, fcnCenter.size(), inPlane.size(), perpDir.size(), gyroidScale.size()});
    SLIC_ASSERT_MSG((fcnCenter.empty() || fcnCenter.size() == ndim) &&
                      (inPlane.empty() || inPlane.size() == ndim) &&
                      (perpDir.empty() || perpDir.size() == ndim) &&
                      (gyroidScale.empty() || gyroidScale.size() == ndim),
                    "fcnCenter, inPlane and perpDir must have consistent sizes "
                    "if specified.");

    // inPlane defaults to origin if omitted.
    if(usingPlanar() && inPlane.empty())
    {
      inPlane.insert(inPlane.begin(), ndim, 0.0);
    }
  }

  bool usingPlanar() { return !perpDir.empty(); }
  bool usingRound() { return !fcnCenter.empty(); }
  bool usingGyroid() { return !gyroidScale.empty(); }

  template <int DIM>
  axom::primal::Point<double, DIM> roundContourCenter() const
  {
    SLIC_ASSERT(fcnCenter.size() == DIM);
    return axom::primal::Point<double, DIM>(fcnCenter.data());
  }

  template <int DIM>
  axom::primal::Point<double, DIM> gyroidScaleFactor() const
  {
    SLIC_ASSERT(gyroidScale.size() == DIM);
    return axom::primal::Point<double, DIM>(gyroidScale.data());
  }

  template <int DIM>
  axom::primal::Point<double, DIM> inplanePoint() const
  {
    SLIC_ASSERT(inPlane.size() == DIM);
    return axom::primal::Point<double, DIM>(inPlane.data());
  }

  template <int DIM>
  axom::primal::Vector<double, DIM> planeNormal() const
  {
    SLIC_ASSERT(perpDir.size() == DIM);
    return axom::primal::Vector<double, DIM>(perpDir.data());
  }
};

//!@brief Our allocator id, based on execution policy.
static int s_allocatorId = axom::INVALID_ALLOCATOR_ID;  // Set in main.

//!@brief Put a conduit::Node array data into the specified memory space.
template <typename T>
void moveConduitDataToNewMemorySpace(conduit::Node& node, const std::string& path, int allocId)
{
  conduit::Node& dataNode = node.fetch_existing(path);
  SLIC_ASSERT(!dataNode.dtype().is_empty() && !dataNode.dtype().is_object() &&
              !dataNode.dtype().is_list());

  std::size_t count = dataNode.dtype().number_of_elements();
  T* oldPtr = static_cast<T*>(dataNode.data_ptr());
  bool deleteOld = dataNode.is_data_external();
  T* newPtr = axom::allocate<T>(count, allocId);
  axom::copy(newPtr, oldPtr, count * sizeof(T));
  dataNode.set_external(newPtr, count);
  if(deleteOld)
  {
    axom::deallocate(oldPtr);
  }
}

void getIntMinMax(int inVal, int& minVal, int& maxVal, int& sumVal)
{
#ifdef AXOM_USE_MPI
  MPI_Allreduce(&inVal, &minVal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&inVal, &maxVal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&inVal, &sumVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  minVal = inVal;
  maxVal = inVal;
  sumVal = inVal;
#endif
}

Input params;

int myRank = -1, numRanks = -1;  // MPI stuff, set in main().

/**
 \brief Generic computational mesh, to hold cell and node data.
*/
struct BlueprintStructuredMesh
{
public:
  explicit BlueprintStructuredMesh(const std::string& meshFile, const std::string& topologyName)
    : _topologyName(topologyName)
    , _topologyPath("topologies/" + topologyName)
  {
    readBlueprintMesh(meshFile);
    for(int d = 0; d < _mdMesh.number_of_children(); ++d)
    {
      auto dl = domainLengths(d);
      SLIC_INFO_IF(params.isVerbose(), axom::fmt::format("dom[{}] size={}", d, dl));
    }
    _maxSpacing = maxSpacing();
  }

  /// Return the blueprint mesh in a conduit::Node
  conduit::Node& asConduitNode() { return _mdMesh; }

  /// Get number of domains in the multidomain mesh
  axom::IndexType domainCount() const { return _domCount; }

  //// Whether mesh is empty locally.
  bool empty() const { return _domCount == 0; }

  /// Get domain group.
  conduit::Node& domain(axom::IndexType domainIdx)
  {
    SLIC_ASSERT(domainIdx >= 0 && domainIdx < _domCount);
    return _mdMesh.child(domainIdx);
  }

  const conduit::Node& domain(axom::IndexType domainIdx) const
  {
    SLIC_ASSERT(domainIdx >= 0 && domainIdx < _domCount);
    return _mdMesh.child(domainIdx);
  }

  template <int DIM>
  axom::quest::MeshViewUtil<DIM> getDomainView(axom::IndexType domainId)
  {
    return axom::quest::MeshViewUtil<DIM>(domain(domainId), _topologyName);
  }

  template <int DIM>
  axom::quest::MeshViewUtil<DIM> getDomainView(axom::IndexType domainId) const
  {
    return axom::quest::MeshViewUtil<DIM>(domain(domainId), _topologyName);
  }

  /*!
    @brief Get the number of cells in each direction of a blueprint single domain.

    @param domId Index of domain
    @lengths Space for dimension() numbers.
  */
  void domainLengths(axom::IndexType domId, axom::IndexType* lengths) const
  {
    const conduit::Node& dom = domain(domId);
    SLIC_ASSERT_MSG(dom.fetch_existing(_coordsetPath + "/type").as_string() == "explicit",
                    axom::fmt::format("Currently only supporting explicit coordinate types."
                                      "  '{}/type' is '{}'",
                                      _coordsetPath,
                                      dom.fetch_existing(_coordsetPath + "/type").as_string()));
    const conduit::Node& dimsNode = dom.fetch_existing(_topologyPath + "/elements/dims");
    for(int i = 0; i < _ndims; ++i)
    {
      lengths[i] = dimsNode[i].as_int();
    }
  }

  axom::Array<axom::IndexType> domainLengths(axom::IndexType domainId) const
  {
    axom::Array<axom::IndexType> rval(_ndims, _ndims);
    domainLengths(domainId, rval.data());
    return rval;
  }

  /// Returns the number of cells in a domain
  int cellCount(axom::IndexType domId) const
  {
    auto shape = domainLengths(domId);
    int rval = 1;
    for(const auto& l : shape)
    {
      rval *= l;
    }
    return rval;
  }

  /// Returns the number of cells in all mesh domains
  int cellCount() const
  {
    int rval = 0;
    for(int domId = 0; domId < _mdMesh.number_of_children(); ++domId)
    {
      rval += cellCount(domId);
    }
    return rval;
  }

  /// Returns the number of nodes in a domain
  int nodeCount(axom::IndexType domId) const
  {
    auto shape = domainLengths(domId);
    int rval = 1;
    for(const auto& l : shape)
    {
      rval *= 1 + l;
    }
    return rval;
  }

  /// Returns the number of nodes in all mesh domains
  int nodeCount() const
  {
    int rval = 0;
    for(int domId = 0; domId < _mdMesh.number_of_children(); ++domId)
    {
      rval += nodeCount(domId);
    }
    return rval;
  }

  int dimension() const { return _ndims; }

  const std::string& coordsetPath() const { return _coordsetPath; }

  /*!
    @return largest mesh spacing.

    Compute only once, because after that, coordinates data may be
    moved to devices.
  */
  double maxSpacing() const
  {
    if(_maxSpacing >= 0)
    {
      // Max spacing has been computed and cached.
      return _maxSpacing;
    }

    double localRval = 0.0;
    for(axom::IndexType domId = 0; domId < domainCount(); ++domId)
    {
      localRval = std::max(localRval, maxSpacing1(domId));
    }

    double rval = localRval;
#ifdef AXOM_USE_MPI
    MPI_Allreduce(&localRval, &rval, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    return rval;
  }

  /*!
    @return largest mesh spacing in a domain.

    This method takes shortcuts by assuming
    the mesh is structured and cartesian, with explicit coordinates.
  */
  double maxSpacing1(axom::IndexType domId) const
  {
    const conduit::Node& dom = domain(domId);
    const conduit::Node& dimsNode = dom.fetch_existing("topologies/mesh/elements/dims");
    axom::Array<axom::IndexType> ls(_ndims);
    for(int d = 0; d < _ndims; ++d)
    {
      ls[d] = 1 + dimsNode[d].as_int();
    }

    double rval = 0.0;

    const conduit::Node& cVals = dom.fetch_existing(_coordsetPath + "/values");
    if(_ndims == 2)
    {
      axom::ArrayView<const double, 2> xs(cVals["x"].as_double_ptr(), ls[1], ls[0]);
      axom::ArrayView<const double, 2> ys(cVals["y"].as_double_ptr(), ls[1], ls[0]);
      rval = std::max(rval, std::abs(xs(0, 0) - xs(0, 1)));
      rval = std::max(rval, std::abs(ys(0, 0) - ys(1, 0)));
    }
    else
    {
      axom::ArrayView<const double, 3> xs(cVals["x"].as_double_ptr(), ls[2], ls[1], ls[0]);
      axom::ArrayView<const double, 3> ys(cVals["y"].as_double_ptr(), ls[2], ls[1], ls[0]);
      axom::ArrayView<const double, 3> zs(cVals["z"].as_double_ptr(), ls[2], ls[1], ls[0]);
      rval = std::max(rval, std::abs(xs(0, 0, 0) - xs(0, 0, 1)));
      rval = std::max(rval, std::abs(ys(0, 0, 0) - ys(0, 1, 0)));
      rval = std::max(rval, std::abs(zs(0, 0, 0) - zs(1, 0, 0)));
    }
    return rval;
  }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    conduit::Node info;
#ifdef AXOM_USE_MPI
    if(!conduit::blueprint::mpi::verify("mesh", _mdMesh, info, MPI_COMM_WORLD))
#else
    if(!conduit::blueprint::verify("mesh", _mdMesh, info))
#endif
    {
      SLIC_INFO("Invalid blueprint for mesh: \n" << info.to_yaml());
      slic::flushStreams();
      return false;
    }
    return true;
  }

  void printMeshInfo() const { _mdMesh.print(); }

  /*!
    @param[in] path Path to existing data in the blueprint mesh,
    relative to each domain in the mesh.
    @param[in] allocId Allocator id for the new memory space.
    @tparam Type of data being moved.  Should be something Conduit
    supports, i.e., not custom user data.
  */
  template <typename T>
  void moveMeshDataToNewMemorySpace(const std::string& path, int allocId)
  {
    AXOM_ANNOTATE_SCOPE("moveMeshDataToNewMemorySpace");
    for(auto& dom : _mdMesh.children())
    {
      moveConduitDataToNewMemorySpace<T>(dom, path, allocId);
    }
  }

private:
  int _ndims {-1};
  conduit::Node _mdMesh;
  axom::IndexType _domCount;
  bool _coordsAreStrided = false;
  const std::string _topologyName;
  const std::string _topologyPath;
  std::string _coordsetPath;
  double _maxSpacing = -1.0;

  /*!
    @brief Read a blueprint mesh into conduit::Node _mdMesh.
  */
  void readBlueprintMesh(const std::string& meshFilename)
  {
    SLIC_ASSERT(!meshFilename.empty());

    _mdMesh.reset();
#ifdef AXOM_USE_MPI
    conduit::relay::mpi::io::blueprint::load_mesh(meshFilename, _mdMesh, MPI_COMM_WORLD);
#else
    conduit::relay::io::blueprint::load_mesh(meshFilename, _mdMesh);
#endif
    SLIC_ASSERT(conduit::blueprint::mesh::is_multi_domain(_mdMesh));
    _domCount = conduit::blueprint::mesh::number_of_domains(_mdMesh);

    if(_domCount > 0)
    {
      SLIC_ASSERT(_mdMesh[0].has_path(_topologyPath));
      auto coordsetName = _mdMesh[0].fetch_existing(_topologyPath + "/coordset").as_string();
      _coordsetPath = axom::fmt::format("coordsets/{}", coordsetName);
      SLIC_ASSERT(_mdMesh[0].has_path(_coordsetPath));

      _coordsAreStrided =
        _mdMesh[0].fetch_existing(_topologyPath + "/elements/dims").has_child("strides");
      const conduit::Node coordsetNode = _mdMesh[0].fetch_existing(_coordsetPath);
      _ndims = conduit::blueprint::mesh::coordset::dims(coordsetNode);
    }
#ifdef AXOM_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &_ndims, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    SLIC_ASSERT(_ndims > 0);

    SLIC_ASSERT(isValid());
  }
};  // BlueprintStructuredMesh

/// Output some timing stats
void printTimingStats(axom::utilities::Timer& t, const std::string& description)
{
  auto getDoubleMinMax = [](double inVal, double& minVal, double& maxVal, double& sumVal) {
#ifdef AXOM_USE_MPI
    MPI_Allreduce(&inVal, &minVal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &maxVal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &sumVal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    minVal = inVal;
    maxVal = inVal;
    sumVal = inVal;
#endif
  };

  {
    double minCompute, maxCompute, sumCompute;
    getDoubleMinMax(t.elapsedTimeInSec(), minCompute, maxCompute, sumCompute);

    const auto count = t.cycleCount();
    if(count > 1)
    {
      SLIC_INFO(
        axom::fmt::format("'{}' took {{avg:{}, min:{}, max:{}}} seconds (avg of {} samples)",
                          description,
                          sumCompute / count / numRanks,
                          minCompute / count,
                          maxCompute / count,
                          count));
    }
    else
    {
      SLIC_INFO(axom::fmt::format("'{}' took {{avg:{}, min:{}, max:{}}} seconds",
                                  description,
                                  sumCompute / numRanks,
                                  minCompute,
                                  maxCompute,
                                  count));
    }
  }
}

/// Write blueprint mesh to disk
void saveMesh(const conduit::Node& mesh, const std::string& filename)
{
  AXOM_ANNOTATE_SCOPE("save mesh (conduit)");

#ifdef AXOM_USE_MPI
  conduit::relay::mpi::io::blueprint::save_mesh(mesh, filename, "hdf5", MPI_COMM_WORLD);
#else
  conduit::relay::io::blueprint::save_mesh(mesh, filename, "hdf5");
#endif
}

#if defined(AXOM_USE_SIDRE)
/// Write blueprint mesh to disk
void saveMesh(const sidre::Group& mesh, const std::string& filename)
{
  AXOM_ANNOTATE_SCOPE("save mesh (sidre)");

  conduit::Node tmpMesh;
  mesh.createNativeLayout(tmpMesh);
  {
    conduit::Node info;
  #ifdef AXOM_USE_MPI
    if(!conduit::blueprint::mpi::verify("mesh", tmpMesh, info, MPI_COMM_WORLD))
  #else
    if(!conduit::blueprint::verify("mesh", tmpMesh, info))
  #endif
    {
      SLIC_INFO("Invalid blueprint for mesh: \n" << info.to_yaml());
      slic::flushStreams();
      assert(false);
    }
    // info.print();
  }
  saveMesh(tmpMesh, filename);
}
#endif

template <typename T, int DIM>
T product(const axom::StackArray<T, DIM>& a)
{
  T rval = a[0];
  for(int d = 1; d < DIM; ++d)
  {
    rval *= a[d];
  }
  return rval;
}

template <typename T, int DIM, typename U>
static void addToStackArray(axom::StackArray<T, DIM>& a, U b)
{
  for(int d = 0; d < DIM; ++d)
  {
    a[d] += b;
  }
}

/*!
  @brief Strategy pattern for supporting a variety of contour types.

  The strategy encapsulates the scalar functions and things related to
  it.
*/
template <int DIM>
struct ContourTestStrategy
{
  using PointType = axom::primal::Point<double, DIM>;

  //!@brief Return test name.
  virtual std::string testName() const = 0;

  //!@brief Return field name for storing nodal function.
  virtual std::string functionName() const = 0;

  //!@brief Return error tolerance for contour mesh accuracy check.
  virtual double errorTolerance() const = 0;

  //!@brief Return the analytical value of the scalar field.
  virtual double valueAt(const PointType& pt) const = 0;

  virtual ~ContourTestStrategy() { }
};

template <int DIM, typename ExecSpace>
struct ContourTestBase
{
  static constexpr auto MemorySpace = axom::execution_space<ExecSpace>::memory_space;
  using PointType = axom::primal::Point<double, DIM>;
  // ContourTestBase(const std::shared_ptr<ContourTestStrategy<DIM>>& testStrategy)
  ContourTestBase()
    : m_testStrategies()
    , m_parentCellIdField("parentCellIds")
    , m_domainIdField("domainIdField")
  { }
  virtual ~ContourTestBase() { }

  void addTestStrategy(const std::shared_ptr<ContourTestStrategy<DIM>>& testStrategy)
  {
    m_testStrategies.push_back(testStrategy);
    SLIC_INFO(axom::fmt::format("Add test {}.", testStrategy->testName()));
  }

  axom::Array<std::shared_ptr<ContourTestStrategy<DIM>>> m_testStrategies;
  //!@brief Prefix sum of facet counts from test strategies.
  axom::Array<axom::IndexType> m_strategyFacetPrefixSum;

  const std::string m_parentCellIdField;
  const std::string m_domainIdField;

  int runTest(BlueprintStructuredMesh& computationalMesh)
  {
    AXOM_ANNOTATE_SCOPE("runTest");

    // Conduit data is in host memory, move to devices for testing.
    if(s_allocatorId != axom::execution_space<axom::SEQ_EXEC>::allocatorID())
    {
      AXOM_ANNOTATE_SCOPE("move mesh to device memory");

      const std::string axes[3] = {"x", "y", "z"};
      for(int d = 0; d < DIM; ++d)
      {
        computationalMesh.moveMeshDataToNewMemorySpace<double>(
          computationalMesh.coordsetPath() + "/values/" + axes[d],
          s_allocatorId);
      }
      for(const auto& strategy : m_testStrategies)
      {
        computationalMesh.moveMeshDataToNewMemorySpace<double>(
          axom::fmt::format("fields/{}/values", strategy->functionName()),
          s_allocatorId);
      }
      computationalMesh.moveMeshDataToNewMemorySpace<int>(
        axom::fmt::format("fields/{}/values", "mask"),
        s_allocatorId);
    }

#if defined(AXOM_USE_UMPIRE)
    /*
      Make sure data is correctly on host or device.
      We don't test with Unified memory because it's too forgiving.
    */
    if(!computationalMesh.empty())
    {
      AXOM_ANNOTATE_SCOPE("move mesh from unified memory");

      std::string resourceName = "HOST";
      umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
      for(const auto& strategy : m_testStrategies)
      {
        const std::string dataPath = axom::fmt::format("fields/{}/values", strategy->functionName());
        void* dataPtr = computationalMesh.domain(0).fetch_existing(dataPath).data_ptr();
        bool dataFromUmpire = rm.hasAllocator(dataPtr);
        if(dataFromUmpire)
        {
          umpire::Allocator allocator = rm.getAllocator(dataPtr);
          resourceName = allocator.getName();
        }
        SLIC_INFO(axom::fmt::format("Testing with policy {} and function data on {}",
                                    params.policy,
                                    resourceName));
        if(params.policy == axom::runtime_policy::Policy::seq)
        {
          SLIC_ASSERT(resourceName == "HOST");
        }
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
        else if(params.policy == axom::runtime_policy::Policy::omp)
        {
          SLIC_ASSERT(resourceName == "HOST");
        }
  #endif
        else
        {
          SLIC_ASSERT(resourceName == "DEVICE");
        }
      }
    }
#endif

    // One-time initializations
    axom::utilities::Timer initializationTimer(false);

    // Entire objectRepCount loop.
    axom::utilities::Timer objectRepLoopTimer(false);

    // All contourGenCount loops.
    axom::utilities::Timer contourGenLoopTimer(false);

    // params.objectRepCount setMesh calls
    axom::utilities::Timer setMeshTimer(false);

    // Time steady-state computeIsocontour calls
    axom::utilities::Timer contourTimer(false);

    // Time first contourIsocontour call after setMesh
    axom::utilities::Timer contourTimerM(false);

    std::unique_ptr<quest::MarchingCubes> mcPtr;
    const auto objectLoopName = axom::fmt::format("objectRepLoop {}", params.objectRepCount);
    AXOM_ANNOTATE_BEGIN(objectLoopName);
    objectRepLoopTimer.start();
    for(int j = 0; j < params.objectRepCount; ++j)
    {
      if(!mcPtr)
      {
        AXOM_ANNOTATE_SCOPE("MCInit");
        initializationTimer.start();
        mcPtr =
          std::make_unique<quest::MarchingCubes>(params.policy, s_allocatorId, params.dataParallelism);
        mcPtr->setMesh(computationalMesh.asConduitNode(), "mesh", "mask");
        initializationTimer.stop();
      }
      auto& mc = *mcPtr;

      // Clear and set MarchingCubes object for a "new" mesh.
      setMeshTimer.start();
      mc.setMesh(computationalMesh.asConduitNode(), "mesh", "mask");
      setMeshTimer.stop();

#ifdef AXOM_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      contourGenLoopTimer.start();
      for(int i = 0; i < params.contourGenCount; ++i)
      {
        SLIC_DEBUG(axom::fmt::format("MarchingCubes object rep {} of {}, contour run {} of {}:",
                                     j,
                                     params.objectRepCount,
                                     i,
                                     params.contourGenCount));
        mc.clearOutput();
        m_strategyFacetPrefixSum.clear();
        m_strategyFacetPrefixSum.push_back(0);
        for(const auto& strategy : m_testStrategies)
        {
          mc.setFunctionField(strategy->functionName());
          for(int iMask = 0; iMask < params.maskCount; ++iMask)
          {
            mc.setMaskValue(iMask);
            if(i == 0)
            {
              contourTimerM.start();
            }
            else
            {
              contourTimer.start();
            }
            mc.computeIsocontour(params.contourVal);
            if(i == 0)
            {
              contourTimerM.stop();
            }
            else
            {
              contourTimer.stop();
            }
          }
          m_strategyFacetPrefixSum.push_back(mc.getContourFacetCount());
        }
      }
      contourGenLoopTimer.stop();
    }
    objectRepLoopTimer.stop();
    AXOM_ANNOTATE_END(objectLoopName);
    SLIC_INFO(axom::fmt::format("Finished {} object reps x {} contour reps",
                                params.objectRepCount,
                                params.contourGenCount));
    printTimingStats(initializationTimer, axom::fmt::format("init"));
    printTimingStats(contourTimerM, axom::fmt::format("setMeshContour"));
    printTimingStats(setMeshTimer, axom::fmt::format("setMesh"));
    printTimingStats(contourTimer, axom::fmt::format("steady-contour"));
    printTimingStats(contourTimerM, axom::fmt::format("first-contour"));
    printTimingStats(contourGenLoopTimer, axom::fmt::format("contourGenLoop"));
    printTimingStats(objectRepLoopTimer, axom::fmt::format("objectRepLoop"));

    auto& mc = *mcPtr;
    printRunStats(mc);

    // Return conduit data to host memory.
    if(s_allocatorId != axom::execution_space<axom::SEQ_EXEC>::allocatorID())
    {
      AXOM_ANNOTATE_SCOPE("copy mesh back to host memory");

      const std::string axes[3] = {"x", "y", "z"};
      for(int d = 0; d < DIM; ++d)
      {
        computationalMesh.moveMeshDataToNewMemorySpace<double>(
          computationalMesh.coordsetPath() + "/values/" + axes[d],
          axom::execution_space<axom::SEQ_EXEC>::allocatorID());
      }
      for(const auto& strategy : m_testStrategies)
      {
        computationalMesh.moveMeshDataToNewMemorySpace<double>(
          axom::fmt::format("fields/{}/values", strategy->functionName()),
          axom::execution_space<axom::SEQ_EXEC>::allocatorID());
      }
      computationalMesh.moveMeshDataToNewMemorySpace<int>(
        axom::fmt::format("fields/{}/values", "mask"),
        axom::execution_space<axom::SEQ_EXEC>::allocatorID());
    }

    // Put contour mesh in a mint object for error checking and output.
    AXOM_ANNOTATE_BEGIN("error checking");

    AXOM_ANNOTATE_BEGIN("convert to mint mesh");

#ifdef AXOM_MINT_USE_SIDRE
    std::string sidreGroupName = "contour_mesh";
    sidre::DataStore objectDS;
    auto* meshGroup = objectDS.getRoot()->createGroup(sidreGroupName);
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> contourMesh(
      DIM,
      DIM == 2 ? mint::CellType::SEGMENT : mint::CellType::TRIANGLE,
      meshGroup);
#else
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> contourMesh(
      DIM,
      DIM == 2 ? mint::CellType::SEGMENT : mint::CellType::TRIANGLE);
#endif
    axom::utilities::Timer extractTimer(false);
    extractTimer.start();
    mc.populateContourMesh(contourMesh, m_parentCellIdField, m_domainIdField);
    extractTimer.stop();
    printTimingStats(extractTimer, "extract");

    {
      axom::Array<axom::IndexType, 2> facetNodeIds;
      axom::Array<double, 2> facetNodeCoords;
      axom::Array<axom::IndexType, 1> facetParentIds;
      axom::Array<axom::IndexType> facetDomainIds;
      mc.relinquishContourData(facetNodeIds, facetNodeCoords, facetParentIds, facetDomainIds);
      SLIC_ASSERT(mc.getContourFacetCount() == 0);
    }
    AXOM_ANNOTATE_END("convert to mint mesh");

    int localErrCount = 0;
    if(params.checkResults)
    {
      localErrCount += checkContourSurface(contourMesh, params.contourVal, "diff");

      localErrCount += checkContourCellLimits(computationalMesh, contourMesh);

      localErrCount += checkCellsContainingContour(computationalMesh, contourMesh);
    }

#if defined(AXOM_MINT_USE_SIDRE)
    if(contourMesh.hasSidreGroup())
    {
      assert(contourMesh.getSidreGroup() == meshGroup);
      // Write contour mesh to file.
      std::string outputName = "contour_mesh";
      saveMesh(*contourMesh.getSidreGroup(), outputName);
      SLIC_INFO(axom::fmt::format("Wrote contour mesh to {}", outputName));
    }
    objectDS.getRoot()->destroyGroupAndData(sidreGroupName);
#endif

    AXOM_ANNOTATE_END("error checking");

    return localErrCount;
  }

  void printRunStats(const quest::MarchingCubes& mc)
  {
    {
      int mn, mx, sum;
      getIntMinMax(mc.getContourCellCount(), mn, mx, sum);
      SLIC_INFO(axom::fmt::format("Contour mesh has {{min:{}, max:{}, sum:{}, avg:{}}} cells",
                                  mn,
                                  mx,
                                  sum,
                                  (double)sum / numRanks));
    }
    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("Contour mesh has locally {} cells, {} nodes.",
                                   mc.getContourCellCount(),
                                   mc.getContourNodeCount()));
  }

  void computeNodalDistance(BlueprintStructuredMesh& bpMesh, ContourTestStrategy<DIM>& strat)
  {
    AXOM_ANNOTATE_SCOPE("computeNodalDistance");

    SLIC_ASSERT(bpMesh.dimension() == DIM);
    for(int domId = 0; domId < bpMesh.domainCount(); ++domId)
    {
      auto domainView = bpMesh.getDomainView<DIM>(domId);

      // Create nodal function data with ghosts like node coords.
      domainView.createField(strat.functionName(),
                             "vertex",
                             conduit::DataType::float64(domainView.getCoordsCountWithGhosts()),
                             domainView.getCoordsStrides(),
                             domainView.getCoordsOffsets());
    }

    for(int domId = 0; domId < bpMesh.domainCount(); ++domId)
    {
      auto domainView = bpMesh.getDomainView<DIM>(domId);
      const auto coordsViews = domainView.getConstCoordsViews(false);
      axom::ArrayView<double, DIM> fieldView =
        domainView.template getFieldView<double>(strat.functionName(), false);
      for(int d = 0; d < DIM; ++d)
      {
        SLIC_ASSERT(coordsViews[d].shape() == fieldView.shape());
      }
      populateNodalDistance(coordsViews, fieldView, strat);
    }
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type populateNodalDistance(
    const axom::StackArray<axom::ArrayView<const double, DIM>, DIM>& coordsViews,
    axom::ArrayView<double, DIM>& fieldView,
    ContourTestStrategy<DIM>& strat)
  {
    AXOM_ANNOTATE_SCOPE("populateNodalDistance 2D");

    const auto& fieldShape = fieldView.shape();
    for(int d = 0; d < DIM; ++d)
    {
      SLIC_ASSERT(coordsViews[d].shape() == fieldShape);
    }

    for(axom::IndexType j = 0; j < fieldShape[1]; ++j)
    {
      for(axom::IndexType i = 0; i < fieldShape[0]; ++i)
      {
        PointType pt;
        for(int d = 0; d < DIM; ++d)
        {
          pt[d] = coordsViews[d](i, j);
        }
        fieldView(i, j) = strat.valueAt(pt);
      }
    }
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type populateNodalDistance(
    const axom::StackArray<axom::ArrayView<const double, DIM>, DIM>& coordsViews,
    axom::ArrayView<double, DIM>& fieldView,
    ContourTestStrategy<DIM>& strat)
  {
    AXOM_ANNOTATE_SCOPE("populateNodalDistance 3D");

    const auto& fieldShape = fieldView.shape();
    for(int d = 0; d < DIM; ++d)
    {
      SLIC_ASSERT(coordsViews[d].shape() == fieldShape);
    }

    for(axom::IndexType k = 0; k < fieldShape[2]; ++k)
    {
      for(axom::IndexType j = 0; j < fieldShape[1]; ++j)
      {
        for(axom::IndexType i = 0; i < fieldShape[0]; ++i)
        {
          PointType pt;
          for(int d = 0; d < DIM; ++d)
          {
            pt[d] = coordsViews[d](i, j, k);
          }
          fieldView(i, j, k) = strat.valueAt(pt);
        }
      }
    }
  }

  void addMaskField(BlueprintStructuredMesh& bpMesh)
  {
    std::string maskFieldName = "mask";
    axom::StackArray<axom::IndexType, DIM> zeros;
    for(int d = 0; d < DIM; ++d)
    {
      zeros[d] = 0;
    }
    for(axom::IndexType domId = 0; domId < bpMesh.domainCount(); ++domId)
    {
      auto domainView = bpMesh.getDomainView<DIM>(domId);
      auto cellCount = domainView.getCellCount();
      auto slowestDirs = domainView.getConstCoordsViews()[0].mapping().slowestDirs();
      axom::StackArray<axom::IndexType, DIM> fastestDirs;
      for(int d = 0; d < DIM; ++d)
      {
        fastestDirs[d] = slowestDirs[DIM - 1 - d];
      }
      domainView.createField(maskFieldName,
                             "element",
                             conduit::DataType::c_int(cellCount),
                             zeros,
                             zeros,
                             fastestDirs);
      auto maskView = domainView.template getFieldView<int>(maskFieldName);
      int maskCount = params.maskCount;
      axom::for_all<axom::SEQ_EXEC>(
        0,
        cellCount,
        AXOM_LAMBDA(axom::IndexType cellId) { maskView.flatIndex(cellId) = (cellId % maskCount); });
    }
  }

  void computeNodalDistance(BlueprintStructuredMesh& bpMesh)
  {
    for(auto& strategy : m_testStrategies)
    {
      computeNodalDistance(bpMesh, *strategy);
    }
  }

  /**
     Check for errors in the surface contour mesh.
     - analytical scalar value at surface points should be
       contourVal, within tolerance zero.
  */
  int checkContourSurface(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh,
                          double contourVal,
                          const std::string& diffField = {})
  {
    AXOM_ANNOTATE_SCOPE("checkContourSurface");
    double* diffPtr = nullptr;
    if(!diffField.empty())
    {
      diffPtr = contourMesh.createField<double>(diffField, axom::mint::NODE_CENTERED);
    }

    int errCount = 0;
    for(axom::IndexType iStrat = 0; iStrat < m_testStrategies.size(); ++iStrat)
    {
      auto contourNodeBegin = DIM * m_strategyFacetPrefixSum[iStrat];
      auto contourNodeEnd = DIM * m_strategyFacetPrefixSum[iStrat + 1];

      auto& strategy = *m_testStrategies[iStrat];
      double tol = strategy.errorTolerance();

      PointType pt;
      for(axom::IndexType iNode = contourNodeBegin; iNode < contourNodeEnd; ++iNode)
      {
        contourMesh.getNode(iNode, pt.data());
        double analyticalVal = strategy.valueAt(pt);
        double diff = std::abs(analyticalVal - contourVal);
        if(diffPtr)
        {
          diffPtr[iNode] = diff;
        }
        if(diff > tol)
        {
          ++errCount;
          SLIC_INFO_IF(
            params.isVerbose(),
            axom::fmt::format("checkContourSurface: node {} at {} has dist {}, off by {}",
                              iNode,
                              pt,
                              analyticalVal,
                              diff));
        }
      }
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("checkContourSurface: found {} errors outside tolerance of {}",
                                     errCount,
                                     tol));
    }
    return errCount;
  }

  //!@brief Get view of output domain id data.
  axom::ArrayView<const axom::quest::MarchingCubes::DomainIdType> getDomainIdView(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh) const
  {
    const auto* ptr =
      contourMesh.getFieldPtr<axom::quest::MarchingCubes::DomainIdType>(m_domainIdField,
                                                                        axom::mint::CELL_CENTERED);
    axom::ArrayView<const axom::quest::MarchingCubes::DomainIdType> view(
      ptr,
      contourMesh.getNumberOfCells());
    return view;
  }

  //!@brief Get view of output parent cell idx data.
  axom::ArrayView<const axom::StackArray<axom::IndexType, DIM>> get_parent_cell_idx_view(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh) const
  {
    axom::IndexType numIdxComponents = -1;
    axom::IndexType* ptr = contourMesh.getFieldPtr<axom::IndexType>(m_parentCellIdField,
                                                                    axom::mint::CELL_CENTERED,
                                                                    numIdxComponents);

    SLIC_ASSERT(numIdxComponents == DIM);

    axom::ArrayView<const axom::StackArray<axom::IndexType, DIM>> view(
      (axom::StackArray<axom::IndexType, DIM>*)ptr,
      contourMesh.getNumberOfCells());
    return view;
  }

  axom::ArrayView<const axom::IndexType> get_parent_cell_id_view(
    const axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh) const
  {
    axom::IndexType numIdxComponents = -1;
    const axom::IndexType* ptr = contourMesh.getFieldPtr<axom::IndexType>(m_parentCellIdField,
                                                                          axom::mint::CELL_CENTERED,
                                                                          numIdxComponents);

    SLIC_ASSERT(numIdxComponents == 1);

    axom::ArrayView<const axom::IndexType> view((const axom::IndexType*)ptr,
                                                contourMesh.getNumberOfCells());
    return view;
  }

  /**
     Check that generated cells fall within their parents.
  */
  int checkContourCellLimits(BlueprintStructuredMesh& computationalMesh,
                             axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh)
  {
    AXOM_ANNOTATE_SCOPE("checkContourCellLimits");

    int errCount = 0;

    auto parentCellIdView = get_parent_cell_id_view(contourMesh);
    auto domainIdView = getDomainIdView(contourMesh);

    const axom::IndexType domainCount = computationalMesh.domainCount();
    axom::Array<typename axom::quest::MeshViewUtil<DIM>::ConstCoordsViewsType> allCoordsViews(
      domainCount);
    for(int iDomain = 0; iDomain < domainCount; ++iDomain)
    {
      auto domainView = computationalMesh.getDomainView<DIM>(iDomain);
      allCoordsViews[iDomain] = domainView.getConstCoordsViews(false);
    }

    std::map<axom::quest::MarchingCubes::DomainIdType, axom::quest::MarchingCubes::DomainIdType>
      domainIdToContiguousId;
    for(int iDomain = 0; iDomain < domainCount; ++iDomain)
    {
      const auto& dom = computationalMesh.domain(iDomain);
      axom::quest::MarchingCubes::DomainIdType domainId = iDomain;
      if(dom.has_path("state/domain_id"))
      {
        domainId = dom.fetch_existing("state/domain_id").to_value();
      }
      domainIdToContiguousId[domainId] = iDomain;
    }

    // Indexers to translate between flat and multidim indices.
    axom::Array<axom::MDMapping<DIM>> mappings(domainCount);
    for(int d = 0; d < domainCount; ++d)
    {
      axom::StackArray<axom::IndexType, DIM> domShape;
      computationalMesh.domainLengths(d, domShape);
      mappings[d].initializeShape(domShape,
                                  axom::MDMapping<DIM>(allCoordsViews[d][0].strides()).slowestDirs());
    }

    auto elementGreaterThan = [](const axom::primal::Vector<double, DIM>& a, double b) {
      bool result(true);
      for(int d = 0; d < DIM; ++d)
      {
        result &= a[d] < b;
      }
      return result;
    };

    for(axom::IndexType iStrat = 0; iStrat < m_testStrategies.size(); ++iStrat)
    {
      auto contourCellBegin = m_strategyFacetPrefixSum[iStrat];
      auto contourCellEnd = m_strategyFacetPrefixSum[iStrat + 1];
      for(axom::IndexType iContourCell = contourCellBegin; iContourCell < contourCellEnd;
          ++iContourCell)
      {
        axom::quest::MarchingCubes::DomainIdType domainId = domainIdView[iContourCell];
        axom::quest::MarchingCubes::DomainIdType contiguousIndex = domainIdToContiguousId[domainId];
        typename axom::quest::MeshViewUtil<DIM>::ConstCoordsViewsType& coordsViews =
          allCoordsViews[contiguousIndex];

        axom::IndexType parentCellId = parentCellIdView[iContourCell];

        axom::StackArray<axom::IndexType, DIM> parentCellIdx =
          mappings[contiguousIndex].toMultiIndex(parentCellId);
        axom::StackArray<axom::IndexType, DIM> upperIdx = parentCellIdx;
        addToStackArray(upperIdx, 1);

        PointType lower, upper;
        for(int d = 0; d < DIM; ++d)
        {
          lower[d] = coordsViews[d][parentCellIdx];
          upper[d] = coordsViews[d][upperIdx];
        }
        axom::primal::BoundingBox<double, DIM> parentCellBox(lower, upper);
        auto tol = axom::numerics::floating_point_limits<double>::epsilon();
        axom::primal::BoundingBox<double, DIM> big(parentCellBox);
        big.expand(tol);
        axom::primal::BoundingBox<double, DIM> small(parentCellBox);
        auto range = parentCellBox.range();
        bool checkSmall = elementGreaterThan(range, tol);
        if(checkSmall)
        {
          small.expand(-tol);
        }

        axom::IndexType* cellNodeIds = contourMesh.getCellNodeIDs(iContourCell);
        const axom::IndexType cellNodeCount = contourMesh.getNumberOfCellNodes(iContourCell);

        for(axom::IndexType nn = 0; nn < cellNodeCount; ++nn)
        {
          PointType nodeCoords;
          contourMesh.getNode(cellNodeIds[nn], nodeCoords.data());

          if(!big.contains(nodeCoords))
          {
            ++errCount;
            SLIC_INFO_IF(params.isVerbose(),
                         axom::fmt::format("checkContourCellLimits: node {} at {} "
                                           "too far outside parent cell boundary.",
                                           cellNodeIds[nn],
                                           nodeCoords));
          }

          if(checkSmall && small.contains(nodeCoords))
          {
            ++errCount;
            SLIC_INFO_IF(params.isVerbose(),
                         axom::fmt::format("checkContourCellLimits: node {} at {} "
                                           "too far inside parent cell boundary.",
                                           cellNodeIds[nn],
                                           nodeCoords));
          }
        }
      }
    }

    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("checkContourCellLimits: found {} "
                                   "nodes not on parent cell boundary.",
                                   errCount));
    return errCount;
  }

  /*!
    Check that computational cells that contain the contour value
    have at least one contour mesh cell.
  */
  int checkCellsContainingContour(BlueprintStructuredMesh& computationalMesh,
                                  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh)
  {
    AXOM_ANNOTATE_SCOPE("checkCellsContainingContour");

    int errCount = 0;

    auto parentCellIdView = get_parent_cell_id_view(contourMesh);
    auto domainIdView = getDomainIdView(contourMesh);

    const axom::IndexType domainCount = computationalMesh.domainCount();

    //
    // Compute mapping to look up domains from the domain id.
    //
    std::map<axom::IndexType, axom::IndexType> domainIdToContiguousId;
    for(axom::quest::MarchingCubes::DomainIdType iDomain = 0; iDomain < domainCount; ++iDomain)
    {
      const auto& dom = computationalMesh.domain(iDomain);
      axom::quest::MarchingCubes::DomainIdType domainId = iDomain;
      if(dom.has_path("state/domain_id"))
      {
        domainId = dom.fetch_existing("state/domain_id").to_value();
      }
      domainIdToContiguousId[domainId] = iDomain;
    }

    /*
      If strategy iStrat creates a contour through cell cellId of
      domain contiguousDomainId, then set bit iStrat in hasContours:
      hasContours[contiguousDomainId][cellId] & (1 < iStrat).
    */
    axom::Array<axom::ArrayView<const double, DIM, MemorySpace>> fcnViews(domainCount);
    axom::Array<axom::MDMapping<DIM>> cellIndexers(domainCount);
    axom::Array<axom::Array<axom::IndexType>> hasContours(domainCount);
    for(axom::IndexType domId = 0; domId < domainCount; ++domId)
    {
      axom::quest::MeshViewUtil<DIM> domainView = computationalMesh.getDomainView<DIM>(domId);

      const axom::IndexType cellCount = domainView.getCellCount();

      axom::Array<axom::IndexType>& hasContour = hasContours[domId];
      hasContour.resize(cellCount, 0);
    }

    for(int iStrat = 0; iStrat < m_testStrategies.size(); ++iStrat)
    {
      auto contourCellBegin = m_strategyFacetPrefixSum[iStrat];
      auto contourCellEnd = m_strategyFacetPrefixSum[iStrat + 1];
      axom::IndexType bitFlag = (1 << iStrat);
      for(axom::IndexType iContourCell = contourCellBegin; iContourCell < contourCellEnd;
          ++iContourCell)
      {
        axom::quest::MarchingCubes::DomainIdType domainId = domainIdView[iContourCell];
        axom::quest::MarchingCubes::DomainIdType contiguousId = domainIdToContiguousId[domainId];
        const axom::IndexType parentCellId = parentCellIdView[iContourCell];
        hasContours[contiguousId][parentCellId] |= bitFlag;
      }
    }

    // Verify that cells marked by hasContours touches the contour
    // and other cells don't.
    for(axom::IndexType domId = 0; domId < domainCount; ++domId)
    {
      auto domainView = computationalMesh.getDomainView<DIM>(domId);

      axom::StackArray<axom::IndexType, DIM> domLengths;
      computationalMesh.domainLengths(domId, domLengths);
      assert(domLengths == domainView.getRealShape("element"));

      const axom::IndexType parentCellCount = domainView.getCellCount();
      // axom::Array<bool> hasContour(parentCellCount, parentCellCount);

      for(axom::IndexType parentCellId = 0; parentCellId < parentCellCount; ++parentCellId)
      {
        const axom::IndexType hasContourBits = hasContours[domId][parentCellId];

        for(axom::IndexType iStrat = 0; iStrat < m_testStrategies.size(); ++iStrat)
        {
          auto& strategy = *m_testStrategies[iStrat];
          const axom::IndexType iStratBit = (1 << iStrat);

          const auto& fcnView =
            domainView.template getConstFieldView<double>(strategy.functionName(), false);

          axom::MDMapping<DIM> cellMDMapper(domLengths, axom::MDMapping<DIM>(fcnView.strides()));

          axom::StackArray<axom::IndexType, DIM> parentCellIdx =
            cellMDMapper.toMultiIndex(parentCellId);

          // Compute min and max function values in the cell.
          double minFcnValue = axom::numeric_limits<double>::max();
          double maxFcnValue = axom::numeric_limits<double>::min();
          constexpr short int cornerCount = (1 << DIM);  // Number of nodes in a cell.
          for(short int cornerId = 0; cornerId < cornerCount; ++cornerId)
          {
            // Compute multidim index of current corner of parent cell.
            axom::StackArray<axom::IndexType, DIM> cornerIdx = parentCellIdx;
            for(int d = 0; d < DIM; ++d)
            {
              if(cornerId & (1 << d))
              {
                ++cornerIdx[d];
              }
            }

            double fcnValue = fcnView[cornerIdx];
            minFcnValue = std::min(minFcnValue, fcnValue);
            maxFcnValue = std::max(maxFcnValue, fcnValue);
          }

          const bool hasContour = hasContourBits & iStratBit;

          bool touchesContour =
            (minFcnValue <= params.contourVal && maxFcnValue >= params.contourVal);
          // If the min or max values in the cell is close to params.contourVal
          // touchesContour and hasCont can go either way.  So give it a pass.
          if(minFcnValue == params.contourVal || maxFcnValue == params.contourVal)
          {
            touchesContour = hasContour;
          }

          if(touchesContour != hasContour)
          {
            ++errCount;
            SLIC_INFO_IF(
              params.isVerbose(),
              axom::fmt::format("checkCellsContainingContour: cell {}: hasContour "
                                "({}) and touchesContour ({}) don't agree for strategy {}.",
                                parentCellIdx,
                                hasContour,
                                touchesContour,
                                strategy.testName()));
          }
        }
      }
    }
    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("checkCellsContainingContour: found {} "
                                   "misrepresented computational cells.",
                                   errCount));
    return errCount;
  }
};

template <int DIM>
struct PlanarTestStrategy : public ContourTestStrategy<DIM>
{
  using PointType = axom::primal::Point<double, DIM>;
  PlanarTestStrategy(const axom::primal::Vector<double, DIM>& perpDir, const PointType& inPlane)
    : ContourTestStrategy<DIM>()
    , _plane(perpDir.unitVector(), inPlane)
    , _errTol(axom::numerics::floating_point_limits<double>::epsilon())
  { }
  virtual std::string testName() const override { return std::string("planar"); }
  virtual std::string functionName() const override { return std::string("dist_to_plane"); }
  double errorTolerance() const override { return _errTol; }
  virtual double valueAt(const PointType& pt) const override { return _plane.signedDistance(pt); }
  const axom::primal::Plane<double, DIM> _plane;
  double _errTol;
};

template <int DIM>
struct RoundTestStrategy : public ContourTestStrategy<DIM>
{
  using PointType = axom::primal::Point<double, DIM>;
  RoundTestStrategy(const PointType& center)
    : ContourTestStrategy<DIM>()
    , _sphere(center, 0.0)
    , _errTol(1e-3)
  { }
  virtual std::string testName() const override { return std::string("round"); }
  virtual std::string functionName() const override { return std::string("dist_to_center"); }
  double errorTolerance() const override { return _errTol; }
  virtual double valueAt(const PointType& pt) const override
  {
    return _sphere.computeSignedDistance(pt);
  }
  void setToleranceByLongestEdge(const BlueprintStructuredMesh& bsm)
  {
    // Heuristic of appropriate error tolerance.
    double maxSpacing = bsm.maxSpacing();
    _errTol = 0.1 * maxSpacing;
  }
  const axom::primal::Sphere<double, DIM> _sphere;
  double _errTol;
};

template <int DIM>
struct GyroidTestStrategy : public ContourTestStrategy<DIM>
{
  using PointType = axom::primal::Point<double, DIM>;
  GyroidTestStrategy(const PointType& scale, double offset)
    : ContourTestStrategy<DIM>()
    , _scale(scale)
    , _offset(offset)
  { }
  virtual std::string testName() const override { return std::string("gyroid"); }
  virtual std::string functionName() const override { return std::string("gyroid_fcn"); }
  double errorTolerance() const override { return _errTol; }
  virtual double valueAt(const PointType& pt) const override
  {
    if(DIM == 3)
    {
      return sin(pt[0] * _scale[0]) * cos(pt[1] * _scale[1]) +
        sin(pt[1] * _scale[1]) * cos(pt[2] * _scale[2]) +
        sin(pt[2] * _scale[2]) * cos(pt[0] * _scale[0]) + _offset;
    }
    else
    {
      // Use the 3D function, with z=0.
      return sin(pt[0] * _scale[0]) * cos(pt[1] * _scale[1]) + sin(pt[1] * _scale[1]) + _offset;
    }
  }
  void setToleranceByLongestEdge(const BlueprintStructuredMesh& bsm)
  {
    // Heuristic of appropriate error tolerance.
    double maxSpacing = bsm.maxSpacing();
    axom::primal::Vector<double, DIM> v(_scale);
    _errTol = 0.1 * v.norm() * maxSpacing;
  }
  const PointType _scale;
  const double _offset;
  double _errTol;
};

///
int allocatorIdToTest(axom::runtime_policy::Policy policy)
{
#if defined(AXOM_USE_UMPIRE)
  //---------------------------------------------------------------------------
  // Memory resource.  For testing, choose device memory if appropriate.
  //---------------------------------------------------------------------------
  int allocatorID =
    policy == RuntimePolicy::seq ? axom::detail::getAllocatorID<axom::MemorySpace::Host>() :
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    policy == RuntimePolicy::omp ? axom::detail::getAllocatorID<axom::MemorySpace::Host>()
    :
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    policy == RuntimePolicy::cuda ? axom::detail::getAllocatorID<axom::MemorySpace::Device>()
    :
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    policy == RuntimePolicy::hip ? axom::detail::getAllocatorID<axom::MemorySpace::Device>()
                                 :
  #endif
                                 axom::INVALID_ALLOCATOR_ID;
#else
  AXOM_UNUSED_VAR(policy);
  int allocatorID = axom::getDefaultAllocatorID();
#endif
  return allocatorID;
}

/// Utility function to initialize the logger
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

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

  conduit::utils::set_error_handler(
    [](auto& msg, auto& file, int line) { slic::logErrorMessage(msg, file, line); });
  conduit::utils::set_warning_handler(
    [](auto& msg, auto& file, int line) { slic::logWarningMessage(msg, file, line); });
  conduit::utils::set_info_handler([](auto& msg, auto& file, int line) {
    slic::logMessage(slic::message::Info, msg, file, line);
  });
}

/// Utility function to finalize the logger
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

/*!
  All the test code that depends on DIM to instantiate.
*/
template <int DIM, typename ExecSpace>
int testNdimInstance(BlueprintStructuredMesh& computationalMesh)
{
  //---------------------------------------------------------------------------
  // params specify which tests to run.
  //---------------------------------------------------------------------------

  std::shared_ptr<PlanarTestStrategy<DIM>> planarStrat;
  std::shared_ptr<RoundTestStrategy<DIM>> roundStrat;
  std::shared_ptr<GyroidTestStrategy<DIM>> gyroidStrat;

  ContourTestBase<DIM, ExecSpace> contourTest;

  if(params.usingPlanar())
  {
    planarStrat = std::make_shared<PlanarTestStrategy<DIM>>(params.planeNormal<DIM>(),
                                                            params.inplanePoint<DIM>());
    contourTest.addTestStrategy(planarStrat);
  }

  if(params.usingRound())
  {
    roundStrat = std::make_shared<RoundTestStrategy<DIM>>(params.roundContourCenter<DIM>());
    roundStrat->setToleranceByLongestEdge(computationalMesh);
    contourTest.addTestStrategy(roundStrat);
  }

  if(params.usingGyroid())
  {
    gyroidStrat =
      std::make_shared<GyroidTestStrategy<DIM>>(params.gyroidScaleFactor<DIM>(), params.contourVal);
    gyroidStrat->setToleranceByLongestEdge(computationalMesh);
    contourTest.addTestStrategy(gyroidStrat);
  }

  contourTest.computeNodalDistance(computationalMesh);

  contourTest.addMaskField(computationalMesh);

  if(params.isVerbose())
  {
    computationalMesh.printMeshInfo();
  }

  // Write computational mesh with contour functions.
  saveMesh(computationalMesh.asConduitNode(), params.fieldsFile);

  int localErrCount = 0;
  localErrCount += contourTest.runTest(computationalMesh);

  // Check results

  int errCount = 0;
  if(params.checkResults)
  {
#ifdef AXOM_USE_MPI
    MPI_Allreduce(&localErrCount, &errCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    errCount = localErrCount;
#endif

    if(errCount)
    {
      SLIC_INFO(axom::fmt::format(" Error exit: {} errors found.", errCount));
    }
    else
    {
      SLIC_INFO(banner("Normal exit."));
    }
  }
  else
  {
    SLIC_INFO("Results not checked.");
  }

  return errCount;
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  myRank = mpi_raii_wrapper.my_rank();
  numRanks = mpi_raii_wrapper.num_ranks();

  initializeLogger();
  //slic::setAbortOnWarning(true);

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  axom::CLI::App app {"Driver/test code for marching cubes algorithm"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(myRank == 0)
    {
      retval = app.exit(e);
    }

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    exit(retval);
  }

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(params.annotationMode);
  axom::utilities::Timer questMarchingCubesExample(false);
  questMarchingCubesExample.start();
  AXOM_ANNOTATE_SCOPE("quest marching cubes example");

  s_allocatorId = allocatorIdToTest(params.policy);

  //---------------------------------------------------------------------------
  // Load computational mesh.
  //---------------------------------------------------------------------------
  AXOM_ANNOTATE_BEGIN("load mesh");
  BlueprintStructuredMesh computationalMesh(params.meshFile, "mesh");
  AXOM_ANNOTATE_END("load mesh");

  SLIC_INFO_IF(params.isVerbose(),
               axom::fmt::format("Computational mesh has {} cells in {} domains locally",
                                 computationalMesh.cellCount(),
                                 computationalMesh.domainCount()));
  slic::flushStreams();

  // Output some global mesh size stats
  {
    int mn, mx, sum;
    getIntMinMax(computationalMesh.cellCount(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format("Computational mesh has {{min:{}, max:{}, sum:{}, avg:{}}} cells",
                                mn,
                                mx,
                                sum,
                                (double)sum / numRanks));
  }
  {
    int mn, mx, sum;
    getIntMinMax(computationalMesh.domainCount(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format("Computational mesh has {{min:{}, max:{}, sum:{}, avg:{}}} domains",
                                mn,
                                mx,
                                sum,
                                (double)sum / numRanks));
  }

  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Run test in the execution space set by command line.
  //---------------------------------------------------------------------------
  int errCount = 0;
  if(params.policy == axom::runtime_policy::Policy::seq)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::SEQ_EXEC>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::SEQ_EXEC>(computationalMesh);
    }
  }
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
  else if(params.policy == axom::runtime_policy::Policy::omp)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::OMP_EXEC>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::OMP_EXEC>(computationalMesh);
    }
  }
  #endif
  #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  else if(params.policy == axom::runtime_policy::Policy::cuda)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::CUDA_EXEC<256>>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::CUDA_EXEC<256>>(computationalMesh);
    }
  }
  #endif
  #if defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  else if(params.policy == axom::runtime_policy::Policy::hip)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::HIP_EXEC<256>>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::HIP_EXEC<256>>(computationalMesh);
    }
  }
  #endif
#endif

  questMarchingCubesExample.stop();
  printTimingStats(questMarchingCubesExample, "questMarchingCubesExample");
  finalizeLogger();

  return errCount != 0;
}
