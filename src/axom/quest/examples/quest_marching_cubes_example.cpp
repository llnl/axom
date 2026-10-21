// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 \file marching_cubes_example.cpp
 \brief Extract an isocontour from a scalar field in a Conduit Blueprint mesh
*/

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesFullParallel.hpp requires conduit"
#endif

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/MeshViewUtil.hpp"
#if defined(AXOM_USE_SIDRE)
  #include "axom/sidre.hpp"
#endif
#include "axom/core/Types.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_relay_io_blueprint.hpp"
#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
  #include "conduit_relay_mpi.hpp"
  #include "conduit_relay_mpi_io_blueprint.hpp"
#endif

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// C/C++ includes
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <map>
#include <limits>
#include <memory>
#include <set>
#include <vector>

namespace quest = axom::quest;
namespace slic = axom::slic;
#if defined(AXOM_USE_SIDRE)
namespace sidre = axom::sidre;
#endif
namespace mint = axom::mint;

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
  std::string fieldName;
  bool listFields {false};

  double contourVal {1.0};

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
      ->description("Path to a Conduit Blueprint computational mesh")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_option("-f,--field", fieldName)
      ->description(
        "Name of the vertex-associated scalar field to contour; "
        "required unless --list-fields is used");

    app.add_flag("--list-fields", listFields)
      ->description("List scalar fields, their centering, and ranges, then exit")
      ->capture_default_str();

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("--contourVal", contourVal)->description("Contour value")->capture_default_str();

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

    if(fieldName.empty() && !listFields)
    {
      throw axom::CLI::RequiredError("--field");
    }

    slic::setLoggingMsgLevel(_verboseOutput ? slic::message::Debug : slic::message::Info);
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

namespace
{

enum class ReductionOperation
{
  Min,
  Max,
  Sum,
  LogicalAnd
};

template <typename T>
T allReduce(T value, ReductionOperation operation)
{
#ifdef AXOM_USE_MPI
  MPI_Op mpiOperation = MPI_OP_NULL;
  switch(operation)
  {
  case ReductionOperation::Min:
    mpiOperation = MPI_MIN;
    break;
  case ReductionOperation::Max:
    mpiOperation = MPI_MAX;
    break;
  case ReductionOperation::Sum:
    mpiOperation = MPI_SUM;
    break;
  case ReductionOperation::LogicalAnd:
    mpiOperation = MPI_LAND;
    break;
  }
  SLIC_ASSERT(mpiOperation != MPI_OP_NULL);

  T result {};
  MPI_Allreduce(&value, &result, 1, axom::mpi_traits<T>::type, mpiOperation, MPI_COMM_WORLD);
  return result;
#else
  AXOM_UNUSED_VAR(operation);
  return value;
#endif
}

}  // namespace

template <typename T>
void getMinMaxSum(T inVal, T& minVal, T& maxVal, T& sumVal)
{
  minVal = allReduce(inVal, ReductionOperation::Min);
  maxVal = allReduce(inVal, ReductionOperation::Max);
  sumVal = allReduce(inVal, ReductionOperation::Sum);
}

Input params;

int myRank = -1, numRanks = -1;  // MPI stuff, set in main().

/// Generic computational mesh, to hold cell and node data.
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
   * @brief Get the number of cells in each direction of a blueprint single domain.
   *
   * @param domId Index of domain
   * @lengths Space for dimension() numbers.
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
    _ndims = allReduce(_ndims, ReductionOperation::Max);
    SLIC_ASSERT(_ndims > 0);

    SLIC_ASSERT(isValid());
  }
};  // BlueprintStructuredMesh

namespace
{

std::vector<std::string> getFieldNames(const BlueprintStructuredMesh& mesh)
{
  std::vector<std::string> fieldNames;
  if(!mesh.empty() && mesh.domain(0).has_child("fields"))
  {
    fieldNames = mesh.domain(0).fetch_existing("fields").child_names();
    std::sort(fieldNames.begin(), fieldNames.end());
  }

#ifdef AXOM_USE_MPI
  conduit::Node rootFieldNamesNode;
  if(myRank == 0)
  {
    rootFieldNamesNode.set(conduit::DataType::list());
    for(const auto& fieldName : fieldNames)
    {
      rootFieldNamesNode.append() = fieldName;
    }
  }
  conduit::relay::mpi::broadcast_using_schema(rootFieldNamesNode, 0, MPI_COMM_WORLD);

  std::vector<std::string> rootFieldNames;
  rootFieldNames.reserve(rootFieldNamesNode.number_of_children());
  for(const conduit::Node& fieldName : rootFieldNamesNode.children())
  {
    rootFieldNames.push_back(fieldName.as_string());
  }
  const bool namesAgree =
    allReduce(fieldNames == rootFieldNames ? 1 : 0, ReductionOperation::LogicalAnd) != 0;
  SLIC_ERROR_IF(!namesAgree, "Field names must agree with the field names on MPI rank 0.");
#endif
  return fieldNames;
}

using ScalarRange = axom::primal::BoundingBox<double, 1>;

struct FieldSummary
{
  int numericDomains {0};
  int hasNan {0};
  ScalarRange valueRange;
  std::set<std::string> topologies;
  std::set<std::string> types;
  std::set<std::string> centerings;
};

FieldSummary summarizeLocalField(const BlueprintStructuredMesh& mesh, const std::string& fieldName)
{
  FieldSummary summary;
  for(axom::IndexType domainIdx = 0; domainIdx < mesh.domainCount(); ++domainIdx)
  {
    const conduit::Node& domain = mesh.domain(domainIdx);
    if(!domain.has_child("fields") || !domain.fetch_existing("fields").has_child(fieldName))
    {
      continue;
    }

    const conduit::Node& field = domain.fetch_existing("fields").child(fieldName);
    const std::string association =
      field.has_child("association") ? field.fetch_existing("association").as_string() : "<missing>";
    const std::string topology =
      field.has_child("topology") ? field.fetch_existing("topology").as_string() : "<missing>";
    summary.topologies.insert(topology);
    summary.centerings.insert(association);

    if(!field.has_child("values"))
    {
      summary.types.insert("<missing>");
      continue;
    }

    const conduit::Node& values = field.fetch_existing("values");
    summary.types.insert(values.dtype().name());
    if(!values.dtype().is_number())
    {
      continue;
    }

    ++summary.numericDomains;
    const auto accessor = values.as_double_accessor();
    for(conduit::index_t valueIdx = 0; valueIdx < accessor.number_of_elements(); ++valueIdx)
    {
      const double value = accessor[valueIdx];
      if(std::isnan(value))
      {
        summary.hasNan = 1;
      }
      else
      {
        summary.valueRange.addPoint(ScalarRange::PointType {value});
      }
    }
  }
  return summary;
}

void printFieldSummary(const BlueprintStructuredMesh& mesh)
{
  const auto fieldNames = getFieldNames(mesh);
  const int globalDomainCount =
    allReduce(static_cast<int>(mesh.domainCount()), ReductionOperation::Sum);

  std::string output = axom::fmt::format("Scalar fields across {} domain{}:\n",
                                         globalDomainCount,
                                         globalDomainCount == 1 ? "" : "s");
  using TableRow = std::array<std::string, 5>;
  std::vector<TableRow> rows {{{"Name", "Topology", "Type", "Centering", "Global range"}}};

  for(const auto& fieldName : fieldNames)
  {
    const FieldSummary local = summarizeLocalField(mesh, fieldName);
    FieldSummary global;
    global.numericDomains = allReduce(local.numericDomains, ReductionOperation::Sum);
    global.hasNan = allReduce(local.hasNan, ReductionOperation::Max);

    const int hasRange = allReduce(local.valueRange.isValid() ? 1 : 0, ReductionOperation::Max);
    if(hasRange != 0)
    {
      global.valueRange.addPoint(
        ScalarRange::PointType {allReduce(local.valueRange.getMin()[0], ReductionOperation::Min)});
      global.valueRange.addPoint(
        ScalarRange::PointType {allReduce(local.valueRange.getMax()[0], ReductionOperation::Max)});
    }

    if(global.numericDomains == 0)
    {
      continue;
    }

    std::string range = "empty";
    if(global.valueRange.isValid())
    {
      range = axom::fmt::format("[{:.3g}, {:.3g}]",
                                global.valueRange.getMin()[0],
                                global.valueRange.getMax()[0]);
      if(global.hasNan != 0)
      {
        range += " (contains NaN)";
      }
    }
    else if(global.hasNan != 0)
    {
      range = "NaN only";
    }

    rows.push_back({{fieldName,
                     axom::fmt::format("{}", axom::fmt::join(local.topologies, ", ")),
                     axom::fmt::format("{}", axom::fmt::join(local.types, ", ")),
                     axom::fmt::format("{}", axom::fmt::join(local.centerings, ", ")),
                     range}});
  }

  std::array<std::size_t, 5> widths {};
  for(const auto& row : rows)
  {
    for(std::size_t column = 0; column < row.size(); ++column)
    {
      widths[column] = std::max(widths[column], row[column].size());
    }
  }
  for(const auto& row : rows)
  {
    for(std::size_t column = 0; column < row.size() - 1; ++column)
    {
      output += axom::fmt::format("{:<{}} | ", row[column], widths[column]);
    }
    output += row.back() + "\n";
  }

  if(rows.size() == 1)
  {
    output += "(none)\n";
  }

  if(myRank == 0)
  {
    SLIC_INFO(output);
  }
}

}  // namespace

/// Output some timing stats
void printTimingStats(axom::utilities::Timer& t, const std::string& description)
{
  {
    double minCompute, maxCompute, sumCompute;
    getMinMaxSum(t.elapsedTimeInSec(), minCompute, maxCompute, sumCompute);

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

template <int DIM, typename ExecSpace>
struct ContourRunner
{
  ContourRunner() : m_parentCellIdField("parentCellIds"), m_domainIdField("domainIdField") { }

  const std::string m_parentCellIdField;
  const std::string m_domainIdField;

  void run(BlueprintStructuredMesh& computationalMesh)
  {
    AXOM_ANNOTATE_SCOPE("run contour");

    // Conduit data is in host memory. Move the contour input to device memory when needed.
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
      computationalMesh.moveMeshDataToNewMemorySpace<double>(
        axom::fmt::format("fields/{}/values", params.fieldName),
        s_allocatorId);
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
      const std::string dataPath = axom::fmt::format("fields/{}/values", params.fieldName);
      void* dataPtr = computationalMesh.domain(0).fetch_existing(dataPath).data_ptr();
      bool dataFromUmpire = rm.hasAllocator(dataPtr);
      if(dataFromUmpire)
      {
        umpire::Allocator allocator = rm.getAllocator(dataPtr);
        resourceName = allocator.getName();
      }
      SLIC_INFO(axom::fmt::format("Running with policy {} and field '{}' on {}",
                                  params.policy,
                                  params.fieldName,
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
        mc.setFunctionField(params.fieldName);
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
      computationalMesh.moveMeshDataToNewMemorySpace<double>(
        axom::fmt::format("fields/{}/values", params.fieldName),
        axom::execution_space<axom::SEQ_EXEC>::allocatorID());
      computationalMesh.moveMeshDataToNewMemorySpace<int>(
        axom::fmt::format("fields/{}/values", "mask"),
        axom::execution_space<axom::SEQ_EXEC>::allocatorID());
    }

    // Put contour mesh in a mint object for output.
    AXOM_ANNOTATE_BEGIN("contour output");

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

    AXOM_ANNOTATE_END("contour output");
  }

  void printRunStats(const quest::MarchingCubes& mc)
  {
    {
      axom::IndexType mn, mx, sum;
      getMinMaxSum(mc.getContourCellCount(), mn, mx, sum);
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
      axom::for_all<axom::SEQ_EXEC>(0, cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
        maskView.flatIndex(cellId) = (cellId % maskCount);
      });
    }
  }
};

///
int allocatorIdToTest(axom::runtime_policy::Policy policy)
{
#if defined(AXOM_USE_UMPIRE)
  //---------------------------------------------------------------------------
  // Memory resource.  For testing, choose device memory if appropriate.
  //---------------------------------------------------------------------------
  int allocatorID = policy == RuntimePolicy::seq ? axom::detail::getDefaultHostAllocatorID() :
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    policy == RuntimePolicy::omp ? axom::detail::getDefaultHostAllocatorID()
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

//! Run contour extraction for one dimension and execution space.
template <int DIM, typename ExecSpace>
void runNdimInstance(BlueprintStructuredMesh& computationalMesh)
{
  ContourRunner<DIM, ExecSpace> contour;
  contour.addMaskField(computationalMesh);

  if(params.isVerbose())
  {
    computationalMesh.printMeshInfo();
  }

  contour.run(computationalMesh);
  SLIC_INFO(banner("Normal exit."));
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
  axom::CLI::App app {"Marching Cubes isocontour example"};

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

    finalizeLogger();
    return retval;
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

  if(params.listFields)
  {
    printFieldSummary(computationalMesh);
    finalizeLogger();
    return 0;
  }

  const int dimension = computationalMesh.dimension();

  SLIC_INFO_IF(params.isVerbose(),
               axom::fmt::format("Computational mesh has {} cells in {} domains locally",
                                 computationalMesh.cellCount(),
                                 computationalMesh.domainCount()));
  slic::flushStreams();

  // Output some global mesh size stats
  {
    int mn, mx, sum;
    getMinMaxSum(computationalMesh.cellCount(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format("Computational mesh has {{min:{}, max:{}, sum:{}, avg:{}}} cells",
                                mn,
                                mx,
                                sum,
                                (double)sum / numRanks));
  }
  {
    axom::IndexType mn, mx, sum;
    getMinMaxSum(computationalMesh.domainCount(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format("Computational mesh has {{min:{}, max:{}, sum:{}, avg:{}}} domains",
                                mn,
                                mx,
                                sum,
                                (double)sum / numRanks));
  }

  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Run contour extraction in the selected execution space.
  //---------------------------------------------------------------------------
  if(params.policy == axom::runtime_policy::Policy::seq)
  {
    if(dimension == 2)
    {
      runNdimInstance<2, axom::SEQ_EXEC>(computationalMesh);
    }
    else if(dimension == 3)
    {
      runNdimInstance<3, axom::SEQ_EXEC>(computationalMesh);
    }
  }
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
  else if(params.policy == axom::runtime_policy::Policy::omp)
  {
    if(dimension == 2)
    {
      runNdimInstance<2, axom::OMP_EXEC>(computationalMesh);
    }
    else if(dimension == 3)
    {
      runNdimInstance<3, axom::OMP_EXEC>(computationalMesh);
    }
  }
  #endif
  #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  else if(params.policy == axom::runtime_policy::Policy::cuda)
  {
    if(dimension == 2)
    {
      runNdimInstance<2, axom::CUDA_EXEC<256>>(computationalMesh);
    }
    else if(dimension == 3)
    {
      runNdimInstance<3, axom::CUDA_EXEC<256>>(computationalMesh);
    }
  }
  #endif
  #if defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  else if(params.policy == axom::runtime_policy::Policy::hip)
  {
    if(dimension == 2)
    {
      runNdimInstance<2, axom::HIP_EXEC<256>>(computationalMesh);
    }
    else if(dimension == 3)
    {
      runNdimInstance<3, axom::HIP_EXEC<256>>(computationalMesh);
    }
  }
  #endif
#endif

  questMarchingCubesExample.stop();
  printTimingStats(questMarchingCubesExample, "questMarchingCubesExample");
  finalizeLogger();

  return 0;
}
