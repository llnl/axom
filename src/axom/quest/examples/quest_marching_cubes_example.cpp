// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file marching_cubes_example.cpp
 * \brief Driver for Marching Cubes isocontour generation
 *
 * Extracts an isocontour from a scalar field in a Conduit Blueprint mesh.
 */

#include "axom/config.hpp"

// This example requires Conduit and Bump
#ifndef AXOM_USE_CONDUIT
  #error "quest_marching_cubes_example.cpp requires Conduit"
#endif
#ifndef AXOM_USE_BUMP
  #error "quest_marching_cubes_example.cpp requires Bump"
#endif

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/MeshViewUtil.hpp"

#if defined(AXOM_USE_SIDRE)
  #include "axom/sidre.hpp"
#endif

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
#include <cmath>
#include <string>
#include <map>
#include <limits>
#include <memory>
#include <set>
#include <type_traits>
#include <variant>
#include <vector>

namespace quest = axom::quest;
namespace slic = axom::slic;
#if defined(AXOM_USE_SIDRE)
namespace sidre = axom::sidre;
#endif
namespace mint = axom::mint;

using RuntimePolicy = axom::runtime_policy::Policy;

//-----------------------------------------------------------------------------
// converts the input string into an 80 character string
// padded on both sides with '=' symbols
//-----------------------------------------------------------------------------
std::string banner(const std::string& str) { return axom::fmt::format("{:=^80}", str); }

//-----------------------------------------------------------------------------
// Struct to parse and store the input parameters
//-----------------------------------------------------------------------------
struct Input
{
public:
  std::string meshFile;
  std::string fieldName;
  bool listFields {false};
  //! @brief Optional file for Bump's welded Blueprint contour.
  std::string blueprintContourFile {};

  double contourVal {1.0};

  RuntimePolicy policy {RuntimePolicy::seq};

  quest::MarchingCubesDataParallelism dataParallelism = quest::MarchingCubesDataParallelism::byPolicy;

  // Use Bump's CutField backend instead of the legacy backend.
  bool useBumpBackend {false};

  // Number of distinct MarchingCubes objects.
  int objectRepCount {1};
  // Number of contour extractions per MarchingCubes object.
  int contourGenCount {1};
  // Number of masking cycles.
  int maskCount {1};

  std::string annotationMode {"none"};

private:
  bool _verboseOutput {false};

  const std::map<std::string, quest::MarchingCubesDataParallelism> s_validImplChoices {
    {"byPolicy", quest::MarchingCubesDataParallelism::byPolicy},
    {"hybridParallel", quest::MarchingCubesDataParallelism::hybridParallel},
    {"fullParallel", quest::MarchingCubesDataParallelism::fullParallel}};

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-p, --policy", policy)
      ->description("Set the runtime policy for contour extraction")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.add_option("--dataParallelism", dataParallelism)
      ->description(
        "Select the scan mode for the legacy backend "
        "(ignored by --useBumpBackend)")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(s_validImplChoices));

    app.add_flag("--useBumpBackend", useBumpBackend)
      ->description("Use the Bump CutField backend instead of the legacy structured-only backend")
      ->capture_default_str();

    app.add_option("-m,--mesh-file", meshFile)
      ->description("Path to a Conduit Blueprint computational mesh")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_option("-f,--field", fieldName)
      ->description(
        "Name of the vertex-associated scalar field to contour; "
        "required unless --list-fields is used");

    app.add_flag("--list-fields", listFields)
      ->description("List scalar vertex fields and exit")
      ->capture_default_str();

    app.add_option("--blueprint-contour-file", blueprintContourFile)
      ->description("Write Bump's welded contour to a Blueprint file; requires --useBumpBackend")
      ->capture_default_str();

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("--contourVal", contourVal)->description("Contour value")->capture_default_str();

    app.add_option("--objectReps", objectRepCount)
      ->description("Number of setMesh and extraction batches to run")
      ->capture_default_str();

    app.add_option("--contourGenReps", contourGenCount)
      ->description("Number of contour extractions after each setMesh call")
      ->capture_default_str();

    app.add_option("--maskCount", maskCount)
      ->description("Assign cells cyclically to this many mask values")
      ->capture_default_str()
      ->check(axom::CLI::Range(1, std::numeric_limits<int>::max()));

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "Caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    app.get_formatter()->column_width(60);

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

void getIntMinMax(int inVal, int& minVal, int& maxVal, int& sumVal)
{
  minVal = allReduce(inVal, ReductionOperation::Min);
  maxVal = allReduce(inVal, ReductionOperation::Max);
  sumVal = allReduce(inVal, ReductionOperation::Sum);
}

void loadBlueprintMesh(const std::string& meshFilename, conduit::Node& mesh)
{
#ifdef AXOM_USE_MPI
  conduit::relay::mpi::io::blueprint::load_mesh(meshFilename, mesh, MPI_COMM_WORLD);
#else
  conduit::relay::io::blueprint::load_mesh(meshFilename, mesh);
#endif
}

bool verifyBlueprintMesh(const conduit::Node& mesh, conduit::Node& info)
{
#ifdef AXOM_USE_MPI
  return conduit::blueprint::mpi::verify("mesh", mesh, info, MPI_COMM_WORLD);
#else
  return conduit::blueprint::verify("mesh", mesh, info);
#endif
}

int myRank = -1, numRanks = -1;  // MPI stuff, set in main().

/// \brief Host-side access to the example's Blueprint mesh and derived sizes.
struct BlueprintStructuredMesh
{
public:
  explicit BlueprintStructuredMesh(const std::string& meshFile,
                                   const std::string& topologyName,
                                   bool verboseOutput = false)
    : _topologyName(topologyName)
    , _topologyPath("topologies/" + topologyName)
  {
    readBlueprintMesh(meshFile);

    if(verboseOutput)
    {
      for(int d = 0; d < _mdMesh.number_of_children(); ++d)
      {
        if(isStructured(d))
        {
          SLIC_INFO(axom::fmt::format("dom[{}] size={}", d, domainLengths(d)));
        }
        else
        {
          SLIC_INFO(axom::fmt::format("dom[{}] cells={}, nodes={}", d, cellCount(d), nodeCount(d)));
        }
      }
    }
  }

  /// Return the Blueprint mesh.
  conduit::Node& asConduitNode() { return _mdMesh; }

  /// Return the number of local domains.
  axom::IndexType domainCount() const { return _domCount; }

  /// Return whether this rank has no domains.
  bool empty() const { return _domCount == 0; }

  /// Return one local domain.
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
   * @brief Return the logical cell dimensions of a structured Blueprint domain.
   *
   * @param[in] domId Local domain index.
   * @param[out] lengths Buffer for dimension() values.
   */
  void domainLengths(axom::IndexType domId, axom::IndexType* lengths) const
  {
    const conduit::Node& dom = domain(domId);
    SLIC_ASSERT_MSG(isStructured(domId), "domainLengths() is only defined for structured domains.");
    SLIC_ASSERT_MSG(dom.fetch_existing(_coordsetPath + "/type").as_string() == "explicit",
                    axom::fmt::format("Currently only supporting explicit coordinate types."
                                      "  '{}/type' is '{}'",
                                      _coordsetPath,
                                      dom.fetch_existing(_coordsetPath + "/type").as_string()));
    const conduit::Node& dimsNode = dom.fetch_existing(_topologyPath + "/elements/dims");
    for(int i = 0; i < _ndims; ++i)
    {
      lengths[i] = static_cast<axom::IndexType>(dimsNode[i].to_int64());
    }
  }

  axom::Array<axom::IndexType> domainLengths(axom::IndexType domainId) const
  {
    axom::Array<axom::IndexType> rval(_ndims, _ndims);
    domainLengths(domainId, rval.data());
    return rval;
  }

  /// Return the number of cells in a domain.
  int cellCount(axom::IndexType domId) const
  {
    if(isStructured(domId))
    {
      const auto shape = domainLengths(domId);
      int rval = 1;
      for(const auto& l : shape)
      {
        rval *= l;
      }
      return rval;
    }
    return static_cast<int>(
      conduit::blueprint::mesh::topology::length(domain(domId).fetch_existing(_topologyPath)));
  }

  /// Return the number of cells in all local domains.
  int cellCount() const
  {
    int rval = 0;
    for(int domId = 0; domId < _mdMesh.number_of_children(); ++domId)
    {
      rval += cellCount(domId);
    }
    return rval;
  }

  /// Return the number of nodes in a domain.
  int nodeCount(axom::IndexType domId) const
  {
    if(isStructured(domId))
    {
      auto shape = domainLengths(domId);
      int rval = 1;
      for(const auto& l : shape)
      {
        rval *= 1 + l;
      }
      return rval;
    }
    return static_cast<int>(
      conduit::blueprint::mesh::coordset::length(domain(domId).fetch_existing(_coordsetPath)));
  }

  /// Return the number of nodes in all local domains.
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

  std::string topologyType(axom::IndexType domId) const
  {
    return domain(domId).fetch_existing(_topologyPath + "/type").as_string();
  }

  bool isStructured(axom::IndexType domId) const { return topologyType(domId) == "structured"; }

  bool isStridedStructured(axom::IndexType domId) const
  {
    return isStructured(domId) &&
      domain(domId).fetch_existing(_topologyPath + "/elements/dims").has_child("strides");
  }

  /*!
   * @brief Whether this domain uses compact field indexing.
   *
   * Strided structured fields occupy a padded window and require their offsets
   * and strides. Other supported fields use flat node indices.
   */
  bool useFlatFields(axom::IndexType domId) const { return !isStridedStructured(domId); }

  /// Check the Blueprint mesh and print diagnostics when validation fails.
  bool isValid() const
  {
    conduit::Node info;
    if(!verifyBlueprintMesh(_mdMesh, info))
    {
      SLIC_INFO("Invalid blueprint for mesh: \n" << info.to_yaml());
      slic::flushStreams();
      return false;
    }
    return true;
  }

  void printMeshInfo() const { _mdMesh.print(); }

  template <typename ExecSpace>
  void copyMeshToMemorySpace(int allocId = axom::execution_space<ExecSpace>::allocatorID())
  {
    AXOM_ANNOTATE_SCOPE("copyMeshToMemorySpace");
    conduit::Node newMesh;
    axom::bump::utilities::copy<ExecSpace>(newMesh, _mdMesh, allocId);
    _mdMesh.swap(newMesh);
  }

private:
  int _ndims {-1};
  conduit::Node _mdMesh;
  axom::IndexType _domCount;
  const std::string _topologyName;
  const std::string _topologyPath;
  std::string _coordsetPath;

  axom::IndexType dimValue(const conduit::Node& node, int dim, axom::IndexType defaultValue = 0) const
  {
    static const char* dimNames[] = {"i", "j", "k"};
    if(node.has_child(dimNames[dim]))
    {
      return static_cast<axom::IndexType>(node.fetch_existing(dimNames[dim]).to_int64());
    }
    if(node.dtype().is_int32())
    {
      return static_cast<axom::IndexType>(node.as_int32_ptr()[dim]);
    }
    if(node.dtype().is_int64())
    {
      return static_cast<axom::IndexType>(node.as_int64_ptr()[dim]);
    }
    if(dim < node.number_of_children())
    {
      return static_cast<axom::IndexType>(node[dim].to_int64());
    }
    return defaultValue;
  }

  //! @brief Read a Blueprint mesh and normalize it to a multi-domain node.
  void readBlueprintMesh(const std::string& meshFilename)
  {
    SLIC_ASSERT(!meshFilename.empty());

    conduit::Node loadedMesh;
    loadBlueprintMesh(meshFilename, loadedMesh);
    // Normalize to a multi-domain node. MarchingCubes::setMesh() performs the
    // equivalent normalization for its input. This wrapper needs its own copy
    // because its size and coordset helpers operate independently of MarchingCubes.
    _mdMesh.reset();
    if(conduit::blueprint::mesh::is_multi_domain(loadedMesh))
    {
      _mdMesh.swap(loadedMesh);
    }
    else
    {
      _mdMesh.append().set(loadedMesh);
    }
    _domCount = conduit::blueprint::mesh::number_of_domains(_mdMesh);

    if(_domCount > 0)
    {
      SLIC_ASSERT(_mdMesh[0].has_path(_topologyPath));
      auto coordsetName = _mdMesh[0].fetch_existing(_topologyPath + "/coordset").as_string();
      _coordsetPath = axom::fmt::format("coordsets/{}", coordsetName);
      SLIC_ASSERT(_mdMesh[0].has_path(_coordsetPath));

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
  // lambda that sends rank 0's fields to the other ranks and checks that they're all present
  // It's a no-op when not using MPI
  const auto fieldNamesAgree = [](const std::vector<std::string>& fieldNames) {
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

    return allReduce(fieldNames == rootFieldNames ? 1 : 0, ReductionOperation::LogicalAnd) != 0;
#else
    AXOM_UNUSED_VAR(fieldNames);
    return true;
#endif
  };

  // collect the sorted list of field names on this rank
  std::vector<std::string> fieldNames;
  if(!mesh.empty() && mesh.domain(0).has_child("fields"))
  {
    fieldNames = mesh.domain(0).fetch_existing("fields").child_names();
    std::sort(fieldNames.begin(), fieldNames.end());
  }

  SLIC_ERROR_IF(!fieldNamesAgree(fieldNames),
                "Field names must agree with the field names on MPI rank 0.");
  return fieldNames;
}

// A one-dimensional bounding box represents the minimum and maximum of a scalar range.
using ScalarRange = axom::primal::BoundingBox<double, 1>;

struct FieldSummary
{
  int presentDomains {0};
  int vertexDomains {0};
  int numericDomains {0};
  int float64Domains {0};
  int topologyDomains {0};
  int hasNan {0};
  ScalarRange valueRange;
  std::set<std::string> topologies;
  std::set<std::string> types;
};

FieldSummary summarizeLocalField(const BlueprintStructuredMesh& mesh,
                                 const std::string& fieldName,
                                 const std::string& topologyName)
{
  FieldSummary summary;
  for(axom::IndexType domainIdx = 0; domainIdx < mesh.domainCount(); ++domainIdx)
  {
    const conduit::Node& domain = mesh.domain(domainIdx);
    if(!domain.has_child("fields") || !domain.fetch_existing("fields").has_child(fieldName))
    {
      continue;
    }

    ++summary.presentDomains;
    const conduit::Node& field = domain.fetch_existing("fields").child(fieldName);
    const std::string association =
      field.has_child("association") ? field.fetch_existing("association").as_string() : "<missing>";
    const std::string topology =
      field.has_child("topology") ? field.fetch_existing("topology").as_string() : "<missing>";
    summary.topologies.insert(topology);
    summary.topologyDomains += topology == topologyName ? 1 : 0;

    if(association != "vertex")
    {
      continue;
    }
    ++summary.vertexDomains;

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
    summary.float64Domains += values.dtype().is_float64() ? 1 : 0;
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

void printFieldSummary(const BlueprintStructuredMesh& mesh, const std::string& topologyName)
{
  const auto fieldNames = getFieldNames(mesh);
  const int globalDomainCount =
    allReduce(static_cast<int>(mesh.domainCount()), ReductionOperation::Sum);

  std::string output = axom::fmt::format("Scalar vertex fields across {} domain{}:\n",
                                         globalDomainCount,
                                         globalDomainCount == 1 ? "" : "s");
  output += "Name | Topology | Type | Domain coverage | Global range | Marching Cubes\n";

  int listedFieldCount = 0;
  for(const auto& fieldName : fieldNames)
  {
    FieldSummary local = summarizeLocalField(mesh, fieldName, topologyName);
    FieldSummary global;
    global.presentDomains = allReduce(local.presentDomains, ReductionOperation::Sum);
    global.vertexDomains = allReduce(local.vertexDomains, ReductionOperation::Sum);
    global.numericDomains = allReduce(local.numericDomains, ReductionOperation::Sum);
    global.float64Domains = allReduce(local.float64Domains, ReductionOperation::Sum);
    global.topologyDomains = allReduce(local.topologyDomains, ReductionOperation::Sum);
    global.hasNan = allReduce(local.hasNan, ReductionOperation::Max);

    const int hasRange = allReduce(local.valueRange.isValid() ? 1 : 0, ReductionOperation::Max);
    if(hasRange != 0)
    {
      global.valueRange.addPoint(
        ScalarRange::PointType {allReduce(local.valueRange.getMin()[0], ReductionOperation::Min)});
      global.valueRange.addPoint(
        ScalarRange::PointType {allReduce(local.valueRange.getMax()[0], ReductionOperation::Max)});
    }

    const auto& topologies = local.topologies;
    const auto& types = local.types;
    if(global.numericDomains == 0)
    {
      continue;
    }

    ++listedFieldCount;
    std::string range = "empty";
    if(global.valueRange.isValid())
    {
      range = axom::fmt::format("[{:.17g}, {:.17g}]",
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

    std::vector<std::string> incompatibilities;
    if(global.presentDomains != globalDomainCount)
    {
      incompatibilities.emplace_back("missing domains");
    }
    if(global.vertexDomains != global.presentDomains)
    {
      incompatibilities.emplace_back("association must be vertex");
    }
    if(global.numericDomains != global.vertexDomains)
    {
      incompatibilities.emplace_back("values must be numeric scalars");
    }
    if(global.float64Domains != global.numericDomains)
    {
      incompatibilities.emplace_back("type must be float64");
    }
    if(global.topologyDomains != global.presentDomains)
    {
      incompatibilities.emplace_back("topology must be " + topologyName);
    }

    const std::string compatibility = incompatibilities.empty()
      ? "compatible"
      : "incompatible: " + axom::fmt::format("{}", axom::fmt::join(incompatibilities, "; "));
    output += axom::fmt::format("{} | {} | {} | {}/{} | {} | {}\n",
                                fieldName,
                                axom::fmt::join(topologies, ", "),
                                axom::fmt::join(types, ", "),
                                global.presentDomains,
                                globalDomainCount,
                                range,
                                compatibility);
  }

  if(listedFieldCount == 0)
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
    const double minCompute = allReduce(t.elapsedTimeInSec(), ReductionOperation::Min);
    const double maxCompute = allReduce(t.elapsedTimeInSec(), ReductionOperation::Max);
    const double sumCompute = allReduce(t.elapsedTimeInSec(), ReductionOperation::Sum);

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
    if(!verifyBlueprintMesh(tmpMesh, info))
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
struct ContourTestBase
{
  explicit ContourTestBase(const Input& params)
    : m_params(params)
    , m_parentCellIdField("parentCellIds")
    , m_domainIdField("domainIdField")
  { }

  const Input& m_params;
  const std::string m_parentCellIdField;
  const std::string m_domainIdField;

  int runTest(BlueprintStructuredMesh& computationalMesh)
  {
    AXOM_ANNOTATE_SCOPE("runTest");

    // Conduit data is in host memory, move to devices for testing.
    if(s_allocatorId != axom::execution_space<axom::SEQ_EXEC>::allocatorID())
    {
      AXOM_ANNOTATE_SCOPE("move mesh to device memory");

      computationalMesh.template copyMeshToMemorySpace<ExecSpace>(s_allocatorId);
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
      const std::string dataPath = axom::fmt::format("fields/{}/values", m_params.fieldName);
      void* dataPtr = computationalMesh.domain(0).fetch_existing(dataPath).data_ptr();
      bool dataFromUmpire = rm.hasAllocator(dataPtr);
      if(dataFromUmpire)
      {
        umpire::Allocator allocator = rm.getAllocator(dataPtr);
        resourceName = allocator.getName();
      }
      SLIC_INFO(axom::fmt::format("Testing with policy {} and function data on {}",
                                  m_params.policy,
                                  resourceName));
      if(m_params.policy == axom::runtime_policy::Policy::seq)
      {
        SLIC_ASSERT(resourceName == "HOST");
      }
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      else if(m_params.policy == axom::runtime_policy::Policy::omp)
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

    // objectRepCount setMesh calls
    axom::utilities::Timer setMeshTimer(false);

    // Time steady-state computeIsocontour calls
    axom::utilities::Timer contourTimer(false);

    // Time first contourIsocontour call after setMesh
    axom::utilities::Timer contourTimerM(false);

    std::unique_ptr<quest::MarchingCubes> mcPtr;
    const auto objectLoopName = axom::fmt::format("objectRepLoop {}", m_params.objectRepCount);
    AXOM_ANNOTATE_BEGIN(objectLoopName);
    objectRepLoopTimer.start();
    for(int j = 0; j < m_params.objectRepCount; ++j)
    {
      if(!mcPtr)
      {
        AXOM_ANNOTATE_SCOPE("MCInit");
        initializationTimer.start();
        mcPtr = std::make_unique<quest::MarchingCubes>(m_params.policy,
                                                       s_allocatorId,
                                                       m_params.dataParallelism);
        mcPtr->setUseBumpBackend(m_params.useBumpBackend);
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
      for(int i = 0; i < m_params.contourGenCount; ++i)
      {
        SLIC_DEBUG(axom::fmt::format("MarchingCubes object rep {} of {}, contour run {} of {}:",
                                     j,
                                     m_params.objectRepCount,
                                     i,
                                     m_params.contourGenCount));
        mc.clearOutput();
        mc.setFunctionField(m_params.fieldName);
        for(int iMask = 0; iMask < m_params.maskCount; ++iMask)
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
          mc.computeIsocontour(m_params.contourVal);
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
                                m_params.objectRepCount,
                                m_params.contourGenCount));
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

      computationalMesh.template copyMeshToMemorySpace<axom::SEQ_EXEC>(
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

    // Optionally write Bump's welded Blueprint contour.
    if(!m_params.blueprintContourFile.empty())
    {
      if(!m_params.useBumpBackend)
      {
        SLIC_WARNING(
          "--blueprint-contour-file requires --useBumpBackend; the legacy kernel has no "
          "Blueprint contour output. Skipping.");
      }
      else
      {
        AXOM_ANNOTATE_SCOPE("write blueprint contour");
        conduit::Node contourBp;
        mc.populateContourMeshBlueprint(contourBp);
        SLIC_INFO(axom::fmt::format("Blueprint contour has {} domains; writing to '{}'",
                                    contourBp.number_of_children(),
                                    m_params.blueprintContourFile));
        saveMesh(contourBp, m_params.blueprintContourFile);
      }
    }

    {
      axom::Array<axom::IndexType, 2> facetNodeIds;
      axom::Array<double, 2> facetNodeCoords;
      axom::Array<axom::IndexType, 1> facetParentIds;
      axom::Array<axom::IndexType> facetDomainIds;
      mc.relinquishContourData(facetNodeIds, facetNodeCoords, facetParentIds, facetDomainIds);
      SLIC_ASSERT(mc.getContourFacetCount() == 0);
    }
    AXOM_ANNOTATE_END("convert to mint mesh");

    // main() reduces this value so all ranks return the same exit code.
    const int localErrCount = 0;

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
    SLIC_INFO_IF(m_params.isVerbose(),
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
      if(bpMesh.useFlatFields(domId))
      {
        addMaskFieldFlat(bpMesh.domain(domId));
        continue;
      }

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
      int maskCount = m_params.maskCount;
      axom::for_all<axom::SEQ_EXEC>(
        0,
        cellCount,
        AXOM_LAMBDA(axom::IndexType cellId) { maskView.flatIndex(cellId) = (cellId % maskCount); });
    }
  }

  void addMaskFieldFlat(conduit::Node& dom)
  {
    const axom::IndexType cellCount = static_cast<axom::IndexType>(
      conduit::blueprint::mesh::topology::length(dom.fetch_existing("topologies/mesh")));

    conduit::Node& mask = dom["fields/mask"];
    mask["association"] = "element";
    mask["topology"] = "mesh";
    mask["values"].set(conduit::DataType::c_int(cellCount));
    auto* maskValues = mask["values"].as_int_ptr();

    const int maskCount = m_params.maskCount;
    for(axom::IndexType cellId = 0; cellId < cellCount; ++cellId)
    {
      maskValues[cellId] = static_cast<int>(cellId % maskCount);
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

// ----------------------------------------------------------------------------
// Utility RAII struct to set up and tear down the example's logger
// ----------------------------------------------------------------------------
struct ParallelLoggerRAII
{
  ParallelLoggerRAII()
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

  void flush() { slic::flushStreams(); }

  /// Utility function to finalize the logger
  ~ParallelLoggerRAII()
  {
    if(slic::isInitialized())
    {
      slic::flushStreams();
      slic::finalize();
    }
  }
};

// ----------------------------------------------------------------------------
// Tag dispatch for choosing the desired execution policy and dimension
// ----------------------------------------------------------------------------
template <typename T>
struct TypeTag
{
  using type = T;
};

template <int DIM_, typename ExecSpace_>
struct TestInstance
{
  static constexpr int DIM = DIM_;
  using ExecSpace = ExecSpace_;
};

using TestInstanceVariant = std::variant<TestInstance<2, axom::SEQ_EXEC>,
                                         TestInstance<3, axom::SEQ_EXEC>
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
                                         ,
                                         TestInstance<2, axom::OMP_EXEC>,
                                         TestInstance<3, axom::OMP_EXEC>
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
                                         ,
                                         TestInstance<2, axom::CUDA_EXEC<256>>,
                                         TestInstance<3, axom::CUDA_EXEC<256>>
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
                                         ,
                                         TestInstance<2, axom::HIP_EXEC<256>>,
                                         TestInstance<3, axom::HIP_EXEC<256>>
#endif
                                         >;

template <typename ExecSpace>
TestInstanceVariant selectTestDimension(TypeTag<ExecSpace>, int dimension)
{
  if(dimension == 2)
  {
    return TestInstance<2, ExecSpace> {};
  }
  if(dimension == 3)
  {
    return TestInstance<3, ExecSpace> {};
  }

  SLIC_ERROR(axom::fmt::format("Unsupported mesh dimension {}", dimension));
  return TestInstance<2, axom::SEQ_EXEC> {};
}

TestInstanceVariant selectTestInstance(const Input& params, int dimension)
{
  if(params.policy == RuntimePolicy::seq)
  {
    return selectTestDimension(TypeTag<axom::SEQ_EXEC> {}, dimension);
  }
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  if(params.policy == RuntimePolicy::omp)
  {
    return selectTestDimension(TypeTag<axom::OMP_EXEC> {}, dimension);
  }
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(params.policy == RuntimePolicy::cuda)
  {
    return selectTestDimension(TypeTag<axom::CUDA_EXEC<256>> {}, dimension);
  }
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  if(params.policy == RuntimePolicy::hip)
  {
    return selectTestDimension(TypeTag<axom::HIP_EXEC<256>> {}, dimension);
  }
#endif

  SLIC_ERROR(axom::fmt::format("Unsupported runtime policy {}", params.policy));
  return TestInstance<2, axom::SEQ_EXEC> {};
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  myRank = mpi_raii_wrapper.my_rank();
  numRanks = mpi_raii_wrapper.num_ranks();

  ParallelLoggerRAII raii_logger;
  //slic::setAbortOnWarning(true);

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  axom::CLI::App app {"Driver/test code for marching cubes algorithm"};
  Input params;

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
  BlueprintStructuredMesh computationalMesh(params.meshFile, "mesh", params.isVerbose());
  AXOM_ANNOTATE_END("load mesh");

  if(params.listFields)
  {
    printFieldSummary(computationalMesh, "mesh");
    raii_logger.flush();
    return 0;
  }

  SLIC_INFO_IF(params.isVerbose(),
               axom::fmt::format("Computational mesh has {} cells in {} domains locally",
                                 computationalMesh.cellCount(),
                                 computationalMesh.domainCount()));
  raii_logger.flush();

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

  raii_logger.flush();

  //---------------------------------------------------------------------------
  // Run test in the execution space set by command line.
  //---------------------------------------------------------------------------
  auto testInstance = selectTestInstance(params, computationalMesh.dimension());
  int errCount = std::visit(
    [&](const auto& instance) {
      AXOM_UNUSED_VAR(instance);
      using Instance = std::decay_t<decltype(instance)>;
      constexpr int DIM = Instance::DIM;
      using ExecSpace = typename Instance::ExecSpace;

      ContourTestBase<DIM, ExecSpace> contourTest(params);
      contourTest.addMaskField(computationalMesh);

      if(params.isVerbose())
      {
        computationalMesh.printMeshInfo();
      }

      int localErrCount = contourTest.runTest(computationalMesh);

      const int globalErrCount = allReduce(localErrCount, ReductionOperation::Sum);

      if(globalErrCount)
      {
        SLIC_INFO(axom::fmt::format(" Error exit: {} errors found.", globalErrCount));
      }
      else
      {
        SLIC_INFO(banner("Normal exit."));
      }

      return globalErrCount;
    },
    testInstance);

  questMarchingCubesExample.stop();
  printTimingStats(questMarchingCubesExample, "questMarchingCubesExample");

  return errCount != 0;
}
