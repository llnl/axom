// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file marching_cubes_example.cpp
 * \brief Driver for Marching Cubes isocontour generation
 *
 * Generates planar, round, and gyroid scalar fields and extracts their contours.
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
#include "axom/primal.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/Shapes.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/MeshViewUtil.hpp"

#if defined(AXOM_USE_SIDRE)
  #include "axom/sidre.hpp"
#endif

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
#include <limits>
#include <functional>
#include <memory>
#include <type_traits>
#include <variant>

namespace quest = axom::quest;
namespace slic = axom::slic;
#if defined(AXOM_USE_SIDRE)
namespace sidre = axom::sidre;
#endif
namespace primal = axom::primal;
namespace bumpviews = axom::bump::views;
namespace mint = axom::mint;
namespace numerics = axom::numerics;

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
  std::string fieldsFile {"fields"};
  //! @brief Optional file for Bump's welded Blueprint contour.
  std::string blueprintContourFile {};

  // Center of round contour function
  std::vector<double> fcnCenter;

  // Scaling factor for gyroid function
  std::vector<double> gyroidScale;

  // Parameters for planar contour function
  std::vector<double> inPlane;
  std::vector<double> perpDir;

  std::size_t ndim {0};

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
      ->check(axom::CLI::ExistingFile);

    app.add_option("--blueprint-contour-file", blueprintContourFile)
      ->description("Write Bump's welded contour to a Blueprint file; requires --useBumpBackend")
      ->capture_default_str();

    app.add_option("-s,--fields-file", fieldsFile)
      ->description("Write the input mesh and generated fields to this file")
      ->capture_default_str();

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    auto* distanceFunctionOption =
      app.add_option_group("distanceFunctionOption", "Scalar-field options");

    auto* distFromPtOption =
      distanceFunctionOption->add_option_group("distFromPtOption", "Distance from a point");
    distFromPtOption->add_option("--center", fcnCenter)
      ->description("Center of the circle or sphere field (x,y[,z])")
      ->expected(2, 3);

    auto* gyroidOption = distanceFunctionOption->add_option_group("gyroidOption", "Gyroid field");
    gyroidOption->add_option("--scale", gyroidScale)
      ->description("Gyroid scale factors (x,y[,z])")
      ->expected(2, 3);

    auto* distFromPlaneOption =
      distanceFunctionOption->add_option_group("distFromPlaneOption", "Distance from a plane");
    auto* perpDirOption = distFromPlaneOption->add_option("--dir", perpDir)
                            ->description("Plane normal direction (x,y[,z])")
                            ->expected(2, 3);
    distFromPlaneOption->add_option("--inPlane", inPlane)
      ->description("Point on the plane (x,y[,z])")
      ->expected(2, 3)
      ->needs(perpDirOption);

    // Require at least one distance function, and allow all three.
    distanceFunctionOption->require_option(1, 3);

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

    slic::setLoggingMsgLevel(_verboseOutput ? slic::message::Debug : slic::message::Info);

    ndim = std::max({ndim, fcnCenter.size(), inPlane.size(), perpDir.size(), gyroidScale.size()});
    SLIC_ASSERT_MSG((fcnCenter.empty() || fcnCenter.size() == ndim) &&
                      (inPlane.empty() || inPlane.size() == ndim) &&
                      (perpDir.empty() || perpDir.size() == ndim) &&
                      (gyroidScale.empty() || gyroidScale.size() == ndim),
                    "--center, --inPlane, --dir, and --scale must have matching dimensions.");

    // The plane passes through the origin when --inPlane is omitted.
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

    _maxSpacing = maxSpacing();
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

  bool isUnstructured(axom::IndexType domId) const { return topologyType(domId) == "unstructured"; }

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

  const std::string& coordsetPath() const { return _coordsetPath; }

  //! @brief Corner node ids of a zone, for an unstructured single-shape topology.
  void unstructuredCellNodeIds(axom::IndexType domId,
                               axom::IndexType cellId,
                               axom::Array<axom::IndexType>& nodeIds) const
  {
    const conduit::Node& elems = domain(domId).fetch_existing(_topologyPath + "/elements");
    const std::string shape = elems.fetch_existing("shape").as_string();
    const axom::IndexType cornersPerCell = shape == "hex" ? 8 : 4;
    const auto conn = elems.fetch_existing("connectivity").as_index_t_accessor();
    nodeIds.resize(cornersPerCell);
    for(axom::IndexType c = 0; c < cornersPerCell; ++c)
    {
      nodeIds[c] = static_cast<axom::IndexType>(conn[cellId * cornersPerCell + c]);
    }
  }

  /*!
   * @return The maximum cell-edge length across all ranks.
   *
   * The constructor caches this value before the coordinate data may move to a device.
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
   * @return The maximum cell-edge length in one domain.
   *
   * The compact structured and unstructured paths inspect every edge.
   * The strided structured path samples one edge per axis and assumes Cartesian coordinates.
   */
  double maxSpacing1(axom::IndexType domId) const
  {
    const conduit::Node& dom = domain(domId);
    if(useFlatFields(domId) && isStructured(domId))
    {
      return maxStructuredEdgeLengthFlat(dom);
    }
    if(isUnstructured(domId))
    {
      return maxUnstructuredEdgeLength(dom);
    }

    const conduit::Node& dimsNode = dom.fetch_existing("topologies/mesh/elements/dims");
    axom::Array<axom::IndexType> ls(_ndims);
    for(int d = 0; d < _ndims; ++d)
    {
      ls[d] = 1 + static_cast<axom::IndexType>(dimsNode[d].to_int64());
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

  double maxStructuredEdgeLengthFlat(const conduit::Node& dom) const
  {
    const conduit::Node& dimsNode = dom.fetch_existing("topologies/mesh/elements/dims");
    axom::StackArray<axom::IndexType, 3> nodeShape {{1, 1, 1}};
    nodeShape[0] = dimsNode.fetch_existing("i").to_int64() + 1;
    nodeShape[1] = dimsNode.fetch_existing("j").to_int64() + 1;
    if(_ndims == 3)
    {
      nodeShape[2] = dimsNode.fetch_existing("k").to_int64() + 1;
    }

    const conduit::Node& coords = dom.fetch_existing(_coordsetPath + "/values");
    const auto xs = coords.fetch_existing("x").as_double_accessor();
    const auto ys = coords.fetch_existing("y").as_double_accessor();
    const bool hasZ = _ndims == 3;
    const auto zs = hasZ ? coords.fetch_existing("z").as_double_accessor()
                         : coords.fetch_existing("x").as_double_accessor();

    auto nodeIndex = [&](axom::IndexType i, axom::IndexType j, axom::IndexType k) {
      return i + j * nodeShape[0] + k * nodeShape[0] * nodeShape[1];
    };

    double maxLen = 0.0;
    for(axom::IndexType k = 0; k < nodeShape[2]; ++k)
    {
      for(axom::IndexType j = 0; j < nodeShape[1]; ++j)
      {
        for(axom::IndexType i = 0; i < nodeShape[0]; ++i)
        {
          const axom::IndexType a = nodeIndex(i, j, k);
          const axom::IndexType maxAxis = hasZ ? 3 : 2;
          for(axom::IndexType axis = 0; axis < maxAxis; ++axis)
          {
            axom::IndexType ni = i, nj = j, nk = k;
            if(axis == 0)
            {
              ++ni;
            }
            else if(axis == 1)
            {
              ++nj;
            }
            else
            {
              ++nk;
            }
            if(ni >= nodeShape[0] || nj >= nodeShape[1] || nk >= nodeShape[2])
            {
              continue;
            }
            const axom::IndexType b = nodeIndex(ni, nj, nk);
            const double dx = xs[a] - xs[b];
            const double dy = ys[a] - ys[b];
            const double dz = hasZ ? zs[a] - zs[b] : 0.0;
            maxLen = std::max(maxLen, std::sqrt(dx * dx + dy * dy + dz * dz));
          }
        }
      }
    }
    return maxLen;
  }

  /*!
   * @brief Longest cell edge over an unstructured single-shape topology.
   *
   * This uses Bump's shape traits for edge connectivity.
   */
  template <typename ShapeTraits>
  double maxEdgeLengthForTraits(const conduit::Node& topo, const conduit::Node& coords) const
  {
    const auto xs = coords.fetch_existing("x").as_double_accessor();
    const auto ys = coords.fetch_existing("y").as_double_accessor();
    const bool hasZ = _ndims == 3;
    const auto zs = hasZ ? coords.fetch_existing("z").as_double_accessor()
                         : coords.fetch_existing("x").as_double_accessor();
    const auto conn = topo.fetch_existing("elements/connectivity").as_index_t_accessor();

    constexpr auto cornersPerCell = ShapeTraits::numberOfNodes();
    constexpr auto edgeCount = ShapeTraits::numberOfEdges();
    const axom::IndexType numCells =
      static_cast<axom::IndexType>(conn.number_of_elements()) / cornersPerCell;

    double maxLen = 0.0;
    for(axom::IndexType cell = 0; cell < numCells; ++cell)
    {
      for(int e = 0; e < edgeCount; ++e)
      {
        const auto edge = ShapeTraits::getEdge(e);
        const auto a = static_cast<axom::IndexType>(conn[cell * cornersPerCell + edge[0]]);
        const auto b = static_cast<axom::IndexType>(conn[cell * cornersPerCell + edge[1]]);
        const double dx = xs[a] - xs[b];
        const double dy = ys[a] - ys[b];
        const double dz = hasZ ? zs[a] - zs[b] : 0.0;
        maxLen = std::max(maxLen, std::sqrt(dx * dx + dy * dy + dz * dz));
      }
    }
    return maxLen;
  }

  double maxUnstructuredEdgeLength(const conduit::Node& dom) const
  {
    const conduit::Node& topo = dom.fetch_existing(_topologyPath);
    const conduit::Node& coords = dom.fetch_existing(_coordsetPath + "/values");
    const std::string shape = topo.fetch_existing("elements/shape").as_string();

    if(shape == "hex")
    {
      return maxEdgeLengthForTraits<bumpviews::HexTraits>(topo, coords);
    }
    if(shape == "quad")
    {
      return maxEdgeLengthForTraits<bumpviews::QuadTraits>(topo, coords);
    }
    SLIC_ERROR(axom::fmt::format("Unsupported unstructured shape '{}'.", shape));
    return 0.0;
  }

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
  double _maxSpacing = -1.0;

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
 * @brief Analytic scalar field and its Blueprint field name.
 *
 * Command-line options select fields at run time, so \c valueAt uses \c std::function.
 * Only host loops call it.
 */
template <int DIM>
struct ContourTestStrategy
{
  using PointType = axom::primal::Point<double, DIM>;

  std::string testName;
  std::string functionName;
  std::function<double(const PointType&)> valueAt;
};

//! @brief Create a signed-distance field for a plane.
template <int DIM>
ContourTestStrategy<DIM> makePlanarStrategy(const axom::primal::Vector<double, DIM>& perpDir,
                                            const axom::primal::Point<double, DIM>& inPlane)
{
  const axom::primal::Plane<double, DIM> plane(perpDir.unitVector(), inPlane);
  return {"planar", "dist_to_plane", [plane](const axom::primal::Point<double, DIM>& pt) {
            return plane.signedDistance(pt);
          }};
}

//! @brief Create the distance-to-center field used for round contours.
template <int DIM>
ContourTestStrategy<DIM> makeRoundStrategy(const axom::primal::Point<double, DIM>& center)
{
  const axom::primal::Sphere<double, DIM> sphere(center, 0.0);
  return {"round", "dist_to_center", [sphere](const axom::primal::Point<double, DIM>& pt) {
            return sphere.computeSignedDistance(pt);
          }};
}

//! @brief Create a gyroid field shifted by @a offset.
template <int DIM>
ContourTestStrategy<DIM> makeGyroidStrategy(const axom::primal::Point<double, DIM>& scale,
                                            double offset)
{
  return {"gyroid", "gyroid_fcn", [scale, offset](const axom::primal::Point<double, DIM>& pt) {
            if(DIM == 3)
            {
              return sin(pt[0] * scale[0]) * cos(pt[1] * scale[1]) +
                sin(pt[1] * scale[1]) * cos(pt[2] * scale[2]) +
                sin(pt[2] * scale[2]) * cos(pt[0] * scale[0]) + offset;
            }
            // Evaluate the 3D function at z = 0.
            return sin(pt[0] * scale[0]) * cos(pt[1] * scale[1]) + sin(pt[1] * scale[1]) + offset;
          }};
}

template <int DIM, typename ExecSpace>
struct ContourTestBase
{
  static constexpr auto MemorySpace = axom::execution_space<ExecSpace>::memory_space;
  using PointType = axom::primal::Point<double, DIM>;
  explicit ContourTestBase(const Input& params)
    : m_params(params)
    , m_testStrategies()
    , m_parentCellIdField("parentCellIds")
    , m_domainIdField("domainIdField")
  { }
  void addTestStrategy(const ContourTestStrategy<DIM>& testStrategy)
  {
    m_testStrategies.push_back(testStrategy);
    SLIC_INFO(axom::fmt::format("Added contour field '{}'.", testStrategy.testName));
  }

  const Input& m_params;
  axom::Array<ContourTestStrategy<DIM>> m_testStrategies;
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
      for(const auto& strategy : m_testStrategies)
      {
        const std::string dataPath = axom::fmt::format("fields/{}/values", strategy.functionName);
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
        for(const auto& strategy : m_testStrategies)
        {
          mc.setFunctionField(strategy.functionName);
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

  void computeNodalDistance(BlueprintStructuredMesh& bpMesh, ContourTestStrategy<DIM>& strat)
  {
    AXOM_ANNOTATE_SCOPE("computeNodalDistance");

    SLIC_ASSERT(bpMesh.dimension() == DIM);
    for(int domId = 0; domId < bpMesh.domainCount(); ++domId)
    {
      if(bpMesh.useFlatFields(domId))
      {
        computeNodalDistanceFlat(bpMesh.domain(domId), strat);
        continue;
      }

      auto domainView = bpMesh.getDomainView<DIM>(domId);

      // Match the field's padded layout to the coordset so node ids remain valid.
      domainView.createField(strat.functionName,
                             "vertex",
                             conduit::DataType::float64(domainView.getCoordsCountWithGhosts()),
                             domainView.getCoordsStrides(),
                             domainView.getCoordsOffsets());
    }

    for(int domId = 0; domId < bpMesh.domainCount(); ++domId)
    {
      if(bpMesh.useFlatFields(domId))
      {
        continue;
      }

      auto domainView = bpMesh.getDomainView<DIM>(domId);
      const auto coordsViews = domainView.getConstCoordsViews(false);
      axom::ArrayView<double, DIM> fieldView =
        domainView.template getFieldView<double>(strat.functionName, false);
      for(int d = 0; d < DIM; ++d)
      {
        SLIC_ASSERT(coordsViews[d].shape() == fieldView.shape());
      }
      populateNodalDistance(coordsViews, fieldView, strat);
    }
  }

  void computeNodalDistanceFlat(conduit::Node& dom, ContourTestStrategy<DIM>& strat)
  {
    conduit::Node& fieldNode = dom["fields/" + strat.functionName];
    fieldNode["association"] = "vertex";
    fieldNode["topology"] = "mesh";

    const conduit::Node& values = dom.fetch_existing("coordsets/coords/values");
    const auto xs = values.fetch_existing("x").as_double_accessor();
    const auto ys = values.fetch_existing("y").as_double_accessor();
    const bool hasZ = DIM == 3;
    const auto zs = hasZ ? values.fetch_existing("z").as_double_accessor()
                         : values.fetch_existing("x").as_double_accessor();

    const conduit::index_t nodeCount = xs.number_of_elements();
    fieldNode["values"].set(conduit::DataType::float64(nodeCount));
    auto* fieldValues = fieldNode["values"].as_double_ptr();

    for(conduit::index_t nodeId = 0; nodeId < nodeCount; ++nodeId)
    {
      PointType pt;
      pt[0] = xs[nodeId];
      pt[1] = ys[nodeId];
      if(DIM == 3)
      {
        pt[2] = zs[nodeId];
      }
      fieldValues[nodeId] = strat.valueAt(pt);
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
      axom::for_all<axom::SEQ_EXEC>(0, cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
        maskView.flatIndex(cellId) = (cellId % maskCount);
      });
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

  void computeNodalDistance(BlueprintStructuredMesh& bpMesh)
  {
    for(auto& strategy : m_testStrategies)
    {
      computeNodalDistance(bpMesh, strategy);
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
TestInstanceVariant selectTestDimension(TypeTag<ExecSpace>, const Input& params)
{
  if(params.ndim == 2)
  {
    return TestInstance<2, ExecSpace> {};
  }
  if(params.ndim == 3)
  {
    return TestInstance<3, ExecSpace> {};
  }

  SLIC_ERROR(axom::fmt::format("Unsupported mesh dimension {}", params.ndim));
  return TestInstance<2, axom::SEQ_EXEC> {};
}

TestInstanceVariant selectTestInstance(const Input& params)
{
  if(params.policy == RuntimePolicy::seq)
  {
    return selectTestDimension(TypeTag<axom::SEQ_EXEC> {}, params);
  }
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  if(params.policy == RuntimePolicy::omp)
  {
    return selectTestDimension(TypeTag<axom::OMP_EXEC> {}, params);
  }
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(params.policy == RuntimePolicy::cuda)
  {
    return selectTestDimension(TypeTag<axom::CUDA_EXEC<256>> {}, params);
  }
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  if(params.policy == RuntimePolicy::hip)
  {
    return selectTestDimension(TypeTag<axom::HIP_EXEC<256>> {}, params);
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

  SLIC_ERROR_IF(
    params.ndim != static_cast<std::size_t>(computationalMesh.dimension()),
    axom::fmt::format(
      "Function parameter dimension {} does not match input mesh dimension {} for '{}'.",
      params.ndim,
      computationalMesh.dimension(),
      params.meshFile));

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
  auto testInstance = selectTestInstance(params);
  int errCount = std::visit(
    [&](const auto& instance) {
      AXOM_UNUSED_VAR(instance);
      using Instance = std::decay_t<decltype(instance)>;
      constexpr int DIM = Instance::DIM;
      using ExecSpace = typename Instance::ExecSpace;

      ContourTestBase<DIM, ExecSpace> contourTest(params);

      if(params.usingPlanar())
      {
        contourTest.addTestStrategy(
          makePlanarStrategy<DIM>(params.planeNormal<DIM>(), params.inplanePoint<DIM>()));
      }

      if(params.usingRound())
      {
        contourTest.addTestStrategy(makeRoundStrategy<DIM>(params.roundContourCenter<DIM>()));
      }

      if(params.usingGyroid())
      {
        contourTest.addTestStrategy(
          makeGyroidStrategy<DIM>(params.gyroidScaleFactor<DIM>(), params.contourVal));
      }

      contourTest.computeNodalDistance(computationalMesh);
      contourTest.addMaskField(computationalMesh);

      if(params.isVerbose())
      {
        computationalMesh.printMeshInfo();
      }

      // Write computational mesh with contour functions.
      saveMesh(computationalMesh.asConduitNode(), params.fieldsFile);

      int localErrCount = contourTest.runTest(computationalMesh);

      int globalErrCount = 0;
#ifdef AXOM_USE_MPI
      MPI_Allreduce(&localErrCount, &globalErrCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
      globalErrCount = localErrCount;
#endif

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
