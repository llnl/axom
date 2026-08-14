// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_distributed_distance_query_example.cpp
 * \brief Driver for a distributed distance query
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/quest.hpp"
#include "axom/slam.hpp"
#include "axom/core/Types.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_io_blueprint.hpp"
#include "conduit_relay_mpi_io_blueprint.hpp"

#include "axom/quest/DistributedClosestPoint.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifndef AXOM_USE_MPI
  #error This example requires Axom to be configured with MPI
#endif
#include "mpi.h"

// C/C++ includes
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <atomic>
#include <chrono>
#include <thread>
#if defined(__GLIBC__)
  #include <malloc.h>  // mallinfo2 / mallinfo / malloc_trim
#endif

namespace quest = axom::quest;
namespace slic = axom::slic;
namespace sidre = axom::sidre;
namespace slam = axom::slam;
namespace primal = axom::primal;

using RuntimePolicy = axom::runtime_policy::Policy;

// MPI stuff, initialized in main()
int my_rank = -1, num_ranks = -1;

// converts the input string into an 80 character string
// padded on both sides with '=' symbols
std::string banner(const std::string& str) { return axom::fmt::format("{:=^80}", str); }

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;
  std::string distanceFile {"cp_coords"};
  std::string objectMeshFile;

  // Memory instrumentation (all off by default).
  bool trackMemory {false};     // report RSS / malloc / umpire at phase boundaries
  bool trimAfterQuery {false};  // malloc_trim(0) after the query, then re-report
  int sampleMemoryMs {0};       // if > 0, background-sample peak RSS during the query

  RuntimePolicy policy {RuntimePolicy::seq};

  double distThreshold {axom::numeric_limits<double>::max()};

  bool dynamicDistanceFiltering {true};

private:
  bool m_verboseOutput {false};

public:
  bool isVerbose() const { return m_verboseOutput; }

  std::string getDCMeshName() const
  {
    using axom::utilities::string::removeSuffix;

    // Remove the parent directories and file suffix
    std::string name = axom::Path(meshFile).baseName();
    name = removeSuffix(name, ".root");

    return name;
  }

  std::string getMdMeshName() const
  {
    using axom::utilities::string::removeSuffix;

    // Remove the parent directories and file suffix
    std::string name = axom::Path(meshFile).baseName();
    name = removeSuffix(name, ".root");

    return name;
  }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-m,--mesh-file", meshFile)
      ->description(
        "Path to multidomain computational mesh following conduit blueprint "
        "convention.")
      ->check(axom::CLI::ExistingFile);

    app.add_option("-s,--distance-file", distanceFile)
      ->description("Name of output mesh file containing closest distance.")
      ->capture_default_str();

    app.add_option("--object-mesh-file", objectMeshFile)
      ->description(
        "Path to a conduit blueprint point mesh root file. "
        "Generate this mesh with src/tools/gen-multidom-point-mesh.py.")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_flag("--track-memory", trackMemory)
      ->description(
        "Report RSS, glibc malloc live/arena bytes (and Umpire high-water when available) "
        "before/after BVH build and before/after the closest-point query, reduced across ranks.")
      ->capture_default_str();

    app.add_flag("--trim-after-query", trimAfterQuery)
      ->description(
        "After the query, call malloc_trim(0) and report memory again. "
        "If RSS drops here (while malloc-live was already back at baseline), "
        "the memory was arena retention, not a leak.  Implies --track-memory.")
      ->capture_default_str();

    app.add_option("--sample-memory-ms", sampleMemoryMs)
      ->description(
        "If > 0, run a background thread sampling RSS at this interval (ms) "
        "during the query and report the peak.  Useful for observing in-flight "
        "send-buffer accumulation.  Implies --track-memory.")
      ->capture_default_str();

    app.add_flag("-v,--verbose,!--no-verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("-d,--dist-threshold", distThreshold)
      ->check(axom::CLI::NonNegativeNumber)
      ->description("Distance threshold to search")
      ->capture_default_str();

    app
      .add_flag("--dynamic-distance-filtering,!--no-dynamic-distance-filtering",
                dynamicDistanceFiltering)
      ->description("Enable/disable dynamic distance filtering")
      ->capture_default_str();

    app.add_option("-p, --policy", policy)
      ->description("Set runtime policy for point query method")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug : slic::message::Info);
  }
};

// Input params set in main()
Input params;

/**
 *  \brief Simple wrapper to a blueprint particle mesh
 *
 *  Given a sidre Group, creates the stubs for a mesh blueptint particle mesh
 *
 *  BlueprintParticleMesh is used by both the object mesh and the query mesh.
 */
struct BlueprintParticleMesh
{
public:
  explicit BlueprintParticleMesh(sidre::Group* group,
                                 const std::string& topology,
                                 const std::string& coordset)
    : m_topologyName(topology)
    , m_coordsetName(coordset)
    , m_group(group)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nranks);
  }

  explicit BlueprintParticleMesh(sidre::Group* group)
    : m_topologyName()
    , m_coordsetName()
    , m_group(group)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nranks);
  }

  /// Gets the root group for this mesh blueprint
  sidre::Group* root_group() const { return m_group; }

  /// Gets number of domains in the multidomain particle mesh
  axom::IndexType domain_count() const { return m_group->getNumGroups(); }

  /// Gets a domain group.
  sidre::Group* domain_group(axom::IndexType groupIdx) const
  {
    SLIC_ASSERT(groupIdx < m_group->getNumGroups());
    return m_group->getGroup(groupIdx);
  }
  /// Gets the parent group for the blueprint coordinate set
  sidre::Group* coords_group(axom::IndexType groupIdx) const
  {
    return domain_group(groupIdx)->getGroup("coordsets")->getGroup(m_coordsetName);
  }
  /// Gets the parent group for the blueprint mesh topology
  sidre::Group* topo_group(axom::IndexType groupIdx) const
  {
    return domain_group(groupIdx)->getGroup("topologies")->getGroup(m_topologyName);
  }
  /// Gets the parent group for the blueprint fields
  sidre::Group* fields_group(axom::IndexType groupIdx) const
  {
    auto* domain = domain_group(groupIdx);
    return domain->hasGroup("fields") ? domain->getGroup("fields") : domain->createGroup("fields");
  }

  const std::string& getTopologyName() const { return m_topologyName; }
  const std::string& getCoordsetName() const { return m_coordsetName; }

  /// Gets the MPI rank for this mesh
  int getRank() const { return m_rank; }
  /// Gets the number of ranks in the problem
  int getNumRanks() const { return m_nranks; }

  /*!
    @brief Returns the number of points in a particle mesh domain
    including ghost points.
  */
  int numPoints(axom::IndexType dIdx) const
  {
    int rval = 0;
    auto* cg = coords_group(dIdx);
    SLIC_ASSERT(cg != nullptr && cg->hasView("values/x"));
    rval = cg->getView("values/x")->getNumElements();
    return rval;
  }
  /// Returns the number of points in the particle mesh
  int numPoints() const
  {
    int rval = 0;
    const axom::IndexType domCount = domain_count();
    for(axom::IndexType dIdx = 0; dIdx < domCount; ++dIdx)
    {
      rval += numPoints(dIdx);
    }
    return rval;
  }

  int dimension() const { return m_dimension; }

  /*!
    @brief Read a blueprint mesh and store it internally in m_group.

    If the topology wasn't specified in the constructor, the first
    topology from the file is used.  The coordset name will be
    replaced with the one corresponding to the topology.
  */
  void read_blueprint_mesh(const std::string& meshFilename)
  {
    SLIC_ASSERT(!meshFilename.empty());

    conduit::Node mdMesh;
    conduit::relay::mpi::io::blueprint::load_mesh(meshFilename, mdMesh, MPI_COMM_WORLD);
    if(!conduit::blueprint::mesh::is_multi_domain(mdMesh) && mdMesh.number_of_children() > 0)
    {
      conduit::Node singleDomainMesh;
      singleDomainMesh.update(mdMesh);
      mdMesh.reset();
      conduit::blueprint::mesh::to_multi_domain(singleDomainMesh, mdMesh);
    }
    conduit::index_t domCount = conduit::blueprint::mesh::number_of_domains(mdMesh);

    if(domCount > 0)
    {
      if(m_topologyName.empty())
      {
        // No topology given.  Pick the first one.
        m_topologyName = mdMesh[0].fetch_existing("topologies")[0].name();
      }
      auto topologyPath = axom::fmt::format("topologies/{}", m_topologyName);

      // Detect strided structured coordinates.
      // Structured topologies have elements/dims, but points and unstructured topology don't.
      // Guard the probe use the resolved topology name rather than a hardcoded "mesh" topology.
      const std::string dimsPath = topologyPath + "/elements/dims";
      m_coordsAreStrided =
        mdMesh[0].has_path(dimsPath) && mdMesh[0].fetch_existing(dimsPath).has_child("strides");
      if(m_coordsAreStrided)
      {
        SLIC_WARNING(
          axom::fmt::format("Mesh '{}' is strided.  Stride support is under development.",
                            meshFilename));
      }

      m_coordsetName = mdMesh[0].fetch_existing(topologyPath + "/coordset").as_string();
      const conduit::Node coordsetNode =
        mdMesh[0].fetch_existing("coordsets").fetch_existing(m_coordsetName);
      m_dimension = conduit::blueprint::mesh::coordset::dims(coordsetNode);
    }

    MPI_Allreduce(MPI_IN_PLACE, &m_dimension, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    SLIC_ASSERT(m_dimension > 0);

    if(domCount > 0)
    {
      // Put mdMesh into sidre Group.
      const bool goodImport = m_group->importConduitTree(mdMesh, false);
      SLIC_ASSERT(goodImport);
      SLIC_ASSERT(m_group->getNumGroups() == domCount);
      AXOM_UNUSED_VAR(goodImport);
    }

    bool valid = isValid();
    SLIC_ASSERT(valid);
    AXOM_UNUSED_VAR(valid);
  }

  /*!  @brief Set the coordinate data from an array of primal Points

    The points are assigned to a new domain (in the multidomain context).

    This method is for manually creating the mesh.  Don't use it if
    the mesh is read in.
  */
  void setPoints(axom::IndexType domainId, const axom::Array<double, 2>& pts)
  {
    axom::IndexType localIdx = createBlueprintStubs();
    SLIC_ASSERT(domain_group(localIdx) != nullptr);
    domain_group(localIdx)->createViewScalar<std::int64_t>("state/domain_id", domainId);

    const int SZ = pts.shape()[0];

    if(m_dimension == -1)
    {
      m_dimension = pts.shape()[1];
    }
    else
    {
      SLIC_ASSERT(pts.shape()[1] == m_dimension);
    }

    // lambda to create a strided view into the buffer
    // uses workaround for empty meshes since apply() requires size > 0
    const auto dimension = m_dimension;
    auto createAndApplyView =
      [=](sidre::Group* grp, const std::string& path, sidre::Buffer* buf, int dim, int sz) {
        if(sz > 0)
        {
          grp->createView(path)->attachBuffer(buf)->apply(sz, dim, dimension);
        }
        else
        {
          grp->createViewAndAllocate(path, sidre::DOUBLE_ID, 0);
        }
      };

    // create views into a shared buffer for the coordinates, with stride m_dimension
    {
      auto* buf = domain_group(localIdx)
                    ->getDataStore()
                    ->createBuffer(sidre::DOUBLE_ID, m_dimension * SZ)
                    ->allocate();

      createAndApplyView(coords_group(localIdx), "values/x", buf, 0, SZ);
      if(m_dimension > 1)
      {
        createAndApplyView(coords_group(localIdx), "values/y", buf, 1, SZ);
      }
      if(m_dimension > 2)
      {
        createAndApplyView(coords_group(localIdx), "values/z", buf, 2, SZ);
      }

      // copy coordinate data into the buffer
      const std::size_t nbytes = sizeof(double) * SZ * m_dimension;
      axom::copy(buf->getVoidPtr(), pts.data(), nbytes);
    }

    // set the default connectivity
    // May be required by an old version of visit.  May not be needed by newer versions of visit.
    sidre::Array<int> arr(topo_group(localIdx)->createView("elements/connectivity"), SZ, SZ);
    for(int i = 0; i < SZ; ++i)
    {
      arr[i] = i;
    }
  }

  template <int DIM>
  axom::Array<primal::Point<double, DIM>> getPoints(int domainIdx)
  {
    auto* cGroup = coords_group(domainIdx);
    auto* xView = cGroup->getView("values/x");
    auto* yView = cGroup->getView("values/y");
    auto* zView = DIM >= 3 ? cGroup->getView("values/z") : nullptr;
    const auto ptCount = xView->getNumElements();
    assert(xView->getStride() == 1);
    assert(yView->getStride() == 1);
    assert(zView == nullptr || zView->getStride() == 1);
    double* xs = xView->getArray();
    double* ys = yView->getArray();
    double* zs = zView ? (double*)(zView->getArray()) : nullptr;

    using PointType = primal::Point<double, DIM>;
    axom::Array<PointType> pts;
    pts.resize(ptCount);
    for(int i = 0; i < ptCount; ++i)
    {
      pts[i][0] = xs[i];
    }
    for(int i = 0; i < ptCount; ++i)
    {
      pts[i][1] = ys[i];
    }
    if(DIM == 3)
    {
      for(int i = 0; i < ptCount; ++i)
      {
        pts[i][2] = zs[i];
      }
    }
    return pts;
  }

  void printMeshSizeStats(const std::string& meshLabel) const
  {
    SLIC_INFO(axom::fmt::format("{} has {} points in {} domains locally",
                                meshLabel,
                                numPoints(),
                                domain_count()));

    auto getIntMinMax = [](int inVal, int& minVal, int& maxVal, int& sumVal) {
      MPI_Allreduce(&inVal, &minVal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &maxVal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &sumVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    };

    // Output some global mesh size stats
    {
      int mn, mx, sum;
      getIntMinMax(numPoints(), mn, mx, sum);
      SLIC_INFO(axom::fmt::format("{} has {{min:{}, max:{}, sum:{}, avg:{}}} points",
                                  meshLabel,
                                  mn,
                                  mx,
                                  sum,
                                  (double)sum / num_ranks));
    }
    {
      int mn, mx, sum;
      getIntMinMax(domain_count(), mn, mx, sum);
      SLIC_INFO(axom::fmt::format("{} has {{min:{}, max:{}, sum:{}, avg:{}}} domains",
                                  meshLabel,
                                  mn,
                                  mx,
                                  sum,
                                  (double)sum / num_ranks));
    }
  }

  template <typename T>
  void registerNodalScalarField(const std::string& fieldName)
  {
    for(axom::IndexType dIdx = 0; dIdx < domain_count(); ++dIdx)
    {
      auto* fld = fields_group(dIdx)->createGroup(fieldName);
      fld->createViewString("association", "vertex");
      fld->createViewString("topology", topo_group(dIdx)->getName());
      if(m_coordsAreStrided)
      {
        auto* offsets = topo_group(dIdx)->getView("elements/dims/offsets");
        auto* strides = topo_group(dIdx)->getView("elements/dims/strides");
        fld->copyView(offsets);
        fld->copyView(strides);
      }
      fld->createViewAndAllocate("values", sidre::detail::SidreTT<T>::id, numPoints(dIdx));
    }
  }

  template <typename T>
  void registerNodalVectorField(const std::string& fieldName)
  {
    const int DIM = dimension();
    for(axom::IndexType dIdx = 0; dIdx < domain_count(); ++dIdx)
    {
      const int SZ = numPoints(dIdx);

      auto* fld = fields_group(dIdx)->createGroup(fieldName);
      fld->createViewString("association", "vertex");
      fld->createViewString("topology", topo_group(dIdx)->getName());
      if(m_coordsAreStrided)
      {
        auto* offsets = topo_group(dIdx)->getView("elements/dims/offsets");
        auto* strides = topo_group(dIdx)->getView("elements/dims/strides");
        fld->copyView(offsets);
        fld->copyView(strides);
      }

      if(SZ == 0)
      {
        fld->createViewAndAllocate("values/x", sidre::detail::SidreTT<T>::id, 0);
        if(DIM > 1)
        {
          fld->createViewAndAllocate("values/y", sidre::detail::SidreTT<T>::id, 0);
        }
        if(DIM > 2)
        {
          fld->createViewAndAllocate("values/z", sidre::detail::SidreTT<T>::id, 0);
        }
        continue;
      }

      // create views into a shared buffer for the coordinates, with stride DIM
      auto* buf = domain_group(dIdx)
                    ->getDataStore()
                    ->createBuffer(sidre::detail::SidreTT<T>::id, DIM * SZ)
                    ->allocate();
      switch(DIM)
      {
      case 3:
        fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
        fld->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, DIM);
        fld->createView("values/z")->attachBuffer(buf)->apply(SZ, 2, DIM);
        break;
      case 2:
        fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
        fld->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, DIM);
        break;
      default:
        fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
        break;
      }
    }
  }

  bool hasScalarField(const std::string& fieldName, int domainIdx = 0) const
  {
    auto* domain = m_group->getGroup(domainIdx);
    auto* fields = domain->getGroup("fields");
    return fields != nullptr && fields->hasGroup(fieldName);
  }

  bool hasVectorField(const std::string& fieldName, int domainIdx = 0) const
  {
    auto* fields = m_group->getGroup(domainIdx)->getGroup("fields");
    return fields != nullptr && fields->hasGroup(fieldName);
  }

  template <typename T>
  axom::ArrayView<T> getNodalScalarField(const std::string& fieldName, int domainIdx)
  {
    SLIC_ASSERT_MSG(domainIdx >= 0 && axom::IndexType(domainIdx) < domain_count(),
                    axom::fmt::format("Rank {} has no domain {}, only {} domains",
                                      m_rank,
                                      domainIdx,
                                      domain_count()));

    auto* domain = m_group->getGroup(domainIdx);
    auto* fields = domain->getGroup("fields");
    auto* field = fields != nullptr ? fields->getGroup(fieldName) : nullptr;
    T* data = field ? static_cast<T*>(field->getView("values")->getVoidPtr()) : nullptr;

    return field ? axom::ArrayView<T>(data, numPoints(domainIdx)) : axom::ArrayView<T>();
  }

  template <typename T>
  axom::ArrayView<T> getNodalVectorField(const std::string& fieldName, int domainIdx)
  {
    SLIC_ASSERT_MSG(domainIdx >= 0 && axom::IndexType(domainIdx) < domain_count(),
                    axom::fmt::format("Rank {} has only {} domains, no domain index {}",
                                      m_rank,
                                      domain_count(),
                                      domainIdx));

    // Note: the implementation currently assumes that the field data is
    // interleaved, so it is safe to get a pointer to the beginning of the
    // x-coordinate's data. This will be relaxed in the future, and we will
    // need to modify this implementation accordingly.
    T* data = nullptr;
    axom::IndexType npts = 0;
    auto* fields = m_group->getGroup(domainIdx)->getGroup("fields");
    bool has = fields != nullptr && fields->hasGroup(fieldName);
    if(has)
    {
      auto xView = fields->getGroup(fieldName)->getView("values/x");
      data = static_cast<T*>(xView->getVoidPtr());
      npts = xView->getNumElements();
    }
    return axom::ArrayView<T>(data, npts);
  }

  /// Returns an array containing the positions of the mesh vertices
  axom::Array<double, 2> getVertexPositions(axom::IndexType domainIdx)
  {
    sidre::Group* cvg =
      getDomain(domainIdx)->getGroup(axom::fmt::format("coordsets/{}/values", getCoordsetName()));
    int ndim = cvg->getNumViews();
    sidre::View* xv = cvg->getView("x");
    sidre::View* yv = cvg->getView("y");
    sidre::View* zv = ndim == 3 ? cvg->getView("z") : nullptr;
    axom::IndexType npts = xv->getNumElements();
    double* xp = xv->getData();
    double* yp = yv->getData();
    double* zp = zv ? (double*)(zv->getData()) : nullptr;
    double* xyzs[3] {xp, yp, zp};
    axom::Array<double, 2> rval(npts, ndim);
    for(int i = 0; i < npts; ++i)
    {
      for(int d = 0; d < ndim; ++d)
      {
        rval[i][d] = xyzs[d][i];
      }
    }
    return rval;
  }

  sidre::Group* getDomain(axom::IndexType domain) { return m_group->getGroup(domain); }
  sidre::Group* getFields(axom::IndexType domainIdx)
  {
    auto* domain = m_group->getGroup(domainIdx);
    return domain->hasGroup("fields") ? domain->getGroup("fields") : domain->createGroup("fields");
  }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    {
      conduit::Node meshNode;
      m_group->createNativeLayout(meshNode);
      conduit::Node info;
      if(!conduit::blueprint::mpi::verify("mesh", meshNode, info, MPI_COMM_WORLD))
      {
        SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
        slic::flushStreams();
        return false;
      }
    }
    return true;
  }

  /// Outputs the particle mesh to disk
  void saveMesh(const std::string& filename)
  {
    conduit::Node meshNode;
    m_group->createNativeLayout(meshNode);
    conduit::relay::mpi::io::blueprint::save_mesh(meshNode, filename, "hdf5", MPI_COMM_WORLD);
  }

  void print_mesh_info() const
  {
    // Copy to conduit::Node.  It's output is easier to read, especially in parallel.
    conduit::Node meshNode;
    m_group->createNativeLayout(meshNode);
    meshNode.print();
  }

private:
  /// Creates blueprint stubs for this mesh
  // for the "coordset", "topologies", "fields" and "state"
  // Return the domain index created.
  axom::IndexType createBlueprintStubs()
  {
    SLIC_ASSERT(m_group != nullptr);

    auto* domainGroup = m_group->createUnnamedGroup();

    auto* coordsGroup = domainGroup->createGroup("coordsets")->createGroup(m_coordsetName);
    coordsGroup->createViewString("type", "explicit");
    coordsGroup->createGroup("values");

    auto* topoGroup = domainGroup->createGroup("topologies")->createGroup(m_topologyName);
    topoGroup->createViewString("coordset", m_coordsetName);
    topoGroup->createViewString("type", "unstructured");
    topoGroup->createViewString("elements/shape", "point");

    domainGroup->createGroup("fields");
    domainGroup->createGroup("state");

    return m_group->getNumGroups() - 1;
  }

private:
  //!@brief Whether stride/offsets are given for blueprint mesh coordinates data.
  bool m_coordsAreStrided = false;
  std::string m_topologyName;
  std::string m_coordsetName;
  /// Parent group for the entire mesh
  sidre::Group* m_group;

  int m_rank;
  int m_nranks;
  int m_dimension {-1};
};  // BlueprintParticleMesh

/**
 * Helper class to generate a mesh blueprint-conforming particle mesh for the input object.
 * The mesh is represented using a Sidre hierarchy
 */
class ObjectMeshWrapper
{
public:
  //!@brief Construct by reading a blueprint object mesh from disk.
  //! Uses the topology/coordset names found in the file, exactly like the
  //! query mesh reader.  (The empty-topology BlueprintParticleMesh ctor lets
  //! read_blueprint_mesh pick the file's actual topology.)
  ObjectMeshWrapper(sidre::Group* group, const std::string& meshFilename) : m_objectMesh(group)
  {
    SLIC_ASSERT(group != nullptr);
    m_objectMesh.read_blueprint_mesh(meshFilename);
  }

  BlueprintParticleMesh& getParticleMesh() { return m_objectMesh; }

  /// Get a pointer to the root group for this mesh
  sidre::Group* getBlueprintGroup() const { return m_objectMesh.root_group(); }

  std::string getTopologyName() const { return m_objectMesh.getTopologyName(); }
  std::string getCoordsetName() const { return m_objectMesh.getCoordsetName(); }

private:
  BlueprintParticleMesh m_objectMesh;
};

class QueryMeshWrapper
{
public:
  //!@brief Construct with blueprint mesh.
  QueryMeshWrapper(sidre::Group* group, const std::string& meshFilename) : m_queryMesh(group)
  {
    // Test reading in multidomain mesh.
    m_queryMesh.read_blueprint_mesh(meshFilename);
    // setupParticleMesh();
  }

  BlueprintParticleMesh& getParticleMesh() { return m_queryMesh; }

  sidre::Group* getBlueprintGroup() const { return m_queryMesh.root_group(); }

  std::string getTopologyName() const { return m_queryMesh.getTopologyName(); }
  std::string getCoordsetName() const { return m_queryMesh.getCoordsetName(); }

  /// Saves the mesh to disk
  void saveMesh(const std::string& filename)
  {
    SLIC_INFO(banner(axom::fmt::format("Saving query mesh '{}' to disk", filename)));

    m_queryMesh.saveMesh(filename);
  }

  void setupParticleMesh()
  {
    {
      m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_rank");
      m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_index");
      m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_domain_index");
      m_queryMesh.registerNodalScalarField<double>("cp_distance");
      m_queryMesh.registerNodalVectorField<double>("cp_coords");
    }

    SLIC_ASSERT(m_queryMesh.isValid());
  }

  /// Prints some info about the mesh
  void print_mesh_info() { m_queryMesh.print_mesh_info(); }

  /*!
    @brief Update results from closest point search.
  */
  void update_closest_points(const conduit::Node& node)
  {
    sidre::Group* dstDomains = m_queryMesh.root_group();
    bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(node);
    if(!isMultidomain)
    {
      SLIC_ASSERT(!isMultidomain || dstDomains->getNumGroups() == node.number_of_children());
    }
    const int domainCount = dstDomains->getNumGroups();
    for(int d = 0; d < domainCount; ++d)
    {
      sidre::Group& domGroup = *dstDomains->getGroup(d);
      const conduit::Node& domNode = isMultidomain ? node.child(d) : node;

      sidre::Group& dstFieldsGroup =
        *(domGroup.hasGroup("fields") ? domGroup.getGroup("fields") : domGroup.createGroup("fields"));
      const conduit::Node& srcFieldsNode = domNode.fetch_existing("fields");
      {
        if(!m_queryMesh.hasScalarField("cp_rank"))
        {
          m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_rank");
        }
        auto dst = dstFieldsGroup.getGroup("cp_rank");
        auto src = srcFieldsNode.fetch_existing("cp_rank");
        bool goodImport = dst->importConduitTree(src);
        SLIC_ASSERT(goodImport);
        AXOM_UNUSED_VAR(goodImport);
      }
      {
        if(!m_queryMesh.hasScalarField("cp_index"))
        {
          m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_index");
        }
        auto dst = dstFieldsGroup.getGroup("cp_index");
        auto src = srcFieldsNode.fetch_existing("cp_index");
        bool goodImport = dst->importConduitTree(src);
        SLIC_ASSERT(goodImport);
        AXOM_UNUSED_VAR(goodImport);
      }
      if(srcFieldsNode.has_child("cp_domain_index"))
      {
        if(!m_queryMesh.hasScalarField("cp_domain_index"))
        {
          m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_domain_index");
        }
        auto src = srcFieldsNode.fetch_existing("cp_domain_index");
        auto dst = dstFieldsGroup.getGroup("cp_domain_index");
        bool goodImport = dst->importConduitTree(src);
        SLIC_ASSERT(goodImport);
        AXOM_UNUSED_VAR(goodImport);
      }
      {
        if(!m_queryMesh.hasVectorField("cp_coords"))
        {
          m_queryMesh.registerNodalVectorField<double>("cp_coords");
        }
        auto dstGroup = dstFieldsGroup.getGroup("cp_coords");
        auto srcNode = srcFieldsNode.fetch_existing("cp_coords");
        int dim = srcNode.fetch_existing("values").number_of_children();
        for(int d = 0; d < dim; ++d)
        {
          auto* dstView = dstGroup->getGroup("values")->getView(d);
          const auto& srcComponent = srcNode.fetch_existing("values").child(d);
          int nPts = srcComponent.dtype().number_of_elements();
          SLIC_ASSERT(nPts == dstView->getNumElements());
          if(nPts == 0)
          {
            continue;
          }

          conduit::float64_array dst = dstView->getArray();
          const conduit::float64_array src = srcComponent.value();
          for(int i = 0; i < nPts; ++i)
          {
            dst[i] = src[i];
          }
        }
      }
    }
  }

private:
  BlueprintParticleMesh m_queryMesh;
};

//---------------------------------------------------------------------------
// Transform closest points to distances and directions
//---------------------------------------------------------------------------
template <int DIM>
void computeDistancesAndDirections(BlueprintParticleMesh& queryMesh,
                                   const std::string& cpCoordsField,
                                   const std::string& cpIndexField,
                                   const std::string& distanceField,
                                   const std::string& directionField)
{
  SLIC_ASSERT(queryMesh.dimension() == DIM);

  using primal::squared_distance;
  using PointType = primal::Point<double, DIM>;
  using IndexSet = slam::PositionSet<>;

  PointType nowhere(axom::numeric_limits<double>::signaling_NaN());
  const double nodist = axom::numeric_limits<double>::signaling_NaN();

  queryMesh.registerNodalScalarField<double>(distanceField);
  queryMesh.registerNodalVectorField<double>(directionField);
  for(axom::IndexType di = 0; di < queryMesh.domain_count(); ++di)
  {
    auto cpCoords = queryMesh.getNodalVectorField<PointType>(cpCoordsField, di);

    auto cpIndices = queryMesh.getNodalScalarField<axom::IndexType>(cpIndexField, di);

    axom::Array<double, 2> qPts = queryMesh.getVertexPositions(di);
    axom::ArrayView<double> distances = queryMesh.getNodalScalarField<double>("distance", di);
    axom::ArrayView<PointType> directions = queryMesh.getNodalVectorField<PointType>("direction", di);
    axom::IndexType ptCount = queryMesh.numPoints(di);
    for(auto ptIdx : IndexSet(ptCount))
    {
      const bool has_cp = cpIndices[ptIdx] >= 0;
      const PointType& cp = has_cp ? cpCoords[ptIdx] : nowhere;
      const PointType& qPt = has_cp ? PointType(&qPts[ptIdx][0]) : nowhere;
      distances[ptIdx] = has_cp ? sqrt(squared_distance(qPt, cp)) : nodist;
      directions[ptIdx] = PointType(has_cp ? (cp - qPt).array() : nowhere.array());
    }
  }
}

void make_coords_contiguous(conduit::Node& coordValues)
{
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  if(isInterleaved)
  {
    conduit::Node oldValues = coordValues;
    conduit::blueprint::mcarray::to_contiguous(oldValues, coordValues);
  }
}

void make_coords_interleaved(conduit::Node& coordValues)
{
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  if(!isInterleaved)
  {
    conduit::Node oldValues = coordValues;
    conduit::blueprint::mcarray::to_interleaved(oldValues, coordValues);
  }
}

//-----------------------------------------------------------------------------
// Optional memory instrumentation
//
// Purpose: distinguish a genuine leak (memory still referenced after the call)
// from arena retention (memory freed at the allocator but not returned to the OS).
// RSS shows what the process actually holds; glibc "malloc live" shows what is still allocated;
// "malloc arena" shows free memory kept in the arena.
//   * RSS after call ~= baseline                -> transient, no problem
//   * RSS high, malloc-live ~= baseline         -> arena retention (malloc_trim should reclaim it)
//   * malloc-live high after call               -> genuine leak
// Everything here is Linux/glibc-specific and degrades to "n/a" elsewhere.
// All of it is opt-in (--track-memory) and off by default.
//-----------------------------------------------------------------------------

/// Read current (VmRSS) and peak (VmHWM) resident set size in bytes; -1 if n/a.
inline void readProcRss(long long& rssBytes, long long& peakRssBytes)
{
  rssBytes = -1;
  peakRssBytes = -1;
#if defined(__linux__)
  std::ifstream status("/proc/self/status");
  std::string line;
  while(std::getline(status, line))
  {
    long long kb = 0;
    if(std::sscanf(line.c_str(), "VmRSS: %lld kB", &kb) == 1)
    {
      rssBytes = kb * 1024;
    }
    else if(std::sscanf(line.c_str(), "VmHWM: %lld kB", &kb) == 1)
    {
      peakRssBytes = kb * 1024;
    }
  }
#endif
}

/// Reset the kernel's peak-RSS high-water mark (VmHWM) to the current RSS, so a
/// later VmHWM read reflects the peak of just the intervening phase.
/// Best-effort: requires Linux clear_refs type 5 (kernel >= 4.0).
inline void resetPeakRss()
{
#if defined(__linux__)
  std::ofstream clear("/proc/self/clear_refs");
  if(clear)
  {
    clear << "5\n";
  }
#endif
}

/// Human-readable byte count.
inline std::string humanBytes(long long b)
{
  if(b < 0)
  {
    return "n/a";
  }
  const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
  double v = static_cast<double>(b);
  int i = 0;
  while(v >= 1024.0 && i < 4)
  {
    v /= 1024.0;
    ++i;
  }
  return axom::fmt::format("{:.2f} {}", v, units[i]);
}

/// MPI-reduce a per-rank byte count to its max (with the rank achieving it) and
/// its sum across ranks.  A negative local value means "unavailable".
inline void reduceBytes(long long local, long long& maxVal, int& maxRank, long long& sumVal)
{
  struct
  {
    long val;
    int rank;
  } in {static_cast<long>(local), my_rank}, out {0, 0};
  MPI_Allreduce(&in, &out, 1, MPI_LONG_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  maxVal = out.val;
  maxRank = out.rank;
  MPI_Allreduce(&local, &sumVal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
}

/*!
 * \brief Opt-in per-run memory probe for the closest-point query.
 *
 * Reports RSS, peak RSS, glibc live/arena bytes (and Umpire high-water when built with Umpire),
 * reduced across ranks (total and hottest rank).  Also offers a background sampler to capture
 * the peak RSS *during* a phase (useful for observing in-flight send-buffer accumulation),
 * and a resetPeak() to make VmHWM phase-local.
 */
class MemoryProbe
{
public:
  MemoryProbe(bool enabled, int sampleMs, int umpireAllocatorId = -1)
    : m_enabled(enabled)
    , m_sampleMs(sampleMs)
    , m_umpireAllocatorId(umpireAllocatorId)
  { }

  bool enabled() const { return m_enabled; }

  /// Reset VmHWM so the next report()'s peak reflects only the next phase.
  void resetPeak()
  {
    if(m_enabled)
    {
      resetPeakRss();
    }
  }

  /// Start a background thread sampling RSS every sampleMs; no-op if disabled or sampleMs <= 0
  /// Records the peak RSS seen until stopSampler().
  void startSampler()
  {
    if(!m_enabled || m_sampleMs <= 0)
    {
      return;
    }
    m_samplerPeak = 0;
    m_stopSampler.store(false, std::memory_order_relaxed);
    m_samplerThread = std::thread([this]() {
      while(!m_stopSampler.load(std::memory_order_relaxed))
      {
        long long rss = -1, peak = -1;
        readProcRss(rss, peak);
        if(rss > m_samplerPeak)
        {
          m_samplerPeak = rss;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(m_sampleMs));
      }
    });
  }

  /// Stop the sampler and report the peak RSS during the sampled phase.
  void stopSampler(const std::string& label)
  {
    if(!m_enabled || m_sampleMs <= 0)
    {
      return;
    }
    m_stopSampler.store(true, std::memory_order_relaxed);
    if(m_samplerThread.joinable())
    {
      m_samplerThread.join();
    }
    long long rss = -1, peak = -1;  // final read to catch the tail
    readProcRss(rss, peak);
    if(rss > m_samplerPeak)
    {
      m_samplerPeak = rss;
    }

    long long maxVal, sumVal;
    int maxRank;
    reduceBytes(m_samplerPeak, maxVal, maxRank, sumVal);
    if(my_rank == 0)
    {
      SLIC_INFO(axom::fmt::format(
        "[mem] {}: peak RSS during phase (sampled @ {} ms): max/rank={} (rank {}), total={}",
        label,
        m_sampleMs,
        humanBytes(maxVal),
        maxRank,
        humanBytes(sumVal)));
    }
  }

  /// Take a snapshot on every rank, reduce it, and print an aggregate (rank 0).
  void report(const std::string& label)
  {
    if(!m_enabled)
    {
      return;
    }

    long long rss = -1, peakRss = -1;
    readProcRss(rss, peakRss);

    long long mallocLive = -1, mallocArena = -1;
#if defined(__GLIBC__)
  #if defined(__GLIBC_PREREQ) && __GLIBC_PREREQ(2, 33)
    struct mallinfo2 mi = mallinfo2();  // size_t fields: safe above 2 GiB
    mallocLive = static_cast<long long>(mi.uordblks) + static_cast<long long>(mi.hblkhd);
    mallocArena = static_cast<long long>(mi.arena);
  #else
    struct mallinfo mi = mallinfo();  // NOTE: int fields saturate above ~2 GiB
    mallocLive = static_cast<long long>(mi.uordblks) + static_cast<long long>(mi.hblkhd);
    mallocArena = static_cast<long long>(mi.arena);
  #endif
#endif

    long long umpireCur = -1, umpireHwm = -1;
#if defined(AXOM_USE_UMPIRE)
    if(m_umpireAllocatorId >= 0)
    {
      auto& rm = umpire::ResourceManager::getInstance();
      umpire::Allocator alloc = rm.getAllocator(m_umpireAllocatorId);
      umpireCur = static_cast<long long>(alloc.getCurrentSize());
      umpireHwm = static_cast<long long>(alloc.getHighWatermark());
    }
#endif

    long long rssMax, rssSum, peakMax, peakSum, liveMax, liveSum, arenaMax, arenaSum;
    long long umpCurMax, umpCurSum, umpHwmMax, umpHwmSum;
    int rssMaxRank, peakMaxRank, liveMaxRank, arenaMaxRank, umpCurMaxRank, umpHwmMaxRank;
    reduceBytes(rss, rssMax, rssMaxRank, rssSum);
    reduceBytes(peakRss, peakMax, peakMaxRank, peakSum);
    reduceBytes(mallocLive, liveMax, liveMaxRank, liveSum);
    reduceBytes(mallocArena, arenaMax, arenaMaxRank, arenaSum);
    reduceBytes(umpireCur, umpCurMax, umpCurMaxRank, umpCurSum);
    reduceBytes(umpireHwm, umpHwmMax, umpHwmMaxRank, umpHwmSum);

    if(my_rank == 0)
    {
      std::string msg = axom::fmt::format(
        "[mem] {}  (total over {} ranks | max on one rank)\n"
        "         RSS           : {:>11} | {:>11} (rank {})\n"
        "         RSS peak      : {:>11} | {:>11} (rank {})\n"
        "         malloc live   : {:>11} | {:>11} (rank {})\n"
        "         malloc arena  : {:>11} | {:>11} (rank {})",
        label,
        num_ranks,
        humanBytes(rssSum),
        humanBytes(rssMax),
        rssMaxRank,
        humanBytes(peakSum),
        humanBytes(peakMax),
        peakMaxRank,
        humanBytes(liveSum),
        humanBytes(liveMax),
        liveMaxRank,
        humanBytes(arenaSum),
        humanBytes(arenaMax),
        arenaMaxRank);
#if defined(AXOM_USE_UMPIRE)
      if(umpHwmMax >= 0)
      {
        msg += axom::fmt::format(
          "\n         umpire current: {:>11} | {:>11} (rank {})"
          "\n         umpire hi-water: {:>10} | {:>11} (rank {})",
          humanBytes(umpCurSum),
          humanBytes(umpCurMax),
          umpCurMaxRank,
          humanBytes(umpHwmSum),
          humanBytes(umpHwmMax),
          umpHwmMaxRank);
      }
#endif
      SLIC_INFO(msg);
    }
  }

private:
  bool m_enabled {false};
  int m_sampleMs {0};
  int m_umpireAllocatorId {-1};
  std::atomic<bool> m_stopSampler {false};
  std::thread m_samplerThread;
  long long m_samplerPeak {0};
};

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

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  initializeLogger();

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  axom::CLI::App app {"Driver for distributed distance query"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(my_rank == 0)
    {
      retval = app.exit(e);
    }

    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();

    exit(retval);
  }

#if defined(AXOM_USE_UMPIRE)
  //---------------------------------------------------------------------------
  // Memory resource.  For testing, choose device memory if appropriate.
  //---------------------------------------------------------------------------
  const std::string umpireResourceName = params.policy == RuntimePolicy::seq ? "HOST" :
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    params.policy == RuntimePolicy::omp ? "HOST"
                                        :
  #endif
  #if defined(UMPIRE_ENABLE_DEVICE)
                                        "DEVICE"
  #elif defined(UMPIRE_ENABLE_UM)
                                                                             "UM"
  #elif defined(UMPIRE_ENABLE_PINNED)
                                                                             "PINNED"
  #else
                                                                             "HOST"
  #endif
    ;
  auto& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator umpireAllocator = rm.getAllocator(umpireResourceName);
#endif

  // Storage for meshes.
  sidre::DataStore dataStore;

  //---------------------------------------------------------------------------
  // Load query mesh and generate a particle mesh over its nodes
  // These will be used to query the closest points on the object mesh(es)
  //---------------------------------------------------------------------------

  QueryMeshWrapper queryMeshWrapper(dataStore.getRoot()->createGroup("queryMesh", true),
                                    params.meshFile);

  if(params.isVerbose())
  {
    queryMeshWrapper.getParticleMesh().printMeshSizeStats("Query mesh");
  }
  slic::flushStreams();

  const size_t spatialDim = queryMeshWrapper.getParticleMesh().dimension();

  //---------------------------------------------------------------------------
  // Object (second) mesh
  //---------------------------------------------------------------------------

  ObjectMeshWrapper objectMeshWrapper(dataStore.getRoot()->createGroup("object_mesh", true),
                                      params.objectMeshFile);

  if(params.isVerbose())
  {
    objectMeshWrapper.getParticleMesh().printMeshSizeStats("Object mesh");
  }
  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Initialize spatial index for object points, and run query
  //---------------------------------------------------------------------------

  int globalObjectPointCount = objectMeshWrapper.getParticleMesh().numPoints();
  MPI_Allreduce(MPI_IN_PLACE, &globalObjectPointCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  auto init_str =
    banner(axom::fmt::format("Initializing BVH tree over {} object points", globalObjectPointCount));

  axom::utilities::Timer initTimer(false);
  axom::utilities::Timer queryTimer(false);

  // Convert blueprint representation from sidre to conduit
  conduit::Node objectMeshNode;
  if(objectMeshWrapper.getParticleMesh().numPoints() > 0)
  {
    objectMeshWrapper.getBlueprintGroup()->createNativeLayout(objectMeshNode);
  }

  // Put sidre data into Conduit Node.
  conduit::Node queryMeshNode;
  queryMeshWrapper.getBlueprintGroup()->createNativeLayout(queryMeshNode);

  // Create distributed closest point query object and set some parameters
  quest::DistributedClosestPoint query;
  query.setRuntimePolicy(params.policy);
#if defined(AXOM_USE_UMPIRE)
  query.setAllocatorID(umpireAllocator.getId());
#endif
  query.setMpiCommunicator(MPI_COMM_WORLD, true);
  query.setVerbosity(params.isVerbose());
  query.setDistanceThreshold(params.distThreshold);
  query.setDynamicDistanceFiltering(params.dynamicDistanceFiltering);
  // To test support for single-domain format, use single-domain when possible.
  query.setObjectMesh(objectMeshNode.number_of_children() == 1 ? objectMeshNode[0] : objectMeshNode,
                      objectMeshWrapper.getTopologyName());

  // Optional memory instrumentation around the index build and the query.
  int memUmpireId = -1;
#if defined(AXOM_USE_UMPIRE)
  memUmpireId = umpireAllocator.getId();
#endif
  const bool trackMem = params.trackMemory || params.trimAfterQuery || params.sampleMemoryMs > 0;
  MemoryProbe memProbe(trackMem, params.sampleMemoryMs, memUmpireId);
  memProbe.report("baseline (meshes read, before BVH)");

  // Build the spatial index over the object on each rank
  SLIC_INFO(init_str);
  slic::flushStreams();
  initTimer.start();
  query.generateBVHTree();
  initTimer.stop();

  memProbe.report("after generateBVHTree (object index built)");

  // Run the distributed closest point query over the nodes of the computational mesh
  // To test support for single-domain format, use single-domain when possible.
  slic::flushStreams();
  memProbe.resetPeak();  // make the post-query VmHWM reflect only the query phase
  memProbe.startSampler();
  queryTimer.start();
  query.computeClosestPoints(
    queryMeshNode.number_of_children() == 1 ? queryMeshNode[0] : queryMeshNode,
    queryMeshWrapper.getTopologyName());
  queryTimer.stop();
  memProbe.stopSampler("computeClosestPoints");

  memProbe.report("after computeClosestPoints");

  // Optionally return free arena memory to the OS and re-measure.
  // If RSS drops here while "malloc live" was already back at baseline,
  // the retained memory was glibc arena, not a leak.
  // (This is a diagnostic; the library-side fix would trim inside computeClosestPoints after releasing transfer buffers.)
  if(params.trimAfterQuery)
  {
#if defined(__GLIBC__)
    ::malloc_trim(0);
    memProbe.report("after malloc_trim(0)");
#else
    SLIC_WARNING("--trim-after-query requested but malloc_trim is glibc-only; skipping.");
#endif
  }

  auto getDoubleMinMax = [](double inVal, double& minVal, double& maxVal, double& sumVal) {
    MPI_Allreduce(&inVal, &minVal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &maxVal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &sumVal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  };

  // Output some timing stats
  {
    double minInit, maxInit, sumInit;
    getDoubleMinMax(initTimer.elapsedTimeInSec(), minInit, maxInit, sumInit);

    double minQuery, maxQuery, sumQuery;
    getDoubleMinMax(queryTimer.elapsedTimeInSec(), minQuery, maxQuery, sumQuery);

    SLIC_INFO(
      axom::fmt::format("Initialization with policy {} took {{avg:{}, min:{}, max:{}}} seconds",
                        axom::runtime_policy::s_policyToName.at(params.policy),
                        sumInit / num_ranks,
                        minInit,
                        maxInit));
    SLIC_INFO(axom::fmt::format("Query with policy {} took {{avg:{}, min:{}, max:{}}} seconds",
                                axom::runtime_policy::s_policyToName.at(params.policy),
                                sumQuery / num_ranks,
                                minQuery,
                                maxQuery));
  }
  slic::flushStreams();
  queryMeshWrapper.update_closest_points(queryMeshNode);

  if(spatialDim == 2)
  {
    computeDistancesAndDirections<2>(queryMeshWrapper.getParticleMesh(),
                                     "cp_coords",
                                     "cp_index",
                                     "distance",
                                     "direction");
  }
  else if(spatialDim == 3)
  {
    computeDistancesAndDirections<3>(queryMeshWrapper.getParticleMesh(),
                                     "cp_coords",
                                     "cp_index",
                                     "distance",
                                     "direction");
  }

  // queryMeshNode.print();
  queryMeshNode.reset();

  queryMeshWrapper.saveMesh(params.distanceFile);

  SLIC_INFO("Normal exit.");

  finalizeLogger();
  MPI_Finalize();

  return 0;
}
