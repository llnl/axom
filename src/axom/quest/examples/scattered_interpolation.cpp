// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file scattered_interpolation.cpp
 * \brief This examples uses a Delaunay mesh to perform scattered data
 * interpolation on fields defined over a point set
 */

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/sidre.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <memory>
#include <cmath>

namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;

namespace internal
{
namespace blueprint
{
/**
 *  \brief Simple wrapper to a blueprint particle mesh
 *
 *  Given a sidre Group, creates the stubs for a mesh blueptint particle mesh
 */
struct PointMesh
{
public:
  explicit PointMesh(sidre::Group* group = nullptr,
                     const std::string& coordset = "coords",
                     const std::string& topology = "mesh")
    : m_group(group)
    , m_meshName(topology)
    , m_coordName(coordset)
  {
    setBlueprintGroup(m_group, coordset, topology);
  }
  /// Gets the root group for this mesh blueprint
  sidre::Group* rootGroup() const { return m_group; }
  /// Gets the parent group for the blueprint coordinate set
  sidre::Group* coordsGroup() const { return m_coordsGroup; }
  /// Gets the parent group for the blueprint mesh topology
  sidre::Group* topoGroup() const { return m_topoGroup; }

  std::string meshName() const { return m_meshName; }
  std::string coordsName() const { return m_coordName; }

  /// Loads the point mesh from a conduit hierarchy
  bool loadFromConduitNode(const conduit::Node& mesh_node, bool verboseOutput = false)
  {
    conduit::Node info;
    if(!conduit::blueprint::verify("mesh", mesh_node, info))
    {
      SLIC_ERROR("Invalid blueprint for particle mesh: \n" << info.to_yaml());
      return false;
    }
    else if(verboseOutput)
    {
      SLIC_INFO("mesh_node info: " << info.to_yaml());
    }

    m_group->destroyViewsAndData();
    m_group->destroyGroups();

    auto domains = conduit::blueprint::mesh::domains(mesh_node);
    SLIC_ERROR_IF(domains.size() != 1,
                  "loadFromConduitNode only currently supports loading "
                  "blueprint files with a single mesh.");

    m_group->importConduitTree(*domains[0]);
    if(!this->isValid())
    {
      return false;
    }

    // Attempt to rename the mesh from the file
    if(domains[0]->has_child("sidre_group_name"))
    {
      m_group->rename(domains[0]->child("sidre_group_name").name());
    }
    else if(!domains[0]->name().empty())
    {
      m_group->rename(domains[0]->name());
    }

    // Log some debug info
    if(verboseOutput)
    {
      SLIC_INFO(axom::fmt::format("Root: {{name:{}, path:{}, num groups:{}}}",
                                  m_group->getName(),
                                  m_group->getPathName(),
                                  m_group->getNumGroups()));

      m_group->printTree(1, std::cout);
    }

    // extract coordset name and group
    {
      auto* coords = m_group->getGroup("coordsets");

      SLIC_ERROR_IF(coords->getNumGroups() != 1,
                    "PointMesh::loadFromConduitNode() only currently supports "
                    "meshes with a single coordset");

      m_coordName = coords->getGroup(0)->getName();

      if(verboseOutput)
      {
        SLIC_INFO("Coords name is: " << m_coordName);
      }
      SLIC_ERROR_IF(!coords->hasGroup(m_coordName), "Missing required coordinates: " << m_coordName);

      m_coordsGroup = coords->getGroup(m_coordName);
    }

    // extract mesh/topology name and group
    {
      auto* topos = m_group->getGroup("topologies");

      SLIC_ERROR_IF(topos->getNumGroups() != 1,
                    "PointMesh::loadFromConduitNode() only currently supports "
                    "meshes with a single topology");

      m_meshName = topos->getGroup(0)->getName();

      if(verboseOutput)
      {
        SLIC_INFO("Mesh name is: " << m_meshName);
      }
      SLIC_ERROR_IF(!topos->hasGroup(m_meshName), "Missing required mesh: " << m_meshName);

      m_topoGroup = topos->getGroup(m_meshName);
    }

    // extract fields group
    {
      m_fieldsGroup =
        m_group->hasGroup("fields") ? m_group->getGroup("fields") : m_group->createGroup("fields");
    }

    return true;
  }

  /// Returns true if points have been added to the particle mesh
  bool hasPoints() const { return m_coordsGroup != nullptr && m_coordsGroup->hasGroup("values"); }

  /// Returns the number of points in the particle mesh
  int numPoints() const
  {
    return hasPoints() ? m_coordsGroup->getView("values/x")->getNumElements() : 0;
  }

  int dimension() const
  {
    if(hasPoints())
    {
      return m_coordsGroup->hasView("values/z") ? 3 : (m_coordsGroup->hasView("values/y") ? 2 : 1);
    }
    return 0;
  }

  /**
   * Sets the parent group for the entire mesh and sets up the blueprint stubs
   * for the "coordset", "topologies", "fields" and "state"
   */
  void setBlueprintGroup(sidre::Group* group,
                         const std::string& coordset = "coords",
                         const std::string& topology = "mesh")
  {
    // TODO: Ensure that we delete previous hierarchy if it existed

    m_group = group;
    m_meshName = topology;
    m_coordName = coordset;

    if(m_group != nullptr)
    {
      createBlueprintStubs(coordset, topology);
    }
  }

  /// Set the coordinate data from an array of primal Points, templated on the dimension
  template <int NDIMS>
  void setPoints(const axom::Array<primal::Point<double, NDIMS>>& pts)
  {
    SLIC_ASSERT_MSG(m_group != nullptr, "Must set blueprint group before setPoints()");

    const int SZ = pts.size();

    // create views into a shared buffer for the coordinates, with stride NDIMS
    auto* buf = m_group->getDataStore()->createBuffer(sidre::DOUBLE_ID, NDIMS * SZ)->allocate();
    switch(NDIMS)
    {
    case 3:
      m_coordsGroup->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, NDIMS);
      m_coordsGroup->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, NDIMS);
      m_coordsGroup->createView("values/z")->attachBuffer(buf)->apply(SZ, 2, NDIMS);
      break;
    case 2:
      m_coordsGroup->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, NDIMS);
      m_coordsGroup->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, NDIMS);
      break;
    default:
      m_coordsGroup->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, NDIMS);
      break;
    }

    // copy coordinate data into the buffer
    const std::size_t nbytes = sizeof(double) * SZ * NDIMS;
    axom::copy(buf->getVoidPtr(), pts.data(), nbytes);

    // set the default connectivity
    sidre::Array<int> arr(m_topoGroup->createView("elements/connectivity"), SZ, SZ);
    for(int i = 0; i < SZ; ++i)
    {
      arr[i] = i;
    }
  }

  template <typename T>
  void registerNodalScalarField(const std::string& fieldName)
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot register a field with the BlueprintParticleMesh "
                    "before adding points");

    auto* fld = m_fieldsGroup->createGroup(fieldName);
    fld->createViewString("association", "vertex");
    fld->createViewString("topology", m_topoGroup->getName());
    fld->createViewAndAllocate("values", sidre::detail::SidreTT<T>::id, numPoints());
  }

  template <typename T>
  void registerNodalVectorField(const std::string& fieldName)
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot register a field with the BlueprintParticleMesh "
                    "before adding points");

    const int SZ = numPoints();
    const int DIM = dimension();

    auto* fld = m_fieldsGroup->createGroup(fieldName);
    fld->createViewString("association", "vertex");
    fld->createViewString("topology", m_topoGroup->getName());

    // create views into a shared buffer for the coordinates, with stride NDIMS
    auto* buf =
      m_group->getDataStore()->createBuffer(sidre::detail::SidreTT<T>::id, DIM * SZ)->allocate();
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

  bool hasField(const std::string& fieldName) const { return m_fieldsGroup->hasGroup(fieldName); }

  template <typename T>
  axom::ArrayView<T> getNodalScalarField(const std::string& fieldName) const
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot extract a field from the BlueprintParticleMesh "
                    "before adding points");

    T* data = hasField(fieldName)
      ? static_cast<T*>(
          m_fieldsGroup->getView(axom::fmt::format("{}/values", fieldName))->getVoidPtr())
      : nullptr;

    return axom::ArrayView<T>(data, numPoints());
  }

  template <typename T>
  axom::ArrayView<T> getNodalVectorField(const std::string& fieldName) const
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot extract a field from the BlueprintParticleMesh "
                    "before adding points");

    T* data = hasField(fieldName)
      ? static_cast<T*>(
          m_fieldsGroup->getView(axom::fmt::format("{}/values/x", fieldName))->getVoidPtr())
      : nullptr;

    return axom::ArrayView<T>(data, numPoints());
  }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    conduit::Node mesh_node;
    m_group->createNativeLayout(mesh_node);

    bool success = true;
    conduit::Node info;
    if(!conduit::blueprint::verify("mesh", mesh_node, info))
    {
      SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
      success = false;
    }

    return success;
  }

  /// Outputs the mesh to disk
  void saveMesh(const std::string& outputMesh, const std::string& protocol)
  {
    auto* ds = m_group->getDataStore();

    const std::string ext = (protocol.find("hdf5") != std::string::npos)
      ? "hdf5"
      : (protocol.find("json") != std::string::npos ? "json" : protocol);
    const std::string fname = axom::fmt::format("{}.{}", outputMesh, ext);

    // generate metadata and write out as "root" file
    {
      auto* rootfile_grp = ds->getRoot()->hasGroup("rootfile_data")
        ? ds->getRoot()->getGroup("rootfile_data")
        : ds->getRoot()->createGroup("rootfile_data");

      rootfile_grp->createViewScalar("number_of_files", 1);
      rootfile_grp->createViewString("file_pattern", fname);
      rootfile_grp->createViewScalar("number_of_trees", 1);
      rootfile_grp->createViewString("tree_pattern", "/");
      rootfile_grp->createViewString("protocol/name", protocol);
      rootfile_grp->createViewString("protocol/version", "0.8.2");

      ds->generateBlueprintIndex(
        m_group->getPathName(),
        m_group->getName(),
        axom::fmt::format("rootfile_data/blueprint_index/{}", m_group->getName()),
        1);

      rootfile_grp->save(axom::fmt::format("{}.root", outputMesh), "json");

      ds->getRoot()->destroyGroup("rootfile_data");
    }

    ds->getRoot()->save(fname, protocol);
  }

  /// Returns a list of registered field names
  std::vector<std::string> getFieldNames() const
  {
    std::vector<std::string> names;

    for(auto idx = m_fieldsGroup->getFirstValidGroupIndex(); sidre::indexIsValid(idx);
        idx = m_fieldsGroup->getNextValidGroupIndex(idx))
    {
      names.push_back(m_fieldsGroup->getGroup(idx)->getName());
    }

    return names;
  }

private:
  /// Creates blueprint stubs for this mesh
  void createBlueprintStubs(const std::string& coords, const std::string& topo)
  {
    SLIC_ASSERT(m_group != nullptr);

    m_coordsGroup = m_group->createGroup("coordsets")->createGroup(coords);
    m_coordsGroup->createViewString("type", "explicit");
    m_coordsGroup->createGroup("values");

    m_topoGroup = m_group->createGroup("topologies")->createGroup(topo);
    m_topoGroup->createViewString("coordset", coords);
    m_topoGroup->createViewString("type", "unstructured");
    m_topoGroup->createViewString("elements/shape", "point");

    m_fieldsGroup = m_group->createGroup("fields");
  }

private:
  sidre::Group* m_group;

  sidre::Group* m_coordsGroup;
  sidre::Group* m_topoGroup;
  sidre::Group* m_fieldsGroup;

  std::string m_meshName;
  std::string m_coordName;
};

}  // namespace blueprint
}  // namespace internal

/// Struct to parse and contain command line arguments
struct Input
{
  std::string outputFile {"scattered_interpolation"};
  std::string inputFile;

  bool verboseOutput {false};
  int numRandPoints {20};
  int numQueryPoints {20};
  int dimension {2};
  std::vector<double> boundsMin;
  std::vector<double> boundsMax;
  std::string outputProtocol = sidre::Group::getDefaultIOProtocol();

  const std::set<std::string> s_validProtocols {"json",
                                                "sidre_json"
#ifdef AXOM_USE_HDF5
                                                ,
                                                "sidre_hdf5"
#endif
  };

public:
  bool hasInputMesh() const { return !inputFile.empty(); }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_flag("-v,--verbose", verboseOutput)
      ->description("Increases the output verbosity while running the application")
      ->capture_default_str();

    // Options for input data
    // Either provide `-n` and `-d`; or `-i` (input mesh)
    auto input_grp = app.add_option_group("Input",
                                          "Parameters associated with input data.\n"
                                          "If an input mesh is provided, it will override the "
                                          "`-n` and `-d` options.");

    input_grp->add_option("-n,--nrandpt", numRandPoints)
      ->description("The number of points to generate for the input mesh")
      ->capture_default_str();

    input_grp->add_option("-d,--dim", dimension)
      ->description(
        "The dimension of the mesh. 2 for triangle mesh in 2D; "
        "3 for tetrahedral mesh in 3D")
      ->capture_default_str();

    input_grp->add_option("-i,--infile", inputFile)
      ->description("Input point mesh with associated scalar fields")
      ->check(axom::CLI::ExistingFile);

    // Options for defining the query data
    auto query_grp =
      app.add_option_group("Query",
                           "Parameters associated with query data. "
                           "Query box defaults to unit square/cube when not provided.\n"
                           "If no input mesh is provided, the input points "
                           "will lie in query bounding box.");

    query_grp->add_option("-q,--nquerypt", numQueryPoints)
      ->description("The number of query points")
      ->capture_default_str();

    auto* minbb = query_grp->add_option("--min", boundsMin)
                    ->description("Min bounds for query box (x,y[,z])")
                    ->expected(2, 3);
    auto* maxbb = query_grp->add_option("--max", boundsMax)
                    ->description("Max bounds for query box (x,y[,z])")
                    ->expected(2, 3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    // Options for outputting the mesh
    auto output_grp = app.add_option_group("Output", "Parameters associated with output");
    output_grp->add_option("-o,--outfile", outputFile)
      ->description("The output file")
      ->capture_default_str();

    output_grp->add_option("-p,--protocol", outputProtocol)
      ->description("Set the output protocol for sidre point meshes")
      ->capture_default_str()
      ->check(axom::CLI::IsMember(s_validProtocols));

    // --

    app.get_formatter()->column_width(45);

    // could throw an exception
    app.parse(argc, argv);

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
      inputFile: '{}'
      nquerypt: {}
      bounding box min: {{{}}}
      bounding box max: {{{}}}
      outfile = '{}'
      output protocol = '{}'
    }})",
                                dimension,
                                numRandPoints,
                                inputFile,
                                numQueryPoints,
                                axom::fmt::join(boundsMin, ", "),
                                axom::fmt::join(boundsMax, ", "),
                                outputFile,
                                outputProtocol));

    axom::slic::setLoggingMsgLevel(verboseOutput ? axom::slic::message::Debug
                                                 : axom::slic::message::Info);
  }
};

template <int DIM>
axom::Array<primal::Point<double, DIM>> generatePts(int numPts,
                                                    const std::vector<double>& bb_min,
                                                    const std::vector<double>& bb_max)
{
  using axom::utilities::random_real;

  using PointType = typename primal::Point<double, DIM>;
  using BoundingBox = typename primal::BoundingBox<double, DIM>;

  axom::Array<PointType> pts(numPts, numPts);

  BoundingBox bbox {PointType(bb_min.data()), PointType(bb_max.data())};

  // generate random points within bounding box
  for(int i = 0; i < numPts; ++i)
  {
    PointType& pt = pts[i];
    for(int d = 0; d < DIM; ++d)
    {
      pt[d] = random_real(bbox.getMin()[d], bbox.getMax()[d]);
    }
  }

  return pts;
}

template <int DIM>
void initializeInputMesh(Input& params, internal::blueprint::PointMesh& inputMesh)
{
  using BBoxType = primal::BoundingBox<double, DIM>;

  if(params.hasInputMesh())
  {
    // Load mesh from file into inputMesh
    conduit::Node mesh;
    conduit::relay::io::blueprint::read_mesh(params.inputFile, mesh);

    bool validMesh = inputMesh.loadFromConduitNode(mesh, params.verboseOutput);
    SLIC_ERROR_IF(!validMesh, axom::fmt::format("Could not read mesh '{}'", params.inputFile));

    // recompute and reset bounding box based on input mesh
    const int nPts = inputMesh.numPoints();
    auto coords = quest::detail::InterleavedOrStridedPoints<double, DIM>(
      inputMesh.coordsGroup()->getGroup("values"));

    BBoxType bbox;
    for(int idx = 0; idx < nPts; ++idx)
    {
      bbox.addPoint(coords[idx]);
    }
    SLIC_INFO("Input mesh bounding box: " << bbox);

    params.boundsMin.clear();
    params.boundsMax.clear();
    for(int d = 0; d < DIM; ++d)
    {
      params.boundsMin.push_back(bbox.getMin()[d]);
      params.boundsMax.push_back(bbox.getMax()[d]);
    }

    // Reset other input parameters
    params.dimension = DIM;
    params.numRandPoints = 0;
  }
  else
  {
    inputMesh.setPoints(generatePts<DIM>(params.numRandPoints, params.boundsMin, params.boundsMax));

    // Extract coordinate positions as scalar fields
    const int nPts = inputMesh.numPoints();
    auto coords = quest::detail::InterleavedOrStridedPoints<double, DIM>(
      inputMesh.coordsGroup()->getGroup("values"));

    inputMesh.registerNodalScalarField<double>("pos_x");
    inputMesh.registerNodalScalarField<double>("pos_y");
    if(DIM > 2)
    {
      inputMesh.registerNodalScalarField<double>("pos_z");
    }

    auto pos_x = inputMesh.getNodalScalarField<double>("pos_x");
    auto pos_y = inputMesh.getNodalScalarField<double>("pos_y");
    auto pos_z = (DIM > 2) ? inputMesh.getNodalScalarField<double>("pos_z")
                           : inputMesh.getNodalScalarField<double>("");

    for(int i = 0; i < nPts; ++i)
    {
      auto pt = coords[i];
      pos_x[i] = pt[0];
      pos_y[i] = pt[1];
      if(DIM > 2)
      {
        pos_z[i] = pt[DIM - 1];
      }
    }
  }
}

template <int DIM>
void initializeQueryMesh(const Input& params,
                         internal::blueprint::PointMesh& queryMesh,
                         std::vector<std::string>&& fieldNames)
{
  queryMesh.setPoints(generatePts<DIM>(params.numQueryPoints, params.boundsMin, params.boundsMax));

  // Register and allocate required fields for query mesh
  queryMesh.registerNodalScalarField<axom::IndexType>("cell_idx");

  for(const auto& fld : fieldNames)
  {
    queryMesh.registerNodalScalarField<double>(fld);
  }
}

/**
 * Compares interpolation using ScatteredInterpolation::getInterpolationWeights()
 * against results from ScattteredInterpolation::interpolateField().
 * The latter stores its results in the blueprint mesh, while the former
 * returns the interpolation weights and indexes for the user to apply themselves.
 */
template <int DIM>
bool checkInterpolation(std::unique_ptr<quest::ScatteredInterpolation<DIM>>& scattered_interp,
                        const internal::blueprint::PointMesh& inputMesh,
                        const internal::blueprint::PointMesh& queryMesh)
{
  using axom::utilities::isNearlyEqual;
  using quest::detail::InterleavedOrStridedPoints;

  using InterpIndices = primal::Point<axom::IndexType, DIM + 1>;
  using InterpWeights = primal::Point<double, DIM + 1>;
  using PointArray = InterleavedOrStridedPoints<double, DIM>;

  constexpr double EPS = 1e-8;

  bool interpolationsAgree = true;

  // Get arrays over the input and query points
  PointArray inputPts(inputMesh.coordsGroup()->getGroup("values"));
  PointArray queryPts(queryMesh.coordsGroup()->getGroup("values"));

  // Create temporary arrays of interpolation variables
  const int nQueryPts = queryMesh.numPoints();
  axom::Array<InterpIndices> interp_indices(nQueryPts, nQueryPts);
  axom::Array<InterpWeights> interp_weights(nQueryPts, nQueryPts);
  axom::Array<bool> interp_valid(nQueryPts, nQueryPts);
  for(int i = 0; i < nQueryPts; ++i)
  {
    interp_valid[i] =
      scattered_interp->getInterpolationWeights(queryPts[i], interp_indices[i], interp_weights[i]);
  }

  // Interpolate points on each field and check against previously interpolated values
  for(const auto& fld : inputMesh.getFieldNames())
  {
    auto inputField = inputMesh.getNodalScalarField<double>(fld);
    auto queryField = queryMesh.getNodalScalarField<double>(fld);

    for(int i = 0; i < nQueryPts; ++i)
    {
      const double query_val = queryField[i];

      if(!interp_valid[i])
      {
        if(!std::isnan(query_val))
        {
          SLIC_WARNING(
            axom::fmt::format("Bad interpolation: Query point {} had value {} for field `{}` "
                              "using `ScatteredInterpolation::interpolateField()` but was not "
                              "found for `ScatteredInterpolation::getInterpolationWeights()`",
                              i,
                              query_val,
                              fld));
          interpolationsAgree = false;
        }
      }
      else
      {
        // perform interpolation
        double val = 0.;
        auto& indices = interp_indices[i];
        auto& weights = interp_weights[i];
        for(int d = 0; d < indices.dimension(); ++d)
        {
          val += inputField[indices[d]] * weights[d];
        }

        // check against value computed via ScatteredInterpolation::interpolateField()
        // stored in query mesh's blueprint fields
        if(!isNearlyEqual(query_val, val, EPS))
        {
          SLIC_WARNING(
            axom::fmt::format("Bad interpolation: Query point {} had value {} for field `{}` "
                              "using `ScatteredInterpolation::interpolateField()` but {} when "
                              "interpolating using "
                              "`ScatteredInterpolation::getInterpolationWeights()`",
                              i,
                              query_val,
                              fld,
                              val));
          interpolationsAgree = false;
        }

        // Interpoate position and check against actual query point
        primal::Point<double, DIM> interp_pt;
        for(int d = 0; d < indices.dimension(); ++d)
        {
          interp_pt.array() += inputPts[indices[d]].array() * weights[d];
        }

        if(!isNearlyEqual(primal::squared_distance(interp_pt, queryPts[i]), EPS))
        {
          SLIC_WARNING(
            axom::fmt::format("Bad interpolation: When interpolating position of query point {} "
                              "@ {}, got {} with distance {}.",
                              i,
                              queryPts[i],
                              interp_pt,
                              sqrt(primal::squared_distance(interp_pt, queryPts[i]))));

          interpolationsAgree = false;
        }
      }
    }
  }

  return interpolationsAgree;
}

int main(int argc, char** argv)
{
  // Initialize the SLIC logger
  axom::slic::SimpleLogger logger(axom::slic::message::Info);
  // axom::slic::enableAbortOnWarning();

  // Initialize default parameters and update with command line arguments:
  Input params;
  axom::CLI::App app {"Scattered interpolation over 2D or 3D point collection"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  /// Set up the input mesh

  sidre::DataStore input_ds;

  const std::string input_coords_name = "input_coords";
  const std::string input_mesh_name = "input_mesh";
  auto* input_pt_mesh_group = input_ds.getRoot()->createGroup("input_point_mesh");
  internal::blueprint::PointMesh inputMesh(input_pt_mesh_group, input_coords_name, input_mesh_name);

  // Initialize the mesh points
  switch(params.dimension)
  {
  case 2:
    initializeInputMesh<2>(params, inputMesh);
    break;
  case 3:
    initializeInputMesh<3>(params, inputMesh);
    break;
  }

  // Write input mesh to file
  if(params.verboseOutput)
  {
    std::string file = params.outputFile + "_input_mesh";
    SLIC_INFO(
      axom::fmt::format("Writing input mesh to '{}/{}'", axom::utilities::filesystem::getCWD(), file));
    inputMesh.saveMesh(file, params.outputProtocol);
  }

  /// Set up the query mesh

  sidre::DataStore query_ds;
  const std::string query_coords_name = "query_coords";
  const std::string query_mesh_name = "query_mesh";
  auto* query_pt_mesh_group = query_ds.getRoot()->createGroup("query_point_mesh");
  internal::blueprint::PointMesh queryMesh(query_pt_mesh_group, query_coords_name, query_mesh_name);

  // Initialize the mesh points
  switch(params.dimension)
  {
  case 2:
    initializeQueryMesh<2>(params, queryMesh, inputMesh.getFieldNames());
    break;
  case 3:
    initializeQueryMesh<3>(params, queryMesh, inputMesh.getFieldNames());
    break;
  }
  const int numQueryPts = queryMesh.numPoints();
  const int numFields = inputMesh.getFieldNames().size();

  // Convert blueprint meshes from Sidre to Conduit
  conduit::Node bp_input, bp_query;
  inputMesh.rootGroup()->createNativeLayout(bp_input);
  queryMesh.rootGroup()->createNativeLayout(bp_query);

  /// Run the query
  std::unique_ptr<quest::ScatteredInterpolation<2>> scattered_2d;
  std::unique_ptr<quest::ScatteredInterpolation<3>> scattered_3d;

  // Initialize the ScatteredInterpolation instance (templated on the dimension)
  switch(params.dimension)
  {
  case 2:
    scattered_2d = std::make_unique<quest::ScatteredInterpolation<2>>();
    break;
  case 3:
    scattered_3d = std::make_unique<quest::ScatteredInterpolation<3>>();
    break;
  }

  // Generate the Delaunay complex
  int numVerts = 0, numSimps = 0;
  std::string bboxStr;
  axom::utilities::Timer timer(true);
  switch(params.dimension)
  {
  case 2:
    scattered_2d->buildTriangulation(bp_input, inputMesh.coordsName());
    numVerts = scattered_2d->numVertices();
    numSimps = scattered_2d->numSimplices();
    bboxStr = axom::fmt::format("{}", scattered_2d->boundingBox());
    break;
  case 3:
    scattered_3d->buildTriangulation(bp_input, inputMesh.coordsName());
    numVerts = scattered_3d->numVertices();
    numSimps = scattered_3d->numSimplices();
    bboxStr = axom::fmt::format("{}", scattered_3d->boundingBox());
    break;
  }
  timer.stop();
  SLIC_INFO(
    axom::fmt::format("It took {} seconds to create a Delaunay complex with {} "
                      "points, {} {} and bounding box {}. "
                      "Insertion rate of {:.1f} points per second.",
                      timer.elapsedTimeInSec(),
                      numVerts,
                      numSimps,
                      params.dimension == 2 ? "triangles" : "tetrahedra",
                      bboxStr,
                      numVerts / timer.elapsedTimeInSec()));

  // Dump the Delaunay complex to disk as a vtk file
  if(params.verboseOutput)
  {
    switch(params.dimension)
    {
    case 2:
      scattered_2d->exportDelaunayComplex(bp_input, "delaunay_2d.vtk");
      break;
    case 3:
      scattered_3d->exportDelaunayComplex(bp_input, "delaunay_3d.vtk");
      break;
    }
  }

  // Find the simplices containing each of the query points
  timer.reset();
  timer.start();
  switch(params.dimension)
  {
  case 2:
    scattered_2d->locatePoints(bp_query, query_coords_name);
    break;
  case 3:
    scattered_3d->locatePoints(bp_query, query_coords_name);
    break;
  }
  timer.stop();
  SLIC_INFO(
    axom::fmt::format("It took {} seconds to locate {} points. Query rate of "
                      "{:.1f} points per second.",
                      timer.elapsedTimeInSec(),
                      numQueryPts,
                      numQueryPts / timer.elapsedTimeInSec()));

  // Perform the interpolation on each of the fields
  timer.reset();
  timer.start();
  switch(params.dimension)
  {
  case 2:
    for(const auto& fld : inputMesh.getFieldNames())
    {
      scattered_2d->interpolateField(bp_query, query_coords_name, bp_input, fld, fld);
    }
    break;
  case 3:
    for(const auto& fld : inputMesh.getFieldNames())
    {
      scattered_3d->interpolateField(bp_query, query_coords_name, bp_input, fld, fld);
    }
    break;
  }
  timer.stop();
  SLIC_INFO(
    axom::fmt::format("It took {} seconds to interpolate the data on {} fields. "
                      "Interpolation rate of {:.1f} points per second.",
                      timer.elapsedTimeInSec(),
                      numFields,
                      numQueryPts * numFields / timer.elapsedTimeInSec()));

  // Check interpolation using local interpolation from weights
  timer.reset();
  timer.start();
  switch(params.dimension)
  {
  case 2:
    checkInterpolation<2>(scattered_2d, inputMesh, queryMesh);
    break;
  case 3:
    checkInterpolation<3>(scattered_3d, inputMesh, queryMesh);
    break;
  }
  timer.stop();
  SLIC_INFO(
    axom::fmt::format("It took {} seconds to check interpolation values by "
                      "interpolating from weights. "
                      "Interpolation rate of {:.1f} points per second.",
                      timer.elapsedTimeInSec(),
                      numQueryPts * numFields / timer.elapsedTimeInSec()));

  // Write query mesh to file
  if(params.verboseOutput)
  {
    std::string file = params.outputFile + "_output_mesh";
    SLIC_INFO(axom::fmt::format("Writing interpolated point mesh to '{}/{}'",
                                axom::utilities::filesystem::getCWD(),
                                file));

    queryMesh.saveMesh(file, params.outputProtocol);
  }

  return 0;
}
