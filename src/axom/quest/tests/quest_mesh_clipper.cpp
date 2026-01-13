// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_mesh_clipper.cpp
 * \brief Test clipping codes built around MeshClipper class.
 * 3D only.  Extensible to 2D when we have 2D clipping.
 */

// Axom includes
#include "axom/config.hpp"

#if !defined(AXOM_USE_CONDUIT)
  #error Shaping functionality requires Axom to be configured with Conduit
#endif

#if !defined(AXOM_USE_SIDRE)
  #error Shaping functionality requires Axom to be have sidre enabled
#endif

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/klee.hpp"
#include "axom/quest.hpp"
#include "axom/quest/detail/clipping/HexClipper.hpp"
#include "axom/quest/detail/clipping/Plane3DClipper.hpp"
#include "axom/quest/detail/clipping/SphereClipper.hpp"
#include "axom/quest/detail/clipping/TetClipper.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_relay_io_blueprint.hpp"
#include "conduit_utils.hpp"

#include <math.h>

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// RAJA
#if !defined(AXOM_USE_RAJA)
  #error quest_mesh_clipper example require RAJA
#endif
#include "RAJA/RAJA.hpp"

// C/C++ includes
#include <string>
#include <vector>
#include <memory>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace sidre = axom::sidre;

//------------------------------------------------------------------------------

using RuntimePolicy = axom::runtime_policy::Policy;

#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
[[maybe_unused]] static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType =
  conduit::DataType::INT64_ID;
#else
[[maybe_unused]] static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType =
  conduit::DataType::INT32_ID;
#endif

/// Struct to parse and store the input parameters
// Some parameters are used to override defaults.
struct Input
{
public:
  std::string outputFile;

  std::vector<double> center;
  double radius {-1.0};
  double radius2 {-0.3};
  double length {-2.0};
  std::vector<double> direction;

  // Shape transformation parameters
  std::vector<double> scaleFactors;

  // Inline mesh parameters
  std::vector<double> boxMins {-2, -2, -2};
  std::vector<double> boxMaxs {2, 2, 2};
  std::vector<int> boxResolution {20, 20, 20};
  int getBoxDim() const
  {
    auto d = boxResolution.size();
    SLIC_ASSERT(boxMins.size() == d);
    SLIC_ASSERT(boxMaxs.size() == d);
    return int(d);
  }
  int getBoxCellCount() const { return boxResolution[0] * boxResolution[1] * boxResolution[2]; }

  // The shape to run.
  std::vector<std::string> testGeom;
  // The shapes this example is set up to run.
  const std::set<std::string> availableShapes {"tet", "hex", "sphere", "plane"};  // More geometries to come.

  RuntimePolicy policy {RuntimePolicy::seq};
  int refinementLevel {7};
  double weldThresh {1e-9};
  std::string annotationMode {"none"};

  std::string backgroundMaterial;

  int screenLevel = -1;

  // clang-format off
  enum class MeshType { bpSidre = 0, bpConduit = 1 };
  const std::map<std::string, MeshType> meshTypeChoices
    { {"bpSidre", MeshType::bpSidre} , {"bpConduit", MeshType::bpConduit} };
  // clang-format on
  MeshType meshType {MeshType::bpSidre};
  bool useBlueprintSidre() { return meshType == MeshType::bpSidre; }
  bool useBlueprintConduit() { return meshType == MeshType::bpConduit; }

private:
  bool m_verboseOutput {false};

public:
  bool isVerbose() const { return m_verboseOutput; }

  /// @brief Return volume of input box mesh
  double boxMeshVolume() const
  {
    primal::Vector<double, 3> x {boxMaxs[0] - boxMins[0], 0, 0};
    primal::Vector<double, 3> y {0, boxMaxs[1] - boxMins[1], 0};
    primal::Vector<double, 3> z {0, 0, boxMaxs[2] - boxMins[2]};
    double volume = primal::Vector<double, 3>::scalar_triple_product(x, y, z);
    return volume;
  }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("--screenLevel", screenLevel)
      ->description("Developer feature for MeshClipper.")
      ->capture_default_str();

    app.add_option("-o,--outputFile", outputFile)->description("Path to output file(s)");

    app.add_flag("-v,--verbose,!--no-verbose", m_verboseOutput)
      ->description("Enable/disable verbose output, including SLIC_DEBUG")
      ->capture_default_str();

    app.add_option("--meshType", meshType)
      ->description("Type of computational mesh to shape on")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(meshTypeChoices));

    app.add_option("-t,--weld-threshold", weldThresh)
      ->description("Threshold for welding")
      ->check(axom::CLI::NonNegativeNumber)
      ->capture_default_str();

    app.add_option("-s,--testGeom", testGeom)
      ->description(
        "The geometry(s) to run.  Specifying multiple shapes will override scaling and "
        "translations to shrink shapes and shift them to individual octants of the mesh.")
      ->check(axom::CLI::IsMember(availableShapes))
      ->delimiter(',')
      ->expected(1, 60);

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    app.add_option("--center", center)
      ->description("Center of sphere or base of cone/cyl/SOR (x,y[,z]) shape")
      ->expected(2, 3);

    app.add_option("--radius", radius)
      ->description("Radius of sphere or cylinder shape")
      ->check(axom::CLI::PositiveNumber);

    app.add_option("--length", length)
      ->description("Length of cone/cyl/SOR shape, avg length of hex.")
      ->check(axom::CLI::PositiveNumber);

    app.add_option("--dir", direction)
      ->description(
        "Direction of axis of rotation (cone/cyl/SOR (x,y[,z])), or rotated "
        "x-axis (hex, tet, tetmesh, and sphere), or positive normal direction "
        "(plane).")
      ->expected(2, 3);

    app.add_option("--radius2", radius2)
      ->description("Second radius of cone shape")
      ->check(axom::CLI::PositiveNumber);

    app.add_option("--scale", scaleFactors)
      ->description("Scale factor to apply to shape (x,y[,z])")
      ->expected(2, 3)
      ->check(axom::CLI::PositiveNumber);

    // use either an input mesh file or a simple inline Cartesian mesh
    {
      auto* inline_mesh_subcommand = app.add_subcommand("inline_mesh")
                                       ->description("Options for setting up a simple inline mesh")
                                       ->fallthrough();

      inline_mesh_subcommand->add_option("--min", boxMins)
        ->description("Min bounds for box mesh (x,y[,z])")
        ->expected(2, 3)
        ->required();
      inline_mesh_subcommand->add_option("--max", boxMaxs)
        ->description("Max bounds for box mesh (x,y[,z])")
        ->expected(2, 3)
        ->required();

      inline_mesh_subcommand->add_option("--res", boxResolution)
        ->description("Resolution of the box mesh (i,j[,k])")
        ->expected(2, 3)
        ->required();
    }

    app.add_option("--background-material", backgroundMaterial)
      ->description("Sets the name of the background material");

    // parameters that only apply to the intersection method
    {
      auto* intersection_options =
        app.add_option_group("intersection", "Options related to intersection-based queries");

      intersection_options->add_option("-r, --refinements", refinementLevel)
        ->description("Number of refinements to perform for revolved contour")
        ->capture_default_str()
        ->check(axom::CLI::NonNegativeNumber);

      std::stringstream pol_sstr;
      pol_sstr << "Set runtime policy for intersection-based sampling method.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
      pol_sstr << "\nSet to 'seq' or 0 to use the RAJA sequential policy.";
  #ifdef AXOM_USE_OPENMP
      pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
  #endif
  #ifdef AXOM_USE_CUDA
      pol_sstr << "\nSet to 'cuda' or 2 to use the RAJA CUDA policy.";
  #endif
  #ifdef AXOM_USE_HIP
      pol_sstr << "\nSet to 'hip' or 3 to use the RAJA HIP policy.";
  #endif
#endif

      intersection_options->add_option("-p, --policy", policy, pol_sstr.str())
        ->capture_default_str()
        ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));
    }
    app.get_formatter()->column_width(50);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug : slic::message::Info);
  }
};  // struct Input
Input params;

/************************************************************
 * Shared variables.
 ************************************************************/

const std::string topoName = "mesh";
const std::string matsetName = "matset";
const std::string coordsetName = "coords";
int cellCount = -1;
// Translation to individual octants (override) when running multiple shapes.
// Exception: the plane always placed at the center of the box mesh
// to facilitate finding its exact overlap volume.
const double tDist = 0.9;  // Bias toward origin to help keep shape inside domain.
std::vector<axom::NumericArray<double, 3>> translations {{tDist, tDist, -tDist},
                                                         {-tDist, tDist, -tDist},
                                                         {-tDist, -tDist, -tDist},
                                                         {tDist, -tDist, -tDist},
                                                         {tDist, tDist, tDist},
                                                         {-tDist, tDist, tDist},
                                                         {-tDist, -tDist, tDist},
                                                         {tDist, -tDist, tDist}};
int translationIdx = 0;  // To track what translations have been used.

std::map<std::string, int> geomReps;  // Repetitions of the geometry.
std::map<std::string, double> exactGeomVols;
std::map<std::string, double> errorToleranceRel;  // Relative error tolerance.
std::map<std::string, double> errorToleranceAbs;  // Absolute error tolerance.
double vScale = 1.0;                              // Volume scale due to geometry scale.

// Start property for all 3D shapes.
axom::klee::TransformableGeometryProperties startProp {axom::klee::Dimensions::Three,
                                                       axom::klee::LengthUnit::unspecified};

// Add scale operator if specified by input parameters.
void addScaleOperator(axom::klee::CompositeOperator& compositeOp)
{
  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  if(!params.scaleFactors.empty())
  {
    std::shared_ptr<axom::klee::Scale> scaleOp =
      std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                          params.scaleFactors[1],
                                          params.scaleFactors[2],
                                          startProp);
    compositeOp.addOperator(scaleOp);
  }
}

// Add translate operator.
void addTranslateOperator(axom::klee::CompositeOperator& compositeOp)
{
  if(params.testGeom.size() > 1)
  {
    const axom::NumericArray<double, 3>& shifts =
      translations[(translationIdx++) % translations.size()];
    primal::Vector3D shift({shifts[0], shifts[1], shifts[2]});
    auto translateOp = std::make_shared<axom::klee::Translation>(shift, startProp);
    compositeOp.addOperator(translateOp);
  }
  else
  {
    // Use zero shift as a smoke test.
    primal::Vector3D shift({0, 0, 0});
    auto translateOp = std::make_shared<axom::klee::Translation>(shift, startProp);
    compositeOp.addOperator(translateOp);
  }
}

// Add operator to rotate x-axis to params.direction, if it is given.
void addRotateOperator(axom::klee::CompositeOperator& compositeOp)
{
  if(!params.direction.empty())
  {
    static const primal::Point3D center {0.0, 0.0, 0.0};
    static const primal::Vector3D x {1.0, 0.0, 0.0};
    primal::Vector3D rotateTo(params.direction.data());
    // Note that the rotation matrix is not unique.
    primal::Vector3D a = rotateTo.unitVector();
    primal::Vector3D u;  // Rotation vector, the cross product of x and a.
    axom::numerics::cross_product(x.data(), a.data(), u.data());
    double angle = asin(u.norm()) * 180 / M_PI;

    auto rotateOp = std::make_shared<axom::klee::Rotation>(angle, center, u, startProp);
    compositeOp.addOperator(rotateOp);
  }
}

// Computational mesh in different forms, initialized in main
axom::sidre::Group* compMeshGrp = nullptr;
axom::sidre::Group* compMeshGrpOnHost = nullptr;
std::shared_ptr<conduit::Node> compMeshNode;

axom::sidre::Group* createBoxMesh(axom::sidre::Group* meshGrp)
{
  using BBox3D = primal::BoundingBox<double, 3>;
  using Pt3D = primal::Point<double, 3>;
  auto res = axom::NumericArray<int, 3>(params.boxResolution.data());
  auto bbox = BBox3D(Pt3D(params.boxMins.data()), Pt3D(params.boxMaxs.data()));
  axom::quest::util::make_unstructured_blueprint_box_mesh_3d(meshGrp,
                                                             bbox,
                                                             res,
                                                             topoName,
                                                             coordsetName,
                                                             params.policy);
#if defined(AXOM_DEBUG)
  conduit::Node meshNode, info;
  meshGrp->createNativeLayout(meshNode);
  SLIC_ASSERT(conduit::blueprint::mesh::verify(meshNode, info));
#endif

  // State group is optional to blueprint, and we don't use it, but mint checks for it.
  meshGrp->createGroup("state");

  return meshGrp;
}

/// \brief Utility function to initialize the logger
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream {nullptr};

#ifdef AXOM_USE_MPI
  int num_ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  SLIC_ERROR_IF(num_ranks > 1, "Sorry, this test is serial.");
  if(num_ranks > 1)
  {
    std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
    const int RLIMIT = 8;
    logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
    logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
  }
  else
#endif  // AXOM_USE_MPI
  {
    std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
    logStream = new slic::GenericOutputStream(&std::cout, fmt);
  }

  slic::addStreamToAllMsgLevels(logStream);
}

/// \brief Utility function to finalize the logger
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

/*
 * For the test shapes, try to get good volume with compact shape
 * that stays in domain when rotated (else volume check is invalid).
 */

axom::klee::Geometry createGeom_Sphere(const std::string& geomName)
{
  Point3D center = params.center.empty() ? Point3D {0, 0, 0} : Point3D {params.center.data()};
  double radius = params.radius < 0 ? 1.0 : params.radius;
  axom::primal::Sphere<double, 3> sphere {center, radius};

  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp);

  const axom::IndexType levelOfRefinement = params.refinementLevel;
  axom::klee::Geometry sphereGeometry(prop, sphere, levelOfRefinement, compositeOp);
  exactGeomVols[geomName] = vScale * 4. / 3 * M_PI * radius * radius * radius;
  errorToleranceRel[geomName] = 1e-3;
  // Tolerance should account for discretization errors.
  errorToleranceRel[geomName] = params.refinementLevel <= 5 ? 0.0015 : 0.0001;
  errorToleranceAbs[geomName] = errorToleranceRel[geomName] * exactGeomVols[geomName];

  return sphereGeometry;
}

axom::klee::Geometry createGeom_Tet(const std::string& geomName)
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  // Tetrahedron at origin.
  const double len = params.length < 0 ? 1.55 : params.length;
  const Point3D a {Point3D::NumericArray {.8, 0., -1.} * len};
  const Point3D b {Point3D::NumericArray {-.8, 1, -1.} * len};
  const Point3D c {Point3D::NumericArray {-.8, -1, -1.} * len};
  const Point3D d {Point3D::NumericArray {0., 0., +1.} * len};
  const primal::Tetrahedron<double, 3> tet {a, b, c, d};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp);
  exactGeomVols[geomName] = vScale * tet.volume();
  errorToleranceRel[geomName] = 1e-12;
  errorToleranceAbs[geomName] = errorToleranceRel[geomName] * exactGeomVols[geomName];

  axom::klee::Geometry tetGeometry(prop, tet, compositeOp);

  return tetGeometry;
}

axom::klee::Geometry createGeom_Hex(const std::string& geomName)
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  const double md = params.length < 0 ? 0.82 : params.length / 2;
  const double lg = 1.2 * md;
  const double sm = 0.8 * md;
  const Point3D p {-lg, -md, -sm};
  const Point3D q {+lg, -md, -sm};
  const Point3D r {+lg, +md, -sm};
  const Point3D s {-lg, +md, -sm};
  const Point3D t {-lg, -md, +sm};
  const Point3D u {+lg, -md, +sm};
  const Point3D v {+lg, +md, +sm};
  const Point3D w {-lg, +md, +sm};
  const primal::Hexahedron<double, 3> hex {p, q, r, s, t, u, v, w};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp);
  exactGeomVols[geomName] = vScale * hex.volume();
  errorToleranceRel[geomName] = 0.000075;
  errorToleranceAbs[geomName] = 0.0003;

  axom::klee::Geometry hexGeometry(prop, hex, compositeOp);

  return hexGeometry;
}

axom::klee::Geometry createGeom_Plane(const std::string& geomName)
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  // Create a plane crossing center of mesh.  No matter the normal,
  // it cuts the mesh in half.
  Point3D center {0.5 *
                  (axom::NumericArray<double, 3>(params.boxMins.data()) +
                   axom::NumericArray<double, 3>(params.boxMaxs.data()))};
  primal::Vector<double, 3> normal = params.direction.empty()
    ? primal::Vector3D {1.0, 0.0, 0.0}
    : primal::Vector3D {params.direction.data()}.unitVector();
  const primal::Plane<double, 3> plane {normal, center, true};

  axom::klee::Geometry planeGeometry(prop, plane, {nullptr});

  // Exact mesh overlap volume, assuming plane passes through center of box mesh.
  using Pt3D = primal::Point<double, 3>;
  Pt3D lower(params.boxMins.data());
  Pt3D upper(params.boxMaxs.data());
  auto diag = upper.array() - lower.array();
  double meshVolume = diag[0] * diag[1] * diag[2];
  exactGeomVols[geomName] = 0.5 * meshVolume;
  errorToleranceRel[geomName] = 1e-6;
  errorToleranceAbs[geomName] = 1e-8;

  return planeGeometry;
}

/*!
  @brief Return the element volumes as a sidre::View containing
  the volumes in an array.

  If it doesn't exist, allocate and compute it.
  \post The volume data is in the blueprint field \c volFieldName.
*/
template <typename ExecSpace>
axom::sidre::View* getElementVolumes(
  sidre::Group* meshGrp,
  const std::string& volFieldName = std::string("elementVolumes"))
{
  using XS = axom::execution_space<ExecSpace>;
  using HexahedronType = axom::primal::Hexahedron<double, 3>;

  axom::sidre::View* volSidreView = nullptr;

  const auto fieldPath = axom::fmt::format("fields/{}", volFieldName);
  if(meshGrp->hasGroup(fieldPath))
  {
    sidre::Group* fieldGrp = meshGrp->getGroup(fieldPath);
    volSidreView = fieldGrp->getView("values");
  }
  else
  {
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(meshGrp, topoName);

    constexpr int NUM_VERTS_PER_HEX = 8;
    constexpr int NUM_COMPS_PER_VERT = 3;
    constexpr double ZERO_THRESHOLD = 1.e-10;

    /*
      Get vertex coordinates.  We use UnstructuredMesh for this,
      so get it on host first then transfer to device if needed.
    */
    auto* connData = meshGrp->getGroup("topologies")
                       ->getGroup(topoName)
                       ->getGroup("elements")
                       ->getView("connectivity");
    SLIC_ASSERT(connData->getNode().dtype().id() == conduitDataIdOfAxomIndexType);

    conduit::Node coordNode;
    meshGrp->getGroup("coordsets")->getGroup(coordsetName)->createNativeLayout(coordNode);
    const conduit::Node& coordValues = coordNode.fetch_existing("values");
    axom::IndexType vertexCount = coordValues["x"].dtype().number_of_elements();
    bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
    int stride = isInterleaved ? NUM_COMPS_PER_VERT : 1;
    axom::StackArray<axom::ArrayView<const double>, 3> coordArrays {
      axom::ArrayView<const double>(coordValues["x"].as_double_ptr(), {vertexCount}, stride),
      axom::ArrayView<const double>(coordValues["y"].as_double_ptr(), {vertexCount}, stride),
      axom::ArrayView<const double>(coordValues["z"].as_double_ptr(), {vertexCount}, stride)};

    const axom::IndexType* connPtr = connData->getArray();
    SLIC_ASSERT(connPtr != nullptr);
    axom::ArrayView<const axom::IndexType, 2> conn(connPtr, cellCount, NUM_VERTS_PER_HEX);
    axom::Array<Point3D> vertCoords(cellCount * NUM_VERTS_PER_HEX,
                                    cellCount * NUM_VERTS_PER_HEX,
                                    XS::allocatorID());
    auto vertCoordsView = vertCoords.view();

    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellIdx) {
        // Get the indices of this element's vertices
        auto verts = conn[cellIdx];

        // Get the coordinates for the vertices
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIdx = cellIdx * NUM_VERTS_PER_HEX + j;
          for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
          {
            vertCoordsView[vertIdx][k] = coordArrays[k][verts[j]];
            // vertCoordsView[vertIdx][k] = mesh.getNodeCoordinate(verts[j], k);
          }
        }
      });

    // Set vertex coords to zero if within threshold.
    // (I don't know why we do this.  I'm following examples.)
    axom::ArrayView<double> flatCoordsView((double*)vertCoords.data(),
                                           vertCoords.size() * Point3D::dimension());
    assert(flatCoordsView.size() == cellCount * NUM_VERTS_PER_HEX * 3);
    axom::for_all<ExecSpace>(
      cellCount * 3,
      AXOM_LAMBDA(axom::IndexType i) {
        if(axom::utilities::isNearlyEqual(flatCoordsView[i], 0.0, ZERO_THRESHOLD))
        {
          flatCoordsView[i] = 0.0;
        }
      });

    // Initialize hexahedral elements.
    axom::Array<HexahedronType> hexes(cellCount, cellCount, meshGrp->getDefaultAllocatorID());
    auto hexesView = hexes.view();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellIdx) {
        // Set each hexahedral element vertices
        hexesView[cellIdx] = HexahedronType();
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIndex = (cellIdx * NUM_VERTS_PER_HEX) + j;
          auto& hex = hexesView[cellIdx];
          hex[j] = vertCoordsView[vertIndex];
        }
      });  // end of loop to initialize hexahedral elements and bounding boxes

    // Allocate and populate cell volumes.
    axom::sidre::Group* fieldGrp = meshGrp->createGroup(fieldPath);
    fieldGrp->createViewString("topology", topoName);
    fieldGrp->createViewString("association", "element");
    fieldGrp->createViewString("volume_dependent", "true");
    volSidreView =
      fieldGrp->createViewAndAllocate("values", axom::sidre::detail::SidreTT<double>::id, cellCount);
    axom::IndexType shape2d[] = {cellCount, 1};
    volSidreView->reshapeArray(2, shape2d);
    axom::ArrayView<double> volView(volSidreView->getData(), volSidreView->getNumElements());
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellIdx) { volView[cellIdx] = hexesView[cellIdx].volume(); });
  }

  return volSidreView;
}

template <typename ExecSpace>
double sumMaterialVolumesImpl(sidre::Group* meshGrp, const std::string& material)
{
  conduit::Node meshNode;
  meshGrp->createNativeLayout(meshNode);
#if defined(AXOM_DEBUG)
  // Conduit can verify Blueprint mesh, but only if data is on host.
  if(axom::execution_space<axom::SEQ_EXEC>::usesAllocId(meshGrp->getDefaultAllocatorID()))
  {
    conduit::Node info;
    conduit::blueprint::mesh::verify(meshNode, info);
    SLIC_ASSERT(conduit::blueprint::mesh::verify(meshNode, info));
  }
#endif
  std::string topoPath = "topologies/" + topoName;
  conduit::Node& topoNode = meshNode.fetch_existing(topoPath);
  const int cellCount = conduit::blueprint::mesh::topology::length(topoNode);

  // Get cell volumes from meshGrp.
  const std::string volsName = "vol_" + material;
  axom::sidre::View* elementVols = getElementVolumes<ExecSpace>(meshGrp, volsName);
  axom::ArrayView<double> elementVolsView(elementVols->getData(), elementVols->getNumElements());

  // Get material volume fractions
  const auto vfFieldName = "vol_frac_" + material;
  const auto vfFieldValuesPath = "fields/" + vfFieldName + "/values";
  axom::sidre::View* volFrac = meshGrp->getView(vfFieldValuesPath);
  axom::ArrayView<double> volFracView(volFrac->getArray(), cellCount);

  using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<ReducePolicy, double> localVol(0);
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType i) { localVol += volFracView[i] * elementVolsView[i]; });

  double globalVol = localVol.get();
#ifdef AXOM_USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &globalVol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return globalVol;
}

double sumMaterialVolumes(sidre::Group* meshGrp, const std::string& material)
{
  double rval = 0.0;
  if(params.policy == RuntimePolicy::seq)
  {
    rval = sumMaterialVolumesImpl<axom::SEQ_EXEC>(meshGrp, material);
  }
#if defined(AXOM_USE_OPENMP)
  if(params.policy == RuntimePolicy::omp)
  {
    rval = sumMaterialVolumesImpl<axom::OMP_EXEC>(meshGrp, material);
  }
#endif
#if defined(AXOM_USE_CUDA)
  if(params.policy == RuntimePolicy::cuda)
  {
    rval = sumMaterialVolumesImpl<axom::CUDA_EXEC<256>>(meshGrp, material);
  }
#endif
#if defined(AXOM_USE_HIP)
  if(params.policy == RuntimePolicy::hip)
  {
    rval = sumMaterialVolumesImpl<axom::HIP_EXEC<256>>(meshGrp, material);
  }
#endif
  return rval;
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

/// Write blueprint mesh to disk
void saveMesh(const sidre::Group& mesh, const std::string& filename)
{
  AXOM_ANNOTATE_SCOPE("save mesh (sidre)");

  axom::sidre::DataStore ds;
  const sidre::Group* meshOnHost = &mesh;
  if(mesh.getDefaultAllocatorID() != axom::execution_space<axom::SEQ_EXEC>::allocatorID())
  {
    meshOnHost =
      ds.getRoot()->deepCopyGroup(&mesh, axom::execution_space<axom::SEQ_EXEC>::allocatorID());
  }
  conduit::Node tmpMesh;
  meshOnHost->createNativeLayout(tmpMesh);
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

//!@brief Fill a sidre array View with a value.
// No error checking.
template <typename T>
void fillSidreViewData(axom::sidre::View* view, const T& value)
{
  double* valuesPtr = view->getData<T*>();
  switch(params.policy)
  {
#if defined(AXOM_USE_CUDA)
  case RuntimePolicy::cuda:
    axom::for_all<axom::CUDA_EXEC<256>>(
      view->getNumElements(),
      AXOM_LAMBDA(axom::IndexType i) { valuesPtr[i] = value; });
    break;
#endif
#if defined(AXOM_USE_HIP)
  case RuntimePolicy::hip:
    axom::for_all<axom::HIP_EXEC<256>>(
      view->getNumElements(),
      AXOM_LAMBDA(axom::IndexType i) { valuesPtr[i] = value; });
    break;
#endif
#if defined(AXOM_USE_OMP)
  case RuntimePolicy::omp:
    axom::for_all<axom::OMP_EXEC>(
      view->getNumElements(),
      AXOM_LAMBDA(axom::IndexType i) { valuesPtr[i] = value; });
    break;
#endif
  case RuntimePolicy::seq:
  default:
    axom::for_all<axom::SEQ_EXEC>(
      view->getNumElements(),
      AXOM_LAMBDA(axom::IndexType i) { valuesPtr[i] = value; });
    break;
  }
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  const int my_rank = mpi_raii_wrapper.my_rank();

  initializeLogger();

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  axom::CLI::App app {"Driver for Klee shaping query"};

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
    finalizeLogger();

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    return retval;
  }

  if(params.testGeom.size() > 1)
  {
    if(params.scaleFactors.empty())
    {
      params.scaleFactors.resize(3, 1.0);
    }
    for(auto& f : params.scaleFactors) f *= 0.5;
    axom::StackArray<double, 3> tmpOutput {params.scaleFactors[0],
                                           params.scaleFactors[1],
                                           params.scaleFactors[2]};
    SLIC_WARNING(
      axom::fmt::format("Multiple test configurations specified.\n"
                        "Adding additional 0.5 scaling to shrink the geometries\n"
                        "and move them to individual octants so they don't overlap\n"
                        "with each other.  Final scaling: {}",
                        tmpOutput));
  }
  for(auto sf : params.scaleFactors)
  {
    vScale *= sf;
  }

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(params.annotationMode);

  /*
    Host allocator is for metadata and arrays that must be on host.
    Data allocator uses host or device, depending on the runtime policy.
  */
  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  int dataAllocId = axom::policyToDefaultAllocatorID(params.policy);

#if defined(AXOM_USE_UMPIRE)
  if(dataAllocId != axom::MALLOC_ALLOCATOR_ID)
  {
    // Use Umpire pool for performance benchmarking.
    constexpr size_t bytesPerCell = 100 * sizeof(double);
    size_t poolSize = params.getBoxCellCount() * bytesPerCell;
    umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
    const umpire::Allocator dataAllocator = rm.getAllocator(dataAllocId);
    const umpire::Allocator dataPoolAllocator =
      rm.makeAllocator<umpire::strategy::QuickPool>("data_pool", dataAllocator, poolSize);
    dataAllocId = dataPoolAllocator.getId();
    const std::string poolName = dataAllocator.getName() + "_pool";
    SLIC_INFO(axom::fmt::format("Using allocator pool id {}, '{}' with size {}",
                                dataAllocId,
                                poolName,
                                poolSize));
  }
#endif

  AXOM_ANNOTATE_BEGIN("quest clipping test");
  AXOM_ANNOTATE_BEGIN("init");

  // Storage for the some geometry meshes.
  sidre::DataStore ds;

  //---------------------------------------------------------------------------
  // Create shapes for the test
  //---------------------------------------------------------------------------
  axom::Array<std::shared_ptr<axom::quest::experimental::MeshClipperStrategy>> geomStrategies;
  geomStrategies.reserve(params.testGeom.size());
  SLIC_ERROR_IF(params.getBoxDim() != 3, "This example is only in 3D.");
  for(const auto& tg : params.testGeom)
  {
    if(geomReps.count(tg) == 0)
    {
      geomReps[tg] = 0;
    }
    std::string name = axom::fmt::format("{}.{}", tg, geomReps[tg]++);

    if(tg == "plane")
    {
      geomStrategies.push_back(
        std::make_shared<axom::quest::experimental::Plane3DClipper>(createGeom_Plane(name), name));
    }
    else if(tg == "hex")
    {
      geomStrategies.push_back(
        std::make_shared<axom::quest::experimental::HexClipper>(createGeom_Hex(name), name));
    }
    else if(tg == "sphere")
    {
      geomStrategies.push_back(
        std::make_shared<axom::quest::experimental::SphereClipper>(createGeom_Sphere(name), name));
    }
    else if(tg == "tet")
    {
      geomStrategies.push_back(
        std::make_shared<axom::quest::experimental::TetClipper>(createGeom_Tet(name), name));
    }
    // More geometries to come.
  }

  {
    SLIC_INFO(axom::fmt::format("{:-^80}", "Generating Blueprint mesh"));
    compMeshGrp = ds.getRoot()->createGroup("compMesh");
    compMeshGrp->setDefaultAllocator(dataAllocId);

    createBoxMesh(compMeshGrp);

    cellCount = params.getBoxCellCount();
  }

  //---------------------------------------------------------------------------
  // Initialize computational mesh.
  //---------------------------------------------------------------------------
  std::shared_ptr<quest::experimental::ShapeMesh> sMeshPtr;
  AXOM_ANNOTATE_BEGIN("setup shaping problem");
  if(params.useBlueprintSidre())
  {
    sMeshPtr = std::make_shared<quest::experimental::ShapeMesh>(params.policy,
                                                                dataAllocId,
                                                                compMeshGrp,
                                                                topoName,
                                                                matsetName);
  }
  if(params.useBlueprintConduit())
  {
    compMeshNode.reset(new conduit::Node);
    compMeshNode->set_allocator(sidre::ConduitMemory::axomAllocIdToConduit(dataAllocId));
    compMeshGrp->createNativeLayout(*compMeshNode);
    compMeshNode->set_allocator(sidre::ConduitMemory::axomAllocIdToConduit(dataAllocId));

    sMeshPtr = std::make_shared<quest::experimental::ShapeMesh>(params.policy,
                                                                dataAllocId,
                                                                *compMeshNode,
                                                                topoName,
                                                                matsetName);
  }
  quest::experimental::ShapeMesh& sMesh = *sMeshPtr;

  AXOM_ANNOTATE_END("setup shaping problem");

  // Compute and cache shared data so they are not associated with the first geometry.
  SLIC_INFO(axom::fmt::format("{:-^80}", "Precomputing mesh data"));
  sMesh.precomputeMeshData();

  AXOM_ANNOTATE_END("init");

  //---------------------------------------------------------------------------
  // Process each of the shapes
  //---------------------------------------------------------------------------

  int failCounts = 0;

  SLIC_INFO(axom::fmt::format("{:=^80}", "Shaping loop"));
  AXOM_ANNOTATE_BEGIN("clipping");
  for(axom::IndexType i = 0; i < geomStrategies.size(); ++i)
  {
    const auto geomName = geomStrategies[i]->name();
    const auto annotationName = "clipping:" + geomName;

    SLIC_INFO(axom::fmt::format("{:-^80}", axom::fmt::format("Processing geometry '{}'", geomName)));

    if(my_rank == 0)
    {
      std::cout << "Info for geometry '" << geomName << "':" << std::endl;
      geomStrategies[i]->info().print();
    }

    quest::experimental::MeshClipper clipper(sMesh, geomStrategies[i]);
    clipper.setVerbose(params.isVerbose());
    if(params.screenLevel >= 0)
    {
      clipper.setScreenLevel(params.screenLevel);
    }
    SLIC_INFO(axom::fmt::format("MeshClipper screen level: {}", clipper.getScreenLevel()));

    axom::Array<double> ovlap;
    AXOM_ANNOTATE_BEGIN(annotationName);
    clipper.clip(ovlap);
    AXOM_ANNOTATE_END(annotationName);

    clipper.logClippingStats();

    // Save volume fractions in mesh, for plotting and checking.
    sMesh.setMatsetFromVolume(geomStrategies[i]->name(), ovlap.view(), false);

    // Correctness check on overlap volume.
    if(!axom::execution_space<axom::SEQ_EXEC>::usesAllocId(ovlap.getAllocatorID()))
    {
      // Move to host for check.
      ovlap = axom::Array<double>(ovlap, hostAllocId);
    }
    auto ovlapView = ovlap.view();
    using reduce_policy = typename axom::execution_space<axom::SEQ_EXEC>::reduce_policy;
    RAJA::ReduceSum<reduce_policy, double> ovlapSumReduce(0.0);
    axom::for_all<axom::SEQ_EXEC>(
      ovlap.size(),
      AXOM_LAMBDA(axom::IndexType i) { ovlapSumReduce += ovlapView[i]; });
    double computedOverlapVol = ovlapSumReduce.get();
    double exactGeomVol = exactGeomVols[geomName];

    bool err = !axom::utilities::isNearlyEqualRelative(computedOverlapVol,
                                                       exactGeomVol,
                                                       errorToleranceRel.at(geomName),
                                                       errorToleranceAbs.at(geomName));
    failCounts += err;

    SLIC_INFO(axom::fmt::format("{:-^80}",
                                axom::fmt::format("Shape '{}' has volume {} vs {}, diff of {}, {}.",
                                                  geomName,
                                                  computedOverlapVol,
                                                  exactGeomVol,
                                                  computedOverlapVol - exactGeomVol,
                                                  (err ? "ERROR" : "OK"))));
  }
  AXOM_ANNOTATE_END("clipping");

  AXOM_ANNOTATE_BEGIN("setFreeVolumeFractions");
  sMesh.setFreeVolumeFractions("free");
  AXOM_ANNOTATE_END("setFreeVolumeFractions");

  /*
    Copy mesh to host check results and plot.
  */
  SLIC_INFO(axom::fmt::format("{:-^80}", "Copying mesh to host and write out"));

  AXOM_ANNOTATE_BEGIN("Copy results to host and write out");

  if(params.useBlueprintConduit())
  {
    compMeshGrpOnHost = ds.getRoot()->createGroup("onHost");
    compMeshGrpOnHost->setDefaultAllocator(hostAllocId);
    compMeshGrpOnHost->importConduitTree(*sMesh.getMeshAsConduit());
  }
  if(params.useBlueprintSidre())
  {
    if(sMesh.getMeshAsSidre()->getDefaultAllocatorID() != hostAllocId)
    {
      compMeshGrpOnHost = ds.getRoot()->createGroup("onHost");
      compMeshGrpOnHost->setDefaultAllocator(hostAllocId);
      compMeshGrpOnHost->deepCopyGroup(sMesh.getMeshAsSidre(), hostAllocId);
    }
    else
    {
      SLIC_ASSERT(sMesh.getMeshAsSidre() == compMeshGrp);
      compMeshGrpOnHost = compMeshGrp;
    }
  }

  compMeshNode.reset(new conduit::Node);
  compMeshGrpOnHost->createNativeLayout(*compMeshNode);

  /*
    Check blueprint validity.
  */

  conduit::Node whyNotValid;
  if(!conduit::blueprint::mesh::verify(*compMeshNode, whyNotValid))
  {
    SLIC_ERROR("Computational mesh is invalid after shaping:\n" + whyNotValid.to_summary_string());
  }

  /*
    Save meshes and fields
  */

  if(!params.outputFile.empty())
  {
    std::string fileName = params.outputFile + ".volfracs";
    saveMesh(*compMeshNode, fileName);
    SLIC_INFO(axom::fmt::format("{:-^80}", "Wrote output mesh " + fileName));
  }

  AXOM_ANNOTATE_END("Copy results to host and write out");

  /*
    Cleanup and exit
  */
  SLIC_INFO(axom::fmt::format("{:-^80}", ""));
  slic::flushStreams();

  AXOM_ANNOTATE_END("quest clipping test");

  SLIC_INFO(axom::fmt::format("exiting with failure count {}", failCounts));

  finalizeLogger();

  return failCounts;
}
