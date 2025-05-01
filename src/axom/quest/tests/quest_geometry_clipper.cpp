// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_shape_clipper.cpp
 * \brief Test clipping codes built around GeometryClipper class.
 * 3D only.  Extensible to 2D when we have 2D clipping.
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/klee.hpp"
#include "axom/quest.hpp"
#include "axom/core/WhereMacro.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#if !defined(AXOM_USE_CONDUIT)
  #error Shaping functionality requires Axom to be configured with Conduit
#endif

#include "conduit_blueprint.hpp"
#include "conduit_relay_io_blueprint.hpp"
#include "conduit_utils.hpp"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// RAJA
#if !defined(AXOM_USE_RAJA)
  #error quest_shape_clipper example require RAJA
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
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType = conduit::DataType::INT64_ID;
#else
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType = conduit::DataType::INT32_ID;
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
  std::vector<std::string> testGeom;  // {"tetmesh"};
  // The shapes this example is set up to run.
  const std::set<std::string>
    availableShapes {"tetmesh", "sphere", "cyl", "cone", "sor", "tet", "hex", "plane", "all"};

  RuntimePolicy policy {RuntimePolicy::seq};
  int outputOrder {2};
  int refinementLevel {7};
  double weldThresh {1e-9};
  double percentError {-1.};
  std::string annotationMode {"none"};

  std::string backgroundMaterial;

  // clang-format off
  enum class MeshType { bpSidre = 0, bpConduit = 1, mfem = 2 };
  const std::map<std::string, MeshType> meshTypeChoices
    { {"bpSidre", MeshType::bpSidre} , {"bpConduit", MeshType::bpConduit}, {"mfem", MeshType::mfem} };
  // clang-format on
  MeshType meshType {MeshType::bpSidre};
  bool useMfem() { return meshType == MeshType::mfem; }
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
    app.add_option("-o,--outputFile", outputFile)->description("Path to output file(s)");

    app.add_flag("-v,--verbose,!--no-verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("--meshType", meshType)
      ->description("Type of computational mesh to shape on")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(meshTypeChoices));

    app.add_option("-t,--weld-threshold", weldThresh)
      ->description("Threshold for welding")
      ->check(axom::CLI::NonNegativeNumber)
      ->capture_default_str();

    app.add_option("-e,--percent-error", percentError)
      ->description(
        "Percent error used for calculating curve refinement and revolved "
        "volume.\n"
        "If this value is provided then dynamic curve refinement will be used\n"
        "instead of segment-based curve refinement.")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();

    app.add_option("-s,--testGeom", testGeom)
      ->description("The shape(s) to run")
      ->check(axom::CLI::IsMember(availableShapes))
      ->delimiter(',')
      ->expected(1, 8);

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
void addTranslateOperator(axom::klee::CompositeOperator& compositeOp,
                          double shiftx,
                          double shifty,
                          double shiftz)
{
  primal::Vector3D shift({shiftx, shifty, shiftz});
  auto translateOp = std::make_shared<axom::klee::Translation>(shift, startProp);
  compositeOp.addOperator(translateOp);
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

const std::string topoName = "mesh";
const std::string matsetName = "matset";
const std::string coordsetName = "coords";
int cellCount = -1;
std::map<std::string, double> exactOverlapVols;
std::map<std::string, double> errorToleranceRel; // Relative error tolerance.
std::map<std::string, double> errorToleranceAbs; // Absolute error tolerance.

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

double volumeOfTetMesh(const axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& tetMesh)
{
  using TetType = axom::primal::Tetrahedron<double, 3>;
  axom::StackArray<axom::IndexType, 1> nodesShape {tetMesh.getNumberOfNodes()};
  axom::ArrayView<const double> x(tetMesh.getCoordinateArray(0), nodesShape);
  axom::ArrayView<const double> y(tetMesh.getCoordinateArray(1), nodesShape);
  axom::ArrayView<const double> z(tetMesh.getCoordinateArray(2), nodesShape);
  const axom::IndexType tetCount = tetMesh.getNumberOfCells();
  axom::Array<double> tetVolumes(tetCount, tetCount);
  double meshVolume = 0.0;
  for(axom::IndexType ic = 0; ic < tetCount; ++ic)
  {
    const axom::IndexType* nodeIds = tetMesh.getCellNodeIDs(ic);
    TetType tet;
    for(int j = 0; j < 4; ++j)
    {
      auto cornerNodeId = nodeIds[j];
      tet[j][0] = x[cornerNodeId];
      tet[j][1] = y[cornerNodeId];
      tet[j][2] = z[cornerNodeId];
    }
    meshVolume += tet.volume();
  }
  return meshVolume;
}

axom::klee::Geometry createGeom_Sphere()
{
  Point3D center = params.center.empty() ? Point3D {0, 0, 0} : Point3D {params.center.data()};
  double radius = params.radius < 0 ? 1.0 : params.radius;
  axom::primal::Sphere<double, 3> sphere {center, radius};

  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp, 1, 1, 1);

  const axom::IndexType levelOfRefinement = params.refinementLevel;
  axom::klee::Geometry sphereGeometry(prop, sphere, levelOfRefinement, compositeOp);
  exactOverlapVols["sphere"] = 4. / 3 * M_PI * radius * radius * radius;
  errorToleranceRel["sphere"] = 1e-3;
  errorToleranceAbs["sphere"] = 1e-5;

  return sphereGeometry;
}

axom::klee::Geometry createGeom_TetMesh(sidre::DataStore& ds)
{
  // Shape a tetrahedal mesh.
  sidre::Group* meshGroup = ds.getRoot()->createGroup("tetMesh");
  AXOM_UNUSED_VAR(meshGroup);  // variable is only referenced in debug configs
  const std::string topo = "mesh";
  const std::string coordset = "coords";
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tetMesh(3,
                                                                 axom::mint::CellType::TET,
                                                                 meshGroup,
                                                                 topo,
                                                                 coordset);

  double lll = params.length < 0 ? 1.17 : params.length;

  // Insert tets around origin.
  tetMesh.appendNode(-lll, -lll, -lll);
  tetMesh.appendNode(+lll, -lll, -lll);
  tetMesh.appendNode(-lll, +lll, -lll);
  tetMesh.appendNode(-lll, -lll, +lll);
  tetMesh.appendNode(+lll, +lll, +lll);
  tetMesh.appendNode(-lll, +lll, +lll);
  tetMesh.appendNode(+lll, +lll, -lll);
  tetMesh.appendNode(+lll, -lll, +lll);
  axom::IndexType conn0[4] = {0, 1, 2, 3};
  tetMesh.appendCell(conn0);
  axom::IndexType conn1[4] = {4, 5, 6, 7};
  tetMesh.appendCell(conn1);

  SLIC_ASSERT(axom::mint::blueprint::isValidRootGroup(meshGroup));

  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp, -1, 1, 1);

  axom::klee::Geometry tetMeshGeometry(prop, tetMesh.getSidreGroup(), topo, compositeOp);

  exactOverlapVols["tetmesh"] = volumeOfTetMesh(tetMesh);
  errorToleranceRel["tetmesh"] = 1e-6;
  errorToleranceAbs["tetmesh"] = 1e-8;

  return tetMeshGeometry;
}

axom::klee::Geometry createGeometry_Sor(axom::primal::Point<double, 3>& sorBase,
                                        axom::primal::Vector<double, 3>& sorDirection,
                                        axom::Array<double, 2>& discreteFunction,
                                        std::shared_ptr<axom::klee::CompositeOperator>& compositeOp)
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  const axom::IndexType levelOfRefinement = params.refinementLevel;
  axom::klee::Geometry sorGeometry(prop,
                                   discreteFunction,
                                   sorBase,
                                   sorDirection,
                                   levelOfRefinement,
                                   compositeOp);
  return sorGeometry;
}

double computeVolume_Sor(axom::Array<double, 2>& discreteFunction)
{
  using ConeType = axom::primal::Cone<double, 3>;
  axom::IndexType segmentCount = discreteFunction.shape()[0];
  double vol = 0.0;
  for(axom::IndexType s = 0; s < segmentCount - 1; ++s)
  {
    ConeType cone(discreteFunction(s, 1),
                  discreteFunction(s + 1, 1),
                  discreteFunction(s + 1, 0) - discreteFunction(s, 0));
    vol += cone.volume();
  }
  return vol;
}

axom::klee::Geometry createGeom_Sor()
{
  Point3D sorBase = params.center.empty() ? Point3D {0.0, 0.0, 0.0} : Point3D {params.center.data()};
  axom::primal::Vector<double, 3> sorDirection = params.direction.empty()
    ? primal::Vector3D {0.1, 0.2, 0.4}
    : primal::Vector3D {params.direction.data()};
  const int numIntervals = 5;
  // discreteFunction are discrete z-r pairs describing the function
  // to be rotated around the z axis.
  axom::Array<double, 2> discreteFunction({numIntervals + 1, 2}, axom::ArrayStrideOrder::ROW);
  double zLen = params.length < 0 ? 2.40: params.length;
  double zShift = -zLen / 2;
  double maxR = params.radius < 0 ? 1.10 : params.radius;
  double dz = zLen / numIntervals;
  discreteFunction(0, 0) = 0 * dz + zShift;
  discreteFunction(0, 1) = 0.0 * maxR;
  discreteFunction(1, 0) = 1 * dz + zShift;
  discreteFunction(1, 1) = 0.8 * maxR;
  discreteFunction(2, 0) = 2 * dz + zShift;
  discreteFunction(2, 1) = 0.4 * maxR;
  discreteFunction(3, 0) = 3 * dz + zShift;
  discreteFunction(3, 1) = 0.5 * maxR;
  discreteFunction(4, 0) = 4 * dz + zShift;
  discreteFunction(4, 1) = 1.0 * maxR;
  discreteFunction(5, 0) = 5 * dz + zShift;
  discreteFunction(5, 1) = 1.0 * maxR;

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addTranslateOperator(*compositeOp, -1, -1, 1);

  axom::klee::Geometry sorGeometry =
    createGeometry_Sor(sorBase, sorDirection, discreteFunction, compositeOp);

  exactOverlapVols["sor"] = computeVolume_Sor(discreteFunction);
  errorToleranceRel["sor"] = 0.011;
  errorToleranceAbs["sor"] = 0.013;

  return sorGeometry;
}

axom::klee::Geometry createGeom_Cylinder()
{
  Point3D sorBase = params.center.empty() ? Point3D {0.0, 0.0, 0.0} : Point3D {params.center.data()};
  axom::primal::Vector<double, 3> sorDirection = params.direction.empty()
    ? primal::Vector3D {0.1, 0.2, 0.4}
    : primal::Vector3D {params.direction.data()};
  // discreteFunction are discrete z-r pairs describing the function
  // to be rotated around the z axis.
  axom::Array<double, 2> discreteFunction({2, 2}, axom::ArrayStrideOrder::ROW);
  double radius = params.radius < 0 ? 0.695 : params.radius;
  double height = params.length < 0 ? 2.78 : params.length;
  discreteFunction(0, 0) = -height / 2;
  discreteFunction(0, 1) = radius;
  discreteFunction(1, 0) = height / 2;
  discreteFunction(1, 1) = radius;

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addTranslateOperator(*compositeOp, 1, -1, 1);

  axom::klee::Geometry sorGeometry =
    createGeometry_Sor(sorBase, sorDirection, discreteFunction, compositeOp);

  exactOverlapVols["cyl"] = computeVolume_Sor(discreteFunction);
  errorToleranceRel["cyl"] = 1e-3;
  errorToleranceAbs["cyl"] = 1e-5;

  return sorGeometry;
}

axom::klee::Geometry createGeom_Cone()
{
  Point3D sorBase = params.center.empty() ? Point3D {0.0, 0.0, 0.0} : Point3D {params.center.data()};
  axom::primal::Vector<double, 3> sorDirection = params.direction.empty()
    ? primal::Vector3D {0.1, 0.2, 0.4}
    : primal::Vector3D {params.direction.data()};
  // discreteFunction are discrete z-r pairs describing the function
  // to be rotated around the z axis.
  axom::Array<double, 2> discreteFunction({2, 2}, axom::ArrayStrideOrder::ROW);
  double baseRadius = params.radius < 0 ? 1.23 : params.radius;
  double topRadius = params.radius2 < 0 ? 0.176 : params.radius2;
  double height = params.length < 0 ? 2.3 : params.length;
  discreteFunction(0, 0) = -height / 2;
  discreteFunction(0, 1) = baseRadius;
  discreteFunction(1, 0) = height / 2;
  discreteFunction(1, 1) = topRadius;

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addTranslateOperator(*compositeOp, 1, 1, -1);

  axom::klee::Geometry sorGeometry =
    createGeometry_Sor(sorBase, sorDirection, discreteFunction, compositeOp);

  exactOverlapVols["cone"] = computeVolume_Sor(discreteFunction);
  errorToleranceRel["cone"] = 1e-4;
  errorToleranceAbs["cone"] = 1e-4;

  return sorGeometry;
}

axom::klee::Geometry createGeom_Tet()
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  // Tetrahedron at origin.
  const double len = params.length < 0 ? 1.5 : params.length;
  const Point3D a {-len, -len, -len};
  const Point3D b {+len, -len, -len};
  const Point3D c {+len, +len, -len};
  const Point3D d {-len, +len, +len};
  const primal::Tetrahedron<double, 3> tet {a, b, c, d};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp, -1, 1, -1);
  exactOverlapVols["tet"] = tet.volume();
  errorToleranceRel["tet"] = 1e-6;
  errorToleranceAbs["tet"] = 1e-8;

  axom::klee::Geometry tetGeometry(prop, tet, compositeOp);

  return tetGeometry;
}

axom::klee::Geometry createGeom_Hex()
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
  addTranslateOperator(*compositeOp, -1, -1, -1);
  exactOverlapVols["hex"] = hex.volume();
  errorToleranceRel["hex"] = 1e-6;
  errorToleranceAbs["hex"] = 1e-8;

  axom::klee::Geometry hexGeometry(prop, hex, compositeOp);

  return hexGeometry;
}

axom::klee::Geometry createGeom_Plane()
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
  exactOverlapVols["plane"] = 0.5 * meshVolume;
  errorToleranceRel["plane"] = 1e-6;
  errorToleranceAbs["plane"] = 1e-8;

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
    MPI_Finalize();
#endif
    exit(retval);
  }

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(params.annotationMode);

  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int allocId = axom::policyToDefaultAllocatorID(params.policy);
#if defined(AXOM_USE_UMPIRE)
  const std::string allocatorName =
    umpire::ResourceManager::getInstance().getAllocator(allocId).getName();
  SLIC_INFO(axom::fmt::format("Using allocator id {}, '{}'", allocId, allocatorName));
#endif

  AXOM_ANNOTATE_BEGIN("quest example for shaping primals");
  AXOM_ANNOTATE_BEGIN("init");

  // Storage for the some geometry meshes.
  sidre::DataStore ds;

  //---------------------------------------------------------------------------
  // Create shapes for the test
  //---------------------------------------------------------------------------
  axom::Array<std::shared_ptr<axom::quest::GeometryClipperStrategy>> geomStrategies;
  geomStrategies.reserve(params.testGeom.size());
  SLIC_ERROR_IF(params.getBoxDim() != 3, "This example is only in 3D.");
  for(const auto& tg : params.testGeom)
  {
    if(tg == "plane")
    {
      geomStrategies.push_back(std::make_shared<axom::quest::Plane3DClipper>(createGeom_Plane(), tg));
    }
    else if(tg == "hex")
    {
      geomStrategies.push_back(std::make_shared<axom::quest::HexClipper>(createGeom_Hex(), tg));
    }
    else if(tg == "sphere")
    {
      geomStrategies.push_back(std::make_shared<axom::quest::SphereClipper>(createGeom_Sphere(), tg));
    }
    else if(tg == "tetmesh")
    {
      geomStrategies.push_back(
        std::make_shared<axom::quest::TetMeshClipper>(createGeom_TetMesh(ds), tg));
    }
    else if(tg == "tet")
    {
      geomStrategies.push_back(std::make_shared<axom::quest::TetClipper>(createGeom_Tet(), tg));
    }
    else if(tg == "sor")
    {
      geomStrategies.push_back(std::make_shared<axom::quest::SorClipper>(createGeom_Sor(), tg));
    }
    else if(tg == "cyl")
    {
      geomStrategies.push_back(std::make_shared<axom::quest::SorClipper>(createGeom_Cylinder(), tg));
    }
    else if(tg == "cone")
    {
      geomStrategies.push_back(std::make_shared<axom::quest::SorClipper>(createGeom_Cone(), tg));
    }
  }

  {
    compMeshGrp = ds.getRoot()->createGroup("compMesh");
    compMeshGrp->setDefaultAllocator(allocId);

    createBoxMesh(compMeshGrp);

    SLIC_INFO(axom::fmt::format("{:-^80}", "Generated Blueprint mesh"));
    cellCount = params.getBoxCellCount();
  }

  //---------------------------------------------------------------------------
  // Initialize computational mesh.
  //---------------------------------------------------------------------------
  std::shared_ptr<quest::ShapeeMesh> sMeshPtr;
  AXOM_ANNOTATE_BEGIN("setup shaping problem");
  if(params.useBlueprintSidre())
  {
    sMeshPtr =
      std::make_shared<quest::ShapeeMesh>(params.policy, allocId, compMeshGrp, topoName, matsetName);
  }
  if(params.useBlueprintConduit())
  {
    compMeshNode.reset(new conduit::Node);
    compMeshNode->set_allocator(axom::ConduitMemory::axomAllocIdToConduit(allocId));
    compMeshGrp->createNativeLayout(*compMeshNode);
    compMeshNode->set_allocator(axom::ConduitMemory::axomAllocIdToConduit(allocId));

    sMeshPtr =
      std::make_shared<quest::ShapeeMesh>(params.policy, allocId, *compMeshNode, topoName, matsetName);
  }
  quest::ShapeeMesh& sMesh = *sMeshPtr;

  AXOM_ANNOTATE_END("setup shaping problem");
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

    quest::GeometryClipper clipper(sMesh, geomStrategies[i]);
    clipper.setVerbose(true);
    axom::Array<double> ovlap;
    AXOM_ANNOTATE_BEGIN(annotationName);
    clipper.clip(ovlap);
    AXOM_ANNOTATE_END(annotationName);

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
    double correctOverlapVol = exactOverlapVols[geomName];

    bool err = !axom::utilities::isNearlyEqualRelative(computedOverlapVol,
                                                       correctOverlapVol,
                                                       errorToleranceRel[geomName],
                                                       errorToleranceAbs[geomName]);
    failCounts += err;

    SLIC_INFO(axom::fmt::format("{:-^80}",
                                axom::fmt::format("Shape '{}' has volume {} vs {}, diff of {}, {}.",
                                                  geomName,
                                                  computedOverlapVol,
                                                  correctOverlapVol,
                                                  computedOverlapVol - correctOverlapVol,
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

  AXOM_ANNOTATE_END("quest example for shaping primals");

  finalizeLogger();

  return failCounts;
}
