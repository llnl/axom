// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_shape_in_memory.cpp
 * \brief Driver application for shaping material volume fractions onto a simulation mesh
 * using a mint::UnstructuredMesh of tets.
 * Modeled after shaping_driver.cpp.
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

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#if !defined(AXOM_USE_CONDUIT)
  #error Shaping functionality requires Axom to be configured with Conduit
#endif

#include "conduit_relay_io_blueprint.hpp"

#if defined(AXOM_USE_MFEM)
  #include "mfem.hpp"
#endif

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

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
  int getBoxCellCount() const
  {
    if(getBoxDim() == 3)
    {
      return boxResolution[0] * boxResolution[1] * boxResolution[2];
    }
    else
    {
      return boxResolution[0] * boxResolution[1];
    }
  }

  // The shape to run.
  std::string testShape {"tetmesh"};
  // The shapes this example is set up to run.
  const std::set<std::string>
    availableShapes {"tetmesh", "tri", "sphere", "cyl", "cone", "sor", "tet", "hex", "plane", "all"};

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

#if defined(AXOM_USE_MFEM)
  /// Generate an mfem Cartesian mesh, scaled to the bounding box range
  mfem::Mesh* createBoxMesh()
  {
    mfem::Mesh* mesh = nullptr;

    switch(getBoxDim())
    {
    case 2:
    {
      using BBox2D = primal::BoundingBox<double, 2>;
      using Pt2D = primal::Point<double, 2>;
      auto res = axom::NumericArray<int, 2>(boxResolution.data());
      auto bbox = BBox2D(Pt2D(boxMins.data()), Pt2D(boxMaxs.data()));

      SLIC_INFO(
        axom::fmt::format("Creating inline box mesh of resolution {} and bounding box {}", res, bbox));

      mesh = quest::util::make_cartesian_mfem_mesh_2D(bbox, res, outputOrder);
    }
    break;
    case 3:
    {
      using BBox3D = primal::BoundingBox<double, 3>;
      using Pt3D = primal::Point<double, 3>;
      auto res = axom::NumericArray<int, 3>(boxResolution.data());
      auto bbox = BBox3D(Pt3D(boxMins.data()), Pt3D(boxMaxs.data()));

      SLIC_INFO(
        axom::fmt::format("Creating inline box mesh of resolution {} and bounding box {}", res, bbox));

      mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, res, outputOrder);
    }
    break;
    default:
      SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
      break;
    }

    // Handle conversion to parallel mfem mesh
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    {
      int* partitioning = nullptr;
      int part_method = 0;
      mfem::Mesh* parallelMesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh, partitioning, part_method);
      delete[] partitioning;
      delete mesh;
      mesh = parallelMesh;
    }
  #endif

    return mesh;
  }

  std::unique_ptr<sidre::MFEMSidreDataCollection> loadComputationalMesh()
  {
    constexpr bool dc_owns_data = true;
    mfem::Mesh* mesh = createBoxMesh();
    std::string name = "mesh";

    auto dc = std::unique_ptr<sidre::MFEMSidreDataCollection>(
      new sidre::MFEMSidreDataCollection(name, mesh, dc_owns_data));
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    dc->SetComm(MPI_COMM_WORLD);
  #endif

    return dc;
  }
#endif

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

    app.add_option("-s,--testShape", testShape)
      ->description("The shape to run")
      ->check(axom::CLI::IsMember(availableShapes));

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

#if defined(AXOM_USE_MFEM)
/**
 * \brief Print some info about the mesh
 *
 * \note In MPI-based configurations, this is a collective call, but
 * only prints on rank 0
 */
void printMfemMeshInfo(mfem::Mesh* mesh, const std::string& prefixMessage = "")
{
  namespace primal = axom::primal;

  int myRank = 0;
  #ifdef AXOM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  #endif

  int numElements = mesh->GetNE();

  mfem::Vector mins, maxs;
  #ifdef MFEM_USE_MPI
  auto* parallelMesh = dynamic_cast<mfem::ParMesh*>(mesh);
  if(parallelMesh != nullptr)
  {
    parallelMesh->GetBoundingBox(mins, maxs);
    numElements = parallelMesh->ReduceInt(numElements);
    myRank = parallelMesh->GetMyRank();
  }
  else
  #endif
  {
    mesh->GetBoundingBox(mins, maxs);
  }

  if(myRank == 0)
  {
    switch(mesh->Dimension())
    {
    case 2:
      SLIC_INFO(axom::fmt::format(
        axom::utilities::locale(),
        "{} mesh has {:L} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 2>(primal::Point<double, 2>(mins.GetData()),
                                       primal::Point<double, 2>(maxs.GetData()))));
      break;
    case 3:
      SLIC_INFO(axom::fmt::format(
        axom::utilities::locale(),
        "{} mesh has {:L} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 3>(primal::Point<double, 3>(mins.GetData()),
                                       primal::Point<double, 3>(maxs.GetData()))));
      break;
    }
  }

  slic::flushStreams();
}
#endif

/************************************************************
 * Shared variables.
 ************************************************************/

const std::string topoName = "mesh";
const std::string coordsetName = "coords";
int cellCount = -1;

const auto hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
int arrayAllocId = axom::INVALID_ALLOCATOR_ID;

// Computational mesh in different forms, initialized in main
#if defined(AXOM_USE_MFEM)
std::shared_ptr<sidre::MFEMSidreDataCollection> shapingDC;
#endif
axom::sidre::Group* compMeshGrp = nullptr;
std::shared_ptr<conduit::Node> compMeshNode;

axom::sidre::Group* createBoxMesh(axom::sidre::Group* meshGrp)
{
  switch(params.getBoxDim())
  {
  case 2:
  {
    using BBox2D = primal::BoundingBox<double, 2>;
    using Pt2D = primal::Point<double, 2>;
    auto res = axom::NumericArray<int, 2>(params.boxResolution.data());
    auto bbox = BBox2D(Pt2D(params.boxMins.data()), Pt2D(params.boxMaxs.data()));

    SLIC_INFO(
      axom::fmt::format("Creating inline box mesh of resolution {} and bounding box {}", res, bbox));

    axom::quest::util::make_unstructured_blueprint_box_mesh_2d(meshGrp,
                                                               bbox,
                                                               res,
                                                               topoName,
                                                               coordsetName,
                                                               params.policy);
  }
  break;
  case 3:
  {
    using BBox3D = primal::BoundingBox<double, 3>;
    using Pt3D = primal::Point<double, 3>;
    auto res = axom::NumericArray<int, 3>(params.boxResolution.data());
    auto bbox = BBox3D(Pt3D(params.boxMins.data()), Pt3D(params.boxMaxs.data()));

    SLIC_INFO(
      axom::fmt::format("Creating inline box mesh of resolution {} and bounding box {}", res, bbox));

    axom::quest::util::make_unstructured_blueprint_box_mesh_3d(meshGrp,
                                                               bbox,
                                                               res,
                                                               topoName,
                                                               coordsetName,
                                                               params.policy);
  }
  break;
  default:
    SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
    break;
  }

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

// Single triangle ShapeSet.
axom::klee::ShapeSet create2DShapeSet(sidre::DataStore& ds)
{
  sidre::Group* meshGroup = ds.getRoot()->createGroup("triangleMesh");
  const std::string topo = "mesh";
  const std::string coordset = "coords";
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> triangleMesh(2,
                                                                      axom::mint::CellType::TRIANGLE,
                                                                      meshGroup,
                                                                      topo,
                                                                      coordset);

  const double lll = 2.0;

  // Insert tet at origin.
  triangleMesh.appendNode(0.0, 0.0);
  triangleMesh.appendNode(lll, 0.0);
  triangleMesh.appendNode(0.0, lll);
  axom::IndexType conn[3] = {0, 1, 2};
  triangleMesh.appendCell(conn);

  SLIC_ERROR_ROOT_IF(!axom::mint::blueprint::isValidRootGroup(meshGroup),
                     "Triangle mesh blueprint is not valid");

  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Two,
                                                    axom::klee::LengthUnit::unspecified};
  axom::klee::Geometry triangleGeom(prop, triangleMesh.getSidreGroup(), topo, nullptr);

  std::vector<axom::klee::Shape> shapes;
  axom::klee::Shape triangleShape(
    "triangle",
    "AL",
    {},
    {},
    axom::klee::Geometry {prop, triangleMesh.getSidreGroup(), topo, nullptr});
  shapes.push_back(
    axom::klee::Shape {"triangle",
                       "AL",
                       {},
                       {},
                       axom::klee::Geometry {prop, triangleMesh.getSidreGroup(), topo, nullptr}});

  axom::klee::ShapeSet shapeSet;
  shapeSet.setShapes(shapes);
  shapeSet.setDimensions(axom::klee::Dimensions::Two);

  return shapeSet;
}

axom::klee::Shape createShape_Sphere()
{
  Point3D center = params.center.empty() ? Point3D {0, 0, 0} : Point3D {params.center.data()};
  double radius = params.radius < 0 ? 0.6 : params.radius;
  axom::primal::Sphere<double, 3> sphere {center, radius};

  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp, 1, 1, 1);

  const axom::IndexType levelOfRefinement = params.refinementLevel;
  axom::klee::Geometry sphereGeometry(prop, sphere, levelOfRefinement, compositeOp);
  axom::klee::Shape sphereShape("sphere", "SPHERE", {}, {}, sphereGeometry);

  return sphereShape;
}

axom::klee::Shape createShape_TetMesh(sidre::DataStore& ds)
{
  // Shape a tetrahedal mesh.
  sidre::Group* meshGroup = ds.getRoot()->createGroup("tetMesh");
  const std::string topo = "mesh";
  const std::string coordset = "coords";
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tetMesh(3,
                                                                 axom::mint::CellType::TET,
                                                                 meshGroup,
                                                                 topo,
                                                                 coordset);

  double lll = params.length < 0 ? 0.7 : params.length;

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

  SLIC_ERROR_ROOT_IF(!axom::mint::blueprint::isValidRootGroup(meshGroup),
                     "Tet mesh blueprint is not valid");

  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp, -1, 1, 1);

  axom::klee::Geometry tetMeshGeometry(prop, tetMesh.getSidreGroup(), topo, compositeOp);
  axom::klee::Shape tetShape("tetmesh", "TETMESH", {}, {}, tetMeshGeometry);

  return tetShape;
}

axom::klee::Geometry createGeometry_Sor(axom::primal::Point<double, 3>& sorBase,
                                        axom::primal::Vector<double, 3>& sorDirection,
                                        axom::Array<double, 2>& discreteFunction,
                                        std::shared_ptr<axom::klee::CompositeOperator>& compositeOp)
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);

  const axom::IndexType levelOfRefinement = params.refinementLevel;
  axom::klee::Geometry sorGeometry(prop,
                                   discreteFunction,
                                   sorBase,
                                   sorDirection,
                                   levelOfRefinement,
                                   compositeOp);
  return sorGeometry;
}

axom::klee::Shape createShape_Sor()
{
  Point3D sorBase = params.center.empty() ? Point3D {0.0, 0.0, 0.0} : Point3D {params.center.data()};
  axom::primal::Vector<double, 3> sorDirection = params.direction.empty()
    ? primal::Vector3D {0.1, 0.2, 0.4}
    : primal::Vector3D {params.direction.data()};
  const int numIntervals = 5;
  // discreteFunction are discrete z-r pairs describing the function
  // to be rotated around the z axis.
  axom::Array<double, 2> discreteFunction({numIntervals + 1, 2}, axom::ArrayStrideOrder::ROW);
  double zLen = params.length < 0 ? 1.6 : params.length;
  double zShift = -zLen / 2;
  double maxR = params.radius < 0 ? 0.75 : params.radius;
  double dz = zLen / numIntervals;
  discreteFunction[0][0] = 0 * dz + zShift;
  discreteFunction[0][1] = 0.0 * maxR;
  discreteFunction[1][0] = 1 * dz + zShift;
  discreteFunction[1][1] = 0.8 * maxR;
  discreteFunction[2][0] = 2 * dz + zShift;
  discreteFunction[2][1] = 0.4 * maxR;
  discreteFunction[3][0] = 3 * dz + zShift;
  discreteFunction[3][1] = 0.5 * maxR;
  discreteFunction[4][0] = 4 * dz + zShift;
  discreteFunction[4][1] = 1.0 * maxR;
  discreteFunction[5][0] = 5 * dz + zShift;
  discreteFunction[5][1] = 0.0;

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addTranslateOperator(*compositeOp, -1, -1, 1);

  axom::klee::Geometry sorGeometry =
    createGeometry_Sor(sorBase, sorDirection, discreteFunction, compositeOp);

  axom::klee::Shape sorShape("sor", "SOR", {}, {}, sorGeometry);

  return sorShape;
}

axom::klee::Shape createShape_Cylinder()
{
  Point3D sorBase = params.center.empty() ? Point3D {0.0, 0.0, 0.0} : Point3D {params.center.data()};
  axom::primal::Vector<double, 3> sorDirection = params.direction.empty()
    ? primal::Vector3D {0.1, 0.2, 0.4}
    : primal::Vector3D {params.direction.data()};
  // discreteFunction are discrete z-r pairs describing the function
  // to be rotated around the z axis.
  axom::Array<double, 2> discreteFunction({2, 2}, axom::ArrayStrideOrder::ROW);
  double radius = params.radius < 0 ? 0.5 : params.radius;
  double height = params.length < 0 ? 1.2 : params.length;
  discreteFunction[0][0] = -height / 2;
  discreteFunction[0][1] = radius;
  discreteFunction[1][0] = height / 2;
  discreteFunction[1][1] = radius;

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addTranslateOperator(*compositeOp, 1, -1, 1);

  axom::klee::Geometry sorGeometry =
    createGeometry_Sor(sorBase, sorDirection, discreteFunction, compositeOp);

  axom::klee::Shape sorShape("cyl", "CYL", {}, {}, sorGeometry);

  return sorShape;
}

axom::klee::Shape createShape_Cone()
{
  Point3D sorBase = params.center.empty() ? Point3D {0.0, 0.0, 0.0} : Point3D {params.center.data()};
  axom::primal::Vector<double, 3> sorDirection = params.direction.empty()
    ? primal::Vector3D {0.1, 0.2, 0.4}
    : primal::Vector3D {params.direction.data()};
  // discreteFunction are discrete z-r pairs describing the function
  // to be rotated around the z axis.
  axom::Array<double, 2> discreteFunction({2, 2}, axom::ArrayStrideOrder::ROW);
  double baseRadius = params.radius < 0 ? 0.7 : params.radius;
  double topRadius = params.radius2 < 0 ? 0.1 : params.radius2;
  double height = params.length < 0 ? 1.3 : params.length;
  discreteFunction[0][0] = -height / 2;
  discreteFunction[0][1] = baseRadius;
  discreteFunction[1][0] = height / 2;
  discreteFunction[1][1] = topRadius;

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addTranslateOperator(*compositeOp, 1, 1, -1);

  axom::klee::Geometry sorGeometry =
    createGeometry_Sor(sorBase, sorDirection, discreteFunction, compositeOp);

  axom::klee::Shape sorShape("cone", "CONE", {}, {}, sorGeometry);

  return sorShape;
}

axom::klee::Shape createShape_Tet()
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);

  // Tetrahedron at origin.
  const double len = params.length < 0 ? 0.8 : params.length;
  const Point3D a {-len, -len, -len};
  const Point3D b {+len, -len, -len};
  const Point3D c {+len, +len, -len};
  const Point3D d {-len, +len, +len};
  const primal::Tetrahedron<double, 3> tet {a, b, c, d};

  auto compositeOp = std::make_shared<axom::klee::CompositeOperator>(startProp);
  addScaleOperator(*compositeOp);
  addRotateOperator(*compositeOp);
  addTranslateOperator(*compositeOp, -1, 1, -1);

  axom::klee::Geometry tetGeometry(prop, tet, compositeOp);
  axom::klee::Shape tetShape("tet", "TET", {}, {}, tetGeometry);

  return tetShape;
}

axom::klee::Shape createShape_Hex()
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);

  const double md = params.length < 0 ? 0.6 : params.length / 2;
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

  axom::klee::Geometry hexGeometry(prop, hex, compositeOp);
  axom::klee::Shape hexShape("hex", "HEX", {}, {}, hexGeometry);

  return hexShape;
}

axom::klee::Shape createShape_Plane()
{
  axom::klee::TransformableGeometryProperties prop {axom::klee::Dimensions::Three,
                                                    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  std::shared_ptr<axom::klee::Scale> scaleOp;
  if(!params.scaleFactors.empty())
  {
    scaleOp = std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                                  params.scaleFactors[1],
                                                  params.scaleFactors[2],
                                                  prop);
  }

  // Create a plane crossing center of mesh.  No matter the normal, it cuts the mesh in half.
  Point3D center {0.5 *
                  (axom::NumericArray<double, 3>(params.boxMins.data()) +
                   axom::NumericArray<double, 3>(params.boxMaxs.data()))};
  primal::Vector<double, 3> normal = params.direction.empty()
    ? primal::Vector3D {1.0, 0.0, 0.0}
    : primal::Vector3D {params.direction.data()}.unitVector();
  const primal::Plane<double, 3> plane {normal, center, true};

  axom::klee::Geometry planeGeometry(prop, plane, scaleOp);
  axom::klee::Shape planeShape("plane", "PLANE", {}, {}, planeGeometry);

  return planeShape;
}

//!@brief Create a ShapeSet with a single shape.
axom::klee::ShapeSet createShapeSet(const axom::klee::Shape& shape)
{
  axom::klee::ShapeSet shapeSet;
  shapeSet.setShapes(std::vector<axom::klee::Shape> {shape});
  shapeSet.setDimensions(axom::klee::Dimensions::Three);

  return shapeSet;
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

double areaOfTriMesh(const axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& triMesh)
{
  using PolygonStaticType = primal::Polygon<double, 2, axom::primal::PolygonArray::Static>;
  using Point2D = primal::Point<double, 2>;
  axom::StackArray<axom::IndexType, 1> nodesShape {triMesh.getNumberOfNodes()};
  axom::ArrayView<const double> x(triMesh.getCoordinateArray(0), nodesShape);
  axom::ArrayView<const double> y(triMesh.getCoordinateArray(1), nodesShape);
  const axom::IndexType triCount = triMesh.getNumberOfCells();
  axom::Array<double> triAreas(triCount, triCount);
  double meshArea = 0.0;
  for(axom::IndexType ic = 0; ic < triCount; ++ic)
  {
    const axom::IndexType* nodeIds = triMesh.getCellNodeIDs(ic);
    PolygonStaticType tri;
    for(int j = 0; j < 3; ++j)
    {
      auto cornerNodeId = nodeIds[j];
      tri.addVertex(Point2D({x[cornerNodeId], y[cornerNodeId]}));
    }
    meshArea += tri.area();
  }
  return meshArea;
}

#if defined(AXOM_USE_MFEM)
/*!
 *  @brief Return the element volumes as a sidre::View.
 *
 *  If it doesn't exist, allocate and compute it.
 *  \post The volume data is in \c dc->GetNamedBuffer(volFieldName).
 *
 *  Most of this is lifted from IntersectionShaper::runShapeQueryImpl.
 */
template <typename ExecSpace>
axom::sidre::View* getElementVolumes(
  sidre::MFEMSidreDataCollection* dc,
  const std::string& volFieldName = std::string("elementVolumes"))
{
  axom::sidre::View* volSidreView = dc->GetNamedBuffer(volFieldName);

  switch(params.getBoxDim())
  {
  case 2:
  {
    using PolygonStaticType = primal::Polygon<double, 2, axom::primal::PolygonArray::Static>;
    using Point2D = primal::Point<double, 2>;

    if(volSidreView == nullptr)
    {
      mfem::Mesh* mesh = dc->GetMesh();

      constexpr int NUM_VERTS_PER_QUAD = 4;
      constexpr int NUM_COMPS_PER_VERT = 2;

      axom::Array<Point2D> vertCoords(cellCount * NUM_VERTS_PER_QUAD, cellCount * NUM_VERTS_PER_QUAD);
      auto vertCoordsView = vertCoords.view();

      // This runs only only on host, because the mfem::Mesh only uses host memory, I think.
      for(axom::IndexType cellIdx = 0; cellIdx < cellCount; ++cellIdx)
      {
        // Get the indices of this element's vertices
        mfem::Array<int> verts;
        mesh->GetElementVertices(cellIdx, verts);
        SLIC_ASSERT(verts.Size() == NUM_VERTS_PER_QUAD);

        // Get the coordinates for the vertices
        for(int j = 0; j < NUM_VERTS_PER_QUAD; ++j)
        {
          int vertIdx = cellIdx * NUM_VERTS_PER_QUAD + j;
          for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
          {
            vertCoordsView[vertIdx][k] = (mesh->GetVertex(verts[j]))[k];
          }
        }
      }

      // Initialize quad elements.
      axom::Array<PolygonStaticType> quads(cellCount, cellCount);
      auto quadsView = quads.view();
      axom::for_all<ExecSpace>(
        cellCount,
        AXOM_LAMBDA(axom::IndexType cellIdx) {
          // Set each quad element vertices
          quadsView[cellIdx] = PolygonStaticType();
          for(int j = 0; j < NUM_VERTS_PER_QUAD; ++j)
          {
            int vertIndex = (cellIdx * NUM_VERTS_PER_QUAD) + j;
            auto& quad = quadsView[cellIdx];
            quad.addVertex(vertCoordsView[vertIndex]);
          }
        });  // end of loop to initialize quad elements

      // Allocate and populate cell volumes.
      volSidreView = dc->AllocNamedBuffer(volFieldName, cellCount);
      axom::ArrayView<double> volView(volSidreView->getData(), volSidreView->getNumElements());
      axom::for_all<ExecSpace>(
        cellCount,
        AXOM_LAMBDA(axom::IndexType cellIdx) { volView[cellIdx] = quadsView[cellIdx].area(); });
    }
  }  // end of case 2
  break;
  case 3:
  {
    using HexahedronType = axom::primal::Hexahedron<double, 3>;

    if(volSidreView == nullptr)
    {
      mfem::Mesh* mesh = dc->GetMesh();

      constexpr int NUM_VERTS_PER_HEX = 8;
      constexpr int NUM_COMPS_PER_VERT = 3;

      axom::Array<Point3D> vertCoords(cellCount * NUM_VERTS_PER_HEX, cellCount * NUM_VERTS_PER_HEX);
      auto vertCoordsView = vertCoords.view();

      // This runs only only on host, because the mfem::Mesh only uses host memory, I think.
      for(axom::IndexType cellIdx = 0; cellIdx < cellCount; ++cellIdx)
      {
        // Get the indices of this element's vertices
        mfem::Array<int> verts;
        mesh->GetElementVertices(cellIdx, verts);
        SLIC_ASSERT(verts.Size() == NUM_VERTS_PER_HEX);

        // Get the coordinates for the vertices
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIdx = cellIdx * NUM_VERTS_PER_HEX + j;
          for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
          {
            vertCoordsView[vertIdx][k] = (mesh->GetVertex(verts[j]))[k];
          }
        }
      }

      // Initialize hexahedral elements.
      axom::Array<HexahedronType> hexes(cellCount, cellCount);
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
      volSidreView = dc->AllocNamedBuffer(volFieldName, cellCount);
      axom::ArrayView<double> volView(volSidreView->getData(), volSidreView->getNumElements());
      axom::for_all<ExecSpace>(
        cellCount,
        AXOM_LAMBDA(axom::IndexType cellIdx) { volView[cellIdx] = hexesView[cellIdx].volume(); });
    }
  }  // end of case 3
  break;
  default:
    SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
    break;
  }

  return volSidreView;
}
#endif

/*!
 *  @brief Return the element volumes as a sidre::View containing the volumes in an array
 *
 *  If it doesn't exist, allocate and compute it.
 *  \post The volume data is in the blueprint field \c volFieldName.
 *
 *  Most of this is lifted from IntersectionShaper::runShapeQueryImpl.
 */
template <typename ExecSpace>
axom::sidre::View* getElementVolumes(
  sidre::Group* meshGrp,
  const std::string& volFieldName = std::string("elementVolumes"))
{
  using XS = axom::execution_space<ExecSpace>;
  axom::sidre::View* volSidreView = nullptr;

  switch(params.getBoxDim())
  {
  case 2:
  {
    using PolygonStaticType = primal::Polygon<double, 2, axom::primal::PolygonArray::Static>;
    using Point2D = primal::Point<double, 2>;

    const auto fieldPath = axom::fmt::format("fields/{}", volFieldName);
    if(meshGrp->hasGroup(fieldPath))
    {
      sidre::Group* fieldGrp = meshGrp->getGroup(fieldPath);
      volSidreView = fieldGrp->getView("values");
    }
    else
    {
      axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(meshGrp, topoName);

      constexpr int NUM_VERTS_PER_QUAD = 4;
      constexpr int NUM_COMPS_PER_VERT = 2;

      /*
        Get vertex coordinates.  We use UnstructuredMesh for this,
        so get it on host first then transfer to device if needed.
      */
      auto* connData = meshGrp->getGroup("topologies")
                         ->getGroup(topoName)
                         ->getGroup("elements")
                         ->getView("connectivity");
      SLIC_ASSERT(connData->getNode().dtype().id() == axom::quest::conduitDataIdOfAxomIndexType);

      conduit::Node coordNode;
      meshGrp->getGroup("coordsets")->getGroup(coordsetName)->createNativeLayout(coordNode);
      const conduit::Node& coordValues = coordNode.fetch_existing("values");
      axom::IndexType vertexCount = coordValues["x"].dtype().number_of_elements();
      bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
      int stride = isInterleaved ? NUM_COMPS_PER_VERT : 1;
      axom::StackArray<axom::ArrayView<const double>, 2> coordArrays {
        axom::ArrayView<const double>(coordValues["x"].as_double_ptr(), {vertexCount}, stride),
        axom::ArrayView<const double>(coordValues["y"].as_double_ptr(), {vertexCount}, stride)};

      const axom::IndexType* connPtr = connData->getArray();
      SLIC_ASSERT(connPtr != nullptr);
      axom::ArrayView<const axom::IndexType, 2> conn(connPtr, cellCount, NUM_VERTS_PER_QUAD);
      axom::Array<Point2D> vertCoords(cellCount * NUM_VERTS_PER_QUAD,
                                      cellCount * NUM_VERTS_PER_QUAD,
                                      XS::allocatorID());
      auto vertCoordsView = vertCoords.view();

      axom::for_all<ExecSpace>(
        cellCount,
        AXOM_LAMBDA(axom::IndexType cellIdx) {
          // Get the indices of this element's vertices
          auto verts = conn[cellIdx];

          // Get the coordinates for the vertices
          for(int j = 0; j < NUM_VERTS_PER_QUAD; ++j)
          {
            int vertIdx = cellIdx * NUM_VERTS_PER_QUAD + j;
            for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
            {
              vertCoordsView[vertIdx][k] = coordArrays[k][verts[j]];
              // vertCoordsView[vertIdx][k] = mesh.getNodeCoordinate(verts[j], k);
            }
          }
        });

      // Initialize quad elements.
      axom::Array<PolygonStaticType> quads(cellCount, cellCount, meshGrp->getDefaultAllocatorID());
      auto quadsView = quads.view();
      axom::for_all<ExecSpace>(
        cellCount,
        AXOM_LAMBDA(axom::IndexType cellIdx) {
          // Set each quad element vertice
          quadsView[cellIdx] = PolygonStaticType();
          for(int j = 0; j < NUM_VERTS_PER_QUAD; ++j)
          {
            int vertIndex = (cellIdx * NUM_VERTS_PER_QUAD) + j;
            auto& quad = quadsView[cellIdx];
            quad.addVertex(vertCoordsView[vertIndex]);
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
        AXOM_LAMBDA(axom::IndexType cellIdx) { volView[cellIdx] = quadsView[cellIdx].area(); });
    }
  }
  break;
  case 3:
  {
    using HexahedronType = axom::primal::Hexahedron<double, 3>;
    using Point3D = primal::Point<double, 3>;

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

      /*
          Get vertex coordinates.  We use UnstructuredMesh for this,
          so get it on host first then transfer to device if needed.
        */
      auto* connData = meshGrp->getGroup("topologies")
                         ->getGroup(topoName)
                         ->getGroup("elements")
                         ->getView("connectivity");
      SLIC_ASSERT(connData->getNode().dtype().id() == axom::quest::conduitDataIdOfAxomIndexType);

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
  }
  break;
  default:
    SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
    break;
  }

  return volSidreView;
}

#if defined(AXOM_USE_MFEM)
/*!
  @brief Return global sum of volume of the given material.
*/
template <typename ExecSpace>
double sumMaterialVolumes(sidre::MFEMSidreDataCollection* dc, const std::string& material)
{
  // Get cell volumes from dc.
  axom::sidre::View* elementVols = getElementVolumes<ExecSpace>(dc);
  axom::ArrayView<double> elementVolsView(elementVols->getData(), elementVols->getNumElements());

  // Get material volume fractions
  const std::string materialFieldName = axom::fmt::format("vol_frac_{}", material);
  mfem::GridFunction* volFracGf = dc->GetField(materialFieldName);
  axom::ArrayView<double> volFracGfArrayView(volFracGf->GetData(), volFracGf->Size());
  axom::Array<double> volFracGfArray(volFracGfArrayView,
                                     axom::execution_space<ExecSpace>::allocatorID());
  auto volFracView = volFracGfArray.view();

  axom::ReduceSum<ExecSpace, double> localVol(0);
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType i) { localVol += volFracView[i] * elementVolsView[i]; });

  double globalVol = localVol.get();
  #ifdef AXOM_USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &globalVol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  return globalVol;
}
#endif

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

  axom::ReduceSum<ExecSpace, double> localVol(0);
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
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  if(params.policy == RuntimePolicy::omp)
  {
    rval = sumMaterialVolumesImpl<axom::OMP_EXEC>(meshGrp, material);
  }
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  if(params.policy == RuntimePolicy::cuda)
  {
    rval = sumMaterialVolumesImpl<axom::CUDA_EXEC<256>>(meshGrp, material);
  }
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_HIP)
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

axom::ArrayView<double> getFieldAsArrayView(const std::string& fieldName)
{
  axom::ArrayView<double> fieldDataView;
  if(params.useBlueprintSidre())
  {
    std::string valuesPath = "fields/" + fieldName + "/values";
    axom::sidre::View* fieldValues = compMeshGrp->getView(valuesPath);
    double* fieldData = fieldValues->getArray();
    fieldDataView = axom::ArrayView<double>(fieldData, fieldValues->getNumElements());
  }
  if(params.useBlueprintConduit())
  {
    std::string valuesPath = "fields/" + fieldName + "/values";
    conduit::Node& fieldValues = compMeshNode->fetch_existing(valuesPath);
    double* fieldData = fieldValues.as_double_ptr();
    fieldDataView = axom::ArrayView<double>(fieldData, fieldValues.dtype().number_of_elements());
  }
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    // auto* mfemMesh = shapingDC->GetMesh();
    mfem::GridFunction* gridFunc = shapingDC->GetField(fieldName);
    fieldDataView = axom::ArrayView<double>(gridFunc->GetData(), gridFunc->Size());
  }
#endif
  return fieldDataView;
}

//!@brief Fill a sidre array View with a value.
// No error checking.
template <typename T>
void fillSidreViewData(axom::sidre::View* view, const T& value)
{
  double* valuesPtr = view->getData<T*>();
  switch(params.policy)
  {
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  case RuntimePolicy::cuda:
    axom::for_all<axom::CUDA_EXEC<256>>(
      view->getNumElements(),
      AXOM_LAMBDA(axom::IndexType i) { valuesPtr[i] = value; });
    break;
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_HIP)
  case RuntimePolicy::hip:
    axom::for_all<axom::HIP_EXEC<256>>(
      view->getNumElements(),
      AXOM_LAMBDA(axom::IndexType i) { valuesPtr[i] = value; });
    break;
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OMP)
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
    exit(retval);
  }

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(params.annotationMode);

  const int arrayAllocId = axom::policyToDefaultAllocatorID(params.policy);

  AXOM_ANNOTATE_SCOPE("quest shaping example");
  AXOM_ANNOTATE_BEGIN("init");

  // Storage for the shape geometry meshes.
  sidre::DataStore ds;

  //---------------------------------------------------------------------------
  // Create simple ShapeSet for the example.
  //---------------------------------------------------------------------------
  axom::klee::ShapeSet shapeSet;

  if(params.testShape == "tetmesh")
  {
    shapeSet = createShapeSet(createShape_TetMesh(ds));
  }
  else if(params.testShape == "tet")
  {
    shapeSet = createShapeSet(createShape_Tet());
  }
  else if(params.testShape == "tri")
  {
    SLIC_ERROR_IF(params.getBoxDim() != 2, "This example is only in 2D.");
    shapeSet = create2DShapeSet(ds);
  }
  else if(params.testShape == "hex")
  {
    shapeSet = createShapeSet(createShape_Hex());
  }
  else if(params.testShape == "sphere")
  {
    shapeSet = createShapeSet(createShape_Sphere());
  }
  else if(params.testShape == "cyl")
  {
    shapeSet = createShapeSet(createShape_Cylinder());
  }
  else if(params.testShape == "cone")
  {
    shapeSet = createShapeSet(createShape_Cone());
  }
  else if(params.testShape == "sor")
  {
    shapeSet = createShapeSet(createShape_Sor());
  }
  else if(params.testShape == "plane")
  {
    shapeSet = createShapeSet(createShape_Plane());
  }
  else if(params.testShape == "all")
  {
    std::vector<axom::klee::Shape> shapesVec;
    shapesVec.push_back(createShape_TetMesh(ds));
    shapesVec.push_back(createShape_Tet());
    shapesVec.push_back(createShape_Hex());
    shapesVec.push_back(createShape_Sphere());
    shapesVec.push_back(createShape_Sor());
    shapesVec.push_back(createShape_Cylinder());
    shapesVec.push_back(createShape_Cone());
    shapeSet.setShapes(shapesVec);
    shapeSet.setDimensions(axom::klee::Dimensions::Three);
  }

  // Save the discrete shapes for viz and testing.
  auto* shapeMeshGroup = ds.getRoot()->createGroup("shapeMeshGroup");
  std::vector<std::shared_ptr<axom::mint::Mesh>> discreteShapeMeshes;
  for(const auto& shape : shapeSet.getShapes())
  {
    axom::quest::DiscreteShape dShape(shape, shapeMeshGroup);
    auto dMesh = std::dynamic_pointer_cast<axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>>(
      dShape.createMeshRepresentation());
    SLIC_INFO(axom::fmt::format("{:-^80}",
                                axom::fmt::format("Shape '{}' discrete geometry has {} cells",
                                                  shape.getName(),
                                                  dMesh->getNumberOfCells())));

    discreteShapeMeshes.push_back(dMesh);

    if(!params.outputFile.empty())
    {
      std::string shapeFileName = params.outputFile + ".shape";
      conduit::Node tmpNode, info;
      dMesh->getSidreGroup()->createNativeLayout(tmpNode);
      conduit::relay::io::blueprint::save_mesh(tmpNode, shapeFileName, "hdf5");
    }
  }

#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    AXOM_ANNOTATE_BEGIN("load mesh");
    //---------------------------------------------------------------------------
    // Load the computational mesh
    // originalMeshDC is the input MFEM mesh
    // It's converted to shapingDC below, and that is used for shaping calls.
    //---------------------------------------------------------------------------
    auto originalMeshDC = params.loadComputationalMesh();

    //---------------------------------------------------------------------------
    // Set up DataCollection for shaping
    // shapingDC is a "copy" of originalMeshDC, the MFEM mesh.
    // It's created empty then populated by the SetMesh call.
    // shapingMesh and parallelMesh are some kind of temporary versions of originalMeshDC.
    //---------------------------------------------------------------------------
    mfem::Mesh* shapingMesh = nullptr;
    constexpr bool dc_owns_data = true;
    shapingDC =
      std::make_shared<sidre::MFEMSidreDataCollection>("shaping", shapingMesh, dc_owns_data);
    {
      shapingDC->SetMeshNodesName("positions");

  // With MPI, loadComputationalMesh returns a parallel mesh.
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      mfem::ParMesh* parallelMesh = dynamic_cast<mfem::ParMesh*>(originalMeshDC->GetMesh());
      shapingMesh = (parallelMesh != nullptr) ? new mfem::ParMesh(*parallelMesh)
                                              : new mfem::Mesh(*originalMeshDC->GetMesh());
  #else
      shapingMesh = new mfem::Mesh(*originalMeshDC->GetMesh());
  #endif
      shapingDC->SetMesh(shapingMesh);
    }
    AXOM_ANNOTATE_END("load mesh");
    printMfemMeshInfo(shapingDC->GetMesh(), "Loaded MFEM mesh");

    cellCount = shapingMesh->GetNE();
  }
#else
  SLIC_ERROR_IF(params.useMfem(),
                "Cannot use MFEM mesh due to Axom configuration.  Please use "
                "Blueprint mesh or configure with MFEM and "
                "-DAXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION.");
#endif

  if(params.useBlueprintSidre() || params.useBlueprintConduit())
  {
    compMeshGrp = ds.getRoot()->createGroup("compMesh");
    compMeshGrp->setDefaultArrayAllocator(arrayAllocId);

    createBoxMesh(compMeshGrp);

    if(params.useBlueprintConduit())
    {
      // Intersection requires conduit mesh to have array data pre-allocated.
      auto makeField = [&](const std::string& fieldName, double initValue) {
        auto fieldGrp = compMeshGrp->createGroup("fields/" + fieldName);
        axom::IndexType shape[] = {params.getBoxCellCount(), 1};
        fieldGrp->createViewString("association", "element");
        fieldGrp->createViewString("topology", topoName);
        fieldGrp->createViewString("volume_dependent", "true");
        axom::sidre::View* valuesView =
          fieldGrp->createViewWithShapeAndAllocate("values",
                                                   axom::sidre::detail::SidreTT<double>::id,
                                                   2,
                                                   shape);
        fillSidreViewData(valuesView, initValue);
      };
      makeField("vol_frac_free", 1.0);
      const auto& shapes = shapeSet.getShapes();
      for(const auto& shape : shapes)
      {
        makeField("vol_frac_" + shape.getMaterial(),
                  0.0);  // Used in volume fraction computation
        makeField("shape_vol_frac_" + shape.getName(),
                  0.0);  // Used in applyReplacementRules
      }
    }

    /*
      Shallow-copy compMeshGrp into compMeshNode,
      so that any change in one is reflected in the other.
    */
    compMeshNode = std::make_shared<conduit::Node>();
    compMeshGrp->createNativeLayout(*compMeshNode);
    SLIC_INFO(axom::fmt::format("{:-^80}", "Generated Blueprint mesh"));
    cellCount = params.getBoxCellCount();
  }

  //---------------------------------------------------------------------------
  // Initialize the shaping query object
  //---------------------------------------------------------------------------
  AXOM_ANNOTATE_BEGIN("setup shaping problem");
  std::shared_ptr<quest::IntersectionShaper> shaper = nullptr;
  if(params.useBlueprintSidre())
  {
    shaper =
      std::make_shared<quest::IntersectionShaper>(params.policy, arrayAllocId, shapeSet, compMeshGrp);
  }
  if(params.useBlueprintConduit())
  {
    shaper =
      std::make_shared<quest::IntersectionShaper>(params.policy, arrayAllocId, shapeSet, *compMeshNode);
  }
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    shaper = std::make_shared<quest::IntersectionShaper>(params.policy,
                                                         arrayAllocId,
                                                         shapeSet,
                                                         shapingDC.get());
  }
#endif
  SLIC_ASSERT(shaper != nullptr);
  // shaper->setExecPolicy(params.policy);

  // Set generic parameters for the base Shaper instance
  shaper->setVertexWeldThreshold(params.weldThresh);
  shaper->setVerbosity(params.isVerbose());
  if(params.percentError > 0.)
  {
    shaper->setPercentError(params.percentError);
    shaper->setRefinementType(quest::DiscreteShape::RefinementDynamic);
  }

  // Associate any fields that begin with "vol_frac" with "material" so when
  // the data collection is written, a matset will be created.
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    shaper->getDC()->AssociateMaterialSet("vol_frac", "material");
  }
#endif

  // Set specific parameters here for IntersectionShaper
  shaper->setLevel(params.refinementLevel);
  SLIC_INFO(axom::fmt::format("{:-^80}",
                              axom::fmt::format("Setting IntersectionShaper policy to '{}'",
                                                axom::runtime_policy::policyToName(params.policy))));

  if(!params.backgroundMaterial.empty())
  {
    shaper->setFreeMaterialName(params.backgroundMaterial);
  }

  AXOM_ANNOTATE_END("setup shaping problem");
  AXOM_ANNOTATE_END("init");

  //---------------------------------------------------------------------------
  // Process each of the shapes
  //---------------------------------------------------------------------------
  SLIC_INFO(axom::fmt::format("{:=^80}", "Shaping loop"));
  AXOM_ANNOTATE_BEGIN("shaping");
  for(const auto& shape : shapeSet.getShapes())
  {
    const std::string shapeFormat = shape.getGeometry().getFormat();
    SLIC_INFO(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Processing shape '{}' of material '{}' (format '{}')",
                                          shape.getName(),
                                          shape.getMaterial(),
                                          shapeFormat)));

    // Load the shape from file. This also applies any transformations.
    shaper->loadShape(shape);
    slic::flushStreams();

    // Generate a spatial index over the shape
    shaper->prepareShapeQuery(shapeSet.getDimensions(), shape);
    slic::flushStreams();

    // Query the mesh against this shape
    shaper->runShapeQuery(shape);
    slic::flushStreams();

    // Apply the replacement rules for this shape against the existing materials
    shaper->applyReplacementRules(shape);
    slic::flushStreams();

    // Finalize data structures associated with this shape and spatial index
    shaper->finalizeShapeQuery();
    slic::flushStreams();
  }
  AXOM_ANNOTATE_END("shaping");

  //---------------------------------------------------------------------------
  // After shaping in all shapes, generate/adjust the material volume fractions
  //---------------------------------------------------------------------------
  AXOM_ANNOTATE_BEGIN("adjust");
  SLIC_INFO(axom::fmt::format("{:=^80}", "Generating volume fraction fields for materials"));

  shaper->adjustVolumeFractions();

  //---------------------------------------------------------------------------
  // Compute and print volumes of each material's volume fraction
  //---------------------------------------------------------------------------
  using axom::utilities::string::startsWith;
  if(params.useBlueprintSidre() || params.useBlueprintConduit())
  {
#if defined(AXOM_DEBUG)
    conduit::Node info;
    SLIC_ASSERT(conduit::blueprint::mesh::verify(*compMeshNode, info));
#endif
    std::vector<std::string> materialNames = shaper->getMaterialNames();
    for(const auto& materialName : materialNames)
    {
      // Compute and print volume of material.
      const double volume = sumMaterialVolumes(compMeshGrp, materialName);
      SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                  "Volume of material '{}' is {:.6Lf}",
                                  materialName,
                                  volume));
    }
  }
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    for(auto& kv : shaper->getDC()->GetFieldMap())
    {
      if(startsWith(kv.first, "vol_frac_"))
      {
        const auto mat_name = kv.first.substr(9);
        auto* gf = kv.second;

        mfem::ConstantCoefficient one(1.0);
        mfem::LinearForm vol_form(gf->FESpace());
        vol_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
        vol_form.Assemble();

        const double volume = shaper->allReduceSum(*gf * vol_form);

        SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                    "Volume of material '{}' is {:.6Lf}",
                                    mat_name,
                                    volume));
      }
    }
  }
#endif
  AXOM_ANNOTATE_END("adjust");

  int failCounts = 0;

  std::vector<std::string> allVfNames(1,
                                      "free");  // All volume fraction names plus "free".
  for(const auto& shape : shapeSet.getShapes())
  {
    allVfNames.push_back(shape.getMaterial());
  }

  // For error checking, work on host.
  using ExecSpace = axom::SEQ_EXEC;

  //---------------------------------------------------------------------------
  // Correctness test: volume fractions should be in [0,1].
  //---------------------------------------------------------------------------
  axom::ReduceSum<ExecSpace, axom::IndexType> rangeViolationCount(0);
  for(const auto& vfName : allVfNames)
  {
    std::string fieldName = "vol_frac_" + vfName;
    axom::ArrayView<double> vfView = getFieldAsArrayView(fieldName);
    axom::Array<double> vfHostArray(vfView, axom::execution_space<ExecSpace>::allocatorID());
    vfView = vfHostArray.view();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        bool bad = vfView[i] < 0.0 || vfView[i] > 1.0;
        rangeViolationCount += bad;
      });
  }

  failCounts += (rangeViolationCount.get() != 0);

  SLIC_INFO(axom::fmt::format("{:-^80}",
                              axom::fmt::format("Count of volume fractions outside of [0,1]: {}.",
                                                rangeViolationCount.get())));
  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Correctness test: volume fractions in each cell should sum to 1.0.
  //---------------------------------------------------------------------------
  axom::Array<double> volSums(cellCount);
  volSums.fill(0.0);
  axom::ArrayView<double> volSumsView = volSums.view();
  axom::ReduceSum<ExecSpace, axom::IndexType> nonUnitSums(0);
  for(const auto& vfName : allVfNames)
  {
    std::string fieldName = "vol_frac_" + vfName;
    axom::ArrayView<double> vfView = getFieldAsArrayView(fieldName);
    axom::Array<double> vfHostArray(vfView, axom::execution_space<ExecSpace>::allocatorID());
    vfView = vfHostArray.view();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) { volSumsView[i] += vfView[i]; });
  }
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType i) {
      bool bad = !axom::utilities::isNearlyEqual(volSumsView[i], 1.0);
      nonUnitSums += bad;
    });

  failCounts += (nonUnitSums.get() != 0);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Count non-unit volume fraction sums: {}.", nonUnitSums.get())));
  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Correctness test: shape volume in shapingMesh should match volume of the
  // shape mesh for closes shape.  As long as the shapes don't overlap, this
  // should be a good correctness check.
  //---------------------------------------------------------------------------
  auto* meshVerificationGroup = ds.getRoot()->createGroup("meshVerification");
  for(const auto& shape : shapeSet.getShapes())
  {
    axom::quest::DiscreteShape dShape(shape, meshVerificationGroup);
    auto shapeMesh =
      std::dynamic_pointer_cast<axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>>(
        dShape.createMeshRepresentation());
    double shapeMeshVol =
      params.getBoxDim() == 3 ? volumeOfTetMesh(*shapeMesh) : areaOfTriMesh(*shapeMesh);
    SLIC_INFO(axom::fmt::format("{:-^80}",
                                axom::fmt::format("Shape '{}' discrete geometry has {} cells",
                                                  shape.getName(),
                                                  shapeMesh->getNumberOfCells())));

    const std::string& materialName = shape.getMaterial();
    double shapeVol = -1;
    if(params.useBlueprintSidre() || params.useBlueprintConduit())
    {
      shapeVol = sumMaterialVolumes(compMeshGrp, materialName);
    }
#if defined(AXOM_USE_MFEM)
    if(params.useMfem())
    {
      shapeVol = sumMaterialVolumes<axom::SEQ_EXEC>(shapingDC.get(), materialName);
    }
#endif
    double correctShapeVol = params.testShape == "plane" ? params.boxMeshVolume() / 2 : shapeMeshVol;
    SLIC_ASSERT(correctShapeVol > 0.0);  // Indicates error in the test setup.
    double diff = shapeVol - correctShapeVol;

    bool err = !axom::utilities::isNearlyEqualRelative(shapeVol, correctShapeVol, 1e-6, 1e-8);
    failCounts += err;

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format("Material '{}' in shape '{}' has volume {} vs {}, diff of {}, {}.",
                        materialName,
                        shape.getName(),
                        shapeVol,
                        correctShapeVol,
                        diff,
                        (err ? "ERROR" : "OK"))));
  }
  slic::flushStreams();
  ds.getRoot()->destroyGroupAndData("meshVerification");

  //---------------------------------------------------------------------------
  // Save meshes and fields
  //---------------------------------------------------------------------------

  if(!params.outputFile.empty())
  {
    std::string fileName = params.outputFile + ".volfracs";
    if(params.useBlueprintSidre())
    {
      saveMesh(*compMeshGrp, fileName);
      SLIC_INFO(axom::fmt::format("{:-^80}", "Wrote output mesh " + fileName));
    }
#if defined(AXOM_USE_MFEM)
    else
    {
      shaper->getDC()->Save(fileName, sidre::Group::getDefaultIOProtocol());
    }
#endif
  }

#if defined(AXOM_USE_MFEM)
  shapingDC.reset();
#endif
  shaper.reset();

  //---------------------------------------------------------------------------
  // Cleanup and exit
  //---------------------------------------------------------------------------
  SLIC_INFO(axom::fmt::format("{:-^80}", ""));
  slic::flushStreams();

  finalizeLogger();

  return failCounts;
}
