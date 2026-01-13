// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/DiscreteShape.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace axom
{
namespace quest
{
namespace internal
{
/*!
 * \brief Implementation of a GeometryOperatorVisitor for processing klee shape operators
 *
 * This class extracts the matrix form of supported operators and marks the operator as unvalid otherwise
 * To use, check the \a isValid() function after visiting and then call the \a getMatrix() function.
 */
class AffineMatrixVisitor : public klee::GeometryOperatorVisitor
{
public:
  AffineMatrixVisitor() : m_matrix(4, 4) { }

  void visit(const klee::Translation& translation) override
  {
    m_matrix = translation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Rotation& rotation) override
  {
    m_matrix = rotation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Scale& scale) override
  {
    m_matrix = scale.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::UnitConverter& converter) override
  {
    m_matrix = converter.toMatrix();
    m_isValid = true;
  }

  void visit(const klee::CompositeOperator&) override
  {
    SLIC_WARNING_ROOT("CompositeOperator not supported for Shaper query");
    m_isValid = false;
  }
  void visit(const klee::SliceOperator&) override
  {
    SLIC_WARNING_ROOT("SliceOperator not yet supported for Shaper query");
    m_isValid = false;
  }

  const numerics::Matrix<double>& getMatrix() const { return m_matrix; }
  bool isValid() const { return m_isValid; }

private:
  bool m_isValid {false};
  numerics::Matrix<double> m_matrix;
};

}  // end namespace internal

// TODO: These were needed for linking - but why? They are constexpr.
constexpr int DiscreteShape::DEFAULT_SAMPLES_PER_KNOT_SPAN;
constexpr double DiscreteShape::MINIMUM_PERCENT_ERROR;
constexpr double DiscreteShape::MAXIMUM_PERCENT_ERROR;
constexpr double DiscreteShape::DEFAULT_VERTEX_WELD_THRESHOLD;

DiscreteShape::DiscreteShape(const axom::klee::Shape& shape,
                             axom::sidre::Group* parentGroup,
                             const std::string& prefixPath)
  : m_shape(shape)
  , m_sidreGroup(nullptr)
  , m_refinementType(DiscreteShape::RefinementUniformSegments)
  , m_percentError(utilities::clampVal(0.0, MINIMUM_PERCENT_ERROR, MAXIMUM_PERCENT_ERROR))
{
  setPrefixPath(prefixPath);
  setParentGroup(parentGroup);

  if(parentGroup == nullptr && m_shape.getGeometry().getFormat() == "blueprint-tets")
  {
    SLIC_ERROR(
      "DiscreteShape: Support for Blueprint-mesh shape format currently "
      "requires a non-null parentGroup in the constructor.  This restriction "
      "can be lifted with additional coding, so please file an issue with the "
      "Axom team if you need this feature.");
  }
}

void DiscreteShape::clearInternalData() { m_meshRep.reset(); }

std::shared_ptr<mint::Mesh> DiscreteShape::createMeshRepresentation()
{
  if(m_meshRep)
  {
    return m_meshRep;
  }

  const axom::klee::Geometry& geometry = m_shape.getGeometry();
  const auto& geometryFormat = geometry.getFormat();

  if(geometryFormat == "blueprint-tets")
  {
    createRepresentationOfBlueprintTets();
  }
  else if(geometryFormat == "tet3D")
  {
    createRepresentationOfTet();
  }
  else if(geometryFormat == "hex3D")
  {
    createRepresentationOfHex();
  }
  else if(geometryFormat == "plane3D")
  {
    createRepresentationOfPlane();
  }
  else if(geometryFormat == "sphere3D")
  {
    createRepresentationOfSphere();
  }
  if(geometryFormat == "sor3D")
  {
    createRepresentationOfSOR();
  }

  if(m_meshRep)
  {
    return m_meshRep;
  }

  if(!m_shape.getGeometry().hasGeometry())
  {
    // If shape has no geometry, there's nothing to discretize.
    SLIC_DEBUG_ROOT(axom::fmt::format("Current shape '{}' of material '{}' has no geometry",
                                      m_shape.getName(),
                                      m_shape.getMaterial()));
    return m_meshRep;
  }

  // We handled all the non-file formats.  The rest are file formats.
  const std::string& file_format = geometryFormat;

  std::string shapePath =
    axom::utilities::filesystem::prefixRelativePath(m_shape.getGeometry().getPath(), m_prefixPath);
  SLIC_INFO_ROOT("Reading file: " << shapePath << "...");

  // Initialize revolved volume.
  m_revolvedVolume = 0.;

  if(utilities::string::endsWith(shapePath, ".stl"))
  {
    SLIC_ERROR_ROOT_IF(file_format != "stl",
                       axom::fmt::format(" '{}' format requires .stl file type", file_format));

    axom::mint::Mesh* meshRep = nullptr;
#ifdef AXOM_USE_MPI
    const int rc = quest::internal::read_stl_mesh(shapePath, meshRep, m_comm);
#else
    const int rc = quest::internal::read_stl_mesh(shapePath, meshRep);
#endif
    SLIC_ERROR_ROOT_IF(rc != quest::internal::READ_SUCCESS,
                       axom::fmt::format("Failed to read STL shape '{}' from file '{}'.",
                                         m_shape.getName(),
                                         shapePath));

    m_meshRep.reset(meshRep);
    // Transform the coordinates of the linearized mesh.
    applyTransforms();
  }
  else if(utilities::string::endsWith(shapePath, ".proe"))
  {
    SLIC_ERROR_ROOT_IF(file_format != "proe",
                       axom::fmt::format(" '{}' format requires .proe file type", file_format));

    axom::mint::Mesh* meshRep = nullptr;
#ifdef AXOM_USE_MPI
    const int rc = quest::internal::read_pro_e_mesh(shapePath, meshRep, m_comm);
#else
    const int rc = quest::internal::read_pro_e_mesh(shapePath, meshRep);
#endif
    SLIC_ERROR_ROOT_IF(rc != quest::internal::READ_SUCCESS,
                       axom::fmt::format("Failed to read Pro/E shape '{}' from file '{}'.",
                                         m_shape.getName(),
                                         shapePath));

    m_meshRep.reset(meshRep);
  }
#ifdef AXOM_USE_C2C
  else if(utilities::string::endsWith(shapePath, ".contour"))
  {
    SLIC_ERROR_ROOT_IF(file_format != "c2c",
                       axom::fmt::format(" '{}' format requires .contour file type", file_format));

    // Get the transforms that are being applied to the mesh. Get them
    // as a single concatenated matrix.
    auto transform = getTransforms();

    // Pass in the transform so any transformations can figure into computing the revolved volume.
    axom::mint::Mesh* meshRep = nullptr;
    const bool uniform = !(m_refinementType == DiscreteShape::RefinementDynamic &&
                           m_percentError > MINIMUM_PERCENT_ERROR);
  #ifdef AXOM_USE_MPI
    int rc = quest::internal::READ_FAILED;
    try
    {
      rc = quest::internal::read_c2c_mesh(shapePath,
                                          uniform,
                                          transform,
                                          m_samplesPerKnotSpan,
                                          m_vertexWeldThreshold,
                                          m_percentError,
                                          meshRep,
                                          m_revolvedVolume,  // output arg
                                          m_comm);
    }
    catch(const std::exception& e)
    {
      SLIC_ERROR_ROOT(
        axom::fmt::format("Failed to read C2C shape '{}' from file '{}'. Exception: {}",
                          m_shape.getName(),
                          shapePath,
                          e.what()));
    }
    catch(...)
    {
      SLIC_ERROR_ROOT(axom::fmt::format("Failed to read C2C shape '{}' from file '{}'.",
                                        m_shape.getName(),
                                        shapePath));
    }
  #else
    int rc = quest::internal::READ_FAILED;
    try
    {
      rc = quest::internal::read_c2c_mesh(shapePath,
                                          uniform,
                                          transform,
                                          m_samplesPerKnotSpan,
                                          m_vertexWeldThreshold,
                                          m_percentError,
                                          meshRep,
                                          m_revolvedVolume);  // output arg
    }
    catch(const std::exception& e)
    {
      SLIC_ERROR_ROOT(
        axom::fmt::format("Failed to read C2C shape '{}' from file '{}'. Exception: {}",
                          m_shape.getName(),
                          shapePath,
                          e.what()));
    }
    catch(...)
    {
      SLIC_ERROR_ROOT(axom::fmt::format("Failed to read C2C shape '{}' from file '{}'.",
                                        m_shape.getName(),
                                        shapePath));
    }
  #endif
    SLIC_ERROR_ROOT_IF(rc != quest::internal::READ_SUCCESS,
                       axom::fmt::format("Failed to read C2C shape '{}' from file '{}'.",
                                         m_shape.getName(),
                                         shapePath));

    m_meshRep.reset(meshRep);

    // Transform the coordinates of the linearized mesh.
    applyTransforms();
  }
#endif
#ifdef AXOM_USE_MFEM
  else if(utilities::string::endsWith(shapePath, ".mesh"))
  {
    SLIC_ERROR_ROOT_IF(file_format != "mfem",
                       axom::fmt::format(" '{}' format requires .mesh file extension", file_format));
    // Get the transforms that are being applied to the mesh. Get them
    // as a single concatenated matrix.
    auto transform = getTransforms();

    // Pass in the transform so any transformations can figure into computing the revolved volume.
    axom::mint::Mesh* meshRep = nullptr;
    const bool uniform = !(m_refinementType == DiscreteShape::RefinementDynamic &&
                           m_percentError > MINIMUM_PERCENT_ERROR);
    const int rc = quest::internal::read_mfem_mesh(shapePath,
                                                   uniform,
                                                   transform,
                                                   m_samplesPerKnotSpan,
                                                   m_vertexWeldThreshold,
                                                   m_percentError,
                                                   meshRep,
                                                   m_revolvedVolume);  // output arg

    SLIC_ERROR_ROOT_IF(rc != quest::internal::READ_SUCCESS,
                       axom::fmt::format("Failed to read MFEM shape '{}' from file '{}'.",
                                         m_shape.getName(),
                                         shapePath));

    m_meshRep.reset(meshRep);

    // Transform the coordinates of the linearized mesh.
    applyTransforms();
  }
#endif
  else
  {
    SLIC_ERROR_ROOT(
      axom::fmt::format("Unsupported filetype for this Axom configuration. "
                        "Provided file was '{}', with format '{}'",
                        shapePath,
                        file_format));
  }

  return m_meshRep;
}

void DiscreteShape::createRepresentationOfBlueprintTets()
{
  SLIC_ASSERT(m_sidreGroup != nullptr);

  const axom::klee::Geometry& geometry = m_shape.getGeometry();

  // Put the in-memory geometry in m_meshRep.
  const axom::sidre::Group* inputGroup = geometry.getBlueprintMesh();
#ifdef AXOM_USE_UMPIRE
  int allocID = inputGroup->getDefaultAllocatorID();
#else
  int allocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
#endif

  std::string modName = inputGroup->getName() + "_modified";
  while(m_sidreGroup->hasGroup(modName))
  {
    modName = modName + "-";
  }

  /*
    We use a sidre::Group::deepCopyGroup as a convenient way to copy
    the mesh.  This requires a non-null m_sidreGroup.  This
    restriction can be lifted if we implement code to copy the
    inputGroup directly into an mint::UnstructuredMesh without putting
    it into another sidre::Group.  We currently don't have that.
  */
  axom::sidre::Group* modGroup = m_sidreGroup->createGroup(modName);
  modGroup->deepCopyGroup(inputGroup, allocID);

  m_meshRep.reset(axom::mint::getMesh(modGroup->getGroup(inputGroup->getName()),
                                      m_shape.getGeometry().getBlueprintTopology()));

  // Transform the coordinates of the linearized mesh.
  applyTransforms();
}

void DiscreteShape::createRepresentationOfTet()
{
  const axom::klee::Geometry& geometry = m_shape.getGeometry();

  const auto& tet = geometry.getTet();

  const axom::IndexType tetCount = 1;
  const axom::IndexType nodeCount = 4;
  axom::Array<double, 2> nodeCoords(nodeCount, 3);
  axom::Array<axom::IndexType, 2> connectivity(tetCount, 4);

  for(int iNode = 0; iNode < 4; ++iNode)
  {
    const auto& coords = tet[iNode];
    nodeCoords[iNode][0] = coords[0];
    nodeCoords[iNode][1] = coords[1];
    nodeCoords[iNode][2] = coords[2];
    connectivity[0][iNode] = iNode;
  }

  TetMesh* tetMesh = nullptr;
  if(m_sidreGroup != nullptr)
  {
    tetMesh = new TetMesh(3,
                          axom::mint::CellType::TET,
                          m_sidreGroup,
                          nodeCoords.shape()[0],
                          connectivity.shape()[0]);
  }
  else
  {
    tetMesh =
      new TetMesh(3, axom::mint::CellType::TET, nodeCoords.shape()[0], connectivity.shape()[0]);
  }
  tetMesh->appendNodes((double*)nodeCoords.data(), nodeCoords.shape()[0]);
  tetMesh->appendCells(connectivity.data(), connectivity.shape()[0]);
  m_meshRep.reset(tetMesh);

  applyTransforms();
}

void DiscreteShape::createRepresentationOfHex()
{
  const axom::klee::Geometry& geometry = m_shape.getGeometry();
  const auto& hex = geometry.getHex();
  axom::StackArray<TetType, HexType::NUM_TRIANGULATE> tets;
  hex.triangulate(tets);

  const axom::IndexType tetCount = HexType::NUM_TRIANGULATE;
  const axom::IndexType nodeCount = HexType::NUM_TRIANGULATE * 4;
  axom::Array<double, 2> nodeCoords(nodeCount, 3);
  auto nodeCoordsView = nodeCoords.view();
  axom::Array<axom::IndexType, 2> connectivity(tetCount, 4);
  auto connectivityView = connectivity.view();
  // NOTE: This is not much computation, so just run on host.
  axom::for_all<axom::SEQ_EXEC>(
    tetCount,
    AXOM_LAMBDA(axom::IndexType iTet) {
      const auto& tet = tets[iTet];
      for(int i = 0; i < 4; ++i)
      {
        axom::IndexType iNode = iTet * 4 + i;
        const auto& coords = tet[i];
        nodeCoordsView[iNode][0] = coords[0];
        nodeCoordsView[iNode][1] = coords[1];
        nodeCoordsView[iNode][2] = coords[2];
        connectivityView[iTet][i] = iNode;
      }
    });

  TetMesh* tetMesh = nullptr;
  if(m_sidreGroup != nullptr)
  {
    tetMesh = new TetMesh(3,
                          axom::mint::CellType::TET,
                          m_sidreGroup,
                          nodeCoords.shape()[0],
                          connectivity.shape()[0]);
  }
  else
  {
    tetMesh =
      new TetMesh(3, axom::mint::CellType::TET, nodeCoords.shape()[0], connectivity.shape()[0]);
  }
  tetMesh->appendNodes((double*)nodeCoords.data(), nodeCoords.shape()[0]);
  tetMesh->appendCells(connectivity.data(), connectivity.shape()[0]);
  m_meshRep.reset(tetMesh);

  applyTransforms();
}

void DiscreteShape::createRepresentationOfPlane()
{
  const axom::klee::Geometry& geometry = m_shape.getGeometry();

  const auto& plane = geometry.getPlane();
  // Generate a big bounding hex on the positive side of the plane.
  /*
    TODO: Use a specialized plane representation.  Plane is very
    simple, but representing it as huge tets has 2 problems.
    1. The shaper runs slow because the huge tets contains huge
       numbers of cells, disabling attempts at divide-and-conquer
       with the BVH.

    2. When the tets are orders of magnitude larger then the
       computational mesh, I've seen weird overlap computation
       errors that affect certain resolutions and don't follow
       any trends.  I think it's due to round-off, but I've not
       verified.
  */
  axom::primal::Hexahedron<double, 3> boundingHex;
  const double len = 1e2;  // Big enough to contain anticipated mesh,
                           // but too big invites weird round-off errors.
  // We should compute based on the mesh.
  // clang-format off
  boundingHex[0] = Point3D{0.0, -len, -len};
  boundingHex[1] = Point3D{len, -len, -len};
  boundingHex[2] = Point3D{len,  len, -len};
  boundingHex[3] = Point3D{0.0,  len, -len};
  boundingHex[4] = Point3D{0.0, -len,  len};
  boundingHex[5] = Point3D{len, -len,  len};
  boundingHex[6] = Point3D{len,  len,  len};
  boundingHex[7] = Point3D{0.0,  len,  len};
  // clang-format on

  // Rotate and translate boundingHex to align with the plane.
  numerics::Matrix<double> rotate = sorAxisRotMatrix(plane.getNormal());
  const auto translate = plane.getNormal() * plane.getOffset();
  for(int i = 0; i < 8; ++i)
  {
    Point3D newCoords;
    numerics::matrix_vector_multiply(rotate, boundingHex[i].data(), newCoords.data());
    newCoords.array() += translate.array();
    boundingHex[i].array() = newCoords.array();
  }

  axom::StackArray<TetType, HexType::NUM_TRIANGULATE> tets;
  boundingHex.triangulate(tets);

  const axom::IndexType tetCount = HexType::NUM_TRIANGULATE;
  const axom::IndexType nodeCount = HexType::NUM_TRIANGULATE * 4;
  axom::Array<double, 2> nodeCoords(nodeCount, 3);
  auto nodeCoordsView = nodeCoords.view();
  axom::Array<axom::IndexType, 2> connectivity(tetCount, 4);
  auto connectivityView = connectivity.view();
  // NOTE: This is not much computation, so just run on host.
  axom::for_all<axom::SEQ_EXEC>(
    tetCount,
    AXOM_LAMBDA(axom::IndexType iTet) {
      const auto& tet = tets[iTet];
      for(int i = 0; i < 4; ++i)
      {
        axom::IndexType iNode = iTet * 4 + i;
        const auto& coords = tet[i];
        nodeCoordsView[iNode][0] = coords[0];
        nodeCoordsView[iNode][1] = coords[1];
        nodeCoordsView[iNode][2] = coords[2];
        connectivityView[iTet][i] = iNode;
      }
    });

  TetMesh* tetMesh = nullptr;
  if(m_sidreGroup != nullptr)
  {
    tetMesh = new TetMesh(3,
                          axom::mint::CellType::TET,
                          m_sidreGroup,
                          nodeCoords.shape()[0],
                          connectivity.shape()[0]);
  }
  else
  {
    tetMesh =
      new TetMesh(3, axom::mint::CellType::TET, nodeCoords.shape()[0], connectivity.shape()[0]);
  }
  tetMesh->appendNodes((double*)nodeCoords.data(), nodeCoords.shape()[0]);
  tetMesh->appendCells(connectivity.data(), connectivity.shape()[0]);
  m_meshRep.reset(tetMesh);

  applyTransforms();
}

void DiscreteShape::createRepresentationOfSphere()
{
  const axom::klee::Geometry& geometry = m_shape.getGeometry();

  /*
    Discretize the sphere into m_meshRep and apply transforms:
    1. Discretize the sphere into a set of Octahedra.
    2. Split each Octahedron into 8 tets.
    3. Insert tets into tet mesh.
    4. Transform the tet mesh.
    We can reduce memory use by eliminating repeated points if it
    is a problem.
  */
  const auto& sphere = geometry.getSphere();
  axom::Array<OctType> octs;
  int octCount = 0;
  axom::quest::discretize(sphere, geometry.getLevelOfRefinement(), octs, octCount);

  // Note: MSVC required `static` to pass TETS_PER_OCT to the axom::for_all w/ C++14
  // this restriction might be lifted w/ C++17
  static constexpr int TETS_PER_OCT = 8;
  static constexpr int NODES_PER_TET = 4;
  const axom::IndexType tetCount = octCount * TETS_PER_OCT;
  const axom::IndexType nodeCount = tetCount * NODES_PER_TET;

  axom::Array<Point3D> nodeCoords(nodeCount, nodeCount);
  axom::Array<axom::IndexType, 2> connectivity(
    axom::StackArray<axom::IndexType, 2> {tetCount, NODES_PER_TET});

  auto nodeCoordsView = nodeCoords.view();
  auto connectivityView = connectivity.view();
  axom::for_all<axom::SEQ_EXEC>(
    octCount,
    AXOM_LAMBDA(axom::IndexType octIdx) {
      TetType tetsInOct[TETS_PER_OCT];
      axom::primal::split(octs[octIdx], tetsInOct);
      for(int iTet = 0; iTet < TETS_PER_OCT; ++iTet)
      {
        axom::IndexType tetIdx = octIdx * TETS_PER_OCT + iTet;
        for(int iNode = 0; iNode < NODES_PER_TET; ++iNode)
        {
          axom::IndexType nodeIdx = tetIdx * NODES_PER_TET + iNode;
          nodeCoordsView[nodeIdx] = tetsInOct[iTet][iNode];
          connectivityView[tetIdx][iNode] = nodeIdx;
        }
      }
    });

  TetMesh* tetMesh = nullptr;
  if(m_sidreGroup != nullptr)
  {
    tetMesh = new TetMesh(3, axom::mint::CellType::TET, m_sidreGroup, nodeCount, tetCount);
  }
  else
  {
    tetMesh = new TetMesh(3, axom::mint::CellType::TET, nodeCount, tetCount);
  }
  tetMesh->appendNodes((double*)nodeCoords.data(), nodeCount);
  tetMesh->appendCells(connectivity.data(), tetCount);
  m_meshRep.reset(tetMesh);

  applyTransforms();
}

void DiscreteShape::createRepresentationOfSOR()
{
  // Construct the tet m_meshRep from the volume-of-revolution.
  auto& sorGeom = m_shape.getGeometry();
  const auto& discreteFcn = sorGeom.getDiscreteFunction();

  // Generate the Octahedra
  axom::Array<OctType> octs;
  int octCount = 0;
  axom::ArrayView<Point2D> polyline((Point2D*)discreteFcn.data(), discreteFcn.shape()[0]);
  const bool good =
    axom::quest::discretize<axom::SEQ_EXEC>(polyline,
                                            int(polyline.size()),
                                            m_shape.getGeometry().getLevelOfRefinement(),
                                            octs,
                                            octCount);
  AXOM_UNUSED_VAR(good);
  SLIC_ASSERT(good);

  // Rotate to the SOR axis direction and translate to the base location.
  numerics::Matrix<double> rotate = sorAxisRotMatrix(sorGeom.getSorDirection());
  const auto& translate = sorGeom.getSorOriginCoords();
  auto octsView = octs.view();
  axom::for_all<axom::SEQ_EXEC>(
    octCount,
    AXOM_LAMBDA(axom::IndexType iOct) {
      auto& oct = octsView[iOct];
      for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
      {
        auto& newCoords = oct[iVert];
        auto oldCoords = newCoords;
        numerics::matrix_vector_multiply(rotate, oldCoords.data(), newCoords.data());
        newCoords.array() += translate.array();
      }
    });

  // Dump discretized octs as a tet mesh
  //
  // NOTE: mesh_from_discretized_polyline returns a tet mesh.
  // IntersectionShaper can handle octahedra, which may be faster
  // because it has 8x fewer elements.  This method returns the
  // discretization as a mesh.  But mint doesn't support octs and
  // Blueprint support octs only as polyhedra.  We can always return
  // an array of primal::Octahedron instead of a mesh.
  axom::mint::Mesh* mesh;
  axom::quest::mesh_from_discretized_polyline(octs.view(), octCount, polyline.size() - 1, mesh);

  if(m_sidreGroup)
  {
    // If using sidre, copy the tetMesh into sidre.
    auto* tetMesh = static_cast<TetMesh*>(mesh);
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* siMesh =
      new axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>(tetMesh->getDimension(),
                                                                 tetMesh->getCellType(),
                                                                 m_sidreGroup,
                                                                 tetMesh->getTopologyName(),
                                                                 tetMesh->getCoordsetName(),
                                                                 tetMesh->getNumberOfNodes(),
                                                                 tetMesh->getNumberOfCells());
    siMesh->appendNodes(tetMesh->getCoordinateArray(0),
                        tetMesh->getCoordinateArray(1),
                        tetMesh->getCoordinateArray(2),
                        tetMesh->getNumberOfNodes());
    siMesh->appendCells(tetMesh->getCellNodesArray(), tetMesh->getNumberOfCells());
    m_meshRep.reset(siMesh);
    delete mesh;
    mesh = nullptr;
  }
  else
  {
    m_meshRep.reset(mesh);
  }

  // Transform the coordinates of the linearized mesh.
  applyTransforms();
}

void DiscreteShape::applyTransforms()
{
  numerics::Matrix<double> transformation = getTransforms();

  // Apply transformation to coordinates of each vertex in mesh
  if(!transformation.isIdentity())
  {
    const int spaceDim = m_meshRep->getDimension();
    const int numSurfaceVertices = m_meshRep->getNumberOfNodes();
    double* x = m_meshRep->getCoordinateArray(mint::X_COORDINATE);
    double* y = m_meshRep->getCoordinateArray(mint::Y_COORDINATE);
    double* z = spaceDim > 2 ? m_meshRep->getCoordinateArray(mint::Z_COORDINATE) : nullptr;

    double xformed[4];
    for(int i = 0; i < numSurfaceVertices; ++i)
    {
      double coords[4] = {x[i], y[i], (z == nullptr ? 0. : z[i]), 1.};
      numerics::matrix_vector_multiply(transformation, coords, xformed);
      x[i] = xformed[0];
      y[i] = xformed[1];
      if(z != nullptr)
      {
        z[i] = xformed[2];
      }
    }
  }
}

numerics::Matrix<double> DiscreteShape::getTransforms() const
{
  const auto identity4x4 = numerics::Matrix<double>::identity(4);
  numerics::Matrix<double> transformation(identity4x4);
  auto& geometryOperator = m_shape.getGeometry().getGeometryOperator();
  if(geometryOperator)
  {
    auto composite = std::dynamic_pointer_cast<const klee::CompositeOperator>(geometryOperator);
    if(composite)
    {
      // Concatenate the transformations

      // Why don't we multiply the matrices in CompositeOperator::addOperator()?
      // Why keep the matrices factored and multiply them here repeatedly?
      // Combining them would also avoid this if-else logic.  BTNG
      for(auto op : composite->getOperators())
      {
        // Use visitor pattern to extract the affine matrix from supported operators
        internal::AffineMatrixVisitor visitor;
        op->accept(visitor);
        if(!visitor.isValid())
        {
          continue;
        }
        const auto& matrix = visitor.getMatrix();
        numerics::Matrix<double> res(identity4x4);
        numerics::matrix_multiply(matrix, transformation, res);
        transformation = res;
      }
    }
    else
    {
      internal::AffineMatrixVisitor visitor;
      geometryOperator->accept(visitor);
      if(visitor.isValid())
      {
        transformation = visitor.getMatrix();
      }
    }
  }
  return transformation;
}

// Return a 3x3 matrix that rotates coordinates from the x-axis to the given direction.
numerics::Matrix<double> DiscreteShape::sorAxisRotMatrix(const Vector3D& dir)
{
  // Note that the rotation matrix is not unique.
  static const Vector3D x {1.0, 0.0, 0.0};
  Vector3D a = dir.unitVector();
  Vector3D u;  // Rotation vector, the cross product of x and a.
  numerics::cross_product(x.data(), a.data(), u.data());
  double sinT = u.norm();
  double cosT = numerics::dot_product(x.data(), a.data(), 3);
  double ccosT = 1 - cosT;

  // Degenerate case with angle near 0 or pi.
  if(utilities::isNearlyEqual(sinT, 0.0))
  {
    if(cosT > 0)
    {
      return numerics::Matrix<double>::identity(3);
    }
    else
    {
      // Give u a tiny component in any non-x direction
      // so we can rotate around it.
      u[1] = 1e-8;
    }
  }

  u = u.unitVector();
  numerics::Matrix<double> rot(3, 3, 0.0);
  rot(0, 0) = u[0] * u[0] * ccosT + cosT;
  rot(0, 1) = u[0] * u[1] * ccosT - u[2] * sinT;
  rot(0, 2) = u[0] * u[2] * ccosT + u[1] * sinT;
  rot(1, 0) = u[1] * u[0] * ccosT + u[2] * sinT;
  rot(1, 1) = u[1] * u[1] * ccosT + cosT;
  rot(1, 2) = u[1] * u[2] * ccosT - u[0] * sinT;
  rot(2, 0) = u[2] * u[0] * ccosT - u[1] * sinT;
  rot(2, 1) = u[2] * u[1] * ccosT + u[0] * sinT;
  rot(2, 2) = u[2] * u[2] * ccosT + cosT;

  return rot;
}

void DiscreteShape::setSamplesPerKnotSpan(int nSamples)
{
  using axom::utilities::clampLower;
  SLIC_WARNING_ROOT_IF(
    nSamples < 1,
    axom::fmt::format("Samples per knot span must be at least 1. Provided value was {}", nSamples));

  m_samplesPerKnotSpan = clampLower(nSamples, 1);
}

void DiscreteShape::setVertexWeldThreshold(double threshold)
{
  SLIC_WARNING_ROOT_IF(
    threshold <= 0.,
    axom::fmt::format("Vertex weld threshold should be positive. Provided value was {}", threshold));

  m_vertexWeldThreshold = threshold;
}

void DiscreteShape::setPercentError(double percent)
{
  using axom::utilities::clampVal;
  SLIC_WARNING_ROOT_IF(percent <= MINIMUM_PERCENT_ERROR,
                       axom::fmt::format("Percent error must be greater than {}. Provided value "
                                         "was {}. Dynamic refinement will not be used.",
                                         MINIMUM_PERCENT_ERROR,
                                         percent));
  SLIC_WARNING_ROOT_IF(
    percent > MAXIMUM_PERCENT_ERROR,
    axom::fmt::format("Percent error must be less than {}. Provided value was {}",
                      MAXIMUM_PERCENT_ERROR,
                      percent));
  if(percent <= MINIMUM_PERCENT_ERROR)
  {
    m_refinementType = DiscreteShape::RefinementUniformSegments;
  }
  m_percentError = clampVal(percent, MINIMUM_PERCENT_ERROR, MAXIMUM_PERCENT_ERROR);
}

void DiscreteShape::setPrefixPath(const std::string& prefixPath)
{
  SLIC_ERROR_ROOT_IF(!prefixPath.empty() && !axom::utilities::filesystem::pathExists(prefixPath),
                     axom::fmt::format("Path '{}' does not exist.", prefixPath));
  m_prefixPath = prefixPath;
}

void DiscreteShape::setParentGroup(axom::sidre::Group* parentGroup)
{
  if(parentGroup)
  {
    // Use object address to create a unique name for sidre group under parent.
    std::string myGroupName;
    int i = 0;
    do
    {
      const std::string str = axom::fmt::format("DiscreteShapeTemp-{}", i);
      if(!parentGroup->hasGroup(str))
      {
        myGroupName = str;
      }
      ++i;
    } while(myGroupName.empty());

    m_sidreGroup = parentGroup->createGroup(myGroupName);
  }
}

}  // namespace quest
}  // namespace axom
