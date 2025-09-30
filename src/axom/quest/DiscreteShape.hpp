// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_DISCRETE_SHAPE_HPP
#define AXOM_QUEST_DISCRETE_SHAPE_HPP

#include <string>
#include <vector>

#include "axom/config.hpp"
#include "axom/klee/Shape.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"

#if defined(AXOM_USE_MPI)
  #include "mpi.h"
#endif

namespace axom
{
namespace quest
{
/*!
  @brief Post-processed klee::Shape, with the geometry discretized and transformed
  according to the Shape's operators.

  TODO: Move this class into internal namespace.
*/
class DiscreteShape
{
public:
  /// Refinement type.
  using RefinementType = enum { RefinementUniformSegments, RefinementDynamic };

  using Point3D = axom::primal::Point<double, 3>;
  using Vector3D = axom::primal::Vector<double, 3>;
  using TetType = axom::primal::Tetrahedron<double, 3>;
  using OctType = axom::primal::Octahedron<double, 3>;
  using HexType = axom::primal::Hexahedron<double, 3>;

  // Common type for the mesh approximation of the shape.
  using TetMesh = axom::mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>;

  static constexpr int DEFAULT_SAMPLES_PER_KNOT_SPAN {25};
  static constexpr double MINIMUM_PERCENT_ERROR {0.};
  static constexpr double MAXIMUM_PERCENT_ERROR {100.};
  static constexpr double DEFAULT_VERTEX_WELD_THRESHOLD {1e-9};

  /*!
    @brief Constructor.

    @param shape The Klee specifications for the shape.
    @param parentGroup Group under which to put the discrete mesh
      and support blueprint-tets shapes.
      If null, don't use sidre and don't support blueprint-tets.
    @param prefixPath Path prefix for shape from a file specified
      with a relative path.

    Refinement type is set to DiscreteShape::RefinementUniformSegments
    and percent error is set to 0.  See setPercentError() and
    setRefinementType().
  */
  DiscreteShape(const axom::klee::Shape& shape,
                axom::sidre::Group* parentGroup,
                const std::string& prefixPath = {});

  virtual ~DiscreteShape() { clearInternalData(); }

  ///@{
  //! @name Functions to get and set shaping parameters

  //! @brief Set prefix for shape files specified as relative path.
  void setPrefixPath(const std::string& prefixPath);

  /*!
    @brief Set the refinement type.
    Refinement type is used for shaping with C2C contours.
  */
  void setRefinementType(RefinementType refinementType) { m_refinementType = refinementType; }

  void setSamplesPerKnotSpan(int nSamples);
  void setVertexWeldThreshold(double threshold);
  /*!
    @brief Set the percentage error tolerance.

    If percent <= MINIMUM_PERCENT_ERROR, the refinement type
    will change to DiscreteShape::RefinementUniformSegments.
  */
  void setPercentError(double percent);

  ///@}

#if defined(AXOM_USE_MPI)
  /**
   @brief Set the MPI communicator used when reading C2C files.
   */
  void setMPICommunicator(MPI_Comm comm) { m_comm = comm; }
#endif

  /*!
    Get the name of this shape.
    \return the shape's name
  */
  const axom::klee::Shape& getShape() const { return m_shape; }

  /*!
    \brief Get the discrete mesh representation.

    If the sidre parent group was used in the constructor, the
    mesh data is stored under that group.

    If the discrete mesh isn't generated yet (for analytical shapes),
    generate it.
  */
  std::shared_ptr<mint::Mesh> createMeshRepresentation();

  //!@brief Get the discrete mesh representation.
  std::shared_ptr<mint::Mesh> getMeshRepresentation() const { return m_meshRep; }

  /*!
     \brief Get the revolved volume for volumes of revolution,
     which is non-zero only for shapes from C2C contours.
  */
  double getRevolvedVolume() const { return m_revolvedVolume; }

private:
  const axom::klee::Shape& m_shape;

  /*!
    \brief Discrete mesh representation.

    This is either an internally generated mesh representing a
    discretized analytical shape or a modifiable copy of the
    discrete input geometry for geometries specified as a
    discrete mesh.
  */
  std::shared_ptr<axom::mint::Mesh> m_meshRep;

  //!@brief Internal DataStore for working space
  axom::sidre::DataStore m_dataStore;

  //!@brief Sidre store for m_meshRep.
  axom::sidre::Group* m_sidreGroup {nullptr};

  //!@brief Prefix for disc files with relative path.
  std::string m_prefixPath;

  ///@{
  //!@name Various parameters for discretization of analytical shapes.
  RefinementType m_refinementType {RefinementUniformSegments};
  double m_percentError {MINIMUM_PERCENT_ERROR};
  int m_samplesPerKnotSpan {DEFAULT_SAMPLES_PER_KNOT_SPAN};
  double m_vertexWeldThreshold {DEFAULT_VERTEX_WELD_THRESHOLD};
  double m_revolvedVolume {0.0};
  ///@}

#if defined(AXOM_USE_MPI)
  MPI_Comm m_comm {MPI_COMM_WORLD};
#endif

  void applyTransforms();
  numerics::Matrix<double> getTransforms() const;

  /*!
    @brief Set the parent group for this object to store data.
  */
  void setParentGroup(axom::sidre::Group* parentGroup);

  //!@brief Return a 3x3 matrix that rotate coordinates from the x-axis to the given direction.
  numerics::Matrix<double> sorAxisRotMatrix(const Vector3D& dir);

  void clearInternalData();

public:
  // These are public only for the CUDA device compiler.
  //!@brief Create the internal mesh representation of the user's tet mesh.
  void createRepresentationOfBlueprintTets();
  //!@brief Create the internal mesh representation of the analytical tetrahedron.
  void createRepresentationOfTet();
  //!@brief Create the internal mesh representation of the analytical hexahedron.
  void createRepresentationOfHex();
  //!@brief Create the internal mesh representation of the analytical plane.
  void createRepresentationOfPlane();
  //!@brief Create the internal mesh representation of the analytical sphere.
  void createRepresentationOfSphere();
  //!@brief Create the internal mesh representation of the analytical SOR.
  void createRepresentationOfSOR();
};

}  // namespace quest
}  // namespace axom

#endif
