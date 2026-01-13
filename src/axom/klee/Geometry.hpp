// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_GEOMETRY_HPP_
#define AXOM_KLEE_GEOMETRY_HPP_

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Units.hpp"

#include "axom/primal.hpp"
#include "axom/sidre/core/Group.hpp"

#include <memory>
#include <string>

namespace axom
{
namespace klee
{
class GeometryOperator;

/// Properties of a geometric object which can be transformed by operators
struct TransformableGeometryProperties
{
  Dimensions dimensions;
  LengthUnit units;
};

/**
 * Compare transformable properties for equality.
 * \param lhs the left-hand-side operand
 * \param rhs the right-hand-side operand
 * \return true if and only if all properties are equal
 */
bool operator==(const TransformableGeometryProperties &lhs,
                const TransformableGeometryProperties &rhs);

/**
 * Compare transformable properties for inequality.
 * \param lhs the left-hand-side operand
 * \param rhs the right-hand-side operand
 * \return false if and only if all properties are equal
 */
inline bool operator!=(const TransformableGeometryProperties &lhs,
                       const TransformableGeometryProperties &rhs)
{
  return !(lhs == rhs);
}

/// Represents the geometry specified in a Shape.
class Geometry
{
public:
  using Point3D = axom::primal::Point<double, 3>;
  using Vector3D = axom::primal::Vector<double, 3>;
  using Sphere3D = axom::primal::Sphere<double, 3>;
  using Cone3D = axom::primal::Cone<double, 3>;
  using Tet3D = axom::primal::Tetrahedron<double, 3>;
  using Hex3D = axom::primal::Hexahedron<double, 3>;
  using Plane3D = axom::primal::Plane<double, 3>;

  /**
   * Create a Geometry object based on a file representation.
   *
   * \param startProperties the transformable properties before any operators are applied
   * \param format the format of the file
   * \param path the path of the file
   * \param operator_ a possibly null operator to apply to the geometry.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           std::string format,
           std::string path,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a Geometry object based on a blueprint tetrahedral mesh.
   *
   * \param startProperties the transformable properties before any operators are applied
   * \param simplexMeshGroup the geometry in blueprint format.
   *   The elements should be segments, triangles or tetrahedra.
   * \param topology The blueprint topology to use.
   * \param operator_ a possibly null operator to apply to the geometry.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           const axom::sidre::Group *simplexMeshGroup,
           const std::string &topology,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a tetrahedron Geometry object.
   *
   * \param startProperties the transformable properties before any operators are applied
   * \param tet Tetrahedron
   * \param operator_ a possibly null operator to apply to the geometry.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           const axom::primal::Tetrahedron<double, 3> &tet,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a hexahedron Geometry object.
   *
   * \param startProperties the transformable properties before any operators are applied
   * \param hex Hexahedron
   * \param operator_ a possibly null operator to apply to the geometry.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           const axom::primal::Hexahedron<double, 3> &hex,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a sphere Geometry object.
   *
   * \param startProperties the transformable properties before any operators are applied
   * \param sphere Analytical sphere specifications
   * \param levelOfRefinement Number of refinement levels to use for discretizing the sphere.
   * \param operator_ a possibly null operator to apply to the geometry.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           const axom::primal::Sphere<double, 3> &sphere,
           axom::IndexType levelOfRefinement,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a surface-of-revolution (SOR) Geometry object.
   *
   * \param startProperties the transformable properties before any operators are applied
   * \param discreteFunction Discrete function describing the surface of revolution.
   * \param sorOrigin 3D coordinates of the point (z=0, r=0).
   * \param sorDirection SOR axis, in the direction of increasing z.
   * \param levelOfRefinement Number of refinement levels to use for discretizing the SOR.
   * \param operator_ a possibly null operator to apply to the geometry.
   *
   * The \c discreteFunction should be an Nx2 array, interpreted as
   * (z,r) pairs, where z is the axial distance and r is the radius.
   *
   * \c sorAxis should point in the direction of increasing z.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           axom::ArrayView<const double, 2> discreteFunction,
           const Point3D &sorOrigin,
           const Vector3D &sorDirection,
           axom::IndexType levelOfRefinement,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a cone Geometry object.
   *
   * \param startProperties the transformable properties before any
   * operators are applied
   * \param sphere Analytical sphere specifications
   * \param levelOfRefinement Number of refinement levels to use for
   *        discretizing the sphere.
   * \param operator_ a possibly null operator to apply to the geometry.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           const axom::primal::Cone<double, 3> &cone,
           axom::IndexType levelOfRefinement,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a planar Geometry object.
   *
   * \param startProperties the transformable properties before any operators are applied
   * \param tet Tetrahedron
   * \param operator_ a possibly null operator to apply to the geometry.
   *
   * The space on the positive normal side of the plane is considered "inside the shape".
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           const axom::primal::Plane<double, 3> &plane,
           std::shared_ptr<GeometryOperator const> operator_);

  /*!
   * @brief Geometry definition in hierarchical format.
   */
  const conduit::Node &asHierarchy() const { return m_geomInfo; }

  conduit::Node &asHierarchy() { return m_geomInfo; }

  /**
   * \brief Get the format in which the geometry was specified.
   *
   * The format is determined by the constructor used.
   * Values are:
   * - "c2c" = C2C file
   * - "proe" = ProE file
   * - "mfem" = MFEM file
   * - "blueprint-tets" = Blueprint tetrahedral mesh in memory
   * - "tet3D" = 3D tetrahedron (4 points)
   * - "sphere3D" = 3D sphere, as \c primal::Sphere<double,3>
   * - "sor3D" = 3D surface of revolution.
   * - "cone3D" = 3D cone, as \c primal::Cone<double,3>
   * - "hex3D" = 3D hexahedron (8 points)
   * - "plane3D" = 3D plane
   *
   * \return the format of the shape
   *
   * \internal TODO: Put all geometry-specific parameters in m_geomInfo, and
   * deprecate geometry-specific interfaces, so new shapes can be added
   * without modifying this code.
   */
  const std::string &getFormat() const { return m_format; }

  /**
   * Get the path at which to find the specification of the geometry,
   * for geometries stored in files.
   *
   * \return the path to the geometry file
   */
  const std::string &getPath() const { return m_path; }

  /**
   * Returns the dimensions of the geometry before applying operators
   * 
   * For file-based inputs, this is the dimension of the input mesh
   */
  Dimensions getInputDimensions() const { return getStartProperties().dimensions; }

  /// Returns the dimensions of the geometry after applying operators
  Dimensions getOutputDimensions() const { return getEndProperties().dimensions; }

  /**
   * \brief Return the blueprint mesh, for formats that are specified
   * by a blueprint mesh or have been converted to a blueprint mesh.
   */
  const axom::sidre::Group *getBlueprintMesh() const;

  /**
   * \brief Return the blueprint mesh topology, for formats that are specified
   * by a blueprint mesh or have been converted to a blueprint mesh.
   */
  const std::string &getBlueprintTopology() const;

  /// \brief Return the SOR axis direction.
  const Vector3D getSorDirection() const { return m_sorDirection; }

  /// \brief Return the 3D coordinates of the point (z=0, r=0)
  const Point3D getSorOriginCoords() const { return m_sorOrigin; }

  /*
   *  \brief Predicate that returns true when the shape has an associated geometry
   *
   *  A false means that this is set up to determine volume fractions without
   *  computing on the geometry.
   *
   *  TODO: We should just create a new format to represent getting
   *  volume fractions without geometries.  Or move this logic into
   *  Shape, because it's confusing to have a Geometry that has no geometry.
  */
  bool hasGeometry() const;

  /**
   * Get a GeometryOperator to apply to this geometry. Can be null.
   *
   * \return a potentially null operator to apply to the geometry
   */
  std::shared_ptr<GeometryOperator const> const &getGeometryOperator() const { return m_operator; }

  /**
   * Get the initial transformable properties of this geometry
   *
   * \return the initial transformable properties of this geometry
   */
  const TransformableGeometryProperties &getStartProperties() const { return m_startProperties; }

  /**
   * Get the final transformable properties of this geometry after operators are applied
   *
   * \return the initial transformable properties of this geometry
   */
  TransformableGeometryProperties getEndProperties() const;

  /**
   * @brief Return the number of levels of refinement for discretization
   * of analytical curves.
   *
   * This number is unused for geometries that are specified in discrete
   * form.
   */
  axom::IndexType getLevelOfRefinement() const { return m_levelOfRefinement; }

  /**
   * @brief Return the tet geometry, when the Geometry
   * represents a tetrahedron.
   */
  const axom::primal::Tetrahedron<double, 3> &getTet() const { return m_tet; }

  /**
   * @brief Return the hex geometry, when the Geometry
   * represents a hexahedron.
   */
  const axom::primal::Hexahedron<double, 3> &getHex() const { return m_hex; }

  /**
   * @brief Return the sphere geometry, when the Geometry
   * represents an alalytical sphere.
   */
  const axom::primal::Sphere<double, 3> &getSphere() const { return m_sphere; }

  /**
   * @brief Return the cone geometry, when the Geometry
   * represents an alalytical cone.
   */
  const axom::primal::Cone<double, 3> &getCone() const { return m_cone; }

  /**
   * @brief Return the plane geometry, when the Geometry
   * represents a plane.
   */
  const axom::primal::Plane<double, 3> &getPlane() const { return m_plane; }

  /**
   * @brief Get the discrete function used in surfaces of revolution.
   */
  axom::ArrayView<const double, 2> getDiscreteFunction() const { return m_discreteFunction.view(); }

private:
  TransformableGeometryProperties m_startProperties;

  /// \brief Geometry info in hierarchical format.
  conduit::Node m_geomInfo;

  /// \brief Geometry format.
  std::string m_format;

  /// \brief Geometry file path, if it's file-based.
  std::string m_path;

  /// \brief Geometry blueprint simplex mesh, when/if it's in memory.
  const axom::sidre::Group *m_meshGroup {nullptr};

  /// \brief Topology of the blueprint simplex mesh, if it's in memory.
  std::string m_topology;

  /// \brief The tetrahedron, if used.
  Tet3D m_tet;

  /// \brief The hexahedron, if used.
  Hex3D m_hex;

  /// \brief The plane, if used.
  Plane3D m_plane;

  /// \brief The analytical sphere, if used.
  Sphere3D m_sphere;

  /// @brief The analytical cone (or cylinder), if used.
  Cone3D m_cone;

  /// \brief The discrete 2D function, as an Nx2 array, if used.
  axom::Array<double, 2> m_discreteFunction;

  /// \brief The point corresponding to z=0 on the SOR axis.
  Point3D m_sorOrigin;

  /// \brief SOR axis in the direction of increasing z.
  Vector3D m_sorDirection;

  /// \brief Level of refinement for discretizing curved analytical shapes and surfaces of revolutions.
  axom::IndexType m_levelOfRefinement {0};

  std::shared_ptr<const GeometryOperator> m_operator;

  /*!
   * @brief Populate m_geomInfo with the geometry definition.
   *
   * Representing geometries in m_geomInfo is a step away from
   * geometry-specific constructors and methods like @c getTet(), @c
   * getHex() and @c getSphere() and toward a uniform interface for
   * providing geometry info.
   */
  void populateGeomInfo();
};

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_GEOMETRY_HPP_
