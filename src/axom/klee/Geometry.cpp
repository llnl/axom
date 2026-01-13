// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"

#include "conduit_blueprint_mesh.hpp"

#include <utility>

namespace axom
{
namespace klee
{
bool operator==(const TransformableGeometryProperties& lhs, const TransformableGeometryProperties& rhs)
{
  return lhs.dimensions == rhs.dimensions && lhs.units == rhs.units;
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   std::string format,
                   std::string path,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format(std::move(format))
  , m_path(std::move(path))
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const axom::sidre::Group* simplexMeshGroup,
                   const std::string& topology,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("blueprint-tets")
  , m_meshGroup(simplexMeshGroup)
  , m_topology(topology)
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const axom::primal::Tetrahedron<double, 3>& tet,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("tet3D")
  , m_tet(tet)
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const axom::primal::Hexahedron<double, 3>& hex,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("hex3D")
  , m_hex(hex)
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const Sphere3D& sphere,
                   axom::IndexType levelOfRefinement,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("sphere3D")
  , m_sphere(sphere)
  , m_levelOfRefinement(levelOfRefinement)
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const axom::primal::Cone<double, 3>& cone,
                   axom::IndexType levelOfRefinement,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("cone3D")
  , m_path()
  , m_meshGroup(nullptr)
  , m_topology()
  , m_cone(cone)
  , m_levelOfRefinement(levelOfRefinement)
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   axom::ArrayView<const double, 2> discreteFunction,
                   const Point3D& sorOrigin,  // surface of revolution.
                   const Vector3D& sorDirection,
                   axom::IndexType levelOfRefinement,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("sor3D")
  , m_discreteFunction(discreteFunction)
  , m_sorOrigin(sorOrigin)
  , m_sorDirection(sorDirection)
  , m_levelOfRefinement(levelOfRefinement)
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const axom::primal::Plane<double, 3>& plane,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("plane3D")
  , m_plane(plane)
  , m_operator(std::move(operator_))
{
  populateGeomInfo();
}

void Geometry::populateGeomInfo()
{
  if(m_format == "blueprint-tets")
  {
    m_meshGroup->deepCopyToConduit(m_geomInfo["klee::Geometry:tetMesh"]);
    m_geomInfo["topologyName"].set(getBlueprintTopology());
  }

  else if(m_format == "tet3D")
  {
    const auto& tet = getTet();
    m_geomInfo["v0"].set(tet[0].data(), 3);
    m_geomInfo["v1"].set(tet[1].data(), 3);
    m_geomInfo["v2"].set(tet[2].data(), 3);
    m_geomInfo["v3"].set(tet[3].data(), 3);
  }

  else if(m_format == "sphere3D")
  {
    const Sphere3D& sphere = getSphere();
    m_geomInfo["center"].set(sphere.getCenter().data(), 3);
    m_geomInfo["radius"].set(sphere.getRadius());
    m_geomInfo["levelOfRefinement"].set(m_levelOfRefinement);
  }

  else if(m_format == "cone3D")
  {
    const Cone3D& cone = getCone();
    m_discreteFunction = axom::Array<double, 2>(2, 2);
    m_discreteFunction(0, 0) = 0.0;
    m_discreteFunction(0, 1) = cone.getBaseRadius();
    m_discreteFunction(1, 0) = cone.getLength();
    m_discreteFunction(1, 1) = cone.getTopRadius();
    m_geomInfo["discreteFunction"].set(m_discreteFunction.data(), m_discreteFunction.size());
    m_geomInfo["sorOrigin"].set(cone.getBaseCenter().data(), 3);
    m_geomInfo["sorDirection"].set(cone.getDirection().data(), 3);
    m_geomInfo["levelOfRefinement"].set(m_levelOfRefinement);
  }

  else if(m_format == "sor3D")
  {
    m_geomInfo["sorOrigin"].set(m_sorOrigin.data(), 3);
    m_geomInfo["sorDirection"].set(m_sorDirection.data(), 3);
    m_geomInfo["discreteFunction"].set(m_discreteFunction.data(), m_discreteFunction.size());
    m_geomInfo["levelOfRefinement"].set(m_levelOfRefinement);
  }

  else if(m_format == "hex3D")
  {
    const auto& hex = getHex();
    m_geomInfo["v0"].set(hex[0].data(), 3);
    m_geomInfo["v1"].set(hex[1].data(), 3);
    m_geomInfo["v2"].set(hex[2].data(), 3);
    m_geomInfo["v3"].set(hex[3].data(), 3);
    m_geomInfo["v4"].set(hex[4].data(), 3);
    m_geomInfo["v5"].set(hex[5].data(), 3);
    m_geomInfo["v6"].set(hex[6].data(), 3);
    m_geomInfo["v7"].set(hex[7].data(), 3);
  }

  else if(m_format == "plane3D")
  {
    const auto& plane = getPlane();
    m_geomInfo["normal"].set(plane.getNormal().data(), 3);
    m_geomInfo["offset"].set(plane.getOffset());
  }

  // TODO: other formats.
}

bool Geometry::hasGeometry() const
{
  bool isInMemory = (m_format == "blueprint-tets" || m_format == "sphere3D" || m_format == "tet3D" ||
                     m_format == "hex3D" || m_format == "plane3D" || m_format == "cone3D");
  if(isInMemory)
  {
    return true;
  }
  return !m_path.empty();
}

TransformableGeometryProperties Geometry::getEndProperties() const
{
  return m_operator ? m_operator->getEndProperties() : m_startProperties;
}

const axom::sidre::Group* Geometry::getBlueprintMesh() const
{
  SLIC_ASSERT_MSG(m_meshGroup,
                  axom::fmt::format("The Geometry format '{}' is not specified "
                                    "as a blueprint mesh and/or has not been converted into one.",
                                    m_format));
  return m_meshGroup;
}

const std::string& Geometry::getBlueprintTopology() const
{
  SLIC_ASSERT_MSG(m_meshGroup,
                  axom::fmt::format("The Geometry format '{}' is not specified "
                                    "as a blueprint mesh and/or has not been converted into one.",
                                    m_format));
  return m_topology;
}

}  // namespace klee
}  // namespace axom
