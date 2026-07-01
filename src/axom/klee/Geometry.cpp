// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/AffineMatrixVisitor.hpp"
#include "axom/klee/KleeError.hpp"

#include "conduit_blueprint_mesh.hpp"

#include <utility>

namespace axom
{
namespace klee
{
namespace
{
Geometry::Point3D applyMatrix(const numerics::Matrix<double>& matrix, const Geometry::Point3D& point)
{
  double coords[4] = {point[0], point[1], point[2], 1.};
  double xformed[4];
  numerics::matrix_vector_multiply(matrix, coords, xformed);
  return Geometry::Point3D {xformed[0], xformed[1], xformed[2]};
}

bool operatorSequenceHasNonAffine(const std::shared_ptr<const GeometryOperator>& op)
{
  if(!op)
  {
    return false;
  }

  if(auto composite = std::dynamic_pointer_cast<const CompositeOperator>(op))
  {
    for(const auto& child : composite->getOperators())
    {
      if(operatorSequenceHasNonAffine(child))
      {
        return true;
      }
    }
    return false;
  }

  return !std::dynamic_pointer_cast<const MatrixOperator>(op);
}

std::string operatorName(const GeometryOperator& op)
{
  if(dynamic_cast<const PointTransform*>(&op))
  {
    return "transform";
  }
  if(dynamic_cast<const Translation*>(&op))
  {
    return "translate";
  }
  if(dynamic_cast<const Rotation*>(&op))
  {
    return "rotate";
  }
  if(dynamic_cast<const Scale*>(&op))
  {
    return "scale";
  }
  if(dynamic_cast<const UnitConverter*>(&op))
  {
    return "convert_units_to";
  }
  if(dynamic_cast<const SliceOperator*>(&op))
  {
    return "slice";
  }
  return "unknown";
}

void requireMatrixOperator(const std::shared_ptr<const GeometryOperator>& op, int index)
{
  if(std::dynamic_pointer_cast<const MatrixOperator>(op))
  {
    return;
  }

  const auto ordinal =
    index >= 0 ? axom::fmt::format("operator {}", index) : std::string {"operator"};
  throw KleeError({Path {"geometry/operators"},
                   axom::fmt::format("Cannot convert geometry to matrix: {} ({}) is non-affine. "
                                     "Use applyTransform(point) for per-point evaluation.",
                                     ordinal,
                                     operatorName(*op))});
}

Geometry::Point3D applyOperator(const std::shared_ptr<const GeometryOperator>& op,
                                const Geometry::Point3D& point)
{
  if(auto composite = std::dynamic_pointer_cast<const CompositeOperator>(op))
  {
    Geometry::Point3D result = point;
    for(const auto& child : composite->getOperators())
    {
      result = applyOperator(child, result);
    }
    return result;
  }

  if(auto matrixOp = std::dynamic_pointer_cast<const MatrixOperator>(op))
  {
    return applyMatrix(matrixOp->toMatrix(), point);
  }

  if(auto pointTransform = std::dynamic_pointer_cast<const PointTransform>(op))
  {
    return pointTransform->apply(point);
  }

  throw KleeError(
    {Path {"geometry/operators"},
     axom::fmt::format("Cannot apply unsupported geometry operator '{}'.", operatorName(*op))});
}
}  // namespace

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
  // Only certain formats actually use Conduit serialization. Check for non-affine
  // operators only for those formats that will actually serialize to m_geomInfo.
  const bool needsSerialization =
    (m_format == "blueprint-tets" || m_format == "tet3D" || m_format == "sphere3D" ||
     m_format == "cone3D" || m_format == "hex3D" || m_format == "plane3D" || m_format == "sor3D");

  if(needsSerialization && hasNonAffineOperators())
  {
    throw KleeError({Path {"geometry"},
                     "Cannot serialize geometry to Conduit: geometry contains non-affine operators "
                     "(e.g., Lua transform functions) that cannot be represented as pure data. "
                     "Retain the original input deck file (.lua) for full reproducibility."});
  }

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

numerics::Matrix<double> Geometry::getTransform() const
{
  const auto identity4x4 = numerics::Matrix<double>::identity(4);
  numerics::Matrix<double> transformation(identity4x4);
  if(m_operator)
  {
    auto composite = std::dynamic_pointer_cast<const CompositeOperator>(m_operator);
    if(composite)
    {
      // Concatenate the transformations

      // Why don't we multiply the matrices in CompositeOperator::addOperator()?
      // Why keep the matrices factored and multiply them here repeatedly?
      // Combining them would also avoid this if-else logic.  BTNG
      int operatorIndex = 0;
      for(auto op : composite->getOperators())
      {
        ++operatorIndex;
        requireMatrixOperator(op, operatorIndex);
        // Use visitor pattern to extract the affine matrix from supported operators
        AffineMatrixVisitor visitor;
        op->accept(visitor);
        if(!visitor.isValid())
        {
          throw KleeError({Path {"geometry/operators"},
                           axom::fmt::format("Cannot convert geometry to matrix: operator {} ({}) "
                                             "is not supported by matrix extraction.",
                                             operatorIndex,
                                             operatorName(*op))});
        }
        const auto& matrix = visitor.getMatrix();
        numerics::Matrix<double> res(identity4x4);
        numerics::matrix_multiply(matrix, transformation, res);
        transformation = res;
      }
    }
    else
    {
      requireMatrixOperator(m_operator, -1);
      AffineMatrixVisitor visitor;
      m_operator->accept(visitor);
      if(visitor.isValid())
      {
        transformation = visitor.getMatrix();
      }
      else
      {
        throw KleeError({Path {"geometry/operators"},
                         axom::fmt::format("Cannot convert geometry to matrix: operator ({}) is "
                                           "not supported by matrix extraction.",
                                           operatorName(*m_operator))});
      }
    }
  }
  return transformation;
}

bool Geometry::hasNonAffineOperators() const { return operatorSequenceHasNonAffine(m_operator); }

Geometry::Point3D Geometry::applyTransform(const Point3D& point) const
{
  if(!m_operator)
  {
    return point;
  }
  return applyOperator(m_operator, point);
}

}  // namespace klee
}  // namespace axom
