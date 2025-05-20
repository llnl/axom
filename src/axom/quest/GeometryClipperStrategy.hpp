// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_GEOMETRYCLIPPERSTRATEGY_HPP
#define AXOM_QUEST_GEOMETRYCLIPPERSTRATEGY_HPP

#include "axom/config.hpp"

#include "axom/core/Array.hpp"
#include "axom/klee/Geometry.hpp"
#include "axom/quest/ShapeeMesh.hpp"
#include "axom/primal.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Strategy base class for clipping operations for specific
  geometry instances.

  Key methods to implement:  (Some combination of these is required.)

  -# @c labelInOut: Label whether the cells in a mesh is inside,
     outside or on the shape boundary.  If a cell cannot be
     determined, you can conservatively label it as on the boundary.

  -# @c getShapesAsTets: Build an array of tetrahedra to approximate
     the shape.

  -# @c getShapesAsOcts: Build an array of octahedra to approximate
     the shape.

  -# @c specializedClip: Use a fast clipping algorithm (if one is
     available) to clip the cells in a mesh.  Implementation should
     use special knowledge of the geometry.  One version of this
     method clips all cells in the mesh and the other clips only
     cells in a provided index list.

  Every method returns true if it fulfilled the request, or
  false if it was a no-op.

  Implementations of this strategy must provide either a
  @c specializedClip method or one of the @c getShapesAs...() methods.
  The former is prefered if the use of geometry-specific information
  can make it faster.  @c labelInOut is optional but if provided,
  it can improve performance by limiting the slower clipping steps
  to a subset of cells.
*/
class GeometryClipperStrategy
{
public:
  /*!
    @brief A type to denote whether something is inside,
    on or outside the boundary of a geometry.
  */
  using LabelType = char;
  //!@brief Denotes something inside a shape boundary.
  static constexpr LabelType LABEL_IN = 0;
  //!@brief Denotes something on a shape boundary.
  static constexpr LabelType LABEL_ON = 1;
  //!@brief Denotes something outside a shape boundary.
  static constexpr LabelType LABEL_OUT = 2;

  using HexahedronType = axom::primal::Hexahedron<double, 3>;
  using OctahedronType = axom::primal::Octahedron<double, 3>;
  using TetrahedronType = axom::primal::Tetrahedron<double, 3>;
  using SphereType = axom::primal::Sphere<double, 3>;
  using Plane3DType = axom::primal::Plane<double, 3>;
  using Point3DType = axom::primal::Point<double, 3>;
  using Vector3DType = axom::primal::Vector<double, 3>;
  using Segment3DType = axom::primal::Segment<double, 3>;

  using CircleType = axom::primal::Sphere<double, 2>;
  using Plane2DType = axom::primal::Plane<double, 2>;
  using Point2DType = axom::primal::Point<double, 2>;
  using Segment2DType = axom::primal::Segment<double, 2>;

  /*!
    @brief Construct a shape clipper

    @param [in] kGeom Describes the shape to place
      into the mesh.
  */
  GeometryClipperStrategy(const klee::Geometry& kGeom);

  /*!
    @brief Optional name for strategy.

    The base implementation returns "UNNAMED".
  */
  virtual const std::string& name() const;

  /*!
    @brief Free-form representation of geometry.

    The exact information depends on the implementation.
  */
  const conduit::Node info() const
  {
    return m_info;
  }

  //@{
  //!@name Geometry-specialized methods
  /*!
    @brief Label the cells in the mesh as inside, outside or
    both/undetermined, if possible.

    @param [in] shapeeMesh Blueprint mesh to shape into.
    @param [out] labels Output

    The output labels are used in optimizing the clipping algorithm.
    Subclasses should implementation this if it's cost-effective, and
    skip if it's not.  It's safe to label cells as on the boundary if
    it can't be possitively determined as inside or outside.

    @return Whether the operation was done.  (A false means
    not done.)

    The cell labels should be set to
    - @c labelIn if the cell is completely inside the shape,
    - @c labelOut if the cell is completely outside, and
    - @c labelOn if the cell is both inside and outside (or
      cannot be easily determined).

    If implemenation returns true, it should ensure these
    post-conditions hold:
    @post labels.size() == shapeeMesh.getCellCount()
    @post labels.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(labels);
    return false;
  }

  /*!
    @brief Clip with a fast geometry-specialized method if
    possible.

    @param [in] shapeeMesh Blueprint mesh to shape into.
    @param ovlap [out] Shape overlap volume of each cell
      in the shapee mesh.

    The default implementation has no specialized method,
    so it's a no-op and returns false.

    If this method returns false, then exactly one of the
    @c getShapesAs...() methods must be provided.

    @return True if clipping was done and false if a no-op.

    This method need not be implemented if labelInOut()
    returns true.

    If implemenation returns true, it should ensure these
    post-conditions hold:
    @post ovlap.size() == shapeeMesh.getCellCount()
    @post ovlap.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool specializedClip(quest::ShapeeMesh& shapeeMesh, axom::ArrayView<double> ovlap)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(ovlap);
    return false;
  }

  /*!
    @brief Clip with a fast geometry-specialized method if
    possible.

    @param [in] shapeeMesh Blueprint mesh to shape into.
    @param [out] ovlap Shape overlap volume of each cell
      in the shapee mesh.
    @param [in] cellIds Limit computation to these cell ids.

    The default implementation has no specialized method,
    so it's a no-op and returns false.

    If this method returns false, then exactly one of the
    shape discretization methods must be provided.

    @return True if clipping was done and false if a no-op.

    This method need not be implemented if labelInOut()
    returns false.

    If implemenation returns true, it should ensure these
    post-conditions hold:
    @post ovlap.size() == shapeeMesh.getCellCount()
    @post ovlap.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool specializedClip(quest::ShapeeMesh& shapeeMesh,
                               axom::ArrayView<double> ovlap,
                               const axom::ArrayView<IndexType>& cellIds)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(ovlap);
    AXOM_UNUSED_VAR(cellIds);
    return false;
  }

  /*!
    @brief Get the fully transformed geometry as discrete tetrahedra,
    or return false.
    @param [in] shapeeMesh Blueprint mesh to shape into.
    @param [out] tets Array of tetrahedra filling the space of the shape.

    All vertex coordinates close to zero should be snapped to zero.

    @return Whether the shape can be represented as tetrahedra.

    If implemenation returns true, it should ensure these
    post-conditions hold:
    @post tets.size() == shapeeMesh.getCellCount()
    @post tets.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool getGeometryAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)

  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(tets);
    return false;
  }

  /*!
    @brief Get the fully transformed geometry as discrete octahedra,
    or return false.
    @param [in] shapeeMesh Blueprint mesh to shape into.
    @param [out] octs Array of octahedra filling the space of the shape.

    All vertex coordinates close to zero should be snapped to zero.

    @return Whether the shape can be represented as octahedra.

    Post-conditions only apply if method returns true.
    @post octs.size() == shapeeMesh.getCellCount()
    @post octs.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool getGeometryAsOcts(quest::ShapeeMesh& shapeeMesh, axom::Array<OctahedronType>& octs)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(octs);
    return false;
  }

  // Note: in 2D, we should have a getGeometryAsSegments().
  //@}

protected:
  /*!
    @brief Free-form representation of the concrete object.

    The constructor initializes this as a deep copy of the source
    klee::Geometry hierarchy data.  Subclasses may use and change this
    data as needed.

    The base class owns this hierachy, but subclasses use it
    for their specific data.

    Use cases:
    - Construct object from Klee input in the form of hierarchy
      data.
    - Print the object.

    Most if not all data should be in host memory.
  */
  conduit::Node m_info;

  /*!
    @brief Transformation due to the GeometryOperator.

    This is a direct result of the klee::Geometry::getGeometryOperator().
  */
  numerics::Matrix<double> m_transMat;

private:
  /*!
    @brief Compute the transformation matrix of a GeometryOperator.

    TODO: The matrix is equivalent to the operator.  This code is
    duplicated in several classes.  Should it be moved to
    GeometryOperator?

    TODO: I've not implemented coordinate transformations for any
    subclasses.  When the transform matrix is not identity, their
    results are wrong until the transformations take place.  Transformations
    need to happen in getGeometryAsTets, getGeometryAsOcts, labelInOut
    and specializedClip.  Basically everything.
  */
  numerics::Matrix<double> computeTransformationMatrix(
    const std::shared_ptr<const axom::klee::GeometryOperator>& op) const;
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_GEOMETRYCLIPPERSTRATEGY_HPP
