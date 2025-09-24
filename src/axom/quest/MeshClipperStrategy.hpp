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

// Requires Conduit for storing hierarchy-form data.
#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Strategy base class for clipping operations for specific
 * geometry instances.

 * Key methods to implement:  (Some combination of these is required.)

 * -# @c getBoundingBox2D or @c getBoundingBox3D: Axis-alligned
 *    bounding box for the geometry.

 * -# @c labelCellsInOut: Label whether the cells in a mesh is inside,
 *    outside or on the shape boundary.  If a cell cannot be
 *    determined, you can conservatively label it as on the boundary.

 * -# @c getShapesAsTets: Build an array of tetrahedra to approximate
 *    the shape.

 * -# @c getShapesAsOcts: Build an array of octahedra to approximate
 *    the shape.

 * -# @c specializedClipCells: Use a fast clipping algorithm (if one is
 *    available) to clip the cells in a mesh.  Implementation should
 *    use special knowledge of the geometry.  One version of this
 *    method clips all cells in the mesh and the other clips only
 *    cells in a provided index list.  The latter works in
 *    conjunction with @c labelCellsInOut.

 * Every method should return true if it fulfilled the request, or
 * false if it was a no-op.

 * Implementations of this strategy must provide either
 * - a @c specializedClipCells method or
 * - one of the @c getShapesAs...() methods.
 * The former is prefered if the use of geometry-specific information
 * can make it faster.  @c labelCellsInOut is optional but if provided,
 * it can improve performance by limiting the slower clipping steps
 * to a subset of cells.  @c getBoundingBox2D or @c getBoundingBox3D
 * can also improve performance by reducing computation.
*/
class MeshClipperStrategy
{
public:
  /*!
   * @brief A type to denote whether something is inside,
   * on or outside the boundary of a geometry.
  */
  using LabelType = char;
  //!@brief Denotes something inside a shape boundary.
  static constexpr LabelType LABEL_IN = 0;
  //!@brief Denotes something on a shape boundary.
  static constexpr LabelType LABEL_ON = 1;
  //!@brief Denotes something outside a shape boundary.
  static constexpr LabelType LABEL_OUT = 2;

  using BoundingBox3DType = axom::primal::BoundingBox<double, 3>;
  using Cone3DType = axom::primal::Cone<double, 3>;
  using HexahedronType = axom::primal::Hexahedron<double, 3>;
  using OctahedronType = axom::primal::Octahedron<double, 3>;
  using Plane3DType = axom::primal::Plane<double, 3>;
  using Point3DType = axom::primal::Point<double, 3>;
  using Ray3DType = axom::primal::Ray<double, 3>;
  using Segment3DType = axom::primal::Segment<double, 3>;
  using SphereType = axom::primal::Sphere<double, 3>;
  using TetrahedronType = axom::primal::Tetrahedron<double, 3>;
  using Triangle3DType = axom::primal::Triangle<double, 3>;
  using Vector3DType = axom::primal::Vector<double, 3>;

  using BoundingBox2DType = axom::primal::BoundingBox<double, 2>;
  using Point2DType = axom::primal::Point<double, 2>;
  using CircleType = axom::primal::Sphere<double, 2>;
  using Plane2DType = axom::primal::Plane<double, 2>;
  using Ray2DType = axom::primal::Ray<double, 2>;
  using Segment2DType = axom::primal::Segment<double, 2>;

  //!@brief Number of tetrahedra per hexahedron decomposes into
  // @internal We could use a more efficient 18-tet decomposition in the future.
  static constexpr axom::IndexType TETS_PER_HEXAHEDRON = HexahedronType::NUM_TRIANGULATE;

  /*!
   * @brief Construct a strategy for the given klee::Geometry object.
   *
   * @param [in] kGeom Describes the shape to place
   *   into the mesh.
  */
  MeshClipperStrategy(const klee::Geometry& kGeom);

  /*!
   * @brief Optional name for strategy.
   *
   * The base implementation returns "UNNAMED".
  */
  virtual const std::string& name() const;

  /*!
   * @brief Free-form information on the geometry.
   *
   * The exact information is provided by the klee::Geometry
   * and should be sufficient to define the geometry.
  */
  const conduit::Node& info() const { return m_info; }

  //@{
  //!@name Geometry-specialized methods

  /*!
   * @brief Get the 2D axis-aligned bounding box for the geometry,
   * if it's applicable and available.
  */
  virtual const axom::primal::BoundingBox<double, 2>& getBoundingBox2D() const;

  /*!
   * @brief Get the 3D axis-alligned bounding box for the geometry,
   * if it's applicable and available.
  */
  virtual const axom::primal::BoundingBox<double, 3>& getBoundingBox3D() const;

  /*!
   * @brief Label the cells in the mesh as inside, outside or
   * both/undetermined, if possible.
   *
   * @param [in] shapeeMesh Mesh to shape into.
   * @param [out] labels Output
   *
   * The cell labels should be set to
   * - @c labelIn if the cell is completely inside the shape,
   * - @c labelOut if the cell is completely outside, and
   * - @c labelOn if the cell is both inside and outside (or
   *   cannot be easily determined).
   *
   * The output labels are used in optimizing the clipping algorithm.
   * Subclasses should implementation this if it's cost-effective, and
   * skip if it's not.  It's safe to label cells as on the boundary if
   * it can't be possitively determined as inside or outside.
   *
   * @return Whether the operation was done.  (A false means
   * not done.)
   *
   * If implemenation returns true, it should ensure these
   * post-conditions hold:
   * @post labels.size() == shapeeMesh.getCellCount()
   * @post labels.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool labelCellsInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& cellLabels)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(cellLabels);
    return false;
  }

  /*!
   * @brief Label the tetrahedra in certain cells, if possible.
   *
   * @param [in] shapeeMesh Blueprint mesh to shape into.
   * @param [in] cellIds Indices of cells whose constituent
   *   tets should be labeled.
   * @param [out] tetLabels Output
   *
   * Only the cells labeled as ON the boundary are subjected to
   * this labeling.
   *
   * Tet indices refer to the @c shapeeMesh.getCellsAsTets() array.
   *
   * If implemenation returns true, it should ensure these
   * post-conditions hold:
   * @post tetLabels.size() == TETS_PER_HEXAHEDRON * cellIds.size()
   * @post labels.getAllocatorID() == shapeeMesh.getAllocatorId()
   * @post \c tetLabels should have \c TETS_PER_HEXAHEDRON labels
   * for each index in \c cellIds.
  */
  virtual bool labelTetsInOut(quest::ShapeeMesh& shapeeMesh,
                              axom::ArrayView<const axom::IndexType> cellIds,
                              axom::Array<LabelType>& tetLabels)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(cellIds);
    AXOM_UNUSED_VAR(tetLabels);
    return false;
  }

  /*!
   * @brief Clip with a fast geometry-specialized method if
   * possible.
   *
   * @param [in] shapeeMesh Blueprint mesh to shape into.
   * @param ovlap [out] Shape overlap volume of each cell
   *   in the shapee mesh.  It's initialized to zeros.
   *
   * The default implementation has no specialized method,
   * so it's a no-op and returns false.
   *
   * If this method returns false, then exactly one of the
   * @c getShapesAs...() methods must be provided.
   *
   * @return True if clipping was done and false if a no-op.
   *
   * This method need not be implemented if labelCellsInOut()
   * returns true.
   *
   * If implemenation returns true, it should ensure these
   * post-conditions hold:
   * @post ovlap.size() == shapeeMesh.getCellCount()
   * @post ovlap.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool specializedClipCells(quest::ShapeeMesh& shapeeMesh, axom::ArrayView<double> ovlap)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(ovlap);
    return false;
  }

  /*!
   * @brief Clip with a fast geometry-specialized method if
   * possible.
   *
   * @param [in] shapeeMesh Blueprint mesh to shape into.
   * @param [out] ovlap Shape overlap volume of each cell
   *   in the shapee mesh, initialized to the cell volumes
   *   for cell inside the shape and zero for other cells.
   * @param [in] cellIds Limit computation to these cell ids.
   *
   * The default implementation has no specialized method,
   * so it's a no-op and returns false.
   *
   * If this method returns false, then exactly one of the
   * shape discretization methods must be provided.
   *
   * @return True if clipping was done and false if a no-op.
   *
   * This method need not be implemented if labelCellsInOut()
   * returns false.
   *
   * @pre @c ovlap is pre-initialized for the implementation
   * to add or subtract partial volumes to individual cells.
   *
   * If implemenation returns true, it should ensure these
   * post-conditions hold:
   * @post ovlap.size() == shapeeMesh.getCellCount()
   * @post ovlap.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool specializedClipCells(quest::ShapeeMesh& shapeeMesh,
                                    axom::ArrayView<double> ovlap,
                                    const axom::ArrayView<IndexType>& cellIds)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(ovlap);
    AXOM_UNUSED_VAR(cellIds);
    return false;
  }

  /*!
   * Clip the tets listed in tetIds.
   * @param [in] shapeeMesh Blueprint mesh to shape into.
   * @param [out] ovlap Shape overlap volume of each cell
   *   in the shapee mesh, initialized to the cell volumes
   *   for cell inside the shape and zero for other cells.
   * @param [in] tetIds Indices of tets to clip, referring to the
   * shapeeMesh.getCellsAsTets() array.  tetIds[i] is the
   * \c (tetIds[i]%TETS_PER_HEXAHEDRON)-th tetrahedron of cell
   * \c tetIds[i]/TETS_PER_HEXAHEDRON.  Its overlap volume should be added
   * to that cell.
   */
  virtual bool specializedClipTets(quest::ShapeeMesh& shapeeMesh,
                                   axom::ArrayView<double> ovlap,
                                   const axom::ArrayView<IndexType>& tetIds)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(ovlap);
    AXOM_UNUSED_VAR(tetIds);
    return false;
  }

  /*!
   * @brief Get the geometry as discrete tetrahedra, or return false.
   *
   * @param [in] shapeeMesh Blueprint mesh to shape into.
   * @param [out] tets Array of tetrahedra filling the space of the shape,
   * fully transformed.
   *
   * All vertex coordinates close to zero should be snapped to zero.
   *
   * @return Whether the shape can be represented as tetrahedra.
   *
   * If implemenation returns true, it should ensure these
   * post-conditions hold:
   * @post tets.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool getGeometryAsTets(quest::ShapeeMesh& shapeeMesh,
                                 axom::Array<TetrahedronType>& tets)

  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(tets);
    return false;
  }

  /*!
   * @brief Get the geometry as discrete octahedra, or return false.
   *
   * @param [in] shapeeMesh Blueprint mesh to shape into.
   * @param [out] octs Array of octahedra filling the space of the shape,
   * fully transformed.
   *
   * All vertex coordinates close to zero should be snapped to zero.
   *
   * @return Whether the shape can be represented as octahedra.
   *
   * If implemenation returns true, it should ensure these
   * post-conditions hold:
   * @post octs.getAllocatorID() == shapeeMesh.getAllocatorId()
   */
  virtual bool getGeometryAsOcts(quest::ShapeeMesh& shapeeMesh,
                                 axom::Array<OctahedronType>& octs)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(octs);
    return false;
  }

  //@}

protected:
  /*!
   * @brief Free-form representation of the concrete object.
   *
   * The constructor initializes this as a deep copy of the source
   * klee::Geometry hierarchy data.  Subclasses may use and change this
   * data as needed.
   */
  conduit::Node m_info;

  /*!
   * @brief External transformation due to the GeometryOperator.
   *
   * This is a direct result of the klee::Geometry::getGeometryOperator().
   * Not to be confused with any geometry's internal transformation
   * (such as a cylinder's orientation and a sphere's center translation),
   * which apply before m_extTrans.
  */
  numerics::Matrix<double> m_extTrans;

private:
  /*!
    @brief Compute the transformation matrix of a GeometryOperator.
  */
  numerics::Matrix<double> computeTransformationMatrix(
    const std::shared_ptr<const axom::klee::GeometryOperator>& op) const;
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_GEOMETRYCLIPPERSTRATEGY_HPP
