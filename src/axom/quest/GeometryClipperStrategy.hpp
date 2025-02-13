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
  @brief Strategy base class for geometry-specific operations
  in clipping.

  Key methods to implement.  Some combination of these is required:

  -# getShapesAsTets: Build an array of tetrahedra from the
     mesh cells.

  -# getShapesAsOcts: Build an array of octahedra from the
     mesh cells.

  -# specializedClip: Use a fast clipping algorithm to clip
     the cells in a mesh.  Implementation may use special
     knowledge of the geometry.  One version of this method
     clips all cells in the mesh and the other clips only
     specified cells.

  -# labelInOut: Label whether the cells in a mesh is inside,
     outside or on the shape boundary.  If a cell cannot be
     determined, you can conservatively label it as on the boundary.

  Every method returns true if it fulfilled the request, or
  false if it was a no-op.

  Subclass must implement either a @c specializedClip method or one of
  the @c getShapesAs...() methods.  The former is prefered if the use
  of geometry-specific information can make it faster.  @c labelInOut
  is optional but if provided, it can improve performance by limiting
  the slower clipping steps to a subset of cells.
*/
class GeometryClipperStrategy
{
public:
  /*!
    @brief A type with possible values of 0, 1 or 2, denoting whether
    something is inside, on or outside the boundary of a geometry.
  */
  using LabelType = char;

  /*!
    @brief Construct a shape clipper

    @param [in] kGeom Describes the shape to place
      into the mesh.

    @c bpMesh must be an unstructured hex mesh.
    That is the only type currently supported.
  */
  GeometryClipperStrategy(const klee::Geometry& kGeom)
  {
    AXOM_UNUSED_VAR(kGeom);
  }

  //!@brief Optional name for strategy.
  virtual std::string name() const { return "UNNAMED"; }

  //@{
  //!@name Interface for geometry-specialized implementations
  /*!
    @brief Label the cells in the mesh as inside, outside or
    both/undetermined, if possible.

    @param [in/out] shapeeMesh Blueprint mesh to shape into.
    @param [out] labels

    The output labels are used in optimizing the clipping algorithm.
    Subclasses should implementation this if it's cost-effective, and
    skip if it's not.  It's safe to label cells as on the boundary if
    it can't be possitively determined as inside or outside.

    @return Whether the operation was done.  (A false means
    not done.)

    The cell labels should be set to
    - 1 if the cell is completely inside the shape,
    - 2 if the cell is completely outside, and
    - 3 if the cell is both inside and outside (or
      cannot be easily determined).

    Post-conditions only apply if method returns true.
    @post labels.size() == shapeeMesh.getCellCount()
    @post labels.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool labelInOut(quest::ShapeeMesh& shapeeMesh,
                          axom::Array<LabelType>& labels)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(labels);
    return false;
  }

  /*!
    @brief Clip with a fast geometry-specialized method if
    possible.

    @param [in/out] shapeeMesh Blueprint mesh to shape into.
    @param ovlap [out] Shape overlap volume of each cell
      in the shapee mesh.

    The default implementation has no specialized method,
    so it's a no-op and returns false.

    If this method returns false, then exactly one of the
    shape discretization methods must be provided.

    @return True if clipping was done and false if a no-op.

    Post-conditions only apply if method returns true.
    @post ovlap.size() == shapeeMesh.getCellCount()
    @post ovlap.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool specializedClip(quest::ShapeeMesh& shapeeMesh,
                               axom::ArrayView<double> ovlap)
  {
    AXOM_UNUSED_VAR(shapeeMesh);
    AXOM_UNUSED_VAR(ovlap);
    return false;
  }

  /*!
    @brief Clip with a fast geometry-specialized method if
    possible.

    @param [in/out] shapeeMesh Blueprint mesh to shape into.
    @param [out] ovlap Shape overlap volume of each cell
      in the shapee mesh.
    @param [in] cellIds Limit computation to these cell ids.

    The default implementation has no specialized method,
    so it's a no-op and returns false.

    If this method returns false, then exactly one of the
    shape discretization methods must be provided.

    @return True if clipping was done and false if a no-op.

    Post-conditions only apply if method returns true.
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
    @brief Get the shape as discrete tetrahedra, or return false.
    @param [out] tets Array of tetrahedra filling the space of the shape.

    All vertex coordinates close to zero should be snapped to zero.

    @return Whether the shape can be represented as tetrahedra.

    Post-conditions only apply if method returns true.
    @post tets.size() == shapeeMesh.getCellCount()
    @post tets.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool getShapeAsTets(axom::Array<axom::primal::Tetrahedron<double, 3>>& tets)
  {
    AXOM_UNUSED_VAR(tets);
    return false;
  }

  /*!
    @brief Get the shape as discrete tetrahedra, or return false.
    @param [out] octs Array of octahedra filling the space of the shape.

    All vertex coordinates close to zero should be snapped to zero.

    @return Whether the shape can be represented as octahedra.

    Post-conditions only apply if method returns true.
    @post octs.size() == shapeeMesh.getCellCount()
    @post octs.getAllocatorID() == shapeeMesh.getAllocatorId()
  */
  virtual bool getShapeAsOcts(axom::Array<axom::primal::Octahedron<double, 3>>& octs)
  {
    AXOM_UNUSED_VAR(octs);
    return false;
  }

  // Note: in 2D, we should have a getShapeAsSegments().
  //@}
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_GEOMETRYCLIPPERSTRATEGY_HPP
