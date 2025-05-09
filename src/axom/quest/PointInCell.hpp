// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_POINT_IN_CELL_HPP_
#define AXOM_QUEST_POINT_IN_CELL_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"

#include "axom/quest/detail/PointFinder.hpp"

namespace axom
{
namespace quest
{
/*!
 * \class PointInCellTraits
 * \brief A traits class for the mesh associated with a PointInCell query.
 *
 * \tparam mesh_tag A tag struct used to identify the mesh
 *
 * An implementation of \a PointInCellTraits on \a mesh_tag must define:
 * \arg A \a MeshType (e.g. mfem::Mesh or axom::mint::Mesh)
 * \arg An \a IndexType (e.g. a signed or unsigned integral type)
 * \arg An \a IndexType variable named \a NO_CELL
 *     (e.g. with value -1 for signed ints), indicating an invalid cell
 *     index in the mesh
 */
template <typename mesh_tag>
struct PointInCellTraits;

namespace detail
{
/*!
 * \class PointInCellMeshWrapper
 * \brief A wrapper class for accessing the mesh instance
 * associated with a PointInCell query
 *
 * \tparam mesh_tag A tag struct used to identify the mesh
 *
 * An implementation of a \a PointInCellMeshWrapper on \a mesh_tag
 * must expose the following API:
 * \arg template<DIM>
 *      void computeBoundingBox(
 *          double, const std::vector< BoundingBox<double, DIM>&,
 *          const BoundingBox<double, DIM>  &) const;
 * \arg void reconstructPoint(IndexType, const double*, double*) const;
 * \arg bool locatePointInCell(IndexType, const double*, double*) const;
 * \arg int numElements() const;
 * \arg int meshDimension() const;
 */
template <typename mesh_tag>
class PointInCellMeshWrapper;

}  // end namespace detail

/*!
 * \class PointInCell
 *
 * \brief A class to accelerate Point-In-Cell queries on a computational mesh.
 *
 * A point in cell query over a computational mesh determines if a given point
 * is contained in one of its cells.  If so, it returns the index of the cell
 * containing the point as well as the isoparametric coordinates of the point
 * within the cell.
 *
 * \tparam mesh_tag A tag type (e.g. an empty struct) for the underlying
 * computational mesh. There must also be corresponding template specializations
 * of \a PointInCellMeshWrapper<mesh_tag> and \a PointInCellTraits<mesh_tag> for
 * the provided mesh_tag in the \a axom::quest::detail namespace
 *
 * \note This class was designed to support point in cell queries against
 * 2D or 3D computational meshes. The queries to the mesh are wrapped in a
 * PointInCellMeshWrapper class templated on a \a mesh_tag type.
 * To extend this design, one must create a new mesh tag (e.g.. an empty struct
 * \a custom_mesh_tag) and add custom template specializations of
 * PointInCellMeshWrapper and PointInCellTraits for this tag in the
 * axom::quest namespace.
 *
 * \sa PointInCellMeshWrapper_mfem.hpp for a specialized implementation
 * for [mfem](http://mfem.org) meshes of arbitrary order.  It uses
 * the mesh_tag \a quest_point_in_cell_mfem_tag
 */
template <typename mesh_tag, typename ExecSpace = axom::SEQ_EXEC>
class PointInCell
{
public:
  using Point2DType = primal::Point<double, 2>;
  using Point3DType = primal::Point<double, 3>;

  using MeshTraits = PointInCellTraits<mesh_tag>;
  using MeshType = typename MeshTraits::MeshType;
  using IndexType = typename MeshTraits::IndexType;

  using MeshWrapperType = detail::PointInCellMeshWrapper<mesh_tag>;
  using PointFinder2D = detail::PointFinder<2, mesh_tag, ExecSpace>;
  using PointFinder3D = detail::PointFinder<3, mesh_tag, ExecSpace>;

  /*!
   * Construct a point in cell query structure over a computational mesh
   *
   * \param [in] mesh A pointer to the computational mesh
   * \param [in] resolution Grid resolution for the spatial index. Default: NULL
   * \param [in] bboxTolerance A tolerance factor by which to expand
   * the bounding boxes. Default: 1e-8
   * \param [in] allocatorId Currently unused. Default value is based on the
   *  allocator ID set for the specified execution space.
   *
   * \note The bboxTolerance should be a small positive number.  It helps avoid
   * numerical issues in the bounding box containment queries by slightly
   * expanding the cell bounding boxes.
   *
   * \note If the resolution is not provided, a heuristic based on the number
   * of cells in the mesh is used to set the resolution.
   * \sa ImplicitGrid
   *
   * \pre \a mesh must not be a NULL pointer
   * \pre If resolution is not NULL, it must have space for at least
   * meshDimension() entries.
   */
  PointInCell(MeshType* mesh,
              int* resolution = nullptr,
              double bboxTolerance = 1e-8,
              int allocatorID = axom::execution_space<ExecSpace>::allocatorID())
    : m_meshWrapper(mesh)
    , m_pointFinder2D(nullptr)
    , m_pointFinder3D(nullptr)
  {
    SLIC_ASSERT(mesh != nullptr);

    const double bboxScaleFactor = 1. + bboxTolerance;

    // Allocate a 2D or 3D PointFinder instance, depending on mesh dimension.
    // Note: Only one of these will be allocated in a PointInCell instance
    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      m_pointFinder2D = new PointFinder2D(&m_meshWrapper, resolution, bboxScaleFactor, allocatorID);
      break;
    case 3:
      m_pointFinder3D = new PointFinder3D(&m_meshWrapper, resolution, bboxScaleFactor, allocatorID);
      break;
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      break;
    }
  }

  /*! Destructor */
  ~PointInCell()
  {
    if(m_pointFinder2D != nullptr)
    {
      delete m_pointFinder2D;
      m_pointFinder2D = nullptr;
    }

    if(m_pointFinder3D != nullptr)
    {
      delete m_pointFinder3D;
      m_pointFinder3D = nullptr;
    }
  }

  /*!
   * Attempt to find the index of the mesh cell containing the given point.
   *
   * If found, and \a isopar is not NULL, \a isopar contains the isoparametric
   * coordinates the point within this cell.
   *
   * \param[in] pos The coordinates of the query point in space
   * \param[out] isopar The isoparametric coordinates of the query pt.
   * Only valid when a cell is found.
   *
   * \return The index of the mesh cell containing the query point.
   * If no cell is found, returns the special value \a MeshTraits::NO_CELL.
   * Otherwise, the result will be between 0 and N, where N is the number
   * of cells in the mesh.
   *
   * \pre \a pos is a non-null array with at least \a meshDimension() coords
   * \pre When not NULL, \a isopar has space for \a meshDimension() coords
   */
  IndexType locatePoint(const double* pos, double* isopar = nullptr) const
  {
    SLIC_ASSERT(pos != nullptr);

    IndexType cellIndex = MeshTraits::NO_CELL;

    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      cellIndex = m_pointFinder2D->locatePoint(pos, isopar);
      break;
    case 3:
      cellIndex = m_pointFinder3D->locatePoint(pos, isopar);
      break;
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      break;
    }

    return cellIndex;
  }

  void locatePoints(axom::ArrayView<const Point2DType> pts,
                    IndexType* outCellIds,
                    Point2DType* outIsopar = nullptr)
  {
    SLIC_ASSERT(pts.size() > 0);
    SLIC_ASSERT(pts.data() != nullptr);
    SLIC_ASSERT(outCellIds != nullptr);
    SLIC_ASSERT(m_pointFinder2D != nullptr);

    m_pointFinder2D->locatePoints(pts, outCellIds, outIsopar);
  }

  void locatePoints(axom::ArrayView<const Point3DType> pts,
                    IndexType* outCellIds,
                    Point3DType* outIsopar = nullptr) const
  {
    SLIC_ASSERT(pts.size() > 0);
    SLIC_ASSERT(pts.data() != nullptr);
    SLIC_ASSERT(outCellIds != nullptr);
    SLIC_ASSERT(m_pointFinder3D != nullptr);

    m_pointFinder3D->locatePoints(pts, outCellIds, outIsopar);
  }

  /*!
   *  Determine if a query point is located within a specified mesh cell
   *
   *  \param [in] cellIdx The index of the cell within the mesh
   *  \param [in] pos The coordinates of the query point in space
   *  \param [out] isopar The isoparametric coordinates of the point
   *  within cell \a cellIdx.  Only valid when the return value is true
   *
   *  \return True, if the query point is located within the specified cell
   *  \pre \a pos is not NULL and has \a meshDimension() entries
   *  \pre \a isopar is not NULL and has space for \a meshDimension() entries
   */
  bool locatePointInCell(IndexType cellIdx, const double* pos, double* isopar) const
  {
    return m_meshWrapper.locatePointInCell(cellIdx, pos, isopar);
  }

  /*!
   * Evaluate the position of a point within a mesh cell at the given
   * isoparametric coordinates.
   *
   * \param [in] cellIdx The index of the cell within the mesh
   * \param [in[ isopar The isoparametric coordinates at which to evaluate
   * \param [out] pos The computed coordinates of the evaluated point
   */
  void reconstructPoint(IndexType cellIdx, const double* isopar, double* pos) const
  {
    m_meshWrapper.reconstructPoint(cellIdx, isopar, pos);
  }

  /*! Returns the dimension of the mesh */
  int meshDimension() const { return m_meshWrapper.meshDimension(); }

  /*!
   * \brief Sets the print verbosity level for the point in cell query
   *
   * \param [in] level The verbosity level (increases with level)
   *  
   * This is useful for debugging the point in cell query
   * 
   *  For the mfem mesh wrapper, the valid options are: 
   *  - -1: never print (default)
   *  -  0: print only errors
   *  -  1: print the first and last iterations
   *  -  2: print every iteration
   *  -  3: print every iteration including point coordinates.
   */
  void setPrintLevel(int level) { m_meshWrapper.setPrintLevel(level); }

  /*!
   * \brief Sets the initial guess type for the element-based point in cell query
   *
   * \param [in] guessType The guess type
   *  
   *  For the mfem mesh wrapper, the valid options are: 
   *  - 0: Use the element center in reference space
   *  - 1: Use the closest physical node on a grid of points in physical space
   *  - 2: Use the closest reference node on a grid of points in reference space
   * 
   *  The grid size is controlled by setInitialGridOrder()
   */
  void setInitialGuessType(int guessType) { m_meshWrapper.setInitialGuessType(guessType); }

  /*!
   * \brief Sets the grid size for the initial guess in the element-based point in cell query
   *
   * \param [in] order The order for the grid size 
   *  
   *  For the mfem mesh wrapper, the number of points in each spatial direction 
   *  is given by `max(trans_order+order,0)+1`, where trans_order is the order 
   *  of the current element.
   * 
   *  \sa setInitialGuessType
   */
  void setInitialGridOrder(int order) { m_meshWrapper.setInitialGridOrder(order); }

  /*!
   * \brief Sets the solution strategy for the element-based point in cell query
   *
   * \param [in] type The strategy type
   *  
   *  For the mfem mesh wrapper, the valid options all use a Newton solve
   *  but differ in their handling of iterates that leave the reference element
   *  - 0: Allow the iterates to leave the reference element
   *  - 1: Project external iterates to the reference space boundary along their current line
   *  - 2: Project external iterates to the closest reference space boundary location
   */
  void setSolverProjectionType(int type) { m_meshWrapper.setSolverProjectionType(type); }

  /*!
   * \brief Gets the candidate IDs for the provided query point \a pt
   *
   * \param [in] pt The 2D query point.
   * \returns A std::vector of the candidates associated with \a pt
   * \note If the query point is 2D and the spatial index is 3D, we use 0 for the z-coordinate
   */
  std::vector<IndexType> getCandidatesForPt(axom::primal::Point<double, 2> const& pt)
  {
    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      return m_pointFinder2D->getCandidates(pt);
    case 3:
      return m_pointFinder3D->getCandidates(axom::primal::Point<double, 3>(pt.data(), 2));
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      return std::vector<IndexType> {};
    }
  }

  /*!
   * \brief Gets the candidate IDs for the provided query point \a pt
   *
   * \param [in] pt The 3D query point.
   * \returns A std::vector of the candidates associated with \a pt
   * \note If the query point is 3D and the spatial index is 2D, we omit the z-coordinate,
   * which effectively projects the point to the XY plane.
   */
  std::vector<IndexType> getCandidatesForPt(axom::primal::Point<double, 3> const& pt)
  {
    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      return m_pointFinder2D->getCandidates(axom::primal::Point<double, 2>(pt.data(), 2));
    case 3:
      return m_pointFinder3D->getCandidates(pt);
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      return std::vector<IndexType> {};
    }
  }

private:
  MeshWrapperType m_meshWrapper;

  PointFinder2D* m_pointFinder2D;
  PointFinder3D* m_pointFinder3D;
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_POINT_IN_CELL_HPP_
