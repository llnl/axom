// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_INTERSECT_BOUNDING_BOX_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_BOUNDING_BOX_IMPL_HPP_

#include "axom/core/Macros.hpp"

#include <type_traits>  // for std::is_floating_point< T >()

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Helper routine for AABB / AABB intersection test.
 *
 * \param [in] min1 the 1st AABB min coordinate along a direction
 * \param [in] max1 the 1st AABB max coordinate along a direction
 * \param [in] min2 the 2nd AABB min coordinate along a direction
 * \param [in] max2 the 2nd AABB max coordinate along a direction
 *
 * \return status true if the AABB's intersect, otherwise false.
 *
 * \note This routine is called by the intersect AABB/AABB methods for each
 *  spatial dimension.
 * \note This routine does not have a tolerance parameter for fuzzy
 *  intersection, but the AABBs can be scaled to achieve the same thing.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_bbox_bbox(const T& min1,
                                                 const T& max1,
                                                 const T& min2,
                                                 const T& max2)
{
  return max1 >= min2 && min1 <= max2;
}

/*!
 * \brief Checks if a specified bounding box intersects with the supplied
 *  bounding box.
 *
 * \param [in] xmin1 x-coordinate of the specified lower bounding box corner
 * \param [in] xmax1 x-coordinate of the specified upper bounding box corner
 * \param [in] ymin1 y-coordinate of the specified lower bounding box corner
 * \param [in] ymax1 y-coordinate of the specified upper bounding box corner
 * \param [in] xmin2 x-coordinate of the supplied lower bounding box corner
 * \param [in] xmax2 x-coordinate of the supplied upper bounding box corner
 * \param [in] ymin2 y-coordinate of the supplied lower bounding box corner
 * \param [in] ymax2 y-coordinate of the supplied upper bounding box corner
 *
 * \return status true if the bounding boxes intersect, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_bounding_box(const T& xmin1,
                                                    const T& xmax1,
                                                    const T& ymin1,
                                                    const T& ymax1,
                                                    const T& xmin2,
                                                    const T& xmax2,
                                                    const T& ymin2,
                                                    const T& ymax2)
{
  bool status = true;
  status = status && intersect_bbox_bbox(xmin1, xmax1, xmin2, xmax2);
  status = status && intersect_bbox_bbox(ymin1, ymax1, ymin2, ymax2);
  return status;
}

/*!
 * \brief Checks if a specified bounding box intersects with the supplied
 *  bounding box.
 *
 * \param [in] xmin1 x-coordinate of the specified lower bounding box corner
 * \param [in] xmax1 x-coordinate of the specified upper bounding box corner
 * \param [in] ymin1 y-coordinate of the specified lower bounding box corner
 * \param [in] ymax1 y-coordinate of the specified upper bounding box corner
 * \param [in] zmin1 z-coordinate of the specified lower bounding box corner
 * \param [in] zmax1 z-coordinate of the specified upper bounding box corner
 * \param [in] xmin2 x-coordinate of the supplied lower bounding box corner
 * \param [in] xmax2 x-coordinate of the supplied upper bounding box corner
 * \param [in] ymin2 y-coordinate of the supplied lower bounding box corner
 * \param [in] ymax2 y-coordinate of the supplied upper bounding box corner
 * \param [in] zmin2 z-coordinate of the supplied lower bounding box corner
 * \param [in] zmax2 z-coordinate of the supplied upper bounding box corner
 *
 * \return status true if the bounding boxes intersect, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_bounding_box(const T& xmin1,
                                                    const T& xmax1,
                                                    const T& ymin1,
                                                    const T& ymax1,
                                                    const T& zmin1,
                                                    const T& zmax1,
                                                    const T& xmin2,
                                                    const T& xmax2,
                                                    const T& ymin2,
                                                    const T& ymax2,
                                                    const T& zmin2,
                                                    const T& zmax2)
{
  bool status = true;
  status = status && intersect_bbox_bbox(xmin1, xmax1, xmin2, xmax2);
  status = status && intersect_bbox_bbox(ymin1, ymax1, ymin2, ymax2);
  status = status && intersect_bbox_bbox(zmin1, zmax1, zmin2, zmax2);
  return status;
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_BOUNDING_BOX_IMPL_HPP_
