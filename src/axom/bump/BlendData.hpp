// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_BLEND_DATA_HPP_
#define AXOM_BUMP_BLEND_DATA_HPP_

#include "axom/core.hpp"

namespace axom
{
namespace bump
{
/*!
 * \brief This class contains views of blend data. BlendData contains the recipe
 *        for taking existing node values (e.g. coordinates, fields) and making
 *        new values that consist of node values that have been blended together
 *        in some way such as when it is necessary to introduce a new node along
 *        an edge or a new interior point.
 *
 *        The BlendData can be created and reused to make new coordsets and fields.
 *
 *        Field data are sampled using m_originalIdsView which is a compact list
 *        of the original node ids that we want to preserve without any blending.
 *
 *        This stream is followed by a second stream of data made using the field
 *        and the blend groups. Each blend group has m_blendGroupSizesView[i] elements,
 *        starts at m_blendGroupStartView[i] and uses values from the m_blendIdsView,
 *        m_blendCoeffView members to blend the data values.
 *
 * Example:
 * \verbatim
 *  Quad 0123    Make new points A,B,C using BlendGroup
 *
 *  3        2   A is 50% along the segment 0,1
 *  *--------*   B is in the center of quad 0123
 *  |        C   C is 75% along the segment 1,2
 *  |    B   |
 *  |        |   There are 3 blend groups (A,B,C)
 *  *----A---*
 *  0        1   blendGroupSizesView = {2, 4, 2};
 *               blendGroupStartView = {0, 2, 6};
 *               blendIdsView = {0, 1,       // A
 *                               0, 1, 2, 3, // B
 *                               1, 2};      // C
 *               blendCoeffView = {0.5, 0.5,               // A
 *                                 0.25, 0.25, 0.25, 0.25, // B
 *                                 0.25, 0.75};            // C
 * \endverbatim
 *   
 * \endcode
 */
struct BlendData
{
  using IndexView = axom::ArrayView<IndexType>;
  using CoeffView = axom::ArrayView<float>;

  IndexView m_originalIdsView {};      //!< Indices of original node ids.
  IndexView m_selectedIndicesView {};  //!< Indices of the selected blend groups.
  IndexView m_blendGroupSizesView {};  //!< Number of ids/weights in each blend group.
  IndexView m_blendGroupStartView {};  //!< Starting offset for a blend group in ids.
  IndexView m_blendIdsView {};         //!< Ids that make up blend groups.
  CoeffView m_blendCoeffView {};       //!< Weights that make up blend groups.
};

}  // end namespace bump
}  // end namespace axom

#endif
