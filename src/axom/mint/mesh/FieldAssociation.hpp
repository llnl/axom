// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FIELDASSOCIATION_HPP_
#define MINT_FIELDASSOCIATION_HPP_

namespace axom
{
namespace mint
{
/*!
 * \brief Enumerates the number of supported associations of field variables
 *  with a corresponding mesh entity, e.g., node, cell-center, face-center, etc.
 */
enum FieldAssociation
{
  ANY_CENTERING = -1,  //!< wild-card used to indicate any centering
  NODE_CENTERED = 0,   //!< used for fields computed at mesh nodes
  CELL_CENTERED,       //!< used for fields computed at cell centers
  FACE_CENTERED,       //!< used for fields computed at face centers
  EDGE_CENTERED,       //!< used for fields computed at edge centers

  NUM_FIELD_ASSOCIATIONS  //!< max number of field associations
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_FIELDASSOCIATION_HPP_ */
