// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_MESHTYPES_HPP_
#define MINT_MESHTYPES_HPP_

namespace axom
{
namespace mint
{
/*!
 * \brief Defines the basic mesh types supported by mint.
 */
enum MeshTypes
{
  UNDEFINED_MESH = -1,  ///< UNDEFINED_MESH

  UNSTRUCTURED_MESH,  ///< UNSTRUCTURED_MESH

  STRUCTURED_CURVILINEAR_MESH,  ///< STRUCTURED_MESH
  STRUCTURED_RECTILINEAR_MESH,  ///< RECTILINEAR_MESH
  STRUCTURED_UNIFORM_MESH,      ///< UNIFORM_MESH

  PARTICLE_MESH,  ///< PARTICLE_MESH

  NUM_MESH_TYPES  ///< NUM_MESH_TYPES
};

} /* namespace mint */
} /* namespace axom */
#endif /* MINT_MESHTYPE_HPP_ */
