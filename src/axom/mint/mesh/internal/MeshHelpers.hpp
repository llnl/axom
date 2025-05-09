// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef MINT_MESH_HELPERS_HPP_
#define MINT_MESH_HELPERS_HPP_

#include "axom/core/Macros.hpp"  // for AXOM_UNUSED_PARAM
#include "axom/core/Types.hpp"   // for nullptr
#include "axom/core/Array.hpp"
#include "axom/mint/config.hpp"          // for mint compile-time type
#include "axom/mint/mesh/CellTypes.hpp"  // for CellType

#include <string>

namespace axom
{
namespace mint
{
class Mesh;  // forward declaration

namespace internal
{
inline int dim(const double* AXOM_UNUSED_PARAM(x), const double* y, const double* z)
{
  return ((z != nullptr) ? 3 : ((y != nullptr) ? 2 : 1));
}

std::string join_ints_into_string(int vcount, IndexType* values, char sep);

std::string make_face_key(int vcount, IndexType* values, char sep);

/*! \brief Record a Mesh's face-to-cell, cell-to-face, and face-to-node
 *         relations.
 *
 * \param [in] m pointer to a Mesh.
 * \param [out] facecount the number of unique faces of m's cells.
 * \param [out] f2c the relation between face f and its two incident cells
 *              with cellIDs at f2c[2*f] and f2c[2*f+1].
 * \param [out] c2f the relation between cell c and its n faces with faceIDs
 *              stored contiguously starting at c2f[c2foffsets[c]].
 * \param [out] c2n the relation between cell c and its n neighbors with
 *              cellIDs stored contiguously starting at c2n[c2foffsets[c]].
 * \param [out] c2foffsets the offset in c2f of the first face of each cell.
 * \param [out] f2n the relation between face f and its nodes stored
 *              contiguously starting at f2n[f2noffsets[f]].
 * \param [out] f2noffsets the offset in f2n of the first node of each face.
 * \param [out] f2ntypes the CellType of each face.
 *
 * \returns success true if each face has one or two incident cells.
 *
 * \note The seven output arrays f2c, c2f, c2n, c2foffsets, f2n, f2noffsets,
 * and f2ntypes are allocated in this routine if the routine is successful.
 * It is the caller's responsibility to free this memory.  If the routine
 * returns false, the output arrays are set to nullptr and facecount is set
 * to 0.
 *
 * This routine visits each of the cells of the mesh.  For each cell face, it
 * retrieves the face's nodes and joins the sorted node IDs to make a unique
 * hash key.  The incident cells are recorded in a list for each face's hash
 * key.  The final face-cell and cell-face relations are constructed from this
 * data structure.
 *
 * This routine is intended to be used in constructing an UnstructuredMesh's
 * face relations, though it will give correct results for any Mesh.
 */
bool initFaces(Mesh* m,
               IndexType& facecount,
               Array<IndexType>& f2c,
               Array<IndexType>& c2f,
               Array<IndexType>& c2n,
               Array<IndexType>& c2foffsets,
               Array<IndexType>& f2n,
               Array<IndexType>& f2noffsets,
               Array<CellType>& f2ntypes);

} /* namespace internal */
} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESH_HELPERS_HPP_ */
