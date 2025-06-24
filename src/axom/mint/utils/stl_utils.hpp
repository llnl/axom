// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_UTILS_STL_UTILS_HPP_
#define MINT_UTILS_STL_UTILS_HPP_

#include <string>  // for std::string

namespace axom
{
namespace mint
{
// Forward Declarations
class Mesh;

/*!
 * \brief Reads an unstructured mesh from an STL Mesh file.
 *
 * \param [in] file path corresponding to STL mesh file.
 * \param [out] mesh pointer to a mesh object where the mesh will be loaded.
 *
 * \return status error code, zero on success
 *
 * \note The STL Mesh file format is documented here:
 *  https://en.wikipedia.org/wiki/STL_(file_format)
 *
 * \note Ownership of the mesh object is passed to the caller. Consequently,
 *  the caller is responsible for properly deallocating the mesh object that
 *  the return mesh pointer points to.
 *
 * \pre file.length() > 0
 * \pre mesh == nullptr
 * \post mesh->isUnstructured()==true
 */
int read_stl(const std::string& file, Mesh*& mesh);

/*!
 * \brief Writes an unstructured mesh to the specified file according to the
 *  STL Mesh file format.
 *
 * \param [in] mesh pointer to the mesh to write.
 * \param [in] filename path to the file where to write the mesh.
 * \param [in] binary Whether the file should be written in binary.
 *
 * \return status error code, zero on success
 *
 * \note The STL Mesh file format is documented here:
 *  https://en.wikipedia.org/wiki/STL_(file_format)
 *
 * \pre mesh != nullptr
 * \pre mesh->isUnstructured()==true
 * \pre file.length() > 0
 */
int write_stl(const mint::Mesh* mesh, const std::string& filename, bool binary = true);

} /* namespace mint */

} /* namespace axom */

#endif /* MINT_UTILS_STL_UTILS_HPP_ */
