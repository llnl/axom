// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_STLWRITER_HPP_
#define QUEST_STLWRITER_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/mint/mesh/Mesh.hpp"

#include <string>

namespace axom
{
namespace quest
{
/*!
 * \class STLWriter
 *
 * \brief A simple STL writer for Mint meshes.
 *
 * STL (STereoLithography) is a common file format for triangle meshes.
 * It encodes a "soup of triangles" by explicitly listing the coordinate
 * positions of the three vertices of each triangle.
 */
class STLWriter
{
public:
  /*!
   * \brief Constructor.
   */
  STLWriter() = default;

  /*!
   * \brief Constructor.
   *
   * \param filename The name of the file to write.
   * \param binary Whether or not to write a binary STL file.
   */
  STLWriter(const std::string &filename, bool binary = false);

  /*!
   * \brief Destructor.
   */
  ~STLWriter() = default;

  /*!
   * \brief Sets the name of the file to write.
   * \param [in] fileName the name of the file to write.
   */
  void setFileName(const std::string& fileName) { m_fileName = fileName; }

  /*!
   * \brief Sets whether to use binary output.
   * \param [in] binary True to write binary output; false to write ASCII.
   */
  void setBinary(bool binary) { m_binary = binary; }

  /*!
   * \brief Writes the mesh into an STL file.
   * \param [in] mesh pointer to the mesh to write.
   * \pre path to input file has been set by calling `setFileName()`
   * \return status set to zero on success; set to a non-zero value otherwise.
   */
  int write(const mint::Mesh* mesh);

protected:
  std::string m_fileName{"output.stl"};
  bool m_binary{false};
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_STLWRITER_HPP_
