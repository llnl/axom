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
 * \brief A simple STL writer for a collection of triangular faces with
 *        limited support for planar polygons
 *
 * STL (STereoLithography) is a common file format for triangle meshes.
 * It encodes a "soup of triangles" by explicitly listing the coordinate
 * positions of the three vertices of each triangle.
 *
 * \note Any faces that contain more than 3 vertices will be treated
 *       like a triangle fan and will result in multiple triangles.
 *       No care is taken to prevent triangle overlaps for non-convex faces.
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
  STLWriter(const std::string& filename, bool binary = false);

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

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  /*!
   * \brief Compute the number of triangles produced for the input mesh.
   *
   * \return The number of triangles.
   */
  IndexType getNumberOfTriangles() const;

  /*!
   * \brief Determines whether the input mesh is topologically 2D.
   *
   * \return True if the mesh is topologically 2D, false otherwise.
   */
  bool isTopologically2D() const;

  const mint::Mesh* m_mesh {nullptr};
  std::string m_fileName {"output.stl"};
  bool m_binary {false};
};

/*!
 * \brief Function interface for STL writer.
 *
 * \param mesh The mesh whose faces will be written.
 * \param filename The name of the file to write.
 * \param binary Whether to write binary or ascii files.
 *
 * \return 0 on success; non-zero otherwise.
 */
int write_stl(const mint::Mesh* mesh, const std::string& filename, bool binary = false);

}  // namespace quest
}  // namespace axom

#endif  // QUEST_STLWRITER_HPP_
