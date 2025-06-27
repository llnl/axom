// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/quest/io/STLWriter.hpp"

#include "axom/core/execution/reductions.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/execution/interface.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <cstring>
#include <fstream>

namespace axom
{
namespace quest
{
namespace internal
{

/*!
 * \brief Write a triangle to an output stream.
 *
 * \param out The stream to use for writing.
 * \param binary Whether to write the triangle in binary format. The file must
 *               have been opened in binary mode.
 * \param coords The triangle coordinates.
 * \param N The triangle normal.
 */
template <typename NormalType>
void writeTriangle(std::ofstream &out, bool binary, double coords[3][3], const NormalType &N)
{
  if(binary)
  {
    using float32 = float;
    float32 n32[3], coords32[3][3];
    for(int comp = 0; comp < 3; comp++)
    {
      n32[comp] = static_cast<float32>(N[comp]);
      coords32[0][comp] = static_cast<float32>(coords[0][comp]);
      coords32[1][comp] = static_cast<float32>(coords[1][comp]);
      coords32[2][comp] = static_cast<float32>(coords[2][comp]);
    }
    // The attribute is sometimes used as colors. Set bits to white.
    // See https://en.wikipedia.org/wiki/STL_(file_format).
    const std::uint16_t attr = 0x7fff;
    out.write(reinterpret_cast<const char *>(n32), 3 * sizeof(float32));
    out.write(reinterpret_cast<const char *>(coords32[0]), 3 * sizeof(float32));
    out.write(reinterpret_cast<const char *>(coords32[1]), 3 * sizeof(float32));
    out.write(reinterpret_cast<const char *>(coords32[2]), 3 * sizeof(float32));
    out.write(reinterpret_cast<const char *>(&attr), sizeof(std::uint16_t));
  }
  else
  {
    out << "\t facet normal " << N[0] << " " << N[1] << " " << N[2] << "\n";
    out << "\t\t outer loop\n";
    out << "\t\t\t vertex " << coords[0][0] << " " << coords[0][1] << " " << coords[0][2] << "\n";
    out << "\t\t\t vertex " << coords[1][0] << " " << coords[1][1] << " " << coords[1][2] << "\n";
    out << "\t\t\t vertex " << coords[2][0] << " " << coords[2][1] << " " << coords[2][2] << "\n";
    out << "\t\t endloop\n";
    out << "\t endfacet\n";
  }
}

}  // end namespace internal

STLWriter::STLWriter(const std::string &filename, bool binary)
  : m_mesh(nullptr)
  , m_fileName(filename)
  , m_binary(binary)
{ }

bool STLWriter::isTopologically2D() const
{
  SLIC_ERROR_IF(m_mesh == nullptr, "mesh pointer is null!");

  bool is2D = false;
  if(m_mesh->getDimension() == 2)
  {
    is2D = true;
  }
  else if(m_mesh->getDimension() == 3)
  {
    // Heuristic for determining whether a m_mesh seems 2D independent of the coordinates.
    if(m_mesh->getNumberOfFaces() == 0 && m_mesh->getNumberOfCells() > 0)
    {
      // 2D shapes in 3D don't seem to have faces in mint. Get the type of the first cell.
      const auto ct = m_mesh->getCellType(0);
      is2D = (ct == mint::CellType::TRIANGLE || ct == mint::CellType::QUAD);
    }
  }
  return is2D;
}

IndexType STLWriter::getNumberOfTriangles() const
{
  SLIC_ERROR_IF(m_mesh == nullptr, "mesh pointer is null!");

  IndexType ntri = 0;
  if(isTopologically2D())
  {
    for(IndexType cellId = 0; cellId < m_mesh->getNumberOfCells(); cellId++)
    {
      ntri += (m_mesh->getNumberOfCellNodes(cellId) - 2);
    }
  }
  else if(m_mesh->getDimension() == 3)
  {
    axom::ReduceSum<axom::SEQ_EXEC, axom::IndexType> ntri_reduce(0);
    axom::mint::for_all_faces<axom::SEQ_EXEC, axom::mint::xargs::nodeids>(
      m_mesh,
      AXOM_LAMBDA(IndexType AXOM_UNUSED_PARAM(faceID),
                  const IndexType *AXOM_UNUSED_PARAM(nodes),
                  IndexType N) { ntri_reduce += (N - 2); });
    ntri = ntri_reduce.get();
  }
  return ntri;
}

int STLWriter::write(const mint::Mesh *mesh)
{
  using VectorType = axom::primal::Vector<double, 3>;

  SLIC_ERROR_IF(mesh == nullptr, "mesh pointer is null!");
  SLIC_ERROR_IF(m_fileName.length() <= 0, "STL filename is empty!");
  SLIC_ERROR_IF(mesh->getDimension() < 2 || mesh->getDimension() > 3, "Input mesh is not 2D/3D.");

  // Save mesh pointer
  m_mesh = mesh;

  std::ofstream out(m_fileName.c_str(),
                    m_binary ? (std::ofstream::out | std::ofstream::binary) : std::ofstream::out);
  if(!out.is_open())
  {
    SLIC_WARNING("cannot write to file [" << m_fileName << "]");
    return -1;
  }

  // Write header
  if(m_binary)
  {
    constexpr int STL_HEADER_SIZE = 80;
    std::uint8_t header[STL_HEADER_SIZE];
    // Fill with spaces
    memset(header, ' ', sizeof(std::uint8_t) * STL_HEADER_SIZE);
    // Copy in string (without terminator)
    const char *msg = "STL Binary File Written By Axom";
    memcpy(header, msg, strlen(msg));
    out.write(reinterpret_cast<const char *>(header), STL_HEADER_SIZE);

    // Write number of triangles
    std::uint32_t ntri = static_cast<std::uint32_t>(getNumberOfTriangles());
    out.write(reinterpret_cast<const char *>(&ntri), sizeof(std::uint32_t));
  }
  else
  {
    out << "solid triangles\n";
  }

  // Write triangle data.
  if(isTopologically2D())
  {
    // For meshes with cells that look 2D, we iterate over the cells.
    axom::Array<axom::IndexType> nodes;
    VectorType N = VectorType::make_vector(0., 0., 1.);
    double coords[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    for(axom::IndexType cellId = 0; cellId < m_mesh->getNumberOfCells(); cellId++)
    {
      // Get nodes for this cell.
      nodes.resize(m_mesh->getNumberOfCellNodes(cellId));
      const auto nnodes = m_mesh->getCellNodeIDs(cellId, nodes.data());

      // Iterate over the face like a triangle fan.
      m_mesh->getNode(nodes[0], coords[0]);
      const IndexType ntri = nnodes - 2;
      for(IndexType ti = 0; ti < ntri; ti++)
      {
        m_mesh->getNode(nodes[ti + 1], coords[1]);
        m_mesh->getNode(nodes[ti + 2], coords[2]);

        // If the cell is in 3D space, make a better normal.
        if(m_mesh->getDimension() == 3)
        {
          const VectorType A(coords[0], 3);
          const VectorType B(coords[1], 3);
          const VectorType C(coords[2], 3);
          N = VectorType::cross_product(B - A, C - A).unitVector();
        }

        internal::writeTriangle(out, m_binary, coords, N);
      }
    }
  }
  else
  {
    // For value capture.
    std::ofstream *out_ptr = &out;
    const bool binary = m_binary;

    axom::mint::for_all_faces<axom::SEQ_EXEC, axom::mint::xargs::nodeids>(
      m_mesh,
      AXOM_LAMBDA(IndexType AXOM_UNUSED_PARAM(faceID), const IndexType *nodes, IndexType nnodes) {
        // NOTE: Here in the lambda, we use "mesh" instead of "m_mesh" so we do
        //       not capture STLWriter's "this" pointer.

        // Iterate over the face like a triangle fan.
        double coords[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
        mesh->getNode(nodes[0], coords[0]);
        const IndexType ntri = nnodes - 2;
        for(IndexType ti = 0; ti < ntri; ti++)
        {
          mesh->getNode(nodes[ti + 1], coords[1]);
          mesh->getNode(nodes[ti + 2], coords[2]);

          // Compute facet normal.
          const VectorType A(coords[0], 3);
          const VectorType B(coords[1], 3);
          const VectorType C(coords[2], 3);
          const VectorType N = VectorType::cross_product(B - A, C - A).unitVector();

          internal::writeTriangle(*out_ptr, binary, coords, N);
        }
      });
  }

  if(!m_binary)
  {
    out << "endsolid triangles\n";
  }
  out.close();

  return 0;
}

int write_stl(const mint::Mesh *mesh, const std::string &filename, bool binary)
{
  STLWriter w(filename, binary);
  return w.write(mesh);
}

}  // namespace quest
}  // namespace axom
