// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/quest/io/STLWriter.hpp"

#include "axom/core/execution/reductions.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/execution/interface.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <fstream>

namespace axom
{
namespace quest
{

STLWriter::STLWriter(const std::string &filename, bool binary) : m_fileName(filename), m_binary(binary)
{
}

int STLWriter::write(const mint::Mesh* mesh)
{
  SLIC_ERROR_IF(mesh == nullptr, "mesh pointer is null!");
  SLIC_ERROR_IF(m_fileName.length() <= 0, "STL filename is empty!");
  SLIC_ERROR_IF(mesh->getDimension() == 2 || mesh->getDimension() == 3, "Input mesh must is not 2D/3D.");

  std::ofstream out(m_fileName.c_str(), m_binary ? std::ofstream::binary : std::ofstream::out);
  if(!out.is_open())
  {
    SLIC_WARNING("cannot write to file [" << m_fileName << "]");
    return -1;
  }

  // Write header
  if(m_binary)
  {
    std::uint8_t header[80] = {};
    strcpy(reinterpret_cast<char *>(header), "STL Binary Header Written by Axom");
    out.write(reinterpret_cast<const char*>(header), 80);

    // Compute the number of triangles we'll write.
    axom::ReduceSum<axom::SEQ_EXEC, std::uint32_t> ntri_reduce(0);
    axom::mint::for_all_faces<axom::SEQ_EXEC, axom::mint::xargs::nodeids>(
      mesh,
      AXOM_LAMBDA(IndexType AXOM_UNUSED_PARAM(faceID), const IndexType* AXOM_UNUSED_PARAM(nodes), IndexType N) {
        ntri_reduce += static_cast<std::uint32_t>(N - 2);
      });

    // Write number of triangles
    std::uint32_t ntri = ntri_reduce.get();
    out.write(reinterpret_cast<const char*>(&ntri), sizeof(std::uint32_t));
  }
  else
  {
    out << "solid triangles\n";
  }

  // For value capture.
  std::ofstream *out_ptr = &out;
  const bool binary = m_binary;

  // Write triangle data.
  if(mesh->getDimension() == 2 && !m_binary)
  {
    axom::mint::for_all_faces<axom::SEQ_EXEC, axom::mint::xargs::nodeids>(
      mesh,
      AXOM_LAMBDA(IndexType AXOM_UNUSED_PARAM(faceID), const IndexType* nodes, IndexType nnodes) {
        // Iterate over the face like a triangle fan.
        double coords[3][2];
        mesh->getNode(nodes[0], coords[0]);
        const IndexType ntri = nnodes - 2;
        for(IndexType ti = 0; ti < ntri; ti++)
        {
          mesh->getNode(nodes[ti + 1], coords[1]);
          mesh->getNode(nodes[ti + 2], coords[2]);

          *out_ptr << "\t facet normal 0. 0. 1.\n";
          *out_ptr << "\t\t outer loop\n";
          *out_ptr << "\t\t\t vertex " << coords[0][0] << " " << coords[0][1] << " 0.\n";
          *out_ptr << "\t\t\t vertex " << coords[1][0] << " " << coords[1][1] << " 0.\n";
          *out_ptr << "\t\t\t vertex " << coords[2][0] << " " << coords[2][1] << " 0.\n";
          *out_ptr << "\t\t endloop\n";
          *out_ptr << "\t endfacet\n";
        }
      });
  }
  else
  {
    using VectorType = axom::primal::Vector<double, 3>;
    using float32 = float;
    using uint16 = unsigned short;

    axom::mint::for_all_faces<axom::SEQ_EXEC, axom::mint::xargs::nodeids>(
      mesh,
      AXOM_LAMBDA(IndexType AXOM_UNUSED_PARAM(faceID), const IndexType* nodes, IndexType nnodes) {
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
          const VectorType N = VectorType::cross_product(B - A, C - A);

          if(binary)
          {
            float32 n32[3], coords32[3][3];
            for(int comp = 0; comp < 3; comp++)
            {
              n32[comp] = static_cast<float32>(N[comp]);
              coords32[0][comp] = static_cast<float32>(coords32[0][comp]);
              coords32[1][comp] = static_cast<float32>(coords32[1][comp]);
              coords32[2][comp] = static_cast<float32>(coords32[2][comp]);
            }
            const uint16 attr = 0x7fff;
            out_ptr->write(reinterpret_cast<const char *>(n32), 3 * sizeof(float32));
            out_ptr->write(reinterpret_cast<const char *>(coords32[0]), 3 * sizeof(float32));
            out_ptr->write(reinterpret_cast<const char *>(coords32[1]), 3 * sizeof(float32));
            out_ptr->write(reinterpret_cast<const char *>(coords32[2]), 3 * sizeof(float32));
            out_ptr->write(reinterpret_cast<const char *>(&attr), sizeof(uint16));
          }
          else
          {
            *out_ptr << "\t facet normal " << N[0] << " " << N[1] << " " << N[2] << "\n";
            *out_ptr << "\t\t outer loop\n";
            *out_ptr << "\t\t\t vertex " << coords[0][0] << " " << coords[0][1] << " " << coords[0][2] << "\n";
            *out_ptr << "\t\t\t vertex " << coords[1][0] << " " << coords[1][1] << " " << coords[1][2] << "\n";
            *out_ptr << "\t\t\t vertex " << coords[2][0] << " " << coords[2][1] << " " << coords[2][2] << "\n";
            *out_ptr << "\t\t endloop\n";
            *out_ptr << "\t endfacet\n";
          }
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

int write_stl(const mint::Mesh* mesh, const std::string &filename, bool binary)
{
   STLWriter w(filename, binary);
   return w.write(mesh);
}

}  // namespace quest
}  // namespace axom
