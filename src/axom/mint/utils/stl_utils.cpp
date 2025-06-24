// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint/utils/stl_utils.hpp"

#include "axom/core/execution/reductions.hpp"
#include "axom/mint/mesh/Mesh.hpp"             /* for Mesh base class */
#include "axom/mint/execution/interface.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <fstream>

namespace axom
{
namespace mint
{
namespace internal
{
} // end namespace internal

//------------------------------------------------------------------------------
int read_stl(const std::string& file, Mesh*& mesh)
{
#if 0
  SLIC_ERROR_IF(file.length() <= 0, "No SU2 file was supplied!");
  SLIC_ERROR_IF(mesh != nullptr, "supplied mesh pointer should be a nullptr");

  std::ifstream ifs(file.c_str());
  if(!ifs.is_open())
  {
    SLIC_WARNING("cannot read from file [" << file << "]");
    return -1;
  }

  // STEP 0: read the raw data
  int ndime = -1;
  int nelem = -1;
  int npoin = -1;

  double* points = nullptr;
  axom::IndexType* connectivity = nullptr;
  mint::CellType* cellTypes = nullptr;
  bool isMixed = false;

  read_data(ifs, ndime, nelem, npoin, isMixed, points, connectivity, cellTypes);
  SLIC_ERROR_IF(ndime < 2 || ndime > 3, "mesh dimension must be 2 or 3!");
  SLIC_ERROR_IF(nelem <= 0, "mesh has zero cells!");
  SLIC_ERROR_IF(npoin <= 0, "mesh has zero nodes!");

  SLIC_ASSERT(points != nullptr);
  SLIC_ASSERT(connectivity != nullptr);
  SLIC_ASSERT(cellTypes != nullptr);

  ifs.close();

  // STEP 1: construct a mint mesh object
  if(isMixed)
  {
    using MeshType = mint::UnstructuredMesh<mint::MIXED_SHAPE>;
    MeshType* m = new MeshType(ndime, npoin, nelem);

    for(int i = 0; i < npoin; ++i)
    {
      m->appendNodes(&points[i * ndime], 1);
    }

    for(int i = 0; i < nelem; ++i)
    {
      m->appendCell(&connectivity[i * mint::MAX_CELL_NODES], cellTypes[i]);
    }

    mesh = m;
  }
  else
  {
    using MeshType = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
    MeshType* m = new MeshType(ndime, cellTypes[0], npoin, nelem);

    for(int i = 0; i < npoin; ++i)
    {
      m->appendNodes(&points[i * ndime], 1);
    }

    for(int i = 0; i < nelem; ++i)
    {
      m->appendCell(&connectivity[i * mint::MAX_CELL_NODES]);
    }

    mesh = m;
  }

  axom::deallocate(points);
  axom::deallocate(connectivity);
  axom::deallocate(cellTypes);

  SLIC_ASSERT(mesh != nullptr);
#endif
  return 0;
}

//------------------------------------------------------------------------------
int write_stl(const mint::Mesh* mesh, const std::string& filename, bool binary)
{
  SLIC_ERROR_IF(mesh == nullptr, "mesh pointer is null!");
  SLIC_ERROR_IF(filename.length() <= 0, "STL filename is empty!");
  SLIC_ERROR_IF(mesh->getDimension() == 2 || mesh->getDimension() == 3, "Input mesh must is not 2D/3D.");

  std::ofstream out(filename.c_str(), binary ? std::ofstream::binary : std::ofstream::out);
  if(!out.is_open())
  {
    SLIC_WARNING("cannot write to file [" << filename << "]");
    return -1;
  }

  // Write header
  if(binary)
  {
    std::uint8_t header[80] = {};
    strcpy(reinterpret_cast<char *>(header), "STL Binary Header Written by Axom");
    out.write(reinterpret_cast<const char*>(header), 80);

    // Compute the number of triangles we'll write.
    axom::ReduceSum<axom::SEQ_EXEC, std::uint32_t> ntri_reduce(0);
    for_all_faces<axom::SEQ_EXEC, xargs::nodeids>(
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
    out << "solid axom\n";
  }

  // Write triangle data.
  std::ofstream *out_ptr = &out;
  if(mesh->getDimension() == 2 && !binary)
  {
    for_all_faces<axom::SEQ_EXEC, xargs::nodeids>(
      mesh,
      AXOM_LAMBDA(IndexType AXOM_UNUSED_PARAM(faceID), const IndexType* nodes, IndexType N) {
        // Iterate over the face like a triangle fan.
        double coords[3][2];
        mesh->getNode(nodes[0], coords[0]);
        const IndexType ntri = N - 2;
        for(IndexType ti = 0; ti < ntri; ti++)
        {
          mesh->getNode(nodes[ti + 1], coords[1]);
          mesh->getNode(nodes[ti + 2], coords[2]);

          *out_ptr << "facet normal 0. 0. 1.\n";
          *out_ptr << "    outer loop\n";
          *out_ptr << "        vertex " << coords[0][0] << " " << coords[0][1] << " 0.\n";
          *out_ptr << "        vertex " << coords[1][0] << " " << coords[1][1] << " 0.\n";
          *out_ptr << "        vertex " << coords[2][0] << " " << coords[2][1] << " 0.\n";
          *out_ptr << "    endloop\n";
          *out_ptr << "endfacet\n";
        }
      });
  }
  else
  {
    using VectorType = axom::primal::Vector<double, 3>;
    using float32 = float;
    using uint16 = unsigned short;

    for_all_faces<axom::SEQ_EXEC, xargs::nodeids>(
      mesh,
      AXOM_LAMBDA(IndexType AXOM_UNUSED_PARAM(faceID), const IndexType* nodes, IndexType N) {
        // Iterate over the face like a triangle fan.
        double coords[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
        mesh->getNode(nodes[0], coords[0]);
        const IndexType ntri = N - 2;
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
            *out_ptr << "facet normal " << N[0] << " " << N[1] << " " << N[2] << "\n";
            *out_ptr << "    *out_ptrer loop\n";
            *out_ptr << "        vertex " << coords[0][0] << " " << coords[0][1] << " " << coords[0][2] << "\n";
            *out_ptr << "        vertex " << coords[1][0] << " " << coords[1][1] << " " << coords[1][2] << "\n";
            *out_ptr << "        vertex " << coords[2][0] << " " << coords[2][1] << " " << coords[2][2] << "\n";
            *out_ptr << "    endloop\n";
            *out_ptr << "endfacet\n";
          }
        }
      });
  }

  if(!binary)
  {
    out << "endsolid axom\n";
  }
  out.close();

  return 0;
}

} // end namespace mint
} // end namespace axom
