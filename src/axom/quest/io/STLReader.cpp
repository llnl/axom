// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/STLReader.hpp"

// Axom includes
#include "axom/core/utilities/Utilities.hpp"
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/slic/interface/slic.hpp"

// C/C++ includes
#include <fstream>

namespace
{
const std::size_t BINARY_HEADER_SIZE = 80;  // bytes
const std::size_t BINARY_TRI_SIZE = 50;     // bytes
}  // namespace

//------------------------------------------------------------------------------
//      STLReader Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace quest
{
STLReader::STLReader() : m_fileName(""), m_num_nodes(0), m_num_faces(0) { }

//------------------------------------------------------------------------------
STLReader::~STLReader() { this->clear(); }

//------------------------------------------------------------------------------
void STLReader::clear()
{
  m_num_nodes = 0;
  m_num_faces = 0;
  m_nodes.clear();
}

//------------------------------------------------------------------------------
bool STLReader::isAsciiFormat() const
{
  // The binary format consists of
  //    a header of size BINARY_HEADER_SIZE==80 bytes
  //    followed by a four byte int encoding the number of triangles
  //    followed by the triangle data (BINARY_TRI_SIZE == 50 bytes per triangle)

  // Open the file
  std::ifstream ifs(m_fileName.c_str(), std::ios::in | std::ios::binary);

  if(!ifs.is_open())
  {
    /* short-circuit */
    SLIC_WARNING("Cannot open the provided STL file [" << m_fileName << "]");
    return false;
  }

  // Find out the file size
  ifs.seekg(0, ifs.end);
  std::int32_t fileSize = static_cast<std::int32_t>(ifs.tellg());

  const int totalHeaderSize = (BINARY_HEADER_SIZE + sizeof(std::int32_t));
  if(fileSize < totalHeaderSize)
  {
    return true;
  }

  // Find the number of triangles (if the file were binary)
  int numTris = 0;
  ifs.seekg(BINARY_HEADER_SIZE, ifs.beg);
  ifs.read((char*)&numTris, sizeof(std::int32_t));

  if(!utilities::isLittleEndian())
  {
    numTris = utilities::byteswap(numTris);
  }

  // Check if the size matches our expectation
  int expectedBinarySize = totalHeaderSize + (numTris * BINARY_TRI_SIZE);

  ifs.close();

  return (fileSize != expectedBinarySize);
}

//------------------------------------------------------------------------------
int STLReader::readAsciiSTL()
{
  std::ifstream ifs(m_fileName.c_str());

  if(!ifs.is_open())
  {
    SLIC_WARNING("Cannot open the provided STL file [" << m_fileName << "]");
    return (-1);
  }

  std::string junk;
  double x, y, z;

  // In an STL  file, we only care about the vertex positions
  // Vertices are strings of the form: "vertex v_x v_y v_z"
  while(true)
  {
    do
    {
      ifs >> junk;
    } while(ifs.good() && junk != "vertex");

    if(ifs.fail())
    {
      break;
    }

    ifs >> x >> y >> z;
    m_nodes.push_back(x);
    m_nodes.push_back(y);
    m_nodes.push_back(z);
  }

  // Set the number of nodes and faces
  m_num_nodes = static_cast<axom::IndexType>(m_nodes.size()) / 3;
  m_num_faces = m_num_nodes / 3;

  ifs.close();
  return (0);
}

//------------------------------------------------------------------------------
int STLReader::readBinarySTL()
{
  // Binary STL format consists of
  //    an 80 byte header (BINARY_HEADER_SIZE)
  //    followed by a 32 bit int encoding the number of faces
  //    followed by the triangles, each of which is 50 bytes (BINARY_TRI_SIZE)

  // A local union data structure for triangles in a binary STL
  union BinarySTLTri
  {
    std::int8_t raw[BINARY_TRI_SIZE];

    struct
    {
      float normal[3];
      float vert[9];
      std::uint16_t attr;
    } data;

  } tri;

  bool const isLittleEndian = axom::utilities::isLittleEndian();

  // Open binary file, skip the header
  std::ifstream ifs(m_fileName.c_str(), std::ios::in | std::ios::binary);
  if(!ifs.is_open())
  {
    SLIC_WARNING("Cannot open the provided STL file [" << m_fileName << "]");
    return (-1);
  }

  ifs.seekg(BINARY_HEADER_SIZE);

  // read the num faces and reserve room for the vertex positions
  ifs.read((char*)&m_num_faces, sizeof(std::int32_t));

  if(!isLittleEndian)
  {
    m_num_faces = utilities::byteswap(m_num_faces);
  }

  m_num_nodes = m_num_faces * 3;
  m_nodes.reserve(m_num_nodes * 3);

  // Read the triangles. Cast to doubles and ignore normals and attributes
  for(axom::IndexType i = 0; i < m_num_faces; ++i)
  {
    ifs.read((char*)tri.raw, BINARY_TRI_SIZE);

    for(int j = 0; j < 9; ++j)
    {
      float coord = isLittleEndian ? tri.data.vert[j] : utilities::byteswap(tri.data.vert[j]);

      m_nodes.push_back(static_cast<double>(coord));
    }
  }

  ifs.close();

  return (0);
}

//------------------------------------------------------------------------------
int STLReader::read()
{
  if(m_fileName.empty())
  {
    return (-1);
  }

  // Clear internal data, check the format and load the data
  this->clear();

  int rc = (isAsciiFormat()) ? readAsciiSTL() : readBinarySTL();
  return (rc);
}

//------------------------------------------------------------------------------
void STLReader::getMesh(axom::mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh)
{
  /* Sanity checks */
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(static_cast<axom::IndexType>(m_nodes.size()) != 3 * m_num_nodes,
                "nodes vector size doesn't match expected size!");
  SLIC_ERROR_IF(mesh->getDimension() != 3, "STL reader expects a 3D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::TRIANGLE, "STL reader expects a triangle mesh!");

  // pre-allocate space to store the mesh
  if(!mesh->isExternal())
  {
    mesh->resize(m_num_nodes, m_num_faces);
  }

  SLIC_ERROR_IF(mesh->getNumberOfNodes() != m_num_nodes,
                "mesh number of nodes does not match the number of nodes in the STL file!");
  SLIC_ERROR_IF(mesh->getNumberOfCells() != m_num_faces,
                "mesh number of cells does not match number of triangles in the STL file!");

  double* x = mesh->getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh->getCoordinateArray(mint::Y_COORDINATE);
  double* z = mesh->getCoordinateArray(mint::Z_COORDINATE);

  // Load the vertices into the mesh
  for(axom::IndexType i = 0; i < m_num_nodes; ++i)
  {
    const axom::IndexType offset = i * 3;
    x[i] = m_nodes[offset];
    y[i] = m_nodes[offset + 1];
    z[i] = m_nodes[offset + 2];
  }

  // Load the triangles.  Note that the indices are implicitly defined.
  axom::IndexType* conn = mesh->getCellNodesArray();
  for(axom::IndexType i = 0; i < m_num_faces; ++i)
  {
    const axom::IndexType offset = i * 3;
    conn[offset] = offset;
    conn[offset + 1] = offset + 1;
    conn[offset + 2] = offset + 2;
  }
}

}  // end namespace quest
}  // end namespace axom
