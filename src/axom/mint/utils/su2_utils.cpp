// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint/utils/su2_utils.hpp"

#include "axom/mint/mesh/Mesh.hpp"             /* for Mesh base class */
#include "axom/mint/mesh/UnstructuredMesh.hpp" /* for UnstructuredMesh */
#include "axom/mint/mesh/CellTypes.hpp"

// C/C++ includes
#include <fstream>  // for std::ifstream

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
constexpr int SU2_LINE = 3;
constexpr int SU2_TRIANGLE = 5;
constexpr int SU2_QUAD = 9;
constexpr int SU2_TET = 10;
constexpr int SU2_HEX = 12;
constexpr int SU2_PRISM = 13;
constexpr int SU2_PYRAMID = 14;

//------------------------------------------------------------------------------
mint::CellType getMintCellType(int su2CellType)
{
  mint::CellType c = mint::UNDEFINED_CELL;

  switch(su2CellType)
  {
  case SU2_LINE:
    c = mint::SEGMENT;
    break;
  case SU2_TRIANGLE:
    c = mint::TRIANGLE;
    break;
  case SU2_QUAD:
    c = mint::QUAD;
    break;
  case SU2_TET:
    c = mint::TET;
    break;
  case SU2_HEX:
    c = mint::HEX;
    break;
  case SU2_PRISM:
    c = mint::PRISM;
    break;
  case SU2_PYRAMID:
    c = mint::PYRAMID;
    break;
  default:
    c = mint::UNDEFINED_CELL;
  }  // END switch

  return (c);
}

//------------------------------------------------------------------------------
int getSU2CellType(mint::CellType mint_type)
{
  int su2type = mint::cellTypeToInt(mint::UNDEFINED_CELL);

  switch(mint_type)
  {
  case mint::SEGMENT:
    su2type = SU2_LINE;
    break;
  case mint::TRIANGLE:
    su2type = SU2_TRIANGLE;
    break;
  case mint::QUAD:
    su2type = SU2_QUAD;
    break;
  case mint::TET:
    su2type = SU2_TET;
    break;
  case mint::HEX:
    su2type = SU2_HEX;
    break;
  case mint::PRISM:
    su2type = SU2_PRISM;
    break;
  case mint::PYRAMID:
    su2type = SU2_PYRAMID;
    break;
  default:
    su2type = mint::cellTypeToInt(mint::UNDEFINED_CELL);
  }  // END switch

  return (su2type);
}

//------------------------------------------------------------------------------
void read_points(double* points, int npoin, int ndime, std::ifstream& ifs)
{
  SLIC_ASSERT(points != nullptr);

  if(npoin <= 0)
  {
    return;
  }

  std::string line = "";
  for(int ipoint = 0; ipoint < npoin; ++ipoint)
  {
    const int offset = ipoint * ndime;
    for(int idim = 0; idim < ndime; ++idim)
    {
      ifs >> points[offset + idim];
    }

    std::getline(ifs, line);

  }  // END for all points
}

//------------------------------------------------------------------------------
void read_connectivity(axom::IndexType* connectivity,
                       mint::CellType* cellTypes,
                       int nelem,
                       bool& isMixed,
                       std::ifstream& ifs)
{
  SLIC_ASSERT(cellTypes != nullptr);
  SLIC_ASSERT(connectivity != nullptr);

  if(nelem <= 0)
  {
    return;
  }

  std::string line = "";
  int ctype;
  for(int icell = 0; icell < nelem; ++icell)
  {
    ifs >> ctype;
    mint::CellType c = getMintCellType(ctype);
    SLIC_ASSERT(c != mint::UNDEFINED_CELL);
    cellTypes[icell] = c;

    if(!isMixed)
    {
      isMixed = ((icell > 0) && (cellTypes[icell - 1] != c)) ? true : false;
    }

    const int offset = icell * mint::MAX_CELL_NODES;
    const int numNodes = mint::getCellInfo(c).num_nodes;
    for(int inode = 0; inode < numNodes; ++inode)
    {
      ifs >> connectivity[offset + inode];
    }

    std::getline(ifs, line);

  }  // END for all cells
}

//------------------------------------------------------------------------------
void read_data(std::ifstream& ifs,
               int& ndime,
               int& nelem,
               int& npoin,
               bool& isMixed,
               double*& points,
               axom::IndexType*& connectivity,
               mint::CellType*& cellTypes)
{
  // sanity checks
  SLIC_ASSERT(points == nullptr);
  SLIC_ASSERT(connectivity == nullptr);
  SLIC_ASSERT(cellTypes == nullptr);

  std::string line = "";
  while(std::getline(ifs, line))
  {
    if(line.substr(0, 1) == "%" || line.length() == 0)
    {
      // skip empty lines and comments
      continue;
    }

    if(line.substr(0, 6) == "NDIME=")
    {
      ndime = std::stoi(line.substr(6, line.length()));
      SLIC_ERROR_IF((ndime < 2) || (ndime > 3), "mesh dimension must be 2 or 3!");
    }
    else if(line.substr(0, 6) == "NPOIN=")
    {
      SLIC_ERROR_IF(ndime == -1, "dimension must be set prior to parsing the mesh points!");

      npoin = std::stoi(line.substr(6, line.length()));
      points = axom::allocate<double>(npoin * ndime);
      read_points(points, npoin, ndime, ifs);
    }
    else if(line.substr(0, 6) == "NELEM=")
    {
      SLIC_ERROR_IF(ndime == -1, "dimension must be set prior to parsing mesh connectivity!");

      nelem = std::stoi(line.substr(6, line.length()));
      connectivity = axom::allocate<axom::IndexType>(nelem * mint::MAX_CELL_NODES);
      cellTypes = axom::allocate<mint::CellType>(nelem);
      read_connectivity(connectivity, cellTypes, nelem, isMixed, ifs);
    }

  }  // END while

  // sanity checks
  SLIC_ASSERT(points != nullptr);
  SLIC_ASSERT(connectivity != nullptr);
  SLIC_ASSERT(cellTypes != nullptr);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
int read_su2(const std::string& file, Mesh*& mesh)
{
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
  return 0;
}

//------------------------------------------------------------------------------
int write_su2(const mint::Mesh* mesh, const std::string& file)
{
  SLIC_ERROR_IF(mesh == nullptr, "mesh pointer is null!");
  SLIC_ERROR_IF(file.length() <= 0, "SU2 filename is empty!");

  if(mesh->isStructured())
  {
    SLIC_WARNING("SU2 format is supported only for unstructured meshes!");
    return -1;
  }

  std::ofstream ofs(file.c_str());
  if(!ofs.is_open())
  {
    SLIC_WARNING("cannot write to file [" << file << "]");
    return -1;
  }

  const int ndims = mesh->getDimension();
  const axom::IndexType nnodes = mesh->getNumberOfNodes();
  const axom::IndexType ncells = mesh->getNumberOfCells();

  ofs << "NDIME= " << ndims << std::endl << std::endl;
  ofs << "NPOIN= " << nnodes << std::endl;

  double coords[3];
  for(axom::IndexType inode = 0; inode < nnodes; ++inode)
  {
    mesh->getNode(inode, coords);

    for(int idim = 0; idim < ndims; ++idim)
    {
      ofs << coords[idim] << " ";
    }  // END for all dimensions
    ofs << std::endl;

  }  // END for all nodes
  ofs << std::endl;

  ofs << "NELEM= " << ncells << std::endl;
  axom::IndexType cell[mint::MAX_CELL_NODES];
  for(axom::IndexType icell = 0; icell < ncells; ++icell)
  {
    const int ctype = getSU2CellType(mesh->getCellType(icell));
    ofs << ctype << " ";

    const axom::IndexType ncnodes = mesh->getNumberOfCellNodes(icell);

    mesh->getCellNodeIDs(icell, cell);
    for(axom::IndexType inode = 0; inode < ncnodes; ++inode)
    {
      ofs << cell[inode] << " ";
    }  // END for all cell nodes
    ofs << std::endl;

  }  // END for all cells
  ofs << std::endl;

  ofs.close();
  return 0;
}

} /* namespace mint */

} /* namespace axom */
