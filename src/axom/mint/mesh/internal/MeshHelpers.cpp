// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint/mesh/internal/MeshHelpers.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace axom
{
namespace mint
{
namespace internal
{
//------------------------------------------------------------------------------
// IMPLEMENTATION OF FREE METHODS
//------------------------------------------------------------------------------
std::string join_ints_into_string(int count, IndexType* values, char sep)
{
  std::stringstream joined;

  for(int i = 0; i < count; ++i)
  {
    if(i > 0)
    {
      joined << sep;
    }
    joined << values[i];
  }

  return joined.str();
}

//------------------------------------------------------------------------------
std::string make_face_key(int count, IndexType* values, char sep)
{
  std::vector<IndexType> locvalues(values, values + count);
  std::sort(locvalues.begin(), locvalues.end());
  return join_ints_into_string(count, locvalues.data(), sep);
}

//------------------------------------------------------------------------------
IndexType otherSide(IndexType* f2c, IndexType thisSide)
{
  return (f2c[0] != thisSide) ? f2c[0] : f2c[1];
}

//------------------------------------------------------------------------------
struct FaceTypeCellsNodes
{
  FaceTypeCellsNodes() : facetype(UNDEFINED_CELL) { }

  FaceTypeCellsNodes(CellType ftype, std::vector<IndexType>& fcells, std::vector<IndexType>& fnodes)
    : facetype(ftype)
    , facecells(fcells)
    , facenodes(fnodes)
  { }

  CellType facetype;
  std::vector<IndexType> facecells;
  std::vector<IndexType> facenodes;
};

//------------------------------------------------------------------------------
bool initFaces(Mesh* mesh,
               IndexType& facecount,
               Array<IndexType>& f2c,
               Array<IndexType>& c2f,
               Array<IndexType>& c2n,
               Array<IndexType>& c2foffsets,
               Array<IndexType>& f2n,
               Array<IndexType>& f2noffsets,
               Array<CellType>& f2ntypes)
{
  bool success = true;

  facecount = 0;

  using IDtoKeyType = std::unordered_map<IndexType, std::vector<IndexType>>;
  using FaceBuilderType =
    std::map<std::vector<IndexType>, FaceTypeCellsNodes, utilities::LexiComparator<IndexType>>;
  FaceBuilderType workface;

  // Iterate over each cell.
  const IndexType cellcount = mesh->getNumberOfCells();
  for(int c = 0; c < cellcount; ++c)
  {
    // Step 1. For every cell, get the nodes.
    IndexType nodes[MAX_CELL_NODES];
    mesh->getCellNodeIDs(c, nodes);
    const CellType celltype = mesh->getCellType(c);
    const CellInfo thisCell = getCellInfo(celltype);

    int base = 0;
    for(int f = 0; f < thisCell.num_faces; ++f)
    {
      // Step 2. The cell nodes will be in "VTK Order," as specified by
      // https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf.
      // This routine selects the face nodes as listed in the registered
      // cell info to make sure that face normals point outward.
      //
      // For every face, sort the node IDs and join them together in a
      // string with '.' delimiters.  This is the key for an associative
      // array whose value is a list of cell IDs.
      int num_face_nodes = thisCell.face_nodecount[f];
      std::vector<IndexType> inorder_facenodes(num_face_nodes);
      const IndexType* face_node_offsets = thisCell.face_nodes + base;
      for(int fn = 0; fn < num_face_nodes; ++fn)
      {
        inorder_facenodes[fn] = nodes[face_node_offsets[fn]];
      }
      base += num_face_nodes;

      std::vector<IndexType> sorted_facenodes(inorder_facenodes);
      std::sort(sorted_facenodes.begin(), sorted_facenodes.end());

      FaceTypeCellsNodes cells_of_face;
      if(workface.count(sorted_facenodes) > 0)
      {
        cells_of_face = workface[sorted_facenodes];
      }
      else
      {
        cells_of_face.facetype = thisCell.face_types[f];
        cells_of_face.facenodes = inorder_facenodes;
      }
      cells_of_face.facecells.push_back(c);
      workface[sorted_facenodes] = cells_of_face;
    }
  }

  // Step 3. For each face, record its incident cells.
  // Here we use ConnectivityArray::reserve() and then append() each
  // face.  This means the kth inserted face gets ID k.
  IndexType faceID = 0;
  IndexType faceNodeTotal = 0;
  IDtoKeyType keys;
  f2c.resize(2 * workface.size());
  for(FaceBuilderType::value_type v : workface)
  {
    FaceTypeCellsNodes theFace = v.second;
    faceNodeTotal += static_cast<IndexType>(theFace.facenodes.size());
    int faceCellCount = static_cast<int>(theFace.facecells.size());
    IndexType* faceCells = theFace.facecells.data();

    if(faceCellCount < 1 || faceCellCount > 2)
    {
      success = false;
    }
    else
    {
      f2c[2 * faceID] = faceCells[0];
      f2c[2 * faceID + 1] = (faceCellCount > 1 ? faceCells[1] : -1);
    }

    keys[faceID] = v.first;

    faceID += 1;
  }

  // If we have any face with less than one or more than two incident cells,
  // clean up and return failure.  We won't do any more work here.
  if(!success)
  {
    f2c.clear();
    return success;
  }

  // Record how many faces we have in this mesh.
  facecount = faceID;

  // Now that we have a count of all the face-nodes, and we know we have no
  // faces with more than two nodes, record the face-node relations.
  f2n.resize(faceNodeTotal);
  f2noffsets.resize(facecount + 1);
  f2ntypes.resize(facecount);
  int faceNodeOffset = 0;
  for(int fidx = 0; fidx < facecount; ++fidx)
  {
    FaceTypeCellsNodes theFace = workface[keys[fidx]];
    std::vector<IndexType>& faceNodes = theFace.facenodes;
    std::copy(faceNodes.begin(), faceNodes.end(), f2n.begin() + faceNodeOffset);
    f2noffsets[fidx] = faceNodeOffset;
    f2ntypes[fidx] = theFace.facetype;
    faceNodeOffset += static_cast<int>(faceNodes.size());
  }
  f2noffsets[facecount] = faceNodeOffset;

  // Step 4. Now that we have face IDs, record cell-to-face relation.
  using CellFaceBuilderType = std::unordered_map<IndexType, std::vector<IndexType>>;

  CellFaceBuilderType cell_to_face;
  int cellFaceCount = 0;
  for(int f = 0; f < facecount; ++f)
  {
    for(IndexType c : workface[keys[f]].facecells)
    {
      std::vector<IndexType> cell_faces;
      if(cell_to_face.count(c) > 0)
      {
        cell_faces = cell_to_face[c];
      }
      cellFaceCount += 1;
      cell_faces.push_back(f);
      cell_to_face[c] = cell_faces;
    }
  }

  // Step 4b. Put cell-to-face relation into output arrays.
  c2f.resize(cellFaceCount);
  c2n.resize(cellFaceCount);
  c2foffsets.resize(cellcount + 1);
  cellFaceCount = 0;

  for(IndexType cellID = 0; cellID < cellcount; ++cellID)
  {
    // For the current cell ID, set its type and value offset
    //  cellTypes[cellID] = m_cell_connectivity->getIDType(cellID);
    c2foffsets[cellID] = cellFaceCount;

    if(cell_to_face.count(cellID) > 0)
    {
      // If we have faces for the current cell (which we certainly
      // *SHOULD* have), copy them in.
      std::vector<IndexType>& theCells = cell_to_face[cellID];
      IndexType* faceIDs = theCells.data();
      int thisCellFaceCount = static_cast<int>(theCells.size());

      // Copy in values
      for(int f = 0; f < thisCellFaceCount; ++f)
      {
        c2f[cellFaceCount + f] = faceIDs[f];
        c2n[cellFaceCount + f] = otherSide(&f2c[2 * faceIDs[f]], cellID);
      }

      // Maintain offset
      cellFaceCount += thisCellFaceCount;
    }
  }
  c2foffsets[cellcount] = cellFaceCount;

  return success;
}

} /* namespace internal */
} /* namespace mint */
} /* namespace axom */
