// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/slic.hpp"
#include "axom/mint/mesh/CellTypes.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cstring>

// namespace aliases
namespace mint = axom::mint;

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_cell_types, check_cell_types)
{
  const int num_nodes[mint::NUM_CELL_TYPES] = {
    1,   // VERTEX
    2,   // SEGMENT
    3,   // TRIANGLE
    4,   // QUAD
    4,   // TET
    8,   // HEX
    6,   // PRISM
    5,   // PYRAMID
    9,   // QUAD9
    27,  // HEX27
  };

  const int num_faces[mint::NUM_CELL_TYPES] = {
    0,  // VERTEX
    0,  // SEGMENT
    3,  // TRIANGLE
    4,  // QUAD
    4,  // TET
    6,  // HEX
    5,  // PRISM
    5,  // PYRAMID
    4,  // QUAD9
    6,  // HEX27
  };

  const int face_nodecount[mint::NUM_CELL_TYPES][mint::MAX_FACE_NODES] = {
    {0},                 // VERTEX
    {0},                 // SEGMENT
    {2, 2, 2},           // TRIANGLE
    {2, 2, 2, 2},        // QUAD
    {3, 3, 3, 3},        // TET
    {4, 4, 4, 4, 4, 4},  // HEX
    {3, 4, 4, 4, 3},     // PRISM
    {4, 3, 3, 3, 3},     // PYRAMID
    {2, 2, 2, 2},        // QUAD9
    {9, 9, 9, 9, 9, 9},  // HEX27
  };

  const mint::CellType face_types[mint::NUM_CELL_TYPES][mint::MAX_FACE_NODES] = {
    {mint::UNDEFINED_CELL},                                                          // VERTEX
    {mint::UNDEFINED_CELL},                                                          // SEGMENT
    {mint::SEGMENT, mint::SEGMENT, mint::SEGMENT},                                   // TRIANGLE
    {mint::SEGMENT, mint::SEGMENT, mint::SEGMENT, mint::SEGMENT},                    // QUAD
    {mint::TRIANGLE, mint::TRIANGLE, mint::TRIANGLE, mint::TRIANGLE},                // TET
    {mint::QUAD, mint::QUAD, mint::QUAD, mint::QUAD, mint::QUAD, mint::QUAD},        // HEX
    {mint::TRIANGLE, mint::QUAD, mint::QUAD, mint::QUAD, mint::TRIANGLE},            // PRISM
    {mint::QUAD, mint::TRIANGLE, mint::TRIANGLE, mint::TRIANGLE, mint::TRIANGLE},    // PYRAMID
    {mint::SEGMENT, mint::SEGMENT, mint::SEGMENT, mint::SEGMENT},                    // QUAD9
    {mint::QUAD9, mint::QUAD9, mint::QUAD9, mint::QUAD9, mint::QUAD9, mint::QUAD9},  // HEX27
  };

  const axom::IndexType face_nodes[mint::NUM_CELL_TYPES][mint::MAX_ALL_FACES_NODES] = {
    {0},  //            VERTEX

    {0},  //            SEGMENT

    {0,
     1,  // face 0 for TRIANGLE
     1,
     2,  // face 1
     2,
     0},  // face 2

    {0,
     1,  // face 0 for QUAD
     1,
     2,  // face 1
     2,
     3,  // face 2
     3,
     0},  // face 3

    {0,
     2,
     1,  // face 0 for TET
     0,
     3,
     2,  // face 1
     0,
     1,
     3,  // face 2
     1,
     2,
     3},  // face 3

    {0, 3, 2, 1,   // face 0 for HEX
     1, 2, 6, 5,   // face 1
     1, 5, 4, 0,   // face 2
     0, 4, 7, 3,   // face 3
     7, 6, 2, 3,   // face 4
     4, 5, 6, 7},  // face 5

    {0,
     1,
     2,  // face 0 for PRISM
     0,
     2,
     5,
     3,  // face 1
     0,
     3,
     4,
     1,  // face 2
     1,
     4,
     5,
     2,  // face 3
     3,
     5,
     4},  // face 4

    {0,
     3,
     2,
     1,  // face 0 for PYRAMID
     0,
     1,
     4,  // face 1
     1,
     2,
     4,  // face 2
     2,
     3,
     4,  // face 3
     3,
     0,
     4},  // face 4

    {0,
     1,  // face 0 for QUAD9
     1,
     2,  // face 1
     2,
     3,  // face 2
     3,
     0},  // face 3

    {0, 3, 2, 1, 11, 10, 9,  8,  24,   // face 0 for HEX27
     1, 2, 6, 5, 9,  18, 13, 17, 21,   // face 1
     1, 5, 4, 0, 17, 12, 16, 8,  22,   // face 2
     0, 4, 7, 3, 16, 15, 19, 11, 20,   // face 3
     7, 6, 2, 3, 14, 18, 10, 19, 23,   // face 4
     4, 5, 6, 7, 12, 13, 14, 15, 25},  // face 5
  };

  for(int i = 0; i < mint::NUM_CELL_TYPES; ++i)
  {
    const int cell_id = cellTypeToInt(mint::cell_info[i].cell_type);
    SLIC_INFO("[" << mint::cell_info[i].name << "] =>"
                  << "type=" << cell_id << " "
                  << "blueprint_name=(" << mint::cell_info[i].blueprint_name << ") "
                  << "vtk_type=" << mint::cell_info[i].vtk_type << " "
                  << "num_nodes=" << mint::cell_info[i].num_nodes);

    EXPECT_EQ(cell_id, i);
    EXPECT_TRUE(strlen(mint::cell_info[i].name) > 0);
    EXPECT_TRUE(strlen(mint::cell_info[i].blueprint_name) > 0);
    EXPECT_TRUE(mint::cell_info[i].vtk_type > 0);
    EXPECT_EQ(mint::cell_info[i].num_nodes, num_nodes[i]);
    EXPECT_EQ(mint::cell_info[i].num_faces, num_faces[i]);
    int face_offset = 0;
    for(int faceidx = 0; faceidx < num_faces[i]; ++faceidx)
    {
      int face_ncount = face_nodecount[i][faceidx];
      EXPECT_EQ(mint::cell_info[i].face_nodecount[faceidx], face_ncount);
      EXPECT_EQ(mint::cell_info[i].face_types[faceidx], face_types[i][faceidx]);
      for(int facenidx = 0; facenidx < face_ncount; ++facenidx)
      {
        EXPECT_EQ(mint::cell_info[i].face_nodes[face_offset + facenidx],
                  face_nodes[i][face_offset + facenidx]);
      }
      face_offset += face_nodecount[i][faceidx];
    }

  }  // END for all types
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
