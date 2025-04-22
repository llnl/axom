// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_CELLTYPES_HPP_
#define MINT_CELLTYPES_HPP_

#include "axom/mint/config.hpp"

namespace axom
{
namespace mint
{
static constexpr int MAX_CELL_NODES = 27;
static constexpr int MAX_CELL_FACES = 6;
static constexpr int MAX_FACE_NODES = 9;
static constexpr int MAX_ALL_FACES_NODES = MAX_CELL_FACES * MAX_FACE_NODES;

/*!
 * \brief Enumerates all cell types supported by Mint
 */
enum class CellType : signed char
{
  UNDEFINED_CELL = -1,  ///< UNDEFINED

  VERTEX,   ///< VERTEX
  SEGMENT,  ///< LINE_SEGMENT

  TRIANGLE,  ///< LINEAR_TRIANGLE
  QUAD,      ///< LINEAR_QUAD
  TET,       ///< LINEAR_TET
  HEX,       ///< LINEAR_HEX
  PRISM,     ///< LINEAR_PRISM
  PYRAMID,   ///< LINEAR_PYRAMID

  QUAD9,  ///< QUADRATIC QUAD
  HEX27,  ///< QUADRATIC HEX

  NUM_CELL_TYPES  ///<  total number of cell types
};

constexpr CellType UNDEFINED_CELL = CellType::UNDEFINED_CELL;
constexpr CellType VERTEX = CellType::VERTEX;
constexpr CellType SEGMENT = CellType::SEGMENT;
constexpr CellType TRIANGLE = CellType::TRIANGLE;
constexpr CellType QUAD = CellType::QUAD;
constexpr CellType TET = CellType::TET;
constexpr CellType HEX = CellType::HEX;
constexpr CellType PRISM = CellType::PRISM;
constexpr CellType PYRAMID = CellType::PYRAMID;
constexpr CellType QUAD9 = CellType::QUAD9;
constexpr CellType HEX27 = CellType::HEX27;
constexpr int NUM_CELL_TYPES = static_cast<int>(CellType::NUM_CELL_TYPES);

/*!
 * \def REGISTER_CELL_INFO( MINT_CELL_TYPE, MINT_NAME, BP_NAME, VTK_TYPE,
 *                          N_NODES, N_FACES, N_FACE_NODES, FACE_CELL_TYPES,
 *                          FACE_NODES )
 *
 * \brief Convenience macro used to register information about a cell type.
 *
 * \param MINT_CELL_TYPE the mint cell type, e.g., mint::QUAD, mint::HEX, etc.
 * \param MINT_NAME the associated mint name for the cell type.
 * \param BP_NAME the associated name in the mesh blueprint.
 * \param VTK_TYPE the corresponding VTK type.
 * \param N_NODES the number of nodes that the cell has.
 * \param N_FACES the number of faces that the cell has.
 * \param N_FACE_NODES an array; the number of nodes that each face has.
 * \param FACE_CELL_TYPES an array; the VTK type of each face.
 * \param FACE_NODES an array; the node offsets specifying each face
 *        (CCW, so normal points out).
 */
#define REGISTER_CELL_INFO(MINT_CELL_TYPE,                        \
                           MINT_NAME,                             \
                           BP_NAME,                               \
                           VTK_TYPE,                              \
                           N_NODES,                               \
                           N_FACES,                               \
                           N_FACE_NODES,                          \
                           FACE_CELL_TYPES,                       \
                           FACE_NODES)                            \
  namespace internal                                              \
  {                                                               \
  static const CellInfo MINT_CELL_TYPE##_INFO = {MINT_CELL_TYPE,  \
                                                 MINT_NAME,       \
                                                 BP_NAME,         \
                                                 VTK_TYPE,        \
                                                 N_NODES,         \
                                                 N_FACES,         \
                                                 N_FACE_NODES,    \
                                                 FACE_CELL_TYPES, \
                                                 FACE_NODES};     \
  }

/*!
 * \def CELL_INFO( MINT_CELL_TYPE )
 *
 * \brief Convenience macro used to access a registered CellInfo struct for
 *  the specified mint cell type.
 *
 * \param MINT_CELL_TYPE the mint cell type, e.g., mint::QUAD, mint::HEX, etc.
 */
#define CELL_INFO(MINT_CELL_TYPE) internal::MINT_CELL_TYPE##_INFO

/*!
 * \struct CellInfo
 *
 * \brief Holds information associated with a given cell type.
 */
struct CellInfo
{
  CellType cell_type;                        /*!< cell type, e.g. mint::QUAD, mint::HEX */
  const char* name;                          /*!< the name associated with the cell */
  const char* blueprint_name;                /*!< corresponding mesh blueprint name */
  int vtk_type;                              /*!< corresponding vtk_type */
  int num_nodes;                             /*!< number of nodes for the given cell */
  int num_faces;                             /*!< number of faces for the given cell */
  int face_nodecount[MAX_CELL_FACES];        /*!< number of nodes for each of cell's faces */
  CellType face_types[MAX_CELL_FACES];       /*!< face type, e.g. mint::SEGMENT, mint::QUAD */
  IndexType face_nodes[MAX_ALL_FACES_NODES]; /*!< nodes for each of cell's faces */
};

// This construct lets us pass literal arrays to function-like macros.
// AR stands for ARray.
#define AR(...) __VA_ARGS__

// Cell Info registration
REGISTER_CELL_INFO(VERTEX, "VERTEX", "point", 1, 1, 0, AR({0}), AR({UNDEFINED_CELL}), AR({0}));

REGISTER_CELL_INFO(SEGMENT,
                   "SEGMENT",
                   "line",
                   3,
                   2,
                   0,
                   AR({0}),
                   AR({UNDEFINED_CELL}),
                   AR({
                     0,  // face 0
                     1   // face 1
                   }));

REGISTER_CELL_INFO(TRIANGLE,
                   "TRIANGLE",
                   "tri",
                   5,
                   3,
                   3,
                   AR({2, 2, 2}),
                   AR({SEGMENT, SEGMENT, SEGMENT}),
                   AR({
                     0,
                     1,  // face 0
                     1,
                     2,  // face 1
                     2,
                     0  // face 2
                   }));

REGISTER_CELL_INFO(QUAD,
                   "QUAD",
                   "quad",
                   9,
                   4,
                   4,
                   AR({2, 2, 2, 2}),
                   AR({SEGMENT, SEGMENT, SEGMENT, SEGMENT}),
                   AR({
                     0,
                     1,  // face 0
                     1,
                     2,  // face 1
                     2,
                     3,  // face 2
                     3,
                     0  // face 3
                   }));

REGISTER_CELL_INFO(TET,
                   "TET",
                   "tet",
                   10,
                   4,
                   4,
                   AR({3, 3, 3, 3}),
                   AR({TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE}),
                   AR({
                     0,
                     2,
                     1,  // face 0
                     0,
                     3,
                     2,  // face 1
                     0,
                     1,
                     3,  // face 2
                     1,
                     2,
                     3  // face 3
                   }));

REGISTER_CELL_INFO(HEX,
                   "HEX",
                   "hex",
                   12,
                   8,
                   6,
                   AR({4, 4, 4, 4, 4, 4}),
                   AR({QUAD, QUAD, QUAD, QUAD, QUAD, QUAD}),
                   AR({
                     0, 3, 2, 1,  // face 0
                     1, 2, 6, 5,  // face 1
                     1, 5, 4, 0,  // face 2
                     0, 4, 7, 3,  // face 3
                     7, 6, 2, 3,  // face 4
                     4, 5, 6, 7   // face 5
                   }));

REGISTER_CELL_INFO(PRISM,
                   "PRISM",
                   "prism-no-bp",
                   13,
                   6,
                   5,
                   AR({3, 4, 4, 4, 3}),
                   AR({TRIANGLE, QUAD, QUAD, QUAD, TRIANGLE}),
                   AR({
                     0,
                     1,
                     2,  // face 0
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
                     4  // face 4
                   }));

REGISTER_CELL_INFO(PYRAMID,
                   "PYRAMID",
                   "pyramid-no-bp",
                   14,
                   5,
                   5,
                   AR({4, 3, 3, 3, 3}),
                   AR({QUAD, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE}),
                   AR({
                     0,
                     3,
                     2,
                     1,  // face 0
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
                     4  // face 4
                   }));

REGISTER_CELL_INFO(QUAD9,
                   "QUAD9",
                   "quad9-no-bp",
                   28,
                   9,
                   4,
                   AR({2, 2, 2, 2}),
                   AR({SEGMENT, SEGMENT, SEGMENT, SEGMENT}),
                   AR({
                     0,
                     1,  // face 0
                     1,
                     2,  // face 1
                     2,
                     3,  // face 2
                     3,
                     0  // face 3
                   }));

REGISTER_CELL_INFO(HEX27,
                   "HEX27",
                   "hex27-no-bp",
                   29,
                   27,
                   6,
                   AR({9, 9, 9, 9, 9, 9}),
                   AR({QUAD9, QUAD9, QUAD9, QUAD9, QUAD9, QUAD9}),
                   AR({
                     0, 3, 2, 1, 11, 10, 9,  8,  24,  // face 0
                     1, 2, 6, 5, 9,  18, 13, 17, 21,  // face 1
                     1, 5, 4, 0, 17, 12, 16, 8,  22,  // face 2
                     0, 4, 7, 3, 16, 15, 19, 11, 20,  // face 3
                     7, 6, 2, 3, 14, 18, 10, 19, 23,  // face 4
                     4, 5, 6, 7, 12, 13, 14, 15, 25   // face 5
                   }));

/*!
 * \brief Array of CellInfo corresponding to each cell type
 * \note The order at which CellInfo for each type is added has to match
 *  the order of the cell types in the CellTypes enum above.
 */
static const CellInfo cell_info[NUM_CELL_TYPES] = {CELL_INFO(VERTEX),
                                                   CELL_INFO(SEGMENT),
                                                   CELL_INFO(TRIANGLE),
                                                   CELL_INFO(QUAD),
                                                   CELL_INFO(TET),
                                                   CELL_INFO(HEX),
                                                   CELL_INFO(PRISM),
                                                   CELL_INFO(PYRAMID),
                                                   CELL_INFO(QUAD9),
                                                   CELL_INFO(HEX27)};

/*!
 * \brief Return the underlying integer associated with the given CellType.
 *
 * \param [in] type the CellType in question.
 */
inline constexpr int cellTypeToInt(CellType type) { return static_cast<int>(type); }

/*!
 * \brief Return the CellInfo struct associated with the given type.
 *
 * \param [in] type the CellType in question.
 */
inline constexpr const CellInfo& getCellInfo(CellType type)
{
  return cell_info[cellTypeToInt(type)];
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_CellTypes_HPP_ */
