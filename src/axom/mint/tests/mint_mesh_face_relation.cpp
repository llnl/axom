// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/config.hpp"                    /* for mint defintions */
#include "axom/mint/mesh/UnstructuredMesh.hpp"     /* for mint::UnstructuredMesh */
#include "axom/mint/mesh/internal/MeshHelpers.hpp" /* for mint::initFaces */

#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic.hpp"

#include "gtest/gtest.h"

#include <string>
#include <sstream>
#include <vector>

using IndexType = axom::IndexType;

namespace axom
{
namespace mint
{
namespace internal
{
struct MeshFaceTest
{
  MeshFaceTest() : name(""), mesh(nullptr), initShouldSucceed(false), totalFaceCount(-1) { }

  MeshFaceTest(std::string thename,
               Mesh *themesh,
               bool itsInitShouldSucceed,
               IndexType theTotalFaceCount,
               std::vector<IndexType> theCellFaceCount,
               std::vector<IndexType> theCellNeighbors,
               std::vector<CellType> theFaceTypes,
               std::vector<IndexType> theFaceNodes)
    : name(thename)
    , mesh(themesh)
    , initShouldSucceed(itsInitShouldSucceed)
    , totalFaceCount(theTotalFaceCount)
  {
    cellFaceCount = theCellFaceCount;
    cellNeighbors = theCellNeighbors;
    faceTypes = theFaceTypes;
    faceNodes = theFaceNodes;
  }

  ~MeshFaceTest() { delete mesh; }

  std::string name;
  Mesh *mesh;
  bool initShouldSucceed;
  IndexType totalFaceCount;
  std::vector<IndexType> cellFaceCount;
  std::vector<IndexType> cellNeighbors;
  std::vector<CellType> faceTypes;
  std::vector<IndexType> faceNodes;
};

} /* end namespace internal */

/*! Generate the tests in mint_mesh_face_relation.svg.
 *  Caller must clean up.
 */
std::vector<internal::MeshFaceTest *> generateFaceTestCases()
{
  std::vector<internal::MeshFaceTest *> tests;

  // Each test mesh in "Face relation tests.svg" is instantiated
  // and put into tests.
  constexpr int TWO_D = 2;
  constexpr int THREE_D = 3;
  constexpr bool EXPECT_INIT_SUCCESS = true;
  constexpr bool EXPECT_INIT_FAILURE = false;

  {
    // 1. tri  ===============================================================
    UnstructuredMesh<SINGLE_SHAPE> *tri = new UnstructuredMesh<SINGLE_SHAPE>(TWO_D, TRIANGLE);
    double trinodes[] = {0, 0, 1, 0, 0, 1};
    IndexType tricells[] = {0, 1, 2};

    tri->appendNodes(trinodes, 3);
    tri->appendCells(tricells, 1);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces, face count,
                    ("tri",
                     tri,
                     EXPECT_INIT_SUCCESS,
                     3,
                     // per-cell face count, neighbor cells for each cell,
                     {3},
                     {-1, -1, -1},
                     // face types, face nodes.
                     {SEGMENT, SEGMENT, SEGMENT},
                     {0, 1, 1, 2, 2, 0}));
  }

  {
    // 2. two tris  ==========================================================
    UnstructuredMesh<SINGLE_SHAPE> *twotris = new UnstructuredMesh<SINGLE_SHAPE>(TWO_D, TRIANGLE);
    double twotrisnodes[] = {0, 0, 1, 0, 0, 1, 0.8, 1.2};
    IndexType twotriscells[] = {0, 1, 2, 1, 3, 2};

    twotris->appendNodes(twotrisnodes, 4);
    twotris->appendCells(twotriscells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("two tris",
                     twotris,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     5,
                     {3, 3},
                     // neighbor cells for each cell.
                     {-1, 1, -1, -1, -1, 0},
                     // face types,
                     {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 1, 1, 2, 2, 0, 1, 3, 3, 2, 2, 0}));
  }

  {
    // 3. three quads and a tri  =============================================
    UnstructuredMesh<MIXED_SHAPE> *thrqtri = new UnstructuredMesh<MIXED_SHAPE>(TWO_D);
    double thrqtrixs[] = {-1, -1, -1, 0, 0, 0, 1, 1};
    double thrqtriys[] = {-1, 0, 1, -1, 0, 1, 0, 1};
    IndexType thrqtricells[] = {0,
                                3,
                                4,
                                1,  // lower-left quad
                                1,
                                4,
                                5,
                                2,  // upper-left quad
                                3,
                                6,
                                4,  // lower-right tri
                                4,
                                6,
                                7,
                                5};  // upper-right quad
    // since thrqtri is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType thrqtrioffs[] = {0, 4, 8, 11, 15};
    CellType thrqtritype[] = {QUAD, QUAD, TRIANGLE, QUAD};

    thrqtri->appendNodes(thrqtrixs, thrqtriys, 8);
    thrqtri->appendCells(thrqtricells, 4, thrqtrioffs, thrqtritype);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("three quads and a tri",
       thrqtri,
       // should init faces,
       EXPECT_INIT_SUCCESS,
       // total face count, per-cell face count,
       11,
       {4, 4, 3, 4},
       // neighbors at each face of the cells
       {-1, 2, 1, -1, 0, 3, -1, -1, -1, 3, 0, 2, -1, -1, 1},
       // face types,
       {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
       // face nodes
       {0, 3, 3, 4, 4, 1, 1, 0, 4, 5, 5, 2, 2, 1, 3, 6, 6, 4, 6, 7, 7, 5}));
  }

  {
    // 4. four tris, with a hole  ============================================
    UnstructuredMesh<SINGLE_SHAPE> *fourtris = new UnstructuredMesh<SINGLE_SHAPE>(TWO_D, TRIANGLE);
    double fourtrisxs[] = {-1, 0, 0, -.2, .2, 1};
    double fourtrisys[] = {-.1, -1, 0, 1, -.2, 0};
    IndexType fourtriscells[] = {0, 1, 2, 0, 2, 3, 1, 5, 4, 2, 5, 3};

    fourtris->appendNodes(fourtrisxs, fourtrisys, 6);
    fourtris->appendCells(fourtriscells, 4);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("four tris, with a hole",
       fourtris,
       // should init faces, total face count,
       EXPECT_INIT_SUCCESS,
       10,
       // per-cell face count,
       {3, 3, 3, 3},
       // neighbors at each face of each cell.
       {-1, -1, 1, 0, 3, -1, -1, -1, -1, -1, -1, 1},
       // face types,
       {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
       // face nodes
       {0, 1, 1, 2, 2, 0, 2, 3, 3, 0, 1, 5, 5, 4, 4, 1, 2, 5, 5, 3}));
  }

  {
    // 5. one 3D tri  ========================================================
    UnstructuredMesh<SINGLE_SHAPE> *threeDtri = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TRIANGLE);
    double threeDtrinodes[] = {-1, 0, 0, 0, 1, 0.5, 1.2, -.2, 3};
    IndexType threeDtricells[] = {0, 1, 2};

    threeDtri->appendNodes(threeDtrinodes, 3);
    threeDtri->appendCells(threeDtricells, 1);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("one 3D tri",
                     threeDtri,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     3,
                     {3},
                     // neighbors at each face of each cell.
                     {-1, -1, -1},
                     // face types,
                     {SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 1, 1, 2, 2, 0}));
  }

  {
    // 6. a tri not far from a quad  =========================================
    UnstructuredMesh<MIXED_SHAPE> *qandtri = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    double qandtrixs[] = {-1, 1, 0, 0.3, 0.9, 2, 1.4};
    double qandtriys[] = {0, 0, 0.8, 0.9, 0.2, 0.5, 1.5};
    double qandtrizs[] = {1, 0, 0, 0, 0, 0, 0};
    IndexType qandtricells[] = {0,
                                1,
                                2,  // left: the triangle
                                3,
                                4,
                                5,
                                6};  // right: the quad
    // since qandtri is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType qandtrioffs[] = {0, 3, 7};
    CellType qandtritype[] = {TRIANGLE, QUAD};

    qandtri->appendNodes(qandtrixs, qandtriys, qandtrizs, 7);
    qandtri->appendCells(qandtricells, 2, qandtrioffs, qandtritype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("tri separated from quad",
                     qandtri,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     7,
                     {3, 4},
                     // neighbors at each face of the cells
                     {-1, -1, -1, -1, -1, -1, -1},
                     // face types
                     {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 1, 1, 2, 2, 0, 3, 4, 4, 5, 5, 6, 6, 3}));
  }

  {
    // 7. tet from four tris  ================================================
    UnstructuredMesh<SINGLE_SHAPE> *tettris = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TRIANGLE);
    double tettrisxs[] = {0, 1, 1, 1};
    double tettrisys[] = {0, 0, 1, 1};
    double tettriszs[] = {0, -.1, 0.2, 1};
    IndexType tettriscells[] = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};

    tettris->appendNodes(tettrisxs, tettrisys, tettriszs, 4);
    tettris->appendCells(tettriscells, 4);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("tet from four tris",
                     tettris,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     6,
                     {3, 3, 3, 3},
                     // neighbors at each face of each cell.
                     {3, 2, 1, 0, 2, 3, 0, 3, 1, 0, 1, 2},
                     // face types,
                     {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 1, 1, 2, 2, 0, 1, 3, 3, 0, 2, 3}));
  }

  {
    // 8. hex from six quads  ================================================
    UnstructuredMesh<SINGLE_SHAPE> *hexquads = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, QUAD);
    double hexquadsxs[] = {0, 1, 1, 0, 0, 1, 1, 0};
    double hexquadsys[] = {0, 0, 1, 1, 0, 0, 1, 1};
    double hexquadszs[] = {0, 0, 0, 0, 1, 1, 1, 1};
    IndexType hexquadscells[] = {0, 3, 2, 1, 0, 1, 5, 4, 1, 2, 6, 5,
                                 2, 3, 7, 6, 3, 0, 4, 7, 4, 5, 6, 7};

    hexquads->appendNodes(hexquadsxs, hexquadsys, hexquadszs, 8);
    hexquads->appendCells(hexquadscells, 6);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("hex from six quads",
       hexquads,
       // should init faces,
       EXPECT_INIT_SUCCESS,
       // total face count, per-cell face count,
       12,
       {4, 4, 4, 4, 4, 4},
       // neighbors at each face of each cell.
       {4, 3, 2, 1, 0, 2, 5, 4, 0, 3, 5, 1, 0, 4, 5, 2, 0, 1, 5, 3, 1, 2, 3, 4},
       // face types
       {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
       // face nodes
       {0, 1, 1, 2, 2, 3, 3, 0, 1, 5, 5, 4, 4, 0, 2, 6, 6, 5, 3, 7, 7, 6, 4, 7}));
  }

  {
    // 9. pyramid from a quad and four tris  =================================
    UnstructuredMesh<MIXED_SHAPE> *pyr = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    double pyrxs[] = {-1, 0, 0, 1, 0};
    double pyrys[] = {0, -1, 1, 0, 0};
    double pyrzs[] = {0, 0, 0, 0, 1};
    IndexType pyrcells[] = {0,
                            1,
                            3,
                            2,  // quad for the base
                            0,
                            1,
                            4,  // four tris for the sides
                            1,
                            3,
                            4,
                            3,
                            2,
                            4,
                            2,
                            0,
                            4};
    // since pyr is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType pyroffs[] = {0, 4, 7, 10, 13, 16};
    CellType pyrtype[] = {QUAD, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE};

    pyr->appendNodes(pyrxs, pyrys, pyrzs, 5);
    pyr->appendCells(pyrcells, 5, pyroffs, pyrtype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("hollow pyramid",
                     pyr,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     8,
                     {4, 3, 3, 3, 3},
                     // neighbors at each face of the cells
                     {1, 2, 3, 4, 0, 2, 4, 0, 3, 1, 0, 4, 2, 0, 1, 3},
                     // face types
                     {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 1, 1, 3, 3, 2, 2, 0, 1, 4, 4, 0, 3, 4, 2, 4}));
  }

  {
    // 10. two tris back to back, forming a closed surface  ==================
    UnstructuredMesh<SINGLE_SHAPE> *b2btris = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TRIANGLE);
    double b2btrisnodes[] = {0, 0, 0, 2, -.3, -.1, 1, 1, 1};
    IndexType b2btriscells[] = {0, 1, 2, 0, 2, 1};

    b2btris->appendNodes(b2btrisnodes, 3);
    b2btris->appendCells(b2btriscells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("back-to-back tris",
                     b2btris,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     3,
                     {3, 3},
                     // neighbors at each face of each cell.
                     {1, 1, 1, 0, 0, 0},
                     // face types
                     {SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 1, 1, 2, 2, 0}));
  }

  {
    // 11. three quads (corner of a box)  ====================================
    UnstructuredMesh<MIXED_SHAPE> *threeq = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    // for variety, you can use a MIXED_SHAPE for a homogeneous mesh---
    // it just means a little more typing.
    double threeqxs[] = {0, 1, 1, 0, 1.4, 1.4, 0.4};
    double threeqys[] = {0, 0, 1, 1, 0.4, 1.4, 1.4};
    double threeqzs[] = {0, 0, 0, 0, 0.4, 0.4, 0.4};
    IndexType threeqcells[] = {0, 1, 2, 3, 1, 4, 5, 2, 2, 5, 6, 3};
    // since threeq is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType threeqoffs[] = {0, 4, 8, 12};
    CellType threeqtype[] = {QUAD, QUAD, QUAD};

    threeq->appendNodes(threeqxs, threeqys, threeqzs, 7);
    threeq->appendCells(threeqcells, 3, threeqoffs, threeqtype);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer, should init faces,
      ("box corner",
       threeq,
       EXPECT_INIT_SUCCESS,
       // total face count, per-cell face count,
       9,
       {4, 4, 4},
       // neighbors at each face of the cells
       {-1, 1, 2, -1, -1, -1, 2, 0, 1, -1, -1, 0},
       // face types
       {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
       // face nodes
       {0, 1, 1, 2, 2, 3, 3, 0, 1, 4, 4, 5, 5, 2, 5, 6, 6, 3}));
  }

  {
    // 12. two quads, two tris  ==============================================
    UnstructuredMesh<MIXED_SHAPE> *twoqtwot = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    double twoqtwotxs[] = {0, 1, 1, 0, 1.4, 1.4, 0.4};
    double twoqtwotys[] = {0, 0, 1, 1, 0.4, 1.4, 1.4};
    double twoqtwotzs[] = {0, 0, 0, 0, 0.4, 0.4, 0.4};
    IndexType twoqtwotcells[] = {0, 1, 2, 1, 4, 5, 2, 2, 5, 6, 3, 3, 0, 2};
    // since twoqtwot is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType twoqtwotoffs[] = {0, 3, 7, 11, 14};
    CellType twoqtwottype[] = {TRIANGLE, QUAD, QUAD, TRIANGLE};

    twoqtwot->appendNodes(twoqtwotxs, twoqtwotys, twoqtwotzs, 7);
    twoqtwot->appendCells(twoqtwotcells, 4, twoqtwotoffs, twoqtwottype);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("box corner (tris and quads)",
       twoqtwot,
       // should init faces,
       EXPECT_INIT_SUCCESS,
       // total face count, per-cell face count,
       10,
       {3, 4, 4, 3},
       // neighbors at each face of the cells
       {-1, 1, 3, -1, -1, 2, 0, 1, -1, -1, 3, -1, 0, 2},
       // face types
       {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
       // face nodes
       {0, 1, 1, 2, 2, 0, 1, 4, 4, 5, 5, 2, 5, 6, 6, 3, 3, 2, 3, 0}));
  }

  {
    // 13. two quads, two tris forming a prism open at both ends  ============
    UnstructuredMesh<MIXED_SHAPE> *oprism = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    double oprismxs[] = {0, 1, 2, 0, 1, 2};
    double oprismys[] = {0, -1, 1, 0, -1, 1};
    double oprismzs[] = {0, 0, 0, 1, 1, 1};
    IndexType oprismcells[] = {0, 1, 4, 3, 1, 2, 4, 2, 5, 4, 5, 2, 0, 3};
    // since oprism is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType oprismoffs[] = {0, 4, 7, 10, 14};
    CellType oprismtype[] = {QUAD, TRIANGLE, TRIANGLE, QUAD};

    oprism->appendNodes(oprismxs, oprismys, oprismzs, 6);
    oprism->appendCells(oprismcells, 4, oprismoffs, oprismtype);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("open prism (quads, tris)",
       oprism,
       // should init faces,
       EXPECT_INIT_SUCCESS,
       // total face count, per-cell face count,
       10,
       {4, 3, 3, 4},
       // neighbors at each face of the cells
       {-1, 1, -1, 3, -1, 2, 0, 3, -1, 1, 2, -1, 0, -1},
       // face types
       {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
       // face nodes
       {0, 1, 1, 4, 4, 3, 3, 0, 1, 2, 2, 4, 2, 5, 5, 4, 2, 0, 3, 5}));
  }

  {
    // 14. cracked tet  ======================================================
    UnstructuredMesh<SINGLE_SHAPE> *crackedtet =
      new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TRIANGLE);
    double crackedtetxs[] = {0, 1, 1, 1, 0.9};
    double crackedtetys[] = {0, 0, 1, 1, 0.9};
    double crackedtetzs[] = {0, -.1, 0.2, 1, 1};
    IndexType crackedtetcells[] = {0, 2, 1, 0, 1, 4, 1, 2, 3, 2, 0, 3};

    crackedtet->appendNodes(crackedtetxs, crackedtetys, crackedtetzs, 5);
    crackedtet->appendCells(crackedtetcells, 4);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("cracked tet",
                     crackedtet,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     8,
                     {3, 3, 3, 3},
                     // neighbors at each face of each cell
                     {3, 2, 1, 0, -1, -1, 0, 3, -1, 0, -1, 2},
                     // face types,
                     {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 2, 2, 1, 1, 0, 1, 4, 4, 0, 2, 3, 3, 1, 0, 3}));
  }

  {
    // 15. cracked pyramid  ==================================================
    UnstructuredMesh<MIXED_SHAPE> *crackedpyr = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    double crackedpyrxs[] = {-1, 0, 0, 1, 0, 0.2};
    double crackedpyrys[] = {0, -1, 1, 0, 0, -0.2};
    double crackedpyrzs[] = {0, 0, 0, 0, 1, 1};
    IndexType crackedpyrcells[] = {0,
                                   1,
                                   3,
                                   2,  // quad for the base
                                   0,
                                   1,
                                   4,  // four tris for the sides
                                   1,
                                   3,
                                   5,
                                   3,
                                   2,
                                   4,
                                   2,
                                   0,
                                   4};
    // since crackedpyr is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType crackedpyroffs[] = {0, 4, 7, 10, 13, 16};
    CellType crackedpyrtype[] = {QUAD, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE};

    crackedpyr->appendNodes(crackedpyrxs, crackedpyrys, crackedpyrzs, 6);
    crackedpyr->appendCells(crackedpyrcells, 5, crackedpyroffs, crackedpyrtype);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("cracked pyramid",
       crackedpyr,
       // should init faces,
       EXPECT_INIT_SUCCESS,
       // total face count, per-cell face count,
       10,
       {4, 3, 3, 3, 3},
       // neighbors at each face of the cells
       {1, 2, 3, 4, 0, -1, 4, 0, -1, -1, 0, 4, -1, 0, 1, 3},
       // face types,
       {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
       // face nodes
       {0, 1, 1, 3, 3, 2, 2, 0, 1, 4, 4, 0, 3, 5, 5, 1, 3, 4, 4, 2}));
  }

  {
    // 16. not a manifold  ===================================================
    UnstructuredMesh<SINGLE_SHAPE> *notmanf = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TRIANGLE);
    double notmanfxs[] = {-1, 0, 0, 0.2, 1};
    double notmanfys[] = {0, 0, 0, -.6, 0};
    double notmanfzs[] = {0.6, 0, 1, 0.4, 0.4};
    IndexType notmanfcells[] = {0, 1, 2, 3, 1, 2, 4, 1, 2};

    notmanf->appendNodes(notmanfxs, notmanfys, notmanfzs, 5);
    notmanf->appendCells(notmanfcells, 3);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("not a manifold",
                     notmanf,
                     // should NOT init faces,
                     EXPECT_INIT_FAILURE,
                     // total face count, per-cell face count,
                     7,
                     {3, 3, 3},
                     // neighbors at each face of each cell
                     // can't be properly expressed.
                     {-1, -1, -1, -1, -1, -1, -1, -1, -1},
                     // face types,
                     {SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT, SEGMENT},
                     // face nodes
                     {0, 1, 1, 2, 2, 0, 3, 1, 2, 3, 4, 1, 2, 4}));
  }

  {
    // 17. egregiously not a manifold  =======================================
    UnstructuredMesh<MIXED_SHAPE> *egreg = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    double egregxs[] = {0, 1, 1, 0, 0, 1, 0, 1, 1};
    double egregys[] = {0, 0, 1, 1, 0.8, 0, 0, -0.4, -0.9};
    double egregzs[] = {0, 0, 0, 0, 0.4, 1, 1, 1.2, 0.5};
    IndexType egregcells[] = {0, 1, 2, 3, 0, 1, 4, 0, 1, 5, 6, 0, 1, 7, 0, 1, 8};
    // since egreg is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType egregoffs[] = {0, 4, 7, 11, 14, 17};
    CellType egregtype[] = {QUAD, TRIANGLE, QUAD, TRIANGLE, TRIANGLE};

    egreg->appendNodes(egregxs, egregys, egregzs, 9);
    egreg->appendCells(egregcells, 5, egregoffs, egregtype);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("egregiously not a manifold",
       egreg,
       //  should NOT init faces,
       EXPECT_INIT_FAILURE,
       // total face count, per-cell face count,
       13,
       {4, 3, 4, 3, 3},
       // neighbors at each face of each cell
       // can't be properly expressed.
       {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       // face types
       {SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT,
        SEGMENT},
       // face nodes
       {0, 1, 1, 2, 2, 3, 3, 4, 1, 4, 4, 0, 1, 5, 5, 6, 6, 0, 1, 7, 7, 0, 1, 8, 8, 0}));
  }

  {
    // 18. 3D tet  ===========================================================
    UnstructuredMesh<SINGLE_SHAPE> *tet = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TET);
    double tetxs[] = {0, 1, 1, 1};
    double tetys[] = {0, 0, 1, 1};
    double tetzs[] = {0, -.1, 0.2, 1};
    IndexType tetcells[] = {0, 1, 2, 3};

    tet->appendNodes(tetxs, tetys, tetzs, 4);
    tet->appendCells(tetcells, 1);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("3D tet",
                     tet,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     4,
                     {4},
                     // neighbors at each face of each cell.
                     {-1, -1, -1, -1},
                     // face types
                     {TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE},
                     // face nodes
                     {0, 2, 1, 0, 3, 2, 0, 1, 3, 1, 2, 3}));
  }

  {
    // 19. two hexs  =========================================================
    UnstructuredMesh<SINGLE_SHAPE> *hexs = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, HEX);
    double hexsxs[] = {0, 1, 2, 0, 1, 3, 0, 1, 2, 0, 1, 3};
    double hexsys[] = {0, 0, 0, 1, 0.8, 1, 0, 0.2, 0, 1, 1, 1};
    double hexszs[] = {0, 0.2, 0, 0, 0, 0, 1, 1, 1, 1, 0.8, 1};
    IndexType hexscells[] = {0, 1, 4, 3, 6, 7, 10, 9, 1, 2, 5, 4, 7, 8, 11, 10};

    hexs->appendNodes(hexsxs, hexsys, hexszs, 12);
    hexs->appendCells(hexscells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("two hexs",
                     hexs,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     11,
                     {6, 6},
                     // neighbors at each face of each cell.
                     {-1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 0, -1},
                     // face types
                     {QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD},
                     // face nodes
                     {0,  3, 4, 1, 1, 4, 10, 7, 1,  7, 6, 0, 0, 6, 9,  3,  9, 10, 4, 3, 6,  7,
                      10, 9, 1, 4, 5, 2, 2,  5, 11, 8, 2, 8, 7, 1, 10, 11, 5, 4,  7, 8, 11, 10}));
  }

  {
    // 20. three hexs  =======================================================
    UnstructuredMesh<SINGLE_SHAPE> *hexs3 = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, HEX);
    double hexs3xs[] = {0, 0, 1, 1, 1, 2, 2, 0, 0, 1, 1, 1, 2, 2};
    double hexs3ys[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1};
    double hexs3zs[] = {0, 2, -1, 1, 3, 0, 2, 0, 2, -1, 1, 3, 0, 2};
    IndexType hexs3cells[] = {0, 3, 10, 7,  1, 4, 11, 8, 3, 5, 12, 10,
                              4, 6, 13, 11, 0, 2, 9,  7, 3, 5, 12, 10};

    hexs3->appendNodes(hexs3xs, hexs3ys, hexs3zs, 14);
    hexs3->appendCells(hexs3cells, 3);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer, should init faces,
      ("three hexs",
       hexs3,
       EXPECT_INIT_SUCCESS,
       // total face count, per-cell face count,
       15,
       {6, 6, 6},
       // neighbors at each face of each cell.
       {2, -1, 1, -1, -1, -1, 2, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, 1},
       // face types
       {QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD, QUAD},
       // face nodes
       {0, 7, 10, 3,  3, 10, 11, 4, 3, 4,  1,  0, 0, 1, 8, 7, 8,  11, 10, 7,
        1, 4, 11, 8,  3, 10, 12, 5, 5, 12, 13, 6, 5, 6, 4, 3, 11, 13, 12, 10,
        4, 6, 13, 11, 0, 7,  9,  2, 2, 9,  12, 5, 2, 5, 3, 0, 10, 12, 9,  7}));
  }

  {
    // 21. two coincident tets, one inside-out, forming a closed manifold  ===
    UnstructuredMesh<SINGLE_SHAPE> *mtet = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TET);
    double mtetxs[] = {0, 1, 1, 1};
    double mtetys[] = {0, 0, 1, 1};
    double mtetzs[] = {0, -.1, 0.2, 1};
    IndexType mtetcells[] = {0, 1, 2, 3, 0, 2, 1, 3};

    mtet->appendNodes(mtetxs, mtetys, mtetzs, 4);
    mtet->appendCells(mtetcells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("two \"back-to-back\" tets",
                     mtet,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     4,
                     {4, 4},
                     // neighbors at each face of each cell.
                     {1, 1, 1, 1, 0, 0, 0, 0},
                     // face types
                     {TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE},
                     // face nodes
                     {0, 2, 1, 0, 3, 2, 0, 1, 3, 1, 2, 3}));
  }

  {
    // 22. bad tet mesh (not a manifold)  ====================================
    UnstructuredMesh<SINGLE_SHAPE> *badtets = new UnstructuredMesh<SINGLE_SHAPE>(THREE_D, TET);
    double badtetsxs[] = {0, 1, 1, 1, 0, 0.3};
    double badtetsys[] = {0, 0, 1, 1, 1, 1.2};
    double badtetszs[] = {0, -.1, 0.2, 1, 0.5, 0.8};
    IndexType badtetscells[] = {0, 1, 2, 3, 0, 3, 2, 4, 0, 3, 2, 5};

    badtets->appendNodes(badtetsxs, badtetsys, badtetszs, 6);
    badtets->appendCells(badtetscells, 3);
    tests.push_back(
      new internal::MeshFaceTest
      // test name, test pointer,
      ("non-manifold tet mesh",
       badtets,
       // should NOT init faces,
       EXPECT_INIT_FAILURE,
       // total face count, per-cell face count,
       10,
       {4, 4, 4},
       // neighbors at each face of each cell can't be properly
       // expressed.
       {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       // face types
       {TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE},
       // face nodes
       {0, 2, 1, 0, 3, 2, 0, 1, 3, 1, 2, 3, 0, 4, 2, 0, 3, 4, 3, 2, 4, 0, 5, 2, 0, 3, 5, 3, 2, 5}));
  }

  {
    // 23. 3D pyramid  =======================================================
    UnstructuredMesh<MIXED_SHAPE> *pyramid = new UnstructuredMesh<MIXED_SHAPE>(THREE_D);
    double pyramidxs[] = {-1, 0, 0, 1, 0};
    double pyramidys[] = {0, -1, 1, 0, 0};
    double pyramidzs[] = {0, 0, 0, 0, 1};
    IndexType pyramidcells[] = {0, 1, 3, 2, 4};

    pyramid->appendNodes(pyramidxs, pyramidys, pyramidzs, 5);
    pyramid->appendCell(pyramidcells, PYRAMID);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("3D pyramid",
                     pyramid,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     5,
                     {5},
                     // neighbors at each face of each cell.
                     {-1, -1, -1, -1, -1},
                     // face types
                     {QUAD, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE},
                     // face nodes
                     {0, 2, 3, 1, 0, 1, 4, 1, 3, 4, 3, 2, 4, 2, 0, 4}));
  }

  return tests;
}

/*! \brief Return true if testnbrs contains the same values as nbrs,
 *         in any order.  Return false otherwise.
 *
 * \param [in] facecount the number of entries to check
 * \param [in] testnbrs an array of values to check against a standard
 * \param [in] nbrs the supplied correct values
 *
 * \returns true if testnbrs contains all values in nbrs and none extra.
 */
bool verifyNeighbors(IndexType facecount, IndexType *testnbrs, IndexType *nbrs)
{
  std::map<IndexType, int> testnbrset, nbrset;

  for(IndexType i = 0; i < facecount; ++i)
  {
    if(testnbrset.count(testnbrs[i]) < 1)
    {
      testnbrset[testnbrs[i]] = 1;
    }
    else
    {
      testnbrset[testnbrs[i]] = testnbrset[testnbrs[i]] + 1;
    }

    if(nbrset.count(nbrs[i]) < 1)
    {
      nbrset[nbrs[i]] = 1;
    }
    else
    {
      nbrset[nbrs[i]] = nbrset[nbrs[i]] + 1;
    }
  }

  for(auto nbrcount : nbrset)
  {
    if(testnbrset.count(nbrcount.first) != 1)
    {
      return false;
    }
    if(testnbrset[nbrcount.first] != nbrcount.second)
    {
      return false;
    }
    testnbrset.erase(nbrcount.first);
  }
  for(auto testnbrcount : testnbrset)
  {
    if(testnbrcount.first != -1)
    {
      return false;
    }
  }

  return true;
}

/*! Helper struct that associates the type of a face with its cells. */
struct FaceTypeNodes
{
  FaceTypeNodes() : facetype(UNDEFINED_CELL) { }

  FaceTypeNodes(CellType ftype, std::vector<IndexType> &fnodes) : facetype(ftype), facenodes(fnodes)
  { }

  CellType facetype;
  std::vector<IndexType> facenodes;
};

/*! If fn is empty, return true.  Otherwise, return false
 *  and compose an error message.
 */
bool checkAndReportFaceNodes(std::map<std::string, FaceTypeNodes> &fn,
                             std::string label,
                             std::stringstream &mesg)
{
  bool success = true;
  if(fn.size() > 0)
  {
    int fcount = static_cast<int>(fn.size());
    mesg << fcount << label << std::endl;
    for(auto fit = fn.begin(), fend = fn.end(); fit != fend; ++fit)
    {
      FaceTypeNodes &ftn = fit->second;
      mesg << "Type " << getCellInfo(ftn.facetype).name << " (";
      mesg << internal::join_ints_into_string(static_cast<int>(ftn.facenodes.size()),
                                              ftn.facenodes.data(),
                                              ' ');
      mesg << ") " << std::endl;
    }
    success = false;
  }

  return success;
}

/*! Return true if a and b contain the same contents, even if shifted. */
template <typename T>
bool matchRotateList(std::vector<T> &a, std::vector<T> &b)
{
  if(a.size() != b.size() || a.size() < 1)
  {
    return false;
  }

  bool matches = false;
  int elts = static_cast<int>(a.size());
  for(int offset = 0; offset < elts && !matches; ++offset)
  {
    matches = true;
    for(int i = 0; i < elts && matches; ++i)
    {
      if(a[(offset + i) % elts] != b[i])
      {
        matches = false;
      }
    }
  }

  return matches;
}

/*! Check face type and (possibly shifted) face nodes for equality. */
bool faceMatches(FaceTypeNodes &a, FaceTypeNodes &b)
{
  return (a.facetype == b.facetype && matchRotateList(a.facenodes, b.facenodes));
}

/*! \brief Check if a list of faces matches a supplied correct answer.
 *
 * \param [in] fcount the number of faces to test
 * \param [in] f2n an array of nodes comprising all the faces to test
 * \param [in] f2noffsets the offset of each face into f2n
 * \param [in] f2ntypes the type of each face
 * \param [in] stdfacecount the number of faces in the correct answer
 * \param [in] stdFaceNodes an array of nodes comprising all faces
 *             in the correct answer
 * \param [in] stdFaceTypes the type of each correct face
 * \param [out] errmesg an explanation of any problems found
 *
 * \returns success true if test faces match correct answer
 *
 * Each face in the test arrays (f2n) must match a face in the correct
 * answer (std).  "Match" means that the face type must be equal, there
 * must not be missing or extra nodes in a face, and the node order
 * (modulo start point) must match.  The order of faces can differ between
 * test and answer data.
 */
bool verifyFaceNodesTypes(int fcount,
                          IndexType *f2n,
                          IndexType *f2noffsets,
                          CellType *f2ntypes,
                          int stdfacecount,
                          IndexType *stdFaceNodes,
                          CellType *stdFaceTypes,
                          std::string &errmesg)
{
  using FaceBuilderType = std::map<std::string, FaceTypeNodes>;

  FaceBuilderType testfnodes, stdfnodes;

  // First, build a map of the test data face nodes.
  for(int f = 0; f < fcount; ++f)
  {
    int nodecount = f2noffsets[f + 1] - f2noffsets[f];
    std::string key = internal::make_face_key(nodecount, &f2n[f2noffsets[f]], '.');
    std::vector<IndexType> facenodes(nodecount);
    std::copy(f2n + f2noffsets[f], f2n + f2noffsets[f] + nodecount, facenodes.begin());
    testfnodes[key] = FaceTypeNodes(f2ntypes[f], facenodes);
  }

  // Then, build a map of the standard face nodes.  Use the same hash scheme.
  int offset = 0;
  for(int f = 0; f < stdfacecount; ++f)
  {
    int nodecount = getCellInfo(stdFaceTypes[f]).num_nodes;
    std::string key = internal::make_face_key(nodecount, &stdFaceNodes[offset], '.');
    std::vector<IndexType> facenodes(nodecount);
    std::copy(stdFaceNodes + offset, stdFaceNodes + offset + nodecount, facenodes.begin());
    stdfnodes[key] = FaceTypeNodes(stdFaceTypes[f], facenodes);
    offset += nodecount;
  }

  // With maps built of both test and standard, compare keys.
  // Iterate over standard, removing each matching face from both standard
  // and test.
  for(auto sit = stdfnodes.begin(), send = stdfnodes.end(); sit != send;)
  {
    if(testfnodes.count(sit->first) == 1 && faceMatches(sit->second, testfnodes[sit->first]))
    {
      testfnodes.erase(sit->first);
      sit = stdfnodes.erase(sit);
    }
    else
    {
      ++sit;
    }
  }

  // If anything remains in stdfnodes, the code missed at least one of
  // the faces.  If anything remains in testfnodes, the code came up with at
  // least one extra face.
  std::stringstream mesg;
  bool allthere = checkAndReportFaceNodes(stdfnodes, " missed faces (from standard):", mesg);
  bool noextra = checkAndReportFaceNodes(testfnodes, " extra faces (not in standard):", mesg);
  errmesg = mesg.str();

  return allthere && noextra;
}

/*! \brief Worker function for the mint face relation generation method.
 *
 * The gtest TEST function drives the test---sets up the list of MeshFaceTest
 * objects, calls this function on each MeshFaceTest, and cleans up.
 * This function does the following:
 *
 * -# calls internal::initFaces with the Mesh pointer supplied in t
 * -# checks if initFaces failed, or succeeded when it shouldn't, reporting
 *    appropriately
 * -# checks the following against answers supplied in t, writing an error
 *    message for each failure
 *    -# number of faces on each cell
 *    -# total number of faces in the mesh
 *    -# each cell's neighbors' IDs (cell by cell check)
 *    -# all faces' type and nodes (all at once, not face by face, because
 *       supplied face data may be in a different order than calculated
 *       faces)
 *
 * This function uses gtest's SCOPED_TRACE to distinguish the tests it runs
 * and uses EXPECT_ predicates to record test success or failure.
 */
void runMeshFaceTest(internal::MeshFaceTest *t)
{
  IndexType facecount = -1;
  Array<IndexType> f2c;
  Array<IndexType> c2f;
  Array<IndexType> c2n;
  Array<IndexType> c2foffsets;
  Array<IndexType> f2n;
  Array<IndexType> f2noffsets;
  Array<CellType> f2ntypes;

  bool initresult =
    internal::initFaces(t->mesh, facecount, f2c, c2f, c2n, c2foffsets, f2n, f2noffsets, f2ntypes);

  if(!initresult && !t->initShouldSucceed)
  {
    SCOPED_TRACE(t->name);
    SUCCEED();
  }
  else if(!initresult && t->initShouldSucceed)
  {
    FAIL() << "test mesh \"" << t->name << "\" call to initFaces() failed but should have succeeded.";
  }
  else if(initresult && !(t->initShouldSucceed))
  {
    FAIL() << "test mesh \"" << t->name << "\" call to initFaces() succeeded but should have failed.";
  }
  else
  {
    // initFaces() succeeded where it should have!
    SCOPED_TRACE(t->name);

    Mesh *m = t->mesh;

    // do we have enough faces on each of the cells?
    IndexType cellcount = m->getNumberOfCells();
    for(IndexType c = 0; c < cellcount; ++c)
    {
      IndexType thisCellFaceCount = c2foffsets[c + 1] - c2foffsets[c];
      EXPECT_EQ(thisCellFaceCount, t->cellFaceCount[c])
        << "thisCellFaceCount " << thisCellFaceCount << " differs from t->cellFaceCount[c] "
        << t->cellFaceCount[c] << " with c == " << c;
    }
    EXPECT_EQ(facecount, t->totalFaceCount);

    // does each cell have the neighbors we expect?
    for(IndexType c = 0; c < cellcount; ++c)
    {
      std::stringstream errmesg;

      IndexType thisCellFaceCount = c2foffsets[c + 1] - c2foffsets[c];
      bool neighborsMatched =
        verifyNeighbors(thisCellFaceCount, &c2n[c2foffsets[c]], &(t->cellNeighbors[c2foffsets[c]]));
      if(!neighborsMatched)
      {
        errmesg << "Cell " << c << " expected neighbors ";

        facecount = t->cellFaceCount[c];
        for(IndexType i = 0; i < facecount; ++i)
        {
          errmesg << t->cellNeighbors[c2foffsets[c] + i] << " ";
        }
        errmesg << std::endl << "found neighbors ";
        for(int i = 0; i < thisCellFaceCount; ++i)
        {
          errmesg << c2n[c2foffsets[c] + i] << " ";
        }
      }
      EXPECT_TRUE(neighborsMatched) << errmesg.str();
    }

    // do all faces have the type and nodes we expect?
    std::string errmesg;
    bool faceNodeTypeMatched = verifyFaceNodesTypes(facecount,
                                                    f2n.data(),
                                                    f2noffsets.data(),
                                                    f2ntypes.data(),
                                                    t->totalFaceCount,
                                                    t->faceNodes.data(),
                                                    t->faceTypes.data(),
                                                    errmesg);
    EXPECT_TRUE(faceNodeTypeMatched) << errmesg;
  }
}

} /* end namespace mint */
} /* end namespace axom */

TEST(mint_mesh_face_relation, tf_verifyFaceNodesTypes)
{
  using namespace axom::mint;

  {
    SCOPED_TRACE("different count fails");

    int fcount = 3;
    IndexType f2n[] = {0, 1, 1, 2, 2, 0};
    IndexType f2noffsets[] = {0, 2, 4, 6};
    CellType f2ntypes[] = {SEGMENT, SEGMENT, SEGMENT};
    int stdfacecount = 2;
    IndexType stdFaceNodes[] = {0, 1, 1, 2};
    CellType stdFaceTypes[] = {SEGMENT, SEGMENT};
    std::string errmesg;

    EXPECT_FALSE(verifyFaceNodesTypes(fcount,
                                      f2n,
                                      f2noffsets,
                                      f2ntypes,
                                      stdfacecount,
                                      stdFaceNodes,
                                      stdFaceTypes,
                                      errmesg));
    EXPECT_TRUE(errmesg.size() > 0);
  }

  {
    SCOPED_TRACE("one face type different fails");

    int fcount = 3;
    IndexType f2n[] = {0, 1, 1, 2, 2, 0};
    IndexType f2noffsets[] = {0, 2, 4, 6};
    CellType f2ntypes[] = {SEGMENT, SEGMENT, SEGMENT};
    int stdfacecount = 3;
    IndexType stdFaceNodes[] = {0, 1, 1, 2, 2, 3, 0};
    CellType stdFaceTypes[] = {SEGMENT, SEGMENT, TRIANGLE};
    std::string errmesg;

    EXPECT_FALSE(verifyFaceNodesTypes(fcount,
                                      f2n,
                                      f2noffsets,
                                      f2ntypes,
                                      stdfacecount,
                                      stdFaceNodes,
                                      stdFaceTypes,
                                      errmesg));
    EXPECT_TRUE(errmesg.size() > 0);
  }

  {
    SCOPED_TRACE("one node different fails");

    int fcount = 3;
    IndexType f2n[] = {0, 1, 1, 2, 2, 0};
    IndexType f2noffsets[] = {0, 2, 4, 6};
    CellType f2ntypes[] = {SEGMENT, SEGMENT, SEGMENT};
    int stdfacecount = 3;
    IndexType stdFaceNodes[] = {0, 1, 1, 2, 2, 3};
    CellType stdFaceTypes[] = {SEGMENT, SEGMENT, SEGMENT};
    std::string errmesg;

    EXPECT_FALSE(verifyFaceNodesTypes(fcount,
                                      f2n,
                                      f2noffsets,
                                      f2ntypes,
                                      stdfacecount,
                                      stdFaceNodes,
                                      stdFaceTypes,
                                      errmesg));
    EXPECT_TRUE(errmesg.size() > 0);
  }

  {
    SCOPED_TRACE("identical input succeeds");

    int fcount = 3;
    IndexType f2n[] = {0, 1, 1, 2, 2, 0};
    IndexType f2noffsets[] = {0, 2, 4, 6};
    CellType f2ntypes[] = {SEGMENT, SEGMENT, SEGMENT};
    int stdfacecount = 3;
    IndexType stdFaceNodes[] = {0, 1, 1, 2, 2, 0};
    CellType stdFaceTypes[] = {SEGMENT, SEGMENT, SEGMENT};
    std::string errmesg;

    EXPECT_TRUE(verifyFaceNodesTypes(fcount,
                                     f2n,
                                     f2noffsets,
                                     f2ntypes,
                                     stdfacecount,
                                     stdFaceNodes,
                                     stdFaceTypes,
                                     errmesg));
    EXPECT_TRUE(errmesg.size() < 1);
  }

  {
    SCOPED_TRACE("different face order succeeds");

    int fcount = 5;
    IndexType f2n[] = {0, 3, 2, 1, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4};
    IndexType f2noffsets[] = {0, 4, 7, 10, 13, 16};
    CellType f2ntypes[] = {QUAD, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE};
    int stdfacecount = 5;
    IndexType stdFaceNodes[] = {0, 1, 4, 0, 3, 2, 1, 2, 3, 4, 1, 2, 4, 3, 0, 4};
    CellType stdFaceTypes[] = {TRIANGLE, QUAD, TRIANGLE, TRIANGLE, TRIANGLE};
    std::string errmesg;

    EXPECT_TRUE(verifyFaceNodesTypes(fcount,
                                     f2n,
                                     f2noffsets,
                                     f2ntypes,
                                     stdfacecount,
                                     stdFaceNodes,
                                     stdFaceTypes,
                                     errmesg));
    EXPECT_TRUE(errmesg.size() < 1);
  }

  {
    SCOPED_TRACE("rolled nodes in a face succeeds");

    int fcount = 4;
    IndexType f2n[] = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
    IndexType f2noffsets[] = {0, 3, 6, 9, 12};
    CellType f2ntypes[] = {TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE};
    int stdfacecount = 4;
    IndexType stdFaceNodes[] = {1, 0, 2, 0, 1, 3, 1, 2, 3, 0, 3, 2};
    CellType stdFaceTypes[] = {TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE};
    std::string errmesg;

    EXPECT_TRUE(verifyFaceNodesTypes(fcount,
                                     f2n,
                                     f2noffsets,
                                     f2ntypes,
                                     stdfacecount,
                                     stdFaceNodes,
                                     stdFaceTypes,
                                     errmesg));
    EXPECT_TRUE(errmesg.size() < 1);
  }

  {
    SCOPED_TRACE("different-order nodes in a face fails");

    int fcount = 4;
    IndexType f2n[] = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
    IndexType f2noffsets[] = {0, 3, 6, 9, 12};
    CellType f2ntypes[] = {TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE};
    int stdfacecount = 4;
    IndexType stdFaceNodes[] = {0, 1, 2, 0, 1, 3, 1, 2, 3, 2, 0, 3};
    CellType stdFaceTypes[] = {TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE};
    std::string errmesg;

    EXPECT_FALSE(verifyFaceNodesTypes(fcount,
                                      f2n,
                                      f2noffsets,
                                      f2ntypes,
                                      stdfacecount,
                                      stdFaceNodes,
                                      stdFaceTypes,
                                      errmesg));
    EXPECT_TRUE(errmesg.size() > 0);
  }
}

TEST(mint_mesh_face_relation, tf_faceMatches)
{
  using namespace axom::mint;

  std::vector<IndexType> list1 {2, 3, 4};
  std::vector<IndexType> list2 {3, 4, 2};
  std::vector<IndexType> list3 {3, 2, 4};
  std::vector<IndexType> list4 {3, 4};

  FaceTypeNodes a(SEGMENT, list1), b(TRIANGLE, list1), c(TRIANGLE, list1), d(SEGMENT, list2),
    e(SEGMENT, list3), f(TRIANGLE, list4);

  EXPECT_FALSE(faceMatches(a, b));  // mismatched type
  EXPECT_TRUE(faceMatches(b, c));   // matching types, identical lists
  EXPECT_TRUE(faceMatches(a, d));   // rotated list
  EXPECT_FALSE(faceMatches(a, e));  // wrong-order list
  EXPECT_FALSE(faceMatches(c, f));  // different-length list
}

/*! \brief Test driver for the mesh face relation construction.
 *
 * This TEST does the following:
 * -# generate mesh test sets, each supplied with correct answers, comprising
 *    -# the mesh itself
 *    -# the number of faces for each cell and the mesh's total number of
 *       faces
 *    -# the mesh's cell-to-cell neighbor relation
 *    -# the mesh's face types and nodes
 * -# for each test set, call runMeshFaceTest()
 * -# clean up.
 *
 * The function std::vector<internal::MeshFaceTest> generateFaceTestCases()
 * will produce mesh test sets corresponding to the figures in the file
 * mint_mesh_face_relation.svg along with the test answers.
 *
 * The function void runMeshFaceTest(internal::MeshFaceTest *) will call
 * internal::initFaces() on the test mesh, generate the neighbor relation,
 * and verify the face counts and neighbor relation as scoped tests.
 */
TEST(mint_mesh_face_relation, correct_construction)
{
  std::vector<axom::mint::internal::MeshFaceTest *> tests = axom::mint::generateFaceTestCases();

  for(auto t : tests)
  {
    axom::mint::runMeshFaceTest(t);
  }

  for(auto t : tests)
  {
    delete t;
  }
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
