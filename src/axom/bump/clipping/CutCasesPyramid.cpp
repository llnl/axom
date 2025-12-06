#include "CutCases.h"

namespace axom
{
namespace bump
{
namespace cutting
{
namespace tables
{

int numCutCasesPyramid = 32;

// clang-format off
unsigned char cutShapesPyramid[] = {
  // Case 0
  // Case 1
  ST_TRI,   COLOR0, EA, EE, ED,
  // Case 2
  ST_TRI,   COLOR0, EA, EB, EF,
  // Case 3
  ST_QUAD,  COLOR0, EB, EF, EE, ED,
  // Case 4
  ST_TRI,   COLOR0, EB, EC, EG,
  // Case 5
  ST_TRI,   COLOR0, EA, EE, ED,
  ST_TRI,   COLOR0, EB, EC, EG,
  // Case 6
  ST_QUAD,  COLOR0, EF, EA, EC, EG,
  // Case 7
  ST_POLY5, COLOR0, EC, EG, EF, EE, ED,
  // Case 8
  ST_TRI,   COLOR0, ED, EH, EC,
  // Case 9
  ST_QUAD,  COLOR0, EE, EH, EC, EA,
  // Case 10
  ST_TRI,   COLOR0, EA, EB, EF,
  ST_TRI,   COLOR0, ED, EH, EC,
  // Case 11
  ST_POLY5, COLOR0, EE, EH, EC, EB, EF,
  // Case 12
  ST_QUAD,  COLOR0, EG, EB, ED, EH,
  // Case 13
  ST_POLY5, COLOR0, EB, EA, EE, EH, EG,
  // Case 14
  ST_POLY5, COLOR0, ED, EH, EG, EF, EA,
  // Case 15
  ST_QUAD,  COLOR0, EF, EE, EH, EG,
  // Case 16
  ST_QUAD,  COLOR0, EH, EE, EF, EG,
  // Case 17
  ST_POLY5, COLOR0, EG, EH, ED, EA, EF,
  // Case 18
  ST_POLY5, COLOR0, EB, EG, EH, EE, EA,
  // Case 19
  ST_QUAD,  COLOR0, ED, EB, EG, EH,
  // Case 20
  ST_POLY5, COLOR0, EC, EH, EE, EF, EB,
  // Case 21
  ST_POLY6, COLOR0, EH, ED, EA, EF, EB, EC,
  // Case 22
  ST_QUAD,  COLOR0, EE, EA, EC, EH,
  // Case 23
  ST_TRI,   COLOR0, ED, EC, EH,
  // Case 24
  ST_POLY5, COLOR0, EF, EG, EC, ED, EE,
  // Case 25
  ST_QUAD,  COLOR0, EC, EA, EF, EG,
  // Case 26
  ST_POLY6, COLOR0, EE, EA, EB, EG, EC, ED,
  // Case 27
  ST_TRI,   COLOR0, EB, EG, EC,
  // Case 28
  ST_QUAD,  COLOR0, EB, ED, EE, EF,
  // Case 29
  ST_TRI,   COLOR0, EA, EF, EB,
  // Case 30
  ST_TRI,   COLOR0, EA, ED, EE,
  // Case 31
};
// clang-format on

unsigned char numCutShapesPyramid[] = {0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1,
                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};

unsigned char startCutShapesPyramid[] = {0,   0,   5,   10,  16,  21,  31,  37,  44,  49,  55,
                                         65,  72,  78,  85,  92,  98,  104, 111, 118, 124, 131,
                                         139, 145, 150, 157, 163, 171, 176, 182, 187, 192};

const size_t cutShapesPyramidSize = sizeof(cutShapesPyramid) / sizeof(unsigned char);

}  // namespace tables
}  // namespace cutting
}  // namespace bump
}  // namespace axom
