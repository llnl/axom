#include "CutCases.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{
namespace tables
{
namespace cutting
{

int numCutCasesWdg = 64;

// clang-format off
unsigned char cutShapesWdg[] = {
  // Case 0
  // Case 1
  ST_TRI,   COLOR0, EA, EG, EC,
  // Case 2
  ST_TRI,   COLOR0, EA, EB, EH,
  // Case 3
  ST_QUA,  COLOR0, EB, EH, EG, EC,
  // Case 4
  ST_TRI,   COLOR0, EB, EC, EI,
  // Case 5
  ST_QUA,  COLOR0, EB, EA, EG, EI,
  // Case 6
  ST_QUA,  COLOR0, EA, EC, EI, EH,
  // Case 7
  ST_TRI,   COLOR0, EH, EG, EI,
  // Case 8
  ST_TRI,   COLOR0, ED, EF, EG,
  // Case 9
  ST_QUA,  COLOR0, EA, ED, EF, EC,
  // Case 10
  ST_TRI,   COLOR0, EA, EB, EH,
  ST_TRI,   COLOR0, EG, ED, EF,
  // Case 11
  ST_POLY5, COLOR0, EF, EC, EB, EH, ED,
  // Case 12
  ST_TRI,   COLOR0, EC, EI, EB,
  ST_TRI,   COLOR0, EG, ED, EF,
  // Case 13
  ST_POLY5, COLOR0, EF, EI, EB, EA, ED,
  // Case 14
  ST_TRI,   COLOR0, EG, ED, EF,
  ST_QUA,  COLOR0, EI, EH, EA, EC,
  // Case 15
  ST_QUA,  COLOR0, EH, ED, EF, EI,
  // Case 16
  ST_TRI,   COLOR0, EH, EE, ED,
  // Case 17
  ST_TRI,   COLOR0, EH, EE, ED,
  ST_TRI,   COLOR0, EA, EG, EC,
  // Case 18
  ST_QUA,  COLOR0, ED, EA, EB, EE,
  // Case 19
  ST_POLY5, COLOR0, EG, EC, EB, EE, ED,
  // Case 20
  ST_TRI,   COLOR0, EH, EE, ED,
  ST_TRI,   COLOR0, EC, EI, EB,
  // Case 21
  ST_TRI,   COLOR0, EH, EE, ED,
  ST_QUA,  COLOR0, EB, EA, EG, EI,
  // Case 22
  ST_POLY5, COLOR0, EA, EC, EI, EE, ED,
  // Case 23
  ST_QUA,  COLOR0, ED, EG, EI, EE,
  // Case 24
  ST_QUA,  COLOR0, EG, EH, EE, EF,
  // Case 25
  ST_POLY5, COLOR0, EF, EC, EA, EH, EE,
  // Case 26
  ST_POLY5, COLOR0, EB, EE, EF, EG, EA,
  // Case 27
  ST_QUA,  COLOR0, EF, EC, EB, EE,
  // Case 28
  ST_TRI,   COLOR0, EC, EI, EB,
  ST_QUA,  COLOR0, EF, EG, EH, EE,
  // Case 29
  ST_POLY6, COLOR0, EA, EH, EE, EF, EI, EB,
  // Case 30
  ST_POLY6, COLOR0, EA, EC, EI, EE, EF, EG,
  // Case 31
  ST_TRI,   COLOR0, EI, EE, EF,
  // Case 32
  ST_TRI,   COLOR0, EE, EI, EF,
  // Case 33
  ST_TRI,   COLOR0, EE, EI, EF,
  ST_TRI,   COLOR0, EA, EG, EC,
  // Case 34
  ST_TRI,   COLOR0, EE, EI, EF,
  ST_TRI,   COLOR0, EA, EB, EH,
  // Case 35
  ST_TRI,   COLOR0, EE, EI, EF,
  ST_QUA,  COLOR0, EB, EH, EG, EC,
  // Case 36
  ST_QUA,  COLOR0, EF, EE, EB, EC,
  // Case 37
  ST_POLY5, COLOR0, EB, EA, EG, EF, EE,
  // Case 38
  ST_POLY5, COLOR0, EA, EC, EF, EE, EH,
  // Case 39
  ST_QUA,  COLOR0, EE, EH, EG, EF,
  // Case 40
  ST_QUA,  COLOR0, EI, EG, ED, EE,
  // Case 41
  ST_POLY5, COLOR0, EI, EC, EA, ED, EE,
  // Case 42
  ST_TRI,   COLOR0, EH, EA, EB,
  ST_QUA,  COLOR0, EG, ED, EE, EI,
  // Case 43
  ST_POLY6, COLOR0, EC, EB, EH, ED, EE, EI,
  // Case 44
  ST_POLY5, COLOR0, ED, EE, EB, EC, EG,
  // Case 45
  ST_QUA,  COLOR0, EB, EA, ED, EE,
  // Case 46
  ST_POLY6, COLOR0, EE, EH, EA, EC, EG, ED,
  // Case 47
  ST_TRI,   COLOR0, EH, ED, EE,
  // Case 48
  ST_QUA,  COLOR0, EH, EI, EF, ED,
  // Case 49
  ST_TRI,   COLOR0, EA, EG, EC,
  ST_QUA,  COLOR0, EH, EI, EF, ED,
  // Case 50
  ST_POLY5, COLOR0, EB, EI, EF, ED, EA,
  // Case 51
  ST_POLY6, COLOR0, ED, EG, EC, EB, EI, EF,
  // Case 52
  ST_POLY5, COLOR0, EB, EC, EF, ED, EH,
  // Case 53
  ST_POLY6, COLOR0, EB, EA, EG, EF, ED, EH,
  // Case 54
  ST_QUA,  COLOR0, EA, EC, EF, ED,
  // Case 55
  ST_TRI,   COLOR0, ED, EG, EF,
  // Case 56
  ST_TRI,   COLOR0, EH, EI, EG,
  // Case 57
  ST_QUA,  COLOR0, EA, EH, EI, EC,
  // Case 58
  ST_QUA,  COLOR0, EG, EA, EB, EI,
  // Case 59
  ST_TRI,   COLOR0, EC, EB, EI,
  // Case 60
  ST_QUA,  COLOR0, EG, EH, EB, EC,
  // Case 61
  ST_TRI,   COLOR0, EA, EH, EB,
  // Case 62
  ST_TRI,   COLOR0, EA, EC, EG
  // Case 63
};
// clang-format on

int numCutShapesWdg[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 2,
                         1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1,
                         1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};

int startCutShapesWdg[] = {0,   0,   5,   10,  16,  21,  27,  33,  38,  43,  49,  59,  66,
                           76,  83,  94,  100, 105, 115, 121, 128, 138, 149, 156, 162, 168,
                           175, 182, 188, 199, 207, 215, 220, 225, 235, 245, 256, 262, 269,
                           276, 282, 288, 295, 306, 314, 321, 327, 335, 340, 346, 357, 364,
                           372, 379, 387, 393, 398, 403, 409, 415, 420, 426, 431, 436};

const size_t cutShapesWdgSize = sizeof(cutShapesWdg) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
