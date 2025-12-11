// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "ClipCases.h"

namespace axom
{
namespace bump
{
namespace extraction
{
namespace tables
{
namespace clipping
{

int numClipCasesPoly6 = 64;

int numClipShapesPoly6[] = {1, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 2, 3, 3, 3, 3, 5,
                            3, 3, 2, 3, 3, 3, 2, 3, 2, 2, 2, 2, 3, 2, 3, 3, 3, 2, 3, 3, 5, 3,
                            3, 3, 3, 2, 2, 2, 3, 2, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 1};

int startClipShapesPoly6[] = {0,   8,   22,  36,  50,  64,  84,   98,   112,  126,  146,  166, 186,
                              200, 220, 234, 248, 262, 282, 302,  322,  342,  372,  392,  412, 426,
                              446, 466, 486, 500, 520, 534, 548,  562,  576,  596,  610,  630, 650,
                              670, 684, 704, 724, 754, 774, 794,  814,  834,  848,  862,  876, 896,
                              910, 930, 950, 970, 984, 998, 1012, 1032, 1046, 1060, 1074, 1088};

// clang-format off
unsigned char clipShapesPoly6[] = {
  // Case #0
  ST_POLY6, COLOR0, P0, P1, P2, P3, P4, P5,
  // Case #1
  ST_POLY7, COLOR0, EA, P1, P2, P3, P4, P5, EF,
  ST_TRI, COLOR1, P0, EA, EF,
  // Case #2
  ST_POLY7, COLOR0, P0, EA, EB, P2, P3, P4, P5,
  ST_TRI, COLOR1, EA, P1, EB,
  // Case #3
  ST_POLY6, COLOR0, EB, P2, P3, P4, P5, EF,
  ST_QUA, COLOR1, P0, P1, EB, EF,
  // Case #4
  ST_POLY7, COLOR0, P0, P1, EB, EC, P3, P4, P5,
  ST_TRI, COLOR1, EB, P2, EC,
  // Case #5
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY5, COLOR0, EC, P3, P4, P5, EF,
  ST_POLY6, COLOR1, P0, EA, EB, P2, EC, EF,
  // Case #6
  ST_POLY6, COLOR0, P0, EA, EC, P3, P4, P5,
  ST_QUA, COLOR1, EA, P1, P2, EC,
  // Case #7
  ST_POLY5, COLOR0, EC, P3, P4, P5, EF,
  ST_POLY5, COLOR1, P0, P1, P2, EC, EF,
  // Case #8
  ST_POLY7, COLOR0, P0, P1, P2, EC, ED, P4, P5,
  ST_TRI, COLOR1, EC, P3, ED,
  // Case #9
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_QUA, COLOR0, ED, P4, P5, EF,
  ST_POLY6, COLOR1, P0, EA, EC, P3, ED, EF,
  // Case #10
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY5, COLOR0, P0, EA, ED, P4, P5,
  ST_POLY6, COLOR1, EA, P1, EB, EC, P3, ED,
  // Case #11
  ST_TRI, COLOR0, EB, P2, EC,
  ST_QUA, COLOR0, ED, P4, P5, EF,
  ST_POLY7, COLOR1, P0, P1, EB, EC, P3, ED, EF,
  // Case #12
  ST_POLY6, COLOR0, P0, P1, EB, ED, P4, P5,
  ST_QUA, COLOR1, EB, P2, P3, ED,
  // Case #13
  ST_TRI, COLOR0, EA, P1, EB,
  ST_QUA, COLOR0, ED, P4, P5, EF,
  ST_POLY7, COLOR1, P0, EA, EB, P2, P3, ED, EF,
  // Case #14
  ST_POLY5, COLOR0, P0, EA, ED, P4, P5,
  ST_POLY5, COLOR1, EA, P1, P2, P3, ED,
  // Case #15
  ST_QUA, COLOR0, ED, P4, P5, EF,
  ST_POLY6, COLOR1, P0, P1, P2, P3, ED, EF,
  // Case #16
  ST_POLY7, COLOR0, P0, P1, P2, P3, ED, EE, P5,
  ST_TRI, COLOR1, ED, P4, EE,
  // Case #17
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY5, COLOR0, EA, P1, P2, P3, ED,
  ST_POLY6, COLOR1, P0, EA, ED, P4, EE, EF,
  // Case #18
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_QUA, COLOR0, P0, EA, EE, P5,
  ST_POLY6, COLOR1, EA, P1, EB, ED, P4, EE,
  // Case #19
  ST_TRI, COLOR0, EE, P5, EF,
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_POLY7, COLOR1, P0, P1, EB, ED, P4, EE, EF,
  // Case #20
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY5, COLOR0, P0, P1, EB, EE, P5,
  ST_POLY6, COLOR1, EB, P2, EC, ED, P4, EE,
  // Case #21
  ST_TRI, COLOR0, EA, P1, EB,
  ST_TRI, COLOR0, EC, P3, ED,
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY6, COLOR1, P0, EA, EB, P2, EC, ED,
  ST_POLY5, COLOR1, P0, ED, P4, EE, EF,
  // Case #22
  ST_TRI, COLOR0, EC, P3, ED,
  ST_QUA, COLOR0, P0, EA, EE, P5,
  ST_POLY7, COLOR1, EA, P1, P2, EC, ED, P4, EE,
  // Case #23
  ST_TRI, COLOR0, EC, P3, ED,
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY8, COLOR1, P0, P1, P2, EC, ED, P4, EE, EF,
  // Case #24
  ST_POLY6, COLOR0, P0, P1, P2, EC, EE, P5,
  ST_QUA, COLOR1, EC, P3, P4, EE,
  // Case #25
  ST_TRI, COLOR0, EE, P5, EF,
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_POLY7, COLOR1, P0, EA, EC, P3, P4, EE, EF,
  // Case #26
  ST_TRI, COLOR0, EB, P2, EC,
  ST_QUA, COLOR0, P0, EA, EE, P5,
  ST_POLY7, COLOR1, EA, P1, EB, EC, P3, P4, EE,
  // Case #27
  ST_TRI, COLOR0, EB, P2, EC,
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY8, COLOR1, P0, P1, EB, EC, P3, P4, EE, EF,
  // Case #28
  ST_POLY5, COLOR0, P0, P1, EB, EE, P5,
  ST_POLY5, COLOR1, EB, P2, P3, P4, EE,
  // Case #29
  ST_TRI, COLOR0, EE, P5, EF,
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY8, COLOR1, P0, EA, EB, P2, P3, P4, EE, EF,
  // Case #30
  ST_QUA, COLOR0, P0, EA, EE, P5,
  ST_POLY6, COLOR1, EA, P1, P2, P3, P4, EE,
  // Case #31
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY7, COLOR1, P0, P1, P2, P3, P4, EE, EF,
  // Case #32
  ST_POLY7, COLOR0, P0, P1, P2, P3, P4, EE, EF,
  ST_TRI, COLOR1, EE, P5, EF,
  // Case #33
  ST_POLY6, COLOR0, EA, P1, P2, P3, P4, EE,
  ST_QUA, COLOR1, P0, EA, EE, P5,
  // Case #34
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY5, COLOR0, EB, P2, P3, P4, EE,
  ST_POLY6, COLOR1, EA, P1, EB, EE, P5, EF,
  // Case #35
  ST_POLY5, COLOR0, EB, P2, P3, P4, EE,
  ST_POLY5, COLOR1, P0, P1, EB, EE, P5,
  // Case #36
  ST_QUA, COLOR0, EC, P3, P4, EE,
  ST_QUA, COLOR0, P0, P1, EB, EF,
  ST_POLY6, COLOR1, EB, P2, EC, EE, P5, EF,
  // Case #37
  ST_TRI, COLOR0, EA, P1, EB,
  ST_QUA, COLOR0, EC, P3, P4, EE,
  ST_POLY7, COLOR1, P0, EA, EB, P2, EC, EE, P5,
  // Case #38
  ST_TRI, COLOR0, P0, EA, EF,
  ST_QUA, COLOR0, EC, P3, P4, EE,
  ST_POLY7, COLOR1, EA, P1, P2, EC, EE, P5, EF,
  // Case #39
  ST_QUA, COLOR0, EC, P3, P4, EE,
  ST_POLY6, COLOR1, P0, P1, P2, EC, EE, P5,
  // Case #40
  ST_TRI, COLOR0, ED, P4, EE,
  ST_POLY5, COLOR0, P0, P1, P2, EC, EF,
  ST_POLY6, COLOR1, EC, P3, ED, EE, P5, EF,
  // Case #41
  ST_TRI, COLOR0, ED, P4, EE,
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_POLY7, COLOR1, P0, EA, EC, P3, ED, EE, P5,
  // Case #42
  ST_TRI, COLOR0, EB, P2, EC,
  ST_TRI, COLOR0, ED, P4, EE,
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY6, COLOR1, P1, EB, EC, P3, ED, EE,
  ST_POLY5, COLOR1, EA, P1, EE, P5, EF,
  // Case #43
  ST_TRI, COLOR0, EB, P2, EC,
  ST_TRI, COLOR0, ED, P4, EE,
  ST_POLY8, COLOR1, P0, P1, EB, EC, P3, ED, EE, P5,
  // Case #44
  ST_TRI, COLOR0, ED, P4, EE,
  ST_QUA, COLOR0, P0, P1, EB, EF,
  ST_POLY7, COLOR1, EB, P2, P3, ED, EE, P5, EF,
  // Case #45
  ST_TRI, COLOR0, ED, P4, EE,
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY8, COLOR1, P0, EA, EB, P2, P3, ED, EE, P5,
  // Case #46
  ST_TRI, COLOR0, ED, P4, EE,
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY8, COLOR1, EA, P1, P2, P3, ED, EE, P5, EF,
  // Case #47
  ST_TRI, COLOR0, ED, P4, EE,
  ST_POLY7, COLOR1, P0, P1, P2, P3, ED, EE, P5,
  // Case #48
  ST_POLY6, COLOR0, P0, P1, P2, P3, ED, EF,
  ST_QUA, COLOR1, ED, P4, P5, EF,
  // Case #49
  ST_POLY5, COLOR0, EA, P1, P2, P3, ED,
  ST_POLY5, COLOR1, P0, EA, ED, P4, P5,
  // Case #50
  ST_TRI, COLOR0, P0, EA, EF,
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_POLY7, COLOR1, EA, P1, EB, ED, P4, P5, EF,
  // Case #51
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_POLY6, COLOR1, P0, P1, EB, ED, P4, P5,
  // Case #52
  ST_TRI, COLOR0, EC, P3, ED,
  ST_QUA, COLOR0, P0, P1, EB, EF,
  ST_POLY7, COLOR1, EB, P2, EC, ED, P4, P5, EF,
  // Case #53
  ST_TRI, COLOR0, EA, P1, EB,
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY8, COLOR1, P0, EA, EB, P2, EC, ED, P4, P5,
  // Case #54
  ST_TRI, COLOR0, EC, P3, ED,
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY8, COLOR1, EA, P1, P2, EC, ED, P4, P5, EF,
  // Case #55
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY7, COLOR1, P0, P1, P2, EC, ED, P4, P5,
  // Case #56
  ST_POLY5, COLOR0, P0, P1, P2, EC, EF,
  ST_POLY5, COLOR1, EC, P3, P4, P5, EF,
  // Case #57
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_POLY6, COLOR1, P0, EA, EC, P3, P4, P5,
  // Case #58
  ST_TRI, COLOR0, P0, EA, EF,
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY8, COLOR1, EA, P1, EB, EC, P3, P4, P5, EF,
  // Case #59
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY7, COLOR1, P0, P1, EB, EC, P3, P4, P5,
  // Case #60
  ST_QUA, COLOR0, P0, P1, EB, EF,
  ST_POLY6, COLOR1, EB, P2, P3, P4, P5, EF,
  // Case #61
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY7, COLOR1, P0, EA, EB, P2, P3, P4, P5,
  // Case #62
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY7, COLOR1, EA, P1, P2, P3, P4, P5, EF,
  // Case #63
  ST_POLY6, COLOR1, P0, P1, P2, P3, P4, P5
};
// clang-format on

const size_t clipShapesPoly6Size = sizeof(clipShapesPoly6) / sizeof(unsigned char);

}  // namespace clipping
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
