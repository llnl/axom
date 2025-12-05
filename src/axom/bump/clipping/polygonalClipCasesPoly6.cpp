// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "ClipCases.h"

namespace axom {
namespace bump {
namespace clipping {
namespace tables {

int numClipCasesPoly6 = 64;

int numClipShapesPoly6[] = {
  1, 3, 3, 2, 3, 3, 2, 2, 3, 3, 3, 4, 2, 4, 2, 2, 3, 3, 3, 4, 3, 5, 4, 4, 2, 4, 4, 4, 2, 4, 2, 3, 3, 2, 3, 2, 3, 4, 4, 2, 3, 4, 5, 4, 4, 4, 4, 3, 2, 2, 4, 2, 4, 4, 4, 3, 2, 2, 4, 3, 2, 3, 3, 1
};

int startClipShapesPoly6[] = {
  0, 8, 26, 44, 58, 76, 96, 110, 124, 142, 162, 182, 206, 220, 244, 258, 272, 290, 310, 330, 354, 374, 404, 428, 452, 466, 490, 514, 538, 552, 576, 590, 608, 626, 640, 660, 674, 694, 718, 742, 756, 776, 800, 830, 854, 878, 902, 926, 944, 958, 972, 996, 1010, 1034, 1058, 1082, 1100, 1114, 1128, 1152, 1170, 1184, 1202, 1220
};

// clang-format off
unsigned char clipShapesPoly6[] = {
  // Case #0
  ST_POLY6, COLOR0, P0, P1, P2, P3, P4, P5,
  // Case #1
  ST_POLY6, COLOR0, EA, P1, P2, P3, P4, P5,
  ST_TRI, COLOR0, EA, P5, EF,
  ST_TRI, COLOR1, P0, EA, EF,
  // Case #2
  ST_POLY6, COLOR0, P0, EB, P2, P3, P4, P5,
  ST_TRI, COLOR0, P0, EA, EB,
  ST_TRI, COLOR1, EA, P1, EB,
  // Case #3
  ST_POLY6, COLOR0, EB, P2, P3, P4, P5, EF,
  ST_QUA, COLOR1, P0, P1, EB, EF,
  // Case #4
  ST_POLY6, COLOR0, P0, P1, EC, P3, P4, P5,
  ST_TRI, COLOR0, P1, EB, EC,
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
  ST_POLY6, COLOR0, P0, P1, P2, ED, P4, P5,
  ST_TRI, COLOR0, P2, EC, ED,
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
  ST_POLY6, COLOR1, P0, P1, EB, EC, P3, ED,
  ST_TRI, COLOR1, P0, ED, EF,
  // Case #12
  ST_POLY6, COLOR0, P0, P1, EB, ED, P4, P5,
  ST_QUA, COLOR1, EB, P2, P3, ED,
  // Case #13
  ST_TRI, COLOR0, EA, P1, EB,
  ST_QUA, COLOR0, ED, P4, P5, EF,
  ST_POLY6, COLOR1, P0, EA, EB, P2, P3, EF,
  ST_TRI, COLOR1, P3, ED, EF,
  // Case #14
  ST_POLY5, COLOR0, P0, EA, ED, P4, P5,
  ST_POLY5, COLOR1, EA, P1, P2, P3, ED,
  // Case #15
  ST_QUA, COLOR0, ED, P4, P5, EF,
  ST_POLY6, COLOR1, P0, P1, P2, P3, ED, EF,
  // Case #16
  ST_POLY6, COLOR0, P0, P1, P2, P3, EE, P5,
  ST_TRI, COLOR0, P3, ED, EE,
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
  ST_POLY6, COLOR1, P0, P1, ED, P4, EE, EF,
  ST_TRI, COLOR1, P1, EB, ED,
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
  ST_POLY6, COLOR1, P1, P2, EC, ED, P4, EE,
  ST_TRI, COLOR1, EA, P1, EE,
  // Case #23
  ST_TRI, COLOR0, EC, P3, ED,
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY6, COLOR1, P0, P1, P2, EC, ED, P4,
  ST_QUA, COLOR1, P0, P4, EE, EF,
  // Case #24
  ST_POLY6, COLOR0, P0, P1, P2, EC, EE, P5,
  ST_QUA, COLOR1, EC, P3, P4, EE,
  // Case #25
  ST_TRI, COLOR0, EE, P5, EF,
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_POLY6, COLOR1, P0, EA, P3, P4, EE, EF,
  ST_TRI, COLOR1, EA, EC, P3,
  // Case #26
  ST_TRI, COLOR0, EB, P2, EC,
  ST_QUA, COLOR0, P0, EA, EE, P5,
  ST_POLY6, COLOR1, EA, P1, EB, EC, P3, P4,
  ST_TRI, COLOR1, EA, P4, EE,
  // Case #27
  ST_TRI, COLOR0, EB, P2, EC,
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY6, COLOR1, P0, P1, EB, EC, P3, P4,
  ST_QUA, COLOR1, P0, P4, EE, EF,
  // Case #28
  ST_POLY5, COLOR0, P0, P1, EB, EE, P5,
  ST_POLY5, COLOR1, EB, P2, P3, P4, EE,
  // Case #29
  ST_TRI, COLOR0, EE, P5, EF,
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY6, COLOR1, P0, P2, P3, P4, EE, EF,
  ST_QUA, COLOR1, P0, EA, EB, P2,
  // Case #30
  ST_QUA, COLOR0, P0, EA, EE, P5,
  ST_POLY6, COLOR1, EA, P1, P2, P3, P4, EE,
  // Case #31
  ST_TRI, COLOR0, EE, P5, EF,
  ST_POLY6, COLOR1, P0, P1, P2, P3, P4, EE,
  ST_TRI, COLOR1, P0, EE, EF,
  // Case #32
  ST_POLY6, COLOR0, P0, P1, P2, P3, P4, EF,
  ST_TRI, COLOR0, P4, EE, EF,
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
  ST_POLY6, COLOR1, P0, EA, EB, P2, EC, P5,
  ST_TRI, COLOR1, EC, EE, P5,
  // Case #38
  ST_TRI, COLOR0, P0, EA, EF,
  ST_QUA, COLOR0, EC, P3, P4, EE,
  ST_POLY6, COLOR1, EA, P1, P2, EE, P5, EF,
  ST_TRI, COLOR1, P2, EC, EE,
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
  ST_POLY6, COLOR1, P0, EC, P3, ED, EE, P5,
  ST_TRI, COLOR1, P0, EA, EC,
  // Case #42
  ST_TRI, COLOR0, EB, P2, EC,
  ST_TRI, COLOR0, ED, P4, EE,
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY6, COLOR1, P1, EB, EC, P3, ED, EE,
  ST_POLY5, COLOR1, EA, P1, EE, P5, EF,
  // Case #43
  ST_TRI, COLOR0, EB, P2, EC,
  ST_TRI, COLOR0, ED, P4, EE,
  ST_POLY6, COLOR1, P0, P1, EB, EC, P3, P5,
  ST_QUA, COLOR1, P3, ED, EE, P5,
  // Case #44
  ST_TRI, COLOR0, ED, P4, EE,
  ST_QUA, COLOR0, P0, P1, EB, EF,
  ST_POLY6, COLOR1, P2, P3, ED, EE, P5, EF,
  ST_TRI, COLOR1, EB, P2, EF,
  // Case #45
  ST_TRI, COLOR0, ED, P4, EE,
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY6, COLOR1, P0, P2, P3, ED, EE, P5,
  ST_QUA, COLOR1, P0, EA, EB, P2,
  // Case #46
  ST_TRI, COLOR0, ED, P4, EE,
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY6, COLOR1, P1, P2, P3, ED, EE, P5,
  ST_QUA, COLOR1, EA, P1, P5, EF,
  // Case #47
  ST_TRI, COLOR0, ED, P4, EE,
  ST_POLY6, COLOR1, P0, P1, P2, P3, ED, P5,
  ST_TRI, COLOR1, ED, EE, P5,
  // Case #48
  ST_POLY6, COLOR0, P0, P1, P2, P3, ED, EF,
  ST_QUA, COLOR1, ED, P4, P5, EF,
  // Case #49
  ST_POLY5, COLOR0, EA, P1, P2, P3, ED,
  ST_POLY5, COLOR1, P0, EA, ED, P4, P5,
  // Case #50
  ST_TRI, COLOR0, P0, EA, EF,
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_POLY6, COLOR1, EA, P1, EB, P4, P5, EF,
  ST_TRI, COLOR1, EB, ED, P4,
  // Case #51
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_POLY6, COLOR1, P0, P1, EB, ED, P4, P5,
  // Case #52
  ST_TRI, COLOR0, EC, P3, ED,
  ST_QUA, COLOR0, P0, P1, EB, EF,
  ST_POLY6, COLOR1, EB, P2, EC, ED, P4, P5,
  ST_TRI, COLOR1, EB, P5, EF,
  // Case #53
  ST_TRI, COLOR0, EA, P1, EB,
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY6, COLOR1, P0, EA, EB, P2, P4, P5,
  ST_QUA, COLOR1, P2, EC, ED, P4,
  // Case #54
  ST_TRI, COLOR0, EC, P3, ED,
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY6, COLOR1, P1, P2, EC, ED, P4, P5,
  ST_QUA, COLOR1, EA, P1, P5, EF,
  // Case #55
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY6, COLOR1, P0, P1, P2, EC, P4, P5,
  ST_TRI, COLOR1, EC, ED, P4,
  // Case #56
  ST_POLY5, COLOR0, P0, P1, P2, EC, EF,
  ST_POLY5, COLOR1, EC, P3, P4, P5, EF,
  // Case #57
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_POLY6, COLOR1, P0, EA, EC, P3, P4, P5,
  // Case #58
  ST_TRI, COLOR0, P0, EA, EF,
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY6, COLOR1, EA, P1, P3, P4, P5, EF,
  ST_QUA, COLOR1, P1, EB, EC, P3,
  // Case #59
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY6, COLOR1, P0, P1, EB, P3, P4, P5,
  ST_TRI, COLOR1, EB, EC, P3,
  // Case #60
  ST_QUA, COLOR0, P0, P1, EB, EF,
  ST_POLY6, COLOR1, EB, P2, P3, P4, P5, EF,
  // Case #61
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY6, COLOR1, P0, EA, P2, P3, P4, P5,
  ST_TRI, COLOR1, EA, EB, P2,
  // Case #62
  ST_TRI, COLOR0, P0, EA, EF,
  ST_POLY6, COLOR1, P1, P2, P3, P4, P5, EF,
  ST_TRI, COLOR1, EA, P1, EF,
  // Case #63
  ST_POLY6, COLOR1, P0, P1, P2, P3, P4, P5
};
// clang-format on

const size_t clipShapesPoly6Size = sizeof(clipShapesPoly6) / sizeof(unsigned char);

} // namespace tables
} // namespace clipping
} // namespace bump
} // namespace axom
