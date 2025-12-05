// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "ClipCases.h"

namespace axom {
namespace bump {
namespace clipping {
namespace tables {

int numClipCasesPoly5 = 32;

int numClipShapesPoly5[] = {
  1, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 4, 2, 4, 2, 2, 2, 2, 3, 2, 3, 4, 4, 2, 2, 2, 4, 2, 2, 2, 2, 1
};

int startClipShapesPoly5[] = {
  0, 7, 20, 33, 46, 59, 78, 91, 104, 117, 136, 155, 178, 191, 214, 227, 240, 253, 266, 285, 298, 317, 340, 363, 376, 389, 402, 425, 438, 451, 464, 477
};

// clang-format off
unsigned char clipShapesPoly5[] = {
  // Case #0
  ST_POLY5, COLOR0, P0, P1, P2, P3, P4,
  // Case #1
  ST_POLY6, COLOR0, EA, P1, P2, P3, P4, EE,
  ST_TRI, COLOR1, P0, EA, EE,
  // Case #2
  ST_POLY6, COLOR0, P0, EA, EB, P2, P3, P4,
  ST_TRI, COLOR1, EA, P1, EB,
  // Case #3
  ST_POLY5, COLOR0, EB, P2, P3, P4, EE,
  ST_QUA, COLOR1, P0, P1, EB, EE,
  // Case #4
  ST_POLY6, COLOR0, P0, P1, EB, EC, P3, P4,
  ST_TRI, COLOR1, EB, P2, EC,
  // Case #5
  ST_TRI, COLOR0, EA, P1, EB,
  ST_QUA, COLOR0, EC, P3, P4, EE,
  ST_POLY6, COLOR1, P0, EA, EB, P2, EC, EE,
  // Case #6
  ST_POLY5, COLOR0, P0, EA, EC, P3, P4,
  ST_QUA, COLOR1, EA, P1, P2, EC,
  // Case #7
  ST_QUA, COLOR0, EC, P3, P4, EE,
  ST_POLY5, COLOR1, P0, P1, P2, EC, EE,
  // Case #8
  ST_POLY6, COLOR0, P0, P1, P2, EC, ED, P4,
  ST_TRI, COLOR1, EC, P3, ED,
  // Case #9
  ST_TRI, COLOR0, ED, P4, EE,
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_POLY6, COLOR1, P0, EA, EC, P3, ED, EE,
  // Case #10
  ST_TRI, COLOR0, EB, P2, EC,
  ST_QUA, COLOR0, P0, EA, ED, P4,
  ST_POLY6, COLOR1, EA, P1, EB, EC, P3, ED,
  // Case #11
  ST_TRI, COLOR0, EB, P2, EC,
  ST_TRI, COLOR0, ED, P4, EE,
  ST_POLY6, COLOR1, P0, P1, EB, EC, P3, ED,
  ST_TRI, COLOR1, P0, ED, EE,
  // Case #12
  ST_POLY5, COLOR0, P0, P1, EB, ED, P4,
  ST_QUA, COLOR1, EB, P2, P3, ED,
  // Case #13
  ST_TRI, COLOR0, ED, P4, EE,
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY6, COLOR1, P0, EA, P2, P3, ED, EE,
  ST_TRI, COLOR1, EA, EB, P2,
  // Case #14
  ST_QUA, COLOR0, P0, EA, ED, P4,
  ST_POLY5, COLOR1, EA, P1, P2, P3, ED,
  // Case #15
  ST_TRI, COLOR0, ED, P4, EE,
  ST_POLY6, COLOR1, P0, P1, P2, P3, ED, EE,
  // Case #16
  ST_POLY6, COLOR0, P0, P1, P2, P3, ED, EE,
  ST_TRI, COLOR1, ED, P4, EE,
  // Case #17
  ST_POLY5, COLOR0, EA, P1, P2, P3, ED,
  ST_QUA, COLOR1, P0, EA, ED, P4,
  // Case #18
  ST_TRI, COLOR0, P0, EA, EE,
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_POLY6, COLOR1, EA, P1, EB, ED, P4, EE,
  // Case #19
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_POLY5, COLOR1, P0, P1, EB, ED, P4,
  // Case #20
  ST_TRI, COLOR0, EC, P3, ED,
  ST_QUA, COLOR0, P0, P1, EB, EE,
  ST_POLY6, COLOR1, EB, P2, EC, ED, P4, EE,
  // Case #21
  ST_TRI, COLOR0, EA, P1, EB,
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY6, COLOR1, P0, EA, EB, P2, EC, P4,
  ST_TRI, COLOR1, EC, ED, P4,
  // Case #22
  ST_TRI, COLOR0, EC, P3, ED,
  ST_TRI, COLOR0, P0, EA, EE,
  ST_POLY6, COLOR1, P1, P2, EC, ED, P4, EE,
  ST_TRI, COLOR1, EA, P1, EE,
  // Case #23
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY6, COLOR1, P0, P1, P2, EC, ED, P4,
  // Case #24
  ST_POLY5, COLOR0, P0, P1, P2, EC, EE,
  ST_QUA, COLOR1, EC, P3, P4, EE,
  // Case #25
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_POLY5, COLOR1, P0, EA, EC, P3, P4,
  // Case #26
  ST_TRI, COLOR0, P0, EA, EE,
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY6, COLOR1, EA, P1, EB, P3, P4, EE,
  ST_TRI, COLOR1, EB, EC, P3,
  // Case #27
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY6, COLOR1, P0, P1, EB, EC, P3, P4,
  // Case #28
  ST_QUA, COLOR0, P0, P1, EB, EE,
  ST_POLY5, COLOR1, EB, P2, P3, P4, EE,
  // Case #29
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY6, COLOR1, P0, EA, EB, P2, P3, P4,
  // Case #30
  ST_TRI, COLOR0, P0, EA, EE,
  ST_POLY6, COLOR1, EA, P1, P2, P3, P4, EE,
  // Case #31
  ST_POLY5, COLOR1, P0, P1, P2, P3, P4
};
// clang-format on

const size_t clipShapesPoly5Size = sizeof(clipShapesPoly5) / sizeof(unsigned char);

} // namespace tables
} // namespace clipping
} // namespace bump
} // namespace axom
