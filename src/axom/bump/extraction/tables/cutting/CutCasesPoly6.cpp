// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
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

int numCutCasesPoly6 = 64;

int numCutShapesPoly6[] = {0, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3,
                           2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2,
                           2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 0};

int startCutShapesPoly6[] = {0,   0,   4,   8,   12,  16,  24,  28,  32,  36,  44,  52,  60,
                             64,  72,  76,  80,  84,  92,  100, 108, 116, 128, 136, 144, 148,
                             156, 164, 172, 176, 184, 188, 192, 196, 200, 208, 212, 220, 228,
                             236, 240, 248, 256, 268, 276, 284, 292, 300, 304, 308, 312, 320,
                             324, 332, 340, 348, 352, 356, 360, 368, 372, 376, 380, 384};

// clang-format off
unsigned char cutShapesPoly6[] = {
  // Case #0
  // Case #1
  ST_LIN, COLOR0, EA, EF,
  // Case #2
  ST_LIN, COLOR0, EA, EB,
  // Case #3
  ST_LIN, COLOR0, EB, EF,
  // Case #4
  ST_LIN, COLOR0, EB, EC,
  // Case #5
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EF,
  // Case #6
  ST_LIN, COLOR0, EA, EC,
  // Case #7
  ST_LIN, COLOR0, EC, EF,
  // Case #8
  ST_LIN, COLOR0, EC, ED,
  // Case #9
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #10
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EB, EC,
  // Case #11
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #12
  ST_LIN, COLOR0, EB, ED,
  // Case #13
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EF,
  // Case #14
  ST_LIN, COLOR0, EA, ED,
  // Case #15
  ST_LIN, COLOR0, ED, EF,
  // Case #16
  ST_LIN, COLOR0, ED, EE,
  // Case #17
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #18
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, ED,
  // Case #19
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #20
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #21
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #22
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #23
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #24
  ST_LIN, COLOR0, EC, EE,
  // Case #25
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #26
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, EC,
  // Case #27
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #28
  ST_LIN, COLOR0, EB, EE,
  // Case #29
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EE, EF,
  // Case #30
  ST_LIN, COLOR0, EA, EE,
  // Case #31
  ST_LIN, COLOR0, EE, EF,
  // Case #32
  ST_LIN, COLOR0, EE, EF,
  // Case #33
  ST_LIN, COLOR0, EA, EE,
  // Case #34
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EE,
  // Case #35
  ST_LIN, COLOR0, EB, EE,
  // Case #36
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, EC, EE,
  // Case #37
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  // Case #38
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EC, EE,
  // Case #39
  ST_LIN, COLOR0, EC, EE,
  // Case #40
  ST_LIN, COLOR0, EC, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #41
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #42
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #43
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #44
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #45
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  // Case #46
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #47
  ST_LIN, COLOR0, ED, EE,
  // Case #48
  ST_LIN, COLOR0, ED, EF,
  // Case #49
  ST_LIN, COLOR0, EA, ED,
  // Case #50
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, ED,
  // Case #51
  ST_LIN, COLOR0, EB, ED,
  // Case #52
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, EC, ED,
  // Case #53
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  // Case #54
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EC, ED,
  // Case #55
  ST_LIN, COLOR0, EC, ED,
  // Case #56
  ST_LIN, COLOR0, EC, EF,
  // Case #57
  ST_LIN, COLOR0, EA, EC,
  // Case #58
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EC,
  // Case #59
  ST_LIN, COLOR0, EB, EC,
  // Case #60
  ST_LIN, COLOR0, EB, EF,
  // Case #61
  ST_LIN, COLOR0, EA, EB,
  // Case #62
  ST_LIN, COLOR0, EA, EF
  // Case #63
};
// clang-format on

const size_t cutShapesPoly6Size = sizeof(cutShapesPoly6) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
