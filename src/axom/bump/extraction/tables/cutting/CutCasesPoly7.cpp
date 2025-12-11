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

int numCutCasesPoly7 = 128;

int numCutShapesPoly7[] = {
  0, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1,
  1, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 1, 2, 2, 2, 2, 3, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1,
  1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2, 2, 2, 2, 1, 2, 2, 3, 2, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 1,
  1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 0};

int startCutShapesPoly7[] = {
  0,   0,   4,   8,   12,  16,  24,  28,  32,  36,  44,  52,  60,  64,  72,  76,  80,  84,  92,
  100, 108, 116, 128, 136, 144, 148, 156, 164, 172, 176, 184, 188, 192, 196, 204, 212, 220, 228,
  240, 248, 256, 264, 276, 288, 300, 308, 320, 328, 336, 340, 348, 356, 364, 372, 384, 392, 400,
  404, 412, 420, 428, 432, 440, 444, 448, 452, 456, 464, 468, 476, 484, 492, 496, 504, 512, 524,
  532, 540, 548, 556, 560, 568, 576, 588, 596, 608, 620, 632, 640, 648, 656, 668, 676, 684, 692,
  700, 704, 708, 712, 720, 724, 732, 740, 748, 752, 760, 768, 780, 788, 796, 804, 812, 816, 820,
  824, 832, 836, 844, 852, 860, 864, 868, 872, 880, 884, 888, 892, 896};

// clang-format off
unsigned char cutShapesPoly7[] = {
  // Case #0
  // Case #1
  ST_LIN, COLOR0, EA, EG,
  // Case #2
  ST_LIN, COLOR0, EA, EB,
  // Case #3
  ST_LIN, COLOR0, EB, EG,
  // Case #4
  ST_LIN, COLOR0, EB, EC,
  // Case #5
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EG,
  // Case #6
  ST_LIN, COLOR0, EA, EC,
  // Case #7
  ST_LIN, COLOR0, EC, EG,
  // Case #8
  ST_LIN, COLOR0, EC, ED,
  // Case #9
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EG,
  // Case #10
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EB, EC,
  // Case #11
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EG,
  // Case #12
  ST_LIN, COLOR0, EB, ED,
  // Case #13
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EG,
  // Case #14
  ST_LIN, COLOR0, EA, ED,
  // Case #15
  ST_LIN, COLOR0, ED, EG,
  // Case #16
  ST_LIN, COLOR0, ED, EE,
  // Case #17
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #18
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, ED,
  // Case #19
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #20
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #21
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #22
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #23
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #24
  ST_LIN, COLOR0, EC, EE,
  // Case #25
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EE, EG,
  // Case #26
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, EC,
  // Case #27
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EG,
  // Case #28
  ST_LIN, COLOR0, EB, EE,
  // Case #29
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EE, EG,
  // Case #30
  ST_LIN, COLOR0, EA, EE,
  // Case #31
  ST_LIN, COLOR0, EE, EG,
  // Case #32
  ST_LIN, COLOR0, EE, EF,
  // Case #33
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #34
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EE,
  // Case #35
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #36
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, EC, EE,
  // Case #37
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #38
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EC, EE,
  // Case #39
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #40
  ST_LIN, COLOR0, EC, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #41
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #42
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #43
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #44
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #45
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #46
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #47
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #48
  ST_LIN, COLOR0, ED, EF,
  // Case #49
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #50
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, ED,
  // Case #51
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #52
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, EC, ED,
  // Case #53
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #54
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EC, ED,
  // Case #55
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #56
  ST_LIN, COLOR0, EC, EF,
  // Case #57
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EF, EG,
  // Case #58
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EC,
  // Case #59
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EF, EG,
  // Case #60
  ST_LIN, COLOR0, EB, EF,
  // Case #61
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EF, EG,
  // Case #62
  ST_LIN, COLOR0, EA, EF,
  // Case #63
  ST_LIN, COLOR0, EF, EG,
  // Case #64
  ST_LIN, COLOR0, EF, EG,
  // Case #65
  ST_LIN, COLOR0, EA, EF,
  // Case #66
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EF,
  // Case #67
  ST_LIN, COLOR0, EB, EF,
  // Case #68
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, EF,
  // Case #69
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EF,
  // Case #70
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, EF,
  // Case #71
  ST_LIN, COLOR0, EC, EF,
  // Case #72
  ST_LIN, COLOR0, EC, EG,
  ST_LIN, COLOR0, ED, EF,
  // Case #73
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #74
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #75
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #76
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, ED, EF,
  // Case #77
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EF,
  // Case #78
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, ED, EF,
  // Case #79
  ST_LIN, COLOR0, ED, EF,
  // Case #80
  ST_LIN, COLOR0, ED, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #81
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #82
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #83
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #84
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #85
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #86
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #87
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #88
  ST_LIN, COLOR0, EC, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #89
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #90
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #91
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #92
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #93
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EE, EF,
  // Case #94
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #95
  ST_LIN, COLOR0, EE, EF,
  // Case #96
  ST_LIN, COLOR0, EE, EG,
  // Case #97
  ST_LIN, COLOR0, EA, EE,
  // Case #98
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EE,
  // Case #99
  ST_LIN, COLOR0, EB, EE,
  // Case #100
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, EE,
  // Case #101
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  // Case #102
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, EE,
  // Case #103
  ST_LIN, COLOR0, EC, EE,
  // Case #104
  ST_LIN, COLOR0, EC, EG,
  ST_LIN, COLOR0, ED, EE,
  // Case #105
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #106
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #107
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #108
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, ED, EE,
  // Case #109
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  // Case #110
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, ED, EE,
  // Case #111
  ST_LIN, COLOR0, ED, EE,
  // Case #112
  ST_LIN, COLOR0, ED, EG,
  // Case #113
  ST_LIN, COLOR0, EA, ED,
  // Case #114
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, ED,
  // Case #115
  ST_LIN, COLOR0, EB, ED,
  // Case #116
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, ED,
  // Case #117
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  // Case #118
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, ED,
  // Case #119
  ST_LIN, COLOR0, EC, ED,
  // Case #120
  ST_LIN, COLOR0, EC, EG,
  // Case #121
  ST_LIN, COLOR0, EA, EC,
  // Case #122
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  // Case #123
  ST_LIN, COLOR0, EB, EC,
  // Case #124
  ST_LIN, COLOR0, EB, EG,
  // Case #125
  ST_LIN, COLOR0, EA, EB,
  // Case #126
  ST_LIN, COLOR0, EA, EG
  // Case #127
};
// clang-format on

const size_t cutShapesPoly7Size = sizeof(cutShapesPoly7) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
