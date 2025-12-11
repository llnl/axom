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

int numCutCasesPoly8 = 256;

int numCutShapesPoly8[] = {
  0, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1,
  1, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 1, 2, 2, 2, 2, 3, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1,
  1, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 2, 3, 3, 3, 3, 4, 3, 3, 2, 3, 3, 3, 2, 3, 2, 2,
  1, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 1, 2, 2, 2, 2, 3, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1,
  1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2, 2, 2, 2, 1, 2, 2, 3, 2, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 1,
  2, 2, 3, 2, 3, 3, 3, 2, 3, 3, 4, 3, 3, 3, 3, 2, 2, 2, 3, 2, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 1,
  1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2, 2, 2, 2, 1, 2, 2, 3, 2, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 1,
  1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 0};

int startCutShapesPoly8[] = {
  0,    0,    4,    8,    12,   16,   24,   28,   32,   36,   44,   52,   60,   64,   72,   76,
  80,   84,   92,   100,  108,  116,  128,  136,  144,  148,  156,  164,  172,  176,  184,  188,
  192,  196,  204,  212,  220,  228,  240,  248,  256,  264,  276,  288,  300,  308,  320,  328,
  336,  340,  348,  356,  364,  372,  384,  392,  400,  404,  412,  420,  428,  432,  440,  444,
  448,  452,  460,  468,  476,  484,  496,  504,  512,  520,  532,  544,  556,  564,  576,  584,
  592,  600,  612,  624,  636,  648,  664,  676,  688,  696,  708,  720,  732,  740,  752,  760,
  768,  772,  780,  788,  796,  804,  816,  824,  832,  840,  852,  864,  876,  884,  896,  904,
  912,  916,  924,  932,  940,  948,  960,  968,  976,  980,  988,  996,  1004, 1008, 1016, 1020,
  1024, 1028, 1032, 1040, 1044, 1052, 1060, 1068, 1072, 1080, 1088, 1100, 1108, 1116, 1124, 1132,
  1136, 1144, 1152, 1164, 1172, 1184, 1196, 1208, 1216, 1224, 1232, 1244, 1252, 1260, 1268, 1276,
  1280, 1288, 1296, 1308, 1316, 1328, 1340, 1352, 1360, 1372, 1384, 1400, 1412, 1424, 1436, 1448,
  1456, 1464, 1472, 1484, 1492, 1504, 1516, 1528, 1536, 1544, 1552, 1564, 1572, 1580, 1588, 1596,
  1600, 1604, 1608, 1616, 1620, 1628, 1636, 1644, 1648, 1656, 1664, 1676, 1684, 1692, 1700, 1708,
  1712, 1720, 1728, 1740, 1748, 1760, 1772, 1784, 1792, 1800, 1808, 1820, 1828, 1836, 1844, 1852,
  1856, 1860, 1864, 1872, 1876, 1884, 1892, 1900, 1904, 1912, 1920, 1932, 1940, 1948, 1956, 1964,
  1968, 1972, 1976, 1984, 1988, 1996, 2004, 2012, 2016, 2020, 2024, 2032, 2036, 2040, 2044, 2048};

// clang-format off
unsigned char cutShapesPoly8[] = {
  // Case #0
  // Case #1
  ST_LIN, COLOR0, EA, EH,
  // Case #2
  ST_LIN, COLOR0, EA, EB,
  // Case #3
  ST_LIN, COLOR0, EB, EH,
  // Case #4
  ST_LIN, COLOR0, EB, EC,
  // Case #5
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EH,
  // Case #6
  ST_LIN, COLOR0, EA, EC,
  // Case #7
  ST_LIN, COLOR0, EC, EH,
  // Case #8
  ST_LIN, COLOR0, EC, ED,
  // Case #9
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EH,
  // Case #10
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EB, EC,
  // Case #11
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EH,
  // Case #12
  ST_LIN, COLOR0, EB, ED,
  // Case #13
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EH,
  // Case #14
  ST_LIN, COLOR0, EA, ED,
  // Case #15
  ST_LIN, COLOR0, ED, EH,
  // Case #16
  ST_LIN, COLOR0, ED, EE,
  // Case #17
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EE, EH,
  // Case #18
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, ED,
  // Case #19
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EH,
  // Case #20
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #21
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EH,
  // Case #22
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #23
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EH,
  // Case #24
  ST_LIN, COLOR0, EC, EE,
  // Case #25
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EE, EH,
  // Case #26
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, EC,
  // Case #27
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EH,
  // Case #28
  ST_LIN, COLOR0, EB, EE,
  // Case #29
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EE, EH,
  // Case #30
  ST_LIN, COLOR0, EA, EE,
  // Case #31
  ST_LIN, COLOR0, EE, EH,
  // Case #32
  ST_LIN, COLOR0, EE, EF,
  // Case #33
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #34
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EE,
  // Case #35
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #36
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, EC, EE,
  // Case #37
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #38
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EC, EE,
  // Case #39
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #40
  ST_LIN, COLOR0, EC, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #41
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #42
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #43
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #44
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #45
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #46
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, ED, EE,
  // Case #47
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EH,
  // Case #48
  ST_LIN, COLOR0, ED, EF,
  // Case #49
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EF, EH,
  // Case #50
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, ED,
  // Case #51
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EF, EH,
  // Case #52
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, EC, ED,
  // Case #53
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EH,
  // Case #54
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EC, ED,
  // Case #55
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EH,
  // Case #56
  ST_LIN, COLOR0, EC, EF,
  // Case #57
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EF, EH,
  // Case #58
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EB, EC,
  // Case #59
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EF, EH,
  // Case #60
  ST_LIN, COLOR0, EB, EF,
  // Case #61
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EF, EH,
  // Case #62
  ST_LIN, COLOR0, EA, EF,
  // Case #63
  ST_LIN, COLOR0, EF, EH,
  // Case #64
  ST_LIN, COLOR0, EF, EG,
  // Case #65
  ST_LIN, COLOR0, EA, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #66
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EF,
  // Case #67
  ST_LIN, COLOR0, EB, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #68
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, EF,
  // Case #69
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #70
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, EF,
  // Case #71
  ST_LIN, COLOR0, EC, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #72
  ST_LIN, COLOR0, EC, EG,
  ST_LIN, COLOR0, ED, EF,
  // Case #73
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #74
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #75
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #76
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, ED, EF,
  // Case #77
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #78
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, ED, EF,
  // Case #79
  ST_LIN, COLOR0, ED, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #80
  ST_LIN, COLOR0, ED, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #81
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #82
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #83
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #84
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #85
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #86
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #87
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #88
  ST_LIN, COLOR0, EC, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #89
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #90
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #91
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #92
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #93
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #94
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EE, EF,
  // Case #95
  ST_LIN, COLOR0, EE, EF,
  ST_LIN, COLOR0, EG, EH,
  // Case #96
  ST_LIN, COLOR0, EE, EG,
  // Case #97
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #98
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EE,
  // Case #99
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #100
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, EE,
  // Case #101
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #102
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, EE,
  // Case #103
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #104
  ST_LIN, COLOR0, EC, EG,
  ST_LIN, COLOR0, ED, EE,
  // Case #105
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #106
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #107
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #108
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, ED, EE,
  // Case #109
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #110
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, ED, EE,
  // Case #111
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EG, EH,
  // Case #112
  ST_LIN, COLOR0, ED, EG,
  // Case #113
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EG, EH,
  // Case #114
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, ED,
  // Case #115
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EG, EH,
  // Case #116
  ST_LIN, COLOR0, EB, EG,
  ST_LIN, COLOR0, EC, ED,
  // Case #117
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EG, EH,
  // Case #118
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EC, ED,
  // Case #119
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EG, EH,
  // Case #120
  ST_LIN, COLOR0, EC, EG,
  // Case #121
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EG, EH,
  // Case #122
  ST_LIN, COLOR0, EA, EG,
  ST_LIN, COLOR0, EB, EC,
  // Case #123
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EG, EH,
  // Case #124
  ST_LIN, COLOR0, EB, EG,
  // Case #125
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EG, EH,
  // Case #126
  ST_LIN, COLOR0, EA, EG,
  // Case #127
  ST_LIN, COLOR0, EG, EH,
  // Case #128
  ST_LIN, COLOR0, EG, EH,
  // Case #129
  ST_LIN, COLOR0, EA, EG,
  // Case #130
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EG,
  // Case #131
  ST_LIN, COLOR0, EB, EG,
  // Case #132
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, EG,
  // Case #133
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EG,
  // Case #134
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, EG,
  // Case #135
  ST_LIN, COLOR0, EC, EG,
  // Case #136
  ST_LIN, COLOR0, EC, EH,
  ST_LIN, COLOR0, ED, EG,
  // Case #137
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EG,
  // Case #138
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EG,
  // Case #139
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EG,
  // Case #140
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, ED, EG,
  // Case #141
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EG,
  // Case #142
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, ED, EG,
  // Case #143
  ST_LIN, COLOR0, ED, EG,
  // Case #144
  ST_LIN, COLOR0, ED, EH,
  ST_LIN, COLOR0, EE, EG,
  // Case #145
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #146
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #147
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #148
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #149
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #150
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #151
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EG,
  // Case #152
  ST_LIN, COLOR0, EC, EH,
  ST_LIN, COLOR0, EE, EG,
  // Case #153
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EE, EG,
  // Case #154
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EG,
  // Case #155
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EG,
  // Case #156
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EE, EG,
  // Case #157
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EE, EG,
  // Case #158
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EE, EG,
  // Case #159
  ST_LIN, COLOR0, EE, EG,
  // Case #160
  ST_LIN, COLOR0, EE, EH,
  ST_LIN, COLOR0, EF, EG,
  // Case #161
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #162
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #163
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #164
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #165
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #166
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #167
  ST_LIN, COLOR0, EC, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #168
  ST_LIN, COLOR0, EC, EH,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #169
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #170
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #171
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #172
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #173
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #174
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #175
  ST_LIN, COLOR0, ED, EE,
  ST_LIN, COLOR0, EF, EG,
  // Case #176
  ST_LIN, COLOR0, ED, EH,
  ST_LIN, COLOR0, EF, EG,
  // Case #177
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #178
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #179
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #180
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #181
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #182
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #183
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EF, EG,
  // Case #184
  ST_LIN, COLOR0, EC, EH,
  ST_LIN, COLOR0, EF, EG,
  // Case #185
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EF, EG,
  // Case #186
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EF, EG,
  // Case #187
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EF, EG,
  // Case #188
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EF, EG,
  // Case #189
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EF, EG,
  // Case #190
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EF, EG,
  // Case #191
  ST_LIN, COLOR0, EF, EG,
  // Case #192
  ST_LIN, COLOR0, EF, EH,
  // Case #193
  ST_LIN, COLOR0, EA, EF,
  // Case #194
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EF,
  // Case #195
  ST_LIN, COLOR0, EB, EF,
  // Case #196
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, EF,
  // Case #197
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EF,
  // Case #198
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, EF,
  // Case #199
  ST_LIN, COLOR0, EC, EF,
  // Case #200
  ST_LIN, COLOR0, EC, EH,
  ST_LIN, COLOR0, ED, EF,
  // Case #201
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #202
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #203
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EF,
  // Case #204
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, ED, EF,
  // Case #205
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EF,
  // Case #206
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, ED, EF,
  // Case #207
  ST_LIN, COLOR0, ED, EF,
  // Case #208
  ST_LIN, COLOR0, ED, EH,
  ST_LIN, COLOR0, EE, EF,
  // Case #209
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #210
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #211
  ST_LIN, COLOR0, EB, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #212
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #213
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #214
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #215
  ST_LIN, COLOR0, EC, ED,
  ST_LIN, COLOR0, EE, EF,
  // Case #216
  ST_LIN, COLOR0, EC, EH,
  ST_LIN, COLOR0, EE, EF,
  // Case #217
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #218
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #219
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, EE, EF,
  // Case #220
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EE, EF,
  // Case #221
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EE, EF,
  // Case #222
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EE, EF,
  // Case #223
  ST_LIN, COLOR0, EE, EF,
  // Case #224
  ST_LIN, COLOR0, EE, EH,
  // Case #225
  ST_LIN, COLOR0, EA, EE,
  // Case #226
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EE,
  // Case #227
  ST_LIN, COLOR0, EB, EE,
  // Case #228
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, EE,
  // Case #229
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  // Case #230
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, EE,
  // Case #231
  ST_LIN, COLOR0, EC, EE,
  // Case #232
  ST_LIN, COLOR0, EC, EH,
  ST_LIN, COLOR0, ED, EE,
  // Case #233
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #234
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #235
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #236
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, ED, EE,
  // Case #237
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  // Case #238
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, ED, EE,
  // Case #239
  ST_LIN, COLOR0, ED, EE,
  // Case #240
  ST_LIN, COLOR0, ED, EH,
  // Case #241
  ST_LIN, COLOR0, EA, ED,
  // Case #242
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, ED,
  // Case #243
  ST_LIN, COLOR0, EB, ED,
  // Case #244
  ST_LIN, COLOR0, EB, EH,
  ST_LIN, COLOR0, EC, ED,
  // Case #245
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  // Case #246
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EC, ED,
  // Case #247
  ST_LIN, COLOR0, EC, ED,
  // Case #248
  ST_LIN, COLOR0, EC, EH,
  // Case #249
  ST_LIN, COLOR0, EA, EC,
  // Case #250
  ST_LIN, COLOR0, EA, EH,
  ST_LIN, COLOR0, EB, EC,
  // Case #251
  ST_LIN, COLOR0, EB, EC,
  // Case #252
  ST_LIN, COLOR0, EB, EH,
  // Case #253
  ST_LIN, COLOR0, EA, EB,
  // Case #254
  ST_LIN, COLOR0, EA, EH
  // Case #255
};
// clang-format on

const size_t cutShapesPoly8Size = sizeof(cutShapesPoly8) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
