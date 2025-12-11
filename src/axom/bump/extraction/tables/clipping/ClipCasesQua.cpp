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

int numClipCasesQua = 16;

int numClipShapesQua[] = {1, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 1};

int startClipShapesQua[] = {0, 6, 18, 30, 42, 54, 72, 84, 96, 108, 120, 138, 150, 162, 174, 186};

// clang-format off
unsigned char clipShapesQua[] = {
  // Case #0
  ST_QUA, COLOR0, P0, P1, P2, P3,
  // Case #1
  ST_POLY5, COLOR0, EA, P1, P2, P3, ED,
  ST_TRI, COLOR1, P0, EA, ED,
  // Case #2
  ST_POLY5, COLOR0, P0, EA, EB, P2, P3,
  ST_TRI, COLOR1, EA, P1, EB,
  // Case #3
  ST_QUA, COLOR0, EB, P2, P3, ED,
  ST_QUA, COLOR1, P0, P1, EB, ED,
  // Case #4
  ST_POLY5, COLOR0, P0, P1, EB, EC, P3,
  ST_TRI, COLOR1, EB, P2, EC,
  // Case #5
  ST_TRI, COLOR0, EC, P3, ED,
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY6, COLOR1, P0, EA, EB, P2, EC, ED,
  // Case #6
  ST_QUA, COLOR0, P0, EA, EC, P3,
  ST_QUA, COLOR1, EA, P1, P2, EC,
  // Case #7
  ST_TRI, COLOR0, EC, P3, ED,
  ST_POLY5, COLOR1, P0, P1, P2, EC, ED,
  // Case #8
  ST_POLY5, COLOR0, P0, P1, P2, EC, ED,
  ST_TRI, COLOR1, EC, P3, ED,
  // Case #9
  ST_QUA, COLOR0, EA, P1, P2, EC,
  ST_QUA, COLOR1, P0, EA, EC, P3,
  // Case #10
  ST_TRI, COLOR0, P0, EA, ED,
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY6, COLOR1, EA, P1, EB, EC, P3, ED,
  // Case #11
  ST_TRI, COLOR0, EB, P2, EC,
  ST_POLY5, COLOR1, P0, P1, EB, EC, P3,
  // Case #12
  ST_QUA, COLOR0, P0, P1, EB, ED,
  ST_QUA, COLOR1, EB, P2, P3, ED,
  // Case #13
  ST_TRI, COLOR0, EA, P1, EB,
  ST_POLY5, COLOR1, P0, EA, EB, P2, P3,
  // Case #14
  ST_TRI, COLOR0, P0, EA, ED,
  ST_POLY5, COLOR1, EA, P1, P2, P3, ED,
  // Case #15
  ST_QUA, COLOR1, P0, P1, P2, P3
};
// clang-format on

const size_t clipShapesQuaSize = sizeof(clipShapesQua) / sizeof(unsigned char);

}  // namespace clipping
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
