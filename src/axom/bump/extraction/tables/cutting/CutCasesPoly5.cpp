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

int numCutCasesPoly5 = 32;

int numCutShapesPoly5[] = {0, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1,
                           1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 0};

int startCutShapesPoly5[] = {0,   0,   4,   8,   12,  16,  24,  28,  32,  36,  44,
                             52,  60,  64,  72,  76,  80,  84,  88,  96,  100, 108,
                             116, 124, 128, 132, 136, 144, 148, 152, 156, 160};

// clang-format off
unsigned char cutShapesPoly5[] = {
  // Case #0
  // Case #1
  ST_LIN, COLOR0, EA, EE,
  // Case #2
  ST_LIN, COLOR0, EA, EB,
  // Case #3
  ST_LIN, COLOR0, EB, EE,
  // Case #4
  ST_LIN, COLOR0, EB, EC,
  // Case #5
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, EE,
  // Case #6
  ST_LIN, COLOR0, EA, EC,
  // Case #7
  ST_LIN, COLOR0, EC, EE,
  // Case #8
  ST_LIN, COLOR0, EC, ED,
  // Case #9
  ST_LIN, COLOR0, EA, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #10
  ST_LIN, COLOR0, EA, ED,
  ST_LIN, COLOR0, EB, EC,
  // Case #11
  ST_LIN, COLOR0, EB, EC,
  ST_LIN, COLOR0, ED, EE,
  // Case #12
  ST_LIN, COLOR0, EB, ED,
  // Case #13
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, ED, EE,
  // Case #14
  ST_LIN, COLOR0, EA, ED,
  // Case #15
  ST_LIN, COLOR0, ED, EE,
  // Case #16
  ST_LIN, COLOR0, ED, EE,
  // Case #17
  ST_LIN, COLOR0, EA, ED,
  // Case #18
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, ED,
  // Case #19
  ST_LIN, COLOR0, EB, ED,
  // Case #20
  ST_LIN, COLOR0, EB, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #21
  ST_LIN, COLOR0, EA, EB,
  ST_LIN, COLOR0, EC, ED,
  // Case #22
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EC, ED,
  // Case #23
  ST_LIN, COLOR0, EC, ED,
  // Case #24
  ST_LIN, COLOR0, EC, EE,
  // Case #25
  ST_LIN, COLOR0, EA, EC,
  // Case #26
  ST_LIN, COLOR0, EA, EE,
  ST_LIN, COLOR0, EB, EC,
  // Case #27
  ST_LIN, COLOR0, EB, EC,
  // Case #28
  ST_LIN, COLOR0, EB, EE,
  // Case #29
  ST_LIN, COLOR0, EA, EB,
  // Case #30
  ST_LIN, COLOR0, EA, EE
  // Case #31
};
// clang-format on

const size_t cutShapesPoly5Size = sizeof(cutShapesPoly5) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
