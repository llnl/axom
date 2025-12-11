// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef AXOM_BUMP_EXTRACTION_CONSTANTS_HPP
#define AXOM_BUMP_EXTRACTION_CONSTANTS_HPP
#include "axom/export/bump.h"

#include <cstdlib>
namespace axom
{
namespace bump
{
namespace extraction
{

// NOTE: These values were derived from VisIt

// Points of original cell (up to 8, for the hex)
// Note: we assume P0 is zero in several places.
// Note: we assume these values are contiguous and monotonic.
constexpr unsigned char P0 = 0;
constexpr unsigned char P1 = 1;
constexpr unsigned char P2 = 2;
constexpr unsigned char P3 = 3;
constexpr unsigned char P4 = 4;
constexpr unsigned char P5 = 5;
constexpr unsigned char P6 = 6;
constexpr unsigned char P7 = 7;

// Edges of original cell (up to 12, for the hex)
// Note: we assume these values are contiguous and monotonic.
constexpr unsigned char EA = 8;
constexpr unsigned char EB = 9;
constexpr unsigned char EC = 10;
constexpr unsigned char ED = 11;
constexpr unsigned char EE = 12;
constexpr unsigned char EF = 13;
constexpr unsigned char EG = 14;
constexpr unsigned char EH = 15;
constexpr unsigned char EI = 16;
constexpr unsigned char EJ = 17;
constexpr unsigned char EK = 18;
constexpr unsigned char EL = 19;

// New interpolated points (ST_PNT outputs)
// Note: we assume these values are contiguous and monotonic.
constexpr unsigned char N0 = 20;
constexpr unsigned char N1 = 21;
constexpr unsigned char N2 = 22;
constexpr unsigned char N3 = 23;

// Shapes
constexpr unsigned char ST_TET = 100;
constexpr unsigned char ST_PYR = 101;
constexpr unsigned char ST_WDG = 102;
constexpr unsigned char ST_HEX = 103;
constexpr unsigned char ST_TRI = 104;
constexpr unsigned char ST_QUA = 105;
constexpr unsigned char ST_POLY5 = 106;
constexpr unsigned char ST_POLY6 = 107;
constexpr unsigned char ST_POLY7 = 108;
constexpr unsigned char ST_POLY8 = 109;
constexpr unsigned char ST_VTX = 110;
constexpr unsigned char ST_LIN = 111;
constexpr unsigned char ST_PNT = 112;

constexpr unsigned char ST_MIN = ST_TET;
constexpr unsigned char ST_MAX = (ST_PNT + 1);

// Colors
constexpr unsigned char COLOR0 = 120;
constexpr unsigned char COLOR1 = 121;
constexpr unsigned char NOCOLOR = 122;

}  // namespace extraction
}  // namespace bump
}  // namespace axom
// clang-format on
//---------------------------------------------------------------------------

#endif
