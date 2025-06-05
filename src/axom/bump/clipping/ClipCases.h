// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef AXOM_VISIT_CLIP_CASES_H
#define AXOM_VISIT_CLIP_CASES_H
//---------------------------------------------------------------------------
// Axom modifications
// NOTE: The values for EA-EL and N0-N3 were reduced.
// NOTE: We're using AXOM_BUMP_EXPORT instead of VISIT_VTK_LIGHT_API througout
// clang-format off

#include "axom/export/bump.h"

#include <cstdlib>
namespace axom {
namespace bump {
namespace clipping {
namespace visit {
//---------------------------------------------------------------------------

// Programmer: Jeremy Meredith
// Date      : August 11, 2003
//
// Modifications:
//    Jeremy Meredith, Mon Sep 15 17:24:15 PDT 2003
//    Added NOCOLOR.
//
//    Jeremy Meredith, Thu Sep 18 11:29:12 PDT 2003
//    Added quad and triangle cases and output shapes.
//
//    Brad Whitlock, Tue Sep 23 09:59:23 PDT 2003
//    Added API so it builds on Windows.
//
//    Jeremy Meredith, Wed Jun 23 15:39:58 PDT 2004
//    Added voxel and pixel cases.  Not output shapes, though.
//
//    Jeremy Meredith, Tue Aug 29 13:52:33 EDT 2006
//    Added line segments and vertexes.
//

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
constexpr unsigned char ST_VTX = 106;
constexpr unsigned char ST_LIN = 107;
constexpr unsigned char ST_PNT = 108;

// Colors
constexpr unsigned char COLOR0 =  120;
constexpr unsigned char COLOR1 =  121;
constexpr unsigned char NOCOLOR = 122;

// Tables
extern AXOM_BUMP_EXPORT int numClipCasesHex;
extern AXOM_BUMP_EXPORT int numClipShapesHex[256];
extern AXOM_BUMP_EXPORT int startClipShapesHex[256];
extern AXOM_BUMP_EXPORT unsigned char clipShapesHex[];

extern AXOM_BUMP_EXPORT int numClipCasesVox;
extern AXOM_BUMP_EXPORT int numClipShapesVox[256];
extern AXOM_BUMP_EXPORT int startClipShapesVox[256];
extern AXOM_BUMP_EXPORT unsigned char clipShapesVox[];

extern AXOM_BUMP_EXPORT int numClipCasesWdg;
extern AXOM_BUMP_EXPORT int numClipShapesWdg[64];
extern AXOM_BUMP_EXPORT int startClipShapesWdg[64];
extern AXOM_BUMP_EXPORT unsigned char clipShapesWdg[];

extern AXOM_BUMP_EXPORT int numClipCasesPyr;
extern AXOM_BUMP_EXPORT int numClipShapesPyr[32];
extern AXOM_BUMP_EXPORT int startClipShapesPyr[32];
extern AXOM_BUMP_EXPORT unsigned char clipShapesPyr[];

extern AXOM_BUMP_EXPORT int numClipCasesTet;
extern AXOM_BUMP_EXPORT int numClipShapesTet[16];
extern AXOM_BUMP_EXPORT int startClipShapesTet[16];
extern AXOM_BUMP_EXPORT unsigned char clipShapesTet[];

extern AXOM_BUMP_EXPORT int numClipCasesQua;
extern AXOM_BUMP_EXPORT int numClipShapesQua[16];
extern AXOM_BUMP_EXPORT int startClipShapesQua[16];
extern AXOM_BUMP_EXPORT unsigned char clipShapesQua[];

extern AXOM_BUMP_EXPORT int numClipCasesPix;
extern AXOM_BUMP_EXPORT int numClipShapesPix[16];
extern AXOM_BUMP_EXPORT int startClipShapesPix[16];
extern AXOM_BUMP_EXPORT unsigned char clipShapesPix[];

extern AXOM_BUMP_EXPORT int numClipCasesTri;
extern AXOM_BUMP_EXPORT int numClipShapesTri[8];
extern AXOM_BUMP_EXPORT int startClipShapesTri[8];
extern AXOM_BUMP_EXPORT unsigned char clipShapesTri[];

extern AXOM_BUMP_EXPORT int numClipCasesLin;
extern AXOM_BUMP_EXPORT int numClipShapesLin[4];
extern AXOM_BUMP_EXPORT int startClipShapesLin[4];
extern AXOM_BUMP_EXPORT unsigned char clipShapesLin[];

extern AXOM_BUMP_EXPORT int numClipCasesVtx;
extern AXOM_BUMP_EXPORT int numClipShapesVtx[2];
extern AXOM_BUMP_EXPORT int startClipShapesVtx[2];
extern AXOM_BUMP_EXPORT unsigned char clipShapesVtx[];

extern AXOM_BUMP_EXPORT int numClipCasesPoly5;
extern AXOM_BUMP_EXPORT int numClipShapesPoly5[32];
extern AXOM_BUMP_EXPORT int startClipShapesPoly5[32];
extern AXOM_BUMP_EXPORT unsigned char clipShapesPoly5[];

extern AXOM_BUMP_EXPORT int numClipCasesPoly6;
extern AXOM_BUMP_EXPORT int numClipShapesPoly6[64];
extern AXOM_BUMP_EXPORT int startClipShapesPoly6[64];
extern AXOM_BUMP_EXPORT unsigned char clipShapesPoly6[];

extern AXOM_BUMP_EXPORT int numClipCasesPoly7;
extern AXOM_BUMP_EXPORT int numClipShapesPoly7[128];
extern AXOM_BUMP_EXPORT int startClipShapesPoly7[128];
extern AXOM_BUMP_EXPORT unsigned char clipShapesPoly7[];

extern AXOM_BUMP_EXPORT int numClipCasesPoly8;
extern AXOM_BUMP_EXPORT int numClipShapesPoly8[256];
extern AXOM_BUMP_EXPORT int startClipShapesPoly8[256];
extern AXOM_BUMP_EXPORT unsigned char clipShapesPoly8[];

//---------------------------------------------------------------------------
// Axom modifications
constexpr unsigned char ST_MIN = ST_TET;
constexpr unsigned char ST_MAX = (ST_PNT + 1);

extern AXOM_BUMP_EXPORT const size_t clipShapesTriSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesQuaSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesTetSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesPyrSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesWdgSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesHexSize;
} // namespace visit
} // namespace clipping
} // namespace bump
} // namespace axom
// clang-format on
//---------------------------------------------------------------------------

#endif
