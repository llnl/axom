// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef AXOM_VISIT_CLIP_CASES_H
#define AXOM_VISIT_CLIP_CASES_H
//---------------------------------------------------------------------------
// Axom modifications
// NOTE: The values for EA-EL and N0-N3 were reduced.
// NOTE: We're using AXOM_MIR_EXPORT instead of VISIT_VTK_LIGHT_API througout
// clang-format off

#include "axom/export/mir.h"

#include <cstdlib>
namespace axom {
namespace mir {
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
#define P0     0
#define P1     1
#define P2     2
#define P3     3
#define P4     4
#define P5     5
#define P6     6
#define P7     7

// Edges of original cell (up to 12, for the hex)
// Note: we assume these values are contiguous and monotonic.
#define EA     8
#define EB     9
#define EC     10
#define ED     11
#define EE     12
#define EF     13
#define EG     14
#define EH     15
#define EI     16
#define EJ     17
#define EK     18
#define EL     19

// New interpolated points (ST_PNT outputs)
// Note: we assume these values are contiguous and monotonic.
#define N0     20
#define N1     21
#define N2     22
#define N3     23

// Shapes
#define ST_TET 100
#define ST_PYR 101
#define ST_WDG 102
#define ST_HEX 103
#define ST_TRI 104
#define ST_QUA 105
#define ST_VTX 106
#define ST_LIN 107
#define ST_PNT 108

// Colors
#define COLOR0  120
#define COLOR1  121
#define NOCOLOR 122

// Tables
extern AXOM_MIR_EXPORT int numClipCasesHex;
extern AXOM_MIR_EXPORT int numClipShapesHex[256];
extern AXOM_MIR_EXPORT int startClipShapesHex[256];
extern AXOM_MIR_EXPORT unsigned char clipShapesHex[];

extern AXOM_MIR_EXPORT int numClipCasesVox;
extern AXOM_MIR_EXPORT int numClipShapesVox[256];
extern AXOM_MIR_EXPORT int startClipShapesVox[256];
extern AXOM_MIR_EXPORT unsigned char clipShapesVox[];

extern AXOM_MIR_EXPORT int numClipCasesWdg;
extern AXOM_MIR_EXPORT int numClipShapesWdg[64];
extern AXOM_MIR_EXPORT int startClipShapesWdg[64];
extern AXOM_MIR_EXPORT unsigned char clipShapesWdg[];

extern AXOM_MIR_EXPORT int numClipCasesPyr;
extern AXOM_MIR_EXPORT int numClipShapesPyr[32];
extern AXOM_MIR_EXPORT int startClipShapesPyr[32];
extern AXOM_MIR_EXPORT unsigned char clipShapesPyr[];

extern AXOM_MIR_EXPORT int numClipCasesTet;
extern AXOM_MIR_EXPORT int numClipShapesTet[16];
extern AXOM_MIR_EXPORT int startClipShapesTet[16];
extern AXOM_MIR_EXPORT unsigned char clipShapesTet[];

extern AXOM_MIR_EXPORT int numClipCasesQua;
extern AXOM_MIR_EXPORT int numClipShapesQua[16];
extern AXOM_MIR_EXPORT int startClipShapesQua[16];
extern AXOM_MIR_EXPORT unsigned char clipShapesQua[];

extern AXOM_MIR_EXPORT int numClipCasesPix;
extern AXOM_MIR_EXPORT int numClipShapesPix[16];
extern AXOM_MIR_EXPORT int startClipShapesPix[16];
extern AXOM_MIR_EXPORT unsigned char clipShapesPix[];

extern AXOM_MIR_EXPORT int numClipCasesTri;
extern AXOM_MIR_EXPORT int numClipShapesTri[8];
extern AXOM_MIR_EXPORT int startClipShapesTri[8];
extern AXOM_MIR_EXPORT unsigned char clipShapesTri[];

extern AXOM_MIR_EXPORT int numClipCasesLin;
extern AXOM_MIR_EXPORT int numClipShapesLin[4];
extern AXOM_MIR_EXPORT int startClipShapesLin[4];
extern AXOM_MIR_EXPORT unsigned char clipShapesLin[];

extern AXOM_MIR_EXPORT int numClipCasesVtx;
extern AXOM_MIR_EXPORT int numClipShapesVtx[2];
extern AXOM_MIR_EXPORT int startClipShapesVtx[2];
extern AXOM_MIR_EXPORT unsigned char clipShapesVtx[];

extern AXOM_MIR_EXPORT int numClipCasesPoly5;
extern AXOM_MIR_EXPORT int numClipShapesPoly5[32];
extern AXOM_MIR_EXPORT int startClipShapesPoly5[32];
extern AXOM_MIR_EXPORT unsigned char clipShapesPoly5[];

extern AXOM_MIR_EXPORT int numClipCasesPoly6;
extern AXOM_MIR_EXPORT int numClipShapesPoly6[64];
extern AXOM_MIR_EXPORT int startClipShapesPoly6[64];
extern AXOM_MIR_EXPORT unsigned char clipShapesPoly6[];

extern AXOM_MIR_EXPORT int numClipCasesPoly7;
extern AXOM_MIR_EXPORT int numClipShapesPoly7[128];
extern AXOM_MIR_EXPORT int startClipShapesPoly7[128];
extern AXOM_MIR_EXPORT unsigned char clipShapesPoly7[];

extern AXOM_MIR_EXPORT int numClipCasesPoly8;
extern AXOM_MIR_EXPORT int numClipShapesPoly8[256];
extern AXOM_MIR_EXPORT int startClipShapesPoly8[256];
extern AXOM_MIR_EXPORT unsigned char clipShapesPoly8[];

//---------------------------------------------------------------------------
// Axom modifications
#define ST_MIN ST_TET
#define ST_MAX (ST_PNT + 1)

extern AXOM_MIR_EXPORT const size_t clipShapesTriSize;
extern AXOM_MIR_EXPORT const size_t clipShapesQuaSize;
extern AXOM_MIR_EXPORT const size_t clipShapesTetSize;
extern AXOM_MIR_EXPORT const size_t clipShapesPyrSize;
extern AXOM_MIR_EXPORT const size_t clipShapesWdgSize;
extern AXOM_MIR_EXPORT const size_t clipShapesHexSize;
} // namespace visit
} // namespace clipping
} // namespace mir
} // namespace axom
// clang-format on
//---------------------------------------------------------------------------

#endif
