// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef AXOM_BUMP_EXTRACTION_CLIP_CASES_H
#define AXOM_BUMP_EXTRACTION_CLIP_CASES_H
#include "axom/export/bump.h"
#include "axom/bump/extraction/ExtractionConstants.hpp"
#include <cstdlib>

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

extern AXOM_BUMP_EXPORT const size_t clipShapesTriSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesQuaSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesPoly5Size;
extern AXOM_BUMP_EXPORT const size_t clipShapesPoly6Size;
extern AXOM_BUMP_EXPORT const size_t clipShapesPoly7Size;
extern AXOM_BUMP_EXPORT const size_t clipShapesPoly8Size;
extern AXOM_BUMP_EXPORT const size_t clipShapesTetSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesPyrSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesWdgSize;
extern AXOM_BUMP_EXPORT const size_t clipShapesHexSize;

}  // namespace clipping
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom

#endif
