// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef AXOM_BUMP_EXTRACTION_CUT_CASES_H
#define AXOM_BUMP_EXTRACTION_CUT_CASES_H
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
namespace cutting
{

// Tables
extern AXOM_BUMP_EXPORT int numCutCasesHex;
extern AXOM_BUMP_EXPORT int numCutShapesHex[256];
extern AXOM_BUMP_EXPORT int startCutShapesHex[256];
extern AXOM_BUMP_EXPORT unsigned char cutShapesHex[];

extern AXOM_BUMP_EXPORT int numCutCasesWdg;
extern AXOM_BUMP_EXPORT int numCutShapesWdg[64];
extern AXOM_BUMP_EXPORT int startCutShapesWdg[64];
extern AXOM_BUMP_EXPORT unsigned char cutShapesWdg[];

extern AXOM_BUMP_EXPORT int numCutCasesPyr;
extern AXOM_BUMP_EXPORT int numCutShapesPyr[32];
extern AXOM_BUMP_EXPORT int startCutShapesPyr[32];
extern AXOM_BUMP_EXPORT unsigned char cutShapesPyr[];

extern AXOM_BUMP_EXPORT int numCutCasesTet;
extern AXOM_BUMP_EXPORT int numCutShapesTet[16];
extern AXOM_BUMP_EXPORT int startCutShapesTet[16];
extern AXOM_BUMP_EXPORT unsigned char cutShapesTet[];

extern AXOM_BUMP_EXPORT int numCutCasesQua;
extern AXOM_BUMP_EXPORT int numCutShapesQua[16];
extern AXOM_BUMP_EXPORT int startCutShapesQua[16];
extern AXOM_BUMP_EXPORT unsigned char cutShapesQua[];

extern AXOM_BUMP_EXPORT int numCutCasesTri;
extern AXOM_BUMP_EXPORT int numCutShapesTri[8];
extern AXOM_BUMP_EXPORT int startCutShapesTri[8];
extern AXOM_BUMP_EXPORT unsigned char cutShapesTri[];

extern AXOM_BUMP_EXPORT int numCutCasesPoly5;
extern AXOM_BUMP_EXPORT int numCutShapesPoly5[32];
extern AXOM_BUMP_EXPORT int startCutShapesPoly5[32];
extern AXOM_BUMP_EXPORT unsigned char cutShapesPoly5[];

extern AXOM_BUMP_EXPORT int numCutCasesPoly6;
extern AXOM_BUMP_EXPORT int numCutShapesPoly6[64];
extern AXOM_BUMP_EXPORT int startCutShapesPoly6[64];
extern AXOM_BUMP_EXPORT unsigned char cutShapesPoly6[];

extern AXOM_BUMP_EXPORT int numCutCasesPoly7;
extern AXOM_BUMP_EXPORT int numCutShapesPoly7[128];
extern AXOM_BUMP_EXPORT int startCutShapesPoly7[128];
extern AXOM_BUMP_EXPORT unsigned char cutShapesPoly7[];

extern AXOM_BUMP_EXPORT int numCutCasesPoly8;
extern AXOM_BUMP_EXPORT int numCutShapesPoly8[256];
extern AXOM_BUMP_EXPORT int startCutShapesPoly8[256];
extern AXOM_BUMP_EXPORT unsigned char cutShapesPoly8[];

extern AXOM_BUMP_EXPORT int numCutCasesLin;
extern AXOM_BUMP_EXPORT int numCutShapesLin[4];
extern AXOM_BUMP_EXPORT int startCutShapesLin[4];
extern AXOM_BUMP_EXPORT unsigned char cutShapesLin[];

extern AXOM_BUMP_EXPORT int numCutCasesVtx;
extern AXOM_BUMP_EXPORT int numCutShapesVtx[2];
extern AXOM_BUMP_EXPORT int startCutShapesVtx[2];
extern AXOM_BUMP_EXPORT unsigned char cutShapesVtx[];

extern AXOM_BUMP_EXPORT const size_t cutShapesTriSize;
extern AXOM_BUMP_EXPORT const size_t cutShapesQuaSize;
extern AXOM_BUMP_EXPORT const size_t cutShapesPoly5Size;
extern AXOM_BUMP_EXPORT const size_t cutShapesPoly6Size;
extern AXOM_BUMP_EXPORT const size_t cutShapesPoly7Size;
extern AXOM_BUMP_EXPORT const size_t cutShapesPoly8Size;
extern AXOM_BUMP_EXPORT const size_t cutShapesTetSize;
extern AXOM_BUMP_EXPORT const size_t cutShapesPyrSize;
extern AXOM_BUMP_EXPORT const size_t cutShapesWdgSize;
extern AXOM_BUMP_EXPORT const size_t cutShapesHexSize;

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom

#endif
