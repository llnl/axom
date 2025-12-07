#include "CutCases.h"

namespace axom
{
namespace bump
{
namespace cutting
{
namespace tables
{

int numCutCasesTri = 8;

// clang-format off
unsigned char cutShapesTri[] = {
  // Case 0
  // Case 1
  ST_LIN,  COLOR0, EA, EC,
  // Case 2
  ST_LIN,  COLOR0, EA, EB,
  // Case 3
  ST_LIN,  COLOR0, EC, EB,
  // Case 4
  ST_LIN,  COLOR0, EC, EB,
  // Case 5
  ST_LIN,  COLOR0, EA, EB,
  // Case 6
  ST_LIN,  COLOR0, EA, EC
  // Case 7
};
// clang-format on

unsigned char numCutShapesTri[] = {0, 1, 1, 1, 1, 1, 1, 0};

unsigned char startCutShapesTri[] = {0, 0, 4, 8, 12, 16, 20, 24};

const size_t cutShapesTriSize = sizeof(cutShapesTri) / sizeof(unsigned char);

}  // namespace tables
}  // namespace cutting
}  // namespace bump
}  // namespace axom
