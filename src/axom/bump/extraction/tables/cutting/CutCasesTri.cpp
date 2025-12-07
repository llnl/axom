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

int numCutShapesTri[] = {0, 1, 1, 1, 1, 1, 1, 0};

int startCutShapesTri[] = {0, 0, 4, 8, 12, 16, 20, 24};

const size_t cutShapesTriSize = sizeof(cutShapesTri) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
