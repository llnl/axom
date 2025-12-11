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

int numCutCasesQua = 16;

// clang-format off
unsigned char cutShapesQua[] = {
  // Case 0
  // Case 1
  ST_LIN,  COLOR0, EA, ED,
  // Case 2
  ST_LIN,  COLOR0, EA, EB,
  // Case 3
  ST_LIN,  COLOR0, EB, ED,
  // Case 4
  ST_LIN,  COLOR0, EB, EC,
  // Case 5
  ST_LIN,  COLOR0, EA, ED,
  ST_LIN,  COLOR0, EB, EC,
  // Case 6
  ST_LIN,  COLOR0, EA, EC,
  // Case 7
  ST_LIN,  COLOR0, ED, EC,
  // Case 8
  ST_LIN,  COLOR0, EC, ED,
  // Case 9
  ST_LIN,  COLOR0, EA, EC,
  // Case 10
  ST_LIN,  COLOR0, EA, EB,
  ST_LIN,  COLOR0, ED, EC,
  // Case 11
  ST_LIN,  COLOR0, EB, EC,
  // Case 12
  ST_LIN,  COLOR0, EB, ED,
  // Case 13
  ST_LIN,  COLOR0, EA, EB,
  // Case 14
  ST_LIN,  COLOR0, EA, ED
  // Case 15
};
// clang-format on

int numCutShapesQua[] = {0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0};

int startCutShapesQua[] = {0, 0, 4, 8, 12, 16, 24, 28, 32, 36, 40, 48, 52, 56, 60, 64};

const size_t cutShapesQuaSize = sizeof(cutShapesQua) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
