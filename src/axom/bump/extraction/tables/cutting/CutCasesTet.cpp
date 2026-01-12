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

int numCutCasesTet = 16;

// clang-format off
unsigned char cutShapesTet[] = {
  // Case 0
  // Case 1
  ST_TRI,   COLOR0, EA, ED, EC,
  // Case 2
  ST_TRI,   COLOR0, EA, EB, EE,
  // Case 3
  ST_QUA,  COLOR0, EE, ED, EC, EB,
  // Case 4
  ST_TRI,   COLOR0, EB, EC, EF,
  // Case 5
  ST_QUA,  COLOR0, ED, EF, EB, EA,
  // Case 6
  ST_QUA,  COLOR0, EA, EC, EF, EE,
  // Case 7
  ST_TRI,   COLOR0, ED, EF, EE,
  // Case 8
  ST_TRI,   COLOR0, ED, EE, EF,
  // Case 9
  ST_QUA,  COLOR0, EA, EE, EF, EC,
  // Case 10
  ST_QUA,  COLOR0, EF, ED, EA, EB,
  // Case 11
  ST_TRI,   COLOR0, EF, EC, EB,
  // Case 12
  ST_QUA,  COLOR0, ED, EE, EB, EC,
  // Case 13
  ST_TRI,   COLOR0, EA, EE, EB,
  // Case 14
  ST_TRI,   COLOR0, EA, EC, ED
  // Case 15
};
// clang-format on

int numCutShapesTet[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};

int startCutShapesTet[] = {0, 0, 5, 10, 16, 21, 27, 33, 38, 43, 49, 55, 60, 66, 71, 76};

const size_t cutShapesTetSize = sizeof(cutShapesTet) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
