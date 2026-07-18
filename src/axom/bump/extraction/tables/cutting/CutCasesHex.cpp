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

int numCutCasesHex = 256;

// clang-format off
unsigned char cutShapesHex[] = {
  // Case 0
  // Case 1
  ST_TRI,   COLOR0, EA, ED, EI,
  // Case 2
  ST_TRI,   COLOR0, EA, EJ, EB,
  // Case 3
  ST_QUA,  COLOR0, EB, ED, EI, EJ,
  // Case 4
  ST_TRI,   COLOR0, EB, EL, EC,
  // Case 5
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_TRI,   COLOR0, EB, EL, EC,
  // Case 6
  ST_QUA,  COLOR0, EJ, EL, EC, EA,
  // Case 7
  ST_POLY5, COLOR0, EI, EJ, EL, EC, ED,
  // Case 8
  ST_TRI,   COLOR0, ED, EC, EK,
  // Case 9
  ST_QUA,  COLOR0, EA, EC, EK, EI,
  // Case 10
  ST_TRI,   COLOR0, EB, EA, EJ,
  ST_TRI,   COLOR0, EC, EK, ED,
  // Case 11
  ST_POLY5, COLOR0, EK, EI, EJ, EB, EC,
  // Case 12
  ST_QUA,  COLOR0, ED, EB, EL, EK,
  // Case 13
  ST_POLY5, COLOR0, EL, EK, EI, EA, EB,
  // Case 14
  ST_POLY5, COLOR0, EJ, EL, EK, ED, EA,
  // Case 15
  ST_QUA,  COLOR0, EI, EJ, EL, EK,
  // Case 16
  ST_TRI,   COLOR0, EE, EI, EH,
  // Case 17
  ST_QUA,  COLOR0, EE, EA, ED, EH,
  // Case 18
  ST_TRI,   COLOR0, EA, EJ, EB,
  ST_TRI,   COLOR0, EI, EH, EE,
  // Case 19
  ST_POLY5, COLOR0, EB, ED, EH, EE, EJ,
  // Case 20
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_TRI,   COLOR0, EI, EH, EE,
  // Case 21
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_QUA,  COLOR0, ED, EH, EE, EA,
  // Case 22
  ST_TRI,   COLOR0, EI, EH, EE,
  ST_QUA,  COLOR0, EJ, EL, EC, EA,
  // Case 23
  ST_POLY6, COLOR0, EH, EE, EJ, EL, EC, ED,
  // Case 24
  ST_TRI,   COLOR0, EI, EH, EE,
  ST_TRI,   COLOR0, ED, EC, EK,
  // Case 25
  ST_POLY5, COLOR0, EE, EA, EC, EK, EH,
  // Case 26
  ST_TRI,   COLOR0, EJ, EB, EA,
  ST_TRI,   COLOR0, EI, EH, EE,
  ST_TRI,   COLOR0, EC, EK, ED,
  // Case 27
  ST_POLY6, COLOR0, EK, EH, EE, EJ, EB, EC,
  // Case 28
  ST_TRI,   COLOR0, EH, EE, EI,
  ST_QUA,  COLOR0, ED, EB, EL, EK,
  // Case 29
  ST_POLY6, COLOR0, EK, EH, EE, EA, EB, EL,
  // Case 30
  ST_TRI,   COLOR0, EE, EI, EH,
  ST_POLY5, COLOR0, EK, ED, EA, EJ, EL,
  // Case 31
  ST_POLY5, COLOR0, EJ, EL, EK, EH, EE,
  // Case 32
  ST_TRI,   COLOR0, EJ, EE, EF,
  // Case 33
  ST_TRI,   COLOR0, EJ, EE, EF,
  ST_TRI,   COLOR0, EA, ED, EI,
  // Case 34
  ST_QUA,  COLOR0, EA, EE, EF, EB,
  // Case 35
  ST_POLY5, COLOR0, EF, EB, ED, EI, EE,
  // Case 36
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_TRI,   COLOR0, EJ, EE, EF,
  // Case 37
  ST_TRI,   COLOR0, ED, EI, EA,
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_TRI,   COLOR0, EE, EF, EJ,
  // Case 38
  ST_POLY5, COLOR0, EC, EA, EE, EF, EL,
  // Case 39
  ST_POLY6, COLOR0, EF, EL, EC, ED, EI, EE,
  // Case 40
  ST_TRI,   COLOR0, EJ, EE, EF,
  ST_TRI,   COLOR0, EC, EK, ED,
  // Case 41
  ST_TRI,   COLOR0, EE, EF, EJ,
  ST_QUA,  COLOR0, EA, EC, EK, EI,
  // Case 42
  ST_TRI,   COLOR0, EC, EK, ED,
  ST_QUA,  COLOR0, EA, EE, EF, EB,
  // Case 43
  ST_POLY6, COLOR0, EI, EE, EF, EB, EC, EK,
  // Case 44
  ST_TRI,   COLOR0, EJ, EE, EF,
  ST_QUA,  COLOR0, EL, EK, ED, EB,
  // Case 45
  ST_TRI,   COLOR0, EE, EF, EJ,
  ST_POLY5, COLOR0, EL, EK, EI, EA, EB,
  // Case 46
  ST_POLY6, COLOR0, EK, ED, EA, EE, EF, EL,
  // Case 47
  ST_POLY5, COLOR0, EL, EK, EI, EE, EF,
  // Case 48
  ST_QUA,  COLOR0, EJ, EI, EH, EF,
  // Case 49
  ST_POLY5, COLOR0, ED, EH, EF, EJ, EA,
  // Case 50
  ST_POLY5, COLOR0, EH, EF, EB, EA, EI,
  // Case 51
  ST_QUA,  COLOR0, EF, EB, ED, EH,
  // Case 52
  ST_TRI,   COLOR0, EL, EC, EB,
  ST_QUA,  COLOR0, EJ, EI, EH, EF,
  // Case 53
  ST_TRI,   COLOR0, EL, EC, EB,
  ST_POLY5, COLOR0, ED, EH, EF, EJ, EA,
  // Case 54
  ST_POLY6, COLOR0, EF, EL, EC, EA, EI, EH,
  // Case 55
  ST_POLY5, COLOR0, ED, EH, EF, EL, EC,
  // Case 56
  ST_TRI,   COLOR0, ED, EC, EK,
  ST_QUA,  COLOR0, EH, EF, EJ, EI,
  // Case 57
  ST_POLY6, COLOR0, EC, EK, EH, EF, EJ, EA,
  // Case 58
  ST_TRI,   COLOR0, EC, EK, ED,
  ST_POLY5, COLOR0, EH, EF, EB, EA, EI,
  // Case 59
  ST_POLY5, COLOR0, EH, EF, EB, EC, EK,
  // Case 60
  ST_QUA,  COLOR0, EF, EJ, EI, EH,
  ST_QUA,  COLOR0, ED, EB, EL, EK,
  // Case 61
  ST_POLY7, COLOR0, EH, EF, EJ, EA, EB, EL, EK,
  // Case 62
  ST_POLY7, COLOR0, EL, EK, ED, EA, EI, EH, EF,
  // Case 63
  ST_QUA,  COLOR0, EF, EL, EK, EH,
  // Case 64
  ST_TRI,   COLOR0, EL, EF, EG,
  // Case 65
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_TRI,   COLOR0, EF, EG, EL,
  // Case 66
  ST_TRI,   COLOR0, EJ, EB, EA,
  ST_TRI,   COLOR0, EF, EG, EL,
  // Case 67
  ST_TRI,   COLOR0, EF, EG, EL,
  ST_QUA,  COLOR0, EB, ED, EI, EJ,
  // Case 68
  ST_QUA,  COLOR0, EB, EF, EG, EC,
  // Case 69
  ST_TRI,   COLOR0, ED, EI, EA,
  ST_QUA,  COLOR0, EB, EF, EG, EC,
  // Case 70
  ST_POLY5, COLOR0, EG, EC, EA, EJ, EF,
  // Case 71
  ST_POLY6, COLOR0, EC, ED, EI, EJ, EF, EG,
  // Case 72
  ST_TRI,   COLOR0, EC, EK, ED,
  ST_TRI,   COLOR0, EL, EF, EG,
  // Case 73
  ST_TRI,   COLOR0, EL, EF, EG,
  ST_QUA,  COLOR0, EK, EI, EA, EC,
  // Case 74
  ST_TRI,   COLOR0, EA, EJ, EB,
  ST_TRI,   COLOR0, EC, EK, ED,
  ST_TRI,   COLOR0, EF, EG, EL,
  // Case 75
  ST_TRI,   COLOR0, EF, EG, EL,
  ST_POLY5, COLOR0, EK, EI, EJ, EB, EC,
  // Case 76
  ST_POLY5, COLOR0, ED, EB, EF, EG, EK,
  // Case 77
  ST_POLY6, COLOR0, EF, EG, EK, EI, EA, EB,
  // Case 78
  ST_POLY6, COLOR0, EG, EK, ED, EA, EJ, EF,
  // Case 79
  ST_POLY5, COLOR0, EK, EI, EJ, EF, EG,
  // Case 80
  ST_TRI,   COLOR0, EF, EG, EL,
  ST_TRI,   COLOR0, EE, EI, EH,
  // Case 81
  ST_TRI,   COLOR0, EG, EL, EF,
  ST_QUA,  COLOR0, EE, EA, ED, EH,
  // Case 82
  ST_TRI,   COLOR0, EB, EA, EJ,
  ST_TRI,   COLOR0, EF, EG, EL,
  ST_TRI,   COLOR0, EI, EH, EE,
  // Case 83
  ST_TRI,   COLOR0, EL, EF, EG,
  ST_POLY5, COLOR0, EH, EE, EJ, EB, ED,
  // Case 84
  ST_TRI,   COLOR0, EE, EI, EH,
  ST_QUA,  COLOR0, EG, EC, EB, EF,
  // Case 85
  ST_QUA,  COLOR0, EC, EB, EF, EG,
  ST_QUA,  COLOR0, EE, EA, ED, EH,
  // Case 86
  ST_TRI,   COLOR0, EI, EH, EE,
  ST_POLY5, COLOR0, EG, EC, EA, EJ, EF,
  // Case 87
  ST_POLY7, COLOR0, ED, EH, EE, EJ, EF, EG, EC,
  // Case 88
  ST_TRI,   COLOR0, ED, EC, EK,
  ST_TRI,   COLOR0, EH, EE, EI,
  ST_TRI,   COLOR0, EL, EF, EG,
  // Case 89
  ST_TRI,   COLOR0, EF, EG, EL,
  ST_POLY5, COLOR0, EC, EK, EH, EE, EA,
  // Case 90
  ST_TRI,   COLOR0, EA, EJ, EB,
  ST_TRI,   COLOR0, EE, EI, EH,
  ST_TRI,   COLOR0, EC, EK, ED,
  ST_TRI,   COLOR0, EF, EG, EL,
  // Case 91
  ST_TRI,   COLOR0, EF, EG, EL,
  ST_POLY6, COLOR0, EJ, EB, EC, EK, EH, EE,
  // Case 92
  ST_TRI,   COLOR0, EI, EH, EE,
  ST_POLY5, COLOR0, EF, EG, EK, ED, EB,
  // Case 93
  ST_POLY7, COLOR0, EB, EF, EG, EK, EH, EE, EA,
  // Case 94
  ST_TRI,   COLOR0, EI, EH, EE,
  ST_POLY6, COLOR0, EA, EJ, EF, EG, EK, ED,
  // Case 95
  ST_POLY6, COLOR0, EJ, EF, EG, EK, EH, EE,
  // Case 96
  ST_QUA,  COLOR0, EL, EJ, EE, EG,
  // Case 97
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_QUA,  COLOR0, EE, EG, EL, EJ,
  // Case 98
  ST_POLY5, COLOR0, EA, EE, EG, EL, EB,
  // Case 99
  ST_POLY6, COLOR0, EG, EL, EB, ED, EI, EE,
  // Case 100
  ST_POLY5, COLOR0, EE, EG, EC, EB, EJ,
  // Case 101
  ST_TRI,   COLOR0, ED, EI, EA,
  ST_POLY5, COLOR0, EE, EG, EC, EB, EJ,
  // Case 102
  ST_QUA,  COLOR0, EC, EA, EE, EG,
  // Case 103
  ST_POLY5, COLOR0, EE, EG, EC, ED, EI,
  // Case 104
  ST_TRI,   COLOR0, EK, ED, EC,
  ST_QUA,  COLOR0, EL, EJ, EE, EG,
  // Case 105
  ST_QUA,  COLOR0, EI, EA, EC, EK,
  ST_QUA,  COLOR0, EL, EJ, EE, EG,
  // Case 106
  ST_TRI,   COLOR0, ED, EC, EK,
  ST_POLY5, COLOR0, EG, EL, EB, EA, EE,
  // Case 107
  ST_POLY7, COLOR0, EE, EG, EL, EB, EC, EK, EI,
  // Case 108
  ST_POLY6, COLOR0, EG, EK, ED, EB, EJ, EE,
  // Case 109
  ST_POLY7, COLOR0, EK, EI, EA, EB, EJ, EE, EG,
  // Case 110
  ST_POLY5, COLOR0, EA, EE, EG, EK, ED,
  // Case 111
  ST_QUA,  COLOR0, EI, EE, EG, EK,
  // Case 112
  ST_POLY5, COLOR0, EL, EJ, EI, EH, EG,
  // Case 113
  ST_POLY6, COLOR0, EH, EG, EL, EJ, EA, ED,
  // Case 114
  ST_POLY6, COLOR0, EH, EG, EL, EB, EA, EI,
  // Case 115
  ST_POLY5, COLOR0, EB, ED, EH, EG, EL,
  // Case 116
  ST_POLY6, COLOR0, EI, EH, EG, EC, EB, EJ,
  // Case 117
  ST_POLY7, COLOR0, EG, EC, EB, EJ, EA, ED, EH,
  // Case 118
  ST_POLY5, COLOR0, EG, EC, EA, EI, EH,
  // Case 119
  ST_QUA,  COLOR0, EC, ED, EH, EG,
  // Case 120
  ST_TRI,   COLOR0, EC, EK, ED,
  ST_POLY5, COLOR0, EI, EH, EG, EL, EJ,
  // Case 121
  ST_POLY7, COLOR0, EA, EC, EK, EH, EG, EL, EJ,
  // Case 122
  ST_TRI,   COLOR0, EC, EK, ED,
  ST_POLY6, COLOR0, EB, EA, EI, EH, EG, EL,
  // Case 123
  ST_POLY6, COLOR0, EB, EC, EK, EH, EG, EL,
  // Case 124
  ST_POLY7, COLOR0, EJ, EI, EH, EG, EK, ED, EB,
  // Case 125
  ST_TRI,   COLOR0, EA, EB, EJ,
  ST_TRI,   COLOR0, EK, EH, EG,
  // Case 126
  ST_POLY6, COLOR0, EA, EI, EH, EG, EK, ED,
  // Case 127
  ST_TRI,   COLOR0, EH, EG, EK,
  // Case 128
  ST_TRI,   COLOR0, EH, EK, EG,
  // Case 129
  ST_TRI,   COLOR0, ED, EI, EA,
  ST_TRI,   COLOR0, EK, EG, EH,
  // Case 130
  ST_TRI,   COLOR0, EA, EJ, EB,
  ST_TRI,   COLOR0, EK, EG, EH,
  // Case 131
  ST_TRI,   COLOR0, EK, EG, EH,
  ST_QUA,  COLOR0, EI, EJ, EB, ED,
  // Case 132
  ST_TRI,   COLOR0, EL, EC, EB,
  ST_TRI,   COLOR0, EG, EH, EK,
  // Case 133
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_TRI,   COLOR0, ED, EI, EA,
  ST_TRI,   COLOR0, EG, EH, EK,
  // Case 134
  ST_TRI,   COLOR0, EG, EH, EK,
  ST_QUA,  COLOR0, EC, EA, EJ, EL,
  // Case 135
  ST_TRI,   COLOR0, EG, EH, EK,
  ST_POLY5, COLOR0, EI, EJ, EL, EC, ED,
  // Case 136
  ST_QUA,  COLOR0, EH, ED, EC, EG,
  // Case 137
  ST_POLY5, COLOR0, EA, EC, EG, EH, EI,
  // Case 138
  ST_TRI,   COLOR0, EA, EJ, EB,
  ST_QUA,  COLOR0, EC, EG, EH, ED,
  // Case 139
  ST_POLY6, COLOR0, EG, EH, EI, EJ, EB, EC,
  // Case 140
  ST_POLY5, COLOR0, EH, ED, EB, EL, EG,
  // Case 141
  ST_POLY6, COLOR0, EB, EL, EG, EH, EI, EA,
  // Case 142
  ST_POLY6, COLOR0, EL, EG, EH, ED, EA, EJ,
  // Case 143
  ST_POLY5, COLOR0, EI, EJ, EL, EG, EH,
  // Case 144
  ST_QUA,  COLOR0, EG, EE, EI, EK,
  // Case 145
  ST_POLY5, COLOR0, EG, EE, EA, ED, EK,
  // Case 146
  ST_TRI,   COLOR0, EJ, EB, EA,
  ST_QUA,  COLOR0, EI, EK, EG, EE,
  // Case 147
  ST_POLY6, COLOR0, ED, EK, EG, EE, EJ, EB,
  // Case 148
  ST_TRI,   COLOR0, EC, EB, EL,
  ST_QUA,  COLOR0, EG, EE, EI, EK,
  // Case 149
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_POLY5, COLOR0, EG, EE, EA, ED, EK,
  // Case 150
  ST_QUA,  COLOR0, EE, EI, EK, EG,
  ST_QUA,  COLOR0, EC, EA, EJ, EL,
  // Case 151
  ST_POLY7, COLOR0, EJ, EL, EC, ED, EK, EG, EE,
  // Case 152
  ST_POLY5, COLOR0, EC, EG, EE, EI, ED,
  // Case 153
  ST_QUA,  COLOR0, EE, EA, EC, EG,
  // Case 154
  ST_TRI,   COLOR0, EB, EA, EJ,
  ST_POLY5, COLOR0, EE, EI, ED, EC, EG,
  // Case 155
  ST_POLY5, COLOR0, EC, EG, EE, EJ, EB,
  // Case 156
  ST_POLY6, COLOR0, EB, EL, EG, EE, EI, ED,
  // Case 157
  ST_POLY5, COLOR0, EG, EE, EA, EB, EL,
  // Case 158
  ST_POLY7, COLOR0, EG, EE, EI, ED, EA, EJ, EL,
  // Case 159
  ST_QUA,  COLOR0, EE, EJ, EL, EG,
  // Case 160
  ST_TRI,   COLOR0, EE, EF, EJ,
  ST_TRI,   COLOR0, EH, EK, EG,
  // Case 161
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_TRI,   COLOR0, EE, EF, EJ,
  ST_TRI,   COLOR0, EK, EG, EH,
  // Case 162
  ST_TRI,   COLOR0, EH, EK, EG,
  ST_QUA,  COLOR0, EF, EB, EA, EE,
  // Case 163
  ST_TRI,   COLOR0, EK, EG, EH,
  ST_POLY5, COLOR0, EF, EB, ED, EI, EE,
  // Case 164
  ST_TRI,   COLOR0, EJ, EE, EF,
  ST_TRI,   COLOR0, EL, EC, EB,
  ST_TRI,   COLOR0, EH, EK, EG,
  // Case 165
  ST_TRI,   COLOR0, EG, EH, EK,
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_TRI,   COLOR0, EE, EF, EJ,
  // Case 166
  ST_TRI,   COLOR0, EH, EK, EG,
  ST_POLY5, COLOR0, EC, EA, EE, EF, EL,
  // Case 167
  ST_TRI,   COLOR0, EK, EG, EH,
  ST_POLY6, COLOR0, ED, EI, EE, EF, EL, EC,
  // Case 168
  ST_TRI,   COLOR0, EF, EJ, EE,
  ST_QUA,  COLOR0, EH, ED, EC, EG,
  // Case 169
  ST_TRI,   COLOR0, EJ, EE, EF,
  ST_POLY5, COLOR0, EG, EH, EI, EA, EC,
  // Case 170
  ST_QUA,  COLOR0, ED, EC, EG, EH,
  ST_QUA,  COLOR0, EF, EB, EA, EE,
  // Case 171
  ST_POLY7, COLOR0, EC, EG, EH, EI, EE, EF, EB,
  // Case 172
  ST_TRI,   COLOR0, EJ, EE, EF,
  ST_POLY5, COLOR0, EH, ED, EB, EL, EG,
  // Case 173
  ST_TRI,   COLOR0, EJ, EE, EF,
  ST_POLY6, COLOR0, EB, EL, EG, EH, EI, EA,
  // Case 174
  ST_POLY7, COLOR0, EA, EE, EF, EL, EG, EH, ED,
  // Case 175
  ST_POLY6, COLOR0, EL, EG, EH, EI, EE, EF,
  // Case 176
  ST_POLY5, COLOR0, EJ, EI, EK, EG, EF,
  // Case 177
  ST_POLY6, COLOR0, EA, ED, EK, EG, EF, EJ,
  // Case 178
  ST_POLY6, COLOR0, EK, EG, EF, EB, EA, EI,
  // Case 179
  ST_POLY5, COLOR0, EF, EB, ED, EK, EG,
  // Case 180
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_POLY5, COLOR0, EK, EG, EF, EJ, EI,
  // Case 181
  ST_TRI,   COLOR0, EB, EL, EC,
  ST_POLY6, COLOR0, EA, ED, EK, EG, EF, EJ,
  // Case 182
  ST_POLY7, COLOR0, EI, EK, EG, EF, EL, EC, EA,
  // Case 183
  ST_POLY6, COLOR0, ED, EK, EG, EF, EL, EC,
  // Case 184
  ST_POLY6, COLOR0, EI, ED, EC, EG, EF, EJ,
  // Case 185
  ST_POLY5, COLOR0, EA, EC, EG, EF, EJ,
  // Case 186
  ST_POLY7, COLOR0, EF, EB, EA, EI, ED, EC, EG,
  // Case 187
  ST_QUA,  COLOR0, EG, EF, EB, EC,
  // Case 188
  ST_POLY7, COLOR0, ED, EB, EL, EG, EF, EJ, EI,
  // Case 189
  ST_POLY6, COLOR0, EA, EB, EL, EG, EF, EJ,
  // Case 190
  ST_TRI,   COLOR0, EA, EI, ED,
  ST_TRI,   COLOR0, EF, EL, EG,
  // Case 191
  ST_TRI,   COLOR0, EL, EG, EF,
  // Case 192
  ST_QUA,  COLOR0, EK, EL, EF, EH,
  // Case 193
  ST_TRI,   COLOR0, EI, EA, ED,
  ST_QUA,  COLOR0, EK, EL, EF, EH,
  // Case 194
  ST_TRI,   COLOR0, EB, EA, EJ,
  ST_QUA,  COLOR0, EF, EH, EK, EL,
  // Case 195
  ST_QUA,  COLOR0, EL, EF, EH, EK,
  ST_QUA,  COLOR0, EI, EJ, EB, ED,
  // Case 196
  ST_POLY5, COLOR0, EB, EF, EH, EK, EC,
  // Case 197
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_POLY5, COLOR0, EH, EK, EC, EB, EF,
  // Case 198
  ST_POLY6, COLOR0, EH, EK, EC, EA, EJ, EF,
  // Case 199
  ST_POLY7, COLOR0, EF, EH, EK, EC, ED, EI, EJ,
  // Case 200
  ST_POLY5, COLOR0, EF, EH, ED, EC, EL,
  // Case 201
  ST_POLY6, COLOR0, EC, EL, EF, EH, EI, EA,
  // Case 202
  ST_TRI,   COLOR0, EJ, EB, EA,
  ST_POLY5, COLOR0, ED, EC, EL, EF, EH,
  // Case 203
  ST_POLY7, COLOR0, EI, EJ, EB, EC, EL, EF, EH,
  // Case 204
  ST_QUA,  COLOR0, ED, EB, EF, EH,
  // Case 205
  ST_POLY5, COLOR0, EB, EF, EH, EI, EA,
  // Case 206
  ST_POLY5, COLOR0, EF, EH, ED, EA, EJ,
  // Case 207
  ST_QUA,  COLOR0, EH, EI, EJ, EF,
  // Case 208
  ST_POLY5, COLOR0, EI, EK, EL, EF, EE,
  // Case 209
  ST_POLY6, COLOR0, EA, ED, EK, EL, EF, EE,
  // Case 210
  ST_TRI,   COLOR0, EA, EJ, EB,
  ST_POLY5, COLOR0, EL, EF, EE, EI, EK,
  // Case 211
  ST_POLY7, COLOR0, EK, EL, EF, EE, EJ, EB, ED,
  // Case 212
  ST_POLY6, COLOR0, EF, EE, EI, EK, EC, EB,
  // Case 213
  ST_POLY7, COLOR0, EE, EA, ED, EK, EC, EB, EF,
  // Case 214
  ST_POLY7, COLOR0, EC, EA, EJ, EF, EE, EI, EK,
  // Case 215
  ST_TRI,   COLOR0, EJ, EF, EE,
  ST_TRI,   COLOR0, EC, ED, EK,
  // Case 216
  ST_POLY6, COLOR0, ED, EC, EL, EF, EE, EI,
  // Case 217
  ST_POLY5, COLOR0, EE, EA, EC, EL, EF,
  // Case 218
  ST_TRI,   COLOR0, EA, EJ, EB,
  ST_POLY6, COLOR0, ED, EC, EL, EF, EE, EI,
  // Case 219
  ST_POLY6, COLOR0, EC, EL, EF, EE, EJ, EB,
  // Case 220
  ST_POLY5, COLOR0, ED, EB, EF, EE, EI,
  // Case 221
  ST_QUA,  COLOR0, EF, EE, EA, EB,
  // Case 222
  ST_POLY6, COLOR0, EF, EE, EI, ED, EA, EJ,
  // Case 223
  ST_TRI,   COLOR0, EJ, EF, EE,
  // Case 224
  ST_POLY5, COLOR0, EK, EL, EJ, EE, EH,
  // Case 225
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_POLY5, COLOR0, EK, EL, EJ, EE, EH,
  // Case 226
  ST_POLY6, COLOR0, EE, EH, EK, EL, EB, EA,
  // Case 227
  ST_POLY7, COLOR0, EB, ED, EI, EE, EH, EK, EL,
  // Case 228
  ST_POLY6, COLOR0, EJ, EE, EH, EK, EC, EB,
  // Case 229
  ST_TRI,   COLOR0, EA, ED, EI,
  ST_POLY6, COLOR0, EJ, EE, EH, EK, EC, EB,
  // Case 230
  ST_POLY5, COLOR0, EC, EA, EE, EH, EK,
  // Case 231
  ST_POLY6, COLOR0, EE, EH, EK, EC, ED, EI,
  // Case 232
  ST_POLY6, COLOR0, EJ, EE, EH, ED, EC, EL,
  // Case 233
  ST_POLY7, COLOR0, EL, EJ, EE, EH, EI, EA, EC,
  // Case 234
  ST_POLY7, COLOR0, EH, ED, EC, EL, EB, EA, EE,
  // Case 235
  ST_TRI,   COLOR0, EB, EC, EL,
  ST_TRI,   COLOR0, EI, EE, EH,
  // Case 236
  ST_POLY5, COLOR0, EH, ED, EB, EJ, EE,
  // Case 237
  ST_POLY6, COLOR0, EB, EJ, EE, EH, EI, EA,
  // Case 238
  ST_QUA,  COLOR0, ED, EA, EE, EH,
  // Case 239
  ST_TRI,   COLOR0, EE, EH, EI,
  // Case 240
  ST_QUA,  COLOR0, EL, EJ, EI, EK,
  // Case 241
  ST_POLY5, COLOR0, EK, EL, EJ, EA, ED,
  // Case 242
  ST_POLY5, COLOR0, EI, EK, EL, EB, EA,
  // Case 243
  ST_QUA,  COLOR0, EL, EB, ED, EK,
  // Case 244
  ST_POLY5, COLOR0, EJ, EI, EK, EC, EB,
  // Case 245
  ST_POLY6, COLOR0, EJ, EA, ED, EK, EC, EB,
  // Case 246
  ST_QUA,  COLOR0, EK, EC, EA, EI,
  // Case 247
  ST_TRI,   COLOR0, ED, EK, EC,
  // Case 248
  ST_POLY5, COLOR0, EL, EJ, EI, ED, EC,
  // Case 249
  ST_QUA,  COLOR0, EC, EL, EJ, EA,
  // Case 250
  ST_POLY6, COLOR0, EI, ED, EC, EL, EB, EA,
  // Case 251
  ST_TRI,   COLOR0, EB, EC, EL,
  // Case 252
  ST_QUA,  COLOR0, EI, ED, EB, EJ,
  // Case 253
  ST_TRI,   COLOR0, EA, EB, EJ,
  // Case 254
  ST_TRI,   COLOR0, EA, EI, ED
  // Case 255
};
// clang-format on

int numCutShapesHex[] = {
  0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 1, 3, 1, 2, 1, 2, 1,
  1, 2, 1, 1, 2, 3, 1, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1,
  1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 3, 2, 1, 1, 1, 1, 2, 2, 3, 2, 2, 2, 2, 1, 3, 2, 4, 2, 2, 1, 2, 1,
  1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1,
  1, 2, 2, 2, 2, 3, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1,
  2, 3, 2, 2, 3, 4, 2, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
  1, 2, 2, 2, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1,
  1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};

int startCutShapesHex[] = {
  0,    0,    5,    10,   16,   21,   31,   37,   44,   49,   55,   65,   72,   78,   85,   92,
  98,   103,  109,  119,  126,  136,  147,  158,  166,  176,  183,  198,  206,  217,  225,  237,
  244,  249,  259,  265,  272,  282,  297,  304,  312,  322,  333,  344,  352,  363,  375,  383,
  390,  396,  403,  410,  416,  427,  439,  447,  454,  465,  473,  485,  492,  504,  513,  522,
  528,  533,  543,  553,  564,  570,  581,  588,  596,  606,  617,  632,  644,  651,  659,  667,
  674,  684,  695,  710,  722,  733,  745,  757,  766,  781,  793,  813,  826,  838,  847,  860,
  868,  874,  885,  892,  900,  907,  919,  925,  932,  943,  955,  967,  976,  984,  993,  1000,
  1006, 1013, 1021, 1029, 1036, 1044, 1053, 1060, 1066, 1078, 1087, 1100, 1108, 1117, 1127, 1135,
  1140, 1145, 1155, 1165, 1176, 1186, 1201, 1212, 1224, 1230, 1237, 1248, 1256, 1263, 1271, 1279,
  1286, 1292, 1299, 1310, 1318, 1329, 1341, 1353, 1362, 1369, 1375, 1387, 1394, 1402, 1409, 1418,
  1424, 1434, 1449, 1460, 1472, 1487, 1507, 1519, 1532, 1543, 1555, 1567, 1576, 1588, 1601, 1610,
  1618, 1625, 1633, 1641, 1648, 1660, 1673, 1682, 1690, 1698, 1705, 1714, 1720, 1729, 1737, 1747,
  1752, 1758, 1769, 1780, 1792, 1799, 1811, 1819, 1828, 1835, 1843, 1855, 1864, 1870, 1877, 1884,
  1890, 1897, 1905, 1917, 1926, 1934, 1943, 1952, 1962, 1970, 1977, 1990, 1998, 2005, 2011, 2019,
  2024, 2031, 2043, 2051, 2060, 2068, 2081, 2088, 2096, 2104, 2113, 2122, 2132, 2139, 2147, 2153,
  2158, 2164, 2171, 2178, 2184, 2191, 2199, 2205, 2210, 2217, 2223, 2231, 2236, 2242, 2247, 2252};

const size_t cutShapesHexSize = sizeof(cutShapesHex) / sizeof(unsigned char);

}  // namespace cutting
}  // namespace tables
}  // namespace extraction
}  // namespace bump
}  // namespace axom
