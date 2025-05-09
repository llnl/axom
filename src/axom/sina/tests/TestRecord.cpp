// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina/tests/TestRecord.hpp"

namespace axom
{
namespace sina
{
namespace testing
{

template <>
TestRecord<std::string>::TestRecord(conduit::Node const &asNode)
  : Record {asNode}
  , value {getRequiredString(TEST_RECORD_VALUE_KEY, asNode, "TestRecord")}
{ }

template <>
TestRecord<int>::TestRecord(conduit::Node const &asNode)
  : Record {asNode}
  , value {getRequiredField(TEST_RECORD_VALUE_KEY, asNode, "TestRecord").as_int()}
{ }

}  // namespace testing
}  // namespace sina
}  // namespace axom
