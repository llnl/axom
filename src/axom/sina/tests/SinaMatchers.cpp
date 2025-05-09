// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina/tests/SinaMatchers.hpp"
#include <string>

namespace axom
{
namespace sina
{
namespace testing
{

conduit::Node parseJsonValue(const std::string& valueAsString)
{
  // If we just try to do node.parse(valueAsString, "json"), then passing
  // in strings does not work. We need to create a document with a key
  // so that valueAsString can be parsed as a value.
  conduit::Node node;
  std::string fullContents = "{\"TEST_KEY\": ";
  fullContents += valueAsString;
  fullContents += "}";
  node.parse(fullContents, "json");
  return node.child("TEST_KEY");
}

// Define the MatchesJson class
MatchesJson::MatchesJson(const std::string& expectedJsonString)
  : expectedJsonString_(expectedJsonString)
{ }

bool MatchesJson::MatchAndExplain(const conduit::Node& node,
                                  ::testing::MatchResultListener* listener) const
{
  conduit::Node expected = parseJsonValue(expectedJsonString_);
  *listener << "Given node is " << node.to_json_default();
  conduit::Node diff;
  bool differ = expected.diff(node, diff);
  return !differ;
}

void MatchesJson::DescribeTo(std::ostream* os) const
{
  *os << "matches JSON: " << expectedJsonString_;
}

void MatchesJson::DescribeNegationTo(std::ostream* os) const
{
  *os << "does not match JSON: " << expectedJsonString_;
}

}  // namespace testing
}  // namespace sina
}  // namespace axom
