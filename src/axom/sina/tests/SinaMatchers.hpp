// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SINAMATCHERS_HPP
#define AXOM_SINAMATCHERS_HPP

#include <gmock/gmock.h>
#include "conduit.hpp"

// Note: LLNL's blueos complained about using GMock's MATCHER_P macros, as originally implemented.
// This required explicitly generating equivalent classes for the Matcher functionality.

namespace axom
{
namespace sina
{
namespace testing
{

// Function to parse JSON value
conduit::Node parseJsonValue(const std::string& valueAsString);

// Matcher class
class MatchesJson
{
public:
  explicit MatchesJson(const std::string& expectedJsonString);

  bool MatchAndExplain(const conduit::Node& node, ::testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* os) const;

  void DescribeNegationTo(std::ostream* os) const;

private:
  const std::string expectedJsonString_;
};

// Helper function to create the matcher
inline ::testing::PolymorphicMatcher<MatchesJson> MatchesJsonMatcher(const std::string& expectedJsonString)
{
  return ::testing::MakePolymorphicMatcher(MatchesJson(expectedJsonString));
}

}  // namespace testing
}  // namespace sina
}  // namespace axom

#endif  //AXOM_SINAMATCHERS_HPP