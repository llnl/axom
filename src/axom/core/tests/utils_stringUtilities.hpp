// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/StringUtilities.hpp"

#include <string>
#include <vector>

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(utils_stringUtilities, startsWith)
{
  {
    std::string testString = "foo.bar";
    std::string testPrefix = "foo";
    EXPECT_TRUE(axom::utilities::string::startsWith(testString, testPrefix));
  }

  {
    std::string testString = "foo.bar";
    EXPECT_TRUE(axom::utilities::string::startsWith(testString, 'f'));
  }

  {
    std::string testString = "foo.bar";
    std::string testPrefix = "boo";
    EXPECT_FALSE(axom::utilities::string::startsWith(testString, testPrefix));
  }

  {
    std::string testString = "foo.bar";
    EXPECT_FALSE(axom::utilities::string::startsWith(testString, 'b'));
  }
}

TEST(utils_stringUtilities, endsWith)
{
  {
    std::string testString = "foo.bar";
    std::string testSuffix = ".bar";
    EXPECT_TRUE(axom::utilities::string::endsWith(testString, testSuffix));
  }

  {
    std::string testString = "foo.bar";
    EXPECT_TRUE(axom::utilities::string::endsWith(testString, 'r'));
  }

  {
    std::string testString = "foo.bar";
    std::string testSuffix = ".baz";
    EXPECT_FALSE(axom::utilities::string::endsWith(testString, testSuffix));
  }

  {
    std::string testString = "foo.bar";
    EXPECT_FALSE(axom::utilities::string::endsWith(testString, 'z'));
  }
}

TEST(utils_stringUtilities, removeSuffix)
{
  // same suffix
  {
    std::string testString = "foo.bar";
    std::string testSuffix = ".bar";
    EXPECT_EQ("foo", axom::utilities::string::removeSuffix(testString, testSuffix));
  }

  // different suffix
  {
    std::string testString = "foo.bar";
    std::string testSuffix = ".baz";
    EXPECT_EQ(testString, axom::utilities::string::removeSuffix(testString, testSuffix));
  }

  // repeated suffix -- only removes one copy
  {
    std::string testString = "foo.bar.bar";
    std::string testSuffix = ".bar";
    EXPECT_EQ("foo.bar", axom::utilities::string::removeSuffix(testString, testSuffix));
  }

  // don't remove from middle
  {
    std::string testString = "foo.bar.baz";
    std::string testSuffix = ".bar";
    EXPECT_EQ(testString, axom::utilities::string::removeSuffix(testString, testSuffix));
  }
}

TEST(utils_stringUtilities, toLower)
{
  // already lower
  {
    std::string testString = "foo.bar";
    std::string expString = "foo.bar";
    axom::utilities::string::toLower(testString);
    EXPECT_EQ(expString, testString);
  }

  // all upper
  {
    std::string testString = "FOO.BAR";
    std::string expString = "foo.bar";
    axom::utilities::string::toLower(testString);
    EXPECT_EQ(expString, testString);
  }

  // mixed
  {
    std::string testString = "FoO.bAr";
    std::string expString = "foo.bar";
    axom::utilities::string::toLower(testString);
    EXPECT_EQ(expString, testString);
  }
}

TEST(utils_stringUtilities, toUpper)
{
  // already upper
  {
    std::string testString = "FOO.BAR";
    std::string expString = "FOO.BAR";
    axom::utilities::string::toUpper(testString);
    EXPECT_EQ(expString, testString);
  }

  // all lower
  {
    std::string testString = "foo.bar";
    std::string expString = "FOO.BAR";
    axom::utilities::string::toUpper(testString);
    EXPECT_EQ(expString, testString);
  }

  // mixed
  {
    std::string testString = "FoO.bAr";
    std::string expString = "FOO.BAR";
    axom::utilities::string::toUpper(testString);
    EXPECT_EQ(expString, testString);
  }
}

TEST(utils_stringUtilities, split)
{
  using StrVec = std::vector<std::string>;

  {
    std::string testString = "foo/bar/baz";
    StrVec exp {"foo", "bar", "baz"};
    StrVec results = axom::utilities::string::split(testString, '/');
    EXPECT_EQ(exp, results);
  }

  // Test empty string
  {
    std::string testString = "";
    StrVec exp;
    StrVec results = axom::utilities::string::split(testString, '/');
    EXPECT_EQ(exp, results);
  }

  // Test single token
  {
    std::string testString = "foo";
    StrVec exp {"foo"};
    StrVec results = axom::utilities::string::split(testString, '/');
    EXPECT_EQ(exp, results);
  }

  // Test other delimiter: '.'
  {
    std::string testString = "foo.bar.baz";
    StrVec exp {"foo", "bar", "baz"};
    StrVec results = axom::utilities::string::split(testString, '.');
    EXPECT_EQ(exp, results);
  }

  // Test different delimiter: '.' in string, ';' in function call
  {
    std::string testString = "foo.bar.baz";
    StrVec exp {"foo.bar.baz"};
    StrVec results = axom::utilities::string::split(testString, ';');
    EXPECT_EQ(exp, results);
  }
}

TEST(utils_stringUtilities, rsplitN)
{
  using StrVec = std::vector<std::string>;

  // Test w/ sufficient tokens
  {
    const std::size_t N = 3;
    std::string testString = "foo/bar/baz";
    StrVec exp {"foo", "bar", "baz"};
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }

  // Test w/ more string tokens than N
  {
    const std::size_t N = 2;
    std::string testString = "foo/bar/baz";
    StrVec exp {"foo/bar", "baz"};
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }

  // Test w/ more string tokens than N
  {
    const std::size_t N = 1;
    std::string testString = "foo/bar/baz";
    StrVec exp {"foo/bar/baz"};
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }

  // Test w/ fewer string tokens than N
  {
    const std::size_t N = 10;
    std::string testString = "foo/bar/baz";
    StrVec exp {"foo", "bar", "baz"};
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }

  // Test w/ N==0
  {
    const std::size_t N = 0;
    std::string testString = "foo/bar/baz";
    StrVec exp;
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_TRUE(results.empty());
    EXPECT_EQ(exp, results);
  }

  // Test w/ empty stringrsplitN
  {
    const std::size_t N = 3;
    std::string testString = "";
    StrVec exp;
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }

  // Test w/ single token in string
  {
    const std::size_t N = 3;
    std::string testString = "foo";
    StrVec exp {"foo"};
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }

  // Test mixed delim: '/'
  {
    const std::size_t N = 3;
    std::string testString = "foo.bar/baz.qux";
    StrVec exp {"foo.bar", "baz.qux"};
    StrVec results = axom::utilities::string::rsplitN(testString, N, '/');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }

  // Test mixed delim: '.'
  {
    const std::size_t N = 3;
    std::string testString = "foo.bar/baz.qux";
    StrVec exp {"foo", "bar/baz", "qux"};
    StrVec results = axom::utilities::string::rsplitN(testString, N, '.');
    EXPECT_LE(results.size(), N);
    EXPECT_EQ(exp, results);
  }
}

TEST(utils_stringUtilities, replaceAllInstances)
{
  // Test single replace
  {
    std::string testString = "baz.bar";
    std::string expString = "foo.bar";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "baz", "foo");
    EXPECT_EQ(expString, results);
  }

  // Test multiple replace
  {
    std::string testString = "baz.baz.baz";
    std::string expString = "foo.foo.foo";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "baz", "foo");
    EXPECT_EQ(expString, results);
  }

  // Test no match
  {
    std::string testString = "foo.bar";
    std::string expString = "foo.bar";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "baz", "foo");
    EXPECT_EQ(expString, results);
  }

  // Test empty find string
  {
    std::string testString = "foo.bar";
    std::string expString = "foo.bar";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "", "foo");
    EXPECT_EQ(expString, results);
  }

  // Test empty replace string
  {
    std::string testString = "foo.bar";
    std::string expString = "foo.";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "bar", "");
    EXPECT_EQ(expString, results);
  }

  // Test empty find and replace string
  {
    std::string testString = "foo.bar";
    std::string expString = "foo.bar";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "", "");
    EXPECT_EQ(expString, results);
  }

  // // Test same find and replace string
  {
    std::string testString = "foo.bar";
    std::string expString = "foo.bar";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "foo", "foo");
    EXPECT_EQ(expString, results);
  }

  // Test find string is substring of replace string
  {
    std::string testString = "foo.bar";
    std::string expString = "fooStr.bar";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "foo", "fooStr");
    EXPECT_EQ(expString, results);
  }

  // Test multiple find string is substring of replace string
  {
    std::string testString = "foo.foo";
    std::string expString = "foofoo.foofoo";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "foo", "foofoo");
    EXPECT_EQ(expString, results);
  }

  // Test replace literal newline with newline
  {
    std::string testString = "foo\\n.bar\\n";
    std::string expString = "foo\n.bar\n";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "\\n", "\n");
    EXPECT_EQ(expString, results);
  }

  // Test replace newline with literal newline
  {
    std::string testString = "foo\n.bar\n";
    std::string expString = "foo\\n.bar\\n";
    std::string results = axom::utilities::string::replaceAllInstances(testString, "\n", "\\n");
    EXPECT_EQ(expString, results);
  }
}
