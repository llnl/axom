// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sina/core/File.hpp"
#include "axom/sina/core/ConduitUtil.hpp"

namespace sina = axom::sina;

char const EXPECTED_MIMETYPE_KEY[] = "mimetype";
char const EXPECTED_TAGS_KEY[] = "tags";

TEST(File, construct_differentType)
{
  sina::File f1 {"from literal"};
  sina::File f2 {std::string {"from std::string"}};
  EXPECT_EQ("from literal", f1.getUri());
  EXPECT_EQ("from std::string", f2.getUri());
}

TEST(File, setMimeType)
{
  sina::File file {"the URI"};
  file.setMimeType("mime");
  EXPECT_EQ("the URI", file.getUri());
  EXPECT_EQ("mime", file.getMimeType());
}

TEST(File, setTags)
{
  std::vector<std::string> tags = {"these", "are", "tags"};
  sina::File file {"the URI"};
  file.setTags(tags);
  EXPECT_EQ("the URI", file.getUri());
  EXPECT_EQ(tags, file.getTags());
}

TEST(File, create_fromNode_basic)
{
  std::string uri = "the URI";
  conduit::Node basic_file(conduit::DataType::object());
  sina::File file {uri, basic_file};
  EXPECT_EQ(uri, file.getUri());
  EXPECT_EQ("", file.getMimeType());
  EXPECT_EQ(0, file.getTags().size());
}

TEST(File, create_fromNode_complete)
{
  std::string uri = "another/uri.txt";
  std::vector<std::string> tags = {"tags", "are", "fun"};
  conduit::Node full_file(conduit::DataType::object());
  full_file[EXPECTED_MIMETYPE_KEY] = "the mime type";
  sina::addStringsToNode(full_file, EXPECTED_TAGS_KEY, tags);
  sina::File file {uri, full_file};
  EXPECT_EQ(uri, file.getUri());
  EXPECT_EQ("the mime type", file.getMimeType());
  EXPECT_EQ(tags, file.getTags());
}

TEST(File, toNode_basic)
{
  sina::File file {"the URI"};
  auto asNode = file.toNode();
  EXPECT_FALSE(asNode.has_child(EXPECTED_MIMETYPE_KEY));
  EXPECT_FALSE(asNode.has_child(EXPECTED_TAGS_KEY));
}

TEST(File, toNode_complete)
{
  std::vector<std::string> tags = {"these", "are", "tags"};
  sina::File file {"the URI"};
  file.setMimeType("the mime type");
  file.setTags(tags);
  auto asNode = file.toNode();
  EXPECT_EQ("the mime type", asNode[EXPECTED_MIMETYPE_KEY].as_string());
  EXPECT_EQ(tags, file.getTags());
}
