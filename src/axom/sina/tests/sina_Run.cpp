// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <stdexcept>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/Run.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::HasSubstr;

char const EXPECTED_TYPE_KEY[] = "type";
char const EXPECTED_LOCAL_ID_KEY[] = "local_id";
char const EXPECTED_GLOBAL_ID_KEY[] = "id";
char const EXPECTED_APPLICATION_KEY[] = "application";
char const EXPECTED_VERSION_KEY[] = "version";
char const EXPECTED_USER_KEY[] = "user";

// Throughout, we have to use "axom::sina::Run" instead of just "Run" due to
// a conflict with the Run() function in gtest

TEST(Run, create_fromnode_valid)
{
  conduit::Node originNode;
  originNode[EXPECTED_TYPE_KEY] = "run";
  originNode[EXPECTED_GLOBAL_ID_KEY] = "the id";
  originNode[EXPECTED_APPLICATION_KEY] = "the app";
  originNode[EXPECTED_VERSION_KEY] = "1.2.3";
  originNode[EXPECTED_USER_KEY] = "jdoe";
  axom::sina::Run run {originNode};
  EXPECT_EQ("run", run.getType());
  EXPECT_EQ("the id", run.getId().getId());
  EXPECT_EQ(IDType::Global, run.getId().getType());
  EXPECT_EQ("the app", run.getApplication());
  EXPECT_EQ("1.2.3", run.getVersion());
  EXPECT_EQ("jdoe", run.getUser());
}

TEST(Run, create_fromNode_missingApplication)
{
  conduit::Node originNode;
  originNode[EXPECTED_TYPE_KEY] = "run";
  originNode[EXPECTED_GLOBAL_ID_KEY] = "the id";
  originNode[EXPECTED_VERSION_KEY] = "1.2.3";
  originNode[EXPECTED_USER_KEY] = "jdoe";
  try
  {
    axom::sina::Run run {originNode};
    FAIL() << "Application should be missing, but is " << run.getApplication();
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_APPLICATION_KEY));
  }
}

TEST(Run, toNode)
{
  ID id {"the id", IDType::Global};
  axom::sina::Run run {id, "the app", "1.2.3", "jdoe"};
  auto asNode = run.toNode();
  EXPECT_TRUE(asNode.dtype().is_object());
  EXPECT_EQ("run", asNode[EXPECTED_TYPE_KEY].as_string());
  EXPECT_EQ("the id", asNode[EXPECTED_GLOBAL_ID_KEY].as_string());
  EXPECT_TRUE(asNode[EXPECTED_LOCAL_ID_KEY].dtype().is_empty());
  EXPECT_EQ("the app", asNode[EXPECTED_APPLICATION_KEY].as_string());
  EXPECT_EQ("1.2.3", asNode[EXPECTED_VERSION_KEY].as_string());
  EXPECT_EQ("jdoe", asNode[EXPECTED_USER_KEY].as_string());
}

TEST(Run, addRunLoader)
{
  conduit::Node originNode;
  originNode[EXPECTED_TYPE_KEY] = "run";
  originNode[EXPECTED_GLOBAL_ID_KEY] = "the id";
  originNode[EXPECTED_APPLICATION_KEY] = "the app";
  originNode[EXPECTED_VERSION_KEY] = "1.2.3";
  originNode[EXPECTED_USER_KEY] = "jdoe";

  RecordLoader loader;
  addRunLoader(loader);

  auto record = loader.load(originNode);
  auto run = dynamic_cast<axom::sina::Run *>(record.get());
  ASSERT_NE(nullptr, run);
  EXPECT_EQ("run", run->getType());
  EXPECT_EQ("the id", run->getId().getId());
  EXPECT_EQ(IDType::Global, run->getId().getType());
  EXPECT_EQ("the app", run->getApplication());
  EXPECT_EQ("1.2.3", run->getVersion());
  EXPECT_EQ("jdoe", run->getUser());
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
