// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <stdexcept>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/Relationship.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

char const EXPECTED_GLOBAL_OBJECT_ID_KEY[] = "object";
char const EXPECTED_LOCAL_OBJECT_ID_KEY[] = "local_object";
char const EXPECTED_GLOBAL_SUBJECT_ID_KEY[] = "subject";
char const EXPECTED_LOCAL_SUBJECT_ID_KEY[] = "local_subject";
char const EXPECTED_PREDICATE_KEY[] = "predicate";

using ::testing::HasSubstr;

TEST(Relationship, create)
{
  std::string subjectID = "the subject";
  std::string objectID = "the object";
  std::string predicate = "is somehow related to";

  Relationship relationship {ID {subjectID, IDType::Global}, predicate, ID {objectID, IDType::Local}};

  EXPECT_EQ(subjectID, relationship.getSubject().getId());
  EXPECT_EQ(IDType::Global, relationship.getSubject().getType());
  EXPECT_EQ(objectID, relationship.getObject().getId());
  EXPECT_EQ(IDType::Local, relationship.getObject().getType());
  EXPECT_EQ(predicate, relationship.getPredicate());
}

TEST(Relationship, create_fromNode_validGlobalIDs)
{
  std::string subjectID = "the subject";
  std::string objectID = "the object";
  std::string predicate = "is somehow related to";

  conduit::Node asNode;
  asNode[EXPECTED_GLOBAL_SUBJECT_ID_KEY] = subjectID;
  asNode[EXPECTED_GLOBAL_OBJECT_ID_KEY] = objectID;
  asNode[EXPECTED_PREDICATE_KEY] = predicate;

  Relationship relationship {asNode};

  EXPECT_EQ(subjectID, relationship.getSubject().getId());
  EXPECT_EQ(IDType::Global, relationship.getSubject().getType());
  EXPECT_EQ(objectID, relationship.getObject().getId());
  EXPECT_EQ(IDType::Global, relationship.getObject().getType());
  EXPECT_EQ(predicate, relationship.getPredicate());
}

TEST(Relationship, create_from_validLocalIDs)
{
  std::string subjectID = "the subject";
  std::string objectID = "the object";
  std::string predicate = "is somehow related to";

  conduit::Node asNode;
  asNode[EXPECTED_LOCAL_SUBJECT_ID_KEY] = subjectID;
  asNode[EXPECTED_LOCAL_OBJECT_ID_KEY] = objectID;
  asNode[EXPECTED_PREDICATE_KEY] = predicate;

  Relationship relationship {asNode};

  EXPECT_EQ(subjectID, relationship.getSubject().getId());
  EXPECT_EQ(IDType::Local, relationship.getSubject().getType());
  EXPECT_EQ(objectID, relationship.getObject().getId());
  EXPECT_EQ(IDType::Local, relationship.getObject().getType());
  EXPECT_EQ(predicate, relationship.getPredicate());
}

TEST(Relationship, create_fromNode_missingSubect)
{
  conduit::Node asNode;
  asNode[EXPECTED_LOCAL_OBJECT_ID_KEY] = "the object";
  asNode[EXPECTED_PREDICATE_KEY] = "some predicate";
  try
  {
    Relationship relationship {asNode};
    FAIL() << "Should have gotten an exception about a missing subject";
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_LOCAL_SUBJECT_ID_KEY));
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_GLOBAL_SUBJECT_ID_KEY));
  }
}

TEST(Relationship, create_fromNode_missingObject)
{
  conduit::Node asNode;
  asNode[EXPECTED_LOCAL_SUBJECT_ID_KEY] = "the subject";
  asNode[EXPECTED_PREDICATE_KEY] = "some predicate";

  try
  {
    Relationship relationship {asNode};
    FAIL() << "Should have gotten an exception about a missing object";
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_LOCAL_OBJECT_ID_KEY));
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_GLOBAL_OBJECT_ID_KEY));
  }
}

TEST(Relationship, create_fromNode_missingPredicate)
{
  conduit::Node asNode;
  asNode[EXPECTED_LOCAL_SUBJECT_ID_KEY] = "the subject";
  asNode[EXPECTED_LOCAL_OBJECT_ID_KEY] = "the object";

  try
  {
    Relationship relationship {asNode};
    FAIL() << "Should have gotten an exception about a missing predicate";
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_PREDICATE_KEY));
    EXPECT_THAT(expected.what(), HasSubstr("Relationship"));
  }
}

TEST(Relationship, toNode_localIds)
{
  std::string subjectID = "the subject";
  std::string objectID = "the object";
  std::string predicate = "is somehow related to";

  Relationship relationship {ID {subjectID, IDType::Local}, predicate, ID {objectID, IDType::Local}};

  conduit::Node asNode = relationship.toNode();

  EXPECT_EQ(subjectID, asNode[EXPECTED_LOCAL_SUBJECT_ID_KEY].as_string());
  EXPECT_EQ(objectID, asNode[EXPECTED_LOCAL_OBJECT_ID_KEY].as_string());
  EXPECT_EQ(predicate, asNode[EXPECTED_PREDICATE_KEY].as_string());
  EXPECT_FALSE(asNode.has_child(EXPECTED_GLOBAL_SUBJECT_ID_KEY));
  EXPECT_FALSE(asNode.has_child(EXPECTED_GLOBAL_OBJECT_ID_KEY));
}

TEST(Relationship, toNode_globalIds)
{
  std::string subjectID = "the subject";
  std::string objectID = "the object";
  std::string predicate = "is somehow related to";

  Relationship relationship {ID {subjectID, IDType::Global}, predicate, ID {objectID, IDType::Global}};

  conduit::Node asNode = relationship.toNode();

  EXPECT_EQ(subjectID, asNode[EXPECTED_GLOBAL_SUBJECT_ID_KEY].as_string());
  EXPECT_EQ(objectID, asNode[EXPECTED_GLOBAL_OBJECT_ID_KEY].as_string());
  EXPECT_EQ(predicate, asNode[EXPECTED_PREDICATE_KEY].as_string());
  EXPECT_FALSE(asNode.has_child(EXPECTED_LOCAL_SUBJECT_ID_KEY));
  EXPECT_FALSE(asNode.has_child(EXPECTED_LOCAL_OBJECT_ID_KEY));
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom