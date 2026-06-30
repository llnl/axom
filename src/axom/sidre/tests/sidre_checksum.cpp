// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

#include "conduit_node.hpp"

#include <string>

namespace
{
void populateConduitNode(conduit::Node& node,
                         const std::string& nestedName = "nested",
                         const char* message = "alpha",
                         int count = 42,
                         double arrayTail = 3.75,
                         bool addExtraChild = false)
{
  const double values[] = {1.5, -2.25, arrayTail};

  node["message"].set(message);
  node["count"].set(count);
  node["values"].set(values, 3);
  node[nestedName + "/pi"].set(3.14159);

  if(addExtraChild)
  {
    node["extra"].set(-7);
  }
}
}  // namespace

TEST(sidre_checksum, conduit_checksum_handles_strided_numeric_arrays)
{
  double contiguousValues[] = {1.0, 2.0, 3.0};
  double interleavedValues[] = {1.0, 99.0, 2.0, 99.0, 3.0, 99.0};

  conduit::Node contiguous;
  contiguous.set_external(contiguousValues, 3);

  conduit::Node strided;
  strided.set_external(
    conduit::DataType::float64(3, 0, 2 * sizeof(double)),
    interleavedValues);

  EXPECT_EQ(axom::sidre::checksum(contiguous), axom::sidre::checksum(strided));
}

TEST(sidre_checksum, conduit_checksum_handles_strided_numeric_array_trees)
{
  int contiguousIds[] = {10, 20, 30};
  int interleavedIds[] = {10, -1, 20, -1, 30, -1};
  double contiguousValues[] = {1.5, 2.5, 3.5};
  double interleavedValues[] = {1.5, -1.0, 2.5, -1.0, 3.5, -1.0};

  conduit::Node contiguous;
  contiguous["fields/ids"].set_external(contiguousIds, 3);
  contiguous["fields/values"].set_external(contiguousValues, 3);

  conduit::Node strided;
  strided["fields/ids"].set_external(
    conduit::DataType::c_int(3, 0, 2 * sizeof(int)),
    interleavedIds);
  strided["fields/values"].set_external(
    conduit::DataType::float64(3, 0, 2 * sizeof(double)),
    interleavedValues);

  EXPECT_EQ(axom::sidre::checksum(contiguous), axom::sidre::checksum(strided));
}

TEST(sidre_checksum, conduit_checksum_changes_for_tree_structure_and_leaf_data)
{
  conduit::Node baseline;
  populateConduitNode(baseline);
  const auto checksum = axom::sidre::checksum(baseline);

  conduit::Node stringChanged;
  populateConduitNode(stringChanged, "nested", "beta");
  EXPECT_NE(checksum, axom::sidre::checksum(stringChanged));

  conduit::Node scalarChanged;
  populateConduitNode(scalarChanged, "nested", "alpha", 7);
  EXPECT_NE(checksum, axom::sidre::checksum(scalarChanged));

  conduit::Node arrayChanged;
  populateConduitNode(arrayChanged, "nested", "alpha", 42, 9.25);
  EXPECT_NE(checksum, axom::sidre::checksum(arrayChanged));

  conduit::Node renamedChild;
  populateConduitNode(renamedChild, "renamed_nested");
  EXPECT_NE(checksum, axom::sidre::checksum(renamedChild));

  conduit::Node addedChild;
  populateConduitNode(addedChild, "nested", "alpha", 42, 3.75, true);
  EXPECT_NE(checksum, axom::sidre::checksum(addedChild));
}

TEST(sidre_checksum, conduit_checksum_supports_empty_object_and_list_nodes)
{
  conduit::Node emptyObject;
  emptyObject.set(conduit::DataType::object());

  conduit::Node emptyList;
  emptyList.set(conduit::DataType::list());

  EXPECT_EQ(0.0L, axom::sidre::checksum(emptyObject, false));
  EXPECT_EQ(0.0L, axom::sidre::checksum(emptyList, false));
  EXPECT_EQ(axom::sidre::checksum(emptyObject, false),
            axom::sidre::checksum(emptyList, false));

  conduit::Node namedEmptyObject;
  namedEmptyObject["empty_object"].set(conduit::DataType::object());

  conduit::Node namedEmptyList;
  namedEmptyList["empty_list"].set(conduit::DataType::list());

  EXPECT_NE(axom::sidre::checksum(namedEmptyObject["empty_object"]),
            axom::sidre::checksum(namedEmptyList["empty_list"]));
}

TEST(sidre_checksum, sidre_view_and_group_checksum_support_strided_external_data)
{
  axom::sidre::DataStore datastore;
  axom::sidre::Group* group = datastore.getRoot()->createGroup("fields");

  int interleavedData[] = {11, -1, 22, -1, 33, -1};
  int contiguousData[] = {11, 22, 33};

  axom::sidre::View* stridedView = group->createView("strided");
  stridedView->setExternalDataPtr(interleavedData);
  stridedView->apply(conduit::DataType::c_int(3, 0, 2 * sizeof(int)));

  conduit::Node contiguousNode;
  contiguousNode.set_external(contiguousData, 3);

  axom::ArrayView<const char> nameView(stridedView->getName().data(),
                                       stridedView->getName().size());
  const auto expectedViewChecksum =
    axom::utilities::checksum(nameView) + axom::sidre::checksum(contiguousNode);

  EXPECT_EQ(expectedViewChecksum, stridedView->checksum());

  const auto groupChecksumBefore = group->checksum();
  interleavedData[2] += 7;
  EXPECT_NE(groupChecksumBefore, group->checksum());
}

TEST(sidre_checksum, view_checksum_changes_on_rename_and_data_mutation)
{
  axom::sidre::DataStore datastore;
  axom::sidre::Group* root = datastore.getRoot();
  axom::sidre::View* view =
    root->createViewAndAllocate("values", axom::sidre::INT_ID, 4);

  int* data = view->getData<int*>();
  data[0] = 2;
  data[1] = 4;
  data[2] = 6;
  data[3] = 8;

  const auto originalChecksum = view->checksum();

  data[2] += 5;
  const auto dataChecksum = view->checksum();
  EXPECT_NE(originalChecksum, dataChecksum);

  ASSERT_TRUE(view->rename("renamed_values"));
  EXPECT_NE(dataChecksum, view->checksum());
}

TEST(sidre_checksum, view_checksum_matches_unnamed_native_layout_path)
{
  axom::sidre::DataStore datastore;
  axom::sidre::View* view =
    datastore.getRoot()->createViewAndAllocate("values", axom::sidre::INT_ID, 3);

  int* data = view->getData<int*>();
  data[0] = 5;
  data[1] = 8;
  data[2] = 13;

  conduit::Node nativeLayout;
  view->createNativeLayout(nativeLayout);

  axom::ArrayView<const char> nameView(view->getName().data(), view->getName().size());
  const auto oldPathChecksum =
    axom::utilities::checksum(nameView) + axom::sidre::checksum(nativeLayout);

  EXPECT_EQ(oldPathChecksum, view->checksum());
}

TEST(sidre_checksum, group_checksum_changes_on_add_remove_rename_and_descendant_data)
{
  axom::sidre::DataStore datastore;
  axom::sidre::Group* group = datastore.getRoot()->createGroup("fields");

  const auto emptyChecksum = group->checksum();

  group->createViewScalar("flag", 7);
  const auto viewAddedChecksum = group->checksum();
  EXPECT_NE(emptyChecksum, viewAddedChecksum);

  group->destroyView("flag");
  EXPECT_EQ(emptyChecksum, group->checksum());

  axom::sidre::Group* child = group->createGroup("child");
  const auto childAddedChecksum = group->checksum();
  EXPECT_NE(emptyChecksum, childAddedChecksum);

  ASSERT_TRUE(child->rename("child_renamed"));
  const auto childRenamedChecksum = group->checksum();
  EXPECT_NE(childAddedChecksum, childRenamedChecksum);

  axom::sidre::View* values =
    child->createViewAndAllocate("values", axom::sidre::INT_ID, 4);
  int* data = values->getData<int*>();
  data[0] = 1;
  data[1] = 3;
  data[2] = 5;
  data[3] = 7;

  const auto beforeMutationChecksum = group->checksum();
  data[1] += 10;
  const auto afterMutationChecksum = group->checksum();
  EXPECT_NE(beforeMutationChecksum, afterMutationChecksum);

  group->destroyGroupAndData("child_renamed");
  EXPECT_NE(afterMutationChecksum, group->checksum());
  EXPECT_EQ(emptyChecksum, group->checksum());
}
