// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina/core/AdiakWriter.hpp"

#ifdef AXOM_USE_ADIAK

  #include <stdexcept>
  #include <string>
  #include <vector>

  #include "gtest/gtest.h"
  #include "gmock/gmock.h"

  #include "adiak.hpp"
extern "C" {
  #include "adiak_tool.h"
  #include "adiak.h"
}

  #include "axom/sina/core/Datum.hpp"
  #include "axom/sina/core/ID.hpp"
  #include "axom/sina/core/Run.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::DoubleEq;
using ::testing::ElementsAre;
using ::testing::HasSubstr;

char const EXPECTED_DATA_KEY[] = "data";
char const EXPECTED_FILES_KEY[] = "files";

class AdiakWriterTest : public ::testing::Test
{
protected:
  static void SetUpTestCase()
  {
    adiak::init(nullptr);
    adiak_register_cb(1, adiak_category_all, AdiakWriterTest::callbackWrapper, 0, &current_test);
  }

  void SetUp() override { current_test = this; }

  static void callbackWrapper(const char *name,
                              adiak_category_t category,
                              const char *subcategory,
                              adiak_value_t *val,
                              adiak_datatype_t *adiak_type,
                              void *adiakwriter)
  {
    auto test = static_cast<AdiakWriterTest **>(adiakwriter);
    adiakSinaCallback(name, category, subcategory, val, adiak_type, &((*test)->record));
  }

  axom::sina::Record record {axom::sina::ID {"test_run", axom::sina::IDType::Local}, "test_type"};
  static AdiakWriterTest *current_test;
};

AdiakWriterTest *AdiakWriterTest::current_test;

TEST_F(AdiakWriterTest, basic_assignment)
{
  //adiak::init(nullptr);
  //adiak_register_cb(1, adiak_category_all, sina::adiakSinaCallback, 0, callback_record_ptr);
  std::string name1 = "name1";
  std::string value1 = "value1";
  std::vector<std::string> tags1 = {"string"};
  std::string name2 = "name2";
  int value2 = 2;
  std::vector<std::string> tags2 = {"int"};
  auto result1 = adiak::value(name1, value1);
  auto result2 = adiak::value(name2, value2);
  EXPECT_TRUE(result1 && result2);
  auto asNode = record.toNode();
  EXPECT_EQ(value1, asNode[EXPECTED_DATA_KEY][name1]["value"].as_string());
  EXPECT_EQ(value2, int(asNode[EXPECTED_DATA_KEY][name2]["value"].as_double()));
  EXPECT_EQ(tags1[0], asNode[EXPECTED_DATA_KEY][name1]["tags"][0].as_string());
  EXPECT_EQ(tags2[0], asNode[EXPECTED_DATA_KEY][name2]["tags"][0].as_string());
}

TEST_F(AdiakWriterTest, scalar_types)
{
  std::string name1 = "my_long";
  long value1 = 0;
  std::string name2 = "my_double";
  double value2 = 3.14;
  auto result1 = adiak::value(name1, value1);
  auto result2 = adiak::value(name2, value2);
  EXPECT_TRUE(result1 && result2);
  auto asNode = record.toNode();
  EXPECT_EQ(value1, asNode[EXPECTED_DATA_KEY][name1]["value"].as_float64());
  EXPECT_EQ(value2, asNode[EXPECTED_DATA_KEY][name2]["value"].as_double());
}

// No extra test for string_types (besides date) as they're handled identically
TEST_F(AdiakWriterTest, date_type)
{
  std::string name1 = "my_date";
  auto result = adiak::value(name1, adiak::date(1568397849));
  EXPECT_TRUE(result);
  auto toNode = record.toNode();
  EXPECT_EQ("Fri, 13 Sep 2019 11:04:09 -0700", toNode[EXPECTED_DATA_KEY][name1]["value"].as_string());
}

TEST_F(AdiakWriterTest, list_types)
{
  std::string name1 = "my_scalar_list";
  std::vector<double> value1 {4.5, 0, 5.12, 42};
  std::string name2 = "my_string_list";
  std::set<std::string> value2 {"spam", "egg and bacon", "egg and spam"};
  auto result1 = adiak::value(name1, value1);
  auto result2 = adiak::value(name2, value2);
  EXPECT_TRUE(result1 && result2);
  auto asNode = record.toNode();
  auto doub_array = asNode[EXPECTED_DATA_KEY][name1]["value"].as_double_ptr();
  std::vector<double> scal_child_vals(
    doub_array,
    doub_array + asNode[EXPECTED_DATA_KEY][name1]["value"].dtype().number_of_elements());
  EXPECT_EQ(value1, scal_child_vals);
  std::set<std::string> node_vals;
  auto val_itr = asNode[EXPECTED_DATA_KEY][name2]["value"].children();
  while(val_itr.has_next()) node_vals.insert(val_itr.next().as_string());
  EXPECT_EQ(value2, node_vals);
}

TEST_F(AdiakWriterTest, files)
{
  std::string name1 = "my_bash";
  std::string value1 = "/bin/bash";
  std::string name2 = "my_cat_pics";
  std::string value2 = "~/pictures/neighbor_cat.png";
  std::vector<std::string> tags2 {name2};
  auto result1 = adiak::value(name1, adiak::path(value1));
  auto result2 = adiak::value(name2, adiak::path(value2));
  EXPECT_TRUE(result1 && result2);
  auto asNode = record.toNode();
  EXPECT_FALSE(asNode[EXPECTED_FILES_KEY].child(value1).dtype().is_empty());
  EXPECT_EQ(1, asNode[EXPECTED_FILES_KEY].child(value2)["tags"].number_of_children());
  EXPECT_EQ(tags2[0], asNode[EXPECTED_FILES_KEY].child(value2)["tags"][0].as_string());
}

TEST_F(AdiakWriterTest, files_list)
{
  std::string fileListName = "my_gecko_pics";
  std::string fileListVal1 = "~/pictures/spike.png";
  std::string fileListVal2 = "~/pictures/sandy.png";
  std::vector<adiak::path> fileListAdiak {adiak::path(fileListVal1), adiak::path(fileListVal2)};
  std::vector<std::string> tags = {"string"};
  EXPECT_TRUE(adiak::value(fileListName, fileListAdiak));
  auto asNode = record.toNode();
  EXPECT_FALSE(asNode[EXPECTED_FILES_KEY].child(fileListVal1).dtype().is_empty());
  EXPECT_EQ(1, asNode[EXPECTED_FILES_KEY].child(fileListVal2)["tags"].number_of_children());
  EXPECT_EQ(fileListName, asNode[EXPECTED_FILES_KEY].child(fileListVal2)["tags"][0].as_string());
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom

#endif  // AXOM_USE_ADIAK
