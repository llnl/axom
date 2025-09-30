// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/mint/mesh/FieldVariable.hpp"  // for mint::FieldVariable
#include "axom/mint/mesh/FieldTypes.hpp"     // for FieldTypes enum
#include "axom/core/numerics/Matrix.hpp"     // for numerics::Matrix
#include "axom/slic.hpp"

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
namespace sidre = axom::sidre;
#endif

// gtest includes
#include "gtest/gtest.h"

// namespace aliases
namespace mint = axom::mint;
namespace numerics = axom::numerics;
namespace utilities = axom::utilities;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
template <typename T>
void populate_array(axom::deprecated::MCArray<T>& data)
{
  const axom::IndexType numTuples = data.size();
  const axom::IndexType numComponents = data.numComponents();

  for(axom::IndexType i = 0; i < numTuples; ++i)
  {
    const double offset = static_cast<T>(i * numComponents);
    for(axom::IndexType j = 0; j < numComponents; ++j)
    {
      data(i, j) = offset + j;
    }  // END for all j
  }  // END for all i
}

//------------------------------------------------------------------------------
template <typename T>
void check_array(axom::deprecated::MCArray<T>& data)
{
  const axom::IndexType numTuples = data.size();
  const axom::IndexType numComponents = data.numComponents();

  for(axom::IndexType i = 0; i < numTuples; ++i)
  {
    const T offset = static_cast<T>(i * numComponents);
    for(axom::IndexType j = 0; j < numComponents; ++j)
    {
      const double expected_value = offset + j;
      EXPECT_DOUBLE_EQ(data(i, j), expected_value);
    }  // END for all j
  }  // END for all i
}

//------------------------------------------------------------------------------
template <typename T>
void populate_field_variable(mint::FieldVariable<T>& fv)
{
  const axom::IndexType numTuples = fv.getNumTuples();
  const axom::IndexType numComponents = fv.getNumComponents();

  T* data = fv.getFieldVariablePtr();

  for(axom::IndexType i = 0; i < numTuples; ++i)
  {
    const T offset = static_cast<T>(i * numComponents);
    for(axom::IndexType j = 0; j < numComponents; ++j)
    {
      const double expected_value = offset + j;
      const int idx = static_cast<int>(offset) + j;
      data[idx] = expected_value;
    }  // END for all components
  }  // END for all tuples
}

//------------------------------------------------------------------------------
template <typename T>
void check_field_variable(mint::FieldVariable<T>& fv)
{
  const axom::IndexType numTuples = fv.getNumTuples();
  const axom::IndexType numComponents = fv.getNumComponents();

  const T* data = fv.getFieldVariablePtr();
  EXPECT_TRUE(data != nullptr);

  for(int i = 0; i < numTuples; ++i)
  {
    const T offset = static_cast<T>(i * numComponents);
    for(int j = 0; j < numComponents; ++j)
    {
      const double expected_value = offset + j;
      const int idx = static_cast<int>(offset) + j;
      EXPECT_DOUBLE_EQ(data[idx], expected_value);
    }  // END for all jj
  }  // END for all ii
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE
void create_sidre_data(sidre::DataStore& ds, int numTuples, int numComponents)
{
  sidre::Group* gp = ds.getRoot();
  EXPECT_TRUE(gp != nullptr);

  // create view to hold field values
  sidre::View* values_view = gp->createView("values");
  sidre::deprecated::MCArray<double> data(values_view, numTuples, numComponents);

  // fill in with some data
  populate_array(data);
}
#endif

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable_DeathTest, invalid_construction)
{
  const char* EMPTY_STRING = "";
  const char* IGNORE_OUTPUT = ".*";
  struct invalid_type
  { };

  EXPECT_EQ(mint::field_traits<invalid_type>::type(), mint::UNDEFINED_FIELD_TYPE);

  EXPECT_DEATH_IF_SUPPORTED(mint::FieldVariable<invalid_type>("foo",
                                                              axom::deprecated::internal::ZERO,
                                                              axom::deprecated::internal::ZERO),
                            IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldVariable<double>(EMPTY_STRING,
                                                        axom::deprecated::internal::ZERO,
                                                        axom::deprecated::internal::ZERO),
                            IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable_DeathTest, invalid_operations)
{
  const char* IGNORE_OUTPUT = ".*";
  double f[] = {1.0, 2.0, 3.0, 4.0};
  mint::FieldVariable<double> field("f", f, 4, 1);

  EXPECT_DEATH_IF_SUPPORTED(field.resize(10), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(field.reserve(10), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, native_constructor)
{
  constexpr axom::IndexType NUM_TUPLES = 4;
  constexpr axom::IndexType NUM_COMPONENTS = 2;
  mint::FieldVariable<double> field("rho", NUM_TUPLES, NUM_COMPONENTS);
  EXPECT_EQ(field.getName(), "rho");
  EXPECT_EQ(field.getNumTuples(), NUM_TUPLES);
  EXPECT_EQ(field.getNumComponents(), NUM_COMPONENTS);
  EXPECT_TRUE(field.getFieldVariablePtr() != nullptr);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, external_constructor)
{
  constexpr axom::IndexType NUM_TUPLES = 4;
  constexpr axom::IndexType NUM_COMPONENTS = 1;

  const double MAGIC_NUM = 42;
  double f[] = {1.0, 2.0, 3.0, 4.0};

  mint::FieldVariable<double>* field =
    new mint::FieldVariable<double>("f", f, NUM_TUPLES, NUM_COMPONENTS);
  EXPECT_TRUE(field != nullptr);
  EXPECT_EQ(field->getName(), "f");
  EXPECT_EQ(field->getNumTuples(), NUM_TUPLES);
  EXPECT_EQ(field->getNumComponents(), NUM_COMPONENTS);
  EXPECT_TRUE(field->getFieldVariablePtr() != nullptr);
  EXPECT_EQ(field->getFieldVariablePtr(), f);

  // set everything to MAGIC_NUM
  double* ptr = field->getFieldVariablePtr();
  for(int i = 0; i < NUM_TUPLES; ++i)
  {
    ptr[i] = MAGIC_NUM;
  }

  delete field;
  field = nullptr;

  // ensure the external buffer remains persistent
  for(int i = 0; i < NUM_TUPLES; ++i)
  {
    EXPECT_DOUBLE_EQ(f[i], MAGIC_NUM);
  }
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, sidre_push_constructor)
{
  constexpr int MAX_TUPLES = 25;
  constexpr int MAX_COMPONENTS = 3;

  for(int i = 1; i <= MAX_TUPLES; ++i)
  {
    const int NUM_TUPLES = i;

    for(int j = 1; j <= MAX_COMPONENTS; ++j)
    {
      const int NUM_COMPONENTS = j;

      sidre::DataStore ds;
      sidre::Group* gp = ds.getRoot();
      sidre::View* values_view = gp->createView("values");

      // BEGIN SCOPE
      {
        mint::FieldVariable<double> fv("f", values_view, NUM_TUPLES, NUM_COMPONENTS);
        EXPECT_EQ(fv.getName(), "f");
        EXPECT_EQ(fv.getNumTuples(), NUM_TUPLES);
        EXPECT_EQ(fv.getNumComponents(), NUM_COMPONENTS);

        // store data
        populate_field_variable(fv);
      }
      // END SCOPE

      // ensure data remains persistent in Sidre
      sidre::deprecated::MCArray<double> dataArray(values_view);
      EXPECT_EQ(dataArray.size(), NUM_TUPLES);
      EXPECT_EQ(dataArray.numComponents(), NUM_COMPONENTS);

      // check data
      check_array(dataArray);

    }  // END for all j

  }  // END for all i
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, sidre_pull_constructor)
{
  constexpr int MAX_TUPLES = 25;
  constexpr int MAX_COMPONENTS = 3;

  for(int i = 1; i <= MAX_TUPLES; ++i)
  {
    const int NUM_TUPLES = i;

    for(int j = 1; j <= MAX_COMPONENTS; ++j)
    {
      const int NUM_COMPONENTS = j;

      sidre::DataStore ds;
      create_sidre_data(ds, NUM_TUPLES, NUM_COMPONENTS);

      sidre::View* values_view = ds.getRoot()->getView("values");
      EXPECT_TRUE(values_view != nullptr);

      // BEGIN SCOPE
      {
        mint::FieldVariable<double> fv("f", values_view);
        EXPECT_EQ(fv.getName(), "f");
        EXPECT_EQ(fv.getNumTuples(), NUM_TUPLES);
        EXPECT_EQ(fv.getNumComponents(), NUM_COMPONENTS);

        double* data = fv.getFieldVariablePtr();
        EXPECT_TRUE(data != nullptr);

        // check data
        check_field_variable(fv);
      }
      // END SCOPE

      // ensure data remains persistent in Sidre
      sidre::deprecated::MCArray<double> dataArray(values_view);
      EXPECT_EQ(dataArray.size(), NUM_TUPLES);
      EXPECT_EQ(dataArray.numComponents(), NUM_COMPONENTS);

      check_array(dataArray);

    }  // END for all components

  }  // END for all tuples
}

#endif

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, field_array_access)
{
  constexpr axom::IndexType NUM_TUPLES = 4;
  constexpr axom::IndexType NUM_COMPONENTS = 3;
  constexpr axom::IndexType TOTAL_SIZE = 12;

  const double EXPECTED_DATA[] = {10.0, 11.0, 12.0, 20.0, 21.0, 22.0, 30.0, 31.0, 32.0, 40.0, 41.0, 42.0};

  // create a field variable
  mint::FieldVariable<double> field("f", NUM_TUPLES, NUM_COMPONENTS);
  double* data = field.getFieldVariablePtr();
  EXPECT_EQ(field.getName(), "f");
  EXPECT_EQ(field.getNumTuples(), NUM_TUPLES);
  EXPECT_EQ(field.getNumComponents(), NUM_COMPONENTS);
  EXPECT_TRUE(data != nullptr);

  // set the field data
  numerics::Matrix<double> A(NUM_COMPONENTS, NUM_TUPLES, data, true);
  for(axom::IndexType i = 0; i < NUM_TUPLES; ++i)
  {
    const double base = (i + 1) * 10;
    for(axom::IndexType j = 0; j < NUM_COMPONENTS; ++j)
    {
      A(j, i) = base + j;
    }  // END for all components
  }  // END for all tuples

  // check the data
  for(int i = 0; i < TOTAL_SIZE; ++i)
  {
    EXPECT_DOUBLE_EQ(data[i], EXPECTED_DATA[i]);
  }
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, resize)
{
  constexpr axom::IndexType SMALL_NUM_TUPLES = 5;
  constexpr axom::IndexType LARGE_NUM_TUPLES = 250;
  constexpr axom::IndexType NUM_COMPONENTS = 3;

  mint::FieldVariable<double> field("f", SMALL_NUM_TUPLES, NUM_COMPONENTS);
  EXPECT_EQ(field.getName(), "f");
  EXPECT_EQ(field.getNumTuples(), SMALL_NUM_TUPLES);
  EXPECT_EQ(field.getNumComponents(), NUM_COMPONENTS);
  EXPECT_TRUE(field.getCapacity() >= field.getNumTuples());

  field.resize(LARGE_NUM_TUPLES);
  EXPECT_EQ(field.getNumTuples(), LARGE_NUM_TUPLES);
  EXPECT_TRUE(field.getCapacity() >= field.getNumTuples());

  field.resize(SMALL_NUM_TUPLES);
  EXPECT_EQ(field.getNumTuples(), SMALL_NUM_TUPLES);
  EXPECT_TRUE(field.getCapacity() >= field.getNumTuples());
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, reserve)
{
  constexpr axom::IndexType SMALL_NUM_TUPLES = 5;
  constexpr axom::IndexType LARGE_NUM_TUPLES = 250;
  constexpr axom::IndexType NUM_COMPONENTS = 3;

  mint::FieldVariable<double> field("f", SMALL_NUM_TUPLES, NUM_COMPONENTS);
  EXPECT_EQ(field.getName(), "f");
  EXPECT_TRUE(field.getCapacity() >= field.getNumTuples());

  axom::IndexType current_capacity = field.getCapacity();

  // reserve space for more tuples -- should increase capacity
  field.reserve(LARGE_NUM_TUPLES);
  EXPECT_TRUE(field.getCapacity() >= field.getNumTuples());
  EXPECT_TRUE(current_capacity < field.getCapacity());

  // test reserving smaller capacity -- should not change capacity
  current_capacity = field.getCapacity();
  field.reserve(SMALL_NUM_TUPLES);
  EXPECT_EQ(current_capacity, field.getCapacity());
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_variable, shrink)
{
  constexpr axom::IndexType SMALL_NUM_TUPLES = 5;
  constexpr axom::IndexType NUM_COMPONENTS = 3;

  mint::FieldVariable<double> field("f", SMALL_NUM_TUPLES, NUM_COMPONENTS);
  EXPECT_EQ(field.getName(), "f");

  const double ratio = field.getResizeRatio();
  axom::IndexType capacity = static_cast<axom::IndexType>(SMALL_NUM_TUPLES * ratio + 0.5);

  if(capacity < axom::deprecated::MCArray<axom::IndexType>::MIN_DEFAULT_CAPACITY)
  {
    capacity = axom::deprecated::MCArray<axom::IndexType>::MIN_DEFAULT_CAPACITY;
  }
  EXPECT_EQ(field.getCapacity(), capacity);

  field.shrink();
  EXPECT_EQ(field.getCapacity(), field.getNumTuples());
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
