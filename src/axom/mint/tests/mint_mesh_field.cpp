// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/Types.hpp"               // for axom Types
#include "axom/mint/mesh/Field.hpp"          // for mint::Field
#include "axom/mint/mesh/FieldTypes.hpp"     // for mint::field_traits
#include "axom/mint/mesh/FieldVariable.hpp"  // for mint::FieldVariable
#include "axom/mint/fem/FEBasis.hpp"         // for FEBasisTypes

#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

// namespace aliases
namespace mint = axom::mint;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template <typename T>
void check_field_instantiation()
{
  constexpr axom::IndexType NUM_TUPLES = 4;
  constexpr axom::IndexType NUM_COMPONENTS = 1;

  mint::Field* f = new mint::FieldVariable<T>("f", NUM_TUPLES, NUM_COMPONENTS);

  EXPECT_TRUE(f != nullptr);
  EXPECT_EQ(f->getName(), "f");
  EXPECT_EQ(f->getNumTuples(), NUM_TUPLES);
  EXPECT_EQ(f->getNumComponents(), NUM_COMPONENTS);
  EXPECT_EQ(f->getBasis(), MINT_UNDEFINED_BASIS);
  EXPECT_EQ(f->getType(), mint::field_traits<T>::type());

  T* data = mint::Field::getDataPtr<T>(f);
  EXPECT_TRUE(data != nullptr);

  delete f;
  f = nullptr;
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(mintMesh_field_DeathTest, invalid_operations)
{
  const char* IGNORE_OUTPUT = ".*";
  constexpr axom::IndexType NUM_TUPLES = 4;
  constexpr axom::IndexType NUM_COMPONENTS = 1;

  mint::Field* f = new mint::FieldVariable<double>("f", NUM_TUPLES, NUM_COMPONENTS);
  EXPECT_DEATH_IF_SUPPORTED(mint::Field::getDataPtr<int>(f), IGNORE_OUTPUT);

  delete f;
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field, instantiate)
{
  check_field_instantiation<double>();
  check_field_instantiation<float>();
  check_field_instantiation<std::int32_t>();
  check_field_instantiation<std::int64_t>();
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field, set_basis)
{
  constexpr axom::IndexType NUM_TUPLES = 4;
  constexpr axom::IndexType NUM_COMPONENTS = 1;

  mint::Field* f = new mint::FieldVariable<double>("f", NUM_TUPLES, NUM_COMPONENTS);
  EXPECT_TRUE(f != nullptr);
  EXPECT_EQ(f->getBasis(), MINT_UNDEFINED_BASIS);

  f->setBasis(MINT_LAGRANGE_BASIS);
  EXPECT_EQ(f->getBasis(), MINT_LAGRANGE_BASIS);

  delete f;
  f = nullptr;
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field, get_dataptr)
{
  constexpr axom::IndexType NUM_TUPLES = 4;
  constexpr axom::IndexType NUM_COMPONENTS = 1;

  double f[] = {1.0, 2.0, 3.0, 4.0};
  mint::Field* field = new mint::FieldVariable<double>("f", f, NUM_TUPLES, NUM_COMPONENTS);

  const double* data = mint::Field::getDataPtr<double>(field);
  EXPECT_TRUE(data != nullptr);
  EXPECT_EQ(data, f);

  for(int i = 0; i < NUM_TUPLES; ++i)
  {
    EXPECT_DOUBLE_EQ(data[i], f[i]);
  }

  delete field;
  field = nullptr;

  for(int i = 0; i < NUM_TUPLES; ++i)
  {
    EXPECT_DOUBLE_EQ(f[i], static_cast<double>(i + 1));
  }
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
