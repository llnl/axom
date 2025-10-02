// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/sidre.hpp"

#include "gtest/gtest.h"

using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::IndexType;
#ifdef AXOM_USE_MFEM
using axom::sidre::MFEMSidreDataCollection;
#endif
using axom::sidre::View;

// This test is meant to mock how Serac stores and loads Quadrature Data.
// The data is stored in Axom Array's of material states and are a user defined
// structure of POD types.
//
// There are two types of tests, one is where the data is managed by Sidre and
// the other is where the data is external to sidre and is saved and loaded to
// disk. The latter is what Serac is doing.

//-----Mock Serac Structs------------------------------------------------------

/**
 * @brief Arbitrary-rank tensor class
 * @tparam T The type stored at each index
 * @tparam n The dimensions of the tensor
 */
template <typename T, int... n>
struct tensor;

template <typename T, int m, int... n>
struct tensor<T, m, n...>
{
  template <typename i_type>
  constexpr auto& operator()(i_type i)
  {
    return data[i];
  }
  template <typename i_type>
  constexpr auto& operator()(i_type i) const
  {
    return data[i];
  }
  template <typename i_type, typename... jklm_type>
  constexpr auto& operator()(i_type i, jklm_type... jklm)
  {
    return data[i](jklm...);
  }
  template <typename i_type, typename... jklm_type>
  constexpr auto& operator()(i_type i, jklm_type... jklm) const
  {
    return data[i](jklm...);
  }

  constexpr auto& operator[](int i) { return data[i]; }
  constexpr const auto& operator[](int i) const { return data[i]; }

  tensor<T, n...> data[m];
};

template <typename T, int m>
struct tensor<T, m>
{
  template <typename i_type>
  constexpr auto& operator()(i_type i)
  {
    return data[i];
  }
  template <typename i_type>
  constexpr auto& operator()(i_type i) const
  {
    return data[i];
  }
  constexpr auto& operator[](int i) { return data[i]; }
  constexpr const auto& operator[](int i) const { return data[i]; }

  template <int last_dimension = m, typename = typename std::enable_if<last_dimension == 1>::type>
  constexpr operator T()
  {
    return data[0];
  }

  template <int last_dimension = m, typename = typename std::enable_if<last_dimension == 1>::type>
  constexpr operator T() const
  {
    return data[0];
  }

  T data[m];
};

//---------------------

template <typename T>
bool check(const T&, double)
{
  SLIC_ERROR("You didn't implement a specialization for check");
  return false;
}

template <typename T>
void fill(T&, double)
{
  SLIC_ERROR("You didn't implement a specialization for fill");
}

//---------------------

template <>
bool check(const double& state, double value)
{
  SLIC_INFO(axom::fmt::format("{} == {}", state, value));
  return state == value;
}

template <>
void fill(double& state, double value)
{
  state = value;
}

//---------------------

struct StateOne
{
  double x;
};

template <>
bool check(const StateOne& state, double value)
{
  return state.x == value;
}

template <>
void fill(StateOne& state, double value)
{
  state.x = value;
}

//---------------------

struct StateTwo
{
  double x;
  double y;
};

template <>
bool check(const StateTwo& state, double value)
{
  return (state.x == value) && (state.y == (value + 1));
}

template <>
void fill(StateTwo& state, double value)
{
  state.x = value;
  state.y = value + 1;
}

//---------------------

struct StateThree
{
  double x;
  double y;
  double z;
};

template <>
bool check(const StateThree& state, double value)
{
  return (state.x == value) && (state.y == (value + 1)) && (state.z == (value + 2));
}

template <>
void fill(StateThree& state, double value)
{
  state.x = value;
  state.y = value + 1;
  state.z = value + 2;
}

//---------------------

struct StateTensorSmall
{
  tensor<double, 2> t;
  double x;
};

template <>
bool check(const StateTensorSmall& state, double value)
{
  return (state.t(0) == value) && (state.t(1) == (value + 1)) && (state.x == (value + 4));
}

template <>
void fill(StateTensorSmall& state, double value)
{
  state.t(0) = value;
  state.t(1) = value + 1;
  state.x = value + 4;
}

//---------------------

struct StateTensorLarge
{
  tensor<double, 2, 2> t;
  double x;
};

template <>
bool check(const StateTensorLarge& state, double value)
{
  return (state.t(0, 0) == value) && (state.t(0, 1) == (value + 1)) &&
    (state.t(1, 0) == (value + 2)) && (state.t(1, 1) == (value + 3)) && (state.x == (value + 4));
}

template <>
void fill(StateTensorLarge& state, double value)
{
  state.t(0, 0) = value;
  state.t(0, 1) = value + 1;
  state.t(1, 0) = value + 2;
  state.t(1, 1) = value + 3;
  state.x = value + 4;
}

//---------------------

template <typename T, int DIM>
void fill_array(T, int, bool)
{
  SLIC_ERROR("You didn't implement the fill array of this type.");
}

template <typename T, int DIM>
void check_array(T)
{
  SLIC_ERROR("You didn't implement the check array of this type.");
}

template <typename T, int DIM>
void fill_array(axom::Array<T, DIM>& states, double fill_value, bool increment)
{
  int count = 0;
  for(auto& state : states)
  {
    fill(state, fill_value);
    if(increment)
    {
      fill_value++;
    }
    count++;
  }
  SLIC_INFO(axom::fmt::format("filled {} states with end value {}", count, fill_value));
}

template <typename T, int DIM>
void check_array(axom::Array<T, DIM> states)
{
  int count = 0;
  double check_value = 0;
  for(auto& state : states)
  {
    EXPECT_TRUE(check(state, check_value++));
    count++;
  }
  SLIC_INFO(axom::fmt::format("checked {} states with end value {}", count, check_value));
}

//------------------------------------------------------------------------------

template <typename T, int DIM>
void test_user_defined_data()
{
  // populate data
  constexpr IndexType size = 10;
  axom::Array<T, DIM> states(size, size);
  fill_array(states, 0, true);

  // Create datastore
  DataStore ds;
  Group* root = ds.getRoot();

  // get size
  auto num_states = static_cast<IndexType>(states.size());
  auto state_size = static_cast<IndexType>(sizeof(*(states.begin())));
  auto total_size = num_states * state_size;

  SLIC_INFO(axom::fmt::format("Num of States={}", num_states));
  SLIC_INFO(axom::fmt::format("State Size={}", state_size));
  SLIC_INFO(axom::fmt::format("Total Size={}", total_size));
  SLIC_INFO(axom::fmt::format("Total Size/State Size={}", total_size / state_size));

  // write shape
  root->createViewScalar("num_states", num_states);
  root->createViewScalar("state_size", state_size);
  root->createViewScalar("total_size", total_size);

  // write data to datastore as bytes
  View* states_view = root->createViewAndAllocate("states", axom::sidre::UINT8_ID, total_size);
  std::uint8_t* sidre_state_data = states_view->getData();
  memcpy(sidre_state_data, states.data(), static_cast<std::size_t>(total_size));

  // mess with data before overriding it back to original in sidre
  fill_array(states, -1, false);

  // Copy original data back over local data
  memcpy(states.data(), sidre_state_data, static_cast<std::size_t>(total_size));

  // Test data is back to original
  check_array(states);
}

template <typename T, int DIM>
void test_external_user_defined_data()
{
#ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::Group::loadExternalData() is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
#endif

  constexpr IndexType size = 10;
  IndexType num_states {0};
  IndexType state_size {0};
  IndexType total_size {0};
  IndexType num_uint8s {0};

  DataStore ds;
  Group* root = ds.getRoot();

  const std::string filename = "sidre_external_states";

  // populate data
  {
    axom::Array<T, DIM> states(size, size);
    fill_array(states, 0, true);

    // get size
    num_states = static_cast<IndexType>(states.size());
    state_size = static_cast<IndexType>(sizeof(*(states.begin())));
    total_size = num_states * state_size;
    num_uint8s = total_size / sizeof(std::uint8_t);

    SLIC_INFO(axom::fmt::format("Num of States={}", num_states));
    SLIC_INFO(axom::fmt::format("State Size={}", state_size));
    SLIC_INFO(axom::fmt::format("Total Size={}", total_size));
    SLIC_INFO(axom::fmt::format("Total Size/State Size={}", total_size / state_size));
    SLIC_INFO(axom::fmt::format("Num of uint8={}", num_uint8s));

    // write shape
    root->createViewScalar("num_states", num_states);
    root->createViewScalar("state_size", state_size);
    root->createViewScalar("total_size", total_size);

    // Add states as an external buffer
    View* states_view = root->createView("states");
    states_view->setExternalDataPtr(axom::sidre::UINT8_ID, num_uint8s, states.data());

    // Save the array data into a file
    root->save(filename);
  }

  // Load data into new array
  {
    Group* loaded_group = root->createGroup("loaded_data");
    loaded_group->load(filename);

    // Verify size correctness
    EXPECT_TRUE(loaded_group->hasView("num_states"));
    auto loaded_num_states = loaded_group->getView("num_states")->getData<axom::IndexType>();
    EXPECT_EQ(num_states, loaded_num_states);

    EXPECT_TRUE(loaded_group->hasView("state_size"));
    auto loaded_state_size = loaded_group->getView("state_size")->getData<axom::IndexType>();
    EXPECT_EQ(state_size, loaded_state_size);

    EXPECT_TRUE(loaded_group->hasView("total_size"));
    auto loaded_total_size = loaded_group->getView("total_size")->getData<axom::IndexType>();
    EXPECT_EQ(total_size, loaded_total_size);

    EXPECT_TRUE(loaded_group->hasView("states"));
    auto loaded_states_size = loaded_group->getView("states")->getNumElements();
    EXPECT_EQ(num_uint8s, loaded_states_size);

    // Create new array to fill with loaded data
    axom::Array<T, DIM> loaded_states(size, size);
    fill_array(loaded_states, -1, false);

    // Load external data
    View* loaded_states_view = loaded_group->getView("states");
    loaded_states_view->setExternalDataPtr(loaded_states.data());
    loaded_group->loadExternalData(filename);

    // Test data was read into the new location and overwrote the bad data
    check_array(loaded_states);
  }
}

#if defined(AXOM_USE_MFEM) && defined(AXOM_USE_HDF5)

template <typename T, int DIM>
void test_MFEMSidreDataCollection_user_defined_data()
{
  std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  constexpr IndexType size = 10;
  IndexType num_states {0};
  IndexType state_size {0};
  IndexType total_size {0};
  IndexType num_uint8s {0};

  // populate data into Data Collection and save
  {
    axom::Array<T, DIM> states(size, size);
    fill_array(states, 0, true);

    // Create MFEMSidreDataCollection
    auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian1D(10));
    const bool owns_mesh_data = true;
    MFEMSidreDataCollection sdc_writer(test_name, mesh, owns_mesh_data);
    Group* bp_group = sdc_writer.GetBPGroup();

    // get size
    num_states = static_cast<IndexType>(states.size());
    state_size = static_cast<IndexType>(sizeof(*(states.begin())));
    total_size = num_states * state_size;
    num_uint8s = total_size / sizeof(std::uint8_t);

    SLIC_INFO(axom::fmt::format("Num of States={}", num_states));
    SLIC_INFO(axom::fmt::format("State Size={}", state_size));
    SLIC_INFO(axom::fmt::format("Total Size={}", total_size));
    SLIC_INFO(axom::fmt::format("Total Size/State Size={}", total_size / state_size));
    SLIC_INFO(axom::fmt::format("Num of uint8={}", num_uint8s));

    // write shape
    bp_group->createViewScalar("num_states", num_states);
    bp_group->createViewScalar("state_size", state_size);
    bp_group->createViewScalar("total_size", total_size);

    // Add states as an external buffer
    View* states_view = bp_group->createView("states");
    states_view->setExternalDataPtr(axom::sidre::UINT8_ID, num_uint8s, states.data());

  // Save the array data into a file
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    sdc_writer.SetComm(MPI_COMM_WORLD);
  #endif
    sdc_writer.SetCycle(0);
    sdc_writer.Save();
  }

  // Load data into new Data Collection
  {
    MFEMSidreDataCollection sdc_reader(test_name);
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    sdc_reader.SetComm(MPI_COMM_WORLD);
  #endif
    sdc_reader.Load();

    Group* loaded_bp_group = sdc_reader.GetBPGroup();

    // Verify size correctness
    EXPECT_TRUE(loaded_bp_group->hasView("num_states"));
    auto loaded_num_states = loaded_bp_group->getView("num_states")->getData<axom::IndexType>();
    EXPECT_EQ(num_states, loaded_num_states);

    EXPECT_TRUE(loaded_bp_group->hasView("state_size"));
    auto loaded_state_size = loaded_bp_group->getView("state_size")->getData<axom::IndexType>();
    EXPECT_EQ(state_size, loaded_state_size);

    EXPECT_TRUE(loaded_bp_group->hasView("total_size"));
    auto loaded_total_size = loaded_bp_group->getView("total_size")->getData<axom::IndexType>();
    EXPECT_EQ(total_size, loaded_total_size);

    EXPECT_TRUE(loaded_bp_group->hasView("states"));
    auto loaded_states_size = loaded_bp_group->getView("states")->getNumElements();
    EXPECT_EQ(num_uint8s, loaded_states_size);

    // Create new array to fill with loaded data
    axom::Array<T, DIM> loaded_states(size, size);
    fill_array(loaded_states, -1, false);

    // Load external data
    EXPECT_TRUE(loaded_bp_group->hasView("states"));
    View* loaded_states_view = loaded_bp_group->getView("states");
    loaded_states_view->setExternalDataPtr(loaded_states.data());
    sdc_reader.LoadExternalData();

    // Test data was read into the new location and overwrote the bad data
    check_array(loaded_states);
  }
}

#endif  // defined(AXOM_USE_MFEM) && defined(AXOM_USE_HDF5)

//------------------------------------------------------------------------------

TEST(sidre, OneD_double_readandwrite) { test_user_defined_data<double, 1>(); }

TEST(sidre, OneD_StateOne_readandwrite) { test_user_defined_data<StateOne, 1>(); }

TEST(sidre, OneD_StateTwo_readandwrite) { test_user_defined_data<StateTwo, 1>(); }

TEST(sidre, OneD_StateThree_readandwrite) { test_user_defined_data<StateThree, 1>(); }

TEST(sidre, OneD_StateTensorSmall_readandwrite) { test_user_defined_data<StateTensorSmall, 1>(); }

TEST(sidre, OneD_StateTensorLarge_readandwrite) { test_user_defined_data<StateTensorLarge, 1>(); }

//-------------------------

TEST(sidre, TwoD_double_readandwrite) { test_user_defined_data<double, 2>(); }

TEST(sidre, TwoD_StateOne_readandwrite) { test_user_defined_data<StateOne, 2>(); }

TEST(sidre, TwoD_StateTwo_readandwrite) { test_user_defined_data<StateTwo, 2>(); }

TEST(sidre, TwoD_StateThree_readandwrite) { test_user_defined_data<StateThree, 2>(); }

TEST(sidre, TwoD_StateTensorSmall_readandwrite) { test_user_defined_data<StateTensorSmall, 2>(); }

TEST(sidre, TwoD_StateTensorLarge_readandwrite) { test_user_defined_data<StateTensorLarge, 2>(); }

//------------------------------------------------------------------------------

TEST(sidre, OneD_double_external_readandwrite) { test_external_user_defined_data<double, 1>(); }

TEST(sidre, OneD_StateOne_external_readandwrite) { test_external_user_defined_data<StateOne, 1>(); }

TEST(sidre, OneD_StateTwo_external_readandwrite) { test_external_user_defined_data<StateTwo, 1>(); }

TEST(sidre, OneD_StateThree_external_readandwrite)
{
  test_external_user_defined_data<StateThree, 1>();
}

TEST(sidre, OneD_StateTensorSmall_external_readandwrite)
{
  test_external_user_defined_data<StateTensorSmall, 1>();
}

TEST(sidre, OneD_StateTensorLarge_external_readandwrite)
{
  test_external_user_defined_data<StateTensorLarge, 1>();
}

//-------------------------

TEST(sidre, TwoD_double_external_readandwrite) { test_external_user_defined_data<double, 2>(); }

TEST(sidre, TwoD_StateOne_external_readandwrite) { test_external_user_defined_data<StateOne, 2>(); }

TEST(sidre, TwoD_StateTwo_external_readandwrite) { test_external_user_defined_data<StateTwo, 2>(); }

TEST(sidre, TwoD_StateThree_external_readandwrite)
{
  test_external_user_defined_data<StateThree, 2>();
}

TEST(sidre, TwoD_StateTensorSmall_external_readandwrite)
{
  test_external_user_defined_data<StateTensorSmall, 2>();
}

TEST(sidre, TwoD_StateTensorLarge_external_readandwrite)
{
  test_external_user_defined_data<StateTensorLarge, 2>();
}

//-------------------------

#if defined(AXOM_USE_MFEM) && defined(AXOM_USE_HDF5)

TEST(sidre, OneD_double_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<double, 1>();
}

TEST(sidre, OneD_StateOne_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateOne, 1>();
}

TEST(sidre, OneD_StateTwo_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateTwo, 1>();
}

TEST(sidre, OneD_StateThree_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateThree, 1>();
}

TEST(sidre, OneD_StateTensorSmall_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateTensorSmall, 1>();
}

TEST(sidre, OneD_StateTensorLarge_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateTensorLarge, 1>();
}

//-------------------------

TEST(sidre, TwoD_double_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<double, 2>();
}

TEST(sidre, TwoD_StateOne_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateOne, 2>();
}

TEST(sidre, TwoD_StateTwo_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateTwo, 2>();
}

TEST(sidre, TwoD_StateThree_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateThree, 2>();
}

TEST(sidre, TwoD_StateTensorSmall_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateTensorSmall, 2>();
}

TEST(sidre, TwoD_StateTensorLarge_MFEMSidreDataCollection_readandwrite)
{
  test_MFEMSidreDataCollection_user_defined_data<StateTensorLarge, 2>();
}

#endif  // defined(AXOM_USE_MFEM) && defined(AXOM_USE_HDF5)

//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
#endif

  result = RUN_ALL_TESTS();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return result;
}
