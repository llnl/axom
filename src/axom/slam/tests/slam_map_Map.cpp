// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_Map.cpp
 *
 * \brief Unit tests for Slam's Map
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/core/execution/runtime_policy.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

using SetType = slam::RangeSet<SetPosition, SetElement>;
using BaseSet = slam::Set<SetPosition, SetElement>;
using IntMap = slam::Map<int, BaseSet>;
using RealMap = slam::Map<double, BaseSet>;

template <typename T>
using VecIndirection = policies::STLVectorIndirection<SetPosition, T>;

constexpr SetPosition MAX_SET_SIZE = 10;

template <int S>
using CompileTimeStrideType = policies::CompileTimeStride<int, S>;

using RunTimeStrideType = policies::RuntimeStride<int>;
using OneStrideType = policies::StrideOne<int>;

}  // end anonymous namespace

TEST(slam_map, construct_empty_map)
{
  IntMap m;

  EXPECT_TRUE(m.isValid(true));
}

template <typename T>
bool constructAndTestMap()
{
  SetType s(MAX_SET_SIZE);

  SLIC_INFO("Creating set of size " << s.size());

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string() << " map on the set ");

  slam::Map<T, BaseSet> m(&s);
  EXPECT_TRUE(m.isValid());

  SLIC_INFO("Setting the elements.");
  double multFac = 1.0001;
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    m[idx] = static_cast<T>(idx * multFac);
  }

  SLIC_INFO("Checking the elements.");
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    EXPECT_EQ(m[idx], static_cast<T>(idx * multFac));
    EXPECT_EQ(m(idx), static_cast<T>(idx * multFac));
  }

  EXPECT_TRUE(m.isValid(true));

  m.print();
  return true;
}

TEST(slam_map, construct_int_map) { EXPECT_TRUE(constructAndTestMap<int>()); }

TEST(slam_map, construct_double_map) { EXPECT_TRUE(constructAndTestMap<double>()); }

TEST(slam_map, out_of_bounds)
{
  int defaultElt = 2;

  SetType s(MAX_SET_SIZE);
  IntMap m(&s, defaultElt);

  SLIC_INFO("Testing Map element access -- in bounds");
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    EXPECT_EQ(defaultElt, m[idx]);
  }

  // Test out of bounds
  SLIC_INFO("Testing Map element access " << "-- out of bounds access; Expecting the test to fail");
#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(m[-1], "") << " Accessed element -1 of Map -- out of bounds";
  EXPECT_DEATH_IF_SUPPORTED(m[m.size()], "")
    << " Accessed element " << m.size() << " of Map -- out of bounds";

#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

TEST(slam_map, map_builder)
{
  SLIC_INFO("Testing construction of Map using MapBuilders");

  using DataType = double;
  using MapSet = slam::Set<>;
  using MapType = slam::Map<DataType, MapSet>;
  using MapBuilder = MapType::MapBuilder;

  MapType m(MapBuilder().set(policies::EmptySetTraits<MapSet>::emptySet()));
  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(m.size(), 0);
  EXPECT_EQ(m.stride(), 1);

  SetType s(MAX_SET_SIZE);
  std::vector<DataType> data_arr(s.size());
  for(auto i = 0u; i < data_arr.size(); ++i)
  {
    data_arr[i] = static_cast<DataType>(i * 1.01);
  }

  MapType m2(MapBuilder().set(&s).data(data_arr.data()));
  EXPECT_TRUE(m2.isValid());
  EXPECT_EQ(m2.size(), s.size());
  EXPECT_EQ(m2.stride(), 1);
  for(auto i = 0u; i < data_arr.size(); ++i)
  {
    EXPECT_EQ(m2[i], data_arr[i]);
  }
}

template <typename T, typename StrideType>
void constructAndTestMapWithStride(int stride)
{
  SetType s(MAX_SET_SIZE);

  SLIC_INFO("\nCreating set of size " << s.size());

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("\nCreating " << slam::util::TypeToString<T>::to_string() << " map with stride "
                          << stride << " on the set ");

  using MapType = slam::Map<T, BaseSet, VecIndirection<T>, StrideType>;
  MapType m(&s, 0, stride);
  EXPECT_TRUE(m.isValid());

  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("\nSetting the elements.");
  double multFac = 100.0001;
  double multFac2 = 1.010;
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    for(auto idx2 = 0; idx2 < stride; ++idx2)
    {
      m(idx, idx2) = static_cast<T>(idx * multFac + idx2 * multFac2);
    }
  }

  EXPECT_TRUE(m.isValid());

  SLIC_INFO("\nChecking the elements.");
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    for(auto idx2 = 0; idx2 < stride; ++idx2)
    {
      EXPECT_DOUBLE_EQ(m(idx, idx2), static_cast<T>(idx * multFac + idx2 * multFac2));
    }
  }

  EXPECT_TRUE(m.isValid(true));
}

TEST(slam_map, construct_double_map_with_stride)
{
  constructAndTestMapWithStride<double, RunTimeStrideType>(1);
  constructAndTestMapWithStride<double, RunTimeStrideType>(2);
  constructAndTestMapWithStride<double, RunTimeStrideType>(3);

  constructAndTestMapWithStride<double, CompileTimeStrideType<1>>(1);
  constructAndTestMapWithStride<double, CompileTimeStrideType<2>>(2);
  constructAndTestMapWithStride<double, CompileTimeStrideType<3>>(3);

  constructAndTestMapWithStride<double, OneStrideType>(1);
}

TEST(slam_map, iterate)
{
  using IterType = RealMap::iterator;

  SLIC_INFO("Testing iterator access");

  SetType s(MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating '" << slam::util::TypeToString<double>::to_string() << "' map on the set ");
  RealMap m(&s);
  EXPECT_TRUE(m.isValid());

  SLIC_INFO("Setting the elements.");
  double multFac = 1.0001;
  {
    int idx = 0;
    for(IterType iter = m.begin(); iter != m.end(); iter++)
    {
      *iter = static_cast<double>(idx * multFac);
      idx++;
    }
  }

  SLIC_INFO("Checking the elements by iterator.");
  EXPECT_EQ(RealMap::iterator(0, &m), m.begin());

  //iter++ access
  {
    int idx = 0;
    for(IterType iter = m.begin(); iter != m.end(); iter++)
    {
      EXPECT_EQ(*iter, static_cast<double>(idx * multFac));
      idx++;
    }
    EXPECT_EQ(idx, m.size());
  }

  //iter+n access
  {
    IterType beginIter = m.begin();
    for(int idx = 0; idx < m.size(); ++idx)
    {
      IterType iter = beginIter + idx;
      EXPECT_EQ(*iter, static_cast<double>(idx * multFac));
    }
  }

  //iter-n access
  {
    IterType endIter = m.end();
    for(int idx = 1; idx <= m.size(); ++idx)
    {
      IterType iter = endIter - idx;
      EXPECT_EQ(*iter, static_cast<double>((m.size() - idx) * multFac));
    }
  }

  //iter+=n access
  {
    for(int idx = 0; idx < m.size(); idx++)
    {
      IterType iter = m.begin();
      iter += idx;
      EXPECT_EQ(*iter, static_cast<double>(idx * multFac));
    }
  }

  //iter-=n access
  {
    for(int idx = 1; idx <= m.size(); ++idx)
    {
      IterType iter = m.end();
      iter -= idx;
      EXPECT_EQ(*iter, static_cast<double>((m.size() - idx) * multFac));
    }
  }

  //iter1 - iter2
  {
    IterType beginIter = m.begin();
    IterType endIter = m.end();
    for(int idx = 0; idx < m.size(); idx++)
    {
      IterType iter = beginIter + idx;
      EXPECT_EQ(idx, iter - beginIter);
      EXPECT_EQ(idx - m.size(), iter - endIter);
    }
  }

  EXPECT_TRUE(m.isValid(true));
}

template <typename StrideType>
void constructAndTestMapIteratorWithStride(int stride)
{
  using RealMap = slam::Map<double, slam::Set<>, VecIndirection<double>, StrideType>;

  SetType s(MAX_SET_SIZE);

  SLIC_INFO("\nCreating set of size " << s.size());

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating double map with stride " << stride << " on the set ");

  RealMap m(&s, 0, stride);
  EXPECT_TRUE(m.isValid());

  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements using iterator.");
  double multFac = 100.0001;
  double multFac2 = 1.010;
  {
    int idx = 0;
    for(auto submap : m.set_elements())
    {
      for(auto idx2 = 0; idx2 < submap.size(); ++idx2)
      {
        submap[idx2] = static_cast<double>(idx * multFac + idx2 * multFac2);
      }
      ++idx;
    }
  }

  EXPECT_TRUE(m.isValid());

  SLIC_INFO("Checking the elements by iterator.");

  //iter++ access
  {
    int idx = 0;
    for(auto submap : m.set_elements())
    {
      for(auto idx2 = 0; idx2 < submap.size(); ++idx2)
      {
        EXPECT_DOUBLE_EQ(submap[idx2], static_cast<double>(idx * multFac + idx2 * multFac2));
      }
      idx++;
    }
    EXPECT_EQ(idx, m.size());
  }
}

TEST(slam_map, iterate_with_stride)
{
  constructAndTestMapIteratorWithStride<OneStrideType>(1);

  constructAndTestMapIteratorWithStride<RunTimeStrideType>(1);
  constructAndTestMapIteratorWithStride<RunTimeStrideType>(2);
  constructAndTestMapIteratorWithStride<RunTimeStrideType>(3);

  constructAndTestMapIteratorWithStride<CompileTimeStrideType<1>>(1);
  constructAndTestMapIteratorWithStride<CompileTimeStrideType<2>>(2);
  constructAndTestMapIteratorWithStride<CompileTimeStrideType<3>>(3);

  SLIC_INFO("Done");
}

//----------------------------------------------------------------------
namespace testing
{
//------------------------------------------------------------------------------
// Define some mappings between execution space and allocator.
//  - Host/OpenMP -> Umpire host allocator/default
//  - CUDA/HIP -> Umpire device/unified allocator
//------------------------------------------------------------------------------
template <typename ExecSpace>
struct ExecTraits
{
  constexpr static bool OnDevice = false;
  static int getAllocatorId()
  {
#ifdef AXOM_USE_UMPIRE
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
#else
    return axom::getDefaultAllocatorID();
#endif
  }

  static int getUnifiedAllocatorId()
  {
#ifdef AXOM_USE_UMPIRE
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
#else
    return axom::getDefaultAllocatorID();
#endif
  }
};

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
template <int BLK_SZ>
struct ExecTraits<axom::CUDA_EXEC<BLK_SZ>>
{
  constexpr static bool OnDevice = true;

  static int getAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
  }

  static int getUnifiedAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
  }
};
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
template <int BLK_SZ>
struct ExecTraits<axom::HIP_EXEC<BLK_SZ>>
{
  constexpr static bool OnDevice = true;

  static int getAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
  }

  static int getUnifiedAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
  }
};
#endif

//------------------------------------------------------------------------------
//  This test harness defines some types that are useful for the tests below
//------------------------------------------------------------------------------
template <typename ExecutionSpace>
class slam_map_templated : public ::testing::Test
{
public:
  using ExecSpace = ExecutionSpace;
  using ConcreteSetType = typename slam::RangeSet<SetPosition, SetElement>::ConcreteSet;

  using RealData = axom::Array<double>;
  using IndirectionPolicy = slam::policies::ArrayViewIndirection<SetPosition, double>;
  using StridePolicy = slam::policies::RuntimeStride<axom::IndexType>;
  using InterfacePolicy = slam::policies::ConcreteInterface;

  using RealMap =
    slam::Map<double, ConcreteSetType, IndirectionPolicy, StridePolicy, InterfacePolicy>;

  template <int Dims>
  using MDStridePolicy = slam::policies::MultiDimStride<axom::IndexType, Dims>;

  template <int Dims>
  using MultiDimMap =
    slam::Map<double, ConcreteSetType, IndirectionPolicy, MDStridePolicy<Dims>, InterfacePolicy>;

  slam_map_templated()
    : m_allocatorId(ExecTraits<ExecSpace>::getAllocatorId())
    , m_unifiedAllocatorId(ExecTraits<ExecSpace>::getUnifiedAllocatorId())
  { }

  void initializeWithStride(int stride)
  {
    // Create associated set.
    m_set = ConcreteSetType(MAX_SET_SIZE);

    SLIC_INFO("\nCreating set of size " << m_set.size());

    EXPECT_EQ(m_set.size(), MAX_SET_SIZE);
    EXPECT_TRUE(m_set.isValid());

    // Create array of elements to back the map.
    m_allocatorId = ExecTraits<ExecSpace>::getAllocatorId();
    axom::IndexType backingSize = m_set.size() * stride;

    m_realBacking = RealData(backingSize, backingSize, m_allocatorId);
  }

  template <int Dims>
  void initializeWithMultiDimStride(axom::StackArray<axom::IndexType, Dims> shape)
  {
    // Create associated set.
    m_set = ConcreteSetType(MAX_SET_SIZE);

    SLIC_INFO("\nCreating set of size " << m_set.size());

    EXPECT_EQ(m_set.size(), MAX_SET_SIZE);
    EXPECT_TRUE(m_set.isValid());

    int stride = 1;
    for(int dim = 0; dim < Dims; dim++)
    {
      stride *= shape[dim];
    }

    // Create array of elements to back the map.
    axom::IndexType backingSize = m_set.size() * stride;

    m_realBacking = RealData(backingSize, backingSize, m_unifiedAllocatorId);
  }

protected:
  int m_allocatorId;
  int m_unifiedAllocatorId;
  ConcreteSetType m_set;
  RealData m_realBacking;
};

using MyTypes = ::testing::Types<
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  axom::OMP_EXEC,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  axom::CUDA_EXEC<256>,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  axom::HIP_EXEC<256>,
#endif
  axom::SEQ_EXEC>;

TYPED_TEST_SUITE(slam_map_templated, MyTypes);

//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_map_templated, constructAndTestStride1)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using MapType = typename TestFixture::RealMap;

  const int stride = 1;

  this->initializeWithStride(stride);

  SLIC_INFO("\nCreating double map with stride 1 on the set ");
  MapType m(this->m_set, this->m_realBacking.view(), stride);

  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("\nSetting the elements.");
  const double multFac = 100.0001;
  axom::for_all<ExecSpace>(this->m_set.size(), AXOM_LAMBDA(int index) { m(index) = index * multFac; });

  SLIC_INFO("\nChecking the elements.");
  int totalSize = this->m_set.size() * stride;
  axom::Array<int> isValid(totalSize, totalSize, this->m_unifiedAllocatorId);
  const auto isValid_view = isValid.view();

  axom::for_all<ExecSpace>(
    this->m_set.size(),
    AXOM_LAMBDA(int index) {
      bool entryValid = (m(index) == index * multFac);
      isValid_view[index] = entryValid;
    });

  for(int validEntry : isValid)
  {
    EXPECT_TRUE(validEntry);
  }
}

//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_map_templated, constructAndTestStride3)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using MapType = typename TestFixture::RealMap;

  const int stride = 3;

  this->initializeWithStride(stride);

  SLIC_INFO("\nCreating double map with stride 3 on the set ");
  MapType m(this->m_set, this->m_realBacking.view(), stride);

  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("\nSetting the elements.");
  const double multFac = 100.0001;
  const double multFac2 = 1.010;
  axom::for_all<ExecSpace>(
    this->m_set.size(),
    AXOM_LAMBDA(int index) {
      for(int comp = 0; comp < stride; comp++)
      {
        m(index, comp) = index * multFac + comp * multFac2;
      }
    });

  SLIC_INFO("\nChecking the elements.");
  int totalSize = this->m_set.size() * stride;
  axom::Array<int> isValid(totalSize, totalSize, this->m_unifiedAllocatorId);
  const auto isValid_view = isValid.data();

  axom::for_all<ExecSpace>(
    this->m_set.size(),
    AXOM_LAMBDA(int index) {
      for(int comp = 0; comp < stride; comp++)
      {
        bool entryValid = (m(index, comp) == index * multFac + comp * multFac2);
        isValid_view[index * stride + comp] = entryValid;
      }
    });

  for(int validEntry : isValid)
  {
    EXPECT_TRUE(validEntry);
  }
}

//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_map_templated, constructAndTest2DStride)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using MapType = typename TestFixture::template MultiDimMap<2>;

  const axom::StackArray<axom::IndexType, 2> shape = {3, 5};
  const axom::StackArray<axom::IndexType, 2> strides = {5, 1};
  int stride = 3 * 5;
  this->initializeWithMultiDimStride(shape);

  SLIC_INFO("\nCreating double map with shape (3, 5) on the set ");
  MapType m(this->m_set, this->m_realBacking.view(), shape);

  EXPECT_EQ(m.stride(), stride);
  EXPECT_EQ(m.shape(), shape);

  SLIC_INFO("\nSetting the elements.");
  const double multFac = 100.0001;
  const double multFac2 = 1.00100;
  const double multFac3 = 0.10010;
  axom::for_all<ExecSpace>(
    this->m_set.size(),
    AXOM_LAMBDA(int index) {
      for(int i = 0; i < shape[0]; i++)
      {
        for(int j = 0; j < shape[1]; j++)
        {
          m(index, i, j) = index * multFac + i * multFac2 + j * multFac3;
        }
      }
    });

  SLIC_INFO("\nChecking the elements.");

  for(int setIdx = 0; setIdx < this->m_set.size(); setIdx++)
  {
    for(int i = 0; i < shape[0]; i++)
    {
      for(int j = 0; j < shape[1]; j++)
      {
        double expectedValue = setIdx * multFac + i * multFac2 + j * multFac3;
        EXPECT_DOUBLE_EQ(m(setIdx, i, j), expectedValue);
        EXPECT_DOUBLE_EQ(m.value(setIdx, i, j), expectedValue);

        int flatIndex = i * strides[0] + j * strides[1];
        EXPECT_DOUBLE_EQ(m[setIdx * stride + flatIndex], expectedValue);
      }
    }
  }

  SLIC_INFO("\nChecking iteration through range iterator.");
  for(auto it = m.set_begin(); it != m.set_end(); ++it)
  {
    int setIdx = it.flatIndex();
    EXPECT_EQ(m.index(setIdx), it.index());
    EXPECT_EQ(it->shape(), shape);
    EXPECT_EQ(it->size(), it.numComp());
    EXPECT_EQ(m.set_begin() + setIdx, it);
    for(int i = 0; i < shape[0]; i++)
    {
      for(int j = 0; j < shape[1]; j++)
      {
        double expectedValue = setIdx * multFac + i * multFac2 + j * multFac3;
        EXPECT_DOUBLE_EQ(expectedValue, (*it)(i, j));
        EXPECT_DOUBLE_EQ(expectedValue, it(i, j));
        EXPECT_DOUBLE_EQ(expectedValue, it.value(i, j));
      }
    }
  }
}
//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_map_templated, constructAndTest3DStride)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using MapType = typename TestFixture::template MultiDimMap<3>;

  const axom::StackArray<axom::IndexType, 3> shape = {2, 3, 4};
  const axom::StackArray<axom::IndexType, 3> strides = {12, 4, 1};
  int stride = 2 * 3 * 4;
  this->initializeWithMultiDimStride(shape);

  SLIC_INFO("\nCreating double map with shape (2, 3, 4) on the set ");
  MapType m(this->m_set, this->m_realBacking.view(), shape);

  EXPECT_EQ(m.stride(), stride);
  EXPECT_EQ(m.shape(), shape);

  SLIC_INFO("\nSetting the elements.");
  const double multFac = 100.0001;
  const double multFac2 = 1.00100;
  const double multFac3 = 0.10010;
  const double multFac4 = 0.01001;
  axom::for_all<ExecSpace>(
    this->m_set.size(),
    AXOM_LAMBDA(int index) {
      for(int i = 0; i < shape[0]; i++)
      {
        for(int j = 0; j < shape[1]; j++)
        {
          for(int k = 0; k < shape[2]; k++)
          {
            m(index, i, j, k) = index * multFac + i * multFac2 + j * multFac3 + k * multFac4;
          }
        }
      }
    });

  SLIC_INFO("\nChecking the elements.");

  for(int setIdx = 0; setIdx < this->m_set.size(); setIdx++)
  {
    for(int i = 0; i < shape[0]; i++)
    {
      for(int j = 0; j < shape[1]; j++)
      {
        for(int k = 0; k < shape[2]; k++)
        {
          double expectedValue = setIdx * multFac + i * multFac2 + j * multFac3 + k * multFac4;
          EXPECT_DOUBLE_EQ(m(setIdx, i, j, k), expectedValue);
          EXPECT_DOUBLE_EQ(m.value(setIdx, i, j, k), expectedValue);

          int flatIndex = i * strides[0] + j * strides[1] + k * strides[2];
          EXPECT_DOUBLE_EQ(m[setIdx * stride + flatIndex], expectedValue);
        }
      }
    }
  }

  SLIC_INFO("\nChecking iteration through range iterator.");
  for(auto it = m.set_begin(); it != m.set_end(); ++it)
  {
    int setIdx = it.flatIndex();
    EXPECT_EQ(m.index(setIdx), it.index());
    EXPECT_EQ(it->shape(), shape);
    EXPECT_EQ(it->size(), it.numComp());
    EXPECT_EQ(m.set_begin() + setIdx, it);
    for(int i = 0; i < shape[0]; i++)
    {
      for(int j = 0; j < shape[1]; j++)
      {
        for(int k = 0; k < shape[2]; k++)
        {
          double expectedValue = setIdx * multFac + i * multFac2 + j * multFac3 + k * multFac4;
          EXPECT_DOUBLE_EQ(expectedValue, (*it)(i, j, k));
          EXPECT_DOUBLE_EQ(expectedValue, it(i, j, k));
          EXPECT_DOUBLE_EQ(expectedValue, it.value(i, j, k));
        }
      }
    }
  }
}

}  // namespace testing

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
