// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_map_BivariateMap.cpp
 *
 * \brief Unit tests for Slam's Bivariate Map
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
namespace traits = axom::slam::traits;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;
using SetType = slam::RangeSet<SetPosition, SetElement>;

using StrideOneType = policies::StrideOne<SetPosition>;

template <unsigned int S>
using CompileTimeStrideType = policies::CompileTimeStride<SetPosition, S>;

using RuntimeStrideType = policies::RuntimeStride<SetPosition>;

template <typename T>
using STLIndirection = policies::STLVectorIndirection<SetPosition, T>;
using VariableCardinality = policies::VariableCardinality<SetPosition, STLIndirection<SetElement>>;

using RelationType =
  slam::StaticRelation<SetPosition, SetElement, VariableCardinality, STLIndirection<SetElement>, SetType, SetType>;

using BivariateSetType = slam::BivariateSet<SetType, SetType>;
using ProductSetType = slam::ProductSet<SetType, SetType>;
using RelationSetType = slam::RelationSet<RelationType>;

template <typename T, typename B, typename I, typename S>
using BivariateMapType = slam::BivariateMap<T, B, I, S>;

constexpr SetPosition MAX_SET_SIZE1 = 10;
constexpr SetPosition MAX_SET_SIZE2 = 15;

constexpr double multFac3 = 0000.1;
constexpr double multFac1 = 1000.0;
constexpr double multFac2 = 0010.0;

}  // end anonymous namespace

TEST(slam_bivariate_map, construct_empty_map)
{
  slam::BivariateMap<int> m;

  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(m.totalSize(), 0);
  EXPECT_EQ(m.firstSetSize(), 0);
  EXPECT_EQ(m.secondSetSize(), 0);
}

template <typename T>
AXOM_HOST_DEVICE inline T getVal(SetPosition idx1, SetPosition idx2, SetPosition idx3 = 0)
{
  return static_cast<T>(idx1 * multFac1 + idx2 * multFac2 + idx3 * multFac3);
}

template <typename T, typename B, typename I, typename S>
void constructAndTestCartesianMap(int stride)
{
  SLIC_INFO("Testing BivariateMap on ProductSet with stride " << stride);

  SLIC_INFO("Creating set");
  using BMapType = BivariateMapType<T, B, I, S>;
  using SubMapType = typename BMapType::SubMapType;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string() << " map on the set ");

  BMapType m(&s, static_cast<T>(0), stride);

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements in the map.");

  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
    {
      for(auto i = 0; i < stride; i++)
      {
        T* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<T>(idx1, idx2, i);
      }
    }
  }

  SLIC_INFO("Checking the elements with findValue().");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
    {
      for(auto i = 0; i < stride; i++)
      {
        T* ptr = m.findValue(idx1, idx2, i);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(*ptr, getVal<T>(idx1, idx2, i));
      }
    }
  }

  SLIC_INFO("Checking the elements with SubMap.");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    SubMapType sm = m(idx1);
    for(auto idx2 = 0; idx2 < sm.size(); ++idx2)
    {
      for(auto i = 0; i < stride; i++)
      {
        T v = sm.value(idx2, i);
        EXPECT_EQ(v, getVal<T>(idx1, idx2, i));
        EXPECT_EQ(sm.index(idx2), idx2);
      }
    }
  }

  EXPECT_TRUE(m.isValid());
}

TEST(slam_bivariate_map, construct_int_map)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<int>;

  constructAndTestCartesianMap<int, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestCartesianMap<int, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestCartesianMap<int, BSet, IndPol, RuntimeStrideType>(3);

  constructAndTestCartesianMap<int, BSet, IndPol, StrideOneType>(1);

  constructAndTestCartesianMap<int, BSet, IndPol, CompileTimeStrideType<1>>(1);
  constructAndTestCartesianMap<int, BSet, IndPol, CompileTimeStrideType<2>>(2);
  constructAndTestCartesianMap<int, BSet, IndPol, CompileTimeStrideType<3>>(3);
}

TEST(slam_bivariate_map, construct_double_map)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<double>;

  constructAndTestCartesianMap<double, BSet, IndPol, StrideOneType>(1);

  constructAndTestCartesianMap<double, BSet, IndPol, CompileTimeStrideType<1>>(1);
  constructAndTestCartesianMap<double, BSet, IndPol, CompileTimeStrideType<2>>(2);
  constructAndTestCartesianMap<double, BSet, IndPol, CompileTimeStrideType<3>>(3);

  constructAndTestCartesianMap<double, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestCartesianMap<double, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestCartesianMap<double, BSet, IndPol, RuntimeStrideType>(3);
}

template <typename T, typename B, typename I, typename S>
void constructAndTestRelationSetMap(int stride)
{
  SLIC_INFO("Testing BivariateMap on RelationSet with stride " << stride);

  SLIC_INFO("Creating set");
  using MapType = BivariateMapType<T, B, I, S>;
  using SubMapType = typename MapType::SubMapType;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);

  RelationType rel(&s1, &s2);

  std::vector<SetPosition> begin_vec(MAX_SET_SIZE1 + 1, 0);
  std::vector<SetPosition> indice_vec;

  auto curIdx = SetPosition();

  for(auto i = 0; i < MAX_SET_SIZE1; ++i)
  {
    begin_vec[i] = curIdx;
    if(MAX_SET_SIZE1 / 4 <= i && i <= MAX_SET_SIZE1 / 4 * 3)
    {
      for(auto j = MAX_SET_SIZE2 / 4; j < MAX_SET_SIZE2 / 4 * 3; ++j)
      {
        indice_vec.push_back(j);
        ++curIdx;
      }
    }
  }
  begin_vec[MAX_SET_SIZE1] = curIdx;

  rel.bindBeginOffsets(MAX_SET_SIZE1, &begin_vec);
  rel.bindIndices(indice_vec.size(), &indice_vec);
  RelationSetType s(&rel);

  const SetPosition indice_size = indice_vec.size();
  EXPECT_EQ(indice_size, s.totalSize());
  EXPECT_TRUE(s.isValid(true));

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string() << " map on the set ");

  MapType m(&s, (T)0, stride);

  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(indice_size, m.totalSize());
  EXPECT_EQ(rel.fromSetSize(), m.firstSetSize());
  EXPECT_EQ(rel.toSetSize(), m.secondSetSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements in the map.");

  for(auto idx1 = 0; idx1 < rel.fromSetSize(); idx1++)
  {
    auto relsubset = rel[idx1];
    for(auto si = 0; si < relsubset.size(); ++si)
    {
      auto idx2 = relsubset[si];
      for(auto i = 0; i < stride; i++)
      {
        T* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<T>(idx1, idx2, i);
      }
    }
  }

  SLIC_INFO("Checking the elements with findValue().");
  for(auto idx1 = 0; idx1 < rel.fromSetSize(); idx1++)
  {
    auto relsubset = rel[idx1];
    auto rel_idx = 0;
    for(auto idx2 = 0; idx2 < rel.toSetSize(); ++idx2)
    {
      bool isInRel = relsubset.size() > rel_idx && relsubset[rel_idx] == idx2;
      for(auto i = 0; i < stride; i++)
      {
        T* ptr = m.findValue(idx1, idx2, i);
        if(isInRel)
        {
          EXPECT_NE(ptr, nullptr);
          EXPECT_EQ(*ptr, getVal<T>(idx1, idx2, i));
        }
        else
        {
          EXPECT_EQ(ptr, nullptr);
        }
      }
      if(isInRel)
      {
        rel_idx++;
      }
    }
  }

  SLIC_INFO("Checking the elements with SubMap.");
  for(auto idx1 = 0; idx1 < rel.fromSetSize(); idx1++)
  {
    auto relsubset = rel[idx1];
    SubMapType sm = m(idx1);
    for(auto idx2 = 0; idx2 < sm.size(); ++idx2)
    {
      ASSERT_EQ(relsubset[idx2], sm.index(idx2));
      for(auto i = 0; i < stride; i++)
      {
        T v = sm.value(idx2, i);
        EXPECT_EQ(v, getVal<T>(idx1, sm.index(idx2), i));
      }
    }
  }

  EXPECT_TRUE(m.isValid());
}

TEST(slam_bivariate_map, construct_int_relset_map)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<int>;

  constructAndTestRelationSetMap<int, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestRelationSetMap<int, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestRelationSetMap<int, BSet, IndPol, RuntimeStrideType>(3);

  constructAndTestRelationSetMap<int, BSet, IndPol, StrideOneType>(1);

  constructAndTestRelationSetMap<int, BSet, IndPol, CompileTimeStrideType<1>>(1);
  constructAndTestRelationSetMap<int, BSet, IndPol, CompileTimeStrideType<2>>(2);
  constructAndTestRelationSetMap<int, BSet, IndPol, CompileTimeStrideType<3>>(3);
}

TEST(slam_bivariate_map, construct_double_relset_map)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<double>;

  constructAndTestRelationSetMap<double, BSet, IndPol, StrideOneType>(1);

  constructAndTestRelationSetMap<double, BSet, IndPol, CompileTimeStrideType<1>>(1);
  constructAndTestRelationSetMap<double, BSet, IndPol, CompileTimeStrideType<2>>(2);
  constructAndTestRelationSetMap<double, BSet, IndPol, CompileTimeStrideType<3>>(3);

  constructAndTestRelationSetMap<double, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestRelationSetMap<double, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestRelationSetMap<double, BSet, IndPol, RuntimeStrideType>(3);
}

template <typename BSet, typename StridePolicy>
void constructAndTestBivariateMapIterator(int stride)
{
  SLIC_INFO("Creating set");
  using DataType = double;
  using IndPol = STLIndirection<DataType>;
  using MapType = BivariateMapType<DataType, BSet, IndPol, StridePolicy>;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<DataType>::to_string() << " map on the set ");
  MapType m(&s, 0.0, stride);
  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements in the map.");
  //currently can't set value using iterator
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
    {
      for(auto i = 0; i < stride; i++)
      {
        DataType* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<DataType>(idx1, idx2, i);
      }
    }
  }

  SLIC_INFO("Checking the elements with SubMap flat iterator.");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    int idx2 = 0;
    int compIdx = 0;
    for(auto iter = m.begin(idx1); iter != m.end(idx1); ++iter)
    {
      EXPECT_EQ(*iter, getVal<DataType>(idx1, idx2, compIdx));
      compIdx++;
      if(compIdx == m.numComp())
      {
        compIdx = 0;
        idx2++;
      }
    }
  }
  SLIC_INFO("Checking the elements with BivariateMap flat iterator.");
  {
    auto iter = m.begin();

    for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    {
      for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      {
        for(auto i = 0; i < stride; i++)
        {
          // Check validity of indexing
          EXPECT_EQ(iter.firstIndex(), idx1);
          EXPECT_EQ(iter.secondIndex(), idx2);
          EXPECT_EQ(iter.compIndex(), i);
          EXPECT_NE(iter, m.end());

          // Check equality of values
          DataType val = getVal<DataType>(idx1, idx2, i);
          EXPECT_EQ(val, *iter);

          iter++;
          // flat_idx++;
        }
      }
    }
  }

  SLIC_INFO("Checking the elements with BivariateMap range iterator.");
  {
    auto iter = m.set_begin();
    for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    {
      for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      {
        EXPECT_EQ((*iter).size(), stride);
        EXPECT_EQ(iter.firstIndex(), idx1);
        EXPECT_EQ(iter.secondIndex(), idx2);
        for(auto i = 0; i < stride; i++)
        {
          DataType val = getVal<DataType>(idx1, idx2, i);
          EXPECT_EQ(iter.value(i), val);
          EXPECT_EQ(iter(i), val);
          // Below disabled because we can't test these with just a forward access iterator
          //EXPECT_EQ((begin_iter + flat_idx).value(i), val);
          //EXPECT_EQ((end_iter - (m.totalSize() - flat_idx)).value(i), val);
          //EXPECT_EQ((inval_iter - (m.totalSize() - flat_idx + 1)).value(i), val);
        }
        iter++;
        // flat_idx++;
      }
    }
  }

  EXPECT_TRUE(m.isValid());
}

TEST(slam_bivariate_map, iterate)
{
  using BSet = BivariateSetType;

  constructAndTestBivariateMapIterator<BSet, RuntimeStrideType>(1);
  constructAndTestBivariateMapIterator<BSet, RuntimeStrideType>(2);
  constructAndTestBivariateMapIterator<BSet, RuntimeStrideType>(3);

  constructAndTestBivariateMapIterator<BSet, CompileTimeStrideType<1>>(1);
  constructAndTestBivariateMapIterator<BSet, CompileTimeStrideType<2>>(2);
  constructAndTestBivariateMapIterator<BSet, CompileTimeStrideType<3>>(3);

  constructAndTestBivariateMapIterator<BSet, StrideOneType>(1);
}

template <typename T, typename B, typename I, typename S>
void testScopedCopyBehavior(int stride)
{
  using BMapType = BivariateMapType<T, B, I, S>;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string() << " map on the set ");

  BMapType m;
  {
    BMapType m_inner(&s, static_cast<T>(0), stride);

    EXPECT_TRUE(m_inner.isValid());
    EXPECT_EQ(s.size(), m_inner.totalSize());
    EXPECT_EQ(m_inner.stride(), stride);

    SLIC_INFO("Setting the elements in the map.");

    for(auto idx1 = 0; idx1 < m_inner.firstSetSize(); ++idx1)
    {
      for(auto idx2 = 0; idx2 < m_inner.secondSetSize(); ++idx2)
      {
        for(auto i = 0; i < stride; i++)
        {
          T* valPtr = m_inner.findValue(idx1, idx2, i);
          EXPECT_NE(valPtr, nullptr);
          *valPtr = getVal<T>(idx1, idx2, i);
        }
      }
    }

    m = m_inner;
  }

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Checking the elements with findValue().");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
    {
      for(auto i = 0; i < stride; i++)
      {
        T* ptr = m.findValue(idx1, idx2, i);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(*ptr, getVal<T>(idx1, idx2, i));
      }
    }
  }
}

TEST(slam_bivariate_map, testScopedMapBehavior)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<double>;

  testScopedCopyBehavior<double, BSet, IndPol, StrideOneType>(1);

  testScopedCopyBehavior<double, BSet, IndPol, CompileTimeStrideType<1>>(1);
  testScopedCopyBehavior<double, BSet, IndPol, CompileTimeStrideType<2>>(2);
  testScopedCopyBehavior<double, BSet, IndPol, CompileTimeStrideType<3>>(3);

  testScopedCopyBehavior<double, BSet, IndPol, RuntimeStrideType>(1);
  testScopedCopyBehavior<double, BSet, IndPol, RuntimeStrideType>(2);
  testScopedCopyBehavior<double, BSet, IndPol, RuntimeStrideType>(3);
}

TEST(slam_bivariate_map, traits)
{
  EXPECT_TRUE(traits::indices_use_indirection<RelationSetType>::value);
  EXPECT_TRUE(traits::indices_use_indirection<BivariateSetType>::value);
  EXPECT_FALSE(traits::indices_use_indirection<ProductSetType>::value);
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
class slam_bivariate_map_templated : public ::testing::Test
{
public:
  using ExecSpace = ExecutionSpace;
  using ConcreteSetType = typename slam::RangeSet<SetPosition, SetElement>::ConcreteSet;

  // StaticRelation template types
  using ElemIndirection = slam::policies::ArrayViewIndirection<SetPosition, SetElement>;
  using VariableCardinality = policies::VariableCardinality<SetPosition, ElemIndirection>;
  using RelationType =
    slam::StaticRelation<SetPosition, SetElement, VariableCardinality, ElemIndirection, ConcreteSetType, ConcreteSetType>;

  // BivariateSet concrete types -- ProductSet and RelationSet
  using ProductSetType = typename slam::ProductSet<ConcreteSetType, ConcreteSetType>::ConcreteSet;
  using RelationSetType = typename slam::RelationSet<RelationType>::ConcreteSet;

  // BivariateMap template types
  using RealData = axom::Array<double>;
  using IndirectionPolicy = slam::policies::ArrayViewIndirection<SetPosition, double>;
  using StridePolicy = slam::policies::RuntimeStride<int>;
  using InterfacePolicy = slam::policies::ConcreteInterface;
  using RelationMapType =
    slam::BivariateMap<double, RelationSetType, IndirectionPolicy, StridePolicy, InterfacePolicy>;
  using CartesianMapType =
    slam::BivariateMap<double, ProductSetType, IndirectionPolicy, StridePolicy, InterfacePolicy>;

  template <int Dims>
  using NDStridePolicy = slam::policies::MultiDimStride<int, Dims>;
  template <int Dims>
  using CartesianNDMapType =
    slam::BivariateMap<double, ProductSetType, IndirectionPolicy, NDStridePolicy<Dims>, InterfacePolicy>;
  template <int Dims>
  using RelationNDMapType =
    slam::BivariateMap<double, RelationSetType, IndirectionPolicy, NDStridePolicy<Dims>, InterfacePolicy>;

  slam_bivariate_map_templated()
    : m_allocatorId(ExecTraits<ExecSpace>::getAllocatorId())
    , m_unifiedAllocatorId(ExecTraits<ExecSpace>::getUnifiedAllocatorId())
  { }

  void initializeAndTestCartesianMap(int stride);

  void initializeAndTestCartesianMap(axom::StackArray<int, 3> shape);

  void initializeAndTestRelationMap(int stride);

  void initializeAndTestRelationMap(axom::StackArray<int, 3> shape);

protected:
  int m_allocatorId;
  int m_unifiedAllocatorId;
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

TYPED_TEST_SUITE(slam_bivariate_map_templated, MyTypes);

//----------------------------------------------------------------------
template <typename ExecutionSpace>
void slam_bivariate_map_templated<ExecutionSpace>::initializeAndTestCartesianMap(int stride)
{
  using MapType = CartesianMapType;

  // Create associated sets.
  axom::Array<ConcreteSetType> sets(2, 2, m_unifiedAllocatorId);
  sets[0] = ConcreteSetType(MAX_SET_SIZE1);
  sets[1] = ConcreteSetType(MAX_SET_SIZE2);

  SLIC_INFO("Creating product set with size (" << MAX_SET_SIZE1 << ", " << MAX_SET_SIZE2 << ")");
  ProductSetType prodSet(&sets[0], &sets[1]);
  EXPECT_EQ(prodSet.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(prodSet.isValid());

  // Create array of elements to back the map.
  m_allocatorId = ExecTraits<ExecSpace>::getAllocatorId();
  axom::IndexType backingSize = prodSet.size() * stride;

  RealData realBacking(backingSize, backingSize, m_allocatorId);

  SLIC_INFO("\nCreating double map with stride 1 on the set ");
  const MapType m(prodSet, realBacking.view(), stride);

  EXPECT_EQ(m.stride(), stride);
  SLIC_INFO("\nSetting the elements.");
  axom::for_all<ExecSpace>(
    m.firstSetSize(),
    AXOM_LAMBDA(int idx1) {
      for(auto idx2 = 0; idx2 < m.secondSetSize(); idx2++)
      {
        for(auto comp = 0; comp < stride; comp++)
        {
          m(idx1, idx2, comp) = getVal<double>(idx1, idx2, comp);
        }
      }
    });

  int totalSize = prodSet.size() * stride;
  axom::Array<int> isValid(totalSize, totalSize, m_unifiedAllocatorId);
  const auto isValid_view = isValid.data();

  SLIC_INFO("\nChecking the elements with findValue().");
  axom::for_all<ExecSpace>(
    m.firstSetSize(),
    AXOM_LAMBDA(int idx1) {
      for(auto idx2 = 0; idx2 < m.secondSetSize(); idx2++)
      {
        for(auto comp = 0; comp < stride; comp++)
        {
          int flatIdx = idx1 * m.secondSetSize() * stride;
          flatIdx += idx2 * stride;
          flatIdx += comp;

          double* ptr = m.findValue(idx1, idx2, comp);
          bool hasValue = (ptr != nullptr);
          hasValue = hasValue && (*ptr == getVal<double>(idx1, idx2, comp));
          isValid_view[flatIdx] = hasValue;
        }
      }
    });

  for(int validEntry : isValid)
  {
    EXPECT_TRUE(validEntry);
  }

  SLIC_INFO("\nChecking the elements with SubMap.");
  axom::for_all<ExecSpace>(
    m.firstSetSize(),
    AXOM_LAMBDA(int idx1) {
      auto submap = m(idx1);
      for(auto idx2 = 0; idx2 < m.secondSetSize(); idx2++)
      {
        for(auto comp = 0; comp < stride; comp++)
        {
          int flatIdx = idx1 * m.secondSetSize() * stride;
          flatIdx += idx2 * stride;
          flatIdx += comp;

          double value = submap(idx2, comp);
          bool hasValue = (value == getVal<double>(idx1, idx2, comp));
          isValid_view[flatIdx] = hasValue;
        }
      }
    });

  for(int validEntry : isValid)
  {
    EXPECT_TRUE(validEntry);
  }
}

//----------------------------------------------------------------------
template <typename ExecutionSpace>
void slam_bivariate_map_templated<ExecutionSpace>::initializeAndTestCartesianMap(
  axom::StackArray<int, 3> shape)
{
  using MapType = CartesianNDMapType<3>;

  // Create associated sets.
  axom::Array<ConcreteSetType> sets(2, 2, m_unifiedAllocatorId);
  sets[0] = ConcreteSetType(MAX_SET_SIZE1);
  sets[1] = ConcreteSetType(MAX_SET_SIZE2);

  SLIC_INFO("Creating product set with size (" << MAX_SET_SIZE1 << ", " << MAX_SET_SIZE2 << ")");
  ProductSetType prodSet(&sets[0], &sets[1]);
  EXPECT_EQ(prodSet.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(prodSet.isValid());

  int flatStride = shape[0] * shape[1] * shape[2];
  int strides[3] = {shape[1] * shape[2], shape[2], 1};

  // Create array of elements to back the map.
  m_allocatorId = ExecTraits<ExecSpace>::getAllocatorId();
  axom::IndexType backingSize = prodSet.size() * flatStride;

  RealData realBacking(backingSize, backingSize, m_unifiedAllocatorId);

  SLIC_INFO(axom::fmt::format("\nCreating double map with shape ({}) on the ProductSet",
                              axom::fmt::join(shape, ", ")));
  const MapType m(prodSet, realBacking.view(), shape);

  EXPECT_EQ(m.stride(), flatStride);
  SLIC_INFO("\nSetting the elements.");
  axom::for_all<ExecSpace>(
    m.firstSetSize(),
    AXOM_LAMBDA(int idx1) {
      for(auto idx2 = 0; idx2 < m.secondSetSize(); idx2++)
      {
        for(int i = 0; i < shape[0]; i++)
        {
          for(int j = 0; j < shape[1]; j++)
          {
            for(int k = 0; k < shape[2]; k++)
            {
              int flatCompIdx = i * strides[0] + j * strides[1] + k * strides[2];
              m(idx1, idx2, i, j, k) = getVal<double>(idx1, idx2, flatCompIdx);
            }
          }
        }
      }
    });

  SLIC_INFO("\nChecking the elements with findValue().");
  for(int idx1 = 0; idx1 < m.firstSetSize(); idx1++)
  {
    for(auto idx2 = 0; idx2 < m.secondSetSize(); idx2++)
    {
      int bsetIndex = idx1 * m.secondSetSize() + idx2;
      for(int i = 0; i < shape[0]; i++)
      {
        for(int j = 0; j < shape[1]; j++)
        {
          for(int k = 0; k < shape[2]; k++)
          {
            int flatCompIdx = i * strides[0] + j * strides[1] + k * strides[2];
            int flatIdx = bsetIndex * flatStride;
            flatIdx += flatCompIdx;

            double* ptr = m.findValue(idx1, idx2, i, j, k);
            EXPECT_NE(ptr, nullptr);
            EXPECT_DOUBLE_EQ(*ptr, getVal<double>(idx1, idx2, flatCompIdx));
            // Test other access methods:
            EXPECT_DOUBLE_EQ(*ptr, m.flatValue(bsetIndex, i, j, k));
            EXPECT_DOUBLE_EQ(*ptr, m(idx1, idx2, i, j, k));
            EXPECT_DOUBLE_EQ(*ptr, m[flatIdx]);
          }
        }
      }
    }
  }

  SLIC_INFO("\nChecking the elements with BivariateMap range iterator.");
  for(auto it = m.set_begin(); it != m.set_end(); ++it)
  {
    int idx1 = it.firstIndex();
    int idx2 = it.secondIndex();
    int flatIdx = it.flatIndex();

    EXPECT_EQ(idx1, m.set()->flatToFirstIndex(flatIdx));
    EXPECT_EQ(idx2, m.set()->flatToSecondIndex(flatIdx));
    EXPECT_EQ(it.numComp(), m.stride());

    for(int i = 0; i < shape[0]; i++)
    {
      for(int j = 0; j < shape[1]; j++)
      {
        for(int k = 0; k < shape[2]; k++)
        {
          int flatCompIdx = i * strides[0] + j * strides[1] + k * strides[2];
          double expected_value = getVal<double>(idx1, idx2, flatCompIdx);

          EXPECT_DOUBLE_EQ(expected_value, (*it)(i, j, k));
          EXPECT_DOUBLE_EQ(expected_value, it(i, j, k));
          EXPECT_DOUBLE_EQ(expected_value, it.value(i, j, k));
        }
      }
    }
  }
}

//----------------------------------------------------------------------
template <typename ExecutionSpace>
void slam_bivariate_map_templated<ExecutionSpace>::initializeAndTestRelationMap(int stride)
{
  using MapType = RelationMapType;

  // Create associated sets.
  axom::Array<ConcreteSetType> sets(2, 2, m_unifiedAllocatorId);
  sets[0] = ConcreteSetType(MAX_SET_SIZE1);
  sets[1] = ConcreteSetType(MAX_SET_SIZE2);

  // Create a relation on the two sets.
  SLIC_INFO("Creating static relation between two sets.");
  axom::Array<RelationType> rel(1, 1, m_unifiedAllocatorId);
  rel[0] = RelationType(&sets[0], &sets[1]);
  axom::Array<SetPosition> begin_vec(MAX_SET_SIZE1 + 1, MAX_SET_SIZE1 + 1, m_unifiedAllocatorId);
  axom::Array<SetPosition> index_vec(0, 0, m_unifiedAllocatorId);

  SetPosition curIdx = 0;

  for(auto i = 0; i < MAX_SET_SIZE1; ++i)
  {
    begin_vec[i] = curIdx;
    if(MAX_SET_SIZE1 / 4 <= i && i <= MAX_SET_SIZE1 / 4 * 3)
    {
      for(auto j = MAX_SET_SIZE2 / 4; j < MAX_SET_SIZE2 / 4 * 3; ++j)
      {
        index_vec.push_back(j);
        ++curIdx;
      }
    }
  }
  begin_vec[MAX_SET_SIZE1] = curIdx;

  rel[0].bindBeginOffsets(MAX_SET_SIZE1, begin_vec.view());
  rel[0].bindIndices(index_vec.size(), index_vec.view());

  RelationType* relPtr = &rel[0];

  RelationSetType relSet(&rel[0]);
  EXPECT_EQ(index_vec.size(), relSet.totalSize());
  EXPECT_TRUE(relSet.isValid());

  // Create array of elements to back the map.
  m_allocatorId = ExecTraits<ExecSpace>::getAllocatorId();
  axom::IndexType backingSize = index_vec.size() * stride;

  RealData realBacking(backingSize, backingSize, m_allocatorId);

  SLIC_INFO("\nCreating double map with stride " << stride << " on the RelationSet ");

  MapType m(relSet, realBacking.view(), stride);

  EXPECT_EQ(m.stride(), stride);
  SLIC_INFO("\nSetting the elements.");
  axom::for_all<ExecSpace>(
    m.firstSetSize(),
    AXOM_LAMBDA(int idx1) {
      auto relSubset = (*relPtr)[idx1];
      for(auto slot = 0; slot < relSubset.size(); slot++)
      {
        auto idx2 = relSubset[slot];
        for(auto comp = 0; comp < stride; comp++)
        {
          double* valPtr = m.findValue(idx1, idx2, comp);
#ifndef AXOM_DEVICE_CODE
          EXPECT_NE(valPtr, nullptr);
#endif
          *valPtr = getVal<double>(idx1, idx2, comp);
        }
      }
    });

  SLIC_INFO("\nChecking the elements with findValue().");
  {
    axom::ReduceSum<ExecSpace, int> numIncorrect(0);

    axom::for_all<ExecSpace>(
      m.firstSetSize(),
      AXOM_LAMBDA(int idx1) {
        auto relSubset = (*relPtr)[idx1];
        auto relIndex = 0;
        for(auto idx2 = 0; idx2 < m.secondSetSize(); idx2++)
        {
          bool inRelation = relSubset.size() > relIndex && relSubset[relIndex] == idx2;
          for(auto comp = 0; comp < stride; comp++)
          {
            double* ptr = m.findValue(idx1, idx2, comp);
            if(inRelation)
            {
              numIncorrect += (ptr == nullptr);
              numIncorrect += (*ptr != getVal<double>(idx1, idx2, comp));
            }
            else
            {
              numIncorrect += (ptr != nullptr);
            }
          }
          if(inRelation)
          {
            relIndex++;
          }
        }
      });

    EXPECT_EQ(numIncorrect.get(), 0);
  }
}
//----------------------------------------------------------------------
template <typename ExecutionSpace>
void slam_bivariate_map_templated<ExecutionSpace>::initializeAndTestRelationMap(
  axom::StackArray<int, 3> shape)
{
  using MapType = RelationNDMapType<3>;

  // Compute strides.
  int flatStride = shape[0] * shape[1] * shape[2];
  int strides[3] = {shape[1] * shape[2], shape[2], 1};

  // Create associated sets.
  axom::Array<ConcreteSetType> sets(2, 2, m_unifiedAllocatorId);
  sets[0] = ConcreteSetType(MAX_SET_SIZE1);
  sets[1] = ConcreteSetType(MAX_SET_SIZE2);

  // Create a relation on the two sets.
  SLIC_INFO("Creating static relation between two sets.");
  axom::Array<RelationType> rel(1, 1, m_unifiedAllocatorId);
  rel[0] = RelationType(&sets[0], &sets[1]);
  axom::Array<SetPosition> begin_vec(MAX_SET_SIZE1 + 1, MAX_SET_SIZE1 + 1, m_unifiedAllocatorId);
  axom::Array<SetPosition> index_vec(0, 0, m_unifiedAllocatorId);

  SetPosition curIdx = 0;

  for(auto i = 0; i < MAX_SET_SIZE1; ++i)
  {
    begin_vec[i] = curIdx;
    if(MAX_SET_SIZE1 / 4 <= i && i <= MAX_SET_SIZE1 / 4 * 3)
    {
      for(auto j = MAX_SET_SIZE2 / 4; j < MAX_SET_SIZE2 / 4 * 3; ++j)
      {
        index_vec.push_back(j);
        ++curIdx;
      }
    }
  }
  begin_vec[MAX_SET_SIZE1] = curIdx;

  rel[0].bindBeginOffsets(MAX_SET_SIZE1, begin_vec.view());
  rel[0].bindIndices(index_vec.size(), index_vec.view());

  RelationSetType relSet(&rel[0]);
  EXPECT_EQ(index_vec.size(), relSet.totalSize());
  EXPECT_TRUE(relSet.isValid());

  // Create array of elements to back the map.
  m_allocatorId = ExecTraits<ExecSpace>::getUnifiedAllocatorId();
  axom::IndexType backingSize = index_vec.size() * flatStride;

  RealData realBacking(backingSize, backingSize, m_allocatorId);

  SLIC_INFO(axom::fmt::format("\nCreating double map with shape ({}) on the ProductSet",
                              axom::fmt::join(shape, ", ")));
  const MapType m(relSet, realBacking.view(), shape);

  EXPECT_EQ(m.stride(), flatStride);
  axom::for_all<ExecSpace>(
    m.firstSetSize(),
    AXOM_LAMBDA(int idx1) {
      auto submap = m(idx1);
      for(auto slot = 0; slot < submap.size(); slot++)
      {
        for(int i = 0; i < shape[0]; i++)
        {
          for(int j = 0; j < shape[1]; j++)
          {
            for(int k = 0; k < shape[2]; k++)
            {
              int idx2 = submap.index(slot);
              int flatCompIdx = i * strides[0] + j * strides[1] + k * strides[2];
              submap(slot, i, j, k) = getVal<double>(idx1, idx2, flatCompIdx);
            }
          }
        }
      }
    });

  SLIC_INFO("\nChecking the elements with findValue().");
  for(int idx1 = 0; idx1 < m.firstSetSize(); idx1++)
  {
    auto submap = m(idx1);
    int slot = 0;
    for(int idx2 = 0; idx2 < m.secondSetSize(); idx2++)
    {
      bool inRelation = submap.size() > slot && submap.index(slot) == idx2;
      for(int i = 0; i < shape[0]; i++)
      {
        for(int j = 0; j < shape[1]; j++)
        {
          for(int k = 0; k < shape[2]; k++)
          {
            int flatCompIdx = i * strides[0] + j * strides[1] + k * strides[2];
            double expected_value = getVal<double>(idx1, idx2, flatCompIdx);

            double* valuePtr = m.findValue(idx1, idx2, i, j, k);
            if(inRelation)
            {
              // Test set-based indexing: (idx1, idx2)
              EXPECT_NE(valuePtr, nullptr);
              EXPECT_DOUBLE_EQ(expected_value, *valuePtr);
              EXPECT_DOUBLE_EQ(expected_value, m(idx1, idx2, i, j, k));
            }
            else
            {
              EXPECT_EQ(valuePtr, nullptr);
            }
          }
        }
      }
      if(inRelation)
      {
        slot++;
      }
    }
  }

  SLIC_INFO("\nChecking the elements with BivariateMap range iterator.");
  for(auto it = m.set_begin(); it != m.set_end(); ++it)
  {
    int idx1 = it.firstIndex();
    int idx2 = it.secondIndex();
    int flatIdx = it.flatIndex();

    EXPECT_EQ(idx1, m.set()->flatToFirstIndex(flatIdx));
    EXPECT_EQ(idx2, m.set()->flatToSecondIndex(flatIdx));
    EXPECT_EQ(it.numComp(), m.stride());

    for(int i = 0; i < shape[0]; i++)
    {
      for(int j = 0; j < shape[1]; j++)
      {
        for(int k = 0; k < shape[2]; k++)
        {
          int flatCompIdx = i * strides[0] + j * strides[1] + k * strides[2];
          double expected_value = getVal<double>(idx1, idx2, flatCompIdx);

          EXPECT_DOUBLE_EQ(expected_value, (*it)(i, j, k));
          EXPECT_DOUBLE_EQ(expected_value, it(i, j, k));
          EXPECT_DOUBLE_EQ(expected_value, it.value(i, j, k));
        }
      }
    }
  }
}

//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_bivariate_map_templated, constructAndTestProductSet)
{
  this->initializeAndTestCartesianMap(1);
  this->initializeAndTestCartesianMap(2);
  this->initializeAndTestCartesianMap(3);
}

//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_bivariate_map_templated, constructAndTestProductSet3D)
{
  this->initializeAndTestCartesianMap({2, 3, 5});
  this->initializeAndTestCartesianMap({3, 5, 7});
}

//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_bivariate_map_templated, constructAndTestRelationSet)
{
  this->initializeAndTestRelationMap(1);
  this->initializeAndTestRelationMap(2);
  this->initializeAndTestRelationMap(3);
}

//----------------------------------------------------------------------
AXOM_TYPED_TEST(slam_bivariate_map_templated, constructAndTestRelationSet3D)
{
  this->initializeAndTestRelationMap({2, 3, 5});
  this->initializeAndTestRelationMap({3, 5, 7});
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
