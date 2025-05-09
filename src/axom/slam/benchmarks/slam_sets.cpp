// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "benchmark/benchmark.h"
#include "axom/slic.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/OrderedSet.hpp"

//------------------------------------------------------------------------------
namespace
{
namespace slam = axom::slam;

// constexpr int STRIDE = 7;
// constexpr int OFFSET = 12;

using IndexType = int;
using IndexArray = IndexType*;

using DataType = double;
using DataArray = DataType*;
/*
    // Generate an array of of size sz of indices in the range of [0,sz)
    // NOTE: Caller must delete the array
    IndexArray generateRandomPermutationArray(int sz, bool shouldPermute =
       false)
    {
        IndexArray indices = new IndexType[sz];

        for(IndexType i=0; i< sz; ++i)
        {
            indices[i] = i;
        }

        if(shouldPermute)
        {
            for(IndexType idx=0; idx< sz; ++idx)
            {
                // find a random position in the array and swap value with
                   current idx
                IndexType otherIdx = idx + rand() % (sz - idx);
                std::swap(indices[idx], indices[otherIdx]);
            }
        }

        return indices;
    }

    // Generate an array of size sz of random doubles in the range [0,1)
    // NOTE: Caller must delete the array
    DataArray generateRandomDataField(int sz)
    {
        const DataType rMaxDouble = static_cast<DataType>(RAND_MAX);

        DataArray data = new DataType[sz];
        for(IndexType i=0; i< sz; ++i)
        {
            data[i] = rand() / rMaxDouble;;
        }
        return data;
    }

    class SetFixture : public ::benchmark::Fixture
    {
    public:
        void SetUp() {

            volatile int str_vol = STRIDE;  // pass through volatile variable so the
            str = str_vol;                  // number is not a compile time constant

            volatile int off_vol = OFFSET;  // pass through volatile variable so the
            off = off_vol;                  // number is not a compile time constant

            ind = nullptr;
            data = nullptr;

        }

        void TearDown() {
            if(ind != nullptr)
                delete[] ind;
            if(data != nullptr)
                delete[] data;
        }

        ~SetFixture() {
            SLIC_ASSERT( ind == nullptr);
            SLIC_ASSERT( data == nullptr);
        }


        int maxIndex(int sz) { return (sz * str + off); }

        int off;
        int str;
        IndexArray ind;
        DataArray data;

    };
 */

enum ArrSizes
{
  S0 = 1 << 3,   // small
  S1 = 1 << 16,  // larger than  32K L1 cache
  S2 = 1 << 19,  // Larger than 256K L2 cache
  S3 = 1 << 25   // Larger than  25M L3 cache
};

void CustomArgs(benchmark::internal::Benchmark* b)
{
  b->Arg(S0);
  b->Arg(S1);
  b->Arg(S2);
  b->Arg(S3);
}
}  // namespace
//------------------------------------------------------------------------------

/// --------------------  Benchmarks for array indexing ---------------------

template <int SZ>
void positionSet_compileTimeSize(benchmark::State& state)
{
  using StaticSetSize = slam::policies::CompileTimeSize<int, SZ>;
  using SetType = slam::OrderedSet<int, int, StaticSetSize>;
  SetType set(SZ);

  while(state.KeepRunning())
  {
    for(int i = 0; i < set.size(); ++i)
    {
      int pos = set[i];
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * set.size());
}
BENCHMARK_TEMPLATE(positionSet_compileTimeSize, S0);
BENCHMARK_TEMPLATE(positionSet_compileTimeSize, S1);
BENCHMARK_TEMPLATE(positionSet_compileTimeSize, S2);
BENCHMARK_TEMPLATE(positionSet_compileTimeSize, S3);

template <int SZ>
void positionSet_runtimeTimeSize_template(benchmark::State& state)
{
  using SetType = slam::OrderedSet<>;
  SetType set(SZ);

  while(state.KeepRunning())
  {
    for(int i = 0; i < set.size(); ++i)
    {
      int pos = set[i];
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * set.size());
}
BENCHMARK_TEMPLATE(positionSet_runtimeTimeSize_template, S0);
BENCHMARK_TEMPLATE(positionSet_runtimeTimeSize_template, S1);
BENCHMARK_TEMPLATE(positionSet_runtimeTimeSize_template, S2);
BENCHMARK_TEMPLATE(positionSet_runtimeTimeSize_template, S3);

void positionSet_runtimeTimeSize_function(benchmark::State& state)
{
  using SetType = slam::OrderedSet<>;
  SetType set(state.range(0));

  while(state.KeepRunning())
  {
    for(int i = 0; i < set.size(); ++i)
    {
      int pos = set[i];
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * set.size());
}
BENCHMARK(positionSet_runtimeTimeSize_function)->Apply(CustomArgs);

void positionSet_runtimeTimeSize_function_sizeOutside(benchmark::State& state)
{
  using SetType = slam::OrderedSet<>;
  SetType set(state.range(0));

  const int sz = set.size();
  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      int pos = set[i];
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * set.size());
}
BENCHMARK(positionSet_runtimeTimeSize_function_sizeOutside)->Apply(CustomArgs);

void positionSet_runtimeTimeSize_function_volatileSizeOutside(benchmark::State& state)
{
  using SetType = slam::OrderedSet<>;
  SetType set(state.range(0));

  volatile int sz = set.size();
  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      int pos = set[i];
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * set.size());
}
BENCHMARK(positionSet_runtimeTimeSize_function_volatileSizeOutside)->Apply(CustomArgs);

void positionSet_runtimeTimeSize_iter(benchmark::State& state)
{
  using SetType = slam::OrderedSet<>;
  using SetIter = SetType::iterator;
  SetType set(state.range(0));

  while(state.KeepRunning())
  {
    for(SetIter it = set.begin(), itEnd = set.end(); it != itEnd; ++it)
    {
      int pos = *it;
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * set.size());
}
BENCHMARK(positionSet_runtimeTimeSize_iter)->Apply(CustomArgs);

int main(int argc, char* argv[])
{
  std::srand(std::time(NULL));
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();

  return 0;
}
