// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "benchmark/benchmark.h"
#include "axom/slic.hpp"

// Uncomment the following line to enable larger test sizes
// #define SLAM_BENCHMARK_LARGE_SIZES

//------------------------------------------------------------------------------
namespace
{
constexpr int STRIDE = 7;
constexpr int OFFSET = 12;

using IndexType = int;
using IndexArray = IndexType*;

using DataType = double;
using DataArray = DataType*;

// Generate an array of of size sz of indices in the range of [0,sz)
// NOTE: Caller must delete the array
IndexArray generateRandomPermutationArray(int sz, bool shouldPermute = false)
{
  IndexArray indices = new IndexType[sz];

  for(IndexType i = 0; i < sz; ++i)
  {
    indices[i] = i;
  }

  if(shouldPermute)
  {
    for(IndexType idx = 0; idx < sz; ++idx)
    {
      // find a random position in the array and swap value with current idx
      IndexType otherIdx = idx + rand() % (sz - idx);
      SLIC_ASSERT(otherIdx >= idx && otherIdx < sz);
      std::swap(indices[idx], indices[otherIdx]);
    }
  }

  for(IndexType i = 0; i < sz; ++i)
  {
    SLIC_ASSERT(indices[i] >= 0 && indices[i] < sz);
  }

  return indices;
}

// Generate an array of size sz of random doubles in the range [0,1)
// NOTE: Caller must delete the array
DataArray generateRandomDataField(int sz)
{
  const DataType rMaxDouble = static_cast<DataType>(RAND_MAX);

  DataArray data = new DataType[sz];

  for(IndexType i = 0; i < sz; ++i)
  {
    data[i] = rand() / rMaxDouble;
  }

  for(IndexType i = 0; i < sz; ++i)
  {
    SLIC_ASSERT(data[i] >= 0.0 && data[i] <= 1.0);
  }

  return data;
}

class SetFixture : public ::benchmark::Fixture
{
public:
  void SetUp(const ::benchmark::State& /*state*/) override
  {
    volatile int str_vol = STRIDE;  // pass through volatile variable so the
    str = str_vol;                  // number is not a compile time constant

    volatile int off_vol = OFFSET;  // pass through volatile variable so the
    off = off_vol;                  // number is not a compile time constant

    ind = nullptr;
    data = nullptr;
  }

  void TearDown(const ::benchmark::State& /*state*/) override
  {
    if(ind != nullptr)
    {
      delete[] ind;
      ind = nullptr;
    }
    if(data != nullptr)
    {
      delete[] data;
      data = nullptr;
    }

    SLIC_ASSERT(ind == nullptr);
    SLIC_ASSERT(data == nullptr);
  }

  ~SetFixture()
  {
    // Note: This destructor is called after the logger is finalized,
    //       so its messages are ignored.
  }

  int maxIndex(int sz) { return (sz * str + off); }

  int off;
  int str;
  IndexArray ind;
  DataArray data;
};

enum ArrSizes
{
  S0 = 1 << 3,   // small
  S1 = 1 << 16,  // larger than  32K L1 cache
  S2 = 1 << 19   // Larger than 256K L2 cache
#if defined(SLAM_BENCHMARK_LARGE_SIZES)
  ,
  S3 = 1 << 25  // Larger than  25M L3 cache
#endif
};

void CustomArgs(benchmark::internal::Benchmark* b)
{
  b->Arg(S0);
  b->Arg(S1);
  b->Arg(S2);
#if defined(SLAM_BENCHMARK_LARGE_SIZES)
  b->Arg(S3);
#endif
}

}  // namespace

/// -------------------

template <int SZ>
void contig_sequence_compileTimeSize(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    for(int i = 0; i < SZ; ++i)
    {
      IndexType pos = i;
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * SZ);
}
BENCHMARK_TEMPLATE(contig_sequence_compileTimeSize, S0);
BENCHMARK_TEMPLATE(contig_sequence_compileTimeSize, S1);
BENCHMARK_TEMPLATE(contig_sequence_compileTimeSize, S2);
#if defined(SLAM_BENCHMARK_LARGE_SIZES)
BENCHMARK_TEMPLATE(contig_sequence_compileTimeSize, S3);
#endif

BENCHMARK_DEFINE_F(SetFixture, contig_sequence)(benchmark::State& state)
{
  const int sz = state.range(0);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i;
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, contig_sequence)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, strided_sequence)(benchmark::State& state)
{
  const int sz = state.range(0);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i * str;
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, strided_sequence)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, offset_sequence)(benchmark::State& state)
{
  const int sz = state.range(0);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i + off;
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, offset_sequence)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, offset_strided_sequence)(benchmark::State& state)
{
  const int sz = state.range(0);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i * str + off;
      benchmark::DoNotOptimize(pos);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, offset_strided_sequence)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_ordered)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), false);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_ordered)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_ordered_strided)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), false);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i * str];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_ordered_strided)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_ordered_offset)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), false);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i + off];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_ordered_offset)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_ordered_strided_offset)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), false);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i * str + off];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_ordered_strided_offset)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_permuted)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), true);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_permuted)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_permuted_strided)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), true);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i * str];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_permuted_strided)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_permuted_offset)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), true);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i + off];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_permuted_offset)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_permuted_strided_offset)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(maxIndex(sz), true);

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i * str + off];
      benchmark::DoNotOptimize(pos);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_permuted_strided_offset)->Apply(CustomArgs);

/// --------------------  Benchmarks for array indexing ---------------------
BENCHMARK_DEFINE_F(SetFixture, contig_sequence_field)(benchmark::State& state)
{
  const int sz = state.range(0);

  data = generateRandomDataField(maxIndex(sz));

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i;
      benchmark::DoNotOptimize(data[pos]);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, contig_sequence_field)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, strided_sequence_field)(benchmark::State& state)
{
  const int sz = state.range(0);

  data = generateRandomDataField(maxIndex(sz));

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i * str;
      benchmark::DoNotOptimize(data[pos]);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, strided_sequence_field)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, offset_sequence_field)(benchmark::State& state)
{
  const int sz = state.range(0);

  data = generateRandomDataField(maxIndex(sz));

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i + off;
      benchmark::DoNotOptimize(data[pos]);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, offset_sequence_field)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, offset_strided_sequence_field)
(benchmark::State& state)
{
  const int sz = state.range(0);

  data = generateRandomDataField(maxIndex(sz));

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = i * str + off;
      benchmark::DoNotOptimize(data[pos]);
    }
  }
  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, offset_strided_sequence_field)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_ordered_field)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(sz, false);
  data = generateRandomDataField(maxIndex(sz));

  //    if(sz == 8)
  //    {
  //        std::cout<<"\n array indices (order)\n\t";
  //        for (int i=0; i < sz; ++i)
  //            std::cout<< "<" << ind[i] << "," << data [ ind[i] ] <<">\t";
  //        std::cout <<std::endl;
  //    }

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i];
      benchmark::DoNotOptimize(data[pos]);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_ordered_field)->Apply(CustomArgs);

BENCHMARK_DEFINE_F(SetFixture, indirection_sequence_permuted_field)
(benchmark::State& state)
{
  const int sz = state.range(0);

  ind = generateRandomPermutationArray(sz, true);
  data = generateRandomDataField(maxIndex(sz));

  //    if(sz == 8)
  //    {
  //        std::cout<<"\n array indices (permute)\n\t";
  //        for (int i=0; i < sz; ++i)
  //            std::cout<< "<" << ind[i] << "," << data [ ind[i] ] <<">\t";
  //        std::cout <<std::endl;
  //    }

  while(state.KeepRunning())
  {
    for(int i = 0; i < sz; ++i)
    {
      IndexType pos = ind[i];
      benchmark::DoNotOptimize(data[pos]);
    }
  }

  state.SetItemsProcessed(state.iterations() * sz);
}
BENCHMARK_REGISTER_F(SetFixture, indirection_sequence_permuted_field)->Apply(CustomArgs);

/// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  std::srand(std::time(NULL));

  ::benchmark::Initialize(&argc, argv);
  axom::slic::SimpleLogger logger;

  ::benchmark::RunSpecifiedBenchmarks();

  return 0;
}
