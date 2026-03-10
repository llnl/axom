// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "benchmark/benchmark.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "axom/core/FlatMap.hpp"
#include "axom/core/FlatMapUtil.hpp"
#include "axom/core/detail/FlatTable.hpp"

#if defined(AXOM_USE_SPARSEHASH)
  #include "axom/sparsehash/sparse_hash_map"
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{

using KeyType = std::int64_t;
using ValueType = std::int64_t;

enum class FlatMapFeatureBenchmarks
{
  None = 0,
  Insertion = 1 << 0,
  Lookup = 1 << 1,
  BatchedInsertion = 1 << 2,

  All = Insertion | Lookup | BatchedInsertion
};

inline FlatMapFeatureBenchmarks operator|(FlatMapFeatureBenchmarks lhs, FlatMapFeatureBenchmarks rhs)
{
  using T = std::underlying_type_t<FlatMapFeatureBenchmarks>;
  return static_cast<FlatMapFeatureBenchmarks>(static_cast<T>(lhs) | static_cast<T>(rhs));
}

inline FlatMapFeatureBenchmarks& operator|=(FlatMapFeatureBenchmarks& lhs,
                                            FlatMapFeatureBenchmarks rhs)
{
  lhs = lhs | rhs;
  return lhs;
}

inline FlatMapFeatureBenchmarks operator&(FlatMapFeatureBenchmarks lhs, FlatMapFeatureBenchmarks rhs)
{
  using T = std::underlying_type_t<FlatMapFeatureBenchmarks>;
  return static_cast<FlatMapFeatureBenchmarks>(static_cast<T>(lhs) & static_cast<T>(rhs));
}

std::vector<int> args_benchmark_sizes;
FlatMapFeatureBenchmarks args_benchmark_features {FlatMapFeatureBenchmarks::None};
int args_batch_size = 1 << 10;
}  // namespace

template <>
struct axom::fmt::formatter<FlatMapFeatureBenchmarks>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(FlatMapFeatureBenchmarks feature, FormatContext& ctx) const
  {
    static const std::map<FlatMapFeatureBenchmarks, std::string> feature_map = {
      {FlatMapFeatureBenchmarks::Insertion, "Insertion"},
      {FlatMapFeatureBenchmarks::Lookup, "Lookup"},
      {FlatMapFeatureBenchmarks::BatchedInsertion, "BatchedInsertion"}};

    if(feature == FlatMapFeatureBenchmarks::None)
    {
      return axom::fmt::format_to(ctx.out(), "None");
    }
    else if(feature == FlatMapFeatureBenchmarks::All)
    {
      return axom::fmt::format_to(ctx.out(), "All");
    }

    std::string name;
    for(const auto& kv : feature_map)
    {
      if((feature & kv.first) != FlatMapFeatureBenchmarks::None)
      {
        name += name.empty() ? kv.second : "|" + kv.second;
      }
    }
    return axom::fmt::format_to(ctx.out(), "{}", name);
  }
};

namespace
{

void CustomArgs(benchmark::internal::Benchmark* b)
{
  for(int sz : ::args_benchmark_sizes)
  {
    b->Arg(sz);
  }
}

std::vector<KeyType> make_shuffled_keys(int n, std::uint64_t seed)
{
  std::vector<KeyType> keys;
  keys.reserve(static_cast<std::size_t>(n));
  for(int i = 0; i < n; ++i)
  {
    keys.push_back(static_cast<KeyType>(i));
  }

  std::mt19937_64 rng(seed);
  std::shuffle(keys.begin(), keys.end(), rng);
  return keys;
}

std::vector<std::pair<KeyType, ValueType>> make_pairs(const std::vector<KeyType>& keys)
{
  std::vector<std::pair<KeyType, ValueType>> pairs;
  pairs.reserve(keys.size());
  for(std::size_t i = 0; i < keys.size(); ++i)
  {
    pairs.emplace_back(keys[i], static_cast<ValueType>(i));
  }
  return pairs;
}

std::vector<KeyType> make_miss_keys(const std::vector<KeyType>& keys, KeyType offset)
{
  std::vector<KeyType> misses;
  misses.reserve(keys.size());
  for(KeyType k : keys)
  {
    misses.push_back(k + offset);
  }
  return misses;
}

template <typename MapType>
struct MapFactory
{
  static MapType make_empty(std::size_t) { return MapType {}; }
  static void reserve(MapType&, std::size_t) { }
};

template <typename Key, typename Value, typename Hash, typename Eq, typename Alloc>
struct MapFactory<std::unordered_map<Key, Value, Hash, Eq, Alloc>>
{
  using MapType = std::unordered_map<Key, Value, Hash, Eq, Alloc>;
  static MapType make_empty(std::size_t) { return MapType {}; }
  static void reserve(MapType& map, std::size_t n) { map.reserve(n); }
};

template <typename Key, typename Value, typename Hash>
struct MapFactory<axom::FlatMap<Key, Value, Hash>>
{
  using MapType = axom::FlatMap<Key, Value, Hash>;
  static MapType make_empty(std::size_t) { return MapType {}; }
  static void reserve(MapType& map, std::size_t n) { map.reserve(static_cast<axom::IndexType>(n)); }
};

#if defined(AXOM_USE_SPARSEHASH)
template <typename Key, typename Value, typename Hash, typename Eq, typename Alloc>
void reserve_sparsehash(axom::google::sparse_hash_map<Key, Value, Hash, Eq, Alloc>& map, std::size_t n)
{
  map.max_load_factor(0.8f);
  const auto buckets_needed =
    static_cast<std::size_t>(static_cast<double>(n) / map.max_load_factor()) + 1;
  map.resize(buckets_needed);
}
#endif

#if defined(AXOM_USE_SPARSEHASH)
template <typename Key, typename Value, typename Hash, typename Eq, typename Alloc>
struct MapFactory<axom::google::sparse_hash_map<Key, Value, Hash, Eq, Alloc>>
{
  using MapType = axom::google::sparse_hash_map<Key, Value, Hash, Eq, Alloc>;
  static MapType make_empty(std::size_t) { return MapType {}; }
  static void reserve(MapType& map, std::size_t n) { reserve_sparsehash(map, n); }
};
#endif

template <typename MapType>
MapType make_reserved_map(std::size_t n)
{
  MapType map = MapFactory<MapType>::make_empty(n);
  MapFactory<MapType>::reserve(map, n);
  return map;
}

template <typename MapType>
MapType make_filled_map(const std::vector<std::pair<KeyType, ValueType>>& pairs)
{
  MapType map = make_reserved_map<MapType>(pairs.size());
  map.insert(pairs.begin(), pairs.end());
  return map;
}

template <typename Key, typename Value, typename Hash>
axom::FlatMap<Key, Value, Hash> make_filled_flatmap_with_target_load_factor(
  const std::vector<std::pair<KeyType, ValueType>>& pairs,
  double target_load_factor)
{
  using MapType = axom::FlatMap<Key, Value, Hash>;
  MapType map;

  const double max_lf = map.max_load_factor();
  const double lf = std::max(1e-3, std::min(target_load_factor, max_lf));
  const double n = static_cast<double>(pairs.size());

  // FlatMap's ctor/rehash argument is scaled internally by max_load_factor.
  // To target load factor `lf` for `n` elements, scale the count accordingly.
  const axom::IndexType rehash_count = static_cast<axom::IndexType>(std::ceil((n * max_lf) / lf));

  map.rehash(rehash_count);
  map.insert(pairs.begin(), pairs.end());
  return map;
}

template <typename MapType>
void BM_Insert_StartEmpty(benchmark::State& state)
{
  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0xA2D5B7C4ULL);
  const auto pairs = make_pairs(keys);

  for(auto _ : state)
  {
    MapType map = MapFactory<MapType>::make_empty(pairs.size());
    map.insert(pairs.begin(), pairs.end());
    benchmark::DoNotOptimize(map);
  }
}

template <typename MapType>
void BM_Insert_Reserved(benchmark::State& state)
{
  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0xA2D5B7C4ULL);
  const auto pairs = make_pairs(keys);

  for(auto _ : state)
  {
    MapType map = make_reserved_map<MapType>(pairs.size());
    map.insert(pairs.begin(), pairs.end());
    benchmark::DoNotOptimize(map);
  }
}

template <typename MapType>
void BM_Find_Hit(benchmark::State& state)
{
  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0xC0FFEEULL);
  const auto pairs = make_pairs(keys);
  const MapType map = make_filled_map<MapType>(pairs);

  for(auto _ : state)
  {
    ValueType sum = 0;
    for(KeyType k : keys)
    {
      auto it = map.find(k);
      if(it != map.end())
      {
        sum += it->second;
      }
    }
    benchmark::DoNotOptimize(sum);
  }
}

template <typename MapType>
void BM_Find_Miss(benchmark::State& state)
{
  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0xC0FFEEULL);
  const auto pairs = make_pairs(keys);
  const MapType map = make_filled_map<MapType>(pairs);
  const auto miss_keys = make_miss_keys(keys, static_cast<KeyType>(n) + 11);

  for(auto _ : state)
  {
    std::int64_t misses = 0;
    for(KeyType k : miss_keys)
    {
      misses += (map.find(k) == map.end()) ? 1 : 0;
    }
    benchmark::DoNotOptimize(misses);
  }
}

template <typename Hash>
void BM_FlatMap_Find_Hit_TargetLoad(benchmark::State& state, double target_load_factor)
{
  using MapType = axom::FlatMap<KeyType, ValueType, Hash>;

  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0xC0FFEEULL);
  const auto pairs = make_pairs(keys);
  const MapType map =
    make_filled_flatmap_with_target_load_factor<KeyType, ValueType, Hash>(pairs, target_load_factor);

  for(auto _ : state)
  {
    ValueType sum = 0;
    for(KeyType k : keys)
    {
      auto it = map.find(k);
      if(it != map.end())
      {
        sum += it->second;
      }
    }
    benchmark::DoNotOptimize(sum);
  }
}

void BM_FlatMap_Find_Hit_Prehashed(benchmark::State& state)
{
  using MapType = axom::FlatMap<KeyType, ValueType>;
  using HashResult = typename MapType::hash_result_type;

  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0xC0FFEEULL);
  const auto pairs = make_pairs(keys);
  const MapType map = make_filled_map<MapType>(pairs);

  std::vector<HashResult> hashes;
  hashes.reserve(keys.size());
  for(KeyType k : keys)
  {
    hashes.push_back(typename MapType::hasher {}(k));
  }

  for(auto _ : state)
  {
    ValueType sum = 0;
    for(std::size_t i = 0; i < keys.size(); ++i)
    {
      auto it = map.find_with_hash(keys[i], hashes[i]);
      if(it != map.end())
      {
        sum += it->second;
      }
    }
    benchmark::DoNotOptimize(sum);
  }
}

void BM_FlatMap_Find_Miss_Prehashed(benchmark::State& state)
{
  using MapType = axom::FlatMap<KeyType, ValueType>;
  using HashResult = typename MapType::hash_result_type;

  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0xC0FFEEULL);
  const auto pairs = make_pairs(keys);
  const MapType map = make_filled_map<MapType>(pairs);
  const auto miss_keys = make_miss_keys(keys, static_cast<KeyType>(n) + 11);

  std::vector<HashResult> hashes;
  hashes.reserve(miss_keys.size());
  for(KeyType k : miss_keys)
  {
    hashes.push_back(typename MapType::hasher {}(k));
  }

  for(auto _ : state)
  {
    std::int64_t misses = 0;
    for(std::size_t i = 0; i < miss_keys.size(); ++i)
    {
      misses += (map.find_with_hash(miss_keys[i], hashes[i]) == map.end()) ? 1 : 0;
    }
    benchmark::DoNotOptimize(misses);
  }
}

template <typename MapType>
void insert_pairs_in_batches(MapType& map,
                             const std::vector<std::pair<KeyType, ValueType>>& pairs,
                             int batch_size)
{
  const std::size_t n = pairs.size();
  const std::size_t bs = static_cast<std::size_t>(std::max(1, batch_size));
  for(std::size_t offset = 0; offset < n; offset += bs)
  {
    const std::size_t count = std::min(bs, n - offset);
    map.insert(pairs.begin() + static_cast<std::ptrdiff_t>(offset),
               pairs.begin() + static_cast<std::ptrdiff_t>(offset + count));
  }
}

template <typename Key, typename Value, typename Hash>
void insert_pairs_in_batches(axom::FlatMap<Key, Value, Hash>& map,
                             const std::vector<std::pair<KeyType, ValueType>>& pairs,
                             int batch_size)
{
  const std::size_t n = pairs.size();
  const std::size_t bs = static_cast<std::size_t>(std::max(1, batch_size));
  for(std::size_t offset = 0; offset < n; offset += bs)
  {
    const std::size_t count = std::min(bs, n - offset);
    map.template insert<axom::SEQ_EXEC>(pairs.begin() + static_cast<std::ptrdiff_t>(offset),
                                        pairs.begin() + static_cast<std::ptrdiff_t>(offset + count));
  }
}

template <typename MapType>
void BM_BatchedInsert_Reserved(benchmark::State& state)
{
  const int n = state.range(0);
  const auto keys = make_shuffled_keys(n, 0x1CEB00DAULL);
  const auto pairs = make_pairs(keys);

  for(auto _ : state)
  {
    MapType map = make_reserved_map<MapType>(pairs.size());
    insert_pairs_in_batches(map, pairs, ::args_batch_size);
    benchmark::DoNotOptimize(map);
  }
}

}  // namespace

//-----------------------------------------------------------------------------
// Register benchmarks
//-----------------------------------------------------------------------------

template <typename MapType>
void RegisterBenchmarksFor(const std::string& map_name)
{
  auto name = [&map_name](const std::string& op) {
    return axom::fmt::format("{}::{}", map_name, op);
  };

  // clang-format off
  if((::args_benchmark_features & FlatMapFeatureBenchmarks::Insertion) != FlatMapFeatureBenchmarks::None)
  {
    benchmark::RegisterBenchmark(name("insert_startEmpty"), &BM_Insert_StartEmpty<MapType>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(name("insert_reserved"), &BM_Insert_Reserved<MapType>)->Apply(CustomArgs);
  }

  if((::args_benchmark_features & FlatMapFeatureBenchmarks::Lookup) != FlatMapFeatureBenchmarks::None)
  {
    benchmark::RegisterBenchmark(name("find_hit"), &BM_Find_Hit<MapType>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(name("find_miss"), &BM_Find_Miss<MapType>)->Apply(CustomArgs);
  }

  if((::args_benchmark_features & FlatMapFeatureBenchmarks::BatchedInsertion) != FlatMapFeatureBenchmarks::None)
  {
    benchmark::RegisterBenchmark(name("insert_batched_reserved"), &BM_BatchedInsert_Reserved<MapType>)->Apply(CustomArgs);
  }
  // clang-format on
}

void RegisterFlatMapPrehashedBenchmarks()
{
  auto name = [](const std::string& op) { return axom::fmt::format("axom::FlatMap::{}", op); };

  if((::args_benchmark_features & FlatMapFeatureBenchmarks::Lookup) != FlatMapFeatureBenchmarks::None)
  {
    benchmark::RegisterBenchmark(name("find_hit_prehashed"), &BM_FlatMap_Find_Hit_Prehashed)
      ->Apply(CustomArgs);
    benchmark::RegisterBenchmark(name("find_miss_prehashed"), &BM_FlatMap_Find_Miss_Prehashed)
      ->Apply(CustomArgs);
  }
}

int main(int argc, char* argv[])
{
  std::vector<int> local_test_sizes;
  FlatMapFeatureBenchmarks local_benchmark_features {FlatMapFeatureBenchmarks::None};
  int local_batch_size = ::args_batch_size;

  axom::CLI::App app {"Axom FlatMap benchmarks"};
  app.add_option("-s,--custom_sizes", local_test_sizes)
    ->description("Adds custom map sizes to benchmark (positive numbers only)")
    ->expected(-1)
    ->default_val(std::vector<int> {1 << 16})
    ->each([](const std::string& num_str) {
      int num = std::stoi(num_str);
      if(num < 0)
      {
        throw axom::CLI::ValidationError("Negative numbers are not allowed");
      }
    });

  app
    .add_flag_callback("--use_cache_related_sizes",
                       [&local_test_sizes]() {
                         local_test_sizes.push_back(1 << 3);   // small
                         local_test_sizes.push_back(1 << 16);  // larger than  32K L1 cache
                         local_test_sizes.push_back(1 << 19);  // larger than 256K L2 cache
                         //local_test_sizes.push_back(1 << 25);  // larger than  25M L3 cache
                       })
    ->description("Test map sizes related to typical cache sizes");

  app.add_option("--batch_size", local_batch_size)
    ->description("Batch size for batched insertion benchmarks")
    ->default_val(local_batch_size)
    ->check(axom::CLI::PositiveNumber);

  std::vector<std::string> feature_strings;
  auto feature_opt =
    app.add_option("-f,--features", feature_strings)
      ->description(
        "Features to benchmark (Insertion, Lookup, BatchedInsertion, All); default is 'All'")
      ->expected(-1)
      ->each([&local_benchmark_features](const std::string& feature) {
        static const std::map<std::string, FlatMapFeatureBenchmarks> feature_map = {
          {"insertion", FlatMapFeatureBenchmarks::Insertion},
          {"lookup", FlatMapFeatureBenchmarks::Lookup},
          {"batchedinsertion", FlatMapFeatureBenchmarks::BatchedInsertion},
          {"all", FlatMapFeatureBenchmarks::All}};

        std::string lower_feature = feature;
        std::transform(lower_feature.begin(), lower_feature.end(), lower_feature.begin(), ::tolower);
        auto it = feature_map.find(lower_feature);
        if(it == feature_map.end())
        {
          throw axom::CLI::ValidationError("Invalid feature: " + feature);
        }

        local_benchmark_features |= it->second;
      });

  app.allow_extras();  // pass additional args to gbenchmark
  CLI11_PARSE(app, argc, argv);

  ::benchmark::Initialize(&argc, argv);
  axom::slic::SimpleLogger logger;

  // process input into global variables
  {
    ::args_benchmark_features =
      feature_opt->count() > 0 ? local_benchmark_features : FlatMapFeatureBenchmarks::All;

    std::sort(local_test_sizes.begin(), local_test_sizes.end());
    auto last = std::unique(local_test_sizes.begin(), local_test_sizes.end());
    local_test_sizes.erase(last, local_test_sizes.end());
    std::swap(::args_benchmark_sizes, local_test_sizes);

    ::args_batch_size = local_batch_size;

    SLIC_INFO("Parsed and processed command line arguments:");
    SLIC_INFO(axom::fmt::format("- Map sizes: {}", axom::fmt::join(::args_benchmark_sizes, ",")));
    SLIC_INFO(axom::fmt::format("- Batch size: {}", ::args_batch_size));
    SLIC_INFO(axom::fmt::format("- Map features to test: {}", ::args_benchmark_features));
  }

  RegisterBenchmarksFor<axom::FlatMap<KeyType, ValueType>>("axom::FlatMap");
  RegisterFlatMapPrehashedBenchmarks();
  using FastHash = axom::detail::flat_map::FastHashMixer64<KeyType, axom::DeviceHash>;
  RegisterBenchmarksFor<axom::FlatMap<KeyType, ValueType, FastHash>>("axom::FlatMapFastHash");

  // Explore the impact of lower load factors on successful lookups.
  // This trades memory for potentially fewer probes and fewer cache misses.
  using DefaultHash = axom::FlatMap<KeyType, ValueType>::hasher;
  benchmark::RegisterBenchmark("axom::FlatMap::find_hit_lf0p50", [](benchmark::State& st) {
    BM_FlatMap_Find_Hit_TargetLoad<DefaultHash>(st, 0.50);
  })->Apply(CustomArgs);
  benchmark::RegisterBenchmark("axom::FlatMap::find_hit_lf0p70", [](benchmark::State& st) {
    BM_FlatMap_Find_Hit_TargetLoad<DefaultHash>(st, 0.70);
  })->Apply(CustomArgs);
  benchmark::RegisterBenchmark("axom::FlatMapFastHash::find_hit_lf0p50", [](benchmark::State& st) {
    BM_FlatMap_Find_Hit_TargetLoad<FastHash>(st, 0.50);
  })->Apply(CustomArgs);
  benchmark::RegisterBenchmark("axom::FlatMapFastHash::find_hit_lf0p70", [](benchmark::State& st) {
    BM_FlatMap_Find_Hit_TargetLoad<FastHash>(st, 0.70);
  })->Apply(CustomArgs);

  RegisterBenchmarksFor<std::unordered_map<KeyType, ValueType>>("std::unordered_map");
  RegisterBenchmarksFor<std::map<KeyType, ValueType>>("std::map");

#if defined(AXOM_USE_SPARSEHASH)
  RegisterBenchmarksFor<axom::google::sparse_hash_map<KeyType, ValueType>>(
    "axom::google::sparse_hash_map");
#endif

  ::benchmark::RunSpecifiedBenchmarks();
  return 0;
}
