// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "benchmark/benchmark.h"
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <tuple>
#include <vector>
#include <type_traits>
#include <cstdint>

namespace
{

enum class ArrayFeatureBenchmarks
{
  None = 0,
  Constructors = 1 << 0,
  Insertion = 1 << 1,
  Iterators = 1 << 2,

  All = Constructors | Insertion | Iterators
};

inline ArrayFeatureBenchmarks operator|(ArrayFeatureBenchmarks lhs, ArrayFeatureBenchmarks rhs)
{
  using T = std::underlying_type_t<ArrayFeatureBenchmarks>;
  return static_cast<ArrayFeatureBenchmarks>(static_cast<T>(lhs) | static_cast<T>(rhs));
}

inline ArrayFeatureBenchmarks& operator|=(ArrayFeatureBenchmarks& lhs, ArrayFeatureBenchmarks rhs)
{
  lhs = lhs | rhs;
  return lhs;
}

inline ArrayFeatureBenchmarks operator&(ArrayFeatureBenchmarks lhs, ArrayFeatureBenchmarks rhs)
{
  using T = std::underlying_type_t<ArrayFeatureBenchmarks>;
  return static_cast<ArrayFeatureBenchmarks>(static_cast<T>(lhs) & static_cast<T>(rhs));
}

std::vector<int> args_benchmark_sizes;
ArrayFeatureBenchmarks args_benchmark_features {ArrayFeatureBenchmarks::None};
}  // namespace

// Custom fmt formatter for ArrayFeatureBenchmarks
template <>
struct axom::fmt::formatter<ArrayFeatureBenchmarks>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(ArrayFeatureBenchmarks feature, FormatContext& ctx) const
  {
    static const std::map<ArrayFeatureBenchmarks, std::string> feature_map = {
      {ArrayFeatureBenchmarks::Constructors, "Constructors"},
      {ArrayFeatureBenchmarks::Insertion, "Insertion"},
      {ArrayFeatureBenchmarks::Iterators, "Iterators"}};

    if(feature == ArrayFeatureBenchmarks::None)
    {
      return axom::fmt::format_to(ctx.out(), "None");
    }
    else if(feature == ArrayFeatureBenchmarks::All)
    {
      return axom::fmt::format_to(ctx.out(), "All");
    }

    std::string name;
    for(const auto& kv : feature_map)
    {
      if((feature & kv.first) != ArrayFeatureBenchmarks::None)
      {
        name += name.empty() ? kv.second : "|" + kv.second;
      }
    }

    return axom::fmt::format_to(ctx.out(), "{}", name);
  }
};

//-----------------------------------------------------------------------------
// Arguments and types to test
//-----------------------------------------------------------------------------

void CustomArgs(benchmark::internal::Benchmark* b)
{
  for(int sz : ::args_benchmark_sizes)
  {
    b->Arg(sz);
  }
}

template <typename T>
struct Wrapper
{
  T j;
  bool operator==(const Wrapper& other) const { return j == other.j; }
};

using Types = std::tuple<int, std::pair<int, int>, Wrapper<int>, std::string>;

//-----------------------------------------------------------------------------
// Helper functions for generating type names and values
//-----------------------------------------------------------------------------

template <typename T>
struct get_value_impl
{
  static T get(int i) { return static_cast<T>(i); }
};

template <typename T>
struct get_value_impl<std::pair<T, T>>
{
  static std::pair<T, T> get(int i)
  {
    return std::make_pair(static_cast<T>(i), static_cast<T>(i + 1));
  }
};

template <typename T>
struct get_value_impl<Wrapper<T>>
{
  static Wrapper<T> get(int i) { return Wrapper<T> {static_cast<T>(i)}; }
};

template <>
struct get_value_impl<std::string>
{
  static std::string get(int i)
  {
    return axom::fmt::format("This is a somewhat long string -- {}", i);
  }
};

template <typename T>
T get_value(int i)
{
  return get_value_impl<T>::get(i);
}

template <typename T>
std::string get_type_name();

template <>
std::string get_type_name<int>()
{
  return "int";
}

template <>
std::string get_type_name<std::pair<int, int>>()
{
  return "std::pair<int, int>";
}

template <>
std::string get_type_name<Wrapper<int>>()
{
  return "Wrapper<int>";
}

template <>
std::string get_type_name<std::string>()
{
  return "std::string";
}

// Helper function to check that an array's size and capacity match expectation
// Since we're benchmarking, let's do our best to ensure it's a no-op in release configs
template <typename Arr>
void assert_size_and_capacity(const Arr& arr, int exp_size, int exp_capacity)
{
#if defined(NDEBUG)
  AXOM_UNUSED_VAR(arr);
  AXOM_UNUSED_VAR(exp_size);
  AXOM_UNUSED_VAR(exp_capacity);
#else
  assert(static_cast<int>(arr.size()) == exp_size);
  assert(static_cast<int>(arr.capacity()) >= exp_capacity);
#endif
}

template <typename Arr>
void assert_size_and_strict_capacity(const Arr& arr, int exp_size, int exp_capacity)
{
#if defined(NDEBUG)
  AXOM_UNUSED_VAR(arr);
  AXOM_UNUSED_VAR(exp_size);
  AXOM_UNUSED_VAR(exp_capacity);
#else
  assert(static_cast<int>(arr.size()) == exp_size);
  assert(static_cast<int>(arr.capacity()) == exp_capacity);
#endif
}

//-----------------------------------------------------------------------------
// Benchmarks for array and vector construction
//-----------------------------------------------------------------------------
template <typename Container>
static void ctor(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    Container arr(size);
    assert_size_and_capacity(arr, size, size);
    benchmark::DoNotOptimize(arr);
  }
}

//-----------------------------------------------------------------------------
// Benchmarks for array and vector push_back and emplace_back
//-----------------------------------------------------------------------------
template <typename Container>
void push_back_startEmpty(benchmark::State& state)
{
  using T = typename Container::value_type;
  const int size = state.range(0);
  for(auto _ : state)
  {
    Container arr;
    assert_size_and_strict_capacity(arr, 0, 0);

    for(int i = 0; i < size; ++i)
    {
      arr.push_back(get_value<T>(i));  // using push_back on initially empty container
    }
    assert_size_and_capacity(arr, size, size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename Container>
void emplace_back_startEmpty(benchmark::State& state)
{
  using T = typename Container::value_type;
  const int size = state.range(0);
  for(auto _ : state)
  {
    Container arr;
    assert_size_and_strict_capacity(arr, 0, 0);
    for(int i = 0; i < size; ++i)
    {
      arr.emplace_back(get_value<T>(i));  // using emplace_back on initially empty container
    }
    assert_size_and_capacity(arr, size, size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename Container>
void push_back_initialReserve(benchmark::State& state)
{
  using T = typename Container::value_type;
  const int size = state.range(0);
  for(auto _ : state)
  {
    Container arr;
    arr.reserve(size);
    assert_size_and_strict_capacity(arr, 0, size);
    for(int i = 0; i < size; ++i)
    {
      arr.push_back(get_value<T>(i));  // using push_back on container w/ reserved size
    }
    assert_size_and_strict_capacity(arr, size, size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename Container>
void emplace_back_initialReserve(benchmark::State& state)
{
  using T = typename Container::value_type;
  const int size = state.range(0);
  for(auto _ : state)
  {
    Container arr;
    arr.reserve(size);
    assert_size_and_strict_capacity(arr, 0, size);
    for(int i = 0; i < size; ++i)
    {
      arr.emplace_back(get_value<T>(i));  // using emplace_back on container w/ reserved size
    }
    assert_size_and_strict_capacity(arr, size, size);
    benchmark::DoNotOptimize(arr);
  }
}

//-----------------------------------------------------------------------------
// Use test fixtures for testing other operations
//-----------------------------------------------------------------------------
template <typename Container>
struct ArrFixture
{
public:
  using T = typename Container::value_type;

  ArrFixture(int size)
  {
    data.reserve(size);
    for(int i = 0; i < size; ++i)
    {
      data.emplace_back(get_value<T>(i));
    }
    assert_size_and_capacity(data, size, size);
  }

  Container data;
};

template <typename Container>
void iterate_range(benchmark::State& state)
{
  using T = typename Container::value_type;
  const int size = state.range(0);
  ArrFixture<Container> fixture(size);

  const T element = get_value<T>(size / 2);
  for(auto _ : state)
  {
    int count = 0;
    for(const T& item : fixture.data)
    {
      if(item == element)
      {
        ++count;
      }
    }
    assert(count == 1);
    benchmark::DoNotOptimize(count);
  }
}

template <typename Container>
void iterate_direct(benchmark::State& state)
{
  using T = typename Container::value_type;
  const int size = state.range(0);
  ArrFixture<Container> fixture(size);

  const T element = get_value<T>(size / 2);
  for(auto _ : state)
  {
    int count = 0;
    for(auto it = fixture.data.begin(); it != fixture.data.end(); ++it)
    {
      if(*it == element)
      {
        ++count;
      }
    }
    assert(count == 1);
    benchmark::DoNotOptimize(count);
  }
}

//-----------------------------------------------------------------------------
// Register all the tests
//-----------------------------------------------------------------------------

template <typename T>
void RegisterBenchmark()
{
  auto tname = [](const std::string& n) {
    return axom::fmt::format("{}<{}>", n, get_type_name<T>());
  };

  // clang-format off
  if((args_benchmark_features & ArrayFeatureBenchmarks::Constructors) != ArrayFeatureBenchmarks::None)
  {
    benchmark::RegisterBenchmark(tname("Array::ctor"), &ctor<axom::Array<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("vector::ctor"), &ctor<std::vector<T>>)->Apply(CustomArgs);
  }

  if((args_benchmark_features & ArrayFeatureBenchmarks::Insertion) != ArrayFeatureBenchmarks::None)
  {
    benchmark::RegisterBenchmark(tname("Array::push_back_startEmpty"), &push_back_startEmpty<axom::Array<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("Array::emplace_back_startEmpty"), &emplace_back_startEmpty<axom::Array<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("Array::push_back_initialReserve"), & push_back_initialReserve<axom::Array<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("Array::emplace_back_initialReserve"), &emplace_back_initialReserve<axom::Array<T>>)->Apply(CustomArgs);

    benchmark::RegisterBenchmark(tname("vector::push_back_startEmpty"), &push_back_startEmpty<std::vector<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("vector::emplace_back_startEmpty"), &emplace_back_startEmpty<std::vector<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("vector::push_back_initialReserve"), &push_back_initialReserve<std::vector<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("vector::emplace_back_initialReserve"), &emplace_back_initialReserve<std::vector<T>>)->Apply(CustomArgs);
  }

  if((args_benchmark_features & ArrayFeatureBenchmarks::Iterators) != ArrayFeatureBenchmarks::None)
  {
    benchmark::RegisterBenchmark(tname("Array::iterate_range"), &iterate_range<axom::Array<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("Array::iterate_direct"), &iterate_direct<axom::Array<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("vector::iterate_range"), &iterate_range<std::vector<T>>)->Apply(CustomArgs);
    benchmark::RegisterBenchmark(tname("vector::iterate_direct"), &iterate_direct<std::vector<T>>)->Apply(CustomArgs);
  }
  // clang-format on
}

//-----------------------------------------------------------------------------
// Main and helper functions to consistently register templated test types
//
// Note: This can likely be simplified/streamlined a bit when we move to C++17,
// e.g. using std::apply and parameter pack fold expressions
//-----------------------------------------------------------------------------

template <typename Tuple, std::size_t... Is>
void RegisterBenchmarksImpl(std::index_sequence<Is...>)
{
  using expander = int[];
  (void)expander {0, (RegisterBenchmark<typename std::tuple_element<Is, Tuple>::type>(), 0)...};
}

template <typename Tuple>
void RegisterBenchmarks()
{
  RegisterBenchmarksImpl<Tuple>(std::make_index_sequence<std::tuple_size<Tuple>::value> {});
}

int main(int argc, char* argv[])
{
  // Parse command line args for array benchmarks
  std::vector<int> local_test_sizes;
  ArrayFeatureBenchmarks local_benchmark_features {ArrayFeatureBenchmarks::None};

  axom::CLI::App app {"Axom array benchmarks"};
  app.add_option("-s,--custom_sizes", local_test_sizes)
    ->description("Adds custom array sizes to benchmark (positive numbers only)")
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
                         local_test_sizes.push_back(1 << 25);  // larger than  25M L3 cache
                       })
    ->description("Test array sizes related to typical cache sizes");

  std::vector<std::string> feature_strings;
  auto feature_opt =
    app.add_option("-f,--features", feature_strings)
      ->description(
        "Features to benchmark (Constructors, Insertion, Iterators, All); "
        "default is 'All'")
      ->expected(-1)
      ->each([&local_benchmark_features](const std::string& feature) {
        static const std::map<std::string, ArrayFeatureBenchmarks> feature_map = {
          {"constructors", ArrayFeatureBenchmarks::Constructors},
          {"insertion", ArrayFeatureBenchmarks::Insertion},
          {"iterators", ArrayFeatureBenchmarks::Iterators},
          {"all", ArrayFeatureBenchmarks::All}};

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
    // copy list of Features to test into global array; by default, test everything
    ::args_benchmark_features =
      feature_opt->count() > 0 ? local_benchmark_features : ArrayFeatureBenchmarks::All;

    // sort and unique-ify the input sizes and copy into the global array variable
    std::sort(local_test_sizes.begin(), local_test_sizes.end());
    auto last = std::unique(local_test_sizes.begin(), local_test_sizes.end());
    local_test_sizes.erase(last, local_test_sizes.end());
    std::swap(::args_benchmark_sizes, local_test_sizes);

    SLIC_INFO("Parsed and processed command line arguments:");
    SLIC_INFO(axom::fmt::format("- Array sizes: {}", axom::fmt::join(::args_benchmark_sizes, ",")));
    SLIC_INFO(axom::fmt::format("- Array features to test: {}", ::args_benchmark_features));
  }

  RegisterBenchmarks<Types>();
  ::benchmark::RunSpecifiedBenchmarks();
  return 0;
}
