// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_hilbert_codec.cpp
 * \brief Benchmarks public and scalar Hilbert conversions in 2D, 3D, and 4D.
 *
 * Use Google Benchmark's --benchmark_filter and --benchmark_repetitions options
 * for timing in Release. Run --verify-only to check the adapters.
 *
 * Each timed iteration processes a batch. Divide its time by
 * `batch * operations_per_item` to get the average time per conversion.
 * Mixed cases encode a point and decode an independent key for each item.
 */

#include "axom/config.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"
#include "axom/slic.hpp"
#include "axom/spin/HilbertIndex.hpp"

#include "benchmark/benchmark.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace
{

// Fixed base seeds keep verification and benchmark inputs reproducible.
constexpr std::uint64_t VERIFICATION_SEED = 42;
// Retain the original benchmark's seed to preserve its input sequences.
constexpr std::uint64_t BENCHMARK_SEED = 0x7281268;

constexpr int SMALL_BATCH_SIZE = 256;
constexpr int LARGE_BATCH_SIZE = 65536;
constexpr int CLUSTER_COUNT = 8;

// Offsets are 0..63 per coordinate and 0..4095 per key, giving clusters of
// 64 grid positions per axis or 4096 consecutive key values. The depth masks
// wrap the resulting coordinates and keys into the valid domain.
constexpr std::uint32_t CLUSTER_COORD_OFFSET_MASK = 63;
constexpr std::uint64_t CLUSTER_KEY_OFFSET_MASK = 4095;

constexpr std::size_t CACHE_PRESSURE_BUFFER_BYTES = 8 * 1024 * 1024;
constexpr std::size_t CACHE_PRESSURE_STRIDE_BYTES = 64;
constexpr int CACHE_PRESSURE_ITERATIONS = 32;

// Public follows Hilbertizer's dispatch. Scalar calls the reference transform
// directly to keep a baseline as dispatch changes. Both allow normal inlining.
enum Method
{
  Public,
  Scalar
};

const char* methodName(Method method) { return method == Public ? "Public" : "Scalar"; }

template <int DIM, typename U, Method METHOD>
struct Codec
{
  using Hilbert = axom::spin::Hilbertizer<std::uint32_t, U, DIM>;
  using ScalarTransform = axom::spin::internal::hilbert::HilbertTransform<std::uint32_t, U, DIM>;
  using Point = typename Hilbert::PointType;

  static U encode(const Point& point, int bits)
  {
    if constexpr(METHOD == Public)
    {
      return Hilbert::toIndex(point, bits);
    }
    else
    {
      return ScalarTransform::toIndex(point, bits);
    }
  }

  static Point decode(U key, int bits)
  {
    if constexpr(METHOD == Public)
    {
      return Hilbert::toPoint(key, bits);
    }
    else
    {
      return ScalarTransform::toPoint(key, bits);
    }
  }
};

enum Direction
{
  Encode,
  Decode,
  Mixed
};
enum Workload
{
  Warm,
  Small,
  Pressure,
  Variable
};
enum Dataset
{
  Uniform,
  Sequential,
  Clustered
};

// The scalar transform defines the expected ordering for these compatibility checks.
template <int DIM, typename U, Method METHOD>
bool verifyCodec()
{
  using CodecType = Codec<DIM, U, METHOD>;
  using Hilbert = typename CodecType::Hilbert;
  using ScalarTransform = typename CodecType::ScalarTransform;
  constexpr int max_bits = Hilbert::MAX_UNIQUE_BITS;
  std::mt19937_64 random(VERIFICATION_SEED + DIM);
  for(int bits = 1; bits <= max_bits; ++bits)
  {
    // Start with all bits set to avoid shifting 1 by the full key width.
    const U key_mask = std::numeric_limits<U>::max() >> (sizeof(U) * 8 - bits * DIM);
    const auto coord_mask = static_cast<std::uint32_t>((std::uint64_t {1} << bits) - 1);
    const bool exhaustive = bits * DIM <= 12;
    const unsigned count = exhaustive ? (1u << (bits * DIM)) : 2048;
    for(unsigned i = 0; i < count; ++i)
    {
      // Enumerate coordinates independently of the decoder in small domains.
      // At larger depths, cover corners and then deterministic random points.
      typename CodecType::Point point;
      for(int axis = 0; axis < DIM; ++axis)
      {
        point[axis] = exhaustive
          ? static_cast<std::uint32_t>((std::uint64_t {i} >> (axis * bits)) & coord_mask)
          : (i < (1u << DIM) ? (((i >> axis) & 1) ? coord_mask : 0)
                             : static_cast<std::uint32_t>(random()) & coord_mask);
      }
      // Independent keys include zero, the maximum, and every single-bit key.
      const U key = exhaustive
        ? static_cast<U>(i)
        : (i == 0                                         ? U {0}
             : i == 1                                     ? key_mask
             : i <= static_cast<unsigned>(bits * DIM + 1) ? U {1} << (i - 2)
                                                          : static_cast<U>(random()) & key_mask);
      const U encoded = CodecType::encode(point, bits);
      const auto decoded = CodecType::decode(key, bits);
      if(!Hilbert::isValidIndex(encoded, bits) || !Hilbert::isValidPoint(decoded, bits) ||
         encoded != ScalarTransform::toIndex(point, bits) ||
         decoded != ScalarTransform::toPoint(key, bits) ||
         CodecType::encode(decoded, bits) != key || CodecType::decode(encoded, bits) != point)
      {
        SLIC_WARNING(
          axom::fmt::format("Verification failed: {} D={} key bits={} depth={} sample={}",
                            methodName(METHOD),
                            DIM,
                            sizeof(U) * 8,
                            bits,
                            i));
        return false;
      }
    }
  }
  SLIC_INFO(
    axom::fmt::format("Verified {}/{}D/u{} at every depth", methodName(METHOD), DIM, sizeof(U) * 8));
  return true;
}

template <int DIM, typename U>
bool verifyCodecsForShape()
{
  return verifyCodec<DIM, U, Public>() && verifyCodec<DIM, U, Scalar>();
}

template <int DIM, typename U, Method METHOD, Direction DIR, Workload WORK, Dataset DATA>
void benchmarkHilbertCodec(benchmark::State& state)
{
  using CodecType = Codec<DIM, U, METHOD>;
  using Hilbert = typename CodecType::Hilbert;
  using Point = typename CodecType::Point;
  constexpr int max_bits = Hilbert::MAX_UNIQUE_BITS;
  const int count = WORK == Small || WORK == Pressure ? SMALL_BATCH_SIZE : LARGE_BATCH_SIZE;
  // Allocation and input generation precede the state loop, so neither is timed.
  std::vector<Point> points(count), decoded(count);
  std::vector<U> keys(count), encoded(count);
  std::vector<int> depths(count);
  // Vary the sequence by dimension. Both codecs use the same seed.
  std::mt19937_64 random(BENCHMARK_SEED + DIM);
  Point centers[CLUSTER_COUNT];
  U key_centers[CLUSTER_COUNT];
  for(int c = 0; c < CLUSTER_COUNT; ++c)
  {
    for(int a = 0; a < DIM; ++a)
    {
      centers[c][a] = static_cast<std::uint32_t>(random());
    }
    key_centers[c] = static_cast<U>(random());
  }
  for(int i = 0; i < count; ++i)
  {
    const int bits = WORK == Variable ? 1 + static_cast<int>(random() % max_bits) : max_bits;
    depths[i] = bits;
    const auto cmask = static_cast<std::uint32_t>((std::uint64_t {1} << bits) - 1);
    const U kmask = std::numeric_limits<U>::max() >> (sizeof(U) * 8 - bits * DIM);
    for(int a = 0; a < DIM; ++a)
    {
      if constexpr(DATA == Sequential)
      {
        points[i][a] =
          static_cast<std::uint32_t>((static_cast<std::uint64_t>(i) >> (a * bits)) & cmask);
      }
      else if constexpr(DATA == Clustered)
      {
        points[i][a] = (centers[i % CLUSTER_COUNT][a] +
                        static_cast<std::uint32_t>(random() & CLUSTER_COORD_OFFSET_MASK)) &
          cmask;
      }
      else
      {
        points[i][a] = static_cast<std::uint32_t>(random()) & cmask;
      }
    }
    if constexpr(DATA == Sequential)
    {
      keys[i] = static_cast<U>(i) & kmask;
    }
    else if constexpr(DATA == Clustered)
    {
      keys[i] =
        (key_centers[i % CLUSTER_COUNT] + static_cast<U>(random() & CLUSTER_KEY_OFFSET_MASK)) & kmask;
    }
    else
    {
      keys[i] = static_cast<U>(random()) & kmask;
    }
  }
  std::vector<std::uint64_t> pressure(
    WORK == Pressure ? CACHE_PRESSURE_BUFFER_BYTES / sizeof(std::uint64_t) : 0);
  constexpr auto pressure_stride = CACHE_PRESSURE_STRIDE_BYTES / sizeof(std::uint64_t);
  // Expose the buffers to the compiler barriers. ClobberMemory below keeps the
  // output stores in each batch without adding a barrier to every conversion.
  benchmark::DoNotOptimize(points.data());
  benchmark::DoNotOptimize(keys.data());
  benchmark::DoNotOptimize(encoded.data());
  benchmark::DoNotOptimize(decoded.data());
  benchmark::DoNotOptimize(pressure.data());
  for(auto _ : state)
  {
    if constexpr(WORK == Pressure)
    {
      // Touch the pressure buffer at a fixed byte stride between timed batches.
      // This adds cache pressure; the amount of eviction depends on the machine.
      state.PauseTiming();
      for(std::size_t i = 0; i < pressure.size(); i += pressure_stride)
      {
        pressure[i] += i + 1;
      }
      benchmark::ClobberMemory();
      state.ResumeTiming();
    }
    for(int i = 0; i < count; ++i)
    {
      // Fixed-depth cases let the compiler specialize the codec. Variable reads
      // a depth per item to exercise calls whose depth is only known at runtime.
      const int bits = WORK == Variable ? depths[i] : max_bits;
      if constexpr(DIR == Encode || DIR == Mixed)
      {
        encoded[i] = CodecType::encode(points[i], bits);
      }
      // Independent keys prevent an encode/decode round trip from being optimized away.
      if constexpr(DIR == Decode || DIR == Mixed)
      {
        decoded[i] = CodecType::decode(keys[i], bits);
      }
    }
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * count * (DIR == Mixed ? 2 : 1));
  state.counters["batch"] = count;
  state.counters["operations_per_item"] = DIR == Mixed ? 2 : 1;
}

template <int DIM, typename U, Method METHOD, Direction DIR, Workload WORK, Dataset DATA>
void registerBenchmark()
{
  const char* directions[] = {"Encode", "Decode", "Mixed"};
  const char* workloads[] = {"Warm", "Small", "Pressure", "Variable"};
  const char* datasets[] = {"Uniform", "Sequential", "Clustered"};
  const auto label = axom::fmt::format("{}/{}D/u{}/{}/{}/{}",
                                       methodName(METHOD),
                                       DIM,
                                       sizeof(U) * 8,
                                       directions[DIR],
                                       workloads[WORK],
                                       datasets[DATA]);
  auto* b = benchmark::RegisterBenchmark(label.c_str(),
                                         &benchmarkHilbertCodec<DIM, U, METHOD, DIR, WORK, DATA>);
  b->UseRealTime()->Unit(benchmark::kNanosecond);
  // Limit the untimed buffer walks per repetition.
  if constexpr(WORK == Pressure)
  {
    b->Iterations(CACHE_PRESSURE_ITERATIONS);
  }
}

template <int DIM, typename U, Method METHOD, Direction DIR>
void registerWorkloads()
{
  registerBenchmark<DIM, U, METHOD, DIR, Warm, Uniform>();
  registerBenchmark<DIM, U, METHOD, DIR, Small, Uniform>();
  registerBenchmark<DIM, U, METHOD, DIR, Pressure, Uniform>();
  registerBenchmark<DIM, U, METHOD, DIR, Variable, Uniform>();
  registerBenchmark<DIM, U, METHOD, DIR, Warm, Sequential>();
  registerBenchmark<DIM, U, METHOD, DIR, Warm, Clustered>();
}

template <int DIM, typename U, Method METHOD>
void registerMethodBenchmarks()
{
  registerWorkloads<DIM, U, METHOD, Encode>();
  registerWorkloads<DIM, U, METHOD, Decode>();
  registerWorkloads<DIM, U, METHOD, Mixed>();
}

template <int DIM, typename U>
void registerBenchmarksForShape()
{
  registerMethodBenchmarks<DIM, U, Public>();
  registerMethodBenchmarks<DIM, U, Scalar>();
}
}  // namespace

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;
  bool verify_only = false;
  axom::CLI::App app {"Axom Hilbert codec benchmarks"};
  app.add_flag("--verify-only",
               verify_only,
               "Check both codecs at every valid depth without timing");
  app.allow_extras();
  CLI11_PARSE(app, argc, argv);

  // Forward unconsumed arguments to Google Benchmark in their original order.
  auto benchmark_args = app.remaining();
  benchmark_args.insert(benchmark_args.begin(), argv[0]);
  std::vector<char*> benchmark_argv;
  for(auto& arg : benchmark_args)
  {
    benchmark_argv.push_back(arg.data());
  }
  int benchmark_argc = static_cast<int>(benchmark_argv.size());
  benchmark_argv.push_back(nullptr);
  benchmark::Initialize(&benchmark_argc, benchmark_argv.data());
  if(benchmark::ReportUnrecognizedArguments(benchmark_argc, benchmark_argv.data()))
  {
    benchmark::Shutdown();
    return 1;
  }

  if(verify_only)
  {
    const bool valid = verifyCodecsForShape<2, std::uint32_t>() &&
      verifyCodecsForShape<2, std::uint64_t>() && verifyCodecsForShape<3, std::uint32_t>() &&
      verifyCodecsForShape<3, std::uint64_t>() && verifyCodecsForShape<4, std::uint32_t>() &&
      verifyCodecsForShape<4, std::uint64_t>();
    benchmark::Shutdown();
    return valid ? 0 : 1;
  }

  registerBenchmarksForShape<2, std::uint32_t>();
  registerBenchmarksForShape<2, std::uint64_t>();
  registerBenchmarksForShape<3, std::uint32_t>();
  registerBenchmarksForShape<3, std::uint64_t>();
  registerBenchmarksForShape<4, std::uint32_t>();
  registerBenchmarksForShape<4, std::uint64_t>();
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
