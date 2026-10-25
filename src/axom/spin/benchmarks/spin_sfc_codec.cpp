// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_sfc_codec.cpp
 * \brief Compares Hilbert and Morton space filling curve conversions
 * in 2D, 3D, and 4D.
 *
 * Use Google Benchmark's --benchmark_filter and --benchmark_repetitions options
 * for timing in Release. Run --verify-only to check the adapters.
 *
 * Each timed iteration processes a batch. Divide its time by
 * `batch * operations_per_item` to get the average time per conversion.
 * Mixed cases encode a point and decode an independent key for each item.
 * Cases identify both the curve and its implementation. All use the same inputs.
 * Morton has no depth argument, so Variable changes its input ranges while
 * Hilbert also receives each item's depth.
 */

#include "axom/config.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/slic.hpp"
#include "axom/spin/HilbertIndex.hpp"
#include "axom/spin/MortonIndex.hpp"
#include "axom/spin/SpaceFillingCurve.hpp"
#include "axom/spin/internal/HilbertFSM.hpp"

#include "benchmark/benchmark.h"

#include <algorithm>
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

using Curve = axom::spin::SpaceFillingCurve;

// Public calls the selected curve's API. Scalar calls the reference Hilbert transform.
// FSM uses the production tables.
// The number in each adapter name gives bits per coordinate per lookup.
// All adapters allow normal inlining.
enum Implementation
{
  Public,
  Scalar,
  FSM1,
  FSM2,
  FSM3
};

const char* curveName(Curve curve) { return curve == Curve::Hilbert ? "Hilbert" : "Morton"; }

const char* implementationName(Implementation impl)
{
  constexpr const char* names[] = {"Public", "Scalar", "FSM1", "FSM2", "FSM3"};
  return names[impl];
}

// Setup bounds inputs to the requested depth, so Morton can ignore bits and
// call its native API without adding masking work to the timed conversions.
template <int DIM, typename U, Curve CURVE, Implementation IMPL>
struct Codec
{
  static_assert(CURVE == Curve::Hilbert || IMPL == Public,
                "Morton has only a public benchmark adapter");
  using Point = axom::primal::Point<std::uint32_t, DIM>;
  static constexpr int MAX_BITS =
    std::min(std::numeric_limits<std::uint32_t>::digits, std::numeric_limits<U>::digits / DIM);

  using HilbertPublicTransform = axom::spin::Hilbertizer<std::uint32_t, U, DIM>;
  using HilbertScalarTransform =
    axom::spin::internal::hilbert::HilbertTransform<std::uint32_t, U, DIM>;
  using MortonPublicTransform = axom::spin::Mortonizer<std::uint32_t, U, DIM>;
  static constexpr int LEVELS = IMPL == FSM3 ? 3 : IMPL == FSM2 ? 2 : 1;
  using HilbertFSMTransform =
    axom::spin::internal::hilbert::HilbertFSMTransform<std::uint32_t, U, DIM, LEVELS>;

  static_assert(MAX_BITS == HilbertPublicTransform::MAX_UNIQUE_BITS &&
                  MAX_BITS == MortonPublicTransform::MAX_UNIQUE_BITS,
                "Codec depths must match the shared input domain");

  static U encode(const Point& point, int bits)
  {
    if constexpr(CURVE == Curve::Morton)
    {
      AXOM_UNUSED_VAR(bits);
      return MortonPublicTransform::mortonize(point);
    }
    else if constexpr(IMPL == Public)
    {
      return HilbertPublicTransform::toIndex(point, bits);
    }
    else if constexpr(IMPL == Scalar)
    {
      return HilbertScalarTransform::toIndex(point, bits);
    }
    else
    {
      return HilbertFSMTransform::toIndex(point, bits);
    }
  }

  static Point decode(U key, int bits)
  {
    if constexpr(CURVE == Curve::Morton)
    {
      AXOM_UNUSED_VAR(bits);
      return MortonPublicTransform::demortonize(key);
    }
    else if constexpr(IMPL == Public)
    {
      return HilbertPublicTransform::toPoint(key, bits);
    }
    else if constexpr(IMPL == Scalar)
    {
      return HilbertScalarTransform::toPoint(key, bits);
    }
    else
    {
      return HilbertFSMTransform::toPoint(key, bits);
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

// Check Hilbert against the scalar transform and Morton against bit interleaving.
template <int DIM, typename U, Curve CURVE, Implementation IMPL>
bool verifyCodec()
{
  using CodecType = Codec<DIM, U, CURVE, IMPL>;
  constexpr int max_bits = CodecType::MAX_BITS;
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
      U expected_key = 0;
      typename CodecType::Point expected_point;
      if constexpr(CURVE == Curve::Morton)
      {
        // Interleave each bit plane explicitly, with x in its lowest bit.
        // This reference is independent of Mortonizer's expand/contract helpers.
        for(int bit = 0; bit < bits; ++bit)
        {
          for(int axis = 0; axis < DIM; ++axis)
          {
            const int key_bit = bit * DIM + axis;
            expected_key |= static_cast<U>((point[axis] >> bit) & 1) << key_bit;
            expected_point[axis] |= static_cast<std::uint32_t>((key >> key_bit) & U {1}) << bit;
          }
        }
      }
      else
      {
        using HilbertScalarTransform =
          axom::spin::internal::hilbert::HilbertTransform<std::uint32_t, U, DIM>;
        expected_key = HilbertScalarTransform::toIndex(point, bits);
        expected_point = HilbertScalarTransform::toPoint(key, bits);
      }
      bool valid_point = true;
      for(int axis = 0; axis < DIM; ++axis)
      {
        valid_point = valid_point && decoded[axis] <= coord_mask;
      }
      if(encoded > key_mask || !valid_point || encoded != expected_key || decoded != expected_point ||
         CodecType::encode(decoded, bits) != key || CodecType::decode(encoded, bits) != point)
      {
        SLIC_WARNING(
          axom::fmt::format("Verification failed: {}/{} D={} key bits={} depth={} sample={}",
                            curveName(CURVE),
                            implementationName(IMPL),
                            DIM,
                            sizeof(U) * 8,
                            bits,
                            i));
        return false;
      }
    }
  }
  SLIC_INFO(axom::fmt::format("Verified {}/{}/{}D/u{} at every depth",
                              curveName(CURVE),
                              implementationName(IMPL),
                              DIM,
                              sizeof(U) * 8));
  return true;
}

template <int DIM, typename U>
bool verifyCodecsForShape()
{
  bool valid = verifyCodec<DIM, U, Curve::Hilbert, Public>() &&
    verifyCodec<DIM, U, Curve::Hilbert, Scalar>() && verifyCodec<DIM, U, Curve::Morton, Public>() &&
    verifyCodec<DIM, U, Curve::Hilbert, FSM1>();
  if constexpr(DIM <= 3)
  {
    valid = verifyCodec<DIM, U, Curve::Hilbert, FSM2>() && valid;
  }
  if constexpr(DIM == 2)
  {
    valid = verifyCodec<DIM, U, Curve::Hilbert, FSM3>() && valid;
  }
  return valid;
}

template <int DIM, typename U, Curve CURVE, Implementation IMPL, Direction DIR, Workload WORK, Dataset DATA>
void benchmarkCodec(benchmark::State& state)
{
  using CodecType = Codec<DIM, U, CURVE, IMPL>;
  using Point = typename CodecType::Point;
  constexpr int max_bits = CodecType::MAX_BITS;
  const int count = WORK == Small || WORK == Pressure ? SMALL_BATCH_SIZE : LARGE_BATCH_SIZE;
  // Allocation and input generation precede the state loop, so neither is timed.
  std::vector<Point> points(count), decoded(count);
  std::vector<U> keys(count), encoded(count);
  std::vector<int> depths(count);
  // Vary the sequence by dimension. All codecs use the same seed and inputs.
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
      // Fixed-depth cases let the compiler specialize Hilbert. Variable supplies
      // a runtime depth per item; Morton uses the same inputs without this argument.
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

template <int DIM, typename U, Curve CURVE, Implementation IMPL, Direction DIR, Workload WORK, Dataset DATA>
void registerBenchmark()
{
  const char* directions[] = {"Encode", "Decode", "Mixed"};
  const char* workloads[] = {"Warm", "Small", "Pressure", "Variable"};
  const char* datasets[] = {"Uniform", "Sequential", "Clustered"};
  const auto label = axom::fmt::format("{}/{}/{}D/u{}/{}/{}/{}",
                                       curveName(CURVE),
                                       implementationName(IMPL),
                                       DIM,
                                       sizeof(U) * 8,
                                       directions[DIR],
                                       workloads[WORK],
                                       datasets[DATA]);
  auto* b = benchmark::RegisterBenchmark(label.c_str(),
                                         &benchmarkCodec<DIM, U, CURVE, IMPL, DIR, WORK, DATA>);
  b->UseRealTime()->Unit(benchmark::kNanosecond);
  // Limit the untimed buffer walks per repetition.
  if constexpr(WORK == Pressure)
  {
    b->Iterations(CACHE_PRESSURE_ITERATIONS);
  }
}

template <int DIM, typename U, Curve CURVE, Implementation IMPL, Direction DIR>
void registerWorkloads()
{
  registerBenchmark<DIM, U, CURVE, IMPL, DIR, Warm, Uniform>();
  registerBenchmark<DIM, U, CURVE, IMPL, DIR, Small, Uniform>();
  registerBenchmark<DIM, U, CURVE, IMPL, DIR, Pressure, Uniform>();
  registerBenchmark<DIM, U, CURVE, IMPL, DIR, Variable, Uniform>();
  registerBenchmark<DIM, U, CURVE, IMPL, DIR, Warm, Sequential>();
  registerBenchmark<DIM, U, CURVE, IMPL, DIR, Warm, Clustered>();
}

template <int DIM, typename U, Curve CURVE, Implementation IMPL>
void registerCodecBenchmarks()
{
  registerWorkloads<DIM, U, CURVE, IMPL, Encode>();
  registerWorkloads<DIM, U, CURVE, IMPL, Decode>();
  registerWorkloads<DIM, U, CURVE, IMPL, Mixed>();
}

template <int DIM, typename U>
void registerBenchmarksForShape()
{
  registerCodecBenchmarks<DIM, U, Curve::Hilbert, Public>();
  registerCodecBenchmarks<DIM, U, Curve::Hilbert, Scalar>();
  registerCodecBenchmarks<DIM, U, Curve::Morton, Public>();
  registerCodecBenchmarks<DIM, U, Curve::Hilbert, FSM1>();
  if constexpr(DIM <= 3)
  {
    registerCodecBenchmarks<DIM, U, Curve::Hilbert, FSM2>();
  }
  if constexpr(DIM == 2)
  {
    registerCodecBenchmarks<DIM, U, Curve::Hilbert, FSM3>();
  }
}
}  // namespace

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;
  bool verify_only = false;
  axom::CLI::App app {"Axom space-filling curve codec benchmarks"};
  app.add_flag("--verify-only", verify_only, "Check all codecs at every valid depth without timing");
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
