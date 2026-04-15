// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "benchmark/benchmark.h"

#include "axom/core/numerics/Determinants.hpp"
#include "axom/primal.hpp"

#include <array>
#include <cstdint>

namespace
{
using axom::numerics::determinant;

template <int DIM>
using PointType = axom::primal::Point<double, DIM>;

template <int DIM>
using VectorType = axom::primal::Vector<double, DIM>;

template <int DIM>
struct CircumsphereEval
{
  PointType<DIM> center {};
  double radius_sq {0.};

  CircumsphereEval() = default;

  CircumsphereEval(const PointType<DIM>& origin, const VectorType<DIM>& center_offset)
    : center(origin + center_offset)
    , radius_sq(center_offset.squared_norm())
  { }
};

template <int DIM>
using SimplexType = std::array<PointType<DIM>, DIM + 1>;

template <int DIM>
using SampleSet = std::array<SimplexType<DIM>, 4096>;

double unitInterval(std::uint64_t bits)
{
  return static_cast<double>(bits & 0xFFFFFFFFu) / static_cast<double>(0x100000000ULL);
}

std::uint64_t stepLcg(std::uint64_t state)
{
  return state * 6364136223846793005ULL + 1442695040888963407ULL;
}

template <int DIM>
SampleSet<DIM> makeSamples();

template <>
SampleSet<2> makeSamples<2>()
{
  SampleSet<2> samples;
  std::uint64_t state = 0x1234abcdu;

  for(auto& simplex : samples)
  {
    state = stepLcg(state);
    const double bx = 0.8 * unitInterval(state);
    state = stepLcg(state);
    const double by = 0.8 * unitInterval(state);
    state = stepLcg(state);
    const double s0 = 0.01 + 0.04 * unitInterval(state);
    state = stepLcg(state);
    const double s1 = 0.01 + 0.04 * unitInterval(state);
    state = stepLcg(state);
    const double t0 = 0.01 + 0.04 * unitInterval(state);
    state = stepLcg(state);
    const double t1 = 0.01 + 0.04 * unitInterval(state);

    simplex[0] = PointType<2> {bx, by};
    simplex[1] = PointType<2> {bx + s0, by + 0.15 * s1};
    simplex[2] = PointType<2> {bx + 0.2 * t0, by + t1};
  }

  return samples;
}

template <>
SampleSet<3> makeSamples<3>()
{
  SampleSet<3> samples;
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;

  for(auto& simplex : samples)
  {
    state = stepLcg(state);
    const double bx = 0.75 * unitInterval(state);
    state = stepLcg(state);
    const double by = 0.75 * unitInterval(state);
    state = stepLcg(state);
    const double bz = 0.75 * unitInterval(state);

    state = stepLcg(state);
    const double a0 = 0.01 + 0.03 * unitInterval(state);
    state = stepLcg(state);
    const double a1 = 0.01 + 0.03 * unitInterval(state);
    state = stepLcg(state);
    const double b0 = 0.01 + 0.03 * unitInterval(state);
    state = stepLcg(state);
    const double b1 = 0.01 + 0.03 * unitInterval(state);
    state = stepLcg(state);
    const double c0 = 0.01 + 0.03 * unitInterval(state);
    state = stepLcg(state);
    const double c1 = 0.01 + 0.03 * unitInterval(state);

    simplex[0] = PointType<3> {bx, by, bz};
    simplex[1] = PointType<3> {bx + a0, by + 0.1 * a1, bz + 0.05 * a1};
    simplex[2] = PointType<3> {bx + 0.15 * b0, by + b1, bz + 0.12 * b0};
    simplex[3] = PointType<3> {bx + 0.08 * c0, by + 0.18 * c0, bz + c1};
  }

  return samples;
}

template <int DIM>
inline double consumeEval(const CircumsphereEval<DIM>& eval)
{
  double value = eval.radius_sq;
  for(int dim = 0; dim < DIM; ++dim)
  {
    value += eval.center[dim];
  }
  return value;
}

inline CircumsphereEval<2> evaluateCircumsphereScalarEdges(const SimplexType<2>& simplex)
{
  const PointType<2>& p0 = simplex[0];
  const PointType<2>& p1 = simplex[1];
  const PointType<2>& p2 = simplex[2];

  const double vx0 = p1[0] - p0[0];
  const double vx1 = p2[0] - p0[0];
  const double vy0 = p1[1] - p0[1];
  const double vy1 = p2[1] - p0[1];

  const double sq0 = vx0 * vx0 + vy0 * vy0;
  const double sq1 = vx1 * vx1 + vy1 * vy1;

  const double a = determinant(vx0, vx1, vy0, vy1);
  const double eps = (a >= 0.) ? axom::primal::PRIMAL_TINY : -axom::primal::PRIMAL_TINY;
  const double ood = 1. / (2. * a + eps);

  const double center_offset_x = determinant(sq0, sq1, vy0, vy1) * ood;
  const double center_offset_y = -determinant(sq0, sq1, vx0, vx1) * ood;

  return CircumsphereEval<2>(p0, VectorType<2> {center_offset_x, center_offset_y});
}

inline CircumsphereEval<2> evaluateCircumsphereVectorEdges(const SimplexType<2>& simplex)
{
  const PointType<2>& p0 = simplex[0];
  const PointType<2>& p1 = simplex[1];
  const PointType<2>& p2 = simplex[2];

  const VectorType<2> v0(p0, p1);
  const VectorType<2> v1(p0, p2);

  const double sq0 = v0.squared_norm();
  const double sq1 = v1.squared_norm();

  const double a = determinant(v0[0], v1[0], v0[1], v1[1]);
  const double eps = (a >= 0.) ? axom::primal::PRIMAL_TINY : -axom::primal::PRIMAL_TINY;
  const double ood = 1. / (2. * a + eps);

  const double center_offset_x = determinant(sq0, sq1, v0[1], v1[1]) * ood;
  const double center_offset_y = -determinant(sq0, sq1, v0[0], v1[0]) * ood;

  return CircumsphereEval<2>(p0, VectorType<2> {center_offset_x, center_offset_y});
}

inline CircumsphereEval<3> evaluateCircumsphereScalarEdges(const SimplexType<3>& simplex)
{
  const PointType<3>& p0 = simplex[0];
  const PointType<3>& p1 = simplex[1];
  const PointType<3>& p2 = simplex[2];
  const PointType<3>& p3 = simplex[3];

  const double vx0 = p1[0] - p0[0];
  const double vx1 = p2[0] - p0[0];
  const double vx2 = p3[0] - p0[0];
  const double vy0 = p1[1] - p0[1];
  const double vy1 = p2[1] - p0[1];
  const double vy2 = p3[1] - p0[1];
  const double vz0 = p1[2] - p0[2];
  const double vz1 = p2[2] - p0[2];
  const double vz2 = p3[2] - p0[2];

  const double sq0 = vx0 * vx0 + vy0 * vy0 + vz0 * vz0;
  const double sq1 = vx1 * vx1 + vy1 * vy1 + vz1 * vz1;
  const double sq2 = vx2 * vx2 + vy2 * vy2 + vz2 * vz2;

  const double a = determinant(vx0, vx1, vx2, vy0, vy1, vy2, vz0, vz1, vz2);
  const double eps = (a >= 0.) ? axom::primal::PRIMAL_TINY : -axom::primal::PRIMAL_TINY;
  const double ood = 1. / (2. * a + eps);

  const double center_offset_x = determinant(sq0, sq1, sq2, vy0, vy1, vy2, vz0, vz1, vz2) * ood;
  const double center_offset_y = determinant(sq0, sq1, sq2, vz0, vz1, vz2, vx0, vx1, vx2) * ood;
  const double center_offset_z = determinant(sq0, sq1, sq2, vx0, vx1, vx2, vy0, vy1, vy2) * ood;

  return CircumsphereEval<3>(p0, VectorType<3> {center_offset_x, center_offset_y, center_offset_z});
}

inline CircumsphereEval<3> evaluateCircumsphereVectorEdges(const SimplexType<3>& simplex)
{
  const PointType<3>& p0 = simplex[0];
  const PointType<3>& p1 = simplex[1];
  const PointType<3>& p2 = simplex[2];
  const PointType<3>& p3 = simplex[3];

  const VectorType<3> v0(p0, p1);
  const VectorType<3> v1(p0, p2);
  const VectorType<3> v2(p0, p3);

  const double sq0 = v0.squared_norm();
  const double sq1 = v1.squared_norm();
  const double sq2 = v2.squared_norm();

  const double a = determinant(v0[0], v1[0], v2[0], v0[1], v1[1], v2[1], v0[2], v1[2], v2[2]);
  const double eps = (a >= 0.) ? axom::primal::PRIMAL_TINY : -axom::primal::PRIMAL_TINY;
  const double ood = 1. / (2. * a + eps);

  const double center_offset_x =
    determinant(sq0, sq1, sq2, v0[1], v1[1], v2[1], v0[2], v1[2], v2[2]) * ood;
  const double center_offset_y =
    determinant(sq0, sq1, sq2, v0[2], v1[2], v2[2], v0[0], v1[0], v2[0]) * ood;
  const double center_offset_z =
    determinant(sq0, sq1, sq2, v0[0], v1[0], v2[0], v0[1], v1[1], v2[1]) * ood;

  return CircumsphereEval<3>(p0, VectorType<3> {center_offset_x, center_offset_y, center_offset_z});
}

template <int DIM, typename Kernel>
void runCircumsphereBenchmark(benchmark::State& state, Kernel&& kernel)
{
  const auto& samples = []() -> const SampleSet<DIM>& {
    static const SampleSet<DIM> value = makeSamples<DIM>();
    return value;
  }();

  const std::size_t mask = samples.size() - 1;
  std::size_t idx = 0;
  double checksum = 0.;

  for(auto _ : state)
  {
    const auto& simplex = samples[idx];
    checksum += consumeEval(kernel(simplex));
    idx = (idx + 1) & mask;
  }

  benchmark::DoNotOptimize(checksum);
  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
}

void benchmark_scalar_edges_2d(benchmark::State& state)
{
  runCircumsphereBenchmark<2>(state, [](const SimplexType<2>& simplex) {
    return evaluateCircumsphereScalarEdges(simplex);
  });
}

void benchmark_vector_edges_2d(benchmark::State& state)
{
  runCircumsphereBenchmark<2>(state, [](const SimplexType<2>& simplex) {
    return evaluateCircumsphereVectorEdges(simplex);
  });
}

void benchmark_scalar_edges_3d(benchmark::State& state)
{
  runCircumsphereBenchmark<3>(state, [](const SimplexType<3>& simplex) {
    return evaluateCircumsphereScalarEdges(simplex);
  });
}

void benchmark_vector_edges_3d(benchmark::State& state)
{
  runCircumsphereBenchmark<3>(state, [](const SimplexType<3>& simplex) {
    return evaluateCircumsphereVectorEdges(simplex);
  });
}

}  // namespace

BENCHMARK(benchmark_scalar_edges_2d);
BENCHMARK(benchmark_vector_edges_2d);
BENCHMARK(benchmark_scalar_edges_3d);
BENCHMARK(benchmark_vector_edges_3d);

BENCHMARK_MAIN();
