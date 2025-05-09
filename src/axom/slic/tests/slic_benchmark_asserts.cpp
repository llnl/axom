// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "benchmark/benchmark.h"
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

#include <cstdint>
#include <iostream>

/*!
 * \file
 *
 * A series of tests of SLIC logging macros (using SLIC_WARNING to avoid 
 * killing the program) in the TPL benchmark framework.
 *
 * These tests were written to see why SLIC macros was not being logged
 * when called from a benchmark fixture's destructor.
 * It turns out that the fixtures are destructed after main exits,
 * so SLIC is already finalized.
 */

namespace
{
static constexpr bool IS_VALID = true;
static constexpr bool SHOULD_PRINT = false;

void printMsg(std::string const& str, bool always_print = false)
{
  if(SHOULD_PRINT || always_print)
  {
    std::cout << str << std::endl;
  }
}

struct AssertCtor
{
  AssertCtor() { SLIC_WARNING_IF(IS_VALID, "Testing warning in .ctor"); };

  struct AssertMethod
  {
    void foo() { SLIC_WARNING_IF(IS_VALID, "Testing warning in class method"); }
  };
};

struct AssertDtor
{
  ~AssertDtor() { SLIC_WARNING_IF(IS_VALID, "Testing warning in .dtor"); }
};

class SetFixtureC : public ::benchmark::Fixture
{
public:
  SetFixtureC()
  {
    static bool override_once = true;
    if(override_once)
    {
      std::cout << "*** Note: cannot use slic macros in fixture constructors since "
                   "they're called before main (i.e. before slic::initialize()) ***"
                << std::endl;
      override_once = false;
    }

    //SLIC_WARNING_IF(IS_VALID, "Testing warning in .ctor");
  }
};

class SetFixtureS : public ::benchmark::Fixture
{
public:
  void SetUp(const ::benchmark::State& /*state*/) override
  {
    SLIC_WARNING_IF(IS_VALID, "Testing warning in setup");
  }
};

class SetFixtureT : public ::benchmark::Fixture
{
public:
  void TearDown(const ::benchmark::State& /*state*/) override
  {
    SLIC_WARNING_IF(IS_VALID, "Testing warning in teardown");
  }
};

class SetFixtureD : public ::benchmark::Fixture
{
public:
  ~SetFixtureD()
  {
    static bool override_once = true;
    if(override_once)
    {
      std::cout << "*** Note: cannot use slic macros in fixture destructors since "
                   "they're called after main (i.e. after slic::finalize()) ***"
                << std::endl;
      override_once = false;
    }

    //SLIC_WARNING_IF(IS_VALID, "Testing warning in .dtor");
  }
};

class SetFixtureOutput : public ::benchmark::Fixture
{
public:
  SetFixtureOutput() { printMsg("  fixure .ctor (called before main)"); }
  void SetUp(const ::benchmark::State& /*state*/) override { printMsg("  fixure SetUp"); }
  void TearDown(const ::benchmark::State& /*state*/) override { printMsg("  fixure TearDown"); }
  ~SetFixtureOutput() { printMsg("  fixure ~.dtor (called after main)"); }
};

static void BM_AssertBefore(benchmark::State& state)
{
  SLIC_WARNING_IF(IS_VALID, "Testing warning before running benchmark");
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
}

static void BM_AssertAfter(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
  SLIC_WARNING_IF(IS_VALID, "Testing warning after running benchmark");
}

static void BM_CallAssertCtorBefore(benchmark::State& state)
{
  AssertCtor();
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
}

static void BM_CallAssertDtorBefore(benchmark::State& state)
{
  AssertDtor();
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
}

static void BM_CallAssertDtorDuring(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    AssertDtor();
  }
}

}  // namespace

BENCHMARK(BM_AssertBefore)->Iterations(1);
BENCHMARK(BM_AssertAfter)->Iterations(1);
BENCHMARK(BM_CallAssertCtorBefore)->Iterations(1);
BENCHMARK(BM_CallAssertDtorBefore)->Iterations(1);
BENCHMARK(BM_CallAssertDtorDuring)->Iterations(1);

BENCHMARK_F(SetFixtureC, assertCtor)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
}
BENCHMARK_REGISTER_F(SetFixtureC, assertCtor)->Iterations(1);

BENCHMARK_F(SetFixtureS, assertSetup)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
}
BENCHMARK_REGISTER_F(SetFixtureS, assertSetup)->Iterations(1);

BENCHMARK_F(SetFixtureT, assertTeardown)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
}
BENCHMARK_REGISTER_F(SetFixtureT, assertTeardown)->Iterations(1);

BENCHMARK_F(SetFixtureD, assertDtor)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
}
BENCHMARK_REGISTER_F(SetFixtureD, assertDtor)->Iterations(1);

// The following two benchmarks are here to show the order of construction,
// setup, teardown and destruction of benchmarks fixtures.
// Note: Need to enable printing (by setting SHOULD_PRINT above to true) to see this.
BENCHMARK_F(SetFixtureOutput, process1)(benchmark::State& state)
{
  printMsg("   P1 (before)");
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
  printMsg("   P1 (after)");
}
BENCHMARK_REGISTER_F(SetFixtureOutput, process1)->Iterations(1);

BENCHMARK_F(SetFixtureOutput, process2)(benchmark::State& state)
{
  printMsg("   P2 (before)");
  while(state.KeepRunning())
  {
    std::int64_t x = 0;
    benchmark::DoNotOptimize(x);
  }
  printMsg("   P2 (after)");
}
BENCHMARK_REGISTER_F(SetFixtureOutput, process2)->Iterations(1);

int main(int argc, char* argv[])
{
  printMsg("at start of main");

  printMsg(" Before init logger");
  axom::slic::SimpleLogger logger;
  printMsg(" After init logger");

  printMsg(" Before init benchmark");
  ::benchmark::Initialize(&argc, argv);
  printMsg(" after init benchmark");

  printMsg(" before run benchmark");
  ::benchmark::RunSpecifiedBenchmarks();
  printMsg(" after run benchmark");

  printMsg("about to return from main");
  return 0;
}
