// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_static_asserts.cpp
 *
 * \brief Compile-time conformance harness for slam's constexpr arithmetic
 *
 * The tests in this file exercise the compile time set/value asserts
 * for sizes, offsets, strides, modular wraparound, and policy composition, so the checks cost
 * nothing at run time and run on every backend on every build. 
 * 
 * Many of the checks in this file are static_asserts. If the file compiles, its invariants hold.
 */

#include "gtest/gtest.h"

#include "axom/slam/ModularInt.hpp"
#include "axom/slam/Optional.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/ValuePolicies.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

//------------------------------------------------------------------------------
// Value policies: compile-time-known values report their value, validity,
// and (for size) emptiness without any runtime state.
//------------------------------------------------------------------------------

// --- Size ---
using Size0 = policies::CompileTimeSize<int, 0>;
using Size5 = policies::CompileTimeSize<int, 5>;

static_assert(Size5().size() == 5, "CompileTimeSize reports its size");
static_assert(Size0().size() == 0, "CompileTimeSize of zero reports zero");
static_assert(!Size5().empty(), "non-zero CompileTimeSize is not empty");
static_assert(Size0().empty(), "zero CompileTimeSize is empty");
static_assert(Size5().isValid(false), "non-negative size is valid");
static_assert(Size5().DEFAULT_VALUE == 5, "DEFAULT_VALUE matches the NTTP");

// --- Offset ---
using Off0 = policies::CompileTimeOffset<int, 0>;
using Off3 = policies::CompileTimeOffset<int, 3>;

static_assert(Off3().offset() == 3, "CompileTimeOffset reports its offset");
static_assert(Off0().offset() == 0, "ZeroOffset-equivalent reports zero");
static_assert(policies::ZeroOffset<int>().offset() == 0, "ZeroOffset is zero");
static_assert(Off3().isValid(false), "all offsets are valid");

// --- Stride ---
using Stride1 = policies::CompileTimeStride<int, 1>;
using Stride4 = policies::CompileTimeStride<int, 4>;

static_assert(Stride4().stride() == 4, "CompileTimeStride reports its stride");
static_assert(Stride4().shape() == 4, "CompileTimeStride shape equals stride");
static_assert(policies::StrideOne<int>().stride() == 1, "StrideOne is one");
static_assert(Stride4().isValid(false), "non-zero stride is valid");
static_assert(Stride1::IS_COMPILE_TIME, "CompileTimeStride is compile-time");
static_assert(Stride1::NumDims == 1, "scalar stride is one-dimensional");

//------------------------------------------------------------------------------
// The unified core: RuntimeValue / CompileTimeValue behave through the tags.
//------------------------------------------------------------------------------
static_assert(policies::CompileTimeValue<7, policies::SizeTag<int>>().value() == 7,
              "CompileTimeValue forwards its NTTP");
static_assert(policies::SizeTag<int>::isValidValue(0), "size 0 is valid");
static_assert(!policies::SizeTag<int>::isValidValue(-1), "negative size is invalid");
static_assert(!policies::StrideTag<int>::isValidValue(0), "stride 0 is invalid");
static_assert(policies::StrideTag<int>::isValidValue(-2), "negative stride is valid");
static_assert(policies::OffsetTag<int>::isValidValue(-5), "any offset is valid");
static_assert(policies::StrideTag<int>::defaultValue() == 1, "default stride is one");
static_assert(policies::SizeTag<int>::defaultValue() == 0, "default size is zero");

//------------------------------------------------------------------------------
// Offset / stride composition: the flat index of element i in a strided,
// offset range is  offset + i * stride.  Verify the policy values compose.
//------------------------------------------------------------------------------
constexpr int flatIndex(int offset, int stride, int i) { return offset + i * stride; }

static_assert(flatIndex(Off3().offset(), Stride4().stride(), 0) == 3,
              "element 0 lands at the offset");
static_assert(flatIndex(Off3().offset(), Stride4().stride(), 2) == 11,
              "offset 3 + 2*stride 4 == 11");
static_assert(flatIndex(policies::ZeroOffset<int>().offset(), policies::StrideOne<int>().stride(), 9) ==
                9,
              "unit stride, zero offset is the identity index");

//------------------------------------------------------------------------------
// ModularInt: cyclic arithmetic is fully constexpr-evaluable.
//------------------------------------------------------------------------------
using Mod5 = slam::ModularInt<policies::CompileTimeSize<int, 5>>;

static_assert(int(Mod5(0)) == 0, "0 mod 5 == 0");
static_assert(int(Mod5(5)) == 0, "5 mod 5 == 0 (wraparound)");
static_assert(int(Mod5(7)) == 2, "7 mod 5 == 2");
static_assert(int(Mod5(-1)) == 4, "-1 mod 5 normalizes to 4");
static_assert(int(Mod5(-5)) == 0, "-5 mod 5 == 0");
static_assert(int(Mod5(13)) == 3, "13 mod 5 == 3");

// arithmetic operators
static_assert(int(Mod5(3) + 4) == 2, "(3+4) mod 5 == 2");
static_assert(int(Mod5(4) + 1) == 0, "(4+1) mod 5 wraps to 0");
static_assert(int(4 + Mod5(3)) == 2, "int + ModularInt commutes");
static_assert(int(Mod5(1) - 3) == 3, "(1-3) mod 5 == 3");
static_assert(int(Mod5(2) * 4) == 3, "(2*4) mod 5 == 3");

// pre/post increment in a constexpr lambda
constexpr int incTwice(int start)
{
  Mod5 m(start);
  ++m;
  ++m;
  return int(m);
}
static_assert(incTwice(4) == 1, "++ twice from 4 mod 5 == 1");

// equality across the modulus
static_assert(Mod5(2) == Mod5(7), "2 and 7 are equal mod 5");
static_assert(Mod5(2) != Mod5(3), "2 and 3 differ mod 5");

//------------------------------------------------------------------------------
// slam::Optional: a device-safe "maybe a value", fully constexpr.
//------------------------------------------------------------------------------
static_assert(!slam::Optional<int>().has_value(), "default Optional is disengaged");
static_assert(slam::Optional<int>(42).has_value(), "value-constructed Optional is engaged");
static_assert(*slam::Optional<int>(42) == 42, "engaged Optional yields its value");
static_assert(slam::Optional<int>(42).value_or(-1) == 42, "value_or returns the value when engaged");
static_assert(slam::Optional<int>().value_or(-1) == -1, "value_or returns fallback when empty");
static_assert(static_cast<bool>(slam::Optional<int>(0)), "engaged-with-zero is still engaged");
static_assert(!static_cast<bool>(slam::Optional<int>()), "disengaged converts to false");
// Trivial copyability allows it to be captured on device
static_assert(std::is_trivially_copyable_v<slam::Optional<int>>,
              "slam::Optional<int> is trivially copyable (device-capturable)");

}  // anonymous namespace

//------------------------------------------------------------------------------
TEST(slam_static_asserts, compile_time_value_holds)
{
  // Checks that a constexpr value survives to run time.

  constexpr int wrapped = int(Mod5(7));
  EXPECT_EQ(wrapped, 2);

  constexpr int idx = flatIndex(Off3().offset(), Stride4().stride(), 2);
  EXPECT_EQ(idx, 11);
}

TEST(slam_optional, engaged_and_disengaged)
{
  // Checks runtime API of slam::Optional

  axom::slam::Optional<double> empty;
  EXPECT_FALSE(empty.has_value());
  EXPECT_FALSE(static_cast<bool>(empty));
  EXPECT_DOUBLE_EQ(empty.value_or(2.5), 2.5);

  axom::slam::Optional<double> full(3.25);
  EXPECT_TRUE(full.has_value());
  EXPECT_TRUE(static_cast<bool>(full));
  EXPECT_DOUBLE_EQ(*full, 3.25);
  EXPECT_DOUBLE_EQ(full.value_or(2.5), 3.25);

  *full = 9.0;
  EXPECT_DOUBLE_EQ(*full, 9.0);
}
