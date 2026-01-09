// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic/interface/slic.hpp"

namespace
{
void custom_abort_function() { }
}  // namespace

TEST(slic_scoped_abort, throws_on_error)
{
  axom::slic::initialize();

  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW({ SLIC_ERROR("testing abort-to-throw"); }, axom::slic::SlicAbortException);

  axom::slic::finalize();
}

TEST(slic_scoped_abort, restores_state_on_stack_unwind)
{
  axom::slic::initialize();

  axom::slic::disableAbortOnError();
  axom::slic::setAbortFunction(custom_abort_function);

  const bool prev_abort_on_error = axom::slic::isAbortOnErrorsEnabled();
  const auto prev_abort_function = axom::slic::getAbortFunction();

  try
  {
    axom::slic::ScopedAbortToThrow abort_guard;
    SLIC_ERROR("testing state restoration");
    FAIL() << "Expected SlicAbortException";
  }
  catch(const axom::slic::SlicAbortException&)
  { }

  EXPECT_EQ(prev_abort_on_error, axom::slic::isAbortOnErrorsEnabled());
  EXPECT_EQ(prev_abort_function, axom::slic::getAbortFunction());

  axom::slic::finalize();
}
