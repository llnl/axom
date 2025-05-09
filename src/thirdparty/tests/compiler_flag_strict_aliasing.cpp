// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * This file tests the ability to disable warnings about strict aliasing
 * on all supported compilers using the AXOM_DISABLE_ALIASING_WARNINGS
 * build variable.
 *
 * The example is modeled after a known example in python's headers
 * (see http://legacy.python.org/dev/peps/pep-3123) and might need
 * to be built with optimization enabled to trigger the warning.
 */

#include <iostream>

struct Foo
{
  int i;
  Foo* ob;
};
struct Bar
{
  int i;
  Foo* ob;
};

int main()
{
  Foo foo = {1, nullptr};
  ((Bar*)(&foo))->i++;  // violates strict aliasing

  std::cout << " foo.i: " << foo.i << std::endl;
  return 0;
}
