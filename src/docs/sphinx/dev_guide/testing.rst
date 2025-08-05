.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _testing-label:

****************************************
Tests and Examples
****************************************

This section describes how to build and organize tests and examples
for Axom components. These live in ``tests`` and ``examples`` directories
within each component top-level directory. 

Comprehensive collections of well-designed unit tests and integration 
tests are important tools for developing quality software. Good tests
help to ensure that software operates as expected as it is being developed
and as it is used in different ways. To maintain a high level of usefulness, 
tests must be maintained and developed along with the rest of project
software. Tests should be expanded and modified as software evolves so that 
they can always be used to verify that software functionality is correct.

Unit tests are most effective for designing and modifying
**individual** software units to make sure they continue work properly. 
Integration tests help to make sure that various sets of software units 
work together. Typically, integration testing is more complex and subtle 
than the sum of the independent unit tests of the parts that comprise the
integrated system. Proving that component X and component Y each work 
independently doesn't prove that X and Y are compatible or that they will 
work together. Also, defects in individual components may bear no relationship 
to issues an end user would see.

When developing software, it is important to keep these thoughts in mind and
to use tests effectively to meet your goals. Exposing issues when designing 
and refactoring individual software components may be best accomplished with 
unit tests, often run manually as you are adding or modifying code. Detecting 
broken regressions (e.g., "this used to work, but something changed and now
it doesn't") may be best done by frequently running automated integration tests.

This section describes how to write and manually run tests in Axom.

==================================
A Few Guidelines for Writing Tests
==================================

Before we describe the mechanics of writing tests in the Axom framework,
we describe some test guidelines. The aim of these guidelines is to help 
ensure that tests are complete, correct, easy to run, and easy to modify 
as needed when software changes.

  * Decompose tests so that each test is independent of the others. For 
    example, a unit test file should test features of a single class 
    and not contain tests of other classes.

  * Each specific behavior should be specified in one and only one test.
    For example, all unit tests for a container class may live in a single
    test file, but you should verify each container operation (e.g., 
    container creation/destruction, item insertion, item removal, container 
    traversal, etc.) in exactly one test. In particular, if a test 
    covers some behavior, checking that same behavior in another test
    is unnecessary.

  * Limit each test to as few logical assertions as possible. Ideally, each
    behavior should require one logical assertion. However, sometimes it 
    makes sense to have more than one check. For example, a test for an
    empty container may assert that its 'empty' method returns *true* and 
    also assert that its 'size' method returns *zero*.

  * Tests should be independent on the order in which they are run.

  * Tests should be independent of the platform (hardware architecture,
    compiler, etc.) on which they are run.

  * Tests should be named clearly and consistently. See 
    :ref:`namestargets-label` for a description of Axom conventions for 
    test names.

===================
Unit Tests
===================

In Axom, we use the 
`Google Test framework <https://github.com/google/googletest>`_
for C and C++ unit tests, we use the 
`Fortran Unit Test Framework (FRUIT) <https://sourceforge.net/projects/fortranxunit/>`_ for Fortran unit tests, and we use the `pytest framework <https://docs.pytest.org/en/stable/>`_ for Python unit tests.

Organization of tests in either language/framework are similar should 
follow the principles summarized in the guidelines above. Each Google Test or 
FRUIT file is compiled into its own executable that can be run directly or 
as the build target ``test``. Each executable may contain multiple tests. So that 
running individual tests as needed is not overly burdensome, such as unit 
tests for a C++ class, we put all tests for distinct software units in files 
separate from those for other units. Tests within each file should be sized
so that too many different behaviors are not executed or verified in a 
single test.

See :ref:`namestargets-label` for test file naming and make target conventions.

Google Test (C++/C Tests)
--------------------------

The contents of a typical Google Test file look like this::

  #include "gtest/gtest.h"

  #include ...    // include Axom headers needed to compiler tests in file

  // ...

  TEST(<test_case_name>, <test_name_1>) 
  {
     // Test 1 code here...
     // EXPECT_EQ(...);
  }

  TEST(<test_case_name>, <test_name_2>) 
  {
     // Test 2 code here...
     // EXPECT_TRUE(...);
  }

  // Etc.

Each unit test is defined by the Google Test ``TEST()`` macro which accepts a 
*test case name* identifier, such as the name of the C++ class being tested, 
and a *test name*, which indicates the functionality being verified by the 
test. For each test, logical assertions are defined using 
Google Test `assertion macros`. Failure of expected values will cause the test 
to fail, but other tests will continue to run. 

Note that the Google Test framework will generate a ``main()`` routine for 
each test file if it is not explicitly provided. However, sometimes it is 
necessary to provide a ``main()`` routine that contains operation to run 
before or after the unit tests in a file; e.g., initialization code or 
pre-/post-processing operations. A ``main()`` routine provided in a test 
file should be placed at the end of the file in which it resides.

Here is an example ``main()`` from an Axom test that sets up a slic logger
object to be used in tests:: 

  int main(int argc, char * argv[])
  {
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);

    SimpleLogger logger;  // create & initialize test logger,
                            // finalized when exiting main scope

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    result = RUN_ALL_TESTS();

    return result;
  }

Note that Google Test is initialized first, followed by initialization of the
slic SimpleLogger object. The ``RUN_ALL_TESTS()`` Google Test macro will 
run all the tests in the file. 

As another example, consider a set of tests that use MPI.  The ``main()`` 
routine will initialize and finalize MPI before and after tests are run,
respectively::

  int main(int argc, char * argv[])
  {
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);

    SimpleLogger logger;  // create & initialize test logger,
                          // finalized when exiting main scope

    MPI_Init(&argc, &argv);

    result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
  }

Note that Google test is initialized before ``MPI_Init()`` is called. 

Other Google Test features, such as *fixtures*, may be used as well. 

See the `Google Test Primer <https://github.com/google/googletest/blob/master/googletest/docs/Primer.md>`_ 
for discussion of Google Test concepts, how to use them, and a listing of 
available assertion macros, etc.

FRUIT (Fortran Tests)
--------------------------

Fortran unit tests using the FRUIT framework are similar in structure to 
the Google Test tests for C and C++ described above.

The contents of a typical FRUIT test file look like this::

  module <test_case_name>
    use iso_c_binding
    use fruit
    use <axom_module_name>
    implicit none

  contains

  subroutine test_name_1
  !  Test 1 code here...
  !  call assert_equals(...)
  end subroutine test_name_1

  subroutine test_name_2
  !  Test 2 code here...
  !  call assert_true(...)
  end subroutine test_name_2

  ! Etc.

The tests in a FRUIT test file are placed in a Fortran *module* named for
the *test case name*, such as the name of the C++ class whose Fortran interface
is being tested. Each unit test is in its own Fortran subroutine named
for the *test name*, which indicates the functionality being verified by the
unit test. Within each unit test, logical assertions are defined using
FRUIT methods. Failure of expected values will cause the test
to fail, but other tests will continue to run.

Note that each FRUIT test file defines an executable Fortran program. The
program is defined at the end of the test file and is organized as follows::

  program fortran_test
    use fruit
    use <axom_component_unit_name>
    implicit none
    logical ok

    ! initialize fruit
    call init_fruit

    ! run tests
    call test_name_1
    call test_name_2

    ! compile summary and finalize fruit
    call fruit_summary
    call fruit_finalize

    call is_all_successful(ok)
    if (.not. ok) then
      call exit(1)
    endif
  end program fortran_test

Please refer to the `FRUIT documentation <https://sourceforge.net/projects/fortranxunit/>`_ for more information.

pytest (Python Tests)
--------------------------

Python unit tests use the pytest framework.

The contents of a typical pytest test file look like this::

  def test_name_1():
    #  Test 1 code here...
    #  call assert foo1 == bar1

  def test_name_2():
    #  Test 2 code here...
    #  call assert foo2 == bar2

  # etc.

The tests in a pytest test file are functions with the prefix ``test_``.
Within each unit test, logical assertions are defined using
``assert`` methods. Failure of expected values will cause the test to stop
at the point of failure.

Please refer to the `pytest documentation <https://docs.pytest.org/en/stable/>`_ for more information.

=======================================
CMake Files and Variables for Tests
=======================================

The `CMakeLists.txt` file in component test directory defines the
following items:

  #. Variables for test source files as needed. Separate variables should
     be used for Fortran, C++, etc. For example, ``gtest_sidre_tests`` for
     C++ tests, ``gtest_sidre_C_tests`` for C tests, and ``fruit_sidre_tests``
     for Fortran tests. Note that we use the *Google Test* framework for C
     and C++ tests and *Fruit* for Fortran tests.

  #. An executable and test variable for each test executable to be
     generated. These variables use the ``blt_add_executable`` and
     ``axom_add_test`` macros, respectively, as described above.

.. note:: Fortran executables and tests should be guarded to prevent
          generation when Fortran is not enabled.

See :ref:`testing-label` for details about writing tests in Axom.


===========================
Examples 
===========================

Examples for Axom components serve to illustrate more realistic usage of
those components. They can also be run as tests if that's appropriate.

The source code for each component test should be contained in the component 
``examples`` directory if it is contained in one file. If it contains multiple
files, these should be placed in a descriptively-named subdirectory 
of the ``examples`` directory.

In addition, each example should be given its own CMake-generated make target.


=======================================
CMake Files and Variables for Examples
=======================================

The ``CMakeLists.txt`` file in each component's ``examples`` directory defines the
following items:

  #. Variables for example source files and header files as needed
     Separate variables should be used for Fortran, C++, etc. For example,
     ``example_sources`` for C++, ``F_example_sources`` for Fortran.

  #. An executable and test variable for each example executable to be
     generated and each executable to be run as a test. These definitions
     use the ``blt_add_executable`` and ``axom_add_test`` macros, respectively.
     For example::

       blt_add_executable(NAME  <example executable name>
                          SOURCES <example source>
                          OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                          DEPENDS_ON <example dependencies>)

     and::

       axom_add_test(NAME <example executable name>
                     COMMAND <example executable name>)

     Fortran executables and tests should be guarded to prevent generation if
     Fortran is not enabled.


.. _namestargets-label:

=============================================================================
Filename and CMake Target Conventions for Axom Tests and Examples
=============================================================================

The conventions in this section are intended to make it easy to tell what
is in a given component test or example file and to make it easy to run
desired test or example. In Axom, we use ``make`` targets to build and run 
tests and examples. Typing ``make help`` will list all available targets. When 
the following conventions are followed, all test and example targets for a 
component will be grouped together in this listing. Also, it will be clear 
from each target name what the target is for.

Test file names and make targets
---------------------------------

The format of a test file name is::

  <component name>_<test name>_<optional language specifier>

Examples::

  sidre_buffer.cpp     ('Buffer' class C++ unit test)
  sidre_buffer_C.cpp   ('Buffer' class C unit test)
  sidre_buffer_F.f     ('Buffer' class Fortran unit test)

When test files are named like this, it is easy to see what they contain.
Additionally, when added to the appropriate ``CMakeLists.txt`` file
(see ``src/components/sidre/tests/CmakeLists.txt file`` for example), the 
extension '_test' will be added to the make target name so that the 
test will appear as follows in the make target listing when ``make help``
is typed::

  sidre_buffer_test
  sidre_buffer_C_test
  sidre_buffer_F_test

.. note:: We should also add a target for each component to run all its tests;
          e.g., 'make sidre_tests'


Example file names and make targets
------------------------------------

The format of an example file name is::

  <component name>_<example name>_<optional language specifier>_ex

Examples::

  sidre_shocktube_ex.cpp    ('shocktube' C++ example)
  sidre_shocktube_F_ex.f    ('shocktube' Fortran example)

============================
Running Tests and Examples
============================

Axom examples and tests can be run in different ways. Examples and tests
can be run individually by manually running their executables. Built examples
are in the build directory under ``examples`` and tests are under ``tests``.

Tests can also be ran as a suite with the build target ``test`` from the
build directory.  For example::

  $ make test 

This will run all tests and examples and report a summary of passes and 
failures. Detailed output on individual tests is suppressed.

If a test fails, you can invoke its executable directly to see the detailed
output of which checks passed or failed. This is especially useful when 
you are modifying or adding code and need to understand how unit test details
are working, for example.
