.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _componentorg-label:

******************************************************
Component Structure
******************************************************

This section describes the structure of directories, files, and their contents
for an Axom component. This section should be used as a guide to identify
tasks to be done when adding a new software component to Axom. These include:

  * Creating the appropriate directory structure
  * Modifying and adding CMake files and variables
  * Generating C and Fortran interfaces
  * Writing documentation
  * Writing tests

.. note:: The discussion here does not contain coding guidelines. Please see
          `Axom Coding Guide <../coding_guide/index.html>`_ 
          for that information.

====================================
Component Directory Structure
====================================

In the ``axom/src/axom`` directory, you will find a subdirectory for
each Axom component. For example::

  $ cd axom/src/axom
  $ ls -1
  CMakeLists.txt
  core
  inlet
  lumberjack
  ...

All files for each component are contained in subdirectories in the
component directory. 

To illustrate, consider the *sidre* component directory::

  $ cd axom/src/axom/sidre
  $ ls -1 -F
  CMakeLists.txt
  core/
  docs/
  examples/
  tests/

Note that, besides directories, the top-level component directory contains
a few files: 

* ``CMakeLists.txt`` contains CMake information for the component in the Axom build system.
    
Components are free to organize their header and source files in whatever 
manner makes sense. For example, in *sidre*, these core header and source files
are in a subdirectory called `core`. As is common practice for C++ libraries, 
associated  header and source files are co-located in the same directories. 

The **docs** directory contains the component documentation. Subdirectories in 
the docs directory are named for each type of documentation. The directories 
`doxygen` and `sphinx` are required. Each Axom component uses *Doxygen* for 
source code documentation and *Sphinx* for user documentation. Other 
documentation directories can be used as needed. For example, *sidre* also 
contains documentation directories `dot` for dot-generated figures, and 
`design` for design documents.

The **interface** directory contains interface files for use by languages 
other than C++. To make it easy for applications written in C and
Fortran, for example, to use Axom directly in their native languages,
Axom components provide APIs in these languages. For information about
how we generate these APIs, see :ref:`shroudfiles-label`.

A **test** directory is required for each component which contains a 
comprehensive set of unit tests. See :ref:`testing-label` for information 
about writing tests and inserting them into our testing framework.

An **examples** directory is optional, but recommended. It contains simple 
code examples illustrating component usage.

.. important:: For consistency, these subdirectory names within the top-level 
               component directory should be the same for each Axom component. 

====================================
CMake Files and Variables
====================================

To properly configure and compile code for a component, and generate 
consistent make targets, existing CMake files and variables need to be
modified in addition to adding CMake files for the new component. In this
section, we describe the sort of changes and additions that are required.
For additional details about our CMake and BLT usage, please look in files
in existing Axom components.

Add CMake macro definitions
------------------------------

The top-level CMake directory ``axom/src/cmake`` contains a file called
`AxomConfig.cmake` that defines macro constants for enabling
Axom components and setting third-party library (TPL) dependencies that 
are used to enforce consistency for conditionally-compiled code. When a new
component or dependency is added, that file must be modified by:

  #. Adding the name of the component to the `COMPS` variable
  #. Adding new TPL dependency to the `TPL_DEPS` variable

The CMake variables are used to generate macro constants in the Axom 
configuration header file. For each new CMake variable added, an associated
``#cmakedefine`` definition must be added to the ``config.hpp.in`` file in the 
``axom/src/include`` directory.

Modify top-level CMakeLists.txt file
----------------------------------------

When adding a new Axom component, the file ``axom/src/components/CMakeLists.txt``
must be modified to hook the component into the CMake build configuration 
system. Specifically:

    #. Add option to enable component. For example,::

         axom_add_component(COMPONENT_NAME sidre DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})

    #. Add component dependency target by adding component name to the ``axom_components`` variable.
    
Add component CMakeLists.txt files
----------------------------------------

There are several ``CMakeLists.txt`` files that must be added in various component
directories. We try to maintain consistent organization and usage across all
Axom components to avoid confusion. To illustrate, we describe the key 
contents of the CMakeLists.txt files in the *sidre* Axom component. See those 
files or those in other components for more details.

Top-level component directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The top-level component directory contains a ``CMakeLists.txt``, e.g., 
``axom/src/components/sidre/CMakeLists.txt``, which contains the following items:

  #. A CMake macro call that checks for any of the required components and/or third-party dependencies
     and errors out with a helpful message, e.g.,::

       axom_component_requires(NAME       Sidre
                               COMPONENTS SLIC
                               TPLS       Conduit )

     .. note:: These dependencies should be limited to the requirements of this singular component.
               Do not list inherited dependencies unless they are used directly in this component.
               Instead ensure that the upstream component has the correct requirements listed.

     .. note:: Optional dependencies should *not* be listed here. Instead toggle their behaviors
               via CMake logic by adding defines, source files, and dependencies.

  #. Subdirectories additions with guards as needed; e.g.,::

       add_subdirectory(src)  

     and::

       if (AXOM_ENABLE_TESTS)
         add_subdirectory(tests)
       endif() 

  #. CMake exports of component targets; e.g.,::

       install(EXPORT <component name>-targets DESTINATION lib/cmake)


Component src directory
^^^^^^^^^^^^^^^^^^^^^^^

The ``CMakeLists.txt`` file in the component ``src`` directory defines:

  #. A variable for component header files named ``<component name>_headers``
  #. A variable for component source files named ``<component name>_sources``
  #. A variable for component dependencies named ``<component name>_depends``

For example, these variables for the *sidre* component are ``sidre_headers``,
``sidre_sources``, and ``sidre_depends``. 

.. note:: It is important to account for all conditional inclusion of items
          in these CMake variable names. For example, a C interface is 
          generated to support a Fortran API, typically. So if Fortran is
          not enabled, it is usually not necessary to include the C header 
          files in `sidre_headers`. Similarly, do not include items in
          the dependency variable if they are not found.

This file also adds source subdirectories as needed (using the CMake 
``add_subdirectory`` command), adds the component as a Axom library, and 
adds target definitions for dependencies. For example, the command to 
add *sidre* as a library is::

  axom_add_library( NAME
                       sidre
                   SOURCES
                       ${sidre_sources}
                       ${sidre_fortran_sources}
                   HEADERS
                       ${sidre_headers}
                   DEPENDS_ON
                       ${sidre_depends}
                   )

All components should follow this format to describe the library information.

Component docs directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A component ``docs`` directory contains a ``sphinx`` that has the
hand-written user documentation that is built and hosted on the Axom's
`ReadTheDocs <https://axom.readthedocs.io/en/develop/index.html>`_
page. These are included by listing them in the Table of Contents 
and in the Documentation section of ``src/index.rst``.

Component tests and examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The content of component ``tests`` and ``examples`` directories, including as
CMake files are discussed in :ref:`testing-label`.

=============================================================================
Filename and CMake Target Conventions for Axom Documentation
=============================================================================

The conventions in this section are intended to make it easy to generate 
a specific piece of documentation for a an Axom component manually. In Axom, 
we use 'make' targets to build documentation. Typing `make help` will list 
all available targets.  When the following conventions are followed, all 
documentation targets for a component will be grouped together in this 
listing. Also, it should be clear from each target name what the target is for.

CMake targets for component user guides and source code docs (i.e., Doxygen) 
are::

  <component name>_user_docs

and ::

  <component name>_doxygen_docs

respectively. For example::

  sidre_user_docs     (sidre component user guide)
  sidre_doxygen_docs  (sidre Doxygen source code docs)


.. _shroudfiles-label:

====================================
C and Fortran Interfaces
====================================

Typically, we use the Shroud tool to generate C and Fortran APIs from our C++ 
interface code. Shroud is a python script that generate code
from a *yaml* file that describes C++ types and their interfaces. It was
developed for the Axom project and has since been generalized and is supported
as a `standalone project <https://github.com/LLNL/shroud>`_.
To illustrate what is needed to generate multi-language API code via a make 
target in the Axom build system, we describe the contents of the *sidre* 
Axom component interface directory ``axom/src/components/sidre/src/interface``
that must be added:

  #. A *yaml* file, named ``sidre_shroud.yaml``, which contains an annotated 
     description of C++ types and their interfaces in *sidre* C++ files.
     This file and its contents are generated manually.

  #. Header files, such as ``sidre.h``, that can be included in C files. Such
     a file includes files containing Shroud-generated 'extern C' prototypes.

  #. Directories to hold the generated files for different languages; e.g.,
     ``c_fortran`` for C and Fortran APIs, ``python`` for python API, etc.

  #. 'Splicer' files containing code snippets that get inserted in the
     generated files.

  #. A ``CMakeLists.txt`` files that contains information for generating CMake
     targets for Shroud to generate the desired interface code. For example::

       add_shroud( YAML_INPUT_FILE sidre_shroud.yaml
            YAML_OUTPUT_DIR yaml
            C_FORTRAN_OUTPUT_DIR c_fortran
            PYTHON_OUTPUT_DIR python
            DEPENDS_SOURCE
                c_fortran/csidresplicer.c c_fortran/fsidresplicer.f
                python/pysidresplicer.c
            DEPENDS_BINARY genfsidresplicer.f
       )

     This tells shroud which *yaml* file to generate code files from, which
     directories to put generated files in, which splicer files to use, etc.

The end result of properly setting up these pieces is a make target called
``generate_sidre_shroud`` that can be invoked to generate *sidre* API code
in other languages Axom supports.

====================================
Python Interfaces
====================================

We use the nanobind library to generate Python APIs from our C++
interface code. Nanobind is a python binding library that generates code
from a *cpp* file that describes C++ functions and their interfaces.

Please refer to the `nanobind documentation <https://nanobind.readthedocs.io/en/latest/>`_ for more information.

The python interpreter can be launched with Axom extension(s) in the PYTHONPATH
by running the convenience script::

  ./bin/run_python_with_axom.sh <optional .py to run>

.. note:: The Python interface requires Axom to be configured with nanobind
          to build and use the interface. This requirement is different from shroud,
          which generates interface files. Once shroud generates the interface
          files, users are not required to configure Axom with shroud to use the
          Fortran interface.

.. warning:: nanobind's numpy interface does not currently support
             `arbitrary Python objects <https://nanobind.readthedocs.io/en/latest/ndarray.html#limitations-related-to-dtypes>`_.

====================================
Documentation
==================================== 

Complete documentation for an Axom component consists of several parts
described in the following sections. All user documentation is accessible 
on `Axom Read The Docs page <https://axom.readthedocs.io>`_.

User Documentation
------------------

Each Axom component uses *Sphinx* for user documentation. This documentation 
is generated by invoking appropriate make targets in our build system.
For example, ``make sidre_docs`` builds *html files* from *Sphinx* user 
documentation for the *sidre* component.

The main goal of good user documentation is to introduce the software to
users so that they can quickly understand what it does and how to use it.
A user guide for an Axom component should enable a new user to get a 
reasonable sense of the capabilities the component provides and what the
API looks like in under 30 minutes. Beyond introductory material, the user
guide should also help users understand all major features and ways the
software may be used. Here is a list of tips to help you write good 
documentation:

  #. Try to limit documentation length and complexity. Using figures,
     diagrams, tables, bulleted lists, etc. can help impart useful 
     information more quickly than text alone.
  #. Use examples. Good examples can help users grasp concepts quickly
     and learn to tackle problems easily.
  #. Place yourself in the shoes of targeted users. Detailed
     instructions may be best for some users, but may be onerous for others
     who can quickly figure things out on their own. Consider providing
     step-by-step instructions for completeness in an appendix, separate
     chapter, via hyperlink, etc. to avoid clutter in sections where you 
     are trying to get the main ideas across.
  #. Try to anticipate user difficulties. When possible, describe workarounds,
     caveats, and places where software is immature to help users set
     expectations and assumptions about the quality and state of your software.
  #. *Test* your documentation. Follow your own instructions completely. 
     If something is unclear or missing, fix your documentation. Working with
     a co-worker who is new to your work, or less informed about it, is
     also a good way to get feedback and improve your documentation.
  #. Make documentation interesting to read. While you are not writing a 
     scintillating novel, you want to engage users with your documentation
     enough so that they don't fall asleep reading it.
  #. Quickly incorporate feedback. When a user provides some useful feedback
     on your documentation, it shows they care enough to help you improve
     it to benefit others. Incorporate their suggestions in a timely fashion
     and ask them if you've addressed their concerns. Hopefully, this will
     encourage them to continue to help.

Speaking of good user documentation, the 
`reStructuredText Primer <http://www.sphinx-doc.org/en/stable/rest.html>`_ 
provides enough information to quickly learn enough to start using the
markdown language for generating sphinx documentation.

Code Documentation
------------------

Each Axom component uses *Doxygen* for code documentation. This documentation 
is generated by invoking appropriate make targets in our build system.
For example, `make sidre_doxygen` builds *html* files from *Doxygen* code 
documentation for the *sidre* component.

The main goal of code documentation is to provide an easily navigable 
reference document of your software interfaces and implementations for
users who need to understand details of your code.

We have a useful discussion of our Doxygen usage conventions in the 
`Documentation Section of the Axom Coding Guide <../coding_guide/sec07_documentation.html>`_. 
The `Doxygen Manual <http://www.doxygen.nl/manual/>`_ contains
a lot more details.

`Axom's code documentation <https://axom.readthedocs.io/en/develop/doxygen/html>`_ 
is published along with our `user documentation. <https://axom.readthedocs.io>`_
