.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _tpls-label:

****************************
Third-party Libraries (TPLs)
****************************

Axom dependencies are grouped into four categories: Git submodules,
built-in TPLs in the Axom source tree, system-level
TPLs, and other TPL libraries. The following sections describe how to
install and update these dependencies for Axom. Specifically, the sections
should provide Axom developers with answers to questions such as:

  * How does one add a new compiler or platform to the mix?
  * How does one build a new set of TPLs for a single platform or compiler
    for testing?
  * What is the procedure for changing versions of one or more TPLs?
  * How do we keep things straight when using different TPL versions for 
    different branches?
  * How to use the scripts for team TPL support vs. local development 
    and experimentation?

===========
Determinism
===========

We strive for as close to deterministic behavior in our builds as possible.
By this, we mean that repeated builds should be the same with respect to the 
following:

* Set of libraries with version of each and compile-time options
* Compilers, versions, and compiler flags
* Structure and permissions of installation directory structure and files


===========================================
Build Scripts and Configuration Files
===========================================

There are three levels of build scripts or programs that drive Axom TPL builds.
As you move up levels, away from Spack, the scripts require less configuration
and can build multiple sets of TPLs and/or Axom configurations with a single
script invocation.

The following sections provide brief descriptions of what each level does and 
what important configuration and input files it uses. The sections appear from 
lowest level to highest.

After these levels are described, we discuss Axom development processes that use them.

Level 1: Spack
--------------

Spack is a multi-platform package manager that builds and installs multiple versions
and configurations of software packages. It has recipes for building each package
with available variants to customize it to your needs. For example, Axom
has variants for Fortran, MPI, and others. The recipes drive
the individual package build systems and manage packages they depend on.
Spack also handles system level packages, so you can describe where they are on your
system instead of building them from scratch. You will need to describe which compilers
are available on your system as well.

* Platform specific configuration files live under ``axom/scripts/spack/configs/<platform name>``.
  There is one file (``spack.yaml``) per platform that has sections describing:

   * ``compilers``: This section contains the compiler specs that describe the location
     and any other required information for each compiler. For example, a compiler spec
     contains compiler and version, build and linker flags, etc.
   * ``packages``: This section describes system level packages. For example, version and
     location on the filesystem. This file is very important
     due to its ability to drastically reduce the number of packages that Spack builds.

* Axom specific Spack package files live under ``axom/scripts/spack/packages``. These override
  the package files in Spack's repository under ``var/spack/repos/builtin/packages``.
  We try to minimize these, but we often have to alter the existing Spack packages to apply fixes
  before pushing them up to Spack proper or alter recipes in ways that are Axom specific.
  Such overrides do not happen at the Spack level, but at the next level, Uberenv, 
  described below.

* More detailed information about Spack can be found in the 
  `Spack GitHub repo <https://github.com/spack/spack>`_
  or in the `Spack documentation <https://spack.readthedocs.io/en/latest/>`_

.. note::
   Spack does not stop at the first error it encounters. It attempts to build as many packages
   as possible. As a result, finding the root cause of an error can be difficult. However, looking
   through the log file, whose name will appear in the screen output, for a large indented section 
   will help. The error will be in that section and a message with a path to the full log
   fill will be printed by Spack afterwards. Searching for ``-build-out.txt`` in your output will
   help.

Level 1: Vcpkg
--------------

Vcpkg is an open-source C++ library manager for Windows, Linux, and MacOS by Microsoft.
For Axom, we use it only for Windows TPL builds.

* Project specific package files live under ``axom/scripts/vcpkg_ports``.  There are
  two different files for each package:

   * ``portfile.cmake``: This file is the recipe on how to build the package. Vcpkg
     has strict rules about how your project is laid out and you can do the conversion
     in this file as well.
   * ``vcpkg.json``: This is the manifest file that describes information about the
     package.  For example, dependencies, license information, and optional features.

* More detailed information about Vcpkg can be found in the 
  Vcpkg GitHub repo <https://github.com/microsoft/vcpkg>`_
  and in the `Vcpkg documentation <https://github.com/microsoft/vcpkg#table-of-contents>`_

Level 2: Uberenv
----------------

Uberenv simplifies the use of the two level 1 package managers, Spack and Vcpkg.
We rely on Uberenv for two important things: to collapse multiple Spack commands into
one, and to add as much determinism as possible to our use of the level 1 package managers.
The basic workflow in Uberenv is the following:

#. Setup necessary paths and directories such as the base directory where the
   package manager will be installed.
#. Clone the package manager to a specific Git commit.
#. Apply patches to the package manager. For example, disable extra config scopes in Spack.
#. Add Axom's package repository to Spack, so our packages take precedence.
#. Clean temporary information from previous runs that may bleed into a new run.
#. Optionally create a package source mirror.
#. Install packages via the selected package manager.

* The information provided to Uberenv to start this workflow is defined in one file in
  the top-level Axom source directory:

   * ``.uberenv_config.json``: This file describes project specific configurations,
     such as where to download the package manager, what git commit to use, and
     the top level package to install.

* More detailed information about Uberenv can be found in the 
  `Uberenv GitHub repo <https://github.com/LLNL/uberenv>`_
  and in the `Uberenv documentation <https://uberenv.readthedocs.io/en/latest/>`_

.. note::
   Uberenv warnings and errors are easy to find by searching the output for ``[ERROR:``
   or ``[Warning:``.  Unlike Spack, Uberenv will stop at the first error it encounters.

Level 3: Build Scripts
----------------------

The file ``axom/scripts/spack/specs.json`` contains a list of all specs
that we share for Axom development and GitLab CI testing on the LC platforms
we use for development and testing. The specs automatically handle
platform differences and contain the full list of compilers and package specs
required to build.

The directory ``axom/scripts/llnl_scripts`` contains three "build" scripts that
are designed to build suites of TPLs via Uberenv and Spack.

* ``build_tpls.py``: First, this script builds a set of TPLs for each of the specs
  listed in the ``specs.json`` file for the platform on which it is run. For each TPL set,
  it will generate a host-config file and copy it to the top-level directory of the local
  copy of the Axom repository in which the script is run. After building all of TPL sets
  for a platform, it will attempt to build Axom and the ``using-with-cmake`` example against
  each set of TPLs. This script stops at the first failed TPL build but
  attempts to build the Axom source will each host-config. It will output a summary at
  the end indicating which Axom build succeeded or failed.
* ``build_src.py``: This script uses the existing host-configs in your local clone of the
  Axom repo, or a specific one you point at, and builds and tests Axom against them. It also
  tests the Axom installation via the ``using-with-cmake`` example and Axom tutorials.
* ``build_devtools.py``: This script builds and installs the developer tools
  listed in the ``axom/scripts/spack/packages/axomdevtools/package.py`` Spack
  package. It uses the set of Spack configs located in the
  ``axom/scripts/spack/devtools_config`` directory, so that the regular Spack configs
  can reuse previously built developer tools.

.. note::
   Due to the large amount of information printed to the screen during a full build, the build scripts
   redirect most build step output to log files. This output will tell you what command is being run,
   i.e., ``[exe: some/command --with-options]``, and will tell you the log file being written
   to before it redirects the output from a command, i.e., ``[[log file: /path/to/log``.


=============
Updating TPLs
=============

Git submodules
--------------

Currently, Axom uses four external packages that appear in the project repo
as Git submodules. These are:

  * `BLT <https://github.com/LLNL/blt.git>`_, the CMake-based build
    system we use. It is located in ``axom/src/cmake/blt``.
  * `Axom Data <https://github.com/LLNL/axom_data.git>`_, a collection
    of data files used in testing Axom. It is located in ``axom/data``.
  * `Uberenv <https://github.com/LLNL/uberenv.git>`_, which contains Python
    scripts we use to help automate building third-party dependencies for
    development and deployment. It is located in ``axom/scripts/uberenv``.
  * `RADIUSS Spack Configs <https://github.com/LLNL/radiuss-spack-configs.git>`_,
    which contains Spack packages for some of our LLNL-developed TPLs. It is
    located in ``axom/scripts/spack/radiuss-spack-configs``.

There is no software installation process for these dependencies in the 
traditional sense. To update one of these packages in Axom, simply go into
the directory where the submodule lives in Axom and check out a new version.
If a version is intended to be changed in the Axom repo, make the version change
on a branch and submit a GitHub pull request as you would do for other software
changes. More info on :ref:`building-axom-label`.

Built-in TPLs
-------------

Axom uses several lightweight, header-only libraries internally, which are
exposed for downstream customers to use if they wish. These are:

  * `CLI11 <https://github.com/CLIUtils/CLI11>`_, a command line parser
    for C++ and beyond that provides a rich feature set with a simple and
    intuitive interface.
  * `fmt <https://github.com/fmtlib/fmt>`_, an open-source formatting
    library providing a fast and safe alternative to C stdio and C++ iostreams.
  * `sol <https://github.com/ThePhD/sol2>`_, a C++ library binding to Lua.
  * `Sparsehash <https://github.com/sparsehash/sparsehash>`_, which contains
    several hash-map implementations.

.. note:: Axom patches its built-in TPLs so that they reside in the ``axom`` namespace
   which prevents symbol collisions with other projects, either our
   dependencies or downstream customers who wish to use their own versions.  For
   example, ``fmt::format("foo")`` is ``axom::fmt::format("foo")``.

These TPLs are located in the directory: ``axom/src/thirdparty/axom``. The basic 
instructions on how to update a built-in TPL are:

#. Download the new release and override the source that is already there.
   This may involve removing files no longer needed.

#. Review and apply the existing patch files in the ``axom/src/thirdparty/axom``
   directory. More than likely, you will not be able to directly apply the patch
   file because the source of the library is different than the current version.
   However, the patch files give the general idea of what needs to be changed.
   For example, inclusion in the ``axom`` namespace mentioned above.

#. Ensure that the code builds and tests pass. For more information, please see :ref:`testing-label`.

#. Follow the normal pull request work flow. For more information, please see :ref:`pullrequest-label`.


CLI11
+++++

CLI11 is a 3rd party builtin library that Axom uses to handle command
line processing. Axom packages the library in a header-only format. The CLI11.hpp
header file must be downloaded from a released version of CLI11 since the source
code repository provides sources in separate header files. The patch file contains
many small changes that can be summarized as follows:

#. Add "namspace axom {" near the top of the file.
#. Add "} // namspace axom" near the bottom of the file.
#. Move "#pragma once" from below the copyright to the top of the file.
#. Replace "CLI::" with "axom::CLI::".
#. Replace "Success" with "CLI11_Success". This avoids a symbol collision with X11.


fmt
+++

The fmt formatting library is used for printing program output. Axom packages the library
in header-only format. To update fmt, modify the ``src/thirdparty/axom/fmt/update_fmt.bash``
script to change the version number and any patches that need to be applied. The script
will clone the fmt library and attempt to apply existing patches and update them, if possible.
If patches cannot be applied, then you will need to update them. This is best done by
commenting out the "apply_patches" step and steps that follow it. Run the script to clone
the fmt repository. Next, make the same edits found in the various "*.patch" files to the
files in the current directory. After editing a file, make a new patch:
``patch -u fmt/include/fmt/format.h format.h > format.h.patch`` (be sure to edit the filenames).
Repeat this process until all of the relevant patches are updated. Note that newer versions
of fmt sometimes remove the need for the patches in Axom. Commit the patch files and run
the script to see that they are able to patch the sources. After making the necessary fmt
changes, build Axom on multiple platforms.

.. _local-tpls-label:


Local Third-party Library Installation
--------------------------------------

It is often useful to build a new set of TPLs, other than what we have for GitLab CI testing
and regular development. For example, you may want to try out a new library or version of an
existing library.

.. important:: Running Spack and building TPLs typically requires much more storage
               than you have available in your home directory on an LC system. To
               avoid surpassing your disk space quota, you should run TPL builds
               in a filesystem location with sufficient space. For example, 
               ``/usr/workspace/<username>`` is usually appropriate for this.

From the top-level Axom directory in a local clone of the repo, run the following command
to build all TPLs for all existing compiler specs on the platform you are currently on::

$ ./scripts/llnl_scripts/build_tpls.py -d local/install/path

where ``local/install/path`` is a directory location where you want the 
libraries to be installed.

The TPL build script will output whether each TPL install succeeded and, 
subsequently, whether an Axom build against the TPL install succeeded.

.. note:: When Spack runs, you may see what looks like an error related to ``axom@develop``
          being unable to download. This is not an actual error, but a "feature" of how 
          Spack reports what it's doing, and can be ignored.

Running the script produces new host-config files (i.e., CMake cache files) 
that you can use to build and test Axom against the installation for development or 
if issues arise. The generated host-config files will be placed in the top-level Axom
directory of your local clone of the repo. If any changes to Axom code are 
needed to work with the TPL update(s), make the changes there and test them.

.. note:: You can build a subset of TPLs for a platform, by using the ``uberenv.py``
           script in the top-level Axom directory. For example:: 

            python3 ./scripts/uberenv/uberenv.py --prefix /my/tpl/path --spec clang@10.0.0~cpp14+devtools+mfem+c2c

          will build the TPLs for the clang 10.0.0 compiler, install them
          to the ``/my/tpl/path`` directory, and generate a host-config file
          that you can use to build Axom and its tests. Please see the
          ``scripts/spack/specs.json`` file for a current list of TPL specs
          we use for GitLab CI testing.


Shared Third-party Library Installation Steps
---------------------------------------------

The following instructions describe how to install copies of Axom TPL builds
on Livermore Computing (LC) platforms and recreate our Docker containers
with a new set of TPLs. Typically, this process is followed when you want to 
update one or more TPLs. After they are built and
the associated changes are merged into develop, they will be available for
other Axom developers to use during development, in Axom GitLab CI testing, etc.

#. **Working on a local branch.** 
   Make a local clone of the Axom repo and create a branch to work on.

#. **Changing versions of system packages or other TPLs.**
   To change a version of a system package, which applies to an LC platforms 
   or a Docker container image we use for CI testing on GitHub, go into
   the directory ``axom/scripts/spack/configs``. There you will find a 
   sub-directory for each supported LC system type. Each sub-directory
   has a ``spack.yaml`` file which contains an entry for each system level
   package we rely on. Find the entry for the library you wish to update and 
   change the version number. Do this for each system you want to test/change,
   including configurations in the ``docker`` subdirectory.

   .. note:: Inside of the ``spack.yaml`` file for each system package directory,
             there is a ``compilers`` section containing compiler and 
             version information for compilers we use for development and 
             testing. If you wish to test and build with a new compiler or 
             version on a system, modify the appropriate ``spack.yaml`` 
             file.

   To change a version of a non-system TPL, go into the 
   ``axom/scripts/spack/configs`` directory. There you will find a 
   sub-directory for each system we test on which contains a Spack
   package file ``package.py``. TPL versions are pinned in those package files.
   Modify the contents of the Spack package file ``package.py`` in each
   package sub-directory as needed to change TPL version numbers.

   .. note:: Before continuing, you should test that the installation works
             on all LC systems with the steps in :ref:`local-tpls-label`.


#. **Install TPLs on all required LC machines.**

   When you are confident that everything is correct and working, you will need
   to perform this step on each of the machines named in Axom's standard host-configs.

   .. important:: To install TPL builds to be shared by all Axom developers and used 
                  in our GitLab CI, you will need to become the Axom service user ``atk``.
                  There is a clone of the Axom repo in the ``/usr/workspace/atk/axom_repo``
                  directory. After becoming ``atk``, you can go into that directory and
                  switch to the branch you made to test your changes. Before running the
                  TPL builds, make sure the branch is updated, including all submodules.

   Become the service user ``atk`` via the following command::

   $ xsu atk

   .. note:: This command requires special access permissions. If you need them, contact the Axom team.

   Run the corresponding command for the system you are on::

     # blueos_3_ppc64le_ib_p9 (default cmake is 3.14, need >=3.21 for using-with-cmake example)
     $ module load cmake/3.29.2
     $ lalloc 1 -W 240 scripts/llnl_scripts/build_tpls.py
     
     # toss_4_x86_64_ib
     $ srun -N1 -n 36 --interactive -t 180 scripts/llnl_scripts/build_tpls.py

     # toss_4_x86_64_ib_cray
     $ flux run -N 1 -t 240 scripts/llnl_scripts/build_tpls.py

   .. note:: You may have to adjust the allocation times you ask for the script to complete.

   The ``build_tpls.py`` script will build all third-party libraries for all compilers specs
   for the machine you are on. These will be installed in the shared LC directory
   ``/usr/workspace/axom/libs/<SYS_TYPE>/<time date>/<compiler>``
   used by Axom developers. When completed, they will produce new host-config
   files for each configuration. These host-configs will be located in the top-level directory
   of the Axom repo clone where the script is run and named with the following pattern:
   ``<machine name>-<SYS_TYPE>-<compiler spec>.cmake``. Give these files to your regular user
   account and log back in to that account. Copy these new host-config files to the
   ``host-configs`` subdirectory and commit them to your branch. Make sure all
   file changes from all previous steps are also committed and pushed upstream.

   .. note:: If this step fails, delete the time date stamped directory that was created.
             If you forget to do this, it will eventually be deleted when it is past a certain
             age and no longer needed.

#. **Update and test new Windows builds.**
   We use uberenv with Vcpkg to manage dependencies for our Windows TPL builds.
   The third-party package files, ``portfile.cmake`` and ``vcpkg.json``, may need to be updated
   to reflect the new dependencies.
   To test the Windows updates, go to our
   `Manual Windows TPL builds <https://github.com/LLNL/axom/actions/workflows/test_windows_tpls.yml>`_
   GitHub Actions page. Click on "Actions" and then on "Manual test for Axom's TPLs on Windows" in the "Workflows" menu.
   Find the "Run Workflow" drop-down menu, select your branch, and click on the "Run workflow"
   button. This will launch the tests for Windows.

#. **Build new Docker images.**
   We use pre-built Docker images containing TPLs in our GitHub CI checks.
   To build these, go to our
   `Docker TPL build <https://github.com/LLNL/axom/actions/workflows/docker_build_tpls.yml>`_
   GitHub Actions page. Click on "Actions" and then on "Docker TPL build" in the "Workflows" menu.
   Find the "Run Workflow" drop-down menu, select your branch, and click on the "Run workflow"
   button. This will launch the build of the docker images.

   When the docker image build completes, click on your build and find the
   "Artifacts" listed at the bottom of the page. These contain host-configs
   for building Axom on the docker images. Download them and copy them to the
   ``axom/host-configs/docker`` subdirectory. Rename them to match the corresponding
   host-config.

#. **Update GitHub Actions to the new Docker images.**
   To complete the setup of the new docker images, the ``Compiler_ImageName``
   entries in ``.github/workflows/ci-tests.yml`` at the top-level directory must be updated
   with the timestamped names of the new images. The new names can be found in
   the log files from the successful GitHub action. On the left of the page for
   the successful action is a "Jobs" menu. Click on each job and then find
   the "Get dockerhub repo name" section of the log. The second line of the
   section there should be an entry of the form ``axom/tpls:clang-10_12-18-20_00h-10m``.
   Copy the name beginning with ``axom/tpls`` to the appropriate locations
   in the ``.github/workflows/ci-tests.yml`` file. Repeat this with the names from each compiler
   job used in the GitHub action. 
   Axom's docker images are hosted on our `DockerHub <https://hub.docker.com/r/axom/tpls/tags>`_ page.

#. Make sure all changes in your branch are committed and pushed, and create
   a pull request for a merge to develop. If everything went well, all checks
   on your GitHub PR should pass.
 
