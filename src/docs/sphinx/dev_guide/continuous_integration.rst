.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _continuous_integration-label:

*******************************
Continuous Integration 
*******************************

The Axom project uses two CI tools,
`GitHub Actions <https://github.com/LLNL/axom/actions>`_
and `GitLab CI <https://docs.gitlab.com/ee/ci/>`_ 
on the LLNL LC Collaboration Zone (CZ).

.. _github_actions-label:

===============
GitHub Actions 
===============

Every Pull Request created on GitHub is automatically run through a series of
CI jobs to ensure that the Axom source builds and passes our unit tests.
These configurations mimic the LC compilers and systems as closely as possible
via Docker containers that have our third-party libraries pre-built on them.

Axom's GitHub project is also configured to require pull requests to pass checks 
from our LC GitLab CI (as described below).


Build and Test Axom in a Docker Container
-----------------------------------------

1. Install Docker Desktop

	* You can look at `dockerdocs <https://docs.docker.com/>`_ for tutorials on installing and using Docker.

2. Ensure the Docker Desktop application is running on your computer
3. Open Terminal/Powershell/Command Line prompt of choice
4. Pull the Docker image

	.. code-block:: bash

		docker pull <image name>

	* Checkout the GitHub Actions `configuration file <https://github.com/LLNL/axom/blob/develop/.github/workflows/ci-tests.yml>`_ for latest images (denoted by ``<Compiler>_IMAGENAME``).

	.. code-block:: bash

		docker images

	* Verify you can see your image

5. Create a Docker volume to save your work

	.. code-block:: bash

		docker volume create <volume name>

	* This allows you to save your progress once you exit the container

	.. code-block:: bash

		docker volume ls

	* Verify you can see your volume

6. Run the Docker container with the volume mounted

	.. code-block:: bash

		docker run -it --mount source=<volume name>,target=/home/axom --name <name you want to give container> <image id>

	* You can get the ``<image id>`` through the ``docker images`` command.
	* This will mount/attach your volume to the Docker image. When you exit your running Docker container, you will want to rerun this command to restart the container with your current progress.

	.. note::

		If you exit the container and attempt to restart, sometimes Docker will return a warning stating the initial container has not fully closed. Just run ``docker container prune`` command to cleanup the container and then you can reissue the run command.

7. In the running Docker container, clone the axom github repository

	.. code-block:: bash

		git clone https://github.com/LLNL/axom.git
		cd axom

8. (optional) Checkout your branch

	.. code-block:: bash

		git checkout <name of your branch>

9. Setup submodules

	.. code-block:: bash

		git submodule update --init

10. Configure Axom

	.. code-block:: bash

		python3 ./config-build.py -hc host-configs/docker/<CMake host-config for container> -bt Release -DENABLE_GTEST_DEATH_TESTS=ON -DBUILD_SHARED_LIBS=ON -DAXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS:BOOL=ON -DENABLE_BENCHMARKS:BOOL=ON

	* The host-config for your Docker container can be found in the ``host-configs/docker/`` directory.
	* The CMake flags (marked with ``-D``) are derived from a job in the GitHub Actions `configuration file <https://github.com/LLNL/axom/blob/develop/.github/workflows/ci-tests.yml>`_. The CMake flags are for mimicking a specific job setup, but are not required to configure Axom.

11. Build Axom

	.. code-block:: bash

		cd <build-* directory>
		make -j 8

12. Run unit tests

You have several potential options:

	.. code-block:: bash

		# Run tests with verbose output using make
		make VERBOSE=1 test

		# Run tests with verbose output using ctest
		ctest -VV

		# Run specific tests with verbose output using ctest, filtering for tests that have a specific prefix
		ctest -VV -R  "<some_test_prefix>*"


.. _gitlab-label:

==========
LC GitLab 
==========

We also maintain a mirror of the `Axom project on LLNL's LC GitLab instance <https://lc.llnl.gov/gitlab/axom/axom>`_
primarily for testing Axom pull requests against the various LC System Types and compilers.

There are two types of GitLab plans.
The first is triggered automatically by pull requests on GitHub,
while the second runs nightly and tests
Axom's ``develop`` branch against a new build of our third-party library stack.

Our GitLab CI configuration also allows manual runs. To initiate a new run, 
navigate to the `CI/CD page, <https://lc.llnl.gov/gitlab/axom/axom/-/pipelines>`_
click on the "Run pipeline" button and select the branch to test.

