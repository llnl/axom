.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

****
Axom
****

Axom is an open source project that provides robust and flexible software 
components that serve as building blocks for high performance scientific 
computing applications. A key goal of the project is to have different 
application teams co-develop and share general core infrastructure software 
across their projects instead of individually developing and maintaining 
capabilities that are similar in functionality but are not easily shared.

An important objective of Axom is to facilitate integration of novel, 
forward-looking computer science capabilities into simulation codes. 
A pillar of Axom design is to enable and simplify the exchange of 
simulation data between applications and tools. Axom developers 
emphasize the following principles in software design and implementation:

  * Start design and implementation based on concrete application use cases and maintain flexibility to meet the needs of a diverse set of applications
  * Develop high-quality, robust, high performance software that has well-designed APIs, good documentation, and solid testing
  * Apply consistent software engineering practices across all Axom components so developers can easily work on them
  * Ensure that components integrate well together and are easy for applications to adopt

The main drivers of Axom capabilities originate in the needs of multiphysics
applications in the `Advanced Simulation and Computing (ASC) Program <https://asc.llnl.gov>`_ at `Lawrence Livermore National Laboratory (LLNL) <https://www.llnl.gov>`_ . However, Axom can be employed in a wide range of applications 
beyond that scope, including research codes, proxy application, etc. Often,
developing these types of applications using Axom can facilitate technology
transfer from research efforts into production applications.

==============
Axom Software
==============

Axom software components are maintained and developed on the
`Axom GitHub Project <https://github.com/LLNL/axom>`_. 

.. note:: While Axom is developed in C++, its components have native 
          interfaces in C and Fortran for straightforward usage in 
          applications developed in those languages. Python interfaces 
          are in development.

Our current collection of components is listed here. The number of 
components and their capabilities will expand over time as new needs
are identified.

   * Bump: Blueprint Utilities for Mesh Processing
   * Inlet: Input file parsing and information storage/retrieval
   * Klee: Shaping specification and implementation
   * Lumberjack: Scalable parallel message logging and filtering
   * Mint: Mesh data model
   * Mir: (Material interface reconstruction)
   * Multimat: Managing multimaterial field data
   * Primal: Computational geometry primitives
   * Quest: Querying on surface tool
   * Sidre: Simulation data repository
   * Sina: Simulation insight and analysis
   * Slam: Set-theoretic lightweight API for meshes
   * Slic: Simple Logging Interface Code
   * Spin: Spatial index structures for managing and accelerating spatial searches

=============
Documentation
=============

User guides and source code documentation are always linked on this site.

  * :doc:`Quickstart Guide <docs/sphinx/quickstart_guide/index>`
  *  `Source documentation <doxygen/html/index.html>`__

.. list-table::
   :align: center

   * - Core
     - :doc:`User Guide <axom/core/docs/sphinx/index>`
     - `Source documentation <doxygen/html/coretop.html>`__
   * - Bump
     - :doc:`User Guide <axom/bump/docs/sphinx/index>`
     - `Source documentation <doxygen/html/bumptop.html>`__
   * - Inlet
     - :doc:`User Guide <axom/inlet/docs/sphinx/index>`
     - `Source documentation <doxygen/html/inlettop.html>`__
   * - Klee
     - :doc:`User Guide <axom/klee/docs/sphinx/index>`
     - `Source documentation <doxygen/html/kleetop.html>`__
   * - Lumberjack
     - :doc:`User Guide <axom/lumberjack/docs/sphinx/index>`
     - `Source documentation <doxygen/html/lumberjacktop.html>`__
   * - Mint
     - :doc:`User Guide <axom/mint/docs/sphinx/index>`
     - `Source documentation <doxygen/html/minttop.html>`__
   * - Mir
     - :doc:`User Guide <axom/mir/docs/sphinx/index>`
     - `Source documentation <doxygen/html/mirtop.html>`__
   * - Multimat
     - :doc:`User Guide <axom/multimat/docs/sphinx/index>`
     - `Source documentation <doxygen/html/multimattop.html>`__
   * - Primal
     - :doc:`User Guide <axom/primal/docs/sphinx/index>`
     - `Source documentation <doxygen/html/primaltop.html>`__
   * - Quest
     - :doc:`User Guide <axom/quest/docs/sphinx/index>`
     - `Source documentation <doxygen/html/questtop.html>`__
   * - Sidre
     - :doc:`User Guide <axom/sidre/docs/sphinx/index>`
     - `Source documentation <doxygen/html/sidretop.html>`__
   * - Sina
     - :doc:`User Guide <axom/sina/docs/sphinx/index>`
     - `Source documentation <doxygen/html/sinatop.html>`__
   * - Slam
     - :doc:`User Guide <axom/slam/docs/sphinx/index>`
     - `Source documentation <doxygen/html/slamtop.html>`__
   * - Slic
     - :doc:`User Guide <axom/slic/docs/sphinx/index>`
     - `Source documentation <doxygen/html/slictop.html>`__
   * - Spin
     - :doc:`User Guide <axom/spin/docs/sphinx/index>`
     - `Source documentation <doxygen/html/spintop.html>`__


============================
Component Level Dependencies
============================

Axom has the following inter-component dependencies:

- Core has no dependencies and the other components depend on Core
- Bump depends on Slic and Primal.
- Slic optionally depends on Lumberjack
- Slam, Spin, Primal, Mint, Quest, and Sidre depend on Slic
- Mint depends on Slam, and optionally Sidre
- Mir depends on Bump, Slic, Slam, and Primal.
- Multimat depends on Slic, and Slam
- Inlet depends on Sidre, Slic, and Primal
- Klee depends on Sidre, Slic, Inlet and Primal
- Quest depends on Slam, Spin, Primal, Mint, and, optionally, Klee

The figure below summarizes these dependencies. Solid links indicate hard 
dependencies; dashed links indicate optional dependencies.

.. graphviz:: docs/dependencies.dot


======================================================
Other Tools Application Developers May Find Useful
======================================================

The Axom team develops and supports other software tools that are useful
for software projects independent of the Axom. These include:

  *  `BLT <https://github.com/LLNL/blt>`_ CMake-based build system developed by the Axom team to simplify CMake usage and development tool integration
  *  `Shroud <https://github.com/LLNL/shroud>`_ Generator for C, Fortran, and Python interfaces to C++ libraries, and Fortran and Python interfaces to C libraries
  *  `Conduit <https://github.com/LLNL/conduit>`_ Library for describing and managing in-memory simulation data

===================
Developer Resources
===================

Folks interested in contributing to Axom may be interested in our developer
resource guides.

  * :doc:`Developer Guide <docs/sphinx/dev_guide/index>`
  * :doc:`Coding Guide <docs/sphinx/coding_guide/index>`

================================
Communicating with the Axom Team
================================

--------------
Mailing Lists
--------------

The most effective way to communicate with the Axom team is by using one
of our email lists:

  * 'axom-users@llnl.gov' is for Axom users to contact developers to ask questions, report issues, etc. 
  * 'axom-dev@llnl.gov' is mainly for communication among Axom team members 


---------
Chat Room
---------

We also have the 'Axom Users' chat room on the LLNL Microsoft Teams server. This 
is open to anyone with access to the LLNL network. Just log onto Teams and 
join the room.


======================================================
Axom Copyright and License Information
======================================================

Please see the :ref:`axom-license`.

Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

LLNL-CODE-741217


.. toctree::
   :hidden:

   docs/sphinx/quickstart_guide/index

.. toctree::
   :hidden:
   :titlesonly:
   :caption: Component User Guides

   Core (Widely useful utilities) <axom/core/docs/sphinx/index>
   Bump (Blueprint mesh utilities) <axom/bump/docs/sphinx/index>
   Inlet (Input files) <axom/inlet/docs/sphinx/index>
   Klee (Shaping specification and implementation) <axom/klee/docs/sphinx/index>
   Lumberjack (Scalable parallel message logging and filtering) <axom/lumberjack/docs/sphinx/index>
   Mint (Mesh data model) <axom/mint/docs/sphinx/index>
   Mir (Material Interface Reconstruction) <axom/mir/docs/sphinx/index>
   Multimat (Multimaterial fields) <axom/multimat/docs/sphinx/index>
   Primal (Computational geometry primitives) <axom/primal/docs/sphinx/index>
   Quest (Querying on surface tool) <axom/quest/docs/sphinx/index>
   Sidre (Simulation data repository) <axom/sidre/docs/sphinx/index>
   Sina (Write data in a common format) <axom/sina/docs/sphinx/index>
   Slam (Set-theoretic lightweight API for meshes) <axom/slam/docs/sphinx/index>
   Slic (Simple Logging Interface Code) <axom/slic/docs/sphinx/index>
   Spin (Spatial indexes) <axom/spin/docs/sphinx/index>

.. toctree::
   :hidden:
   :caption: Developer Resources

   docs/sphinx/dev_guide/index
   docs/sphinx/coding_guide/index
   docs/licenses
