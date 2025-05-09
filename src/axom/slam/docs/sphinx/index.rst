.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Slam User Guide
===============

Axom's Set-theoretic Lightweight API for Meshes (SLAM) component provides high performance
building blocks for distributed-memory mesh data structures in HPC simulation codes.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/slamtop.html>`_


Introduction
------------

Simulation codes have a broad range of requirements for their mesh data structures,
spanning the complexity gamut from structured Cartesian grids to fully unstructured
polyhedral meshes. Codes also need to support features like dynamic topology changes,
adaptive mesh refinement (AMR), submesh elements and ghost/halo layers, in
addition to other custom features.

Slam targets the low level implementation of these distributed mesh data structures and is
aimed at developers who implement mesh data structures within HPC applications.


Set-theoretic abstraction
-------------------------

Slam's design is motivated by the observation that despite vast differences in the high
level features of such mesh data structures, many of the core concepts are shared at a
lower level, where we need to define and process mesh entities and their associated data
and relationships.

Slam provides a simple, intuitive, API centered around a set-theoretic abstraction for
meshes and associated data. Specifically, it models three core set-theoretic concepts:

* **Sets** of entities (e.g. vertices, cells, domains)
* **Relations** among a pair of sets (e.g. incidence, adjacency and containment relations)
* **Maps** defining fields and attributes on the elements of a given set.

The goal is for users to program against Slam's interface without having to be aware of
different design choices, such as the memory layout and underlying data containers. The
exposed API is intended to feel natural to end users (e.g. application developers and
domain scientists) who operate on the meshes that are built up from Slam's abstractions.

See :ref:`srm-label` for more details.


Policy-based design
-------------------

There is considerable variability in how these abstractions can be implemented and user
codes make many different design choices.  For example, we often need different data
structures to support dynamic meshes than we do for static meshes. Similarly, codes
might choose different container types for their arrays (e.g. STL vectors vs. raw C-arrays
vs. custom array types).

Performance considerations can also come in to play. For example, in some cases, a code
has knowledge of some fixed parameters (e.g. the stride between elements in a relation).
Specifying this information at compile-time allows the compiler to better optimize the
generated code than specifying it at runtime.

.. Recognizing that iteration over the mesh entities is often a performance critical
   operation in mesh processing algorithms, Slam attempts to balance the tension between
   generality to allow sharing mesh data and performance.

Slam uses a Policy-based design to orthogonally decompose the feature space without
sacrificing performance. This makes it easier to customize the behavior of Slam's sets,
relations and maps and to extend support for custom features extend the basic interface.

See :ref:`policy-label` for more details.


Current limitations
-------------------

* Slam is under active development with many features planned.
* Support for GPUs in Slam is under development.
* Slam's policy-based design enable highly configurable classes  which are explicitly
  defined via type aliases. We are investigating ways to simplify this set up
  using *Generator* classes where enumerated strings can define related types 
  within a mesh configuration.


.. toctree::
   :caption: Contents
   :maxdepth: 2

   first_example
   core_concepts
   implementation_details
