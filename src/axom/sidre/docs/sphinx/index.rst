.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Sidre User Guide
================

The Sidre (Simulation data repository) component of Axom provides tools to
centralize data management in HPC applications: data description, allocation, 
access, and so forth. The goal of Sidre is efficient coordination and sharing 
of data across physics packages in integrated applications, and between 
applications and tools that provide capabilities such as file I/O, *in situ* 
visualization, and analysis.

The design of Sidre is based on substantial production application experience
and requirements for codes to run on advanced architectures. Such codes must 
carefully manage data allocation and placement to run efficiently. Sidre is 
intended to be used in a variety of diverse applications.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/sidretop.html>`_


Introduction
-------------

Sidre provides simple application-level semantics to describe, 
allocate, deallocate, and provide access to data. Currently supported
capabilities include:

* **Separate data description and allocation operations.** This allows 
  applications to describe their data and then decide how best to place 
  the data in memory.
* **Multiple different "views" into a chunk of (shared) data.** A Sidre view
  includes description semantics to define data type, number of elements,
  offset, stride, etc. Thus, a chunk of data in memory can be interpreted
  conceptually in different ways using different views into it.
* **Externally-owned "opaque" or described data.** Sidre can accept a pointer 
  to externally-allocated data and provide access to it by name. When
  external data is described to Sidre, it can be processed in the same ways
  as data that Sidre owns. When data is not described (i.e., it is "opaque"), 
  Sidre can provide access to the data via a pointer, but the consumer of the 
  data must know type information to do anything substantial with the data.
* **Support for different memory spaces and allocation schemes.** Default
  memory space and allocation support is basic. More advanced features are
  available via the `Umpire <https://github.com/LLNL/Umpire>`_ library when
  Axom is configured to be built with it.
* **Attributes or metadata associated with a Sidre view.**  This metadata is
  available to user code to facilitate program logic and is also used
  in Axom to enable selective writing of data subsets to files.
* **Tree-structured data hierarchies.** Many mesh-based application codes 
  organize data into hierarchies of contexts (e.g., domains, regions, blocks, 
  mesh centerings, subsets of elements containing different materials, etc.). 
  Sidre supports hierarchical, tree-based organizations in a simple, flexible 
  way that aligns with the data organization in many applications.
* **APIs for C++, C, and Fortran along with mechanisms to ensure inter-language
  data consistency.**

So far, Sidre development has focused on designing and building flexible and
powerful concepts to build on. The Sidre API includes five main concepts:

* **DataStore.** The main access point to data managed by Sidre; it contains 
  a collection of Buffers, a collection of default Attributes, and a tree 
  structure of Groups.
* **Buffer.** Describes and holds a chunk of data in memory owned by Sidre.
* **Group.** Defines a tree structure like a filesystem, where Groups are like
  folders and Views are like files.
* **View.** Describes a conceptual layout of data in memory. Each View
  has a collection of its explicitly-set Attributes.
* **Attribute.** Provides a metadata descriptor of a View.

These concepts will be described in more detail in later sections.

At this point, Sidre supports simple data types such as scalars, strings, and 
(multidimensional) arrays of scalars. Fundamentally, Sidre does not preclude 
the use of more complex data structures, but does not currently support them
directly. Future Sidre development will expand core functionality to support
additional features such as:

* Support for more complex data types.
* Complex queries and actions involving Attributes.

Support for these enhancements and others will be added based on application
needs and use cases.


.. toctree::
   :caption: Contents
   :maxdepth: 2

   first_example
   core_concepts
   file_io
   parallel_io_concepts
   sidre_conduit
   mfem_sidre_datacollection
