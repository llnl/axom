.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core concepts
******************************************************

As noted earlier and described in the previous example, Sidre provides five 
main classes: ``DataStore``, ``Buffer``, ``Group``, ``View``, and ``Attribute``.
In combination, these classes implement a data store with a tree structure 
to organize data in a hierarchy. Before delving into details of these 
classes/concepts, we summarize their basic intent:

* **DataStore** is the main interface to access a data hierarchy.
* **Buffer** describes and holds a contiguous chunk of data in memory.
* **Group** defines parent-child relationships in a hierarchical tree data structure and provides access to file I/O operations.
* **View** provides a virtual description of data and access to it.
* **Attribute** allows a program to attach metadata to View objects for processing data selectively. 

The following sections summarize the main interface features and functionality
of these Sidre classes.

.. note:: Interfaces for each of these classes are provided natively in C++, 
          C, and Fortran. 

.. toctree::
   :maxdepth: 2

   datastore
   buffer
   group
   view 
   attribute
   data_vs_metadata
   query_data
