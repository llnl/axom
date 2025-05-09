.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Axom Core User Guide
====================

The Axom Core component provides fundamental data structures and operations used
throughout Axom.  Different compilers and platforms support these essential 
building blocks in various ways; the Core library provides a uniform interface.

Axom Core contains numeric limits and a matrix representation with operators in
the ``axom::numerics`` namespace.  It contains various utilities including a 
timer class, lexical comparator class, ``processAbort()`` function, and numeric,
filesystem, and string manipulation utilities in the ``axom::utilities``
namespace.  Axom Core also contains the ``axom::Array`` and ``axom::StackArray``
container classes.

.. note:: The Core library does not depend on any other Axom components/libraries. In particular, Core does not use macros from the Axom Slic library for output.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/coretop.html>`_


.. toctree::
   :caption: Contents
   :maxdepth: 1

   core_numerics
   core_utilities
   core_containers
   core_acceleration
