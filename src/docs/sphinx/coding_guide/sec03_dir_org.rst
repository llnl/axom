.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _dirorgsec-label:

=====================================
3 Directory Organization
=====================================

The goal of the guidelines in this section is to make it easy to locate a file
easily and quickly. Make it easy for your fellow developers to find stuff they 
need.

------------------------------------------
Limit scope of directory contents
------------------------------------------

3.1 The contents of each directory and file **must** be well-defined and
limited so that the directory can be named to clearly indicate its contents. 
The goal is to prevent directories and files from becoming bloated with 
too many divergent concepts.


----------------------------------------------------------
Put files where it's easy to find them
----------------------------------------------------------

3.2 Header files and associated implementation files **should** reside in 
the same directory unless there is a good reason to do otherwise. This is 
common practice for C++ libraries.

3.3 Each file **must** reside in the directory that corresponds to (and named
for) the code functionality supported by the contents of the file.

