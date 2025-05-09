.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Coding Guidelines
==============================

These guidelines define code style conventions for the Axom project. Most of the
items were taken from the cited references, sometimes with modifications and 
simplifications; see :ref:`codingrefs-label`.

The guidelines emphasize code readability, correctness, portability, and
interoperability. Agreement on coding style and following common idioms
and patterns provides many benefits to a project with multiple developers.
A uniform "look and feel" makes it easier to read and understand source code,
which increases team productivity and reduces confusion and coding errors
when developers work with code they did not write. Also, guidelines
facilitate code reviews by enforcing consistency and focusing developers on
common concerns. Some of these guidelines are arbitrary, but all are based
on practical experience and widely accepted sound practices. For brevity,
most guidelines contain little detailed explanation or justification.

Each guideline is qualified by one of three auxiliary verbs:
"must", "should", or "may" (or "must not", "should not", "may not").

* A "must" item is an absolute requirement.
* A "should" item is a strong recommendation.
* A "may" item is a potentially beneficial stylistic suggestion.

How to apply "should" and "may" items often depends on the particular code 
situation. It is best to use these in ways that enhance code readability 
and help reduce user and developer errors.

.. important:: * Variations in coding style for different Axom components 
                 is permitted. However, coding style *within each* Axom 
                 component **must** be consistent.
               * Deviations from these guidelines **must** be 
                 agreed upon by the Axom team.
               * When the team agrees to change the guidelines, this
                 guide **must** be updated.


**Contents:**

.. toctree::
   :maxdepth: 3

   sec01_changing_code
   sec02_names
   sec03_dir_org
   sec04_header_org
   sec05_source_org
   sec06_scope
   sec07_documentation
   sec08_design_implement
   sec09_format
   sec10_dev_macros
   sec11_portability
   references

