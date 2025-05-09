.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Multimat User Guide
====================

Axom's MultiMat component is a data management library for multimaterial field data
within multiphysics simulation codes. Simulation codes use materials to overlay extra
parts and details onto a mesh without requiring those features to be modeled
conformally. Instead of using cells to model the part geometry, the geometry is
instead represented using materials and volume fractions. The method for adding
such details is often called *"shaping"* and it is covered in Axom's
:doc:`Klee <../../../../axom/klee/docs/sphinx/index>` and :doc:`Quest <../../../../axom/quest/docs/sphinx/index>`
components.

In addition to representing materials on a mesh, MultiMat is used to define
fields on the mesh, and over material subsets. This enables fields to contain multiple
values where needed for mixed-material cells. MultiMat supports flexible data mappings,
layouts, and dense vs sparse field storage, allowing field data to occupy less memory
than would otherwise be necessary.



API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/multimattop.html>`_


.. toctree::
   :caption: Contents
   :maxdepth: 1

   multimat_materials
   multimat_field_concepts
   multimat_using_fields
   multimat_dynamic_mode
