.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Field Concepts
******************************************************

Fields are data that are defined over a mesh, typically with one or more values
per cell. Fields can be scalar, indicating 1-component per cell - or they can
contain multiple components as with vector data (2+ components per cell). This
section discusses important MultiMat field concepts that determine where fields
live on the mesh and how their data are organized in memory.

#######################
Field Mapping
#######################

MultiMat associates materials with cells in a mesh, possibly subdividing cells.
MultiMat includes the concept of field *mapping*, which is where on the mesh the
field data live. Fields can be defined over the cells, which is how most simulations
think about cell-centered fields. With MultiMat, fields can also be defined over
the materials, allowing for compact storage of material-level data. Fields can
also be defined over the cells/material pairs from the Cell-Material Relation (CMR),
allowing fields to have data values for each material in a cell.

.. figure:: figures/mapping.png
   :figwidth: 700px

   Diagram showing field mapping concept.

+--------------------+-----------------------------------------------------------+
| FieldMapping       | Meaning                                                   |
+====================+===========================================================+
| PER_CELL           | The field contains up to ncells * ncomponents values (for |
|                    | dense storage) and there are ncomponents values per cell. |
|                    | For scalars *ncomponents* is 1 so the field length is     |
|                    | ncells.                                                   |
+--------------------+-----------------------------------------------------------+
| PER_MAT            | The field contains nmats * ncomponents values and there   |
|                    | are ncomponents values per material. This mapping allows  |
|                    | fields to be defined over the entire material region and  |
|                    | any cell that uses the material inherits the per-material |
|                    | field value, allowing for very compact storage of         |
|                    | material-level properties.                                |
+--------------------+-----------------------------------------------------------+
| PER_CELL_MAT       | The field contains up to ncells * nmats * ncomponents (for|
|                    | dense storage). This mapping allows materials within a    |
|                    | cell to have their own unique values, which makes them    |
|                    | useful for tracking data at a sub-cell level.             |
+--------------------+-----------------------------------------------------------+

#######################
Data Layout
#######################

Simulation codes contain a variety of algorithms that may have a preference for how
data are arranged in memory to ensure good performance. MultiMat supports
fields with a ``PER_CELL_MAT`` mapping and there are two ways to organize such data.
Fields are said to be **Cell-Dominant** (``CELL_DOM``) if they are stored such that
each cell stores all of its material data to memory before proceeding to data for
the next cell. Fields are **Material-Dominant** (``MAT_DOM``) if the data for all
cells that use the material is stored before proceeding to the next material.
The data layout for multi-material data can be thought of as 2 nested for-loops where
the outer loop is the dominant loop. For example, if iterating over materials and
then cells, the data are stored using ``MAT_DOM`` layout.

+--------------------+----------------------------------------------------------+
| DataLayout         | Meaning                                                  |
+====================+==========================================================+
| CELL_DOM           | Data are stored for each cell and then for each material |
|                    | like this *(c=cell, m=material)*:                        |
|                    |                                                          |
|                    | ``{c0m0, c0m1, c0m2, ..., c1m0, c1m1, c1m2, ...}``       |
+--------------------+----------------------------------------------------------+
| MAT_DOM            | Data are stored for each material and then for each cell |
|                    | like this *(m=material, c=cell)*:                        |
|                    |                                                          |
|                    | ``{m0c0, m0c1, m0c2, m0c3, ... , m1c0, m1c1, m1c2, ...}``|
+--------------------+----------------------------------------------------------+

#######################
Sparsity Layout
#######################

Sparsity primarily concerns fields with ``PER_CELL_MAT`` mapping. When initializing
the MultiMat object, the CMR indicates how materials are distributed
over the mesh. It is completely acceptable for materials to skip over certain cells,
which makes sense if we think about materials as a way to divide up the mesh into
various regions or parts. There are ncells * nmats pairs of data that could be entered
for MultiMat fields. For ``DENSE`` fields, the field must contain ncells * nmats values,
with values present for cell/material pairs whether materials are present or not.
This is an easy way to specify the data but it wastes memory by including extra
values whose only purpose is to keep the rectangular shape of the data array.

For large meshes, compressing out unnecessary values can save a lot of memory. MultiMat
supports a ``SPARSE`` layout that does not include any unnecessary values. If we
regard the CMR as a matrix of true/false values, a user must only provide field values
for ``SPARSE`` data where the CMR contains true values.


.. figure:: figures/sparsity.png
   :figwidth: 700px

   Mixed-material volume fraction field with both DENSE and SPARSE representations.



+--------------------+----------------------------------------------------------+
| SparsityLayout     | Meaning                                                  |
+====================+==========================================================+
| DENSE              | Data are provided for all ncells * nmats pairs, even if  |
|                    | there is no cell/material that is valid.                 |
+--------------------+----------------------------------------------------------+
| SPARSE             | Data are provided for only the cell/material pairs that  |
|                    | are valid according to the CMR.                          |
+--------------------+----------------------------------------------------------+
