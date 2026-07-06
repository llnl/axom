.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _srm-label:

=============
Core concepts
=============

This section desribes Slam concepts: what they mean and how they are used.

.. figure:: figs/set_relation_map.png
   :figwidth: 400px
   :alt: Sets, relations and maps in slam
   :align: center

   A **relation** (blue lines) between two **sets** (ovals with red and green dots, as elements)
   and a **map** of scalar values (brown) on the second set.

.. _set-concept-label:

Set
===

Sets model mesh entities such as vertices, zones, materials or refinement levels.
Ordered sets associate each entity with a position so Slam code can iterate
and index efficiently.

Use ``RangeSet`` for contiguous ranges. Use ``ArraySet`` or ``ArrayViewSet``
when set elements are stored in Axom buffers.

.. Future
   Discuss different indexing schemes for ProductSets


.. _relation-concept-label:

Relation
========

A relation connects elements from a `from-set` to those of a `to-set`
and can be used to encode mesh incidence and adjacency relations, such 
as those from cells to vertices, from vertices to cells, from zones to materials
or from elements to neighboring elements.

Slam classifies relations along a few independent axes:

* cardinality: constant per from-set entity or variable per from-set entity
* mutability: static after construction or dynamically editable
* storage: implicit, such as a product set, or explicit, such as index buffers

Use ``ConstantRelation`` / ``ConstantRelationView`` for static fixed-cardinality
relations and ``VariableRelation`` / ``VariableRelationView`` for static
CSR-shaped relations. Use ``DynamicConstantRelation`` or
``DynamicVariableRelation`` when the connectivity needs to be edited.


.. _map-concept-label:

Map
===

Maps attach values to the members of a set. In mesh terms, a map is the Slam
abstraction for scalar, vector or tensor data associated with vertices, cells,
materials or other entity sets.

``ArrayMap`` is the canonical alias for maps that own their data 
and store values in ``axom::Array``. ``ArrayViewMap`` is the canonical alias
for maps that don't own their data.
``ArrayBivariateMap`` and ``ArrayViewBivariateMap`` attach values to a
bivariate set, such as a product set or relation set.

