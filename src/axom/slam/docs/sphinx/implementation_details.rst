.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

======================
Implementation details
======================

.. _policy-label:

Policy-based design
===================

Slam uses policy-based design to combine the storage and indexing behavior a
type needs without paying for features it does not use.

Common policies include:

* SizePolicy, StridePolicy, OffsetPolicy (compile time vs. runtime)
* IndirectionPolicy, which chooses how set elements, relation indices and map
  values are reached through backing storage:

  * ``NoIndirection`` for implicit ``element == position`` storage
  * ``ArrayIndirection`` for the canonical ``axom::Array``-backed policy used
    when Slam owns or manages an Axom buffer
  * ``ArrayViewIndirection`` for the canonical non-owning ``axom::ArrayView``
    policy used to view externally owned storage and to build lightweight,
    device-capturable Slam objects
  * ``CArrayIndirection`` and ``STLVectorIndirection`` as compatibility and
    reference examples
  * custom policies, e.g. ``mfem::Array`` adapters
* SubsettingPolicy (none, virtual parent, concrete parent)
* OwnershipPolicy (local, sidre, other repository)


Feature diagram of OrderedSet policies (subset).

.. figure:: figs/orderedset_feature_diagram.png
   :figwidth: 100%
   :alt: Feature diagram for slam's ordered set

The figure shows how these policies interact with the subscript operator.


.. _setup-label:

Simplifying mesh setup
======================

* Builder classes
    * Chained initialization using named-parameter idiom
* Generator classes to simplify types
