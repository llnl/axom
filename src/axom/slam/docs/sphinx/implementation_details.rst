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

Slam defines several aliases for commonly used configurations (see :doc:`aliases`), 
and its types are extensible to custom storage, legacy storage 
and less common combinations that the alias layer does not intentionally name.

Policies include:

* SizePolicy, StridePolicy, OffsetPolicy (compile time vs. runtime)
* IndirectionPolicy, which chooses how set elements, relation indices and map
  values are reached through backing storage:

  * ``NoIndirection`` for implicit ``element == position`` storage
  * ``ArrayIndirection`` (the default), backed by an ``axom::Array``
    whose lifetime the Slam object manages
  * ``ArrayViewIndirection``, backed by an ``axom::ArrayView`` of a buffer managed elsewhere.
    These are small and trivially copyable, so it is the form used for lightweight, device-capturable Slam objects
  * ``CArrayIndirection`` and ``STLVectorIndirection`` for interoperation 
    with raw-pointer and ``std::vector`` storage
  * custom policies, e.g. ``mfem::Array`` adapters
* SubsettingPolicy (none, virtual parent, concrete parent)
* OwnershipPolicy (local, sidre, other repository)

``Map`` and ``BivariateMap`` default to ``ArrayIndirection``.

The following feature diagram of ``slam::OrderedSet`` policies shows how these policies
interact with the subscript operator:

.. figure:: figs/orderedset_feature_diagram.png
   :figwidth: 100%
   :alt: Feature diagram for slam's ordered set


.. _setup-label:

Simplifying mesh setup
======================

* Builder classes
    * Chained initialization using named-parameter idiom
* Generator classes to simplify types
