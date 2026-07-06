.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

============================
Common aliases and migration
============================

Slam's core concepts are **sets**, **relations** and **maps**.
While these are highly configurable, most concrete Slam types are built
from a common set of policy template parameters provided as aliases 
in ``axom/slam/Aliases.hpp``.

Use the aliases when code is meant to communicate Slam intent. 
Use explicit policies when adapting old storage, third-party storage
or a less common policy stack.

Choosing an alias
=================

Sets
----

``RangeSet<P,E>``
  A contiguous set without explicit storage. Use this for dense ranges of mesh
  entities such as cells, nodes, materials or levels.

``ArraySet<P,E>``
  A set whose element ids are stored in an external ``axom::Array``. 
  The array object is owned elsewhere and must outlive the set.

``ArrayViewSet<P,E>``
  A set whose element ids are stored in an ``axom::ArrayView`` held by value.
  Prefer this for non-owning, lightweight views over storage owned elsewhere.

Relations
---------

``ConstantRelation<FromSet, ToSet, N>``
  A static relation where each from-set entity has exactly ``N`` related
  to-set entities. The relation indices are stored in an ``axom::Array``.

``RuntimeConstantRelation<FromSet, ToSet>``
  A static relation with constant cardinality known at runtime.

``VariableRelation<FromSet, ToSet>``
  A static, CSR-shaped relation whose cardinality can vary per from-set entity.
  It stores begin offsets and relation indices in ``axom::Array`` buffers.

The ``*View`` forms use ``axom::ArrayView`` buffers instead of ``axom::Array`` buffers.
They are the preferred form for lightweight, non-owning traversal and kernel-capturable data.

The dynamic relation classes, ``DynamicConstantRelation`` and ``DynamicVariableRelation``, 
keep their existing names because their cardinality or connectivity can be edited 
after construction.

Maps
----

``ArrayMap<S,T,STRIDE>``
  An owning map from set ``S`` to values of type ``T``. Its values are stored in
  an ``axom::Array``. Use this for registry-owned fields, scratch fields and
  other host-constructed data owned by the map.

``ArrayViewMap<S,T,STRIDE>``
  A non-owning map view over an ``axom::ArrayView``. Use this when storage is
  owned elsewhere, or when passing a map-shaped view to generic algorithms or kernels.

``RuntimeArrayMap<S,T>`` and ``RuntimeArrayViewMap<S,T>``
  Use these when the number of components per set element is known only at runtime.

``ArrayBivariateMap<BSet,T,STRIDE>`` and ``ArrayViewBivariateMap<BSet,T,STRIDE>``
  Map values over a bivariate set, such as a product set or relation set. Use these for
  dense/sparse two-axis data, for example material-by-cell or cell-by-material values.

``RuntimeArrayBivariateMap<BSet,T>`` and ``RuntimeArrayViewBivariateMap<BSet,T>``
  Runtime-stride forms of the bivariate map aliases.

Default map storage
===================

``slam::Map`` and ``slam::BivariateMap`` now default to ``policies::ArrayIndirection``. 
That means a map declared without an explicit indirection policy owns an ``axom::Array`` buffer:

.. code-block:: C++

   using Cells = slam::RangeSet<>;
   using Mats  = slam::RangeSet<>;
   using CellMatSet = slam::ProductSet<Cells, Mats>;

   Cells cells(numCells);
   Mats materials(numMaterials);
   CellMatSet cm(&cells, &materials);
   slam::Map<double, Cells> density(&cells);            // axom::Array-backed
   slam::BivariateMap<double, CellMatSet> volfrac(&cm); // axom::Array-backed

The maps in the last two lines can use the following aliases:

.. code-block:: C++

   using Density = slam::ArrayMap<Cells, double>;
   using Volfrac = slam::ArrayBivariateMap<CellMatSet, double>;

Code that requires ``std::vector`` backing can still use those explicitly:

.. code-block:: C++

   using VecIndirection =
     slam::policies::STLVectorIndirection<Cells::PositionType, double>;

   using VectorBackedMap =
     slam::Map<double, Cells, VecIndirection>;


FieldRegistry migration
=======================

``FieldRegistry`` is still a host-side convenience class, but its owning field
and buffer types now follow the canonical map aliases:

.. code-block:: C++

   using CellSet = slam::RangeSet<>;
   using Registry = slam::FieldRegistry<CellSet, double>;

   CellSet cells(numCells);
   Registry fields;
   auto& density = fields.addField("density", &cells); // Registry::MapType
   auto& buffer = fields.addBuffer("tmp", cells.size()); // Registry::BufferType

``Registry::BufferType`` is an ``axom::Array<double>``.

Guidance for migrating from ``std::vector``:

* Use ``auto&`` or ``Registry::BufferType&`` instead of ``std::vector<T>&``.
* Use ``buffer.data()``, ``buffer.size()`` and ``buffer.resize(n)`` for the
  common vector-like operations.
* Use ``buffer.view()`` when an ``axom::ArrayView`` is the right interface.
* Use ``addFieldView()`` or ``addBufferView()`` for externally owned storage.


When explicit policies are still appropriate
============================================

Use the full policy stack instead of these aliases when code needs to:

* Preserve legacy ``std::vector`` storage with ``STLVectorIndirection``.
* Bind raw pointer storage with ``CArrayIndirection``.
* Adapt a third-party container such as an application-specific array type.
* Combine non-default size, offset, stride, subsetting, indirection or interface
  policies that are not named by the alias layer.
* Demonstrate or test the policy contract itself.

