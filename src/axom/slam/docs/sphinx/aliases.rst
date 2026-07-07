.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _aliases-label:

===============================
Choosing policies and aliases
===============================

A Slam set, relation or map is assembled from policy template parameters including:
cardinality, stride, indirection, offset, size, subsetting and interface.
Choosing those policies is how you describe the data structure, i.e. how its
connectivity is shaped, where its data lives, and what is fixed at compile time.
This page explains the choices that come up most often, and a small set of
aliases in ``axom/slam/Aliases.hpp`` that name the most common relation configurations.

The aliases are a convenience for these common types and are preferred when applicable.
Use the policies directly when needed, e.g. for compile-time strides, indirection buffers
backed by ``std::vector`` or C-arrays, or specialized cardinalities.

Where data lives: managed buffers and views
===========================================

The most common choice is the indirection policy, which selects the
buffer that a set, relation or map indexes through. Slam's default is
``policies::ArrayIndirection``, backed by an ``axom::Array``:

.. code-block:: C++

   using Cells = slam::RangeSet<>;
   Cells cells(numCells);

   slam::Map<double, Cells> density(&cells);   // axom::Array indirection (the default)

The two most common indirection policies differ in who manages the buffer's lifetime:

``policies::ArrayIndirection`` (``axom::Array``)
  The Slam object allocates the buffer and frees it when the object is destroyed.
  Lifetime is tied to the Slam object.

``policies::ArrayViewIndirection`` (``axom::ArrayView``)
  The Slam object refers to a buffer that something else (e.g an application array,
  a registry, another Slam object) allocated and will free.
  That buffer must outlive the Slam object. An ``ArrayView`` is small and trivially
  copyable, so it can be handed to a generic algorithm and can be captured in a device kernel.

.. note:: "Manages the buffer" vs. "views a buffer" is a statement about
   lifetime, not about which piece of code logically owns the data.
   A map with an ``axom::Array`` indirection frees its buffer as part of its destruction.
   A map with an ``axom::ArrayView`` indirection leaves that to whoever created the buffer.
   When a Slam object merely wraps an ``axom::Array`` that lives elsewhere,
   an ``ArrayView`` indirection is the accurate description 
   even though the underlying ``axom::Array`` owns its memory.

Two other indirections are supplied for interoperability and to test the interface:
``policies::STLVectorIndirection`` indexes an ``std::vector``, 
and ``policies::CArrayIndirection`` indexes a raw pointer.
Use them when adapting existing storage of that kind.


Fixing parameters at compile time
==================================

Where a quantity is known at compile time, encoding it in a policy lets the
compiler specialize the generated code and data structures.

The stride of a relation or map is the common case. A quad mesh's
element-to-vertex relation has exactly four vertices per element,
so its stride can be a compile-time constant:

.. code-block:: C++

   using ElemVertRelation =
     slam::ConstantRelation<ElemSet, VertSet, /* stride */ 4>;   // CompileTimeStride

When the same count is only known at runtime, use the runtime form
(``RuntimeConstantRelation``, or ``policies::RuntimeStride``) directly.
The same distinction applies to a set's size (``policies::CompileTimeSize`` vs. a runtime size) and offset.

The set and relation aliases
============================

A handful of relation configurations recur across mesh code and require spelling out all
six ``StaticRelation`` policy parameters, so a named shorthand can be helpful.
The aliases below cover them, together with the two indirection set types.
Each has a ``*View`` form that uses an ``axom::ArrayView`` (a buffer managed elsewhere) 
in place of the ``axom::Array`` (a buffer the object manages).

Sets:

``RangeSet<P,E>``
  A contiguous range of positions, with no separate storage. Use it for dense
  ranges of mesh entities such as cells, nodes, materials or levels.

``ArraySet<P,E>`` / ``ArrayViewSet<P,E>``
  A set whose element ids are read from an ``axom::Array`` it manages (``ArraySet``) 
  or from an ``axom::ArrayView`` of a buffer managed elsewhere (``ArrayViewSet``).

Relations:

``ConstantRelation<FromSet, ToSet, N>``
  A static relation with exactly ``N`` to-set entities per from-set entity,
  with ``N`` fixed at compile time.

``RuntimeConstantRelation<FromSet, ToSet>``
  As above, but with the (constant) cardinality supplied at runtime.

``VariableRelation<FromSet, ToSet>``
  A static, CSR-shaped relation whose cardinality varies per from-set entity;
  it holds begin offsets and relation indices.

The dynamic relation classes, ``DynamicConstantRelation`` and
``DynamicVariableRelation``, keep their own names because their connectivity can
be edited after construction.

There are intentionally no map aliases. 
A ``Map`` or ``BivariateMap`` already defaults to an ``axom::Array`` indirection with stride one, 
so the common cases read clearly on their own:

.. code-block:: C++

   slam::Map<double, Cells> density(&cells);            // one value per cell
   slam::BivariateMap<double, CellMatSet> volfrac(&cm); // one value per (cell, material)

The cases that vary a map, e.g. a view into a buffer managed elsewhere, or a runtime stride,
tend to be the cases where a fixed alias name would not capture the configuration cleanly:

.. code-block:: C++

   using CellPos = Cells::PositionType;

   // A field that views a buffer managed elsewhere, with a runtime component count:
   using ViewedField =
     slam::Map<double,
               Cells,
               slam::policies::ArrayViewIndirection<CellPos, double>,
               slam::policies::RuntimeStride<CellPos>>;

   // An std::vector-backed field, for interoperation with existing storage:
   using VectorField =
     slam::Map<double,
               Cells,
               slam::policies::STLVectorIndirection<CellPos, double>>;

FieldRegistry
=============

``FieldRegistry`` is a host-side convenience that keeps a set of named fields and buffers.
Its field type (``Registry::MapType``) is a ``slam::Map`` with the default ``axom::Array`` indirection,
and its buffer type (``Registry::BufferType``) is an ``axom::Array``:

.. code-block:: C++

   using CellSet = slam::RangeSet<>;
   using Registry = slam::FieldRegistry<CellSet, double>;

   CellSet cells(numCells);
   Registry fields;
   auto& density = fields.addField("density", &cells);    // Registry::MapType
   auto& buffer  = fields.addBuffer("tmp", cells.size()); // Registry::BufferType (axom::Array<double>)

Because the registry's buffers are now ``axom::Array`` rather than
``std::vector``, code that consumes them should:

* use ``auto&`` or ``Registry::BufferType&`` rather than ``std::vector<T>&``;
* use ``buffer.data()``, ``buffer.size()`` and ``buffer.resize(n)`` for the
  common vector-like operations;
* call ``buffer.view()`` when an ``axom::ArrayView`` is the right interface;
* register externally-managed storage with ``addFieldView()`` or ``addBufferView()``.
