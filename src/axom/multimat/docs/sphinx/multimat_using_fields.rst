.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Using Fields
******************************************************

Fields store useful simulation quantities. MultiMat supports defining fields
on the mesh and material subsets of the mesh. This section describes how to
create fields and access their data.

#######################
Adding a Field
#######################

The ``addField()`` method adds a field to a MultiMat object. The method
accepts arguments that indicate the mapping, layout, and sparsity for the supplied
data, which are given using an ``axom::ArrayView``. The data given in the view are
copied into new memory managed by MultiMat.

The field mapping argument indicates the space where the data live: the mesh cells,
the materials, or the cells/material regions defined over the mesh. The
data layout argument indicates how the data are organized with respect to cells and
materials. For data that have 1 value per cell, pass ``PER_CELL``. For data that have
1 value per material *(ignoring how many cells use the material)*, pass ``PER_MAT``. For
data that have a unique value per material within a cell, pass ``PER_CELL_MAT``.
For ``PER_CELL_MAT`` data, it is important to know the data layout. Pass ``CELL_DOM``
if all of the material values for a cell are sequential neighbors in memory; otherwise
pass ``MAT_DOM``.

Sparsity layout indicates whether the data array contains the maximum number of
values (numMaterials * numCells for ``DENSE`` fields) or whether it instead contains
only the subset of elements where materials are defined. The length of ``SPARSE``
fields is determined by the number of true values in the Cell-Mesh Relation (CMR).

.. literalinclude:: ../../examples/basic.cpp
   :start-after: _multimat_using_fields_addfields_begin
   :end-before: _multimat_using_fields_addfields_end
   :language: C++

^^^^^^^^^^^^^^^^^^^^^^^
Multi-Component Data
^^^^^^^^^^^^^^^^^^^^^^^

MultiMat can store fields with multiple components (vector data) by passing a non-unity
stride in the last argument when adding a field. Multi-component data are arranged
in memory as a contiguous block where the components of the first element (cell or material)
exist sequentially in memory, followed immediately by the components for the next element,
and so on.

.. literalinclude:: ../../examples/basic.cpp
   :start-after: _multimat_using_fields_multicomponent_begin
   :end-before: _multimat_using_fields_multicomponent_end
   :language: C++

^^^^^^^^^^^
Allocators
^^^^^^^^^^^

MultiMat supports allocating data through allocators. There are 2 separate allocators.
The "Slam" allocator allocates data for internal data structures. The field allocator
is used to allocate field bulk data, which is useful to override when writing GPU
algorithms. Both allocators can be set at once using the ``setAllocatorID()`` method.

* setAllocatorID()
* setSlamAllocatorID()
* setFieldAllocatorID()

#######################
External Field Data
#######################

MultiMat normally allocates its own memory for fields, however fields that point to
externally-allocated memory can also be added using the ``addExternalField()`` method.
This method has the same arguments as the ``addField()`` method, except that the
the supplied view is used as the field's actual data instead of being used to initialize
additional memory. The ``addExternalField()`` method is useful when MultiMat is
managing fields allocated and initialized externally, such as through Sidre, Conduit,
MFEM, etc.

#######################
Removing a Field
#######################

Removing a field is done by calling the ``removeField()`` method on the MultiMat object
and passing the name of the field to be removed. MultiMat will remove the field from
its list of fields and deallocate memory, as needed. For external fields, deallocating the
field's bulk data is the responsibility of the caller.

.. code-block:: cpp

     mm.removeField("myField");

###############
Introspection
###############

The MultiMat object provides methods that permit host codes to determine the number
of fields, their names, and their properties. The ``getNumberOfFields()`` method returns
the number of fields. The ``getFieldName()`` method takes a field index and returns the
name of the field. The ``getFieldIdx()`` method returns the field index for a given field
name.

.. literalinclude:: ../../examples/basic.cpp
   :start-after: _multimat_using_fields_introspection_begin
   :end-before: _multimat_using_fields_introspection_end
   :language: C++

#######################
Accessing Field Data
#######################

Accessing fields and their data is best done when the field's properties are known. The
field mapping determines the mesh subset where the field is defined. Fields with either
``PER_CELL`` or ``PER_MAT`` mappings are defined along one dimension of the numMaterials * numCells
grid so they are *"1D"* fields. Fields with ``PER_CELL_MAT`` field mapping are defined using
both axes of the numMaterials * numCells grid so they are *"2D"* fields.

MultiMat provides separate field access functions for 1D/2D fields. In addition,
there are specific 2D methods to access fields according to whether their data
are dense or sparse. Each of these templated methods returns a field object 
specific to the field's stored data and layout. The field object is used to read/write
the field's data.

* get1dField()
* get2dField()
* getDense2dField()
* getSparse2dField()

^^^^^^^^^^^^^^
Indexing Sets
^^^^^^^^^^^^^^

The data layout for a 2D field determines how it should be traversed. The
MultiMat object provides methods that access the Cell-Material Relation (CMR) and return
indexing sets that are useful for specific materials or cells. For example, if data
use a ``CELL_DOM`` layout then all of the material values for a cell are contiguous in
memory, even though a given cell might not use all possible materials. To write loops
over sparse data that focus on only the valid cell-material pairs from the
CMR, the ``getIndexingSetOfCell()`` and ``getIndexingSetOfMat()`` methods can be called.

.. literalinclude:: ../../examples/basic.cpp
   :start-after: _multimat_using_fields_index_sets_begin
   :end-before: _multimat_using_fields_index_sets_end
   :language: C++

^^^^^^^^^^
1D Fields
^^^^^^^^^^

1D Fields are those with a field mapping of ``PER_CELL`` or ``PER_MAT``. 1D fields can be
retrieved from MultiMat using the ``get1dField()`` method, which returns an object
that can access the field data. The ``get1dField()`` method takes a template argument
for the type of data stored in the field so if double-precision data are stored in
MultiMat then ``get1dField<double>()`` should be called to access the field.

.. literalinclude:: ../../examples/basic.cpp
   :start-after: _multimat_using_fields_1d_start
   :end-before: _multimat_using_fields_1d_end
   :language: C++

1D fields can store multi-component values as well, which adds a small amount of
complexity. The field provides a ``numComp()`` method that returns the number of
components. A component for a given cell is retrieved using the 2-argument
call ``operator()`` by passing the cell index and then the desired component index.

.. literalinclude:: ../../examples/basic.cpp
   :start-after: _multimat_using_fields_1dmc_start
   :end-before: _multimat_using_fields_1dmc_end
   :language: C++

^^^^^^^^^^
2D Fields
^^^^^^^^^^

2D fields are those with a ``PER_CELL_MAT`` field mapping. Since the fields can vary over
materials and cells (in either order) and they can be dense or sparse, there are multiple
ways to iterate over the field data.

Fields can be iterated using access patterns suitable for ``DENSE`` sparsity, even
when the data may be ``SPARSE``. The field's ``findValue()`` function is useful
in this case since it allows 2 indices to be passed in addition to a component index.
The first index is the cell number for ``CELL_DOM`` fields, making the second index the
material number. For ``MAT_DOM`` fields, the order is reversed. This approach to locating
the data is general and can be used to traverse the data in the opposite order of the
native data layout, if desired. However, each call to ``findValue()`` includes a
short search and the method can return ``nullptr`` if no valid cell/material pair
is located for the field.

.. literalinclude:: ../../examples/traversal.cpp
   :start-after: _multimat_using_fields_dense_start
   :end-before: _multimat_using_fields_dense_end
   :language: C++

For algorithms where sparse data traversal is desired, the MultiMat indexing sets
can be used directly as an alternative to dense traversal patterns. Cells are iterated
first in this example since the field has a ``CELL_DOM`` data layout. The materials for
the current cell are queried are used to to compute an index into the field data.

.. literalinclude:: ../../examples/traversal.cpp
   :start-after: _multimat_using_fields_indexset_start
   :end-before: _multimat_using_fields_indexset_end
   :language: C++

MultiMat 2D fields provide iterators as a means for writing simpler code. The iterators
can be used to write sparse data traversal algorithms without the added complexity
of using index sets.

.. literalinclude:: ../../examples/traversal.cpp
   :start-after: _multimat_using_fields_iterator_start
   :end-before: _multimat_using_fields_iterator_end
   :language: C++

Flat iterators can be used to traverse all data in a field. This is useful when
computing values over the field and the algorithm does not need to know the cell
or material associated with the data.

.. literalinclude:: ../../examples/traversal.cpp
   :start-after: _multimat_using_fields_flatiter_start
   :end-before: _multimat_using_fields_flatiter_end
   :language: C++

#############
Conversions
#############

MultiMat can store fields with various data mappings, layouts, and sparsity values.
This allows for great flexibility in how fields are represented in memory. MultiMat
provides conversion routines that allow fields to be converted internally between
various representations.

Field conversion can be done for a variety of reasons. Perhaps fields are converted from
sparse to dense to expose them to an external library that needs dense data. Or, perhaps
dense fields are retrieved from I/O routines and then made sparse in MultiMat during
simulation execution. Other times, depending on the needs of an algorithm, it can make
sense to transpose the data, changing ``CELL_DOM`` to ``MAT_DOM`` or vice-versa. The following
methods perform these data conversions and they take a field index as an argument.

Convert field sparsity:

* convertFieldToSparse()
* convertFieldToDense()

Convert field data layout:

* transposeField()
* convertFieldToMatDom()
* convertFieldToCellDom()
