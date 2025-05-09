.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Dynamic Mode
******************************************************

The distribution of materials in MultiMat is controlled by the
Cell-Material Relation (CMR). For many MultiMat use cases, this is set
once (static mode) and then fields are defined on the MultiMat object. MultiMat
also supports a dynamic mode that permits materials to move around in the mesh.

When creating a MultiMat object using the default constructor, it will default
to static mode using ``CELL_DOM`` data layout with sparse data. The data layout
argument later impacts the ``addEntry()`` and ``removeEntry()`` methods that modify
the MultiMat object's CMR in dynamic mode. For those methods, when MultiMat is
created with a ``CELL_DOM`` data layout, it means that the first argument to
``addEntry()`` will be a cell number and the second will be a material number.

To convert the MultiMat object to dynamic mode, call the ``convertToDynamic()``
method. This method changes some internal representations (including field 
organization) to better support dynamic modifications of the CMR. For example,
when changing to dynamic mode, ``SPARSE`` fields are converted to ``DENSE``
so further changes do not require field data to be reallocated/reorganized
again. The CMR is modified using calls to the ``addEntry()`` and ``removeEntry()``
methods.

.. literalinclude:: ../../examples/basic.cpp
   :start-after: _multimat_dynamic_mode_begin
   :end-before: _multimat_dynamic_mode_end
   :language: C++
