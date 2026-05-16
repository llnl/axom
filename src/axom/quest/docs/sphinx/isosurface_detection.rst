.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _isosurface-detection:

********************
Isosurface Detection
********************

Quest provides the ``quest::MarchingCubes`` class for generating contours from
node-centered scalar fields on Conduit Blueprint meshes. In 2D, the algorithm
produces line segments; in 3D, it produces triangles. The output contour mesh
also records which input cell and domain generated each contour facet, which is
useful for analysis, debugging, and downstream reconstruction workflows.

``MarchingCubes`` operates on Blueprint meshes in *multi-domain* form. The
parent topology must be structured and the scalar field must be node-centered.
The class supports both 2D and 3D inputs.

.. Note::

   The current implementation is for the original algorithm:

   Lorensen, William E.; Cline, Harvey E. (1 August 1987).
   "Marching cubes: A high resolution 3D surface construction algorithm".
   *ACM SIGGRAPH Computer Graphics*. 21 (**4**): 163-169

   Other similar or improved algorithms could be added in the future.

.. Note::

   If an input mesh cell contains an isosurface saddle point, the
   isocontour topology is ambiguous.  The algorithm will choose
   the topology arbitrarily but consistently.

.. figure:: figs/planar_and_spherical_isosurfaces.png
   :width: 400px

   Planar isocontour generated using the field :math:`f(\mathbf{r}) =
   f_0 + \mathbf{r} \cdot \mathbf{n}` and spherical contour generated using
   the field :math:`g(\mathbf{r}) = |\textbf{r} - \textbf{r}_0|`. Colors
   denote the domain index in the multi-domain mesh.

Basic workflow
--------------

The workflow is:

#. Construct a ``quest::MarchingCubes`` object with a runtime policy,
   allocator, and data-parallel implementation choice.
#. Call ``setMesh()`` with the Blueprint mesh and topology name. A cell mask
   field name may also be supplied.
#. Select the node-centered scalar field with ``setFunctionField()``.
#. Call ``computeIsocontour()`` for each isovalue of interest.
#. Export the contour to either a ``mint::UnstructuredMesh`` or the raw output
   arrays.

The example application in
``<axom>/src/axom/quest/examples/quest_marching_cubes_example.cpp`` uses the
following setup:

.. literalinclude:: ../../examples/quest_marching_cubes_example.cpp
   :start-after: _quest_marching_cubes_init_start
   :end-before: _quest_marching_cubes_init_end
   :language: C++

After the object is configured, the caller selects a scalar field and computes
the contour:

.. literalinclude:: ../../examples/quest_marching_cubes_example.cpp
   :start-after: _quest_marching_cubes_usage_start
   :end-before: _quest_marching_cubes_usage_end
   :language: C++

Output
------

``MarchingCubes`` stores its output internally until it is exported or cleared.
The simplest output path is to populate a ``mint::UnstructuredMesh``:

.. literalinclude:: ../../examples/quest_marching_cubes_example.cpp
   :start-after: _quest_marching_cubes_output_start
   :end-before: _quest_marching_cubes_output_end
   :language: C++

The generated mesh can optionally contain:

* a field with the parent cell ID for each contour facet
* a field with the parent domain ID for each contour facet

If host-side ``mint`` output is not desired, the class also exposes array-based
accessors for connectivity, node coordinates, parent cell IDs, and parent
domain IDs. Those arrays remain in the allocator space associated with the
``MarchingCubes`` object.

Runtime policies and implementation choices
-------------------------------------------

``MarchingCubes`` accepts an Axom runtime policy, so the contour generation can
run on the CPU or on supported GPU backends. The
``MarchingCubesDataParallelism`` enum controls which implementation is used:

* ``byPolicy`` chooses the implementation that best matches the runtime policy.
* ``hybridParallel`` uses a partially parallel implementation that performs
  well on CPUs.
* ``fullParallel`` uses a more fully parallel implementation that is intended
  for highly parallel devices.

Masking and repeated use
------------------------

The optional mask argument to ``setMesh()`` names a cell-centered integer field
used to restrict contour generation. After a mask field is supplied, the caller
can select which mask value to process by calling ``setMaskValue()`` before
``computeIsocontour()``.

The same ``MarchingCubes`` object can be reused for multiple fields, masks, or
isovalues. ``computeIsocontour()`` appends new facets to the existing output,
while ``clearOutput()`` discards the accumulated contour data.

MPI-parallel runs
-----------------

For MPI-parallel runs, the input mesh may contain local and remote domains. The
algorithm itself is local to each rank and does not require communication. The
generated contour mesh uses locally unique node and cell numbering, so callers
that need globally unique IDs must renumber the output.
