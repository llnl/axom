.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _winding-number:

****************
Winding Numbers
****************

Quest provides generalized winding number workflows for querying curved and
linearized geometry on MFEM sample meshes. These workflows are useful when the
goal is not just a binary in/out test at quadrature points, but a
winding-number field and its derived in/out classification over a
user-defined query mesh.

At a high level, there are two paths:

* A direct path that evaluates winding number on curved geometry, such as NURBS
  curves in 2D or NURBS patches in 3D. This is the most natural choice when
  the curved representation is the source of truth and preserving that geometry
  in the query is more important than maximizing throughput.
* A linearized path that first replaces the curved shape with segments or
  triangles and then evaluates the same winding-number quantity on that
  discretized geometry. This path is useful when the input is already discrete
  or when a linearized representation is acceptable for the query.

Fast Approximate Methods
------------------------

The "fast" winding-number methods build on the linearized path. They keep the
same basic query quantity, but add preprocessing over the segment or triangle
representation so that large batches of queries can be evaluated more
efficiently. In the current implementation, that acceleration uses a hierarchy
over the linearized geometry together with precomputed moment data to
approximate the contribution of well-separated clusters.

In practice:

* use the direct workflow when curved geometry fidelity matters most
* use the linearized workflow when a segment or triangle approximation is
  already available or acceptable
* use the fast linearized workflow when query counts are high enough that extra
  preprocessing is worth the reduction in per-query cost

Query Outputs
-------------

The example programs ``quest_winding_number_2d.cpp`` and
``quest_winding_number_3d.cpp`` show these workflows end to end, including
query-mesh setup, preprocessing, winding-number evaluation, and creation of the
derived ``inout`` field by rounding the winding-number result.
