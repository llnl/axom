.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

=======================
MIR User Documentation
=======================

Axom's Material Interface Reconstruction (MIR) component provides algorithms for
reconstructing the interfaces between different materials in multimaterial
meshes. The algorithms take Blueprint meshes
containing a coordset, topology, and matset as input and they output a new Blueprint
node with a new coordset, topology, and matset that contains at most 1 material per zone.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/coretop.html>`_


.. toctree::
   :caption: Contents
   :maxdepth: 1

   mir_algorithms
