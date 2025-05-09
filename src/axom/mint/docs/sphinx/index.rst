.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Mint User Guide
===============

Mint provides a *comprehensive mesh data model* and a mesh-aware, fine-grain,
parallel execution model that underpins the development of
computational tools and numerical discretization methods. Thereby, enable
implementations that are born *parallel* and *portable* to new and emerging
architectures.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/minttop.html>`_


Key Features
------------

* Support for 1D/2D/3D mesh geometry.

* Efficient data structures to represent :ref:`ParticleMesh`,
  :ref:`StructuredMesh` and :ref:`UnstructuredMesh` types, including
  unstructured meshes with :ref:`MixedCellTopology`.

* Native support for a variety of commonly employed :ref:`CellTypes`.

* A flexible :ref:`MeshStorageManagement` system, which can optionally
  inter-operate with `Sidre <../../../sidre/docs/sphinx/index.html>`_ as the underlying, in-memory, hierarchical
  datastore, facilitating the integration across packages.

* Basic support for :ref:`sections/fem`, consisting of
  commonly employed *shape functions* and *quadratures*.

* A Mesh-Aware :ref:`sections/execution_model`, based on the `RAJA <https://raja.readthedocs.io/en/main/index.html>`_ programming
  model abstraction layer that supports on-node parallelism for mesh-traversals,
  enabling the implementation of computational kernels that are born parallel
  and portable across different processor architectures.


Requirements
------------

Mint is designed to be *light-weight* and *self-contained*.
The only requirement for using Mint is a C++11 compliant compiler.
However, to realize the full spectrum of capabilities, support for
the following third-party libraries is provided:

* `RAJA <https://raja.readthedocs.io/en/main/index.html>`_, used for the parallel execution and portability layer.

* `Conduit <https://llnl-conduit.readthedocs.io/en/latest/index.html>`_ , for using `Sidre <../../../sidre/docs/sphinx/index.html>`_ as the :ref:`MeshStorageManagement` system.

* `Umpire <https://umpire.readthedocs.io/en/develop/>`_, for memory management on next-generation architectures.

For further information on how to build the `Axom Toolkit <../../../../index.html>`_ using these
third-party libraries, consult the `Axom Quick Start Guide <../../../../docs/sphinx/quickstart_guide/index.html>`_.


About this Guide
----------------

This guide discusses the basic concepts and architecture of Mint.

* The :ref:`sections/mint/getting_started` section provides a quick introduction
  to Mint, designed to illustrate high-level concepts and key capabilities, in
  the context of a small working example.

* The :ref:`sections/tutorial` section provides code snippets that
  demonstrate specific topics in a structured and simple format.

* For complete documentation of the interfaces of the various classes and
  functions in Mint consult the `Mint Doxygen API Documentation <../../../../doxygen/html/minttop.html>`_.

* Complete examples and code walk-throughs of mini-apps using Mint are
  provided in the :ref:`sections/examples` section.

Additional questions, feature requests or bug reports on Mint can be submitted
by `creating a new issue on Github <https://github.com/LLNL/axom/issues>`_
or by sending e-mail to the Axom Developers mailing list at axom-dev@llnl.gov.

.. toctree::
   :caption: Contents

   sections/getting_started.rst
   sections/preliminary_concepts.rst
   sections/mesh_representation.rst
   sections/mesh_types.rst
   sections/architecture.rst
   sections/execution_model.rst
   sections/fem.rst
   sections/tutorial.rst
   sections/examples.rst
   sections/faq.rst
   sections/appendix.rst
