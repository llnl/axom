.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Slic User Guide
===============

Slic provides a *light-weight*, *modular* and *extensible* logging
infrastructure that simplifies logging application messages.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/slictop.html>`_


Key Features
------------

* Interoperability across the constituent libraries of an
  application. Messages logged by an application and any of its libraries
  using Slic have a unified format and routed to a centralized output
  destination.

* Customizable :ref:`logMessageFormat` to suit application requirements.

* Customizable handling and filtering of log messages by extending the
  :ref:`logStream` base class.

* :ref:`BuiltInLogStreams` to support common logging use cases, e.g., log to
  a file or console.

* Native integration with `Lumberjack <../../../lumberjack/docs/sphinx/index.html>`_ for logging and filtering of messages
  at scale.

* Fortran bindings that provide an idiomatic API for Fortran applications.

Requirements
------------

Slic is designed to be *light-weight* and *self-contained*. The only requirement
for using Slic is a C++11 compliant compiler. However, to use Slic in the
context of a distributed parallel application and in conjunction with
`Lumberjack <../../../lumberjack/docs/sphinx/index.html>`_, support for building with MPI is provided.

For further information on how to build the `Axom Toolkit <../../../../index.html>`_,
consult the `Axom Quick Start Guide <../../../../docs/sphinx/quickstart_guide/index.html>`_.

About this Guide
----------------

This guide presents core concepts, key capabilities, and guiding design
principles of Slic's :ref:`sections/architecture`. To get started with
using Slic quickly within an application, see the
:ref:`sections/getting_started` section. For more detailed information on
the interfaces of the various classes and functions in Slic, developers
are advised to consult the `Slic Doxygen API <../../../../doxygen/html/slictop.html>`_.

Additional questions, feature requests or bug reports on Slic can be submitted
by `creating a new issue on Github <https://github.com/LLNL/axom/issues>`_
or by sending e-mail to the Axom Developers mailing list at axom-dev@llnl.gov.

.. toctree::
   :caption: Contents
   :maxdepth: 3

   sections/getting_started.rst
   sections/architecture.rst
   sections/wrapping_slic_in_macros.rst
   sections/appendix.rst
