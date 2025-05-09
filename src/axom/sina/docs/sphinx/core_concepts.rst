.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

=============
Core Concepts
=============

Sina provides four main classes:

    - ``Document`` - represents the top-level object of a JSON file conforming to the Sina schema.
    - ``Record`` - represents the data to be stored.
    - ``Relationship`` - represents a way to define the relationship between two ``Record`` objects.
    - ``CurveSet`` - a class to store related independent and dependent ``Curve`` objects. \
      ``Curve`` and ``CurveSet`` objects are just two more types of data that ``Record`` objects can store.

More details on each class can be found in their respective pages below.

.. toctree::
   :maxdepth: 2

   documents
   records
   relationships
   curve_sets