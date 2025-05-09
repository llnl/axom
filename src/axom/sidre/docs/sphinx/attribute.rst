.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _attribute-label:

==========
Attribute
==========

Sidre ``Attribute`` objects enable users to attach metadata (strings and
values) to Sidre views to support queries (e.g., search for views with a given
attribute name), outputting data for a subset of views to files, and other
ways an application may need to selectively process views in a Sidre data
hierarchy.

An attribute is created with a string name and a default scalar or string value.
A default value can be changed later as needed.

.. note:: * ``Attribute`` objects can only be created and destroyed using
            ``DataStore`` class methods. The ``Attribute`` class constructor
            and destructor are private (see :ref:`datastore-label`).
          * Each ``Attribute`` object has a unique name and integer identifier.
            Either can be used to retrieve an ``Attribute`` object from a
            ``DataStore`` object..

Each Sidre view inherits all attributes contained in a datastore with their
default strings or values. Then, an application may explicitly set any
attribute on a view. The application may also query the value of a view
attribute, query whether the attribute was explicitly set, or set the
attribute back to its default value. See :ref:`view-label`
for more information about ``View`` objects.

The ``Attribute`` class provides the following operations:

 * Retrieve the name and unique id of an attribute
 * Set the scalar or string value of an attribute
 * Get the type of an attribute's scalar value

One can iterate through the attributes in the datastore using either "range-for" or iterator syntax:

.. code-block:: C++

    // 'range-for' syntax
    for(auto& attr: datastore->attributes()) { /* ... */ }

    // 'iterator' syntax:
    for(auto it = datastore->attributes().begin(),
          itEnd = datastore->attributes().end(); it != itEnd; ++it)
    {
      auto& attr = *it;
      /* ... */
    }
