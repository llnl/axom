.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _records-label:

=======
Records
=======

Sina ``Record`` objects are used to represent the data to be stored for your
study. Some examples for the natural scope of ``Record`` objects include things
like:

    - a single run of an application
    - a `Maestro <https://maestrowf.readthedocs.io/en/stable/>`_ step
    - a cluster of runs that has some metadata attached to the cluster (this
      ``Record`` might have a "contains" :doc:`Relationship <relationships>` for all
      the runs within it)

Each ``Record`` must have an ID and a type that you define. Additionally, Each
``Record`` can have a list of ``File`` objects and a list of ``Datum`` objects.

.. contents:: Topics Covered in this Page
   :depth: 2
   :local:

------
Datums
------

Sina ``Datum`` objects help track the value and (optionally) tags and/or units of a
value associated with a ``Record``. In the Sina schema, a ``Datum`` always belongs
to a ``Record`` or one of ``Record``'s inheriting types.

Some examples of potential ``Datum`` values would be:

    - a scalar
    - a piece of metadata
    - an input parameter

The value of a ``Datum`` may be a string, a double, an array of strings, or an array
of doubles.

Below showcases an example of creating an instance of ``Datum`` with an array of
strings and adding it to a ``Record``:

.. literalinclude:: ../../examples/sina_create_datum.cpp
   :language: cpp

Once executed, this code will output:

.. code:: json

    {
        "data": 
        {
            "my_scalar": 
            {
                "value": 
                [
                    "input"
                ]
            }
        },
        "type": "my_type",
        "local_id": "my_record"
    }

.. _datum-type-label:

+++++++++++++++++++
Checking Datum Type
+++++++++++++++++++

It's possible to check the type of a ``Datum`` with the ``getType()`` method. Types
are tracked in an enumeration called ``ValueType``. The enumeration is as follows:

    - 0: string
    - 1: scalar
    - 2: array of strings
    - 3: array of scalars

Below is an example of this in action:

.. literalinclude:: ../../examples/sina_check_datum_type.cpp
   :language: cpp

++++++++++++++++++++++
Setting Units and Tags
++++++++++++++++++++++

For certain ``Datum`` instances it may be helpful to assign them units and/or tags.
This can be accomplished with the ``setUnits()`` and ``setTags()`` methods respectively.

Below is an example of this functionality:

.. literalinclude:: ../../examples/sina_set_datum_units_tags.cpp
   :language: cpp

+++++++++++++++++++++++++++++++++++++
Viewing Datum From an existing Record
+++++++++++++++++++++++++++++++++++++

Sometimes it's necessary to obtain the current ``Datum`` instances from an existing
``Record``. To do this, you can utilize the ``Record`` object's ``getData`` method.
This method will return an unordered map of ``Datum`` instances.

Below is an example of this process:

.. literalinclude:: ../../examples/sina_view_datum_types.cpp
   :language: cpp

Executing this code will print out:

.. code:: bash

    datum1 is type: 1
    datum2 is type: 0
    datum3 is type: 3

Which, we know from `Checking Datum Type <#datum-type-label>`_, signifies that
datum1 is a scalar, datum2 is a string, and datum3 is an array of scalars.

Using this knowledge we can modify our code to show us the current datum values:

.. literalinclude:: ../../examples/sina_view_datum_values.cpp
   :language: cpp

This will provide the following output:

.. code:: bash

    datum1: 12.34
    datum2: foobar
    datum3: 1 2 20 

-----
Files
-----

Sina ``File`` objects help track the location (URI) and mimetype of a file on the
file system, plus any tags. In the Sina schema, a ``File`` always belongs to a ``Record``
or one of ``Record``'s inheriting types.

Every ``File`` must have a URI, while mimetype and tags are optional.

Below is an example showcasing how to create a file and add it to a record:

.. literalinclude:: ../../examples/sina_file_object_creation.cpp
   :language: cpp

This code will produce the following output:

.. code:: json

    {
        "type": "my_type",
        "local_id": "my_record",
        "files": 
        {
            "/path/to/other/file.txt":
            {
                "tags": 
                [
                    "these",
                    "are",
                    "tags"
                ]
            },
            "/path/to/file.png":
            {
                "mimetype": "image/png"
            }
        }
    }

Similarly, files can be removed from a ``Record`` with the ``remove()`` method:

.. literalinclude:: ../../examples/sina_file_object_removal.cpp
   :language: cpp

As we see from the output, the contents of ``myFile`` are no longer in the
``Record`` instance:

.. code:: json

    {
        "type": "my_type",
        "local_id": "my_record",
        "files": 
        {
            "/path/to/other/file.txt": 
            {
                "tags": 
                [
                    "these",
                    "are",
                    "tags"
                ]
            }
        }
    }

+++++++++++++++++++++++++++++++++++++
Viewing Files From an Existing Record
+++++++++++++++++++++++++++++++++++++

Sometimes it's necessary to view the current ``File`` instances that are stored in
a ``Record``. This can be accomplished by using the ``Record`` object's ``getFiles()``
method which returns an unordered map of ``File`` instances.

Below is an expansion of the previous example where we query the ``Record`` instance
for files:

.. literalinclude:: ../../examples/sina_query_record_for_files.cpp
   :language: cpp

The above code will output:

.. code:: bash

    File with URI '/path/to/file.png' has mimetype 'image/png' and tags ''
    File with URI '/path/to/other/file.txt' has mimetype '' and tags 'these are tags '

----
Runs
----

Sina ``Run`` objects are a subtype of Sina ``Record`` objects corresponding to
a single run of an application as specified in the Sina schema. Similar to ``Record``
objects, ``Run`` objects also require an ID; however, they do *not* require a type as
it will automatically be set to "run". In addition to IDs, there are a few other
required fields:

    - application: the application/code used to create the ``Run``
    - version: the version of the application used to create the ``Run``
    - user: the username of the person who ran the application that generated this ``Run``

The below code shows an example of how a ``Run`` can be created:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create the ID first
        axom::sina::ID run1ID{"run1", axom::sina::IDType::Local};

        // Create a run with:
        // ID: run1ID
        // application: "My Sim Code"
        // version: "1.2.3"
        // user: "jdoe"
        std::unique_ptr<sina::Record> run1{new sina::Run{run1ID, "My Sim Code", "1.2.3", "jdoe"}};
    }
