.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _documents-label:

=========
Documents
=========

Sina ``Document`` objects are a way to represent the top-level object of a
JSON file that conforms to the Sina schema. When serialized, these documents
can be ingested into a Sina database and used with the Sina tool.

``Document`` objects follow a very basic JSON layout consisting of two entries:
``records`` and ``relationships``. Each of these entries will store a list of their
respective objects. An example of an empty document is shown below:

.. code:: json

    {
        "records": [],
        "relationships": []
    }

The ``records`` list can contain ``Record`` objects and their inheritying types,
such as ``Run`` objects. The ``relationships`` list can contain ``Relationship``
objects. For more information on these objects, see :doc:`Records <records>`
and :doc:`Relationships <relationships>`.

--------------------
Assembling Documents
--------------------

``Document`` objects can be assembled programatically. To accomplish this:

1. Create a new instance of the ``Document`` class
2. Create a ``Record``
3. Add the instance of the ``Record`` with the ``add`` method

On the `Sina C++ User Guide <./index>`_ page, you can see an example of this
process. Below we will expand on this example to add a ``Relationship``:

.. literalinclude:: ../../examples/sina_document_assembly.cpp
   :language: cpp

After executing the above code, the resulting ``MySinaData.json`` file will
look like so:

.. code:: json

    {
        "records": [
            {
                "type": "run",
                "local_id": "run1",
                "application": "My Sim Code",
                "version": "1.2.3",
                "user": "jdoe"
            },
            {
                "type": "UQ study",
                "local_id": "study1"
            }
        ],
        "relationships": [
            {
                "predicate": "contains",
                "local_subject": "study1",
                "local_object": "run1"
            }
        ]
    }

------------------------------
Generating Documents From JSON
------------------------------

Alternatively to assembling ``Document`` instances programatically, it is
also possible to generate ``Document`` objects from existing JSON files
or JSON strings.

Using our same example from the previous section, if we instead had the
``MySinaData.json`` file prior to executing our code, we could generate
the document using Sina's ``loadDocument()`` function:

.. code:: cpp

    #include "axom/sina.hpp"

    int main (void) {
        axom::sina::Document myDocument = axom::sina::loadDocument("MySinaData.json");
    }

Similarly, if we had JSON in string format we could also load an instance
of the ``Document`` that way:

.. code:: cpp

    #include "axom/sina.hpp"

    int main (void) {
        std::string my_json = "{\"records\":[{\"type\":\"run\",\"id\":\"test\"}],\"relationships\":[]}";
        axom::sina::Document myDocument = axom::sina::Document(my_json, axom::sina::createRecordLoaderWithAllKnownTypes());
        std::cout << myDocument.toJson() << std::endl;
    }


------------------------------
Generating Documents From HDF5
------------------------------

In addition to assembling ``Document`` instances from existing JSON files, it
is possible to generate ``Document`` objects from existing HDF5 files using
Conduit.

When Axom is configured with HDF5 support, Sina's ``saveDocument()`` and ``loadDocument()`` 
functions support HDF5 assembly through the `Protocol::HDF5` argument.  The functions will 
throw a runtime error with the list of available types in Axom configurations if `Protocol::HDF5` 
is attempted without HDF5 support.

.. code:: cpp

    #include "axom/sina.hpp"

    int main (void) {
        axom::sina::Document myDocument = axom::sina::loadDocument("MySinaData.hdf5", axom::sina::Protocol::HDF5);
    }


---------------------------------------------------------
Obtaining Records & Relationships from Existing Documents
---------------------------------------------------------

Sina provides an easy way to query for both ``Record`` and ``Relationship``
objects that are associated with a ``Document``. The ``getRecords()`` and
``getRelationships()`` methods will handle this respectively.

Below is an example showcasing their usage:

.. literalinclude:: ../../examples/sina_query_records_relationships.cpp
   :language: cpp

Running this will show that both records and the one relationship were
properly queried:

.. code:: bash

    Number of Records: 2
    Number of Relationships: 1

------------------------------
Filetype Comparisons
------------------------------

.. toctree::
   :maxdepth: 2

   hdf5_vs_json
