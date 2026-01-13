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

------------------------------------------
Appending Documents to Existing Sina Files
------------------------------------------

It's normal during a simulation to need to dump data multiple times, for example, to collect quantities
at certain milestones, write timeseries once they reach a certain length, or add new sets
of curves. The append methods such as ``appendDocument()`` cover this case. Simply write your
first document to the filesystem, and when you reach a point where you would like to write another,
use an append function matching the filetype (JSON is recommended for small files, namely anything
not involving timeseries, and HDF5 for anything beefier). There's some nuance to how appending works:

- If your new document contains records not present in the file on disk, it will add them.
  This is generally how you'll want to do snapshots; it determines whether records match by checking the ID,
  so as long as each snapshot has a unique ID, they'll accumulate in the target file.

- If your new document contains records that ARE present in the file on disk, it will attempt to merge
  the data they contain. This is most useful when you have a longer-running simulation where you want
  to accumulate timeseries values. Simply ensure that the record containing new values has the same
  ID as the one containing what you have so far.

- For data, files, user defined, and anything else that isn't curve sets or libraries, it will
  go through field-by-field. If there's a field not already present, it will be added.
  If it IS already present, an optional argument allows you to define the behavior. By default, newest wins,
  but you can also have oldest win or refuse the write.

- For curve sets, it will append new values to the existing dependents/independents. There is some
  checking to ensure all timeseries WITHIN A CURVE SET end up the same length, so while it is possible
  to add new dependents (or independents) to a set of curves while a simulation is running, ensure that
  you're backfilling the required number of values (ex: you dump every 10 cycles, you've dumped 20 times, each
  curve has 200 values. If on your next dump you also add brand_new_curve, it MUST have 210 values.) You'll
  almost always want to define curves up front and add in new values.

- For library_data, append will recurse, following the above rules.

In general, appending is very powerful, but a bit complicated; if Sina encounters any issues, it will
write them to the returned conduit node for troubleshooting.


------------------------------
Filetype Comparisons
------------------------------

.. toctree::
   :maxdepth: 2

   hdf5_vs_json
