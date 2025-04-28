.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _tutorial-label:

========
Tutorial
========

This short tutorial walks you through the basic usage of the Sina library.
For more in-depth details, see the documentation for the individual classes,
such as Record, Relationship, and Document.

.. contents:: Tutorial Contents
   :depth: 2

------------------------------
Creating Documents and Records
------------------------------

The basic working units in Sina are the Document, Record, and Relationship.
A Document is a collection of Records and Relationships. A Record contains
information about a particular entity, such as the run of an application,
or a description of a uncertainty quantification (UQ) study. A Relationship
describes how two records relate to each user (e.g. UQ studies contain runs).

This first example shows how to create a record:

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin create record]
   :end-before: //! [end create record]

The record has an ID "some_record_id", which is unique to the enclosing
document (it will be replaced by a global ID upon ingestion). The only
required field for records is their type, which is "my_record_type" in this
case. Once created, a record can be added to a Document.

We can create Runs. Runs are special types of records. They have the required
fields of application ("My Sim Code"), version ("1.2.3"), and user ("jdoe").
The type is automatically set to "run".

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin create run]
   :end-before: //! [end create run]

-----------
Adding Data
-----------

Once we have a Record, we can add different types of data to it. Any Datum
object that is added will end up in the "data" section of the record in
the output file.

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin adding data]
   :end-before: //! [end adding data]

-----------------
Adding Curve Sets
-----------------

While you can add data items that are vectors of numbers, sometimes you want
to express relationships between them. For example, you may want to express
the fact that a timeplot captures the fact that there is an independent
variable (e.g. "time"), and possibly multiple dependent variables (e.g.
"temperature" and "energy").

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin curve sets]
   :end-before: //! [end curve sets]

------------
Adding Files
------------

It is also useful to add to a record a set of files that it relates to.
For example your application generated some external data, or you want to
point to a license file.

Conversely, at times it may be necessary to remove a file from the record's file list.
For example if the file was deleted or renamed.

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin file add_and_remove]
   :end-before: //! [end file add_and_remove]

-----------------------------
Relationships Between Records
-----------------------------

Relationships between objects can be captured via the Relationship class.
This relates two records via a user-defined predicate. In the example below,
a new relashionship is created between two records: a UQ study and a run. The
study is said to "contain" the run. As a best practice, predicates should
be active verbs, such as "contains" in "the study contains the run", rather
than "is part of", as in "the run is part of the study".

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin relationships]
   :end-before: //! [end relationships]

---------------------
Library-Specific Data
---------------------

Oftentimes, simulation codes are composed of multiple libraries. If those
offer a capability to collect data in a Sina document, you can leverage that
to expose this additional data in your records.

For example, suppose you are using libraries named ``foo`` and ``bar``.
library ``foo`` defines ``foo_collectData()`` like this:

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin library data foo]
   :end-before: //! [end library data foo]

Library ``bar`` defines ``bar_gatherData()`` like this:

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin library data bar]
   :end-before: //! [end library data bar]

In your host application, you can define sections for ``foo`` and ``bar``
to add their own data.

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin library data host]
   :end-before: //! [end library data host]

In the example above, once the record is ingested into a Sina datastore,
users will be able to search for "temperature" (value = 450),
"foo/temperature" (value = 500), and "bar/temperature" (value = 400).

----------------
Input and Output
----------------

Once you have a document, it is easy to save it to a file. To save to a JSON, we
run the saveDocument() with the optional argument Protocol set to JSON or set as
nothing. Alternatively if you wish to save the document to an HDF5 file: Configure 
axom for HDF5 support then you can set saveDocument()'s optional Protocol parameter 
to HDF5. After executing the below, you will output a file named "my_output.json" 
and a file named "my_output.hdf5", both of which you can ingest into a Sina datastore.

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin io write]
   :end-before: //! [end io write]

If needed, you can also load a document from a file. This can be useful,
for example, if you wrote a document when writing a restart and you want to
continue from where you left off.  To load from a JSON file simply run loadDocument()
with the optional argument Protocol set to JSON or set as nothing. If you've configured 
for, and wish to load from an HDF5 simply set the Protocol to HDF5.

Note that due to HDF5's handling of '/' as indicators for nested structures,
parent nodes will have '/' changed to the ``slashSubstitute`` variable located in 
``axom/sina/core/Document.hpp`` as an HDF5 with saveDocument(). loadDocument() 
will restore them to normal:


.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin io read]
   :end-before: //! [end io read]

---------------------------------
Non-Conforming, User-Defined Data
---------------------------------

While the Sina format is capable of expressing and indexing lots of different
types of data, there may be special cases it doesn't cover well. If you want
to add extra data to a record but don't care if it doesn't get indexed, you
can add it to the "user_defined" section of records (or libraries of
a record). This is a JSON object that will be ignored by Sina for
processing purposes, but will be brought back with your record if you
retrieve it from a database.

Sina uses `Conduit <https://llnl-conduit.readthedocs.io/>`_ to
convert to and from JSON. The user-defined section is exposed as a
`Conduit Node <https://llnl-conduit.readthedocs.io/en/latest/tutorial_cpp_basics.html#node-basics>`_.

.. literalinclude:: ../../examples/sina_tutorial.cpp
   :language: cpp
   :start-after: //! [begin user defined]
   :end-before: //! [end user defined]

