.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _relationships-label:

=============
Relationships
=============

Sina ``Relationship`` objects are used as a way to describe the correlation between
two ``Record`` instances (and/or ``Record`` inheritors, e.g. ``Run``). A ``Relationship``
consists of three parts: a subject, an object, and a predicate. The subject and object
must be IDs referring to valid instances of ``Records``, while the predicate may be
any string.

In describing the connection between objects, a ``Relationship`` is read as "<subject>
<predicate> <object>". For example, in the relationship "Alice knows Bob", "Alice" is
the subject, "knows" is the predicate, and "Bob" is the object. Below are some additional
examples:

    - Task_22 contains Run_1024
    - msub_1_1 describes out_j_1_1
    - Carlos sends an email to Dani
    - local_task_12 runs before local_run_14

A ``Relationship`` should be described in the active voice. **Using active voice in predicates
is recommended** to maintain a clear direction in the relationship. For example, instead of the
passive construction "Dani is emailed by Carlos," use the active form "Carlos emails Dani."

Below is an example showcasing how to construct a ``Relationship`` programmatically. Here,
we assemble a ``Relationship`` showing that "Task_22 contains Run_1024":

.. literalinclude:: ../../examples/sina_relationship_assembly.cpp
   :language: cpp

If executed, the above code will output:

.. code:: json

    {
        "predicate": "contains",
        "subject": "Task_22",
        "object": "Run_1024"
    }

As with any other Sina ID, the subject or object may be either local (uniquely refer to
one object in a Sina file) or global (uniquely refer to one object in a database). Local
IDs are replaced with global ones upon ingestion; all ``Relationship`` instances referring
to that local ID (as well as the ``Record`` possessing that ID) will be updated to use the
same global ID.

Let's add on to our previous example to demonstrate this:

.. literalinclude:: ../../examples/sina_local_id_relationship.cpp
   :language: cpp

In the above code, the "my_local_run" ID would be replaced by a global ID on ingestion.
If this new global ID was, for example, "5Aed-BCds-23G1", then "my_local_run" would
automatically be replaced by "5Aed-BCds-23G1" in both the ``Record`` and ``Relationship``
entries.
