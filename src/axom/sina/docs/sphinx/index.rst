.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

===================
Sina C++ User Guide
===================

The Sina ([S]imulation [In]sight and [A]nalysis) C++ library can read and write
JSON files in the Sina schema. It can be used by simulation applications to summarize
run data to be ingested into a database using the Sina tool suite.

The top-level object in the Sina schema is the :doc:`Document <documents>`.
It contains lists of :doc:`Record <records>` and :doc:`Relationship <relationships>`
objects. The example below shows the basics. For more details, see the
:doc:`Tutorial <tutorial>`.

.. literalinclude:: ../../examples/sina_basic.cpp
   :language: cpp

After running the above, the file "MySinaData.json" will contain the
following:

.. code:: json

    {
        "records": [
            {
                "application": "My Sim Code",
                "local_id": "run1",
                "type": "run",
                "user": "jdoe",
                "version": "1.2.3"
            }
        ],
        "relationships": []
    }

.. toctree::
   :caption: Contents
   :maxdepth: 2

   tutorial
   core_concepts
