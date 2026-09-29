.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _sections/mint/appendix:

Appendix
---------

.. _MintApplicationCodeExample:

Mint Application Code Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is the complete :ref:`MintApplicationCodeExample` presented in
the :ref:`sections/mint/getting_started` section. The code can be found in the Axom
source code under ``src/axom/mint/examples/user_guide/mint_getting_started.cpp``.

.. literalinclude:: ../../../examples/user_guide/mint_getting_started.cpp
   :start-after: sphinx_tutorial_basic_example_start
   :end-before: sphinx_tutorial_basic_example_end
   :language: C++
   :linenos:

.. _axomLambdaMacro:

Host/Device Lambdas
^^^^^^^^^^^^^^^^^^^

Use explicit capture and Axom's host/device annotation for portable kernels:

.. code-block:: C++

   [=] AXOM_HOST_DEVICE(axom::IndexType idx)
   {
     // kernel body
   }

.. _rawSidreData:

Raw Sidre Data
^^^^^^^^^^^^^^

.. literalinclude:: raw_sidre_data.txt
   :language: json
   :linenos:

.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst
