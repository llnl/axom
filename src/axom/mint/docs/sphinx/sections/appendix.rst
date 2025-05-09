.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
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

AXOM_LAMBDA Macro
^^^^^^^^^^^^^^^^^

The ``AXOM_LAMBDA`` convenience macro expands to:

 * ``[=]`` capture by value when the `Axom Toolkit`_ is compiled without CUDA.
 * ``[=] __host__ __device__`` when the `Axom Toolkit`_ is compiled with CUDA

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
