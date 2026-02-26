
Lumberjack User Guide
=====================

Lumberjack, named because it cuts down logs, is a C++ library that
provides scalable logging while reducing the amount of messages
written out the screen or file system.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/lumberjacktop.html>`_


Introduction
------------

Lumberjack was created to provide scalable logging with a simple programming
model while allowing developers to customize its behavior. It is named Lumberjack
because it cuts down logs. It uses MPI and a scalable binary tree reduction
scheme to combine duplicate messages and limit output to only the root node.


Requirements
------------

* MPI - MPI is fundamental to Lumberjack and without MPI, Lumberjack is not useful.


Code Guarding
-------------

You tell if Axom was built with Lumberjack enabled by using the following
include and compiler define:

.. code-block:: c

    #include "axom/config.hpp"
    #ifdef AXOM_USE_LUMBERJACK
        // Lumberjack code
    #endif


.. toctree::
   :caption: Contents
   :maxdepth: 2

   quick_start
   core_concepts
   lumberjack_classes
   performance_notes
