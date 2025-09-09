
Lumberjack Classes
==================

Basic
-----

* :ref:`Lumberjack <lumberjack_class_label>` - Performs all high level functionality for the Lumberjack library.
* :ref:`Message <message_class_label>` - Holds all information pertaining to a Message.


Combiners
---------

Handles Message combination and tests whether Message classes should be combined.

* :ref:`Combiner <combiner_class_label>` - Abstract base class that all Combiners must inherit from.
* :ref:`TextTagCombiner <texttagcombiner_class_label>` - Combines Message classes that have equal Text and Tag member variables.
* :ref:`TextEqualityCombiner <textequalitycombiner_class_label>` - Combines Message classes that have equal Text member variables.
* :ref:`LineFileTagCombiner <linefiletagcombiner_class_label>` - Combines Message classes that have equal line number, file name, and tag member variables.


Communicators
-------------

Handles all node-to-node Message passing.

* :ref:`Communicator <communicator_class_label>` - Abstract base class that all Communicators must inherit from.
* :ref:`BinaryTreeCommunicator <binarytreecommunicator_class_label>` - Main Communicator that is implemented with a scalable Binary Tree scheme.
* :ref:`RootCommunicator <rootcommunicator_class_label>` - non-scalable communication scheme that all nodes connect to the root node.  This is given for diagnostic purposes only.
* :ref:`NonCollectiveRootCommunicator <noncollectiverootcommunicator_class_label>` - Communication scheme that all nodes connect to the root node.  Messages are pushed non-collectively.


.. toctree::
   :maxdepth: 1

   lumberjack_class
   combiner_class
   communicator_class
   message_class
