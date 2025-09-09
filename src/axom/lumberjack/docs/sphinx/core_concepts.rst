.. _core_concepts_label:

Core Concepts
=============

The following are core concepts required to understand how Lumberjack works.


.. _combine_label:

Combining
---------

Combining Messages is how Lumberjack cuts down on the number of Messages output from
your program.  It does so by giving the currently held Messages at the current node,
two at a time, to the Combiner classes that are currently registered to Lumberjack
when a Push happens.

Lumberjack provides one Combiner as default, the TextTagCombiner. You can write your own
Combiners and register them with Lumberjack.  The idea is that each Combiner would have
its own criteria for whether a Message should be combined and how to combine that specific
Message with another of the same type.

Combiner's have two main functions, shouldMessagesBeCombined and combine.

The function shouldMessagesBeCombined, returns True if the pair of messages satisfy the associated criteria.  For example in the TextTagCombiner,
if the Text strings and tag strings are exactly equal, it signals they should be combined.

The function combine, takes two Messages and combines them in the way that is specific
to that Combiner class.  For example in the TextTagCombiner, the only thing
that happens is the second Message's ranks gets added to the first and the message count
is increased.  This is because the text strings were equal.  This may not be the case
for all Combiners that you write yourself.


.. _communication_label:

Communication
-------------

Communicating Messages between nodes in an intelligent way is how Lumberjack scales
logging Messages.  The :ref:`Communicator <communicator_class_label>`
class instance handles the specifics on how the communication is implemented.  For
example, it handles where a specific node passes its Messages and which nodes
are allowed to output messages.  As of now, there are three implemented Communicators:
:ref:`BinaryTreeCommunicator <binarytreecommunicator_class_label>`,
:ref:`RootCommunicator <rootcommunicator_class_label>`, and
:ref:`NonCollectiveRootCommunicator <noncollectiverootcommunicator_class_label>`.

BinaryTreeCommunicator, as the name implies, utilizes a standard Binary Tree
algorithm to define how the nodes are connected.  Children pass their
Messages to their parent and the root node is the only node allowed to output Messages.

RootCommunicator has a very simple communication scheme that does not scale well
but is useful in some cases for its simplicity.  All nodes connect directly
to the root node which is also the only node allowed to output Messages.

NonCollectiveRootCommunicator is similar to the RootCommunicator, but messages are
pushed by individual ranks rather than collectively by all ranks.

.. _push_label:

Pushing
-------

A push has three steps: combining, sending, and receiving Messages. When you queue
a Message into Lumberjack, it is held at the node that generated the Message until
you indicate the Lumberjack to push, either once or fully.  If you do not push,
then only the Messages generated at the root node will be outputed.

In a single push, nodes send their currently held Messages to the nodes their are connected
to based on the Communicator's communcation scheme.  For example in the BinaryTreeCommunicator,
children nodes send their Messages to their parent. While the root node only receives Messages.
After a single push, it is not guaranteed that all Messages will be ready to be outputted.

A full push is a number of single pushes until all currently held Messages.  The Communicator
tells the Lumberjack class how many single pushes it takes to fully flush the system of
Messages.  For example in the BinaryTreeCommunicator, it is the log of the number of nodes.
