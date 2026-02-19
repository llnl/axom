.. _combiner_class_label:

Combiner Class
==============

The Combiner class is an abstract base class that defines the interface for
all Combiner classes.  Concrete instances need to inherit from this class and
implement these functions to be used when Message classes are combined by the
Lumberjack class.


Functions
---------

========================= ===================
Name                      Description
========================= ===================
id                        Returns the unique differentiating identifier for the class instance.
shouldMessagesBeCombined  Indicates if two messages should be combined.
combine                   Combines the second message into the first.
========================= ===================

Concrete Instances
------------------

.. _texttagcombiner_class_label:

TextTagCombiner
^^^^^^^^^^^^^^^

This Combiner combines the two given Messages if the Message text strings and tag strings are equal.
It does so by adding the second Message's ranks to the first Message (if not past
the ranksLimit) and incrementing the Message's count as well.  This is handled by
Message.addRanks().  The message creation time of the combined message is assigned as the smaller
of the two creation times.

.. note:: This is the only Combiner automatically added to Lumberjack for you.  You can remove it by calling Lumberjack::removeCombiner("TextTagCombiner").

.. _textequalitycombiner_class_label:

TextEqualityCombiner
^^^^^^^^^^^^^^^^^^^^

This Combiner combines the two given Messages if the Message text strings are equal.
It does so by adding the second Message's ranks to the first Message (if not past
the ranksLimit) and incrementing the Message's count as well.  This is handled by
Message.addRanks().  The message creation time of the combined message is assigned as
the smaller of the two creation times.

.. note:: You can add this Combiner by calling Lumberjack::addCombiner(new TextEqualityCombiner).

.. _linefiletagcombiner_class_label:

LineFileTagCombiner
^^^^^^^^^^^^^^^^^^^

This Combiner combines the two given Messages if the Message line number, file name,
and tag are equal.  It does so by adding the second Message's ranks to the first
Message (if not past the ranksLimit) and incrementing the Message's count as well.
This is handled by Message.addRanks().  The message creation time of the combined
message is assigned as the smaller of the two creation times.

.. note:: You can add this Combiner by calling Lumberjack::addCombiner(new LineFileTagCombiner).  The combiner does not compare the text of the two messages passed in.  Therefore, the message texts do not get combined, and the first message text remains unchanged.
