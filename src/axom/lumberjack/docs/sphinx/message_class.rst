.. _message_class_label:

Message Class
=============

The Message class holds the information about a single message or multiple messages
that were combined via a Combiner instance.

Information
-----------

The Message class contains the following information. All fields have their respective 
getters and setters.

============  ===================
Name          Description
============  ===================
text          Contents of the message
ranks         Truncated list of where the message originated
count         Total count of how many individual messages occurred
fileName      File name that generated the message
lineNumber    Line number that generated the message
level         Message severity (error, warning, debug, etc.)
creationTime  Creation time of this message
tag           Tag for showing what part of the code generated the message
============  ===================

Functions
---------

The Message class also contains the following helper functions to ease use of the class.

============== ===================
Name           Description
============== ===================
stringOfRanks  Returns a string of the ranks
pack           Returns a packed version of all the message's information
unpack         Takes a packed message and overwrites all the message's information
addRank        Add a rank to the message to the given limit
addRanks       Add ranks to the message to the given limit
============== ===================

