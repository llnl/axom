.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _docsec-label: 

========================================
7 Code Documentation
========================================

This section contains guidelines for content and formatting of code
documentation mentioned in earlier sections. The aims of these 
guidelines are to:

   * Document files, data types, functions, etc. consistently.
   * Promote good documentation practices so that essential information is 
     presented clearly and lucidly, and which do not over-burden developers.
   * Generate source code documentation using the Doxygen system.


-----------------------------------------
Document only what's needed
-----------------------------------------

7.1 Documentation **should** only include what is essential for users and 
other developers to easily understand code. Comments **should** be limited to 
describing constraints, pre- and post-conditions, and other issues that 
are important, but not obvious. Extraneous comments (e.g., documenting 
"the obvious") **should** be avoided.

      Code that uses clear, descriptive names (functions, variables, etc.) and 
      clear logical structure is preferable to code that relies on a lot of 
      comments for understanding. To be useful, comments must be understood by 
      others and kept current with the actual code. Generally, maintenance 
      and understanding are better served by rewriting tricky, unclear code 
      than by adding comments to it.


-----------------------------------------
Documenting new code vs. existing code
-----------------------------------------

7.2 New source code **must** be documented following the guidelines in this 
section. Documentation of existing code **should** be modified to conform to 
these guidelines when appropriate. 

7.3 Existing code documentation **should** be improved when its inadequate,
incorrect, or unclear.

.. note:: When code is modified, documentation **must** be changed to reflect 
          the changes.


-----------------------------------------
Write clear documentation
-----------------------------------------

7.4 To make comment text clear and reduce confusion, code comments 
**should** be written in grammatically-correct complete sentences or 
easily understood sentence fragments.


-----------------------------------------
Documentation should be easy to spot
-----------------------------------------

7.5 End-of-line comments **should not** be used to document code logic, 
since they tend to be less visible than other comment forms and may be 
difficult to format cleanly. 

      Short end-of-line comments **may** be useful for labeling closing braces 
      associated with nested loops, conditionals, for scope in general, and 
      for documenting local variable declarations.

7.6 All comments, except end-of-line comments, **should** be indented to 
match the indentation of the code they are documenting. Multiple line comment 
blocks **should** be aligned vertically on the left.

7.7 Comments **should** be clearly delimited from executable code with blank 
lines and "blocking characters" (see examples below) to make them stand out 
and, thus, improve the chances they will be read.

7.8 White space, such as blank lines, indentation, and vertical alignment 
**should** be used in comment blocks to enhance readability, emphasize 
important information, etc.


--------------------------------------------------------------------
General Doxygen usage
--------------------------------------------------------------------

The Doxygen code documentation system uses C or C++ style comment sections 
with special markings and Doxygen-specific commands to extract documentation 
from source and header files. Although Doxygen provides many sophisticated 
documentation capabilities and can generate a source code manual in a variety 
of formats such as LaTeX, PDF, and HTML, these guidelines address only a small 
subset of Doxygen syntax. The goal of adhering to a small, simple set of 
documentation commands is that developers will be encouraged to build useful 
documentation when they are writing code.


Brief vs. detailed comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Doxygen system interprets each documentation comment as either "brief" 
or "detailed". 

 - A brief comment is a concise statement of purpose for an item (usually no 
   more than one line) and starts with the Doxygen command "\\brief" 
   (or "@brief"). Brief comments appear in summary sections of the generated 
   documentation. They are typically seen before detailed comments when 
   scanning the documentation; thus good brief comments make it easier to 
   can or navigate a source code manual.

 - A detailed comment is any comment that is not identified as 'brief'.

7.9 A "brief" description **should** be provided in the Doxygen comment 
section for each of the following items: 

      * A type definition (i.e., class, struct, typedef, enum, etc.) 
      * A macro definition
      * A struct field or class data member
      * A class member function declaration (in the header file class 
        definition) 
      * An unbound function signature (in a header file)
      * A function implementation (when there is no description in the 
        associated header file)

7.10 Important information of a more lengthy nature (e.g., usage examples
spanning multiple lines) **should** be provided for files, major data types 
and definitions, functions, etc. when needed. A detailed comment **must** be 
separated from a brief comment in the same comment block with a line containing
no documentation text.


Doxygen comment blocks
^^^^^^^^^^^^^^^^^^^^^^^

7.11 Doxygen comment blocks **may** use either JavaDoc, Qt style, or one 
of the C++ comment forms described below.

      JavaDoc style comments consist of a C-style comment block starting with
      two \*'s, like this::

         /**
          * ...comment text...
          */

      Qt style comments add an exclamation mark (!) after the opening of a
      C-style comment block,like this::

         /*!
          * ...comment text...
          */

      For JavaDoc or Qt style comments, the asterisk characters ("\*") on
      intermediate lines are optional, but encouraged.

      C++ comment block forms start each line with an additional slash::

         ///
         /// ...comment text...
         ///

      or an exclamation mark::

         //!
         //! ...comment text...
         //!

      For these C++ style comment forms, the comment delimiter is required on
      each line.

7.12 A consistent Doxygen comment block style **must** be used within a component.

7.13 Doxygen comment blocks **must** appear immediately before the items they 
describe; i.e., no blank lines between comment and documented item. This
insures that Doxygen will properly associate the comment with the item.


Doxygen inline comments
^^^^^^^^^^^^^^^^^^^^^^^^

7.14 Inline Doxygen comments **may** be used for class/struct data members, 
enum values, function arguments, etc. 

      When inline comments are used, they **must** appear after the item 
      **on the same line** and **must** use the following syntax::

          /*!< ...comment text... */

      Note that the "<" character must appear immediately after the opening of
      the Doxygen comment (with no space before). This tells Doxygen that the
      comment applies to the item immediately preceding the comment. See
      examples in later sections.

7.15 When an item is documented using the inline form, the comment 
**should not** span multiple lines.


--------------------------------------------------------------------
Copyright and release statement
--------------------------------------------------------------------

7.16 Each file **must** contain a comment section that includes the project
software release information (using whichever comment characters are 
appropriate for the language the file is written in). In the interest of 
brevity, the complete release statement is summarized here to show the 
essential information. The full version can be found in any of the project 
files.

.. note:: Change this when we release the code.

.. code-block:: cpp

   /*
    * Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC.
    * Produced at the Lawrence Livermore National Laboratory.
    *
    * All rights reserved.
    *
    * This source code cannot be distributed without permission and
    * further review from Lawrence Livermore National Laboratory.
    */

See :ref:`headerlayout-label` and :ref:`sourcelayout-label` for guidelines
on placement of copyright and release statement in header and source files,
respectively.


--------------------------------------------------------------------
File documentation
--------------------------------------------------------------------

7.17 Each header file that declares a global type, method, etc. **must** 
have a Doxygen file prologue similar to the following:

.. code-block:: cpp

   /*!
    ***************************************************************************
    *
    * \file ...optional name of file...
    *
    * \brief A brief statement describing the file contents/purpose. (optional)
    *
    * Optional detailed explanatory notes about the file.
    *
    ****************************************************************************
    */

      The "\\file" command **must** appear first in the file prologue. It 
      identifies the comment section as documentation for the file.

      The file name **may** include (part of) the path if the file name is not 
      unique. If the file name is omitted on the line after the "\\file" 
      command, then any documentation in the comment block will belong to 
      the file in which it is located instead of the summary documentation 
      in the listing of documented files.

.. note:: Doxygen requires that a file itself must be documented for
          documentation to be generated for any global item (global function,
          typedef, enum, etc.) defined in the file.

See :ref:`headerlayout-label` and :ref:`sourcelayout-label` for guidelines
on placement of file prologue in header and source files, respectively.


Brief and detailed comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.18 A brief statement of purpose for the file **should** appear as the first 
comment after the file command. If included, the brief statement, **must** be 
preceded by the "\\brief" command.

7.19 Any detailed notes about the file **may** be included after the brief 
comment. If this is done, the detailed comments **must** be separated from 
the brief statement by a line containing no documentation text.


--------------------------------------------------------------------
Type documentation
--------------------------------------------------------------------

7.20 Each type and macro definition appearing in a header file **must** have 
a Doxygen type definition comment prologue immediately before it. For example

.. code-block:: cpp

   /*!
    ****************************************************************************
    *
    * \brief A brief statement of purpose of the type or macro.
    *
    * Optional detailed information that is helpful in understanding the
    * purpose, usage, etc. of the type/macro ...
    *
    * \sa optional cross-reference to other types, functions, etc...
    * \sa etc...
    *
    * \warning This class is only partially functional.
    *
    ****************************************************************************
    */

.. note:: Doxygen requires that a compound entity, such as a class, struct, 
          etc. be documented in order to document any of its members.


Brief and detailed comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.21 A brief statement describing the type **must** appear as the first text 
comment using the Doxygen command "\\brief".

7.22 Important details about the item **should** be included after the brief 
comment and, if included, **must** be separated from the brief comment by a 
blank line.


Cross-references and caveats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.23 Cross-references to other items, such as other related types **should** 
be included in the prologue to enhance the navigability of the documentation. 

      The Doxygen command "\\sa" (for "see also") **should** appear before each
      such cross-reference so that links are generated in the documentation.

7.24 Caveats or limitations about the documented type **should** be noted 
using the "\\warning" Doxygen command as shown above.


--------------------------------------------------------------------
Function documentation
--------------------------------------------------------------------

7.25 Each unbound function **should** be be documented with a function 
prologue in the header file where its prototype appears or in a source file 
immediately preceding its implementation.

7.26 Since C++ class member functions define the class interface, they 
**should** be documented with a function prologue immediately preceding 
their declaration in the class definition.


Example function documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following examples show two function prologue variations that may 
be used to document a method in a class definition. The first shows how
to document the function arguments in the function prologue.

.. code-block:: cpp

      /*!
       *************************************************************************
       *
       * \brief Initialize a Foo object with given operation mode.
       *
       * The "read" mode means one thing, while "write" mode means another.
       *
       * \return bool indicating success or failure of initialization.
       *              Success returns true, failure returns false.
       *
       * \param[in] mode OpMode enum value specifying initialization mode.
       *                 ReadMode and WriteMode are valid options.
       *                 Any other value generates a warning message and the
       *                 failure value ("false") is returned.
       *
       *************************************************************************
       */
       bool initMode(OpMode mode);

The second example shows how to document the function argument inline.

.. code-block:: cpp

      /*!
       ************************************************************************
       *
       * @brief Initialize a Foo object to given operation mode.
       *
       * The "read" mode means one thing, while "write" mode means another.
       *
       * @return bool value indicating success or failure of initialization.
       *             Success returns true, failure returns false.
       *
       *************************************************************************
       */
       bool initMode(OpMode mode /*!< [in] ReadMode, WriteMode are valid options */ );

Note that the first example uses the "\\" character to identify Doxygen 
commands; the second uses "@". 


Brief and detailed comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.27 A brief statement of purpose for a function must appear as the first 
text comment after the Doxygen command "\\brief" (or "@brief"). 

7.28 Any detailed function description, when included, **must** appear 
after the brief comment and **must** be separated from the brief comment by 
a line containing no text.


Return values
^^^^^^^^^^^^^

7.29 If the function has a non-void return type, the return value **should** 
be documented in the prologue using the Doxygen command "\\return" 
(or "@return") preceding a description of the return value. 

      Functions with "void" return type and C++ class constructors and 
      destructors **should not** have such documentation.

Arguments
^^^^^^^^^^^^^^^^^^^

7.30 Function arguments **should** be documented in the function prologue 
or inline (as shown above) when the intent or usage of the arguments is not 
obvious. 

      The inline form of the comment may be preferable when the argument 
      documentation is short. When a longer description is provided (such as 
      when noting the range of valid values, error conditions, etc.) the 
      description **should** be placed within the function prologue for 
      readability. However, the two alternatives for documenting function 
      arguments **must not** be mixed within the documentation of a single 
      function to reduce confusion. 

      In any case, superfluous documentation should be avoided. For example, 
      when there are one or two arguments and their meaning is obvious from 
      their names or the description of the function, providing no comments is 
      better than cluttering the code by documenting the obvious. Comments 
      that impart no useful information are distracting and less helpful than 
      no comment at all.

7.31 When a function argument is documented in the prologue comment section, 
the Doxygen command "\param" **should** appear before the comment as in the 
first example above.

7.32 The "in/out" status of each function argument **should** be documented.

       The Doxygen "\param" command supports this directly by allowing such an
       attribute to be specified as "\param[in]", "\param[out]", or 
       "\param[in,out]". Although the inline comment form does not support 
       this, such a description **should** be included; e.g., by using "[in]", 
       "[out]", or "[in,out]" in the comment.


Grouping small functions
^^^^^^^^^^^^^^^^^^^^^^^^^

7.33 Short, simple functions (e.g., inline methods) **may** be grouped 
together and documented with a single descriptive comment when this is 
sufficient.

      An example of Doxygen syntax for such a grouping is::

         ///@{
         //! @name Setters for data members

         void setMember1(int arg1) { m_member1 = arg1; }
         void setMember2(int arg2) { m_member2 = arg2; }

         ///@}


Header file vs. source file documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.34 Important implementation details (vs. usage detailed) about a function 
**should** be documented in the source file where the function is implemented,
rather than the header file where the function is declared.

      Header file documentation **should** include only purpose and usage 
      information germane to an interface. When a function has separate 
      implementation documentation, the comments **must not** contain Doxygen 
      syntax. Using Doxygen syntax to document an item in more than one location 
      (e.g., header file and source file) can cause undesired Doxygen 
      formatting issues and potentially confusing documentation.

      A member of a class may be documented as follows in the source file 
      for the class as follows (i.e., no Doxygen comments)::

        /*
         ***********************************************************************
         *
         * Set operation mode for a Foo object.
         *
         * Important detailed information about what the function does...
         *
         ***********************************************************************
         */
         bool Foo::initMode(OpMode mode)
         {
            ...function body...
         }


--------------------------------------------------------------------
Data member documentation
--------------------------------------------------------------------

7.35 Each struct field or class data member **should** have a descriptive 
comment indicating its purpose. 

     This comment may as appear as a prologue before the item, such as::

        /*!
         * \brief Brief statement describing the input mode...
         *
         * Optional detailed information about the input mode...
         */
        int m_input_mode;

     or, it may appear after the item as an inline comment such as::

        int m_input_mode; /*!< \brief Brief statement describing the input mode.... */


Brief and detailed comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.36 Regardless of which documentation form is used, a brief description 
**must** be included using the Doxygen command "\\brief" (or "@brief").

7.37 Any detailed description of an item, if included, **must** appear after 
the brief comment and be separated from the brief comment with a line
containing no documentation text.

     When a detailed comment is provided, or the brief statement requires 
     more than one line, the prologue comment form **should** be used instead 
     of the inline form to make the documentation easier to read.


Grouping data members
^^^^^^^^^^^^^^^^^^^^^^^^^

7.38 If the names of data members are sufficiently clear that their meaning 
and purpose are obvious to other developers (which should be determined in 
a code review), then the members **may** be grouped together and documented 
with a single descriptive comment.

      An example of Doxygen syntax for such a grouping is::

         ///@{
         //!  @name Data member description...

         int m_member1;
         int m_member2;
         ...
         ///@}


--------------------------------------------------------------------
Summary of common Doxygen commands
--------------------------------------------------------------------

This Section provides an overview of commonly used Doxygen commands.
Please see the `Doxygen guide <http://www.doxygen.nl/manual/>`_ for more details and information about other commands.

Note that to be processed properly, Doxygen commands **must** be preceded with 
either "\\" or "\@" character. For brevity, we use "\\" for all commands 
described here.

   \\brief 
     The "brief" command is used to begin a brief description of 
     a documented item. The brief description ends at the next blank line.
   \\file 
     The "file" command is used to document a file. Doxygen requires
     that to document any global item (function, typedef, enum, etc.), the file
     in which it is defined **must be** documented. 
   \\name
     The "name" command, followed by a name containing no blank 
     spaces, is used to define a name that can be referred to elsewhere 
     in the documentation (via a link).
   \\param 
     The "param" command documents a function parameter/argument.
     It is followed by the parameter name and description. The "\\param" 
     command can be given an optional attribute to indicate usage of the 
     function argument; possible values are "[in]", "[out]", and "[in,out]".
   \\return 
     The "return" command is used to describe the return value 
     of a function.
   \\sa 
     The "sa" command (i.e., "see also") is used to refer (and 
     provide a link to) another documented item. It is followed by the target 
     of the reference (e.g., class/struct name, function name, documentation 
     page, etc.).
   \@{,\@}  
     These two-character sequences begin and end a 
     grouping of documented items. Optionally, the group can be given a name 
     using the "name" command. Groups are useful for providing additional 
     organization in the documentation, and also when several items can be 
     documented with a single description (e.g., a set of simple, related 
     functions). 
   \\verbatim, \\endverbatim
     The "verbatim/endverbatim" commands are 
     used to start/stop a block of text that is to appear exactly as it is 
     typed, without additional formatting, in the generated documentation.
   -, -#
     The "-" and "-#" symbols begin an item in a bulleted 
     list or numbered list, respectively. In either case, the item ends at 
     the next blank line or next item.
   \\b, \\e 
     These symbols are used to make the next word bold or 
     emphasized/italicized, respectively, in the generated documentation.
   

