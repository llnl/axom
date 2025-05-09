.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _namesec-label:

===========================
2 Names
===========================

Good names are essential to sound software design. This section contains 
guidelines for naming files, types, functions, class members, variables, etc. 
The main goal is to use clear and unambiguous names. Also, we want naming 
conventions for different entities so that, when applied, the role of each 
is obvious from the form of its name.

------------------------------------
Good names are clear and meaningful
------------------------------------

2.1 Every name **must** be meaningful. In particular, its meaning **must** 
be clear to other code developers and users, not just the author of the name.

      A substantial benefit of good name selection is that it can greatly
      reduce the amount of developer debate to define a concept. A good name
      also tends to reduce the amount of documentation required for others to
      understand it. For example, when the name of a function clearly indicates
      what it does and the meaning and purpose of each argument is clear from
      its name, then code comments may be unnecessary. Documentation 
      can be a substantial part of software and requires maintenance. 
      Minimizing the amount of required documentation reduces this burden.

------------------------------------
Avoid cryptic names
------------------------------------

2.2 Tersely abbreviated or cryptic names **should** be avoided. However, 
common acronyms and jargon that are well understood by team members and
users **may** be used.

------------------------------------
Use terminology consistently
------------------------------------

2.3 Terminology **must** be used consistently; i.e., for names and concepts
in the code and in documentation. Multiple terms **should not** be used to 
refer to the same concept and a concept **should not** be referred to by 
multiple terms.

      Using a clear, limited set of terminology in a software project helps
      maintain the consistency and integrity of the software, and it makes
      the code easier to understand for developers and users.

2.4 Each name **must** be consistent with other similar names in the code.

      For example, if getter/setter methods follow the convention "getFoo"
      and "setFoo" respectively, then adding a new setter method called
      "putBar" is clearly inconsistent.


----------------------------------------------------------
Name directories so it's easy to know what's in them
----------------------------------------------------------

2.5 Each directory **must** be named so that the collective purpose 
of the files it contains is clear. All directory names **should** follow
the same style conventions. 

      All directory names **should** use all lower case letters and consist 
      of a single word in most cases. A directory name with more than one 
      word **should** use an 'underscore' to separate words. 

      For example, use::

         cool_stuff

      not ::

         cool-stuff


--------------------------------------------------------
Follow file extension conventions
--------------------------------------------------------

2.6 C++ header and source file extensions **must** be: \*.hpp and \*.cpp, 
respectively.

2.7 C header and source files (e.g., tests, examples, and generated API code)
**must** have extensions \*.h and \*.c, respectively. 

2.8 Fortran source files (e.g., tests and examples, and generated API code) 
**must** have the extension \*.f or \*.F . \*.F **must** be used if the 
preprocessor is needed to compile the source file.


---------------------------------------------------------
Associated source and header file names should match
---------------------------------------------------------

2.9 The names of associated header and source files **should** match, apart from
the file extension, to make their association clear.

      For example, the header and source files for a class "Foo" **should** 
      be named "Foo.hpp" and "Foo.cpp", respectively.

      Also, files that are closely related in other ways, such as a header file
      containing prototypes for a set of methods that are not class members and 
      a source file containing implementations of those methods, **should** be 
      named the same or sufficiently similar so that their relationship is 
      clear. 
 

-------------------------------------------------
File contents should be clear from file name
-------------------------------------------------

2.10 The name of each file **must** clearly indicate its contents.

      For example, the header and source file containing the definition and
      implementation of a major type, such as a class **must** include the 
      type name of the type in the file name. For example, the header and
      implementation file for a class called "MyClass" should be named 
      "MyClass.hpp" and "MyClass.cpp", respectively.

      Files that are not associated with a single type, but which contain 
      closely related functionality or concepts, **must** be named so that
      the functionality or concepts are clear from the name. For example,
      files that define and implement methods that handle file I/O **should** 
      be named "FileIO.hpp" and "FileUtils.cpp", or similar.


-------------------------------------------
File names should not differ only by case 
-------------------------------------------

2.11 File names that differ only in letter case **must not** be used.

      Since we must support Windows platforms, which have limited case
      sensitivity for file names, having files with names "MyClass.hpp" 
      and "myclass.hpp", for example, is not acceptable. 



------------------------
Namespace name format
------------------------

2.12 All namespaces defined **must** use all lowercase letters for 
consistency and to avoid user confusion.


--------------------------
Type name format
--------------------------

2.13 Type names (i.e., classes, structs, typedefs, enums, etc.) **must** be 
nouns and **should** be in mixed case with each word starting with 
an upper case letter and all other letters in lower cases.

      For example, these are preferred type names::

         DataStore, MyCollection, TypeUtils

      These type names should not be used::

         dataStore, mycollection, TYPEUTILS

2.14 Separating characters, such as underscores, **should not** be used 
between words in a type name.

      For example, these names are not preferred type names::

         Data_store, My_Collection

.. note:: **Exceptions to the guidelines above** include cases where types
          play a similar role to those in common use elsewhere. For example, 
          naming an iterator class "base_iterator" would be acceptable if 
          it is conceptually similar with the C++ standard library class.

2.15 Suffixes that may be used by compilers for name mangling, or 
which are used in the C++ standard library, such as "\_t", **must not** 
be used in type names.


-------------------------------------------------------
Function name format
-------------------------------------------------------

2.16 Function names **must** use "camelCase" or "pot_hole" style. camelCase 
is preferred. 

      **camelCase style:** The first word has all lower case letters.
      If multiple words are used, each word after the first starts with
      an upper case letter and all other letters in the word are lower case.
      Underscores must not be used in camelCase names, but numbers may be used.

      For example, these are proper camelCase names::

         getLength(), createView2()

      **pot_hole style:** All letters are lower case. If multiple
      words are used, they are separated by a single underscore. Numbers
      may be used in pothole style names.

      For example, these are acceptable pothole style variable names::

         push_front(), push_back_2()

2.17 Names of related functions, such as methods for a class, **should** 
follow the same style.

.. note:: **Exception:**  While consistency is important, name style may be 
          mixed when it makes sense to do so. While camelCase style is 
          preferred for class member functions, a class may also contain 
          methods that follow pot_hole style if those methods perform 
          operations that are similar to C++ standard library functions, 
          for example.

          For example, the following method names are acceptable for a class 
          with camelCase style names::
          
              push_back(), push_front()

          if those methods are similar in behavior to C++ standard methods.


-------------------------------------------------------
Function names should indicate behavior
-------------------------------------------------------

2.18 Each function name **must** indicate clearly indicate what the 
function does. 

      For example::

        calculateDensity(), getDensity()

      are good function names because they distinguish the fact that the
      first performs a calculation and the second returns a value. If a
      function were named::

        density()

      what it actually does is murky; i.e., folks would have to read its 
      documentation or look at its implementation to see what it actually does.

2.19 Function names **should** begin with a verb because they perform an action.

2.20 Verbs such as "is", "has", "can", etc. **should** be used for functions 
with a boolean return type.

      For example, the following names are preferred::

         isInitialized(), isAllocated()


-------------------------------------------------------
Related functions should have similar names
-------------------------------------------------------

2.21 Complementary verbs such as  "get/set", "add/remove" and "create/destroy"
**must** be used for routines that perform complementary operations.

      Such symmetry prevents confusion and makes interfaces easier to use.


-------------------------------------------
Data member and variable name format
-------------------------------------------

2.22 All variables (class/struct members, function-scoped variables, function
arguments, etc.) **must** use either "camelCase" style or "pot_hole" style. 
Pot_hole style is preferred since it distinguishes variable names from 
method names.

       For example, these are acceptable variable names::

         myAverage, person_name, pressure2

2.23 Non-static class and struct data member names **must** have the 
prefix "m\_".

      This convention makes it obvious which variables are class 
      members/struct fields and which are other local variables. For 
      example, the following are acceptable names for class data members using
      camelCase style::

         m_myAverage, m_personName

      and acceptable pothole style::

         m_my_average, m_person_name

2.24 Static class/struct data member names and static file scope variables
**must** have the prefix "s\_".

      Similar to the guideline above, this makes it obvious that the variable
      is static.


-------------------------------------------
Variable names should indicate type
-------------------------------------------

2.25 Verbs, such as "is", "has", "can", etc., **should** be used for boolean 
variables (i.e., either type bool or integer that indicates true/false).

      For example, these names are preferred::

         m_is_initialized, has_license

      to these names::

         m_initialized, license

2.26 A variable that refers to a non-fundamental type **should** give an 
indication of its type.

      For example,::

         Topic* my_topic;

      is clearer than::

         Topic* my_value;


------------------------------------
Macro and enumeration name format
------------------------------------

2.27 Preprocessor macro constants **must** be named using all uppercase 
letters and underscores **should** be used between words.

      For example, these are acceptable macro names::

         MAX_ITERATIONS, READ_MODE

      These are not acceptable::

         maxiterations, readMode

2.28 The name of each enumeration value **should** start with a capital letter
and use an underscore between words when multiple words are used.

       For example,::

          enum Orange
          {
             Navel,
             Valencia,
             Num_Orange_Types
          };
