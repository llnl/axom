.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _headerguide-label:

=====================================
4 Header File Organization
=====================================

The goal of these guidelines is to make it easy to find essential information 
in header files easily and quickly. Header files define software interfaces so
consistently-applied conventions for file organization can significantly 
improve user understanding and developer productivity. 

---------------------------------------------------------
All header file contents should be related
---------------------------------------------------------

4.1 A header file **may** contain multiple type definitions (e.g., structs, 
classes, enums, etc.). However, type definitions and function declarations in 
a header file **must** be related closely and/or support the primary type for 
which the file is named.

4.2 A "typedef" statement, when used, **should** appear in the header file 
where the type is defined. 

      This practice helps ensure that all names associated with a given type
      are available when the appropriate header file is used and eliminates
      potentially inconsistent type names.


-----------------------------------------------------------------------
Include in a header file only what's needed to compile it
-----------------------------------------------------------------------

4.3 A header file **must** be self-contained and self-sufficient.

    Specifically, each header file
    
      * **Must** have proper header file include guards 
        (see :ref:`headerlayout-label`) to prevent multiple inclusion. The 
        macro symbol name for each guard must be chosen to guarantee uniqueness 
        within *every* compilation unit in which it appears.
      * **Must** include all other headers and/or forward declarations it 
        needs to be compiled standalone. In addition, a file **should not** 
        rely on symbols defined in other header files it includes; the 
        other files **should** be included explicitly.
      * **Must** contain implementations of all generic templates and inline
        methods defined in it. A compiler will require the full definitions of
        these constructs to be seen in every source file that uses them.

.. note:: Function templates or class template members whose implementations 
          are fully specialized with all template arguments **must** be 
          defined in an associated source file to avoid linker errors 
          (e.g., multiply-defined symbols). Fully specialized templates 
          *are not* templates and are treated just like regular functions.

4.4 Extraneous header files or forward declarations (i.e., those not 
required for standalone compilation) **must not** be included in header files.

      Spurious header file inclusions, in particular, introduce spurious file
      dependencies, which can increase compilation time unnecessarily.


---------------------------------------------------------
Use forward declarations when you can
---------------------------------------------------------

4.5 Header files **should** use forward declarations instead of header file 
inclusions when possible. This may speed up compilation, especially when 
recompiling after header file changes.

.. note:: **Exceptions to this guideline:**

    * Header files that define external APIs for the Axom  
      project **must** include all header files for all types that 
      appear in the API. This makes use of the API much easier.
    
    * When using a function, such as an inline method or template, that 
      is implemented in a header file, the header file containing the
      implementation **must** be included.
    
    * When using C++ standard library types in a header file, it **may** be 
      preferable to include the actual headers rather than forward reference 
      headers, such as 'iosfwd', to make the header file easier to use. This 
      prevents users from having to explicitly include standard headers 
      wherever your header file is used.

4.6 A forward type declaration **must** be used in a header file when an 
include statement would result in a circular dependency among header files. 

.. note:: Forward references, or C++ standard 'fwd' headers, are preferred
          over header file inclusions when they are sufficient.


.. _headerincludeorder-label:

---------------------------------------------------------
Organize header file contents for easy understanding
---------------------------------------------------------

4.7 Header file include statements **should** use the same ordering pattern 
for all files.

      This improves code readability, helps to avoid misunderstood
      dependencies, and insures successful compilation regardless of
      dependencies in other files. A common, recommended header file 
      inclusion ordering scheme is (only some of these may be needed):

      #. Headers in the same Axom component
      #. Other headers within the project
      #. TPL headers; e.g., MPI, OpenMP, HDF5, etc.
      #. C++ and C standard library headers
      #. System headers

      Also, code is easier to understand when include files are ordered
      alphabetically within each of these sections and a blank line is
      inserted between sections. Adding comments that describe the
      header file categories can be helpful as well.  For example,

.. code-block:: cpp

         // Headers from this component
         #include "OtherClassInThisComponent.hpp"

         // "other" component headers
         #include "other/SomeOtherClass.hpp"

         // C standard library 
         #include <stdio.h>

         // C++ standard library
         #include <unordered_map>
         #include <vector>

         // Non-std system header
         #include <unistd.h>

.. note:: Ideally, header file inclusion ordering should not matter. 
          Inevitably, this will not always be the case. Following the
          ordering prescription above helps to avoid problems when others'
          header files are not constructed following best practices.


4.8 Routines **should** be ordered and grouped in a header file so that
code readability and understanding are enhanced.

      For example, all related methods should be grouped together. Also,
      public methods, which are part of an interface, should appear before 
      private methods.


---------------------------------------------------------
All function arguments should have names
---------------------------------------------------------

4.9 The name of each function argument **must** be specified in a header 
file declaration. Also, names in function declarations and definitions 
**must** match.

       For example, this is not an acceptable function declaration::

          void doSomething(int, int, int);

       Without argument names, the only way to tell what the arguments mean is
       to look at the implementation or hope that the method is documented 
       well.


.. _headerlayout-label:

---------------------------------------------------------
Header file layout details
---------------------------------------------------------

Content **must** be organized consistently in all header files. 
This section summarizes the recommended header file layout using numbers 
and text to illustrate the basic structure. Details about individual items 
are contained in the guidelines after the summary.

.. code-block:: cpp

   // (1) Axom copyright and release statement

   // (2) Doxygen file prologue

   // (3a) Header file include guard, e.g.,
   #ifndef MYCLASS_HPP
   #define MYCLASS_HPP

   // (4) Header file inclusions (when NEEDED in lieu of forward declarations)
   #include "myHeader.hpp"

   // (5) Forward declarations NEEDED in header file (outside of project namespace)
   class ForwardDeclaredClass;

   // (6a) Axom project namespace declaration
   namespace axom {

   // (7a) Internal namespace (if used); e.g.,
   namespace awesome {

   // (8) Forward declarations NEEDED in header file (in project namespace(s))
   class AnotherForwardDeclaredClass;

   // (9) Type definitions (class, enum, etc.) with Doxygen comments e.g.,
   /*!
    * \brief Brief ...summary comment text...
    *
    * ...detailed comment text...
    */
   class MyClass {
      int m_classMember;
   };

   // (7b) Internal namespace closing brace (if needed)
   } // awesome namespace closing brace

   // (6b) Project namespace closing brace
   } // axom namespace closing brace

   // (3b) Header file include guard closing endif */
   #endif // closing endif for header file include guard


4.10 **(Item 1)** Each header file **must** contain a comment section that 
includes the Axom copyright and release statement.

      See :ref:`docsec-label` for details.

4.11 **(Item 2)** Each header file **must** begin with a Doxygen file prologue.

      See :ref:`docsec-label` for details.

4.12 **(Items 3a,3b)** The contents of each header file **must** be guarded 
using a preprocessor directive that defines a unique "guard name" for the file.

      The guard must appear immediately after the file prologue and use the
      '#ifndef' directive (item 2a); this requires a closing '#endif' 
      statement at the end of the file (item 2b). 

      The preprocessor constant must use the file name followed by "_HPP" for
      C++ header files; e.g., "MYCLASS_HPP" as above.

      The preprocessor constant must use the file name followed by "_H" for
      C header files.

4.13 **(Item 4)** All necessary header file inclusion statements **must** 
appear immediately after copyright and release statement and before any 
forward declarations, type definitions, etc.

4.14 **(Item 5)** Any necessary forward declarations for types defined outside 
the project namespace **must** appear after the header include statements
and before the Axom project namespace statement.

4.15 **(Items 6a, 6b, 7a, 7b)** All types defined and methods defined in a 
header file **must** be included in a namespace.

      Either the project "axom" namespace (item 6a) or a namespace
      nested within the project namespace (item 7a) may be used, or 
      both may be used. A closing brace ( "}" ) is required to close each
      namespace declaration (items 6b and 7b) before the closing '#endif' 
      for the header file include guard.

4.16 **(Item 8)** Forward declarations needed **must** appear in the 
appropriate  namespace before any other statements (item 8).

4.17 **(Item 9)** All class and other type definitions **must** appear 
after header file inclusions and forward declarations. A proper class 
prologue **must** appear before the class definition. See :ref:`docsec-label`
for details.
