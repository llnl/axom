.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _codemacros-label:

=======================================================
10 Common Code Development Macros, Types, etc.
=======================================================

This section provides guidelines for consistent use of macros and types
defined in Axom components, such as "axom_utils" and "slic", and in our build 
system that we use in day-to-day code development.

.. important:: Code that is guarded with macros described in this section 
               **must not** change the *externally-observable* execution 
               behavior of the code.

               The macros are intended to help users and developers avoid
               unintended or potentially erroneous usage, etc. not confuse them. 

------------------------------------
Unused variables
------------------------------------

10.2 To silence compiler warnings and express variable usage intent more 
clearly, macros in the `AxomMacros.hpp` header file in the source include 
directory **must** be used when appropriate. For example,::

    void my_function(int x, int AXOM_DEBUG_PARAM(y))
    {
      // use variable y only for debug compilation
    }

Here, the `AXOM_DEBUG_PARAM` macro indicates that the variable 'y' is only
used when the code is compiled in debug mode. It also removes the variable
name in the argument list in non-debug compilation to prevent unwanted
compiler warnings.

Please see the `AxomMacros.hpp` header file for other available macros and 
usage examples.


------------------------------------
Disabling compiler-generated methods
------------------------------------

10.3 To disable compiler-generated class/struct methods when this is desired 
and to clearly express this intent, the `AXOMMacros.hpp` header file in 
source include directory contains macros that **should** be used for this 
purpose. See :ref:`compilergenmethods-label` for more information about 
compiler-generated methods.

Please see the `AXOMMacros.hpp` header file for other available macros and 
usage examples.


.. _codemacros-conditional-label:

------------------------------------
Conditionally compiled code
------------------------------------

10.4 Macros defined by Axom build system **must** be used to 
control conditional code compilation. 

For example, complex or multi-line code that is intended to be exposed only
for a debug build **must** be guarded using the `AXOM_DEBUG` macro::

   void MyMethod(...) 
   {
   #if defined(AXOM_DEBUG)
     // Code that performs debugging checks on object state, method args,
     // reports diagnostic messages, etc. goes here 
   #endif 

      // rest of method implementation
   }

Axom provides various other macro constants for conditionally-compiled code 
which reflect which built-in and third-party libraries are being used and 
which Axom components are enabled. The macro constants are defined in the 
``config.hpp.in`` file in the top-level Axom source directory. Each of these 
macro constants has the form ``AXOM_USE_<FOO>``, where ``FOO`` is the name of 
an Axom library dependency or the name of an Axom component.

When CMake is run to configure an Axom build, the macro constants are set 
based on CMake options and Axom dependencies and this file is converted to 
the ``config.hpp`` header file in the Axom build space. The ``config.hpp`` 
header file is included in all Axom source and header files for consistent 
application of the macro constants throughout the code.


------------------------------------
Error handling
------------------------------------

10.5 Macros provided in the "slic" component **should** be used to provide 
runtime checks for incorrect or questionable usage and informative messages 
for developers and users.

Runtime checks for incorrect or questionable usage and generation of 
informative warning, error, notice messages can be a tremendous help to 
users and developers. This is an excellent way to provide run-time debugging 
capabilities in code. Using the "slic" macros ensures that syntax and meaning
are consistent and that output information is handled similarly throughout 
the code. 

When certain conditions are encountered, the macros can emit failed boolean 
expressions and descriptive messages that help to understand potentially
problematic usage. Here's an example of common *SLIC* macro usage in AXOM::

   Bar* myCoolFunction(int in_val, Foo* in_foo)
   {
     if ( in_val < 0 || in_foo == nullptr )  
     {
       SLIC_CHECK_MSG( in_val >= 0, "in_val must be non-negative" );
       SLIC_CHECK( in_foo != nullptr );
       return nullptr;
     } else if ( !in_foo->isSet() ) {
       SLIC_CHECK_MSG( in_foo->isSet(), 
                       "in_foo is not set, will use default settings");
       const int val = in_val >= 0 ? in_val : DEFAULT_VAL;
       in_foo->setDefaults( val );
     }
    
     Bar* bar = new Bar(in_foo);

     return bar;
   }

This example uses slic macros that are only active when the code is compiled
in debug mode. When compiled in release mode, for example, the macros are 
empty and so do nothing. Also, when a condition is encountered that is 
problematic, such as 'in_val < 0' or 'in_foo == nullptr', the code will
emit the condition and an optional message and not halt. This allows calling
code to catch the issue (in this case a null return value) and react. There
are other macros (e.g., SLIC_ASSERT) that will halt the code if that is 
desired.

Slic macros operate in one of two compilation-defined modes. Some macros are 
active only in for a debug compile. Others are active for any build type.
Macros provided for each of these modes can be used to halt the code or not 
after describing the condition that triggered them. The following table
summarizes the SLIC macros.

============== ================ ====================
  Macro type     When active?     Halts code?
============== ================ ====================
  ERROR          Always           Yes
  WARNING        Always           No
  ASSERT         Debug only       Yes
  CHECK          Debug only       No
============== ================ ====================

Typically, we use macros ERROR/WARNING macros rarely. They are used primarily
to catch cases that are obvious programming errors or would put an application 
in a state where continuing is seriously in doubt. CHECK macros are used most
often, since they provide useful debugging information and do not halt the 
code -- they allow users to catch cases from which they can recover. ASSERT
macros are used in cases where halting the code is desired, but only in 
debug mode.

Please see the `slic.hpp` header file to see which macros are available and 
how to use them. 

.. important:: It is important to apply these macros judiciously so that they
               benefit users and other developers. We want to help folks use 
               our software correctly and not "spam" them with too much 
               information.

