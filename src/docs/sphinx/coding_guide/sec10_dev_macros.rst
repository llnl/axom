.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
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
clearly, macros in the ``axom/core/Macros.hpp`` header file **must** be used
when appropriate. For example,::

    void my_function(AXOM_MAYBE_UNUSED int y)
    {
      SLIC_ASSERT(y > 0);
    }

Here, ``AXOM_MAYBE_UNUSED`` indicates that the parameter ``y`` might be unused
in configurations where ``SLIC_ASSERT`` is disabled. The parameter remains
named in every configuration, while the annotation prevents unwanted compiler
warnings. Prefer this macro for variables and parameters whose use depends on
debug or SLIC macros.

Please see the ``axom/core/Macros.hpp`` header file for other available macros
and usage examples.


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

The ``AXOM_DEBUG_DEFINE`` CMake setting controls when ``AXOM_DEBUG`` is
defined. Its supported values are:

* ``DEFAULT``: define it for ``Debug`` and ``RelWithDebInfo`` configurations.
* ``ON``: define it for every configuration.
* ``OFF``: do not define it in any configuration.

Thus, ``AXOM_DEBUG`` normally follows the build configuration, but
``AXOM_DEBUG_DEFINE`` can explicitly override that behavior.

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

This example uses Slic macros controlled by the
``AXOM_ENABLE_SLIC_DEBUG_MACROS`` CMake setting. ``DEFAULT`` follows
``AXOM_DEBUG``: the macros are active when the ``AXOM_DEBUG`` compiler
definition is visible and are empty otherwise. With
``AXOM_DEBUG_DEFINE=DEFAULT``, Axom defines
``AXOM_DEBUG`` in ``Debug`` and ``RelWithDebInfo`` configurations and does not
define it in ``Release`` or ``MinSizeRel``. An installed Axom can export that
definition to downstream targets, independently of the downstream build type.
``ON`` or ``OFF`` explicitly enables or disables the Slic macros without
changing other code guarded by ``AXOM_DEBUG``. A downstream CMake project can
set this variable before calling ``find_package(axom)`` to override the setting
exported by Axom.

Choose a preprocessor guard according to the behavior it controls:

* Use ``AXOM_DEBUG`` for general debug-only behavior that is independent of
  whether Slic debug macros are active.
* Use ``AXOM_SLIC_DEBUG_MACROS_ENABLED`` for work performed exclusively to
  support ``SLIC_ASSERT``, ``SLIC_CHECK``, or ``SLIC_DEBUG``. Enabling Slic
  debug macros must not enable unrelated ``AXOM_DEBUG`` behavior.

For example, guard supporting work used only by Slic debug macros as follows::

   #if AXOM_SLIC_DEBUG_MACROS_ENABLED
     // Compute values or perform validation used only by Slic debug macros.
   #endif

This read-only macro evaluates to ``0`` or ``1`` and reflects the effective
state after applying ``AXOM_ENABLE_SLIC_DEBUG_MACROS`` and its default
``AXOM_DEBUG`` behavior. Do not define it directly. In particular, do not use
``AXOM_ENABLE_SLIC_DEBUG_MACROS`` as a source-code condition because it is
undefined in ``DEFAULT`` mode.

When disabled, the Slic debug macros are empty and do nothing. When a
condition is encountered that is problematic, such as 'in_val < 0' or
'in_foo == nullptr', the code will
emit the condition and an optional message and not halt. This allows calling
code to catch the issue (in this case a null return value) and react. There
are other macros (e.g., SLIC_ASSERT) that will halt the code if that is 
desired.

Slic macros operate in one of two compilation-defined modes. Some macros are
controlled by the debug-macro settings described above. Others are active for
any build type.
Macros provided for each of these modes can be used to halt the code or not 
after describing the condition that triggered them. The following table
summarizes the SLIC macros.

============== ========================================== ====================
  Macro type     When active?                               Halts code?
============== ========================================== ====================
  ERROR          Always                                     Yes
  WARNING        Always                                     No
  ASSERT         When Slic debug macros are enabled         Yes
  CHECK          When Slic debug macros are enabled         No
  DEBUG          When Slic debug macros are enabled         No
============== ========================================== ====================

Typically, we use macros ERROR/WARNING macros rarely. They are used primarily
to catch cases that are obvious programming errors or would put an application 
in a state where continuing is seriously in doubt. CHECK macros are used most
often, since they provide useful debugging information and do not halt the 
code -- they allow users to catch cases from which they can recover. ASSERT
macros are used in cases where halting the code is desired, but only in 
builds where Slic debug macros are enabled.

Please see the `slic.hpp` header file to see which macros are available and 
how to use them. 

.. important:: It is important to apply these macros judiciously so that they
               benefit users and other developers. We want to help folks use 
               our software correctly and not "spam" them with too much 
               information.
