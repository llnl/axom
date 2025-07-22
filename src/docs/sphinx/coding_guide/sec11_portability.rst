.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _portsec-label: 

===================================================
11 Portability, Compilation, and Dependencies
===================================================

C++ is a huge language with many advanced and powerful features. To avoid
over-indulgence and obfuscation, we would like to avoid C++ feature bloat.
By constraining or even banning the use of certain language features and
libraries, we hope to keep our code simple and portable. We also hope to 
avoid errors and problems that may occur when language features are not 
completely understood or not used consistently. This section lists such 
restrictions and explains why use of certain features is constrained or 
restricted.


--------------------------------------------------------------------
Portability
--------------------------------------------------------------------

Nothing beyond C++17
^^^^^^^^^^^^^^^^^^^^

11.1 C++ language features beyond standard C++17 **must not** be used unless
reviewed by the team and verified that the features are supported by all 
compilers we need to support.

      Changing this guideline requires full consensus of all team members.


No non-standard language constructs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

11.2 Special non-standard language constructs, such as GNU extensions, 
**must not** be used if they hinder portability.


.. note:: Any deviation from these C++ usage requirements must be 
          agreed on by all members of the team and vetted with our
          main application users.


--------------------------------------------------------------------
Compilation
--------------------------------------------------------------------

Avoid conditional compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

11.3 Excessive use of the preprocessor for conditional compilation at a 
fine granularity (e.g., selectively including or removing individual source 
lines) **should** be avoided. 

      While it may seem convenient, this practice typically produces confusing 
      and error-prone code. Often, it is better to refactor the code into 
      separate routines or large code blocks subject to conditional compilation
      where it is more clear. 

      Code reviews by team members will dictate what is/is not acceptable.


The compiler is your friend
^^^^^^^^^^^^^^^^^^^^^^^^^^^

11.4 Developers **should** rely on compile-time and link-time errors to 
check for code correctness and invariants. 

      Errors that occur at run-time and which depend on specific control flow 
      and program state are inconvenient for users and can be difficult to 
      detect and fix.

.. important::  Add specific guidance on how this should be done...


--------------------------------------------------------------------
Minimize dependencies on third-party libraries (TPLs)
--------------------------------------------------------------------

11.5 While it is generally desirable to avoid recreating functionality that
others have already implemented, we **should** limit third-party library
dependencies for Axom to make it easier for users. We are a library,
and everything we necessarily depend on will become a dependency for our
user.  

      **Before introducing any significant TPL dependency on Axom
      (e.g., hdf5), it must be agreed on by the development team and vetted
      with our main users.**

11.6 Unless absolutely necessary, any TPL we depend on **must not** be 
exposed through any public interface in Axom.
