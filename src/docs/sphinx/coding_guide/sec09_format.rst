.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

=====================================
9 Code Formatting
=====================================

--------------------------------------------------------------------
Conditional statements and loops
--------------------------------------------------------------------

9.1 Curly braces **should** be used in all conditionals, loops, etc. 
even when the content inside the braces is a "one-liner". 

       This helps prevent coding errors and misinterpretation of intent. 
       For example, this::

          if (done) { ... }

       is preferable to this::

          if (done) ...

9.2 One-liners **may** be used for "if" conditionals with 
"else/else if"  clauses when the resulting code is clear. 

       For example, either of the following styles **may** be used::

          if (done) {
             id = 3;
          } else {
             id = 0;
          }

       or::

          if (done) { id = 3; } else { id = 0; }

9.3 Complex "if/else if" conditionals with many "else if" clauses 
**should** be avoided.

      Such statements can always be refactored using local boolean variables 
      or "switch" statements. Doing so often makes code easier to read and 
      understand and may improve performance.

9.4 An explicit test for zero/nonzero **must** be used in a conditional 
unless the tested quantity is a boolean or pointer type. 

      For example, a conditional based on an integer value should use::

         if (num_lines != 0) { ... }

      not::

         if (num_lines) { ... }


--------------------------------------------------------------------
White space and code alignment
--------------------------------------------------------------------

Most conventions for indentation, spacing and code alignment 
preferred by the team are enforced by using the `clang-format` tool. 

There are several build system targets related to code formatting grouped 
under the `check` and `style` targets. The former verify that the code is 
properly formatted, while the latter modify source files to conform to 
axom's rules.

.. important:: Axom's `style` targets modify source files. Please ensure
   that you've committed/staged all your changes before running them.

.. tip:: When axom is configured to use the `make`-based generator, the 
   entire codebase can be formatted by running ``make clangformat_style`` 
   from the build directory. Simiarly, one can verify that the code if 
   properly formatted by running ``make clangformat_check``. There are 
   also component-specific variants for these targets, e.g. for axom's
   `core` component, we have ``core_clangformat_style`` and 
   ``core_clangformat_check``.

Not all preferred formatting conventions are supported by `clang-format`.
The following guidelines provide additional recommendations to make
code easier to read and understand.


White space enhances code readability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

9.5 Blank lines and indentation **should** be used throughout code to 
enhance readability. 

      Examples of helpful white space include:

         * Between operands and operators in arithmetic expressions.
         * After reserved words, such as "while", "for", "if", "switch", etc. 
           and before the parenthesis or curly brace that follows.
         * After commas separating arguments in functions.
         * After semi-colons in for-loop expressions.
         * Before and after curly braces in almost all cases.

9.6 White space **must not** appear between a function name and the opening 
parenthesis to the argument list. In particular, if a function call is broken 
across source lines, the break **must not** come between the function name and 
the opening parenthesis.

9.7 Tabs **must not** be used for indentation since this can be problematic 
for developers with different text editor settings.


Vertical alignment helps to show scope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

9.8 When function arguments (in either a declaration or implementation)
appear on multiple lines, the arguments **should** be vertically aligned
for readability.

9.9 All statements within a function body **should** be indented within 
the surrounding curly braces.

9.10 All source lines in the same scope **should** be vertically aligned.
Continuation of previous lines **may** be indented if it make the code easier
to read.


Break lines where it makes sense
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

9.11 When a line is broken at a comma or semi-colon, it **should** be broken 
after the comma or semi-colon, not before. This helps make it clear that 
the statement continues on the next line.

9.12 When a source line is broken at an arithmetic operator 
(i.e., +, -, etc.), it **should** be broken after the operator, not before. 


Use parentheses for clarity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

9.13 Parentheses **should** be used in non-trivial mathematical and logical 
expressions to clearly indicate structure and intended order of operations. 
Do not assume everyone who reads the code knows all the rules for operator 
precedence.
