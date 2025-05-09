.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _scopesec-label:

=====================================
6 Scope
=====================================

---------------------------------------------------------
Use namespaces to avoid name collisions
---------------------------------------------------------

6.1 All Axom code **must** be included in the project namespace 
'axom'; e.g.,::

         namespace axom {
              // . . .
         }

.. note::We will change the top-level namespace at some point to shorten it.

6.2 Each Axom component **must** define its own unique namespace within
the "axom" namespace. All contents of each component **must** reside
within that namespace.


---------------------------------------------------------
Use namespaces to hide non-API code in header files
---------------------------------------------------------

6.3 Code that must be appear in header files (e.g., templates) that is not 
intended to be part of a public interface, such as helper classes/structs 
and methods, **should** be placed in an internal namespace. 

      Common names for such namespaces include 'internal' (for implementations
      used only internally) and 'detailed' (for types, etc. used only 
      internally). Any reasonable choice is acceptable; however, the choice 
      **must** be the same within each Axom component.

      Note that declaring helper classes/structs private within a class 
      definition is another good option. See :ref:`scopenestedclass-label`
      for details.


---------------------------------------------------------
Use 'unnamed' namespace for hiding code in source files
---------------------------------------------------------

6.4 Classes/structs and methods that are meant to be used only internally to a 
single source file **should** be placed in the 'unnamed' namespace to make
them invisible outside the file.

      This guarantees link-time name conflicts will not occur. For example::

         namespace {
            void myInternalFunction();
         }


---------------------------------------------------------
Apply the 'using directive' carefully
---------------------------------------------------------

6.5 The 'using directive' **must not** be used in any header file.

      Applying this directive in a header file leverages a bad decision to
      circumvent the namespace across every file that directly or indirectly
      includes that header file. 

.. note:: This guideline implies that each type name appearing in a header 
          file **must be fully-qualified** (i.e., using the namespace 
          identifier and scope operator) if it resides in a different 
          namespace than the contents of the file.

6.6 The 'using directive' **may** be used in a source file to avoid using a 
fully-qualified type name at each declaration. Using directives **must** appear
after all "#include" directives in a source file.

6.7 When only parts of a namespace are used in an implementation file, only 
those parts **should** be included with a using directive instead of the 
entire namespace contents.

      For example, if you only need the standard library vector container form
      the "std" namespace, it is preferable to use::

         using std::vector;

      rather than::

         using namespace std;


---------------------------------------------------------
Use access qualifiers to control class interfaces
---------------------------------------------------------

6.8 Class members **must** be declared in the following order: 

      #. "public"
      #. "protected"
      #. "private"

      That is, order members using these access qualifiers in terms of 
      "decreasing scope of visibility".

.. note:: Declaring methods before data members is preferred because methods 
           are more commonly considered part of a class interface. Also,
           separating methods and data into their own access qualified 
           sections usually helps make a class definition clearer.

6.9 Class data members **should** be "private". The choice to use "public" 
or "protected" data members **must** be scrutinized by other team members.

      Information hiding is an essential part of good software engineering 
      and private data is the best means for a class to preserve its 
      invariants. Specifically, a class should maintain control of how object 
      state can be modified to minimize side effects. In addition, restricting
      direct access to class data enforces encapsulation and facilitates 
      design changes through refactoring.

6.10 Guard protected/private access qualifiers when hiding methods that contain
kernel launches.

      When building algorithms that are templated for execution on the GPU,
      methods that instantiate kernels via `axom::for_all()` cannot be marked
      as protected or private when building for certain backends. In these
      situations, add conditional compilation around the access qualifiers.

      For example::

         class Algorithm
         {
         #if !defined(__CUDACC__)
         private:
         #endif
           void helperMethod()
           {
             axom::for_all<ExecSpace>(100,
               AXOM_LAMBDA(axom::IndexType index)
               {
                 // do something
               });
           }
         };

---------------------------------------------------------
Use 'friend' and 'static' rarely
---------------------------------------------------------

6.11 "Friend" declarations **should** be used rarely. When used, they 
**must** appear within the body of a class definition before any class 
member declarations. This helps make the friend relationship obvious.

      Note that placing "friend" declarations before the "public:" keyword 
      makes them private, which preserves encapsulation.

6.12 Static class members (methods or data) **must** be used rarely. In 
every case, their usage **should** be carefully reviewed by the team.

      When it is determined that a static member is needed, it **must** appear 
      first in the appropriate member section. Typically, static member 
      functions **should** be "public" and static data members **should** be 
      "private".


.. _scopenestedclass-label:

---------------------------------------------------------
Hide nested classes when possible
---------------------------------------------------------

6.13 Nested classes **should** be private unless they are part of the 
enclosing class interface.

      For example::

         class Outer
         {
            // ...
         private:
            class Inner
            {
               // ...
            };
         };

      When only the enclosing class uses a nested class, making it private
      does not pollute the enclosing scope needlessly. Furthermore, nested 
      classes may be forward declared within the enclosing class definition 
      and then defined in the implementation file of the enclosing class. 
      For example::

         class Outer
         {
            class Inner; // forward declaration

            // use name 'Inner' in Outer class definition
         };

         // In Outer.cpp implementation file...
         class Outer::Inner
         {
            // Inner class definition
         }

      This makes it clear that the nested class is only needed in the
      implementation and does not clutter the class definition.


---------------------------------------------------------
Limit scope of local variables
---------------------------------------------------------

6.14 Local variables **should** be declared in the narrowest scope possible 
and as close to first use as possible.

      Minimizing variable scope makes source code easier to comprehend and
      may have performance and other benefits. For example, declaring a loop 
      index inside a for-loop statement such as::

         for (int ii = 0; ...) {

      is preferable to::

         int ii;
         ...
         for (ii = 0; ...) {

      Beyond readability, this rule has benefits for thread safety, etc.

.. note:: **Exception:** When a local variable is an object, its constructor 
          and destructor may be invoked every time a scope (such as a loop) 
          is entered and exited, respectively. 

      Thus, instead of this::

         for (int ii = 0; ii < 1000000; ++ii) {
            Foo f;
            f.doSomethingCool(ii);
         }

      it may be more efficient to do this::

         Foo f;
         for (int ii = 0; ii < 1000000; ++ii) {
            f.doSomethingCool(ii);
         }

6.15 A local reference to any item in the global namespace (which should be 
rare if needed at all) **should** use the scope operator ("::") to make 
the fact that it resides in the global namespace clear.

      For example::

         int local_val = ::global_val;
