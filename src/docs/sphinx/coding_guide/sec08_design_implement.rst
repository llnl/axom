.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _designsec-label:

=======================================================
8 Design and Implement for Correctness and Robustness
=======================================================

The guidelines in this section describe various software design and
implementation practices that help enforce correctness and robustness 
and avoid mis-interpretation or confusion by others.


--------------------------------------------------------------------
Keep it simple...
--------------------------------------------------------------------

8.1 Simplicity, clarity, ease of modification and extension **should** 
always be a main goal when writing new code or changing existing code. 

8.2 Each entity (class, struct, variable, function, etc.) **should** embody 
one clear, well-defined concept. 

      The responsibilities of an entity may increase as it is used in new and 
      different ways. However, changes that divert it from its original intent 
      **should** be avoided. Also, large, monolithic entities that provide too 
      much functionality or which include too many concepts tend to increase 
      code coupling and complexity and introduce undesirable side effects. 
      Smaller, clearly constrained objects are easier to write, test, maintain,
      and use correctly. Also, small, simple objects tend to get used more 
      often and reduce code redundancy.


--------------------------------------------------------------------
Avoid global and static data
--------------------------------------------------------------------

8.3 Global, complex, or opaque data sharing **should not** be used. Shared 
data increases coupling and contention between different parts of a code base, 
which makes maintenance and modification difficult.

8.4 Static or global variables of class type **must not** be used.

      Due to indeterminate order of construction, their use may cause bugs
      that are very hard to find. Static or global variables that are pointers
      to class types **may** be used and must be initialized properly in a
      single source file.


--------------------------------------------------------------------
Avoid macros and magic numbers for constants
--------------------------------------------------------------------

8.5 Preprocessor macros **should not** be used when there is a better 
alternative, such as an inline function or a constant variable definition.

      For example, this::

         const double PI = 3.1415926535897932384626433832;

      is preferable to this::

         #define PI (3.1415926535897932384626433832)

      Macros circumvent the ability of a compiler to enforce beneficial
      language concepts such as scope and type safety. Macros are also
      context-specific and can produce errors that cannot be understood
      easily in a debugger. Macros **should be used only** when there is
      no better choice for a particular situation.

8.6 Hard-coded numerical constants and other "magic numbers" **must not** 
be used directly in source code. When such values are needed, they **should** 
be declared as named constants to enhance code readability and consistency.


.. _compilergenmethods-label:

------------------------------------------------------
Avoid issues with compiler-generated class methods
------------------------------------------------------

The guidelines in this section apply to class methods that may be 
*automatically generated* by a compiler, including constructors, destructors,
copy, and move methods. Developers should be aware of the conditions under
which compilers will and will not generate these methods. Developers should
also be aware of when compiler-generated methods suffice and when they do not.
After providing some guidelines, we discuss standard C++ rules that compilers
follow for generating class methods when they are not explicitly defined. 
See :ref:`automethods-label`.

The most important cases to pay attention to involve the destructor, copy
constructor, and copy-assignment operator. Classes that provide these methods,
either explicitly or compiler-generated, are referred to as *copyable*. Failing 
to follow the rules for these methods can be damaging due to errors or 
unexpected behavior. Rules involving the move constructor and move-assignment 
operator are less important since they mostly affect efficiency and not 
correctness. Copy operations can be used to accomplish the same end result
as move operations, just less efficiently. Move semantics are an important
optimization feature of C++. The C++11 standard requires compilers to use 
move operations instead of copy operations when certain conditions are 
fulfilled. Classes that provide move operations, either explicitly or 
compiler-generated, are referred to as *movable*.


Rule of three
^^^^^^^^^^^^^^

8.7 Each class **must** follow the *Rule of Three* which states: if the 
destructor, copy constructor, or copy-assignment operator is explicitly 
defined, then the others **must** be defined.

      Compiler-generated and explicit versions of these methods **must not**
      be mixed. If a class requires one of these methods to be implemented, 
      it almost certainly requires all three to be implemented. 

      This rule helps guarantee that class resources are managed properly. 
      C++ copies and copy-assigns objects of user-defined types in various 
      situations (e.g., passing/returning by value, container manipulations, 
      etc.). These special member functions will be called, if accessible. 
      If they are not user-defined, they are implicitly-defined by the compiler.

      Compiler-generated special member functions can be incorrect 
      if a class manages a resource whose handle is an object of 
      non-class type. Consider a class data member which is a bare pointer to 
      an object. The compiler-generated class destructor will not free the 
      object. Also, the compiler-generated copy constructor and copy-assignment
      operator will perform a "shallow copy"; i.e., they will copy the value 
      of the pointer without duplicating the underlying resource.

      Neglecting to free a pointer or perform a deep copy when those operations
      are expected can result in serious logic errors. Following the Rule of 
      Three guards against such errors. On the rare occasion these actions are 
      intentional, a programmer-written destructor, copy constructor, and 
      copy-assignment operator are ideal places to document intent of
      design decisions.


Restrict copying of non-copyable resources
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.8 A class that manages non-copyable resources through non-copyable handles, 
such as pointers, **should** declare the copy methods private and and leave 
them unimplemented.

      When the intent is that such methods should never be called, this is a 
      good way to help a compiler to catch unintended usage. For example::

	   class MyClass
	   {
	      // ...

	   private:
              DISABLE_DEFAULT_CTOR(MyClass);
              DISABLE_COPY_AND_ASSIGNMENT(MyClass);

	      // ...
	   };

      When code does not have access to the private members of a class tries 
      to use such a method, a compile-time error will result. If a class does 
      have private access and tries to use one of these methods an link-time 
      error will result. 

      This is another application of the "Rule of Three".

      Please see :ref:`codemacros-label` for more information about the 
      macros used in this example to disable compiler-generated methods.

.. note::  **Exception:** If a class inherits from a base class that declares
           these methods private, the subclass need not declare the methods
           private. Including comments in the derived class header indicating 
           that the the parent class enforces the non-copyable properties of 
           the class is helpful.


Rely on compiler-generated methods when appropriate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.9 When the compiler-generated methods are appropriate (i.e.,
correct and sufficiently fast), the default constructor, copy constructor, 
destructor, and copy assignment **may** be left undeclared. In this case, 
it is often helpful to add comments to the class header file indicating that 
the compiler-generated versions of these methods will be used.

8.10 If a class is default-constructable and has POD ("plain old data") or 
pointer data members, a default constructor **should** be provided explicitly 
and its data members **must** be initialized explicitly if a default 
constructor is provided. A compiler-generated default constructor will not 
initialize such members, in general, and so will leave a constructed object 
in an undefined state.

      For example, the following class should provide a default constructor
      and initialize its data members in it::

	   class MyClass
	   {
	      MyClass();

	      // ...

	   private:
              double* m_dvals;
              int[]   m_ivals;
              
	   };


Functors should always be copyable 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.11 By convention, a functor class **should** have a copy constructor and 
copy-assignment operator. 

      Typically, the compiler-generated versions are sufficient when the class 
      has no state or non-POD data members. Since such classes are usually 
      small and simple, the compiler-generated versions of these methods 
      **may** be used without documenting the use of default value semantics 
      in the functor definition.

      For example::

	   class MyFunctor
	   {
	      // Compiler-generated copy ctor and copy assignment sufficient 

	   private:
	      DIABLE_DEFAULT_CTOR(MyFunctor); // prevent default construction

	      // ...
	   };

Note that in this example, the default constructor is disabled to prevent
default construction. This can help prevent programming errors when 
object state must be fully initialialized on construction. For more 
information about common Axom macros, see :ref:`codemacros-label`.


.. _automethods-label:

--------------------------------------------------------
Understand standard rules for compiler-generated methods
--------------------------------------------------------

This section provides some background information related to the guidelines
in the previous section. There, we provide guidelines that help to decide 
when to define class methods that may be generated automatically by a compiler 
and when relying on compiler-generated versions suffices.  Here, we describe
the conditions under which compilers generate methods automatically.

Consider the following simple class::

   class MyClass
   {
   public:
      int x;
   };

How many methods does it have? None?

Actually, MyClass may have as many as **six** methods depending on how it is 
used: a default constructor, destructor, copy constructor, copy-assignment 
operator, move constructor, and move-assignment operator. Any of these may 
be generated by a compiler.

.. note:: See :ref:`portsec-label` for discussion about using C++11 features
          such as *move semantics*.

C++ compiler rules for generating class member functions are:

   * The parameter-less default constructor is generated if a class does
     not define *any* constructor and all base classes and data members
     are default-constructable. This means that once you declare a copy
     constructor (perhaps to disable the automatically provided one),
     the compiler will not supply a default constructor.
   * The destructor is automatically supplied if possible, based on the
     members and the base classes.
   * A copy constructor is generated if all base classes and members are
     copy-constructable. Note that reference members are copy-constructable.
   * The copy-assignment operator is generated if all base classes and members
     are copy-assignable. For this purpose, reference members are not
     considered copy-assignable.
   * A move constructor is supplied unless the class has any of the following: 
     a user-defined copy constructor, copy-assignment operator, 
     move-assignment operator, or destructor. If the move constructor cannot
     be implemented because not all base classes or members are
     move-constructable, the supplied move constructor will be defined
     as deleted.
   * A move-assignment operator is generated under the same conditions as 
     the move constructor.

The importance of understanding these rules and applying the guidelines in 
the previous section is underscored by the fact that compiler-generated 
methods may have different behaviors depending on how they are used. Here 
we provide some examples based on MyClass defined above.

If MyClass has a user-defined constructor, then

.. code-block:: cpp

    MyClass item1;

and

.. code-block:: cpp

    MyClass item2 = MyClass();

will both call the user-defined default constructor "MyClass()" and there is
only one behavior.

However, if MyClass relies on the compiler-generated constructor

.. code-block:: cpp

    MyClass item1;

performs *default initialization*, while

.. code-block:: cpp

    MyClass item2 = MyClass();

performs *value initialization*.

Default initialization calls the constructors of any base classes, and nothing
else. Since constructors for intrinsic types do not do anything, that means
all member variables will have garbage values; specifically, whatever values 
happen to reside in the corresponding addresses.

Value initialization also calls the constructors of any base classes. Then,
one of two things happens:

   * If MyClass is a POD class (all member variables are either intrinsic
     types or classes that only contain intrinsic types and have no
     user-defined constructor/destructor), all data is initialized to 0.
   * If MyClass is not a POD class, the constructor does not touch any data,
     which is the same as default initialization (so member variables have
     garbage values unless explicitly constructed otherwise).

Other points worth noting:

   * Intrinsic types, such as int, float, bool, pointers, etc. have
     constructors that do nothing (not even initialize to zero), destructors
     that do nothing, and copy constructors and copy assignment-ers that
     blindly copy bytes.
   * Comparison operators, such as "==" or "!=" are never automatically
     generated by a compiler, even if all base classes and members are
     comparable.


---------------------------------------------------
Initializing and copying class members
---------------------------------------------------

Initialize all members at construction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.12 Class type variables **should** be defined using direct initialization 
instead of copy initialization to avoid unwanted and spurious type conversions 
and constructor calls that may be generated by compilers.

      For example, use::

         std::string name("Bill");

      instead of::

         std::string name = "Bill";

      or::

         std::string name = std::string("Bill");

8.13 Each class data member **must** be initialized (using default values 
when appropriate) in every class constructor. That is, an initializer or
initialization **must** be provided for each class data member so that 
every object is in a well-defined state upon construction. 

      Generally, this requires a user-defined default constructor when a class 
      has POD members. Do not assume that a compiler-generated default 
      constructor will leave any member variable in a well-defined state.

.. note::  **Exception:** A class that has no data members, including one that
           is derived from a base class with a default constructor that provides 
           full member initialization, does not require a user-defined default 
           constructor since the compiler-generated version will suffice.


Know when to use initialization vs. assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.14 Data member initialization **should** be used instead of assignment in 
constructors, especially for small classes. Initialization prevents needless 
run-time work and is often faster.

8.15 When using initialization instead of assignment to set data member 
values in a constructor, data members **should** always be initialized 
in the order in which they appear in the class definition. 

      Compilers adhere to this order regardless of the order that members 
      appear in the class initialization list. So you may as well agree with 
      the compiler rules and avoid potential errors that could result when
      one member depends on the state of another.

8.16 For classes with complex data members, assignment within the body of 
the constructor **may** be preferable.

      If the initialization process is sufficiently complex, it **may** be
      better to initialize (i.e., assign) member objects in a method that 
      is called after object creation, such as "init()".


Use the copy-and-swap idiom
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.17 A user-supplied implementation of a class copy-assignment operator 
**should** check for assignment to self, **must** copy all data members 
from the object passed to operator, and **must** return a reference to "\*this".

      The *copy-and-swap* idiom **should** be used. 


Initializing, copying, and inheritance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.18 A constructor **must not** call a virtual function on any data member 
object since an overridden method defined in a subclass cannot be called 
until the object is fully constructed. 

      There is no general guarantee that data members are fully-created 
      before a constructor exits.

8.19 All constructors and copy operations for a derived class **must** call 
the necessary constructors and copy operations for each of its base classes 
to insure that each object is properly allocated and initialized.


---------------------------------------------------
Prefer composition to inheritance
---------------------------------------------------

8.20 Class composition **should** be used instead of inheritance to extend 
behavior.

      Looser coupling between objects is typically more flexible and easier
      to maintain and refactor.


---------------------------------------------------
Keep inheritance relationships simple
---------------------------------------------------

8.21 Class hierarchies **should** be designed so that subclasses inherit 
from abstract interfaces; i.e., pure virtual base classes.

      Inheritance is often done to reuse code that exists in a base class.
      However, there are usually better design choices to achieve reuse.
      Good object-oriented use of inheritance is to reuse existing *calling*
      code by exploiting base class interfaces using polymorphism. Put another
      way, "interface inheritance" should be used instead of "implementation
      inheritance".

8.22 Deep inheritance hierarchies; i.e., more than 2 or 3 levels, **should**
be avoided.

8.23 Multiple inheritance **should** be restricted so that only one base 
class contains methods that are not "pure virtual".

8.24 "Private" and "protected" inheritance **must not** be used unless you 
absolutely understand the ramifications of such a choice and are sure that 
it will not create design and implementation problems.

      Such a choice **must** be reviewed with team members. There almost
      always exist better alternatives.


---------------------------------------------------
Design for/against inheritance
---------------------------------------------------

8.25 One **should not** inherit from a class that was not designed to be a 
base class; e.g., if it does not have a virtual destructor.

      Doing so is bad practice and can cause problems that may not be reported 
      by a compiler; e.g., hiding base class members. To add functionality, 
      one **should** employ class composition rather than by "tweaking" an 
      existing class.

8.26 The destructor of a class that is designed to be a base class **must** 
be declared "virtual". 

      However, sometimes a destructor should not be declared virtual, such as 
      when deletion through a pointer to a base class object should be 
      disallowed.


---------------------------------------------------
Use virtual functions responsibly
---------------------------------------------------

8.27 Virtual functions **should** be overridden responsibly. That is, the 
pre- and post-conditions, default arguments, etc. of the virtual functions 
should be preserved.

      Also, the behavior of an overridden virtual function **should not**
      deviate from the intent of the base class. Remember that derived classes
      are subsets, not supersets, of their base classes.

8.28 Inherited non-virtual methods **must not** be overloaded or hidden.

8.29 A virtual function in a base class **should only** be implemented in
the base class if its behavior is always valid default behavior for *any* 
derived class.

8.30 If a method in a base class is not expected to be overridden in any 
derived class, then the method **should not** be declared virtual.

8.31 If each derived class has to provide specific behavior for a base class 
virtual function, then it **should** be declared *pure virtual*.

8.32 Virtual functions **must not** be called in a class constructor or 
destructor. Doing so is undefined behavior. Even if it seems to work 
correctly, it is fragile and potentially non-portable.


--------------------------------------------------------------------
Inline functions
--------------------------------------------------------------------

Function inlining is a compile time operation and the full definition of an 
inline function must be seen wherever it is called. Thus, the implementation
of every function to be inlined must be provided in a header file. 

Whether or not a function implemented in a header file is explicitly declared
inline using the "inline" keyword, the compiler decides if the function will 
be inlined. A compiler will not inline a function that it considers too 
long or too complex (e.g., if it contains complicated conditional logic). 
When a compiler inlines a function, it replaces the function call with the 
body of the function. Most modern compilers do a good job of deciding when 
inlining is a good choice.

It is possible to specify function attributes and compiler flags that can
force a compiler to inline a function. Such options should be applied with 
care to prevent excessive inlining that may cause executable code bloat and/or 
may make debugging difficult.

.. note:: **When in doubt, don't use the "inline" keyword and let the compiler 
          decide whether to inline a function.**


Inline short, simple functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.33 Simple, short frequently called functions, such as accessors, that will
almost certainly be inlined by most compilers **should** be implemented inline 
in header files.


Only inline a class constructor when it makes sense
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.34 Class constructors **should not** be inlined in most cases.

      A class constructor implicitly calls the constructors for its base 
      classes and initializes some or all of its data members, potentially 
      calling more constructors. If a constructor is inlined, the construction 
      and initialization needed for its members and bases will appear at every 
      object declaration.

.. note::  **Exception:** A class/struct that has only POD members, is not 
           a subclass, and does not explicitly declare a destructor, can 
           have its constructor safely inlined in most cases.


Do not inline virtual methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.35 Virtual functions **must not** be inlined due to polymorphism. 

      For example, do not declare a virtual class member function as::

         inline virtual void foo( ) { }

      In most circumstances, a virtual method cannot be inlined because a
      compiler must do runtime dispatch on a virtual method when it doesn't 
      know the complete type at compile time.

.. note:: **Exception:** It is safe to define an empty destructor inline in an
          abstract base class with no data members.

.. important:: Should we add something about C++11 'final' keyword???


--------------------------------------------------------------------
Function and operator overloading
--------------------------------------------------------------------

There's a fine line between clever and...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.36 Operator overloading **must not** be used to be clever to the point of 
obfuscation and cause others to think too hard about an operation. 
Specifically, an overloaded operator must preserve "natural" semantics 
by appealing to common conventions and **must** have meaning similar 
to non-overloaded operators of the same name.

      Overloading operators can be beneficial, but **should not** be overused 
      or abused. Operator overloading is essentially "syntactic sugar" and an
      overloaded operator is just a function like any other function. An 
      important benefit of overloading is that it often allows more 
      appropriate syntax that more easily communicates the meaning of an 
      operation. The resulting code can be easier to write, maintain, and 
      understand, and it may be more efficient since it may allow the compiler
      to take advantage of longer expressions than it could otherwise.


Overload consistently
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.37 Function overloading **must not** be used to define functions that 
do conceptually different things. 

      Someone reading declarations of overloaded functions should be able to 
      assume (and rightfully so!) that functions with the same name do 
      something very similar.

8.38 If an overloaded virtual method in a base class is overridden in a 
derived class, all overloaded methods with the same name in the base class 
**must** be overridden in the derived class. 

      This prevents unexpected behavior when calling such member functions. 
      Remember that when a virtual function is overridden, the overloads of 
      that function in the base class **are not visible** to the derived class.


Common operators
^^^^^^^^^^^^^^^^^

8.39 Both boolean operators "==" and "!=" **should** be implemented if one 
of them is. 

      For consistency and correctness, the "!=" operator **should** be 
      implemented using the "==" operator implementation. For example::

         bool MyClass::operator!= (const MyClass& rhs)
         {
            return !(this == rhs);
         }

8.40 Standard operators, such as "&&", "||", and "," (i.e., comma), 
**must not** be overloaded.

      Built-in versions of these operators are typically treated specially 
      by a compiler. Thus, programmers cannot implement their full semantics. 
      This can cause confusion. For example, the order of operand evaluation 
      cannot be guaranteed when overloading operators "&&" or "||". This may 
      cause problems as someone may write code that assumes that evaluation 
      order is the same as the built-in versions.


--------------------------------------
Function arguments
--------------------------------------

Consistent argument order makes interfaces easier to use
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.41 Function arguments **must** be ordered similarly for all routines 
in an Axom component.

      Common conventions are either to put all input arguments first, then
      outputs, or vice versa. Input and output arguments **must not** be mixed 
      in a function signature. Parameters that are both input and output can 
      make the best choice unclear. Conventions consistent with related 
      functions **must** always be followed. When adding a new parameter to an 
      existing method, the established ordering convention **must** be followed.

.. note:: When adding an argument to an existing method, do not just stick it
          at the end of the argument list.


Pointer and reference arguments and const
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.42 Each function argument that is not a built-in type (i.e., int, double, 
char, etc.) **should** be passed either by reference or as a pointer to avoid 
unnecessary copies.

8.43 Each function reference or pointer argument that is not changed by
the function **must** be declared "const".


Always name function arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.44 Each argument in a function declaration **must** be given a name that 
exactly matches the function implementation.

      For example, use::

         void computeSomething(int op_count, int mode);

      not::

         void computeSomething(int, int);


--------------------------------------
Function return points
--------------------------------------

8.45 Each function **should** have exactly one return point to make 
control logic clear.

      Functions with multiple return points tend to be a source of errors when
      trying to understand or modify code, especially if there are multiple 
      return points within a scope. Such code can always be refactored to 
      have a single return point by using local scope boolean variables and/or 
      different control logic.

      A function **may** have two return points if the first return statement
      is associated with error condition check, for example. In this case,
      the error check **should** be performed at the start of the function body
      before other statements are reached. For example, the following is a
      reasonable use of two function return points because the error condition
      check and the return value for successful completion are clearly visible::

         int computeSomething(int in_val)
         {
            if (in_val < 0) { return -1; }

            // ...rest of function implementation...

            return 0;
         }

.. note:: **Exception.** If multiple return points actually fit well into the
          logical structure of some code, they **may** be used. For example, 
          a routine may contain extended if/else conditional logic with 
          several "if-else" clauses. If needed, the code may be more clear if
          each clause contains a return point.


--------------------
Proper type usage
--------------------

8.46 The "bool" type **should** be used instead of "int" for boolean 
true/false values.

8.47 The "string" type **should** be used instead of "char\*".

      The string type supports and optimizes many character string manipulation
      operations which can be error-prone and less efficient if implemented
      explicitly using "char\*" and standard C library functions. Note that
      "string" and "char\*" types are easily interchangeable, which allows C++
      string data to be used when interacting with C routines.

8.48 An enumeration type **should** be used instead of macro definitions 
or "int" data for sets of related constant values. 

      Since C++ enums are distinct types with a compile-time specified set of 
      values, there values cannot be implicitly cast to integers or 
      vice versa -- a "static_cast" operator must be used to make the 
      conversion explicit. Thus, enums provide type and value safety and 
      scoping benefits.

      In many cases, the C++11 `enum class` construct **should** be used 
      since it provides stronger type safety and better scoping than regular
      enum types.


---------------
Templates
---------------

8.49 A class or function **should** only be made a template when its 
implementation is independent of the template type parameter.

       Note that class member templates (e.g., member functions that are
       templates of a class that is not a template) are often useful to
       reduce code redundancy.

8.50 Generic templates that have external linkage **must** be defined in the 
header file where they are declared since template instantiation is a compile 
time operation. Implementations of class templates and member templates that
are non-trivial **should** be placed in the class header file after the class 
definition.


--------------------------------------------------------------------
Use const to enforce correct usage
--------------------------------------------------------------------

8.51 The "const" qualifier **should** be used for variables and methods 
when appropriate to clearly indicate usage and to take advantage of 
compiler-based error-checking. For example, any class member function 
that does not change the state of the object on which it is called 
**should** be declared "const"

      Constant declarations can make code safer and less error-prone since they 
      enforce intent at compile time. They also improve code understanding
      because a constant declaration clearly indicates that the state
      of a variable or object will not change in the scope in which the 
      declaration appears.

8.52 Any class member function that does not change a data member of the 
associated class **must** be declared "const".

      This enables the compiler to detect unintended usage.

8.53 Any class member function that returns a class data member that 
should not be changed by the caller **must** be declared "const" and 
**must** return the data member as a "const" reference or pointer.

       Often, both "const" and non-"const" versions of member access functions
       are needed so that callers may declare the variable that holds the
       return value with the appropriate "const-ness".


--------------------------------------------------------------------
Casts and type conversions
--------------------------------------------------------------------

Avoid C-style casts, const_cast, and reinterpret_cast
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.54 C-style casts **must not** be used.

      All type conversions **must** be done explicitly using the named C++ 
      casting operators; i.e., "static_cast", "const_cast", "dynamic_cast", 
      "reinterpret_cast".

8.55 The "const_cast" operator **should** be avoided. 

       Casting away "const-ness" is usually a poor programming decision and can 
       introduce errors.

.. note :: **Exception:** It may be necessary in some circumstances to cast 
           away const-ness, such as when calling const-incorrect APIs.

8.56 The "reinterpret_cast" **must not** be used unless absolutely necessary.

       This operator was designed to perform a low-level reinterpretation of 
       the bit pattern of an operand. This is needed only in special 
       circumstances and circumvents type safety.

Use the explicit qualifier to avoid unwanted conversions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.57  A class constructor that takes a single *non-default* argument, or a 
single argument with a *default* value, **must** be declared "explicit".

       This prevents compilers from performing unexpected (and, in many
       cases, unwanted!) implicit type conversions. For example::

          class MyClass
          {
          public:
             explicit MyClass(int i, double x = 0.0);
          };

       Note that, without the explicit declaration, an implicit conversion 
       from an integer to an object of type MyClass could be allowed. For 
       example::

          MyClass mc = 2;

       Clearly, this is confusing. The "explicit" keyword forces the 
       following usage pattern::

          MyClass mc(2);

       to get the same result, which is much more clear.


-----------------------------
Memory management
-----------------------------

Allocate and deallocate memory in the same scope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.58 Memory **should** be deallocated in the same scope in which it is 
allocated.

8.59 All memory allocated in a class constructor **should** be deallocated 
in the class destructor. 

      Note that the intent of constructors is to acquire resources and the 
      intent of destructors is to free those resources.

8.60 Pointers **should** be set to null explicitly when memory is deallocated.
This makes it easy to check pointers for "null-ness" when needed.


Use new/delete consistently
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8.61 Data managed exclusively within C++ code **must** be allocated and 
deallocated using the "new" and "delete" operators.

      The operator "new" is type-safe, simpler to use, and less error-prone
      than the "malloc" family of C functions.  C++ new/delete operators
      **must not** be combined with C malloc/free functions.

8.62 Every C++ array deallocation statement **must** include "[ ]" 
(i.e., "delete[ ]") to avoid memory leaks.

      The rule of thumb is: when "[ ]" appears in the allocation, then "[ ]"
      **must** appear in the corresponding deallocation statement.


