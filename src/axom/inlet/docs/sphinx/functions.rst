##################
Function Callbacks
##################

For input file types that support functions, e.g., Lua, functions can also be read from the input file
into a ``std::function``, the wrapper for callables provided by the C++ standard library.

Defining And Storing
--------------------

This is accomplished by calling ``addFunction`` on an Inlet or Container object.

Consider the following Lua function that accepts a vector in **R**\ :sup:`2` or **R**\ :sup:`3` and returns a double:

.. code-block:: Lua

  coef = function (v)
    if v.dim == 2 then
      return v.x + (v.y * 0.5)
    else
      return v.x + (v.y * 0.5) + (v.z * 0.25)
    end
  end

The schema for this function would be defined as follows:

.. literalinclude:: ../../examples/mfem_coefficient.cpp
   :start-after: _inlet_mfem_func_coef_start
   :end-before: _inlet_mfem_func_coef_end
   :language: C++

The return type and argument types are described with the ``inlet::FunctionTag`` enumeration, which has the following members:

  * ``Double`` - corresponds to a C++ ``double``
  * ``String`` - corresponds to a C++ ``std::string``
  * ``Vector`` - corresponds to a C++ ``inlet::InletVector``
  * ``Void`` - corresponds to C++ ``void``, should only be used for functions that don't return a value

Note that a single type tag is passed for the return type, while a vector of tags is passed
for the argument types.  Currently a maximum of two arguments are supported.
To declare a function with no arguments, simply leave the list of argument types empty.

.. note::  The ``InletVector`` type (and its Lua representation) are statically-sized vectors with
  a maximum dimension of three.  That is, they can also be used to represent two-dimensional vectors.

A Lua callback declared with a ``Vector`` return type may return either a ``Vector.new(...)``
value or an ordinary Lua table containing one to three numeric components. Ordinary table
returns must use contiguous integer indices starting at one; sparse tables, named entries,
and non-numeric components are rejected when the callback is called.

In Lua, the following operations on the ``Vector`` type are supported (for ``Vector`` s ``u``, ``v``, and ``w``):

1. Construction of a 3D vector: ``u = Vector.new(1, 2, 3)``
#. Construction of a 2D vector: ``u = Vector.new(1, 2)``
#. Construction of an empty vector (default dimension is 3): ``u = Vector.new()``
#. Vector addition and subtraction: ``w = u + v``, ``w = u - v``
#. Vector negation: ``v = -u``
#. Scalar multiplication: ``v = u * 0.5``, ``v = 0.5 * u``
#. Indexing (1-indexed for consistency with Lua): ``d = u[1]``, ``u[1] = 0.5``
#. L2 norm and its square: ``d = u:norm()``, ``d = u:squared_norm()``
#. Normalization: ``v = u:unitVector()``
#. Dot and cross products: ``d = u:dot(v)``, ``w = u:cross(v)``
#. Dimension retrieval: ``d = u.dim``
#. Component retrieval: ``d = u.x``, ``d = u.y``, ``d = u.z``

Functions as value alternatives
-------------------------------

Some schemas allow an input to be either a concrete value or a function that computes it.
Declare the concrete entry normally, then associate a function with the same input name:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_function_value_alternative_schema_start
   :end-before: _inlet_function_value_alternative_schema_end
   :language: C++
   :dedent: 2

Declare the alternative before the concrete entry it applies to. Both use the same relative
and slash-delimited paths as other ``Container`` methods, and either may be declared through
a parent or a child ``Container``. The example accepts either of these Lua inputs:

.. code-block:: Lua

  scale = 2.0

  -- or
  scale = function() return 3.0 end

After verification, query which representation was supplied before retrieving it:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_function_value_alternative_access_start
   :end-before: _inlet_function_value_alternative_access_end
   :language: C++
   :dedent: 2

For the shared input name, ``contains`` reports only the concrete representation and
``containsFunctionValueAlternative`` reports only the function representation.
Either form counts as user-provided input, and both are recognized by strict Containers.
An unrelated input type fails verification.

An alternative is an ordinary ``Function``, so it is returned by ``addFunctionAsValueAlternative``
as a ``Verifiable<Function>`` and can carry the usual schema constraints, such as ``required()``.
Validation attached to the *concrete* schema entry does not apply to the function result,
so applications that constrain both forms should validate after resolving the representation.

Generated Sphinx and JSON Schema documentation show only the concrete entry.
Application documentation should describe the function form when it is part of the Lua interface.

Accessing
---------

To retrieve a function, both the implicit conversion and ``get<T>`` syntax is supported.
For example, a function can be retrieved as follows:

.. literalinclude:: ../../examples/mfem_coefficient.cpp
   :start-after: _inlet_mfem_coef_simple_retrieve_start
   :end-before: _inlet_mfem_coef_simple_retrieve_end
   :language: C++

It can also be assigned directly to a ``std::function`` without the need to use ``get<T>``:

.. code-block:: C++

  std::function<double(FunctionType::Vector)> coef = inlet["coef"];

Additionally, if a function does not need to be stored, the overhead of a copy can be eliminated
by calling it directly:

.. code-block:: C++

  double result = inlet["coef"].call<double>(axom::inlet::FunctionType::Vector{3, 5, 7});

.. note::  Using ``call<ReturnType>(ArgType1, ArgType2, ...)`` requires both that the return type
  be explicitly specified and that argument types be passed with the exact type as used in the
  signature defined as part of the schema.  This is because the arguments do not participate in
  overload resolution.

Callbacks retrieved from Inlet keep their Lua state alive, so they remain callable after the
Inlet and Reader are destroyed. Callbacks from one ``LuaReader`` share mutable interpreter
state and must not be invoked concurrently without synchronization. Lua execution errors and
invalid callback return values throw ``std::runtime_error`` at the call site.
