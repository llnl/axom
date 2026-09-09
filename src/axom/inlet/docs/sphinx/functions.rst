##################
Function callbacks
##################

For input file types that support functions, e.g., Lua, functions can also be read from the input file
into a ``std::function``, the wrapper for callables provided by the C++ standard library.

Defining and storing
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

Some schemas accept either a concrete value or a function that computes it.
Declare the function alternative before the concrete entry, using the same input name:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_function_value_alternative_schema_start
   :end-before: _inlet_function_value_alternative_schema_end
   :language: C++
   :dedent: 2

Declaring the alternative afterwards is an error because Inlet has already read the
concrete entry. Both declarations use paths relative to their ``Container``,
with slashes separating nested names. They may use different parent or child Containers
as long as they refer to the same input path. The example accepts either Lua input:

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

For the shared input name, ``contains`` reports the concrete value, including a default,
and ``containsFunctionValueAlternative`` reports a supplied function. Check the function
first because a concrete default can exist alongside it. Either supplied form counts as
user-provided input, and strict Containers recognize both. An unrelated input type fails
verification.

``addFunctionAsValueAlternative`` returns a ``Verifiable<Function>``.
Calling ``required()`` on it requires a function, and a concrete value does not satisfy the reqquirement.
Similarly, a function alone does not satisfy ``required()`` on the concrete entry.
To require either form, leave both entries optional and register a verifier on the
root Container that checks whether either representation exists.

Inlet does not evaluate the function automatically. Constraints on the concrete entry,
such as a range, do not apply to the function result, so applications should validate the
value after resolving and evaluating the supplied representation.

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

Lua callbacks copied into a ``std::function`` keep their Lua state alive and remain callable
after the Inlet and Reader are destroyed. A reference to an Inlet-owned ``Function`` does
not extend its lifetime. Callbacks from one ``LuaReader`` share mutable interpreter state
and must not be invoked concurrently without synchronization.

Lua execution errors and invalid callback return values throw ``axom::inlet::InletError``
when the callback is invoked. Inlet uses different reporting mechanisms for schema and
input validation:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Kind of problem
     - Reported through
     - Examples
   * - API or schema misuse
     - SLIC diagnostics
     - an empty or malformed key, a lookup for an entry that was never defined,
       a name that is ambiguous between a container, field, and function
   * - Contents of the input file
     - ``verify()`` and ``VerificationError``
     - a required entry is missing, a value has the wrong type or fails a
       registered verifier
   * - Failure while calling an input function
     - ``InletError`` (derived from ``std::runtime_error``)
     - the Lua function raises an error, or returns something that cannot be
       converted to the declared return type

Inlet does not evaluate callbacks during ``verify()`` unless a custom verifier calls them.
Successful verification does not guarantee that a later callback invocation will succeed.
If a callback throws ``InletError`` inside a custom verifier, the exception propagates
out of ``verify()`` unless the verifier catches it.

Applications can catch ``InletError`` to report where the callback failed. For example,
Klee wraps it in a ``KleeError`` that identifies the owning shape or named operator,
the operator location, and the callback field.
