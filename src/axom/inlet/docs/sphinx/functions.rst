##################
Function Callbacks
##################

Inlet can read functions from input formats such as Lua into a C++ ``std::function``.

Defining And Storing
--------------------

Call ``addFunction`` on an Inlet or Container object to declare an input function.

Consider the following Lua function that accepts a vector in **R**\ :sup:`2` or **R**\ :sup:`3` and returns a double:

.. code-block:: Lua

  coef = function (v)
    if v.dim == 2 then
      return v.x + (v.y * 0.5)
    else
      return v.x + (v.y * 0.5) + (v.z * 0.25)
    end
  end

Declare its signature in the schema:

.. literalinclude:: ../../examples/mfem_coefficient.cpp
   :start-after: _inlet_mfem_func_coef_start
   :end-before: _inlet_mfem_func_coef_end
   :language: C++

Use ``inlet::FunctionTag`` to specify the return and argument types:

  * ``Double`` - corresponds to a C++ ``double``
  * ``String`` - corresponds to a C++ ``std::string``
  * ``Vector`` - corresponds to a C++ ``inlet::InletVector``
  * ``Void`` - corresponds to C++ ``void``, should only be used for functions that don't return a value

Pass one tag for the return type and a vector of tags for up to two arguments.
Leave the argument list empty for a function with no arguments.

.. note::  ``InletVector`` stores up to three components, and its ``dim`` member records
   the vector's dimension. In Lua, ``Vector.new`` creates 2D or 3D vectors.

A Lua callback with a ``Vector`` return type may return ``Vector.new(...)`` or a Lua
table with one to three numeric components. Table keys must be contiguous integers
starting at one. When converting a callback's return value to a vector,
Inlet reports a SLIC error if the table has gaps, named entries, or non-numeric components.

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

Values or Functions
-------------------

Some entries accept a value or a function at the same input path.
Call ``orFunction`` on the entry to allow both forms:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_or_function_schema_start
   :end-before: _inlet_or_function_schema_end
   :language: C++
   :dedent: 2

For ``addDouble``, the function returns ``Double``, for ``addString``, it returns ``String``,
and for ``addDoubleArray``, it returns ``Vector``.
Declare function arguments with ``FunctionTag`` values, as for ``addFunction``.
Other entry types do not support ``orFunction``.

The example accepts either of the following Lua inputs:

.. code-block:: Lua

  -- a value
  scale = 2.0

  -- or a function
  scale = function(t) return 2.0 * t end

Use ``evaluate`` for either form. When a function was supplied in the input,
it evaluates it with the given arguments, and when a value was supplied, it ignores the arguments:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_or_function_access_start
   :end-before: _inlet_or_function_access_end
   :language: C++
   :dedent: 2

``evaluate`` passes numeric arguments as ``double`` and character strings as ``std::string``.
For an array of doubles, use ``evaluate<std::vector<double>>``.
Its function may return one to three components.

To tell which form the input supplied, check ``type()``, which reports ``InletType::Function``
for a function. Reading a supplied function with ``get`` is an error.

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_or_function_type_start
   :end-before: _inlet_or_function_type_end
   :language: C++
   :dedent: 2

Both forms make the entry present and user-provided. They also satisfy ``required`` and
``strict`` Containers. A Field's default applies only when the input supplies neither form.
Field constraints such as ``range`` apply only to a value because Inlet does not call
the function during ``verify()``. Readers that cannot supply functions,
such as the YAML and JSON readers, read only the value.
In a collection of structs, each element may use either form.

Accessing
---------

Retrieve a function with ``get<T>``:

.. literalinclude:: ../../examples/mfem_coefficient.cpp
   :start-after: _inlet_mfem_coef_simple_retrieve_start
   :end-before: _inlet_mfem_coef_simple_retrieve_end
   :language: C++

Or assign it directly to a ``std::function``:

.. code-block:: C++

  std::function<double(FunctionType::Vector)> coef = inlet["coef"];

Use ``call`` to invoke it without copying it:

.. code-block:: C++

  double result = inlet["coef"].call<double>(axom::inlet::FunctionType::Vector{3, 5, 7});

.. note::  ``call<ReturnType>(ArgType1, ArgType2, ...)`` requires an explicit return type.
   Argument types must match the schema signature exactly because they do not
   participate in overload resolution.

Lua callbacks copied into a ``std::function`` keep their Lua state alive and remain callable
after the Inlet and Reader are destroyed. A reference to an Inlet-owned ``Function`` does
not extend its lifetime. Callbacks from one ``LuaReader`` share mutable interpreter state
and must not be invoked concurrently without synchronization.

Inlet reports errors according to when they occur:

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - Problem
     - Reported through
   * - Invalid keys, undefined entries, or ambiguous schema names
     - SLIC diagnostics
   * - Missing required input, wrong types, or failed verifiers
     - ``verify()`` and ``VerificationError``
   * - Lua execution errors or callback results that cannot convert to the declared type
     - SLIC errors, which include Lua's error message for execution errors

Inlet does not evaluate callbacks during ``verify()`` unless a custom verifier calls them.
Successful verification does not guarantee that a later callback invocation will succeed.
A callback that fails inside a custom verifier reports its SLIC error during ``verify()``.
