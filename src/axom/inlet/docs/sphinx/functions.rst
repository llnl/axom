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

Some schemas accept either a concrete value or a function that computes that value.
Use ``addFunctionAsValueAlternative`` to declare this relationship explicitly.
Inlet owns the callback's internal storage name and associates it with the concrete
field's public name:

.. code-block:: C++

  inlet.addFunctionAsValueAlternative(
    "scale",
    axom::inlet::FunctionTag::Vector,
    {});
  inlet.addDoubleArray("scale");

With this schema, ``scale = {2.0, 3.0, 4.0}`` populates ``scale``,
while ``scale = function() return {2.0, 3.0, 4.0} end`` supplies the function
alternative associated with ``scale``. Use ``containsFunctionValueAlternative("scale")``
and ``getFunctionValueAlternative("scale")`` on the containing ``Container`` to
query and retrieve it without depending on an internal schema name.

Inside a struct collection, Inlet resolves the public value name against each concrete
element. For example, this Lua input supplies ``bar`` directly for one ``foo`` element
and computes it with a callback for another:

.. code-block:: Lua

  foo = {
    [7] = { bar = 2 },
    [12] = { bar = function() return 3 end }
  }

The corresponding schema is compiled as part of Inlet's function callback example:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_nested_callback_alternative_start
   :end-before: _inlet_nested_callback_alternative_end
   :language: C++
   :dedent: 2

After verification, each concrete element contains only the representation that was provided.
A ``FromInlet`` specialization can normalize both representations to a scalar:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_function_value_alternative_access_start
   :end-before: _inlet_function_value_alternative_access_end
   :language: C++

This example evaluates a supplied callback during conversion. An application that needs
deferred evaluation can instead call ``get<double()>()`` on the result
of ``getFunctionValueAlternative("bar")`` and store the returned ``std::function<double()>``.

The two schema entries may be added in either order. A function encountered at a normal
field path remains a type error unless this alternative has been declared.
The narrow API permits one callback alternative for a public value name in a given container
and declaring a second callback alternative for that name is an error.
Since one Lua object cannot simultaneously be a concrete value and a function,
at most one of the two supported representations can match.

Only the selected representation exists: ``contains`` reports the concrete field when a
value was supplied, and ``containsFunctionValueAlternative`` reports the callback when a
function was supplied. A value with an unrelated type matches neither representation and
fails verification. The shared public name is recognized by strict containers and is not
reported as unexpected.

The returned function and the concrete field remain independently verifiable schema entries.
Consequently, ``required()`` and registered verifiers apply to the entry on
which they are configured. The narrow value-alternative API does not currently provide
a group-level annotation meaning "either representation is required."

In Sidre, the callback entry retains the same signature metadata as an ordinary Inlet function
and is tagged as an internal value alternative. The live Lua callable remains in memory.
Inlet restart reconstructs containers and fields, not functions, so the tag prevents
the internal callback group from being mistaken for a field.
A persisted concrete value is reconstructed normally, while a Lua callback is not.

Generated Sphinx and JSON Schema documentation describe the concrete value form.
The internal callback entry is omitted: its storage name is not an input-file path,
and JSON inputs cannot provide a Lua function. Document the callback form separately when exposing
it as part of an application's Lua interface.

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

Callbacks copied out of Inlet keep their Lua state alive and remain callable after the Inlet
object is destroyed. Lua execution errors and invalid callback return values are reported as
``std::runtime_error`` at the call site.

The lifetime behavior is also demonstrated by the compiled function callback example:

.. literalinclude:: ../../examples/functions.cpp
   :start-after: _inlet_callback_after_destruction_start
   :end-before: _inlet_callback_after_destruction_end
   :language: C++
   :dedent: 2

All callbacks obtained from one ``LuaReader`` retain and share that reader's Lua state.
Copying a returned ``std::function`` does not clone the state, and mutations to Lua globals
or captured Lua tables made by one callback are therefore visible to the others.
Callbacks from the same reader must not be invoked concurrently:
even apparently read-only calls use the same Lua interpreter state.
Applications that need concurrent callback evaluation must serialize access to each state
or create independent ``LuaReader`` instances and parse the input separately for each thread.
