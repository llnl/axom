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

Aliased input paths
-------------------

A function's schema name and its path in the input do not need to match.
Inlet provides an ``InputPath`` descriptor to allow a path to retain the same meaning
when a schema is expanded across a struct array or dictionary:

.. code-block:: C++

  using axom::inlet::InputPath;

  // Every element reads the same root-level callback.
  shapes.addFunction(
    "transform_callback",
    axom::inlet::FunctionTag::Vector,
    {axom::inlet::FunctionTag::Vector},
    InputPath::exact("shared_transform"));

  // Every element reads its own "transform" callback.
  shapes.addFunction(
    "transform_callback",
    axom::inlet::FunctionTag::Vector,
    {axom::inlet::FunctionTag::Vector},
    InputPath::relativeToCollectionElement("transform"));

An exact path is used unchanged, regardless of where the function is stored in the schema.
A collection-relative path is joined to each concrete collection element.
Outside a collection, it is relative to the current ``Container``.

The ``addFunction(name, returnType, argumentTypes, description, pathOverride)`` overload
remains available. Its string override retains the rule that
it is exact outside a struct collection and relative to each element inside one.
Prefer ``InputPath`` when adding new aliases so this behavior is explicit at the call site.

Functions as value alternatives
-------------------------------

Some schemas accept either a concrete value or a function that computes that value.
Use ``addFunctionAsValueAlternative`` to declare this relationship explicitly.
The recommended overload lets Inlet own the callback's internal storage name and
associates it with the concrete field's input path:

.. code-block:: C++

  inlet.addFunctionAsValueAlternative(
    axom::inlet::FunctionTag::Vector,
    {},
    "scale");
  inlet.addDoubleArray("scale");

With this schema, ``scale = {2.0, 3.0, 4.0}`` populates ``scale``,
while ``scale = function() return {2.0, 3.0, 4.0} end`` supplies the function
alternative associated with ``scale``. Use ``containsFunctionValueAlternative("scale")``
and ``getFunctionValueAlternative("scale")`` on the containing ``Container`` to
query and retrieve it without depending on an internal schema name.

An overload that takes a schema name first remains available when the callback needs
an independently addressable schema entry. Both forms accept an ``InputPath`` descriptor
when exact or collection-relative resolution needs to be explicit.

The two schema entries may be added in either order. A function encountered at a normal
field path remains a type error unless this alternative has been declared.

Only the selected representation exists: ``contains`` reports the concrete field when a
value was supplied, and ``containsFunctionValueAlternative`` reports the callback when a
function was supplied. A value with an unrelated type matches neither representation and
fails verification. The shared input path is recognized by strict containers and is not
reported as unexpected.

The returned function and the concrete field remain independently verifiable schema entries.
Consequently, ``required()`` and registered verifiers apply to the entry on
which they are configured. The narrow value-alternative API does not currently provide
a group-level annotation meaning "either representation is required."

Accessing
---------

To retrieve a function, both the implicit conversion and ``get<T>`` syntax is supported.  For example,
a function can be retrieved as follows:

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
