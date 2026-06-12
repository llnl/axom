############
Simple Types
############

To help structure your input file, Inlet categorizes the information into
two types: Fields and Containers.

Fields refer to the individual scalar values that are either at the global level or
that are contained inside of a Container.

Containers can contain multiple Fields, other sub-Containers, as well as a single array or a single dictionary.

.. note::  There is a global Container that holds all top-level Fields.  This can be
  accessed via your `Inlet` class instance.


*******
Fields
*******

In Inlet, Fields represent an individual scalar value of primitive type.
There are four supported field types: ``bool``, ``int``, ``double``, and ``string``.

In this example we will be using the following part of an input file:

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_input_start
   :end-before: _inlet_simple_types_fields_input_end
   :language: lua


Defining And Storing
--------------------

This example shows how to add the four simple field types with descriptions to the
input file schema and add their values, if present in the input file, to the Sidre
DataStore to be accessed later.

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_add_start
   :end-before: _inlet_simple_types_fields_add_end
   :language: C++

You can also add default values to Fields to fall back to if they are not defined
in your input file. The last added Field was intentionally not present in the input file.  Not all
fields need to be present, unless they are marked required, like ``a_simple_double``.

Accessing
---------

Accessing field values stored in Inlet can be accessed via their name with the ``[]`` operator
or through the templated ``get<T>`` function.  The ``[]`` operator is more streamlined but
can lead to compile time ambiquity depending on how it is used.  The example below shows
an example of this.

Prior to accessing optional fields, you should verify they were provided by the user via
the ``contains`` function.

The ``contains`` function returns ``true`` if the field was *either*
provided by the user or via a default.  To check if the field was provided by the user (and not
via a default), you can use the ``isUserProvided`` method, which returns ``true`` if the value
was provided by the user in the input file.

Accessing a value that was not provided by the user, or 
a default value, will result in a runtime error.

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_access_start
   :end-before: _inlet_simple_types_fields_access_end
   :language: C++

.. note:: The field ``does_not_exist`` was purposefully left this out of the
   user-provided input file to show no warnings/errors are thrown during runtime
   for defining optional fields in the schema.


**********
Containers
**********

Containers help with grouping associated data together into a single collection. Containers can
contain multiple individually named Fields, multiple sub-Containers, as well as a single
array or a single dictionary.

In this example, we will be using the following part of an input file:

.. literalinclude:: ../../examples/containers.cpp
   :start-after: _inlet_simple_types_containers_input_start
   :end-before: _inlet_simple_types_containers_input_end
   :language: lua


Defining And Storing
--------------------

This example shows how to add a Container with a nested Container to the
input file schema and add the underlying field values to the Sidre
DataStore to be accessed later.

.. literalinclude:: ../../examples/containers.cpp
   :start-after: _inlet_simple_types_containers_add_start
   :end-before: _inlet_simple_types_containers_add_end
   :language: C++

This example also shows that the ``color`` Field that was not given in the
input file but used the default value that was specified in the schema.

.. note:: Inlet also has an ``addStruct`` member for defining more complex types,
   such as nested structures. See :ref:`Advanced Types <inlet_advanced_types_label>`
   for more details


Accessing
---------

Field values stored inside a container can be accessed via their name with the ``[]`` operator.
They can be accessed from the Inlet class instance with their fully qualified name or you
can get the Container instance first, then access it with the relative name.

.. literalinclude:: ../../examples/containers.cpp
   :start-after: _inlet_simple_types_containers_access_start
   :end-before: _inlet_simple_types_containers_access_end
   :language: C++


******
Arrays
******

Primitive arrays store a collection of values under integer keys.  The ``addBoolArray``,
``addIntArray``, ``addDoubleArray``, and ``addStringArray`` schema methods expect
all values in the array to have the requested primitive type.

In this example, both arrays contain only integer values.  The first input array
is contiguous and the second uses explicit integer keys:

.. literalinclude:: ../../examples/homogeneous_collections.cpp
   :start-after: _inlet_simple_types_homogeneous_arrays_input_start
   :end-before: _inlet_simple_types_homogeneous_arrays_input_end
   :language: lua

For arrays whose values can be any supported primitive type, use ``addVariantArray``.
Variant arrays store values as ``inlet::VariantValue``, which is a
``std::variant<bool, int, double, std::string>``.

In this example, the arrays contain mixed primitive value types:

.. literalinclude:: ../../examples/variant_collections.cpp
   :start-after: _inlet_simple_types_variant_arrays_input_start
   :end-before: _inlet_simple_types_variant_arrays_input_end
   :language: lua

Defining And Storing
--------------------

Homogeneous primitive arrays are added with the type-specific ``add*Array`` method:

.. literalinclude:: ../../examples/homogeneous_collections.cpp
   :start-after: _inlet_simple_types_homogeneous_collections_add_start
   :end-before: _inlet_simple_types_homogeneous_collections_add_end
   :language: C++

Variant arrays are added with ``addVariantArray``:

.. literalinclude:: ../../examples/variant_collections.cpp
   :start-after: _inlet_simple_types_variant_collections_add_start
   :end-before: _inlet_simple_types_variant_collections_add_end
   :language: C++

Accessing
---------

Contiguous homogeneous arrays can be retrieved as ``std::vector<T>``:

.. literalinclude:: ../../examples/homogeneous_collections.cpp
   :start-after: _inlet_simple_types_homogeneous_arrays_access_vector_start
   :end-before: _inlet_simple_types_homogeneous_arrays_access_vector_end
   :language: C++

Integer-keyed homogeneous arrays can be retrieved as
``std::unordered_map<int, T>`` when the original indices are needed:

.. literalinclude:: ../../examples/homogeneous_collections.cpp
   :start-after: _inlet_simple_types_homogeneous_arrays_access_map_start
   :end-before: _inlet_simple_types_homogeneous_arrays_access_map_end
   :language: C++

Contiguous variant arrays can be retrieved as ``std::vector<inlet::VariantValue>``.
Integer-keyed variant arrays can be retrieved as
``std::unordered_map<int, inlet::VariantValue>`` when the original indices are
needed.

.. literalinclude:: ../../examples/variant_collections.cpp
   :start-after: _inlet_simple_types_variant_arrays_access_start
   :end-before: _inlet_simple_types_variant_arrays_access_end
   :language: C++

************
Dictionaries
************

Dictionaries store a collection of values under arbitrary string keys or a mix of
string and integer keys.  The ``addBoolDictionary``, ``addIntDictionary``,
``addDoubleDictionary``, and ``addStringDictionary`` schema methods expect all
values in the dictionary to have the requested primitive type.

In this example, all dictionary values are integers:

.. literalinclude:: ../../examples/homogeneous_collections.cpp
   :start-after: _inlet_simple_types_homogeneous_dictionary_input_start
   :end-before: _inlet_simple_types_homogeneous_dictionary_input_end
   :language: lua

For dictionaries whose values can be any supported primitive type, use
``addVariantDictionary``.  Mixed string and integer keys are represented with
``inlet::VariantKey``.

In this example, the dictionary has mixed primitive values and mixed key types:

.. literalinclude:: ../../examples/variant_collections.cpp
   :start-after: _inlet_simple_types_variant_dictionary_input_start
   :end-before: _inlet_simple_types_variant_dictionary_input_end
   :language: lua

Defining And Storing
--------------------

Homogeneous primitive dictionaries are added with the type-specific
``add*Dictionary`` method:

.. literalinclude:: ../../examples/homogeneous_collections.cpp
   :start-after: _inlet_simple_types_homogeneous_collections_add_start
   :end-before: _inlet_simple_types_homogeneous_collections_add_end
   :language: C++

Variant dictionaries are added with ``addVariantDictionary``:

.. literalinclude:: ../../examples/variant_collections.cpp
   :start-after: _inlet_simple_types_variant_collections_add_start
   :end-before: _inlet_simple_types_variant_collections_add_end
   :language: C++

Accessing
---------

Mixed-key homogeneous dictionaries can be retrieved as
``std::unordered_map<inlet::VariantKey, T>``.

.. literalinclude:: ../../examples/homogeneous_collections.cpp
   :start-after: _inlet_simple_types_homogeneous_dictionary_access_start
   :end-before: _inlet_simple_types_homogeneous_dictionary_access_end
   :language: C++

Mixed-key variant dictionaries can be retrieved as
``std::unordered_map<inlet::VariantKey, inlet::VariantValue>``.

.. literalinclude:: ../../examples/variant_collections.cpp
   :start-after: _inlet_simple_types_variant_dictionary_access_start
   :end-before: _inlet_simple_types_variant_dictionary_access_end
   :language: C++
