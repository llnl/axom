.. _inlet_advanced_types_label:

##############
Advanced Types
##############

In addition to Inlet's primitive types (bool, int, double, string), user-defined types
and functions can also be defined as part of an input file.

In this section, we first describe how `Individual Structs`_ can be added to our schemas.
We then extend it to `Arrays and Dictionaries of Structs`_.

.. _indiv structs label:

******************
Individual Structs
******************

Defining And Storing
--------------------

To add a single (i.e., not array) user-defined type to the input file, use the ``addStruct``
method of the Inlet or Container classes to add a Container (collection of Fields and sub-Containers)
that will represent the fields of the struct.

Consider a simple Lua table that contains only primitive types, whose definition might look like:

.. code-block:: Lua

  car = {
      make = "BestCompany",
      seats = 2,
      horsepower = 200,
      color = "blue"
  }

or in YAML, something like:

.. code-block:: YAML
  
  car:
    make: BestCompany
    seats: 2
    horsepower: 200
    color: blue

Its Inlet schema can be defined as follows:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_start
   :end-before: _inlet_userdef_simple_end
   :language: C++

This would be used by creating an ``inlet::Container`` for the ``Car`` instance and then defining the struct
schema on that subcontainer, e.g.:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_usage_start
   :end-before: _inlet_userdef_simple_usage_end
   :language: C++

.. note:: 
  The definition of a static ``defineSchema`` member function is not required, and is just used
  for convenience.  The schema definition for a class or struct could also be implemented as a
  free function for third-party types, or even in the same place as the sub-container declaration.

Accessing
---------

In order to convert from Inlet's internal representation to a custom C++ ``struct``, you must provide 
deserialization logic.  This is accomplished by a specialization of the ``FromInlet<T>`` functor for your
type ``T``, which implements a member function with the signature ``T operator()(const inlet::Container&)``.
This function should return an instance of type ``T`` with its members populated with the corresponding
fields in the input file. For the simple ``Car`` example whose schema is defined above, the specialization
might look like:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_frominlet_start
   :end-before: _inlet_userdef_simple_frominlet_end
   :language: C++

In the above snippet, ``Container::operator[]`` is used to extract data from Inlet's internal representation
which is automatically converted to the correct member variable types when the function's return value is
constructed.  This conversion does not happen automatically for user-defined types.
If a ``Car`` object as defined above is located at the path "car" within the input file, it can be retrieved as follows:

.. code-block:: C++

  Car car = inlet["car"].get<Car>();

**********************************
Arrays and Dictionaries of Structs
**********************************

Arrays of user-defined types are also supported in Inlet.  

Defining And Storing
--------------------

Consider a collection of cars, described in Lua as:

.. code-block:: Lua

  fleet = {
    {
      make = "Globex Corp",
      seats = 3,
      horsepower = 155,
      color = "green"
    },
    {
      make = "Initech",
      seats = 4,
      horsepower = 370
      -- uses default value for color
    },
    {
      make = "Cyberdyne",
      seats = 1,
      horsepower = 101,
      color = "silver"
    }
  }

or in YAML, as:

.. code-block:: YAML

  fleet:
    - make: Globex Corp
      seats: 3
      horsepower: 155
      color: green
    - make: Initech
      seats: 4
      horsepower: 370
      # uses default value for color
    - make: Cyberdyne
      seats: 1
      horsepower: 101
      color: silver

First, use the ``addStructArray`` function to create a subcontainer, then define the schema on that container
using the same ``Car::defineSchema`` used above:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_collection_usage_start
   :end-before: _inlet_userdef_collection_usage_end
   :language: C++

.. note::
  The schema definition logic for a struct is identical between individual instances of structs
  and arrays of structs.  The distinction is made by ``Container`` on which the struct schema is
  defined - specifically, whether it is obtained via ``addStruct`` or ``addStructArray``.

Associative arrays are also supported, using string keys or a mixture of string and integer keys.
The ``addStructDictionary`` function can be used analogously to the ``addStructArray`` function
for these associative arrays.

.. note::
  Although many of Inlet's input file languages do not distinguish between a "dictionary" type
  and a "record" type, Inlet treats them differently for type safety reasons:

  *Dictionaries* use arbitrary strings or integers for their keys, and their values (entries)
  can only be retreived as a homogeneous type.  In other words, dictionaries must map to
  ``std::unordered_map<Key, Value>`` for fixed key and value types.

  *Structs* contain a fixed set of named fields, but these fields can be of any type.
  As the name suggests, these map to ``structs`` in C++.

Accessing
---------

As with the schema definition, the ``FromInlet`` specialization for a user-defined type will work for both
single instances of the type and arrays of the type.

To retrieve an array of structs as a contiguous array of user-defined type, use ``std::vector``:

.. code-block:: C++

  auto fleet = base["fleet"].get<std::vector<Car>>();

Some input file languages support non-contiguous array indexing, so you can also retrieve arrays 
as ``std::unordered_map<int, T>``:

.. code-block:: C++

  auto fleet = inlet["fleet"].get<std::unordered_map<int, Car>>();

.. note::
  If a non-contiguous array is retrieved as a (contiguous) ``std::vector``, the elements will be
  ordered by increasing index.

String-keyed dictionaries are implemented as ``std::unordered_map<std::string, T>`` and can be retrieved
in the same way as the array above.  For dictionaries with a mix of string and integer keys, the
``inlet::VariantKey`` type can be used, namely, by retrieving a ``std::unordered_map<inlet::VariantKey, T>``.

Variant Struct Collections
--------------------------

Collections whose entries can be selected from a finite set of user-defined struct types can be
defined with ``addVariantStructArray`` or ``addVariantStructDictionary``.  Each entry must contain
a string discriminator field, and each alternative still uses its normal ``FromInlet`` specialization
for retrieval:

.. code-block:: C++

  using Shape = std::variant<Circle, Box>;

  auto shapes = inlet.addVariantStructArray<Shape>("shapes", "kind");
  shapes.addAlternative<Circle>("circle", [](inlet::Container& circle) {
    circle.addDouble("radius").required();
  });
  shapes.addAlternative<Box>("box", [](inlet::Container& box) {
    box.addDouble("width").required();
    box.addDouble("height").required();
  });

  auto values = inlet["shapes"].get<std::vector<Shape>>();

Verification fails if an entry omits the discriminator or uses a discriminator value that was not
registered as an alternative.
