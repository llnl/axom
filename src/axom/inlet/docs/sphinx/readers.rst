#######
Readers
#######

Inlet has built-in support for three input file languages: JSON, Lua, and YAML.
Due to language features, not all readers support all Inlet features.
Below is a table that lists supported features:


.. list-table:: Supported Language Features
   :header-rows: 1

   * - 
     - JSON
     - Lua
     - YAML
   * - Primitive Types
     - bool, double, int, string
     - bool, double, int, string
     - bool, double, int, string
   * - Dictionaries
     - X
     - X
     - X
   * - Arrays
     - X
     - X
     - X
   * - Non-contiguous Arrays
     - 
     - X
     - 
   * - Mixed-typed key Arrays
     - 
     - X
     - 
   * - Callback Functions
     - 
     - X
     - 

***********************
Extra Lua Functionality
***********************

The `LuaReader` class has the ability to access the entire Lua State via the protected member function
``LuaReader::solState()``.  This allows you to fully utilize the Sol library, documented in
`Sol's documentation <https://sol2.readthedocs.io/en/v2.20.6/index.html>`_. This is an advanced feature
and not recommended unless there is a good reason.  We provide an example on how to create a derived
reader class here:

.. literalinclude:: ../../examples/lua_library.cpp
   :start-after: _inlet_sol_state_start
   :end-before: _inlet_sol_state_end
   :language: C++

Inlet opens four Lua libraries by default: ``base``, ``math``, ``string``, ``package``.
All libraries are documented in `Sol's open_library documentation <https://sol2.readthedocs.io/en/v2.20.6/api/state.html?highlight=open_libraries#enumerations>`_.

.. warning::

   Lua input is executable code, and Inlet does not sandbox it. The ``package`` library can
   load additional Lua or native modules, and Inlet imposes no CPU, memory, recursion, or
   execution-time limits. Only parse Lua input from trusted sources. Exposing more libraries
   or modifying ``solState()`` can grant the input more capabilities and change values after
   Inlet has read or verified them.

For example, you can add the ``io`` library by doing this:

.. literalinclude:: ../../examples/lua_library.cpp
   :start-after: _inlet_io_library_add_start
   :end-before: _inlet_io_library_add_end
   :language: C++
