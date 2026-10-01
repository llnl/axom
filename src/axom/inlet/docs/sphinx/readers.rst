.. _inlet_readers_label:

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
Extra Lua functionality
***********************

Derive from ``LuaReader`` to access its Lua state through the protected
``solState()`` method. Use `Sol <https://sol2.readthedocs.io/en/v2.20.6/index.html>`_
to add libraries or C++ bindings, as this example shows:

.. literalinclude:: ../../examples/lua_library.cpp
   :start-after: _inlet_sol_state_start
   :end-before: _inlet_sol_state_end
   :language: C++

Inlet opens four Lua libraries by default: ``base``, ``math``, ``string``, ``package``.
All libraries are documented in `Sol's open_library documentation <https://sol2.readthedocs.io/en/v2.20.6/api/state.html?highlight=open_libraries#enumerations>`_.

.. warning::

   Only parse Lua input from trusted sources. Inlet does not sandbox Lua or limit its
   resource use or execution time. The ``package`` library can load Lua or native modules.
   Adding libraries or bindings can grant input code more capabilities. Changing the Lua
   state can also affect callbacks after Inlet has verified the input.

For example, you can add the ``io`` library by doing this:

.. literalinclude:: ../../examples/lua_library.cpp
   :start-after: _inlet_io_library_add_start
   :end-before: _inlet_io_library_add_end
   :language: C++
