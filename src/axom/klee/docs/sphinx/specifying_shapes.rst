Specifying Shapes
=================

"Shaping", or "painting", is the process of adding non-conformal material
regions to a mesh. Traditionally, this has been done in code-specific formats
by each code that provides such a capability. Axom's Klee library provides
a way to read shape specifications in YAML or Lua files and apply the specified
geometry to a mesh.

Basics
------
Shapes in Klee are specified as a list. Each shape specifies its :code:`name`
and :code:`material`, as well as a description of its geometry.

In addition to the shapes themselves, a file must specify the number of
dimensions of the shapes (only 2 and 3 are allowed). This will be important
later when specifying transformations and slices.

.. code-block:: yaml

    dimensions: 3

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: stl
          path: wheel.stl
          units: cm
      - name: windshield
        material: glass
        geometry:
          format: stl
          path: windshield.stl
          units: in


The above example describes a series of 3D shapes. The first is a wheel
made of steel. Its geometry is specified in an STL file named :code:`wheel.stl`.
Since STL files don't have units embedded in them, we must specify them.
The second shape is named "windshield", is made of "glass", and its geometry
is specified in :code:`windshield.stl`. Note that Klee does not specify
what a particular material means. A material is simply a label which can
be used by a host code to apply properties to part of a mesh.

Lua Input Files
***************
Klee can also read Lua input files when Axom is configured with :code:`AXOM_ENABLE_LUA=ON`
and Sol support. Lua files use the same Klee schema as YAML files.
Klee evaluates the Lua once, when the input file is read, and then passes the
resulting values to Inlet. It does not re-evaluate the Lua during later loops,
timesteps, or shape-processing operations.

File-based reads use an explicit :code:`axom::klee::InputFormat` when provided.
Otherwise they infer the format from the file extension 
(:code:`.yaml`, :code:`.yml`, or :code:`.lua`) when present, and default to YAML.
Similarly, stream-based reads are YAML by default unless the caller passes :code:`axom::klee::InputFormat::Lua`.

The following YAML and Lua inputs are equivalent:

.. code-block:: yaml

    dimensions: 3

    shapes:
      - name: windshield
        material: glass
        geometry:
          format: stl
          path: windshield.stl
          units: cm
          operators:
            - rotate: 90
              axis: [0, 1, 0]
              center: [0, 0, -10]
            - translate: [10, 20, 30]

.. code-block:: lua

    dimensions = 3

    shapes = {
      {
        name = "windshield",
        material = "glass",
        geometry = {
          format = "stl",
          path = "windshield.stl",
          units = "cm",
          operators = {
            { rotate = 90, axis = {0, 1, 0}, center = {0, 0, -10} },
            { translate = {10, 20, 30} }
          }
        }
      }
    }

Because Lua is evaluated before Klee reads the input file,
ordinary table values can be generated programmatically:

.. code-block:: lua

    local dim = 2
    local r = 4.0
    local z = 8.0
    local x = 1.0
    local y = 2.0

    dimensions = dim

    shapes = {
      {
        name = "part",
        material = "steel",
        geometry = {
          format = "stl",
          path = "part.stl",
          units = "cm",
          operators = {
            { translate = (dim == 2) and {r, z} or {x, y, z} }
          }
        }
      }
    }

Caller-provided primitive values can be installed as initial Lua globals before a deck is
evaluated. This is useful when an application wants one deck to select between 2D and 3D
geometry, dimensions, or operator values at run time:

.. code-block:: c++

    axom::klee::LuaInputOptions options;
    options.initialGlobals = {
      {"dimensions", axom::klee::LuaGlobalValue {2}},
      {"shape_suffix", axom::klee::LuaGlobalValue {std::string {"2d"}}}
    };
    auto shapeSet = axom::klee::readShapeSet("shape.lua", options);

.. code-block:: lua

    local function shape_path()
      return "part_" .. shape_suffix .. ".stl"
    end

    shapes = {
      {
        name = "part",
        material = "steel",
        geometry = {
          format = "stl",
          path = shape_path(),
          units = "cm",
          operators = {
            { translate = (dimensions == 2) and {1.0, 2.0} or {1.0, 2.0, 3.0} }
          }
        }
      }
    }

Initial globals are Lua-only and may be booleans, integers, doubles, or strings.
Their names must be non-keyword ASCII Lua identifiers. They are ordinary mutable
globals—not a read-only context—and are allowed by Klee's unexpected-global check.
Deck code can reassign or delete them, so applications should treat them as initial
values rather than controls. Other helper values in the deck should still be
declared :code:`local`. Initial globals may not replace standard Lua globals such
as :code:`math` or :code:`package`.

Applications that need richer runtime customization can also provide a Lua
initialization chunk. Klee evaluates the chunk after installing
:code:`initialGlobals` and before parsing the deck. The chunk must return a table;
those table entries are then installed as initial globals while unrelated
unexpected globals in the deck remain errors. This allows host code to provide
helper functions, tables, and local closures without recompiling the application:

.. code-block:: c++

    axom::klee::LuaInputOptions options;
    options.initialization = axom::klee::LuaInitializationChunk {
      R"(
        local dim = 2
        local lift = 3.0

        local function offset(y)
          return function()
            return {0.0, y}
          end
        end

        return {
          dimensions = dim,
          lift = lift,
          offset = offset
        }
      )",
      "runtime_initialization"
    };
    auto shapeSet = axom::klee::readShapeSet("shape.lua", options);

.. code-block:: lua

    shapes = {
      {
        name = "part",
        material = "steel",
        geometry = {
          format = "stl",
          path = "part.stl",
          units = "cm",
          operators = {
            { translate = offset(lift) }
          }
        }
      }
    }

Initialization chunks must return a table whose exported keys are non-keyword
ASCII Lua identifiers. Exported values may be booleans, numbers, strings, tables,
or functions and retain their original Lua representation. An exported Lua integer,
for example, is not converted through a C++ floating-point value. Export names may
not collide with standard Lua globals or :code:`initialGlobals`.

Klee evaluates the chunk in an isolated Lua environment. Preloaded Lua
libraries and caller-provided initial globals remain visible, but global names
assigned by the chunk do not leak into the input deck unless they are returned
in the export table. Exported functions retain access to the chunk's private
environment. Exported globals are mutable: deck code can replace or delete them
and can mutate exported tables.

The environment isolation is shallow. Inherited objects such as :code:`math` and
:code:`package` are shared, so mutating a member of an inherited table is visible
to the deck. Initialization chunks and decks are trusted code: this mechanism is
not a security sandbox, :code:`package` may load additional code, and Klee imposes
no CPU, memory, or recursion limits.

Use :code:`local` helper functions and constants for intermediate values so the
global namespace contains only the Klee schema fields that Inlet should read.
For Lua input, a one-value scale is written as a one-entry table, for example
:code:`{ scale = {2.0} }`.

Selected operator fields may also be written as zero-argument Lua callbacks.
Klee evaluates each callback exactly once while reading the deck; the resulting
shape still contains ordinary affine or slice operators, not runtime Lua
functions. Callbacks should be pure functions of local deck variables.
Callbacks in a named operator are evaluated when that named operator is
constructed. Each :code:`ref` reuses the resulting concrete operator rather
than evaluating its callbacks again for the referring shape.

Callback evaluation order is deterministic. Klee processes
:code:`named_operators` before :code:`shapes`, entries in each list in source
order, and each geometry's operators in source order. Within a multi-field
operator, fields are evaluated in this order:

* :code:`rotate`, then :code:`center`, then :code:`axis` (when present)
* :code:`scale`, then :code:`center` (when present)
* perpendicular slice :code:`x`, :code:`y`, or :code:`z`, then
  :code:`origin`, :code:`normal`, and :code:`up` (when present)
* arbitrary slice :code:`origin`, then :code:`normal`, then :code:`up`

Klee does not coordinate Lua evaluation across MPI ranks. If an application
calls :code:`readShapeSet` on every rank, each rank reads and evaluates the deck
and initialization chunk independently. 

.. code-block:: lua

    local dim = 2
    local r = 4.0
    local z = 8.0
    local x = 1.0
    local y = 2.0

    dimensions = dim

    shapes = {
      {
        name = "part",
        material = "steel",
        geometry = {
          format = "stl",
          path = "part.stl",
          units = "cm",
          operators = {
            {
              translate = function()
                if dim == 2 then
                  return {r, z}
                end
                return {x, y, z}
              end
            }
          }
        }
      }
    }

Vector-valued callbacks return raw numeric Lua tables such as :code:`{x, y}` or :code:`{x, y, z}`.
The typed :code:`Vector.new(...)` object is also accepted.
Scalar-valued callbacks return a number.
String-valued callbacks return a string.
Supported callback fields are
:code:`translate`, :code:`axis`, :code:`center`, :code:`scale`, :code:`slice.origin`,
:code:`slice.normal`, :code:`slice.up`, :code:`rotate`, :code:`slice.x`,
:code:`slice.y`, :code:`slice.z`, :code:`convert_units_to`, and :code:`ref`.
For :code:`scale`, a one-entry table means uniform scaling and a multi-entry table means per-axis scaling.

Common Lua input errors are reported as Klee parsing errors.
A Lua input file read without Lua support reports:

.. code-block:: text

    Lua input files require Axom configured with AXOM_ENABLE_LUA=ON and Sol library support. Rebuild Axom with Lua enabled or convert the file to YAML.

Unsupported file extensions are rejected before parsing, and unexpected top-level globals
or syntax errors are reported during Inlet verification or Lua evaluation.
Unknown nested fields follow the same Inlet schema strictness rules as YAML input.

Error Reporting
***************
Klee validates user-provided YAML and Lua input files while reading them.
When validation fails, Klee reports the problem by throwing an exception,
usually :code:`axom::klee::KleeError`.

Using exceptions lets Klee return detailed feedback at the point where the invalid input is detected, 
including the input path reported by Inlet.

Callers that read input files should catch :code:`axom::klee::KleeError` and display :code:`what()`
or inspect :code:`getErrors()` when multiple verification errors are available.
Klee may still throw standard exceptions such as :code:`std::logic_error` or :code:`std::invalid_argument`
for programming errors or inconsistent manually constructed objects.

Callback failures include the field, owning shape or named operator,
and operator location, for example:

.. code-block:: text

    Error evaluating callback for 'translate' in shape 'part' operator 1: [Inlet] Lua function call failed: ...

Paths
*****
The paths specified in shapes are specified either as absolute paths
(begin with a :code:`/`), or as relative paths. Absolute paths are evaluated
as absolute paths in the file systems. Relative paths are evaluated relative
to the Klee input file (not the current working directory). For example, in the
file :code:`/path/to/my_shapes.yaml`, the table below illustrates how
different paths would be specified.

+---------------------------------+-----------------------------+
| Path in /path/to/my_shapes.yaml | Path in file system         |
+=================================+=============================+
| /an/absolute/path.stl           | /an/absolute/path.stl       |
+---------------------------------+-----------------------------+
| just_a_file.stl                 | /path/to/just_a_file.stl    |
+---------------------------------+-----------------------------+
| dir/and/a/file.stl              | /path/to/dir/and/a/file.stl |
+---------------------------------+-----------------------------+

Changing Dimensions on a Per-Shape Basis
****************************************
Some problem setups can require operating on shapes of different dimensions.

In addition to the global :code:`dimensions` field, Klee provides a per-shape override 
mechanism for the dimension of a shape via the :code:`geometry/dimensions` field

.. code-block:: yaml

    dimensions: 2

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: stl
          path: wheel.stl
          dimensions: 3
          units: cm

In the above snippet, the overall dimension of the problem is `2`, but the `wheel`
shape is 3-dimensional.
In such cases, the user code is responsible for conversions between
the shape dimension (e.g. the dimension in the specified file format) and the problem dimension
(i.e. the dimension of the computational mesh).

Alternatively, some :code:`operators` (described below), can handle the conversion 
from the shape's dimension and the problem dimension.
For example, the input to the :code:`slice`` operator is a three-dimensional shape
and the output is a two-dimensional shape. This is specified via the 
:code:`start_dimensions` field on the :code:`geometry` of a :code:`shape`.
After the `operators` have been applied, the (end) dimension of the shape
will match that of the (global or per-shape) `dimensions`.

.. code-block:: yaml

    dimensions: 2

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: stl
          path: wheel.stl
          start_dimensions: 3
          units: cm
          operators:
            - slice:
                x: 10


Overlay Rules
-------------
Shapes are added to meshes in the order in which they appear in the YAML
file. By default, each one replaces all materials that occupy the space
specified by its geometry file. This can be overridden by using the
:code:`replaces` and :code:`does_not_replace` properties.

.. code-block:: yaml

    dimensions: 3
    units: cm

    shapes:
      - name: wheel
        material: steel
        replaces: [rubber, air]
        geometry:
          format: stl
          path: wheel.stl
      - name: windshield
        does_not_replace: [steel]
        material: glass
        geometry:
          format: stl
          path: windshield.stl

In the example above, the wheel would only replace rubber and air. Any other
materials that happen to be in the same space as it would be left untouched.
The windshield would replace everything except steel. 

.. warning::

   It is an error to specify both :code:`replaces` and :code:`does_not_replace` 
   for the same shape.


Operators
---------

When assembling complex geometries, it is often the case that different parts
are specified in different coordinate systems. For example, a description
of the wheel of a car might be specified around its center, not its position
relative to the rest of the car. To help with this, Klee provides a mechanism
to apply transformations to shapes.

.. code-block:: yaml

    dimensions: 3

    shapes:
      - name: windshield
        material: glass
        geometry:
          format: stl
          path: windshield.stl
          units: cm
          operators:
            - rotate: 90
              axis: [0, 1, 0]
              center: [0, 0, -10]
            - translate: [10, 20, 30]

In the example above, the wheel is rotated 90 degrees counterclockwise
around an axis centered at the point :code:`(0, 0, -10)` and pointing in the
direction of the vector :code:`(0, 1, 0)`. It is then translated by the
vector :code:`(10, 20, 30)`.

Regardless of whether the geometry file has embedded units and what those may
be, units must be specified whenever specifying operators. These are the
units that will be used to interpret any lengths and points specified in
operators. Units may be specified in one of two ways: by specifying
:code:`units`, or by specifying :code:`start_units` and :code:`end_units`.
Specifying :code:`units` is the same as giving the same value for
:code:`start_units` and :code:`end_units`. Being able to change units is
useful for situations where your geometry file is in one set units, but
you're thinking about your larger assembly in another set of units.
For example:

.. code-block:: yaml

    dimensions: 3

    shapes:
      - name: windshield
        material: glass
        geometry:
          format: stl
          path: windshield.stl
          start_units: in
          end_units: ft
          operators:
            # Orient the windshield about its own coordinate system,
            # working in its native units (inches)
            - translate: [10, 20, 30]  # inches
            - rotate: 90
              axis: [0, 1, 0]
              center: [0, 0, -10]  # inches
            # switch to feet to put in the right place while thinking of
            # of the car in car in different units
            - convert_units_to: ft
            - translate: [2, 3, 4]  # feet

It is an error if the :code:`end_units` do not match the units after the
last operator.

Supported Operators
*******************
The supported operators are listed below. Unless otherwise specified,
the only difference between the 2D and 3D versions are that whenever points
or vectors are expected, the points and vectors must be of the dimensionality
specified by the shape file.

Operators take the form of :code:`operator_name: value`, where
:code:`operator_name` is the name of the operator, and
:code:`value` is the value specifying the parameters of the operation.
Operators may also have additional required or optional parameters.

* Translations

  :description: Translate the shape by a given vector.
  :name: :code:`translate`
  :value: a vector specifying the amount by which to translate the shape
  :example:
    ::

        # Translate by vector (1, 2, 3)
        translate: [1, 2, 3]

* Rotations

  :description: Rotate the shape by a given amount around a specified axis
  :name: :code:`rotate`
  :value: an angle, in degrees by which the shape will be rotated
    counterclockwise.
  :additional required parameters:
    :axis: (3D only) the nonzero axis of rotation
  :optional arguments:
    :center: a point specifying the center of rotation
  :example:
    ::

        # Rotate 45 degrees counterclockwise around the ray passing through
        # the point (1, 2, 3) and pointing in the direction of the vector
        # (4, 5, 6)
        rotate: 45
        center: [1, 2, 3]
        axis: [4, 5, 6]

* Scaling

  :description: Scale the shape by a specified amount
  :name: :code:`scale`
  :value: a vector specifying the amount by which to scale in each dimension,
    or a single value specifying by which to scale in all dimensions
  :optional arguments:
    :center: a point specifying the center relative to which to scale.
      If omitted, scaling is performed relative to the origin.
  :example:
    ::

        # Scale by 2x in the x direction 0.5x in y, and 1.5x in z
        scale: [2.0, 0.5, 1.5]

        # Scale by 2x in every direction relative to the point (1, 2, 3)
        scale: 2.0
        center: [1, 2, 3]

* Changing Units

  :description: Change the units in which subsequent operators are expressed.
    This is the same as scaling by the appropriate factor.
  :name: :code:`convert_units_to`
  :value: the name of the units to convert to. Must be one of the named units.
  :example:
    ::

      geometry:
        ...
        start_units: in
        end_units: cm
        operators:
          - translate: [2, 3]  # in inches
          - convert_units_to: cm  # same as scale: 2.54
          - translate: [5, 6]  # in centimeters

* Slices

  :description: Slice a 3D object and convert it into a 2D object. This is
    accomplished by defining a cut plane which will be used to determine
    what slice of the geometry to take. In addition, a point on the plane
    is picked as the new origin, and a vector is used to specify how the
    plane should be oriented with the 2D axes.
  :name: :code:`slice`
  :value: an object with the the following properties

    :origin: the point to use as the origin of the new coordinate system
    :normal: a vector normal to the slice plane
    :up: a vector which will be mapped to the positive Y direction on the cut plane.
  :optional arguments:
    :x: a single value specifying that the cut plane perpendicular to the
      x-axis at this value. See defaults table below.
    :y: a single value specifying that the cut plane perpendicular to the
      y-axis at this value. See defaults table below.
    :z: a single value specifying that the cut plane perpendicular to the
      z-axis at this value. See defaults table below.

    If a plane is specified by just giving "x", "y", or "z", then the origin,
    normal, and up vectors are given the default values specified
    in the table below. They can be overridden so long as the origin is still
    on the plane, and the normal is a multiple of the default normal.

    +------------------+-----------------------+-------------------+-------------------+
    | Usage            | :code:`origin`        | :code:`normal`    | :code:`up`        |
    +==================+=======================+===================+===================+
    | :code:`x: <val>` | :code:`(<val>, 0, 0)` | :code:`(1, 0, 0)` | :code:`(0, 0, 1)` |
    +------------------+-----------------------+-------------------+-------------------+
    | :code:`y: <val>` | :code:`(0, <val>, 0)` | :code:`(0, 1, 0)` | :code:`(1, 0, 0)` |
    +------------------+-----------------------+-------------------+-------------------+
    | :code:`z: <val>` | :code:`(0, 0, <val>)` | :code:`(0, 0, 1)` | :code:`(0, 1, 0)` |
    +------------------+-----------------------+-------------------+-------------------+

  :example:
    ::

        # Cut a 3D object with a plane that passes through the point
        # [10, 20, 30] and is normal to the vector [4, 5, 6]. The vector
        # [-5, 4, 0] will be mapped to the positive Y axis. [10, 20, 30] will
        # be mapped to the origin.
        slice:
          origin: [10, 20, 30]
          normal: [4, 5, 6]
          up: [-5, 4, 0]

Named Operators
***************
It can often be useful to name and reuse operators. For example, you may
have several parts of an assembly specified in one coordinate system
that you then need to transform to another. To enable reuse, we provide
support for named operators.

Named operators are specified via the top-level :code:`named_operators`
object. This is a list where each entry has the following values:

:name (required): the name of the operator. This is how it is referenced later.
:value (required): A list of operators. This is identical to the
  :code:`operators` entry in the :code:`geometry` object of a :code:`shape`.
:start_dimensions (optional): the number of initial dimensions of the
  operator. Must be 2 or 3. If not specified, the number of dimensions of the
  document is used.
:units (required, must specify this or start_units and end_units): the units in
  which the operator is specified
:start_units (optional, must specify this or units): the units in which the
  first operator is specified
:end_units (optional, must specify this or units): the units in which the
  last operator is specified. It is an error if the units aren't properly
  converted to `end_units` after applying all operations.

For Lua input, named-operator values support the same callback-capable fields
as shape operators. Klee constructs named operators before shapes, evaluates
each of their callbacks once, and shares that concrete result through every :code:`ref`.
A named-operator callback therefore cannot depend on the identity of a shape
that later refers to it.

The example below demonstrates how to create and then use a named operator.
Notice how we can use multiple :code:`ref` entries in the list of
operators and we can intermix these with other operators as needed.

.. code-block:: yaml

    dimensions: 2

    named_operators:
      - name: MyFirstOperator
        units: cm
        value:
          - translate: [10, 20]
          - scale: 1.5
      - name: AnotherOperator
        unit: cm
        value:
          - translate: [30, 40]

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: c2c
          path: wheel.contour
          units: cm
          operators:
            - ref: MyFirstOperator
            - rotate: 90
            - ref: AnotherOperator

An important thing to note is that the units of the named operator and
the shape that uses it do not have to match. The appropriate conversions
will be done automatically if needed. This allows you to not worry about how
the transformation was defined when you use it.

.. code-block:: yaml

    dimensions: 2

    named_operators:
      - name: MySampleOperator
        start_units: cm
        end_units: mm
        value:
          - translate: [10, 20]  # cm
          - convert_units_to: mm
          - translate: [30, 40]  # mm
    shapes:
      - name: wheel
        material: steel
        geometry:
          format: c2c
          path: wheel.contour
          units: in
          operators:
            # Automatic conversion from in to cm
            - ref: MySampleOperator
            # Automatic conversion from mm to in
            - translate: [50, 60]  # in



In addition to using :code:`ref` in an individual shape's operators, you
can also use it in other named operators. The only restriction is that it
be defined in the list before it is used. For Lua input, this restriction also
applies when a :code:`ref` callback returns the operator name.

.. code-block:: yaml

    dimensions: 2
    units: cm

    named_operators:
      - name: SomeOperator
        units: cm
        value:
          - translate: [10, 20]
          - scale: 1.5
      - name: AnotherOperator
        units: cm
        value:
          - rotate: 90
          - ref: SomeOperator
