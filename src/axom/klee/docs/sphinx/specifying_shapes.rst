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

Lua Input Decks
***************
Klee can also read Lua input decks when Axom is configured with
:code:`AXOM_ENABLE_LUA=ON` and Sol support. Lua decks use the same Klee schema
as YAML decks in this stage; the difference is that Lua evaluates first and
then Inlet reads the resulting global tables. File-based reads infer the input
format from :code:`.yaml`, :code:`.yml`, or :code:`.lua` extensions. Stream-based
reads are YAML by default unless the caller passes
:code:`axom::klee::InputFormat::Lua`.

Configuration
^^^^^^^^^^^^^
Lua deck support is optional. Downstream code can check for it with the
:code:`AXOM_USE_LUA` preprocessor macro from :code:`axom/config.hpp`. When Lua
support is unavailable, Klee rejects Lua input decks with an actionable error
instead of falling back to YAML parsing.

Migration from YAML
^^^^^^^^^^^^^^^^^^^
There is no automated YAML-to-Lua converter. A practical migration is to copy the
YAML structure directly, replace :code:`key: value` syntax with Lua table fields,
and then introduce :code:`local` constants or helper functions only where they
remove duplication or express a real computed choice.

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

Because Lua is evaluated before Klee reads the deck, ordinary table values can
be generated programmatically:

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

Caller-provided input variables can also be injected into a Lua deck before it
is evaluated. This is useful when an application wants one deck to select between
2D and 3D geometry, dimensions, or transforms at run time:

.. code-block:: c++

    axom::klee::InputVariables vars {
      {"dimensions", axom::klee::InputVariableValue {2}},
      {"shape_suffix", axom::klee::InputVariableValue {std::string {"2d"}}}
    };
    auto shapeSet = axom::klee::readShapeSet("shape.lua", vars);

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

Input variables are Lua-only and may be booleans, integers, doubles, or strings.
Their names must be Lua identifiers. They are globals by construction and are
allowed by Klee's unexpected-global check; other helper values in the deck should
still be declared :code:`local`.

Use :code:`local` helper functions and constants for intermediate values so the
global namespace contains only the Klee schema fields that Inlet should read.
For Lua input, a one-value scale is written as a one-entry table, for example
:code:`{ scale = {2.0} }`.

Parse-Time Callbacks
^^^^^^^^^^^^^^^^^^^^
Selected operator fields may also be written as zero-argument Lua callbacks.
Klee evaluates each callback exactly once while reading the deck; the resulting
shape still contains ordinary affine or slice operators, not runtime Lua
functions. Callbacks should be pure functions of local deck variables.

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

Vector-valued callbacks return raw numeric Lua tables such as :code:`{x, y}` or
:code:`{x, y, z}`. The typed :code:`Vector.new(...)` object is also accepted.
Scalar-valued callbacks return a number. Supported callback fields are
:code:`translate`, :code:`axis`, :code:`center`, :code:`scale`,
:code:`slice.origin`, :code:`slice.normal`, :code:`slice.up`, :code:`rotate`,
:code:`slice.x`, :code:`slice.y`, and :code:`slice.z`. For :code:`scale`, a
number means uniform scaling and a table means per-axis scaling.

Runtime Point Transforms
^^^^^^^^^^^^^^^^^^^^^^^^
Lua decks can also use the runtime :code:`transform` operator for non-affine
point mapping. Unlike zero-argument callbacks, :code:`transform` is stored with
the geometry and evaluated for each point when a consumer applies shape
operators. The transform function receives a point with :code:`x`, :code:`y`,
and :code:`z` fields and returns a raw numeric table with the same dimensionality
as the input geometry.

.. code-block:: lua

    local function spherical_to_cartesian(p)
      local r = p.x
      local theta = p.y
      local phi = p.z
      return {
        r * math.sin(phi) * math.cos(theta),
        r * math.sin(phi) * math.sin(theta),
        r * math.cos(phi)
      }
    end

    dimensions = 3

    shapes = {
      {
        name = "shell",
        material = "steel",
        geometry = {
          format = "stl",
          path = "shell_in_spherical_coords.stl",
          units = "cm",
          operators = {
            { transform = spherical_to_cartesian }
          }
        }
      }
    }

The :code:`transform` operator is forward-only and cannot be represented as a
matrix. Calls to matrix-only APIs such as :code:`Geometry::getTransform()` reject
geometries that contain :code:`transform`; use :code:`Geometry::applyTransform()`
for point-wise evaluation. Quest's discrete-shape path applies mixed affine and
Lua transforms in input order, while Quest algorithms that require affine
matrices reject non-affine transforms with a diagnostic.

Lua API Reference
^^^^^^^^^^^^^^^^^
Lua decks may use ordinary Lua variables and functions, but only Klee schema
fields should be global. Klee/Inlet registers the typed :code:`Vector` object as
an alternative to raw tables. For callbacks and runtime transforms, raw numeric
tables are the recommended return form. Runtime transform arguments expose
:code:`p.x`, :code:`p.y`, :code:`p.z`, and :code:`p.dim`.

Limitations
^^^^^^^^^^^
Runtime :code:`transform` operators are forward-only and do not provide inverse transforms.

**Lua State Lifetime:** Geometry objects with Lua :code:`transform` operators
hold shared ownership of the Lua state used during parsing. The Lua state is
automatically kept alive for the lifetime of these geometries through
:code:`std::function` capture semantics, so no manual lifetime management is
required. Geometries with Lua transforms can be safely used after the original
:code:`readShapeSet` call returns and any parsing artifacts are destroyed.

**Serialization:** The original Lua deck remains the portable representation;
runtime Lua transforms cannot be serialized to Sidre/Conduit as pure data.
Attempting to serialize a geometry with non-affine operators will result in
an error. For reproducibility, retain the original :code:`.lua` input file.

Error Messages
^^^^^^^^^^^^^^
Common Lua input errors are reported as Klee parsing errors. A Lua deck read
without Lua support reports:

.. code-block:: text

    Lua input decks require Axom configured with AXOM_ENABLE_LUA=ON and Sol library support. Rebuild Axom with Lua enabled or convert deck to YAML.

Unsupported file extensions are rejected before parsing, and unexpected globals
or syntax errors are reported during Inlet verification or Lua evaluation.
Callback failures include the field, shape name when available, and operator
location, for example:

.. code-block:: text

    Error evaluating callback for 'translate' in shape 'part' operator 1: [Inlet] Lua function call failed: ...
Runtime transform failures are reported when the transform is evaluated and carry
the same shape/operator context. Matrix-only consumers report that
:code:`transform` is non-affine and point-wise evaluation is required.

Paths
*****
The paths specified in shapes are specified either as absolute paths
(begin with a :code:`/`), or as relative paths. Absolute paths are evaluated
as absolute paths in the file systems. Relative paths are evaluated relative
to the Klee input deck (not the current working directory). For example, in the
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
    :axis: (3D only) the axis of rotation
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
be defined in the list before it is used.

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
