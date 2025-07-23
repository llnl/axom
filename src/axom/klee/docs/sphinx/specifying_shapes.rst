Specifying Shapes
=================

"Shaping", or "painting", is the process of adding non-conformal material
regions to a mesh. Traditionally, this has been done in code-specific formats
by each code that provides such a capability. Axom's Klee library provides
a way to read shape specifications in YAML files and apply the specified
geometry to a mesh.

Basics
------
Shapes in Klee are specified in YAML. A basic file consists of a list of
shapes. Each one specifies its :code:`name` and :code:`material`,
as well as a description of its geometry.

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

Paths
*****
The paths specified in shapes are specified either as absolute paths
(begin with a :code:`/`), or as relative paths. Absolute paths are evaluated
as absolute paths in the file systems. Relative paths are evaluated relative
to the YAML file (not the current working directory). For example, in the
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
  :example:
    ::

        # Scale by 2x in the x direction 0.5x in y, and 1.5x in z
        scale: [2.0, 0.5, 1.5]

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
