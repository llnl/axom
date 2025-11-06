[comment]: # (#################################################################)
[comment]: # (Copyright 2017-2025, Lawrence Livermore National Security, LLC)
[comment]: # (and Axom Project Developers. See the top-level LICENSE file)
[comment]: # (for details.)
[comment]: #
[comment]: # (# SPDX-License-Identifier: BSD-3-Clause)
[comment]: # (#################################################################)


<style>
    figcaption {
        text-align: center;
        font-style: italic;
        color: #555555;
    }
</style>

# Lesson 03: Defining geometry setup with Klee

In this lesson, we temporarily set aside our mesh metadata example and use Axom's `Klee` component to define the geometric setup for a multimaterial simulation.  (Don't worry, we'll return to the mesh metadata in the next lesson).

`Klee` is built on top of `Inlet` to define the schema for geometry setup. A Klee input consists of a list of shapes; each shape specifies its material and geometry, and may optionally include *replacement rules* that describe which previously "shaped in" materials the current shape will replace or preserve.

## Introducing Klee

<div style="text-align: center;">

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#e6f0ff",
    "primaryBorderColor": "#3366cc",
    "primaryTextColor": "#000000",
    "lineColor": "#3366cc",
    "background": "#ffffff",
    "fontFamily": "Verdana"
  }
}}%%
graph LR
  classDef highlight fill: #ffe082, stroke: #e65100, stroke-width: 3px;

  bump --> spin
  inlet --> primal
  inlet --> sidre
  klee:::highlight --> inlet
  lumberjack --> core
  mint --> slam
  mir --> bump
  multimat --> slam
  sidre --> slic
  primal --> slic
  quest --> spin
  quest --> mint
  sina --> core
  slam --> core
  spin --> slam
  spin --> primal

  %% dashed (soft deps)
  mint -.-> sidre
  quest -.-> klee
  slic -.-> lumberjack

  primal --> core
  sidre --> core
  slam --> slic

  %% links can be highlighted as follows:
  %% linkStyle 5 stroke: #d32f2f,color: #d32f2f, stroke-width:3px;
```
<figcaption>Figure: Axom components, highlighting Klee</figcaption>
</div>

Klee provides a flexible framework for defining complex geometric configurations through its key components:

**ShapeSet**
: The top-level container that holds all shapes. It manages:
   - The global dimensionality for shapes (2D or 3D); individual shapes may override it if needed.
   - The ordered collection of shapes

**Shape**
: Represents a single geometric entity with an assigned material. Each shape includes:
   - A unique name identifier
   - A material designation
   - A geometry definition
   - Replacement rules (how it interacts with previously defined materials)

**Geometry**
: Defines the geometry of a shape, often a path to a file, as well as operators defining transformations that can be applied to the shape, e.g. scaling, rotating, unit conversions.

Shapes are processed in order, with later shapes potentially replacing or preserving materials from earlier shapes according to the replacement rules. This approach allows for building complex material layouts through a combination of simple shapes and clearly defined interaction rules.

> :information_source: **Why "Klee"?**
> The Klee component was named after the Swiss-German artist Paul Klee (1879-1940), who was known for his distinctive style featuring geometric shapes, abstract forms, and colorful compositions. The pronunciation of Klee (/kleɪ/, often pronounced “clay”) evokes modeling or shaping materials.
>
> <div style="text-align: center;">
>   <img src="klee_paintings.png" alt="Examples of Paul Klee's paintings featuring geometric shapes and patterns" />
>   <figcaption>Figure: Examples of Paul Klee's paintings featuring geometric shapes and patterns. <br />
>               Public domain images from: https://en.wikipedia.org/wiki/Paul_Klee
>   </figcaption>
> </div>




### Separation of geometric data and replacement rules
Separating the Klee input from the geometric files creates a clear division between the "business logic" of shape ordering and replacement rules and the actual geometric descriptions of shapes.

This separation provides several key benefits:
  - Geometric shapes can be reused across different simulations with different material assignments
  - Material replacement logic and geometric definitions can evolve independently
  - Simulation setup can be quickly iterated by changing replacement rules without modifying complex geometry files
  - Domain experts can focus on creating accurate geometries while simulation scientists define how materials interact
  - Version control becomes more manageable with separate geometry files and material assignment rules

<div style="text-align: center;">
  <p>
    <img src="jack_rects.svg" width="20%" alt="iteration 1: two rectangles">
    <img src="jack_rects_and_circs.svg" width="20%" alt="iteration 2: add circular caps">
    <img src="jack_capped_rects.svg" width="20%" alt="iteration 3: two circular caps">
    <img src="jack_single_path.svg" width="20%" alt="iteration 4: a single path">
  </p>
  <figcaption>Figure: Different design iterations of a "jack" shape that could be used in Klee geometry definitions. These changes can be applied without changing the Klee input file.</figcaption>
</div>


### Geometry definition: 

#### File formats
- Klee currently supports the following input mesh formats: 
    - `.stl` triangle mesh (ascii and binary)
    - `.c2c` contour file
    - `.mfem` contours
    - `.proe` tetrahedral meshes

> :warning: **Note:** Support for `.c2c` files is only available on LLNL's LC systems. For other environments, please use alternative geometry file formats.

> :information_source: Our ascii-based `.proe` format is an indexed representation for the vertices and tetrahedra within a tetrahedral mesh.

> :information_source: Axom provides a Python script to convert SVG files to the MFEM format, making it easy to use vector graphics as geometry inputs for Klee.

<details>
  <summary>STL Format Example (four triangles bounding a tetrahedron)</summary>

<div style="text-align: center;">

<div style="text-align: center;">
<table style="margin: auto;">
  <tr>
    <td style="text-align: left; vertical-align: top;">

```stl
solid Mesh
facet normal 0.0 0.0 -1.0
  outer loop
    vertex 0.0 0.0 0.0
    vertex 1.0 0.0 0.0
    vertex 0.0 1.0 0.0
  endloop
endfacet
facet normal -1.0 0.0 0.0
  outer loop
    vertex 0.0 0.0 0.0
    vertex 0.0 1.0 0.0
    vertex 0.0 0.0 1.0
  endloop
endfacet
facet normal 0.0 -1.0 0.0
  outer loop
    vertex 0.0 0.0 0.0
    vertex 0.0 0.0 1.0
    vertex 1.0 0.0 0.0
  endloop
endfacet
facet normal 0.577 0.577 0.577
  outer loop
    vertex 1.0 0.0 0.0
    vertex 0.0 0.0 1.0
    vertex 0.0 1.0 0.0
  endloop
endfacet
endsolid Mesh
```
</td>
<td>

[Open interactive 3D viewer on GitHub ↗](https://github.com/LLNL/axom_data/blob/main/quest/tetrahedron.stl)

</td>
</tr>
</table>
</div>

</div>
</details>

#### Basic shape file

A minimal shape file includes a ``dimensions`` field and a ``shapes`` list.

The following example has a single shape representing the boundary a tetrahedron.
The name of the shape -- in this case ``my_tetrahedron`` -- is used internally, and it will get shaped into the ``steel`` material:

```yaml
dimensions: 3

shapes:
  - name: my_tetrahedron
    material: steel
    geometry:
      format: stl
      path: tetrahedron.stl
      units: cm
```


#### Operators

The geometry definition can also include `operators` defining affine transformations (scaling, rotations, translations, ...) and unit conversions.

<div style="text-align: center;">
<table style="border-collapse: collapse; border: none;">
  <tr style="border: none;">
    <td style="width:40%; border: none;">

```yaml
dimensions: 2

shapes:
  - name: outer_shell
    material: steel
    geometry:
      format: c2c
      path: ../contours/unit_circle.contour
      units: cm
      operators:
        - scale: 5
  - name: inner_ball
    material: void
    geometry:
      format: c2c
      path: ../contours/unit_circle.contour
      start_units: mm
      end_units: cm
      operators:
        - scale: 25
        - convert_units_to: cm
```

  </td>
  <td style="width:40%; border: none;">

```yaml
dimensions: 2

shapes:
  - name: front_left_wheel
    material: steel
    geometry:
      format: stl
      path: wheel.stl
      units: cm
      operators:
        - rotate: 90
        - translate: [100, -80, 0]
  - name: front_right_wheel
    material: steel
    geometry:
      format: stl
      path: wheel.stl
      units: cm
      operators:
        - rotate: -90
        - translate: [100, 80, 0]
```

  </td>
  </tr>
</table>
  <figcaption style="font-style: italic;">Figure: Example Klee inputs showing `scale`, `translate`, `rotate` and unit conversion operators.</figcaption>
</div>


### Replacement Rules
Replacement rules give users some extra control in how shapes get overlaid. By default, a new shape of a given material will replace all other shapes.

<div style="text-align: center;">
  <img src="klee_replacement_use_cases.png" width="50%" alt="Several use cases for replacement rules" />
  <figcaption>Figure: Replacement rules have many uses, including when there are overlapping parts (top), when we need to expand a shape to close a gap (middle), or when we need to fill a void (bottom).</figcaption>
</div>

<br />

If desired, users can either add an explicit list of materials to replace via the `replaces` entry, or an explicit list of materials to preserve via the `does_not_replace` entry (but not both).

<div style="text-align: center;">
  <img src="klee_replacement_rules.png" width="75%" alt="Visualization of Klee replacement rules">
  <figcaption style="font-style: italic;">Figure: Illustration of Klee replacement rules (default, explicit "replaces" list, and explicit "does_not_replace" list) and how later shapes interact with earlier materials when shaping "wood" boats on a mesh with "water", "mud" and "grass".</figcaption>
</div>

<br />

<details>
<summary>More details on input structure</summary>

- Top-level
  - dimensions: 2 | 3
  - shapes: array of shape entries
- Shape
  - name: unique string
  - material: string
  - geometry
    - format: none | stl | c2c | mfem | proe
    - path: relative path to geometry file
    - start_units/end_units/units: unit metadata
    - operators: ordered list of transforms/conversions
      - translate: [dx, dy, dz]
      - rotate: ang
        axis: [ax, ay, az]
        center: val
      - scale: scalar or [sx, sy, sz]
      - convert_units_to: target_unit
  - replaces: \<list of material names>
  - does_not_replace: \<list of material names>
</details>

## Examples: 

### Simple 2D contours with units and scaling

```yaml
dimensions: 2

shapes:
  - name: outer_shell
    material: steel
    geometry:
      format: c2c
      path: ../contours/unit_circle.contour
      units: cm
      operators:
        - scale: 5
  - name: inner_ball
    material: void
    geometry:
      format: c2c
      path: ../contours/unit_circle.contour
      start_units: mm
      end_units: cm
      operators:
        - scale: 25
        - convert_units_to: cm
```

<details>
<summary> Another example w/ STL meshes </summary>

```yaml
dimensions: 3
units: m
shapes:
  - name: fluid_volume
    material: water
    geometry:
      format: stl
      path: ../meshes/volume.stl
      operators:
        - scale: 0.1
  - name: rotor_blade
    material: composite
    geometry:
      format: stl
      path: ../meshes/blade.stl
      operators:
        - rotate: 45
            axis: [0, 1, 0]
           angle: 30
        - translate: [2.0, 0.0, 0.0]
```
</details>

## Let's see some code!

The code example for this lesson loads a Klee file, performs some validation and then prints out details about the geometric setup

### Load and validate the Klee input

```cpp
axom::klee::ShapeSet shapeset
try
{
  shapeSet = axom::klee::readShapeSet(inputFilename);
}
catch(axom::klee::KleeError& error)
{
  std::vector<std::string> errs;
  for(auto verificationError : error.getErrors())
  {
    errs.push_back(axom::fmt::format(" - '{}': {}",
                                      static_cast<std::string>(verificationError.path),
                                      verificationError.message));
  }

  SLIC_WARNING(
    axom::fmt::format("Error during parsing klee input. Found the following errors:\n{}",
                      axom::fmt::join(errs, "\n")));
}
```

Next, we loop through the shapes and print out information about each shape. We're using an `fmt::memory_buffer` (similar to a `std::stringstream`) to write everything in a single log statement:
```cpp
  axom::fmt::memory_buffer buffer;
  axom::fmt::format_to(std::back_inserter(buffer), "Klee ShapeSet Information:\n");

  axom::fmt::format_to(std::back_inserter(buffer),
                       "  Overall dimensions: {}\n",
                       dimensionsToString(shapeSet.getDimensions()));

  // Step 1: Collect and print sorted list of material names
  std::vector<std::string> materials = [&shapeSet]() {
    std::set<std::string> materialSet;
    for(const auto& shape : shapeSet.getShapes())         // <-- 1(a)
    {
      materialSet.insert(shape.getMaterial());            // <-- 1(b)
    }
    return std::vector<std::string>(materialSet.begin(), materialSet.end());
  }();
  axom::fmt::format_to(std::back_inserter(buffer),
                       "\n  Unique materials ({}): {}\n",
                       materials.size(),
                       materials);

  // Step 2: Print information about each shape
  fmt::format_to(std::back_inserter(buffer),
                 "\n  Details for the {} shapes:\n",
                 shapeSet.getShapes().size());
  for(const auto& shape : shapeSet.getShapes())
  {
    const auto& geom = shape.getGeometry();               // <-- 2(a)

    axom::fmt::format_to(
      std::back_inserter(buffer),
      "  - name: '{}'\n"
      "    material: '{}'\n"
      "    format: '{}'\n"
      "    units: {}\n"
      "    dimensions: {}\n"
      "    replaces materials: {}\n\n",
      shape.getName(),                                   // <-- 2(b)
      shape.getMaterial(),
      geom.getFormat(),
      formattedUnits(geom.getStartProperties().units, geom.getEndProperties().units),
      formattedDimensions(geom.getInputDimensions(), geom.getOutputDimensions()),
      getMaterialsReplacedBy(shape, materials));
  }

  SLIC_INFO(axom::fmt::to_string(buffer));
```

## Ice cream challenge

Let's create a setup for an ice cream cone, which will consist of a cone, and ice cream scoop and a bunch of sprinkles.

<div style="text-align: center;">
<p>
  <img src="ice_cream.svg" width="20%" alt="Ice cream cone illustration used for the challenge">
  <img src="cone.svg" width="20%" alt="Cone">
  <img src="scoop.svg" width="20%" alt="Scoop">
  <img src="sprinkles.svg" width="20%" alt="sprinkles">
  <figcaption style="font-style: italic;">Figure: Ice cream cone geometry for the challenge.</figcaption>
</div>

<br />

**Task:** Given the "cone", "scoop" and "sprinkle" geometries, create a Klee setup such that:
* the cone covers the background
* the scoop covers the cone
* the sprinkles cover the ice cream scoop, but not the cone or the background

**Solution:**
```yaml
dimensions: 2

shapes:
  - name: background
    material: background
    geometry:
      format: none

  - name: vanilla_scoop
    material: icecream
    geometry:
      format: mfem
      path: ice_cream_scoop.mesh
      units: cm

  - name: colorful_sprinkles
    material: sprinkles
    geometry:
      format: mfem
      path: ice_cream_sprinkles.mesh
      units: cm
    replaces: [icecream]

  - name: cone
    material: batter
    geometry:
      format: mfem
      path: ice_cream_cone.mesh
      units: cm
    does_not_replace: [icecream, sprinkles]
```

<details>
  <summary>Sidebar: Converting SVG to MFEM format</summary>

  Axom has an `svg2contours` python script to convert SVG files to the [MFEM NURBS format](https://mfem.org/mesh-format-v1.0/#nurbs-meshes) using the [svgpathtools](https://github.com/mathandy/svgpathtools) python package. We use a virtual environment to set up the proper dependencies (and add a small patch for svgpathtools):

  ```bash
  # create a virtual environment with svgpathtools and dependencies
  > cd <axom>/src/tools/svg2contours/
  > python3 -m venv venv
  > source venv/bin/activate

  (venv)> pip3 install -r requirements.txt
  # Note: We apply a small patch that was upstreamed after v1.7.1
  (venv)> patch  -p1 venv/lib/python3.9/site-packages/svgpathtools/path.py -i svgpathtools-1.7.1-eigenvec-fix.patch --verbose

  # run the script on our svg images: sprinkles.svg, scoop.svg and cone.svg
  (venv)> cd <build_directory>
  (venv)> python3 ../src/tools/svg2contours/svg2contours.py -i ../src/examples/shaping_tutorial/lesson_03/sprinkles.svg -o sprinkles.mesh

    SVG dimensions: width='416.53091' height='506.883' viewBox='0 0 416.53091 506.88299'
    Wrote 'sprinkles.mesh' with 608 vertices and NURBS 304 elements
  
  (venv)> python3 ../src/tools/svg2contours/svg2contours.py -i ../src/examples/shaping_tutorial/lesson_03/cone.svg -o cone.mesh

    SVG dimensions: width='416.53091' height='506.883' viewBox='0 0 416.53091 506.88299'
    Wrote 'cone2.mesh' with 12 vertices and NURBS 6 elements
  
  (venv)> python3 ../src/tools/svg2contours/svg2contours.py -i ../src/examples/shaping_tutorial/lesson_03/scoop.svg -o scoop.mesh

    SVG dimensions: width='416.53091' height='506.883' viewBox='0 0 416.53091 506.88299'
    Wrote 'scoop2.mesh' with 22 vertices and NURBS 11 elements
  ```

  <div style="text-align: center;">
    <img src="ice_cream_visit.png" width="60%" alt="Viewing the generated MFEM files for the ice cream cone in VisIt using the MultiresControl">
    <figcaption>Figure: Generated MFEM meshes containing the contours for our ice cream example (cone, scoop, sprinkles) visualized in VisIt. We use VisIt's `MultiresControl` operator to see the curvature of the contours.</figcaption>
  </div>

</details>

> :clapper: Let's run this example:
>
> ```bash
> > ./bin/lesson_03_klee_operators_and_validation ../lesson_03/ice_cream.yaml
>```
> <details>
> <summary> Output </summary>
> 
> ```bash
> [INFO] Klee ShapeSet Information:
>   Overall dimensions: 2D
> 
>   Unique materials (4): ["background", "batter", "ice_cream", "sprinkles"]
> 
>   Details for the 4 shapes:
>   - name: 'background'
>     material: 'background'
>     format: 'none'
>     units: '<unspecified>'
>     dimensions: 2D
>     replaces materials: ["background", "batter", "ice_cream", "sprinkles"]
> 
>   - name: 'vanilla_scoop'
>     material: 'ice_cream'
>     format: 'mfem'
>     units: 'cm'
>     dimensions: 2D
>     replaces materials: ["background", "batter", "ice_cream", "sprinkles"]
> 
>   - name: 'colorful_sprinkles'
>     material: 'sprinkles'
>     format: 'mfem'
>     units: 'cm'
>     dimensions: 2D
>     replaces materials: ["ice_cream"]
> 
>   - name: 'cone'
>     material: 'batter'
>     format: 'mfem'
>     units: 'cm'
>     dimensions: 2D
>     replaces materials: ["background", "batter"]
> ```
> </details>


<br />

  > :crystal_ball: We'll revisit what this will look like when "shaped-in" to a computational mesh in the next lesson

## Summary and next steps

In this lesson, we learned about geometric setup through Klee. 

In the next lesson, we'll explore Axom's shaping algorithms, which will generate the volume fractions on a computational mesh associated with a user-provided Klee input.
