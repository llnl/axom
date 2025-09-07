# Lesson 02: Inlet Metadata Example

This example demonstrates how to use Axom's Inlet library to define, parse, and validate structured mesh metadata from a YAML input file. The metadata includes a bounding box with minimum and maximum coordinates and a resolution in x and y directions.

## MeshMetadata Struct

- Contains nested structs `BoundingBox` and `Resolution`.
- `BoundingBox` defines the minimum `(min_x, min_y)` and maximum `(max_x, max_y)` coordinates.
- `Resolution` defines the number of elements in the x and y directions.

## Schema Definition

The `defineSchema` static method of `MeshMetadata` defines the Inlet schema, specifying:

- `bounding_box` struct with nested `min` and `max` structs containing `x` and `y` double parameters with default values.
- `resolution` struct containing integer parameters `x` and `y` with default values.

## FromInlet Specialization

A template specialization of `FromInlet` for `MeshMetadata` enables direct deserialization from Inlet’s Container.

## Input YAML File Format

```yaml
mesh:
  bounding_box:
    min:
      x: 0.0
      y: 0.0
    max:
      x: 1.0
      y: 1.5
  resolution:
    x: 15
    y: 25
```

## Building and Running the Example

Build the example using your usual build system, then run:

```bash
./inlet_metadata input.yaml
```

The program will parse the input YAML, validate the parameters, and print the mesh metadata to the console.

## Output Example

```
Bounding Box Min: (0.0, 0.0)
Bounding Box Max: (1.0, 1.5)
Resolution: (15, 25)
```
