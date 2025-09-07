# Setting up a 2D Cartesian Mesh with Axom Sidre and Inlet

## Introduction

This tutorial introduces how to use Axom's Sidre and Inlet components to set up and manage metadata for a simple 2D Cartesian mesh. Sidre provides an efficient way to store hierarchical mesh metadata, while Inlet facilitates flexible metadata input through parameter parsing. We will focus on defining the spatial bounding box and resolution of the mesh, closely following the specifications in the [Conduit mesh blueprint](https://llnl-conduit.readthedocs.io/en/latest/mesh.html).

A basic understanding of mesh concepts and the Conduit mesh blueprint is helpful for this tutorial.

Mesh metadata defines key properties that describe the geometry and discretization of a mesh. For Cartesian meshes, two main pieces of metadata are essential:

- **Bounding Box:** Defines the spatial extent of the mesh, described by minimum and maximum coordinates in each dimension.
- **Resolution:** Specifies the number of discretization points or cells along each dimension, dictating the mesh granularity.

<figure style="text-align: center;">
  <img src="cartesian_mesh.svg" alt="Cartesian Mesh" style="display: inline-block;" />
  <figcaption>Figure: This shows the resolution and bounding box for a 2D Cartesian mesh.</figcaption>
</figure>

In this tutorial, we will demonstrate how to represent these metadata elements using Sidre data structures and configure their input with Inlet.

## Sidre Basics

Sidre is an Axom component designed for managing hierarchical data structures stored in memory efficiently. It is well-suited for storing and organizing mesh metadata and simulation data. Sidre's core concept revolves around a hierarchical tree of data containers managed by a central `DataStore`.

### Key Concepts in Sidre

- **DataStore:** The top-level container that owns the entire hierarchical data structure. All groups and views ultimately belong to a DataStore instance.

- **Groups:** Nodes in the hierarchy that act like directories or folders. Groups can contain other groups or views, helping to organize data logically.

- **Views:** Leaf nodes containing metadata or raw data. Views provide access to actual data buffers.

- **Buffers:** Memory blocks allocated to hold the data referenced by views. Views use buffers to read and write actual data values.

- **Attributes:** Metadata about the views or groups, such as type information or external identifiers, which provide additional context.

Sidre allows flexible and efficient memory management, making it ideal to store structured mesh metadata such as bounding box coordinates and resolution parameters in an accessible and modifiable way.

In the following sections, we will create groups and views within a Sidre DataStore to store the bounding box and resolution data for a Cartesian mesh.


## Defining Mesh Metadata with Sidre

Now, let's set up the Sidre groups and views to store the mesh metadata for a 2D Cartesian mesh. We'll create a root group called `"mesh"` and add two subgroups: `"bounding_box"` and `"resolution"`.

```cpp
#include "axom/sidre.hpp"

void setup_mesh(axom::sidre::DataStore* datastore)
{
  // Create a root group for the mesh metadata
  axom::sidre::Group* meshGroup = datastore.getRoot()->createGroup("mesh");

  // Create bounding box groups and views
  axom::sidre::Group* minGroup = meshGroup->createGroup("bounding_box/min");
  minGroup->createViewScalar("x", input.min_x);
  minGroup->createViewScalar("y", input.min_y);

  axom::sidre::Group* maxGroup = meshGroup->createGroup("bounding_box/max");
  maxGroup->createViewScalar("x", input.max_x);
  maxGroup->createViewScalar("y", input.max_y);

  // Create resolution group and views
  axom::sidre::Group* resGroup = meshGroup->createGroup("resolution");
  resGroup->createViewScalar("x", input.res_x);
  resGroup->createViewScalar("y", input.res_y);
}
```

The hierarchy looks as follows:

<figure style="text-align: center;">
  <img src="hierarchy_diagram.svg" alt="Sidre Group and View Hierarchy" style="display: inline-block;" />
  <figcaption>Figure: Sidre hierarchical structure showing the mesh metadata groups and views.</figcaption>
</figure>


> :clapper: You can try running this example code in the [mesh_metadata_sidre](https://github.com/LLNL/axom/tree/develop/examples/mesh_metadata_sidre) example provided in Axom's GitHub repository to see how Sidre manages mesh metadata in practice. This example allows you to enter the bounding box and resolution parameters.


