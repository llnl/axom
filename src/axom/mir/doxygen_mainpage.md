MIR {#mirtop}
=============

Axom's [Mir](@ref axom::mir), (Material Interface Reconstruction) component
provides classes that can process Blueprint meshes with mixed materials into new
Blueprint meshes such that there is a single material per zone. In addition, the
MIR component provides useful building blocks for developing algorithms using
Blueprint meshes. There are views that simplify dealing with Conduit data and
utility algorithms for processing and constructing meshes.

# Design goals {#mirgoals}

This component's algorithms are mainly delivered as classes that are templated on
an execution space, allowing them to operate on a variety of computing backends.
The algorithms take Conduit nodes (containing Blueprint data) as input and
output new Blueprint data in an output Conduit node. Where possible, algorithms
have been broken out into classes to promote reuse.

# Views {#mirviews}

Blueprint defines protocols for representing various mesh and data constructs in
a hierarchical form inside Conduit nodes. There are objects defined for coordinate
sets (coordsets), mesh topologies, mesh fields, material sets (matsets), and there
are various flavors of each type of object. This can make it difficult to write
algorithms against Blueprint data since the data live in Conduit nodes with different
names and they may have different formats. Conduit can also use multiple data types
for any of the data arrays that represent objects. Views were developed to simplify
some of these challenges by providing common templated interfaces that let different
types of objects be accessed in a uniform way. Templating helps to deal with the
data types. The views also provide functions that can wrap a Conduit node
in a suitable view type and then dispatch that view to a generic user-provided lambda,
enabling algorithms to be instantiated for multiple data types with a compact amount
of code.

# MIR {#mir}

The MIR component provides an implementation of the Equi-Z MIR algorithm, though
additional algorithms are planned.

# Utilities {#mirutilities}

The MIR component provides algorithms for performing useful mesh operations such as
extracting sub-meshes, merging meshes, clipping meshes, and creating relations.
These building blocks can be reused to ease the process of writing additional algorithms
that operate on Blueprint meshes.

