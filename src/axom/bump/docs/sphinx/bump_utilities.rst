.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
BUMP Utilities
******************************************************

The BUMP component contains several useful building blocks for writing algorithms
for Blueprint meshes.

 * Structured as classes with an ``execute()`` method
 * Often templated on execution space and views

Classes are used so algorithms can contain state such as view objects and so various
algorithm stages or utility functions can be written as helper methods. Classes also
allow algorithms to be modified using inheritance and overriding methods. Finally,
classes allow the BUMP component to provide specialized versions of templated algorithms.

Blueprint meshes consist of various parts such as coordsets, topologies, fields, and matsets.
These constructs are organized as various paths within a root Conduit node. Elements such
as strings and scalar values can be stored as usual within a Conduit node in host memory.
Conduit nodes that contain bulk data such as coordinate or field data should point to
memory blocks that are valid for the execution environment in which algorithms will run.
To achieve this, entire Conduit node hierarchies can be copied to the proper memory space
or they can be constructed by forcing Conduit to allocate memory in the proper memory space.

 * `NOTE: Some code examples use a utils:: namespace prefix as shorthand for axom::bump::utilities.`

#######################
Copying Blueprint Data
#######################

If a ``conduit::Node`` containing Blueprint data is not in a memory space appropriate
for the target execution space, the data can be copied to a suitable memory space using
the ``axom::bump::utilities::copy<ExecSpace>()`` function. The target execution
space (and thus its memory space) is specified using the ``copy`` function's ``ExecSpace``
template argument. The ``copy`` function copies the source ``conduit::Node`` to
the destination ``conduit::Node``, making sure to use the appropriate Axom allocator for
non-string bulk arrays (e.g. arrays of ints, floats, doubles, etc.). Data small enough to
fit in a ``conduit::Node`` and strings are left in the host memory space, which lets
algorithms on the host side query them. For data that have been moved to the device,
their sizes and data types can still be queried using normal Conduit mechanisms such
as the ``conduit::Node::dtype()`` method.

.. code-block:: cpp

    conduit::Node hostMesh, deviceMesh, hostMesh2;
    // host->device
    axom::bump::utilities::copy<axom::HIP_EXEC<256>>(deviceMesh, hostMesh);
    // device->host
    axom::bump::utilities::copy<axom::SEQ_EXEC>(hostMesh2, deviceMesh);

############################
ConduitAllocateThroughAxom
############################

When writing algorithms that construct Blueprint data, it is helpful to force Conduit
to allocate its memory through Axom's allocation routines and then make an ``axom::ArrayView``
of the data in the Conduit node. This prevents data from having to be copied from an Axom
data structure into a Conduit node since it can be constructed from the start inside the
Conduit node. The size of the array must be known.

The ``axom::bump::utilities::ConduitAllocateThroughAxom``
class is a template class that takes an execution space as a template argument. The class
installs an allocation routine in Conduit that can be used to allocate data through
Axom. The Conduit allocator is set on each ``conduit::Node`` before setting data into
the object.

.. literalinclude:: ../../extraction/TableBasedExtractor.hpp
   :start-after: _bump_utilities_c2a_begin
   :end-before: _bump_utilities_c2a_end
   :language: C++

##########
ClipField
##########

The ``axom::bump::extraction::ClipField`` class intersects all the zones in the input Blueprint
mesh with an implicit surface where the selected input field equals zero and produces a new
Blueprint mesh based on the selected zone fragments produced by the intersection. This can be thought
of as an isosurface algorithm but with a volumetric output mesh where the mesh is either inside or
outside of the selected isovalue. The ``ClipField`` class has multiple template arguments to
select the execution space, the type of topology view, the type of coordset view, and the
type of intersector used to determine intersections. The default intersection uses an isosurface-
based intersection method, though other intersectors could be created to perform plane
or sphere intersections.

.. literalinclude:: ../../tests/bump_clipfield.cpp
   :start-after: _bump_utilities_clipfield_begin
   :end-before: _bump_utilities_clipfield_end
   :language: C++

################
CoordsetBlender
################

The ``axom::bump::CoordsetBlender`` class takes a ``BlendData`` and makes
a new explicit coordset where each new point corresponds to one blend group. A "BlendData" is
an object that groups several array views that describe a set of blend groups. Each blend group
is formed from a list of node ids and weight values. A new coordinate is formed by looking
up the points in the blend group in the source coordset and multiplying them by their weights
and summing them together to produce the new point for the output coordset. Classes such as
``ClipField`` use ``CoordsetBlender`` to make new coordsets that contain points that were a
combination of multiple points in the input coordset.

.. literalinclude:: ../../extraction/TableBasedExtractor.hpp
   :start-after: _bump_utilities_coordsetblender_begin
   :end-before: _bump_utilities_coordsetblender_end
   :language: C++

################
CoordsetSlicer
################

The ``axom::bump::CoordsetSlicer`` class takes ``SliceData`` and makes a
new explicit coordset where each point corresponds to a single index from the node indices
stored in SliceData. This class can be used to select a subset of a coordset, reorder nodes
in a coordset, or repeat nodes in a coordset.

.. literalinclude:: ../../ExtractZones.hpp
   :start-after: _bump_utilities_coordsetslicer_begin
   :end-before: _bump_utilities_coordsetslicer_end
   :language: C++

##########
CutField
##########

The ``axom::bump::extraction::CutField`` class intersects all the zones in the input Blueprint
mesh with an implicit surface where the selected input field equals zero. This is an isosurface
algorithm, though it could be used with other intersection policies. The default algorithm
produces a new Blueprint mesh containing topologically 2D or 1D zone fragments, based on the 
input mesh's dimension. The ``CutField`` class has multiple template arguments to
select the execution space, the type of topology view, the type of coordset view, and the
type of intersector used to determine intersections. The default intersection uses an isosurface-
based intersection method, though other intersectors could be created to perform plane
or sphere intersections.

.. literalinclude:: ../../tests/bump_cutfield.cpp
   :start-after: _bump_utilities_cutfield_begin
   :end-before: _bump_utilities_cutfield_end
   :language: C++

##################
ExtractZones
##################

The ``axom::bump::ExtractZones`` class takes a list of selected zone ids and extracts
a new mesh from a source mesh that includes only the selected zones. There is a derived class
``ExtractZonesAndMatset`` that also extracts a matset, if present.

.. literalinclude:: ../../ExtractZones.hpp
   :start-after: _bump_utilities_extractzones_begin
   :end-before: _bump_utilities_extractzones_end
   :language: C++

##################
ExtrudeMesh
##################

The ``axom::bump::ExtrudeMesh`` class extrudes a 2D Blueprint mesh composed
of triangles and quad shapes *(polygons are not yet supported)* and produces 3D zones repeated some
number of times in the Z direction. Fields and matsets are also extruded.

.. literalinclude:: ../../tests/bump_topology_mapper.cpp
   :start-after: _bump_utilities_extrudemesh_begin
   :end-before: _bump_utilities_extrudemesh_end
   :language: C++

#############
FieldBlender
#############

The ``axom::bump::FieldBlender`` class is similar to the ``CoordsetBlender``
class, except that it operates on a field instead of coordsets. The class is used to create a
new field that includes values derived from multiple weighted source values.

############
FieldSlicer
############

The ``axom::bump::FieldSlicer`` class selects specific indices from a
field and makes a new field.

.. literalinclude:: ../../tests/bump_slicers.cpp
   :start-after: _bump_utilities_fieldslicer_begin
   :end-before: _bump_utilities_fieldslicer_end
   :language: C++

##############
MakePointMesh
##############

The ``axom::bump::MakePointMesh`` class generates a point at the center
of each zone *(or selected set of zones)* in an input topology and generates a new unstructured
topology consisting of points located at those zone centers.

.. literalinclude:: ../../../mir/ElviraAlgorithm.hpp
   :start-after: _bump_utilities_makepointmesh_begin
   :end-before: _bump_utilities_makepointmesh_end
   :language: C++

#######################
MakePolyhedralTopology
#######################

The ``axom::bump::MakePolyhedralTopology`` class transforms an
input topology from its native form to an unstructured polyhedral topology. The output
topology uses the same coordset as the input topology. The faces produced from each zone
in the source topology will not be unique. The ``MergePolyhedralFaces`` class can be used
to merge polyhedral faces so they are unique.

.. literalinclude:: ../../tests/bump_make_polyhedral_topology.cpp
   :start-after: _bump_utilities_makepolyhedraltopology_begin
   :end-before: _bump_utilities_makepolyhedraltopology_end
   :language: C++

##################
MakeUnstructured
##################

The ``axom::bump::MakeUnstructured`` class takes a structured topology
and creates a new unstructured topology. This class does not need views to wrap the input
structured topology.

.. literalinclude:: ../../tests/bump_mesh_operations.cpp
   :start-after: _bump_utilities_makeunstructured_begin
   :end-before: _bump_utilities_makeunstructured_end
   :language: C++

##################
MakeZoneCenters
##################

The ``axom::bump::MakeZoneCenters`` class takes an input Blueprint
topology and produces a new element-associated Blueprint vector field that contains the zone
centers. The number of components in the vector will match the number of components for the
topology's coordset. The zone center is computed as the average of the node coordinates used
in the zone. Likewise, the type *(e.g. float, double)* used to compute and represent the zone
centers will match the type of the values that define the coordset.

.. literalinclude:: ../../../mir/ElviraAlgorithm.hpp
   :start-after: _bump_utilities_makezonecenters_begin
   :end-before: _bump_utilities_makezonecenters_end
   :language: C++

##################
MakeZoneVolumes
##################

The ``axom::bump::MakeZoneVolumes`` class takes an input Blueprint
topology and produces a new element-associated Blueprint vector field that contains the zone
volumes for 3D, or areas for 2D.

##################
MatsetSlicer
##################

The ``axom::bump::MatsetSlicer`` class is similar to the ``FieldSlicer``
class except it slices matsets instead of fields. The same ``SliceData`` can be passed to
MatsetSlicer to pull out and assemble a new matset data for a specific list of zones.

.. literalinclude:: ../../ExtractZones.hpp
   :start-after: _bump_utilities_matsetslicer_begin
   :end-before: _bump_utilities_matsetslicer_end
   :language: C++

####################
MergeCoordsetPoints
####################

The ``axom::bump::MergeCoordsetPoints`` class merges duplicate
coordinates in an input coordset, within a given tolerance.
The tolerance is passed via an options node with a key value called "tolerance". Points
are merged by first rounding off extra precision in a temporary point copy that is used to
make a hashed name for the point. Points within the tolerance get the same hashed
name and points are made unique using the ``Unique`` class.

The selected point for each unique name is copied into the coordset, overwriting the old
coordset values. The number of points in the coordset will change if any merging is done.
This means that any topologies that reference the coordset will need to be updated.
The class passes out a selectedIds array containing indices from the original coordset
that are used in the new coordset. This information can be used with ``FieldSlicer`` to
slice/update vertex fields. An ``old2new`` array is also returned, which is a
map to convert old node indices to new indices in the updated coordset. This map can be
used to update connectivity node numbers.

.. literalinclude:: ../../../mir/detail/elvira_detail.hpp
   :start-after: _bump_utilities_mergecoordsetpoints_begin
   :end-before: _bump_utilities_mergecoordsetpoints_end
   :language: C++

##################
MergeMeshes
##################

The ``axom::bump::MergeMeshes`` class merges data for coordsets,
topology, and fields from multiple input meshes into a new combined mesh. The class also
supports renaming nodes using a map that converts a local mesh's node ids to the final
output node numbering, enabling meshes to be merged such that some nodes get combined.
A derived class can also merge matsets.

.. literalinclude:: ../../tests/bump_mergemeshes.cpp
   :start-after: _bump_utilities_mergemeshes_begin
   :end-before: _bump_utilities_mergemeshes_end
   :language: C++

#####################
MergePolyhedralFaces
#####################

The ``axom::bump::MergePolyhedralFaces`` class takes an input
Blueprint topology, which may have duplicated faces, and makes the face definitions
in the subelements unique and rewrites the subelement and element connectivity. For
faces to be merged successfully, the faces must reference the same coordinate indices
in the coordset. The ``MergePolyhedralFaces`` class modifies the Conduit node that
contains the input polyhedral topology.

.. literalinclude:: ../../tests/bump_make_polyhedral_topology.cpp
   :start-after: _bump_utilities_makepolyhedraltopology_begin
   :end-before: _bump_utilities_makepolyhedraltopology_end
   :language: C++

###########################
NodeToZoneRelationBuilder
###########################

The ``axom::bump::NodeToZoneRelationBuilder`` class creates a Blueprint
O2M (one to many) relation that relates node numbers to the zones that contain them. This mapping
is akin to inverting the normal mesh connectivity which is a map of zones to node ids. The O2M
relation is useful for recentering data from the zones to the nodes.

.. literalinclude:: ../../tests/bump_node_to_zone_relation.cpp
   :start-after: _bump_utilities_n2zrel_begin
   :end-before: _bump_utilities_n2zrel_end
   :language: C++

###########
PlaneSlice
###########

The ``axom::bump::extraction::PlaneSlice`` class slices input Blueprint geometry using
a slice plane and produces a new Blueprint output mesh. This algorithm is a close cousin
to ``axom::bump::extraction::CutField``, except that it uses a plane intersector that
accepts "origin" and "normal" parameters to specify the slice plane.

.. literalinclude:: ../../tests/bump_planeslice.cpp
   :start-after: _bump_utilities_planeslice_begin
   :end-before: _bump_utilities_planeslice_end
   :language: C++

###############
PrimalAdaptor
############### 

The ``axom::bump::PrimalAdaptor`` class takes a topology view and a
coordset view and makes it possible to retrieve a zone as a shape from Axom's Primal
component. For example, the ``PrimalAdaptor`` class can wrap a topology view that contains 2D
shapes such as triangles, quads, polygons and allow them to be accessed as an
``axom::primal::Polygon``. For 3D, Primal shapes are returned for meshes that contain
tetrahedra, hexahedra, or polyhedra. For unstructured meshes that contain pyramids or
wedges, or mixed shapes, a ``VariableShape`` is returned that allows those shapes to be
represented using one or more primal shapes.

When the class is instantiated with
``axom::bump::views::UnstructuredTopologyPolyhedralView`` as its topology view, the ``getShape()``
method will normally return ``axom::primal::Polyhedron``. Converting between Blueprint polyhedron
zones and Axom ``Polyhedron`` objects is sometimes overkill so the ``PrimalAdaptor`` class can
also expose polyhedra as a special ``PolyhedralFaces`` representation that represents the
polyhedron as a collection of ``axom::primal::Plane`` objects. This mode is selected by
instantiating ``PrimalAdaptor`` with the ``makeFaces`` template parameter set to true.

.. literalinclude:: ../../MakeZoneVolumes.hpp
   :start-after: _bump_utilities_makezonevolumes_begin
   :end-before: _bump_utilities_makezonevolumes_end
   :language: C++


###############
RecenterField
############### 

The ``axom::bump::RecenterField`` class uses an O2M relation to average
field data from multiple values to an averaged value. In Axom, this is used to convert a field
associated with the elements to a new field associated with the nodes.

.. literalinclude:: ../../tests/bump_mesh_operations.cpp
   :start-after: _bump_utilities_recenterfield_begin
   :end-before: _bump_utilities_recenterfield_end
   :language: C++

##################
SelectedZones
##################

The ``axom::bump::SelectedZones`` class creates an array view that
represents selected zone ids. The zone ids are obtained either from a Conduit options
node containing a *"selectedZones"* array, if the array is present. If the "selectedZones"
array is not present, the class makes an array of zone ids that selects all zones in the
associated topology.

.. literalinclude:: ../../../mir/ElviraAlgorithm.hpp
   :start-after: _bump_utilities_selectedzones_begin
   :end-before: _bump_utilities_selectedzones_end
   :language: C++


##################
TopologyMapper
##################

The ``axom::bump::TopologyMapper`` class intersects a source mesh
with a target mesh and maps materials from the source mesh onto a new matset on the 
target mesh. The source mesh must contain a "clean" matset, which is a matset where there
are no mixed-material zones. The matset identifies the unique material for each zone in
the source mesh. The source mesh could be the output of one of the BUMP algorithms.

The source and target meshes should overlap spatially. The zones in the source mesh are
intersected with the zones in the target mesh and their overlaps are determined and are
used to build a new matset on the target mesh. Each zone in the target mesh may recieve
contributions from multiple zones and materials in the source mesh.

.. literalinclude:: ../../tests/bump_topology_mapper.cpp
   :start-after: _bump_utilities_topologymapper_begin
   :end-before: _bump_utilities_topologymapper_end
   :language: C++


##################
Unique
##################

The ``axom::bump::Unique`` class can take an unsorted list of values and produce a
sorted list of unique outputs, along with a list of offsets into the original values to identify
one representative value in the original list for each unique value. This class is used to help
merge points.

.. literalinclude:: ../../tests/bump_clipfield.cpp
   :start-after: _bump_utilities_unique_begin
   :end-before: _bump_utilities_unique_end
   :language: C++

##################
VariableShape
##################

The ``axom::bump::VariableShape`` class behaves like a primal shape
but it can represent various 3D shapes, some not present in primal.

##################
ZoneListBuilder
##################

The ``axom::bump::ZoneListBuilder`` class takes a matset view and a list
of selected zone ids and makes two output lists of zone ids that correspond to clean zones and
mixed zones (more than 1 material in the zone). There are also methods that take into consideration
how zones are connected through their nodes so algorithms that operate on node-centered volume
fractions can operate on adjacent zones that may not be mixed but must participate in BUMP.

.. literalinclude:: ../../../mir/EquiZAlgorithm.hpp
   :start-after: _bump_utilities_zlb_begin
   :end-before: _bump_utilities_zlb_end
   :language: C++
