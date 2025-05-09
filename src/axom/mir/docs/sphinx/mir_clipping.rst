.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

*************
Clipping
*************

The MIR component provides a clipping algorithm that can perform isosurface-based
clipping and return volumetric output for zones and partial zones that are "inside"
or "outside" the clip boundary. The clipping algorithm is implemented in the
``axom::mir::clipping::ClipField`` class. The class can be instantiated with several
template arguments that govern where it will execute, which coordset and topology
types it supports, and how it performs intersection. The input to the algorithm is
a Blueprint mesh. When instantiated with coordset and topology views appropriate
for the input data, the algorithm can operate on a wide variety of mesh types. This
includes 2D/3D structured and unstructured topologies that can be represented using
finite elements.

By default, the algorithm will clip using a field but other intersection routines
can be substituted via a template argument to facilitate creation of clipping using
planes, spheres, surfaces of revolution, etc. The Equi-Z algorithm uses ClipField
with an intersector that examines material volume fractions to determine the clipped geometry.

#######
Inputs
#######

Like the MIR algorithms, the clipping algorithm is designed to accept a Conduit node
containing various options that influence how the algorithm operates. The clipping
algorithm copies the options node to the memory space where it will be used.

+---------------------------------+------------------------------------------------------+
| Option                          | Description                                          |
+=================================+======================================================+
| ``clipField: name``             | A required string argument that specifies the name   |
|                                 | of the field that is used for clipping. At present,  |
|                                 | the field must be a vertex-associated field.         |
+---------------------------------+------------------------------------------------------+
| ``clipValue: value``            | An optional numeric argument that specifies the      |
|                                 | value in the field at which the clip boundary is     |
|                                 | defined. The default is 0.                           |
+---------------------------------+------------------------------------------------------+
| ``colorField: name``            | If inside=1 and outside=1 then a color field is      |
|                                 | generated so it is possible to tell apart regions of |
|                                 | the clip output that were inside or outside the clip |
|                                 | boundary. This field permits the user to change the  |
|                                 | name of the color field, which is called "color" by  |
|                                 | default.                                             |
+---------------------------------+------------------------------------------------------+
| ``coordsetName: name``          | The name of the new coordset in the output mesh. If  |
|                                 | it is not provided, the output coordset will have the|
|                                 | same name as the input coordset.                     |
+---------------------------------+------------------------------------------------------+
|``fields:``                      | The fields node lets the caller provide a list of    |
|                                 | field names that will be processed and added to the  |
|                                 | output mesh. The form is *currentName:newName*. If   |
|                                 | the *fields* node is not given, the algorithm will   |
|                                 | process all input fields. If the fields node is empty|
|                                 | then no fields will be processed.                    |
+---------------------------------+------------------------------------------------------+
| ``inside: number``              | Indicates to the clipping algorithm that it should   |
|                                 | preserve zone fragments that were "inside" the clip  |
|                                 | boundary. Set to 1 to enable, 0 to disable. The      |
|                                 | algorithm will generate these fragments by default.  |
+---------------------------------+------------------------------------------------------+
| ``originalElementsField: name`` | The name of the field in which to store the original |
|                                 | elements map. The default is "originalElements".     |
+---------------------------------+------------------------------------------------------+
| ``outside: number``             | Indicates to the clipping algorithm that it should   |
|                                 | preserve zone fragments "outside" the clip boundary. |
|                                 | Set to 1 to enable, 0 to disable. These fragments are|
|                                 | not on by default.                                   |
+---------------------------------+------------------------------------------------------+
| ``selectedZones: [zone list]``  | An optional argument that provides a list of zone ids|
|                                 | on which to operate. The output mesh will only have  |
|                                 | contributions from zone numbers in this list, if it  |
|                                 | is given.                                            |
+---------------------------------+------------------------------------------------------+
| ``topologyName: name``          | The name of the new topology in the output mesh. If  |
|                                 | it is not provided, the output topology will have the|
|                                 | same name as the input topology.                     |
+---------------------------------+------------------------------------------------------+

##########
ClipField
##########

To use the ``ClipField`` class, one must have Blueprint data with at least one vertex-associated
field. Views for the coordset and topology are created and their types are used to instantiate
a ``ClipField`` object. The ``ClipField`` constructor takes a Conduit node for the input Blueprint mesh, a Conduit
node that contains the options, and a 3rd output Conduit node that will contain the clipped
mesh and fields. The input mesh node needs to contain data arrays for coordinates, mesh
topology, and fields. These data must exist in the memory space of the targeted device.
Other Conduit nodes that contain strings or single numbers that can fit within a Conduit
node are safe remaining in host memory. If the mesh is not in the desired memory space, it
can be moved using ``axom::mir::utilities::blueprint::copy()``.

.. code-block:: cpp

    #include "axom/mir.hpp"

    // Set up views for the mesh in deviceRoot node.
    auto coordsetView = axom::mir::views::make_rectilinear_coordset<float, 3>::view(deviceRoot["coordsets/coords"]);
    auto topologyView = axom::mir::views::make_rectilinear<3>::view(deviceRoot["topologies/Mesh"]);

    // Make a clipper.
    using CoordsetView = decltype(coordsetView);
    using TopologyView = decltype(topologyView);
    using Clip = axom::mir::clipping::ClipField<axom::SEQ_EXEC, TopologyView, CoordsetView>;
    Clip clipper(topologyView, coordsetView);

    // Run the clip algorithm
    conduit::Node options;
    options["clipField"] = "data";
    options["clipValue"] = 3.5;
    options["outside"] = 1;
    options["inside"] = 0;
    clipper.execute(deviceRoot, options, clipOutput);


.. figure:: figures/clipfield.png
   :figwidth: 800px

   Diagram showing original mesh colored by clipping field (left), original mesh colored by a radial field (middle), and the clipped mesh colored by the radial field (right).


^^^^^^^^^^^^^
Intersectors
^^^^^^^^^^^^^

An intersector is a policy class that is passed as a template argument to ``ClipField``. The
intersector determines how the ``ClipField`` algorithm will generate intersection cases, for
each zone in the mesh. The ``ClipField`` algorithm default intersector uses a field to determine clip
cases, resulting in isosurface behavior for the geometry intersections. Alternative intersectors
can be provided to achieve other types of intersections.
