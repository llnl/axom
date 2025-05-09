.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _sidre-conduit:

******************************************************
Sidre Interaction with Conduit
******************************************************

Internally, Sidre uses the in-memory data description capabilities of the
`Conduit <https://github.com/LLNL/conduit>`_ library. Sidre also leverages 
Conduit to facilitate data exchange, shown here applied to visualization. 
The following discussion gives a basic overview of Sidre capabilities
when combined with Conduit.

Mesh Blueprint
--------------

The `Mesh Blueprint <http://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html>`_
is a data exchange protocol supported by Conduit, consisting of a
properly-structured Sidre datastore hierarchy saved as an HDF5 or JSON file 
and a Conduit index file. The Blueprint can accomodate structured or 
unstructured meshes, with node- or element-centered fields. The following 
example shows how to create a Blueprint-conforming datastore hierarchy 
containing two unstructured adjacent hexahedrons with one node-centered field 
and one element-centered field. In the diagram, nodes are labeled in black, 
the node-centered field values are in blue, and the element-centered field 
values are in green.

.. image:: figs/tiny_mesh.png
   :width: 300px

A simulation organizes Sidre data as the code design dictates.
Here is a simple example data hierarchy.

.. image:: figs/ds.png
   :width: 650px

Here is the code to create that example data hierarchy.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _tiny_create_start
   :end-before: _tiny_create_end
   :language: C++

To use the Conduit Mesh Blueprint, make a group hierarchy ``tinymesh`` 
conforming to the Mesh Blueprint protocol. The structure of the group 
hierarchy is shown below (summarizing
the `Mesh Blueprint <http://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html>`_
documentation).

First build top-level groups required by the Mesh Blueprint.

.. image:: figs/cds.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_toplevel_start
   :end-before: _blueprint_restructure_toplevel_end
   :language: C++

Add the node coordinates. The views under ``tinymesh`` will point to the
same buffers that were created for the views under ``nodes`` so that 
``tinymesh`` can use the data without new data allocations or data copying.

.. image:: figs/cdscoords.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_coords_start
   :end-before: _blueprint_restructure_coords_end
   :language: C++

Arrange the nodes into elements. Each simulation has its own knowledge of
its mesh topology. This tiny example didn't previously encode the topology,
so we must explicitly specify it.

.. image:: figs/cdstopo.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_topo_start
   :end-before: _blueprint_restructure_topo_end
   :language: C++

Link the fields into ``tinymesh``. As with the node positions, the views
point to the existing buffers containing the field data.

.. image:: figs/cdsfields.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_field_start
   :end-before: _blueprint_restructure_field_end
   :language: C++

Conduit includes a ``verify`` method to test if the structure of the 
``tinymesh`` conforms to the Mesh Blueprint. This is valuable for writing and 
debugging data adapters. Once the datastore hierarchy is properly structured, 
save it, then use Conduit to save the index file (ending with `.root`). This 
toy data set is small enough that we can choose to save it in a JSON format.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_save_start
   :end-before: _blueprint_restructure_save_end
   :language: C++

The code listed above produces the files `tinymesh.json` and `tinymesh.root`.
Any code that uses Mesh Blueprint can open and use this pair of files.

The DataStore also contains a method that can automatically generate the
Blueprint index within a Sidre group rather than calling directly into
Conduit. Set up a mesh similarly to the example above.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_generate_toplevel_start
   :end-before: _blueprint_generate_toplevel_end
   :language: C++

Then call the ``DataStore::generateBlueprintIndex()`` method to generate the 
index within a group in the datastore. Additional data needed in the root file
can be added and saved using Sidre I/O calls.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_generate_save_start
   :end-before: _blueprint_generate_save_end
   :language: C++

Additionally, the Sidre Parallel I/O (SPIO) class ``IOManager`` provides a
method that both generates a Blueprint index and adds it to a root file.
Using the same mesh data from the last example, first write out all of the
parallel data using the ``IOManager::write()`` method. This will output to 
files all of the data for all domains, and will also create a basic root file.
Then the ``IOManager::writeBlueprintIndexToRootFile()`` methods can be called 
to generate the Blueprint index and add it to the root file. This is currently 
only implemented to work with the ``sidre_hdf5`` I/O protocol.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_generate_spio_start
   :end-before: _blueprint_generate_spio_end
   :language: C++

Data Visualization
------------------

The `VisIt <http://visitusers.org/>`_ tool can read in a Conduit Mesh 
Blueprint conforming file, interpret the index file, and sensibly display the 
data contained in the data file. Starting from version 2.13.1, VisIt can open 
a `.root` file just like any other data file. VisIt produced the following 
image from the Mesh Blueprint file saved above.

.. image:: figs/tiny_mesh_rendered.png
   :width: 600px

Conduit is also a foundational building block for the 
`Ascent <http://ascent.readthedocs.io/>`_ project, which provides a powerful
*in situ* data analytics and visualization facility (without copying memory) to 
distributed-memory simulation codes.
