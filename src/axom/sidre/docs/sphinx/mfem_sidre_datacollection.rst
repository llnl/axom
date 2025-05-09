.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _mfem-sidre-datacollection:

*********************
Using Sidre with MFEM
*********************

.. note::
   The functionality described in this page is only available if Axom is 
   configured with MFEM and if the CMake variable 
   ``AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION`` is set to ``ON``

The ``MFEMSidreDataCollection`` class implements `MFEM <https://mfem.org>`_'s 
``DataCollection`` interface for recording simulation data. Specifically, it 
knows about fields (``mfem::GridFunction`` objects) and the mesh on which 
these fields are defined.

The ``MFEMSidreDataCollection`` class internally organizes its data according 
to the Mesh Blueprint, a hierarchical schema for describing mesh data.

See the :ref:`Conduit page <sidre-conduit>` for more information on the Mesh 
Blueprint.

In this page, we first discuss how MFEM objects can 
be associated with the ``MFEMSidreDataCollection``. We then explain the process
and options available when `Saving Data to a File`_. The workflow for reading 
saved data back in is discussed in `Restarting a Simulation`_.

Associating MFEM Objects with an MFEMSidreDataCollection
--------------------------------------------------------

We begin to describe the data in an ``MFEMSidreDataCollection`` object by 
*registering* an MFEM mesh with an instance of the class when we create the
data collection:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_vis.cpp
   :start-after: _sidredc_vis_construct_start
   :end-before: _sidredc_vis_construct_end
   :language: C++

The MFEM mesh can also be registered after construction:

.. code-block:: C++

  dc.SetMesh(/* mfem::Mesh* */ mesh);

.. note::
   There is a 1-1 relationship between ``MFEMSidreDataCollection`` objects and 
   MFEM meshes. That is, multiple meshes cannot be associated with a single 
   ``MFEMSidreDataCollection``.

After a mesh has been registered, fields can be association with the 
``MFEMSidreDataCollection`` object:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_vis.cpp
   :start-after: _sidredc_vis_field_start
   :end-before: _sidredc_vis_field_end
   :language: C++

Special kinds of fields (material-dependent fields, material set volume 
fraction fields, and species set values fields) require additional setup before 
they can be registered - see `Mixed-material Fields`_ for more informantion.

Saving Data to a File
---------------------

The data in an instance of the ``MFEMSidreDataCollection`` class can be saved 
to a file using a variety of protocols. These files can be visualized with 
tools like `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit>`_
or used to restart a simulation by loading them back into an instance of the 
``MFEMSidreDataCollection`` class (see `Restarting a Simulation`_).

Current options for output protocols/formats include:

   - ``sidre_hdf5``
   - ``sidre_conduit_json``
   - ``sidre_json``
   - ``conduit_hdf5``

.. Note::
   The ``AXOM_USE_HDF5`` CMake build option must be enabled to use
   HDF5-based formats.

Before saving, simulation state metadata should be updated. Currently, this 
metadata consists of:

- ``cycle``, the current iteration number for the simulation
- ``time``, the current simulation time
- ``time_step``, the current simulation time step (sometimes called ``dt``)

Each of these variables has a corresponding *setter* function:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_vis.cpp
   :start-after: _sidredc_vis_state_start
   :end-before: _sidredc_vis_state_end
   :language: C++

.. note::
   There are also corresponding accessors for these state quantities (``GetCycle``, ``GetTime``, ``GetTimeStep``).

Once state information has been updated, the complete simulation state can 
be written to a file:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_vis.cpp
   :start-after: _sidredc_vis_save_start
   :end-before: _sidredc_vis_save_end
   :language: C++

.. note::
   By default, saved files will be written to the current working directory. 
   To write to/read from a different directory, use ``SetPrefixPath`` to 
   change the directory location of where the ``MFEMSidreDataCollection`` 
   object's data will be written.

See the ``sidre_mfem_datacollection_vis`` example for a more thorough example 
of the above functionality.

Restarting a Simulation
-----------------------

Experimental support for complete reconstruction of a simulation mesh, fields, 
and qfields is also provided by the ``MFEMSidreDataCollection`` class. That is,
when an output file is read in using the ``MFEMSidreDataCollection::Load()``
method, the data read in will be used to reconstruct MFEM objects than can be 
accessed with the ``GetField``, ``GetQField``, and ``GetMesh`` methods.

.. warning::
   Currently, an ``MFEMSidreDataCollection`` object must own the mesh and field
   data in order to completely reconstruct simulation state. Mesh data 
   ownership can be configured with the ``owns_mesh_data`` constructor option 
   (should be set to ``true``), and field data ownership requires that each
   ``GridFunction`` object be unallocated when passed to the 
   ``RegisterField()`` method, which performs the allocation within Sidre-owned
   memory. After registration, a ``GridFunction`` object can be used normally.
   The same conditions apply for ``QuadratureFunction`` objects.

A complete demonstration of functionality is provided in the 
``sidre_mfem_datacollection_restart`` example, which is a stripped-down
example of how a simulation code might utilize the automatic reconstruction logic when loading in a datastore.

.. Note::
  The mesh/field reconstruction logic requires that the saved file was created 
  with an ``MFEMSidreDataCollection`` object. In `Mesh Blueprint <http://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html>`_ terms, the following 
  constraints are imposed on the structure of the data:

  * There must be a coordinate set named ``coords``
  * There must be a topology named ``mesh`` with corresponding attributes stored in a field named ``mesh_material_attribute``
  * There must be a topology named ``boundary`` with corresponding attributes stored in a field named ``boundary_material_attribute``

Mixed-material Fields
---------------------

The Mesh Blueprint provides support for mixed-material simulations through 
its *material set* construct. Material metadata (stored in 
``mfem::GridFunction`` objects) can be registered with the 
``MFEMSidreDataCollection`` like any other field (with the 
``RegisterField()`` method), given that some additional setup is performed 
and the field names match a specific naming convention.

Currently, there are three kinds of *special* fields that can be associated 
with mixed-material metadata:

- Volume fraction fields in a material set (the per-element volume fraction of a given material)
- Material-dependent fields (a different set of values for each material)
- Species set fields (a different set of values for each material and for each dimension)

`Material sets <https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#material-sets>`_ are
defined by associating a volume fraction field name with a material set name:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_materials.cpp
   :start-after: _sidredc_material_matset_start
   :end-before: _sidredc_material_matset_end
   :language: C++

Once this association has been made, corresponding fields can be registered.
Their names must follow the format ``<associated volume fraction field name>_<material id>``, for example:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_materials.cpp
   :start-after: _sidredc_material_matset_register_start
   :end-before: _sidredc_material_matset_register_end
   :language: C++

`Material-dependent fields <https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#fields>`_ are
defined by associating a field name with a material set name:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_materials.cpp
   :start-after: _sidredc_material_depfield_start
   :end-before: _sidredc_material_depfield_end
   :language: C++

Material-specific values can be registered with a name in the format ``<associated dependent field name>_<material id>``:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_materials.cpp
   :start-after: _sidredc_material_depfield_register_start
   :end-before: _sidredc_material_depfield_register_end
   :language: C++

`Species sets <https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#species-sets>`_ are
defined by associating a field name with a species set name and corresponding material set name:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_materials.cpp
   :start-after: _sidredc_material_specset_start
   :end-before: _sidredc_material_specset_end
   :language: C++

Species set values can be registered with the name in the format ``<associated field name>_<material id>_<component>``:

.. literalinclude:: ../../examples/sidre_mfem_datacollection_materials.cpp
   :start-after: _sidredc_material_specset_register_start
   :end-before: _sidredc_material_specset_register_end
   :language: C++
