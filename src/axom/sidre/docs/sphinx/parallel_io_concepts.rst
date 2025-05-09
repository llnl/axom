.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _spio-core-concepts:

******************************************************
Parallel File I/O
******************************************************


The Sidre ``IOManager`` class provides an interface to manage parallel I/O of 
hierarchical data managed by Sidre. It enables the writing of data from 
parallel runs and can be used for the purposes of restart or visualization.

Introduction
-------------

The ``IOManager`` class provides parallel I/O services to Sidre. It relies on 
the fact that Sidre ``Group`` and ``View`` objects are capable of saving and 
loading themselves. The I/O operations provided by those two classes are 
inherently serial. So the ``IOManager`` class coordinates I/O operations 
of multiple Sidre objects that exist across the MPI ranks of a parallel run.

* The internal details of the I/O of individual Sidre objects are opaque to
  the ``IOManager`` class, which needs only to make calls to I/O methods in 
  ``Group`` and ``View`` class methods.
* Sidre data can be written from M ranks to N files (M >= N), and the files can
  be read to restart a run on M ranks.
* When saving output, a root file is created that contains bookkeeping
  information that is used to coordinate a subsequent restart read.
* The calling code can also add extra data to the root file to provide 
  metadata that gives necessary instructions to visualization tools.

Parallel I/O using the IOManager class
---------------------------------------

A Sidre ``IOManager`` object is created with an MPI communicator and provides 
various overloads of ``write()`` and ``read()`` methods.  These methods save 
a ``Group`` object in parallel to a set of files and read a ``Group`` object
from existing files. The I/O manager can optionally use the 
`SCR <https://github.com/LLNL/scr>`_ library for 
scalable I/O management (such as using burst buffers when available).

In typical usage, a run that calls ``read()`` on a certain set of files
should be executed on the same number of MPI ranks as the run that created
those files with a ``write()`` call.  However, if using the ``sidre_hdf5``
protocol, there are some usage patterns that do not have this limitation.

A ``read()`` call using the ``sidre_hdf5`` protocol will work when called from 
a greater number of MPI ranks. If ``write()`` was performed on N ranks and 
``read()`` is called while running on M ranks (M > N), then data will be read 
into ranks 0 to N-1, and all ranks higher than N-1 will receive no data.

If the ``read()`` method is called using the ``sidre_hdf5`` protocol to read 
data that was created on a larger number of ranks, this will work only in 
the case that the data was written in a file-per-processor mode (M ranks to 
M files). In this case, the data in the group being filled with file input 
will look a bit different than in other usage patterns, since a group on one 
rank will end up with data from multiple ranks. An integer scalar view named 
``reduced_input_ranks`` will be added to the group with the value being the 
number of ranks that wrote the files.  The data from each output rank will 
be read into subgroups located at ``rank_{%07d}/sidre_input`` in the input 
group's data hierarchy.

.. warning::
   If the ``read()`` is called to read data that was created on a larger
   number of ranks than the current run with files produced in M-to-N
   mode (M > N), an error will occur. Support for this type of usage is
   intended to be added in the future.

In the following example, an ``IOManager`` object is created and used to write 
the contents of the "root" group in parallel.

First include needed headers.

.. literalinclude:: ../../tests/spio/spio_parallel.hpp
   :start-after: _parallel_io_headers_start
   :end-before: _parallel_io_headers_end
   :language: C++

Then use IOManager to save in parallel.

.. literalinclude:: ../../tests/spio/spio_parallel.hpp
   :start-after: _parallel_io_save_start
   :end-before: _parallel_io_save_end
   :language: C++

Loading data in parallel is easy:

.. literalinclude:: ../../tests/spio/spio_parallel.hpp
   :start-after: _parallel_io_load_start
   :end-before: _parallel_io_load_end
   :language: C++

IOManager class use
-------------------

An ``IOManager`` object is constructed with an MPI communicator and does I/O 
operations on all ranks associated with that communicator.

The core functionality of the ``IOManager`` class is contained in its 
``write()`` and ``read()`` methods.

.. code-block:: cpp

  void write(sidre::DataGroup * group,
             int num_files,
             const std::string& file_string,
             const std::string& protocol);

The ``write()`` method is called on each rank and passed a group pointer
for its local data. The calling code specifies the number of output files,
and the I/O manager organizes the output so that each file receives data from a
roughly equal number of ranks.  The files containing the data from the group
will have names of the format ``file_string/file_string_*******.suffix``, with a
7-digit integer value identifying the files from 0 to num_files-1, and the
suffix indicating the file format according to the protocol argument.
Additionally the ``write()`` method will produce a root file with the name 
``file_string.root`` that holds bookkeeping data about other files and can 
also receive extra user-specified data.

.. code-block:: cpp

  void read(sidre::DataGroup * group,
            const std::string& root_file);


The ``read()`` methods is called in parallel with the root file as an argument.
It must be called on a run with the same processor count as the run that 
called the ``write()`` method. The first argument is a pointer to a group that 
contains no child groups and owns no views. The second argument is the name of
the root file; the information in the root file is used to identify the files 
that each processor will read to load data into the argument group.

The ``write()`` and ``read()`` methods above are sufficient to do a restart 
save/load when the data in the group hierarchy is completely owned by the 
Sidre data structures. If Sidre is used to manage data that is externally 
allocated (i.e., the hiearchy contains *external* views), the loading
procedure requires additional steps to restore data in the same
externally-allocated state.

First the ``read()`` method is called, and the full hierarchy structure of the
group is loaded, but no data is allocated for views identified as external. 
Then, the calling code can examine the group and allocate data for the external
views, as needed.  The ``View`` class method ``setExternalDataPtr()`` is used
to associate the pointer with the view.  Once this is done, the I/O manager
method ``loadExternalData()`` can be used to load the data from the file into 
the user-allocated arrays.

Below is a code example for loading external data. We assume that this code
somehow has knowledge that the root group contains a single external view at
the location "fields/external_array" describing an array of doubles. See the
`Sidre API Documentation <../../../../doxygen/html/sidretop.html>`_ for details
about ``Group`` and ``View`` class methods to query the Sidre data structures 
for this type of information when the code does not have a priori knowledge.

.. code-block:: cpp

  // Construct a DataStore with an empty root group.
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  // Read from file into the root group.  The full Sidre hierarchy is built,
  // but the external view is created without allocating a data buffer.
  IOManager reader(MPI_COMM_WORLD);
  reader.read(root, "checkpoint.root");

  // Get a pointer to the external view. 
  DataView * external_view = root->getView("fields/external_array");

  // Allocate storage for the array and associate it with the view.
  double * external_array = new double[external_view->getNumElements()];
  external_view->setExternalDataPtr(external_array);

  // Load the data values from file into the external view.
  reader.loadExternalData(root, "checkpoint.root");


User-specified data in the root file
------------------------------------

The root file is automatically created to provide the ``IOManager`` object 
with bookkeeping information that is used when reading data, but it can also
be used to store additional data that may be useful to the calling code or
is needed to allow other tools to interact with the data in the output files,
such as for visualization.  For example, a `Conduit Mesh Blueprint <https://llnl-conduit.readthedocs.io/en/latest/blueprint.html>`_ index can be 
:ref:`stored in a DataGroup <sidre-conduit>` written to the root file to 
provide metadata about the mesh layout and data fields that can be visualized 
from the output files.

.. code-block:: cpp

  void writeGroupToRootFile(sidre::DataGroup * group,
                            const std::string& file_name);

  void writeGroupToRootFileAtPath(sidre::DataGroup * group,
                                  const std::string& file_name,
                                  const std::string& group_path);

  void writeViewToRootFileAtPath(sidre::DataView * view,
                                 const std::string& file_name,
                                 const std::string& group_path);

The above methods are used to write this extra data to the root file. The
first simply writes data from the given group to the top of the root file,
while the latter two methods write their Sidre objects to a path that must
already exist in the root file.
