.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _dataconcepts-label:

===========================
Data vs. Metadata Concepts
===========================

This section includes examples that illustrate some common Sidre 
usage patterns involving buffers, groups, and views. In particular, it shows
how data associated with a view is related to a buffer, how the lifetimes
of a view, its associated data, and the buffer owning the data are related,
and how the description of data in a view (metadata) is independent of the
the allocation state of the data in a view and a corresponding buffer. 

.. note:: Object and data creation and destruction is very explicit in Sidre
          to allow maximum flexibility to compose complex operations from 
          simpler ones. This should become apparent by following the examples. 

The code snippets shown and described here exist in the file 
``axom/src/axom/sidre/examples/sidre_data_vs_metadata.cpp``, which can be 
built and run to experiment with if you wish.

The starting point for the first three examples is a simple Sidre datastore 
group hierarchy in which the root group has two child groups named "A" and "B". 
The generation of this structure is shown in the following code along with
some print statements to verify the state is what we expect.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex_datastore_initial_start
   :end-before: _ex_datastore_initial_end
   :language: C++

As expected, the datastore has no buffers and the groups have no views,
which we can see from the output of the code::

  Datastore start state....
        Num buffers in datastore: 0
        Num views in group A: 0
        Num views in group B: 0

--------------------------------
Example 1: One View, One Buffer 
--------------------------------

The first example shows a common Sidre usage pattern in which a view is created in a group and the data associated with the view is managed via the view. This 
case represents a one-to-one relationship between a view and a buffer.

We begin by creating a view named "aview" in group "A" describing an integer 
array of length 10 and allocate it, all in one method call. This creates a 
buffer in the datastore which holds the array data. Then, we get a pointer to 
the start of the array from the view and initialize the array values. To give 
some insight into how Sidre works, we access and print various pieces of 
information about the state of the group, view, and buffer. We also print the 
value of the buffer array at slot 5 to confirm that the buffer is indeed 
holding the view's data. Here's the relevant code section.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_create_start
   :end-before: _ex1_oneview_onebuffer_create_end
   :language: C++

The output printed by the code is::

  After A_grp->createViewAndAllocate() call
        Num views in group A: 1
        Num elements in view: 10
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10

  After initialization of view array
        Value of elt 5 in buffer array (expect 7): 7

Next, we deallocate the view and show that its description remains intact; 
for example, the number of elements is still 10. The view is still attached
to the buffer. The view is no longer allocated since the buffer data was 
deallocated. The buffer's data description remains intact. 

Then, we allocate the view again. Since the buffer is still attached to the
view, a new buffer is not created in the datastore and the existing buffer is 
re-allocated. The same data description as before is used for the allocation
since we haven't changed it. We verify this in the code output which follows
the code snippet (i.e., the datastore still has one buffer).

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_deallocalloc_start
   :end-before: _ex1_oneview_onebuffer_deallocalloc_end
   :language: C++

The output of the code is::

  After view deallocate call, the data no longer exists,
  but the view description remains.
        Num elements in view: 10
        View has buffer? 1
        Is view allocated? 0
        Num buffers in datastore: 1
        Is buffer allocated? 0
        Num views attached to buffer: 1
        Num elements in buffer array: 10

  After allocating view again...
        Num buffers in datastore: 1
        Is buffer allocated? 1
        Is view allocated? 1
        Num elements in view: 10 

Lastly, we destroy the view and its data with a single method call and verify
that the view and associated buffer no longer exist.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_destroy_start
   :end-before: _ex1_oneview_onebuffer_destroy_end
   :language: C++

The output of this portion of code is::

  After destroyViewAndData() call
        Num views in group A: 0
        Num buffers in datastore: 0

---------------------------------
Example 2: Two Views, One Buffer 
---------------------------------

The second example illustrates a Sidre usage pattern in which multiple views 
are created to describe distinct portions of data held in a shared buffer. 
This case represents a many-to-one relationship between views and a buffer. 
Before we start, we verify that the datastore contains no buffers and that 
the groups we created earlier contain no views.

We start by creating and allocating a *buffer* holding an array of doubles of 
length 10. We initialize the array so that each element has a value matching 
its position in the array; i.e., the values 0 through 9. Then, we create two 
views in group "A" and *attach* each to the buffer. Next, we apply a data 
description containing an offset and stride to each view so that one view 
is associated with the even values in the buffer and the other is associated
with the odd values. Accessing the data pointer in each view and printing
the values shows that this is indeed the case.

We call a method to destroy the first view and its data, similar to the
last part of the first example. The view is destroyed. However, since the 
buffer that held its data is shared by the other view, the buffer and its
data remain intact. In particular, the data associated with the remaining view
is also untouched.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex2_twoviews_onebuffer_start
   :end-before: _ex2_twoviews_onebuffer_end
   :language: C++

The output of the code is::

  Datastore start state
        Num buffers in datastore: 0
       	Num views in group A: 0
        Num views in group B: 0

  After buffer allocation and attaching to views
       	Num buffers in datastore: 1
       	Buffer num elements: 10

       	aview1 data has even values:
       	aview1 num elements: 5
       	aview1 offset: 0
       	aview1 stride: 2
       	aview1 data:	0   2   4   6   8   

       	aview2 data has odd values:
       	aview2 num elements: 5
       	aview2 offset: 1
       	aview2 stride: 2
       	aview2 data:	1   3   5   7   9   

  After destroyViewAndData(aview1) call
       	Num views in group A: 1
       	Num buffers in datastore: 1
       	Buffer num elements: 10
       	aview2 data still has its odd values:	1   3   5   7   9 

-----------------------------------------------
Example 3: Two Views and One Buffer (View Copy)
-----------------------------------------------

The third example illustrates another Sidre usage pattern in which multiple 
views share the same data in a single buffer. Before we start, we verify that 
the datastore contains one buffer and that the "A" group contains one view. 
This is the Sidre buffer, group, and view state that exists at the end of the 
previous example.

We begin by making a copy of the group "A" view "aview2" in group "B". The new 
view is identical in name, description, and data associated with it as the 
original view. In particular, the data in the new view *is the same data* 
associated with the original view. Recall that Sidre copy operations for 
groups and views are **shallow copy operations**. This means that a copy of 
a group or view is made in the destination group, but the data associated
with the copy is the same as in the original. We verify this by printing 
the values of the array associated with each view and also the base address 
of the array for each view.

Next, we destroy the "A" group which owned the original view. We verify that
the view copy remains in the "B" group and its data is still intact. When
we destroyed the "A" group, its view is also destroyed. So we can no longer 
access it with the usual method calls. If we maintained a handle (e.g., pointer)
to it, it would no longer be valid. 

Lastly, we destroy the "B" group. Similar to the destruction of the "A" group,
the view in the "B" group is destroyed. However, since we did not explicitly
delete (i.e., destroy) the data, we see from the code output below that the
buffer still exists in the datastore and is allocated. However, it has no 
attached views. Here is the complete example source code.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex3_twoviews_onebuffer_copy_start 
   :end-before: _ex3_twoviews_onebuffer_copy_end
   :language: C++

The output of the code is::

  Datastore start state
       	Num buffers in datastore: 1
       	Num views in group A: 1
       	Num views in group B: 0

  After copying aview2 to group B
       	Num buffers in datastore: 1
       	Num views in group A: 1
       	Num views in group B: 1
       	aview2 in A group has values:	1   3   5   7   9   
       	aview2 in B group has values:	1   3   5   7   9   

       	Base address of array in A group: 0xc595b8
       	Base address of array in A group: 0xc595b8

  After destroyGroup(A) call:
       	Num views in group B: 1
       	aview2 in B group has values:	1   3   5   7   9
       	Num buffers in datastore: 1
        Is buffer allocated? 1 
        Num views attached to buffer: 1

  After destroyGroup(B) call:
        Num buffers in datastore: 1
        Is buffer allocated? 1
        Num views attached to buffer: 0
  
The last operation in this example is intended to emphasize the explicit 
nature of Sidre methods. In particular, when a group is destroyed, its views 
are also destroyed, but their *data may remain intact*. At first, this may
seem unexpected. Such behavior was a design choice for Sidre to provide 
maximum flexibility in defining and manipulating views and data independently
(e.g., describing data and allocating it in separate phases of code execution)
while keeping Sidre's internal bookkeeping implementations reasonably simple.
The example following this one continues on this point.

.. note:: Object and data creation and destruction is very explicit in Sidre
          to allow maximum flexibility to compose complex operations from
          simpler ones. **Specific methods must be called to destroy views
          and deallocate their data.** 

-------------------------------------
Example 4: More Basic Mechanics
-------------------------------------

The last example should help to make clear the point made at the beginning of
this section and at the end of previous example about Sidre operations 
being explicit with respect to whether data associated with a group or view
is destroyed.

We start with a fresh datastore and create one group "A" in the root group. We
create a view "aview" in the group and allocate it. The output of the 
code below shows that everything is working as expected. 

We destroy the view and we see that the view is gone from the group, but
that the buffer is still in the datastore and allocated. This is so because
we did not explicitly destroy or deallocate the data. This is in contrast to 
the first example, where we destroyed the view and its data with a single 
method call designed for this purpose.

Next, we create the view again and attach the buffer to it. Then, we apply the
data description. This restores the state of everything before we destroyed 
the view.

Then, we deallocate the buffer and verify that the buffer and view are no 
longer allocated, but are still described. Also, the view is still attached 
to the buffer.

Lastly, we destroy the buffer. It has been deallocated and is gone from the 
datastore. The view remains, is deallocated, but is still described.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex4_more_mechanics_start
   :end-before: _ex4_more_mechanics_end
   :language: C++

The output of the code is::

  After A_grp->createViewAndAllocate() call:
        Num views in group A: 1
        Num elements in view: 10
        Is view allocated? 1
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10
        Is buffer allocated? 1

  After A_grp->destroyView() call:
        Num views in group A: 0
        Num buffers in datastore: 1
        Num views attached to buffer: 0
        Num elements in buffer array: 10
        Is buffer allocated? 1

  After recreating view, attaching buffer, and describing data:
        Num views in group A: 1
        Num elements in view: 10
        Is view allocated? 1
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10
        Is buffer allocated? 1

  After buffer deallocate call:
        Num elements in view: 10
        Is view allocated? 0
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10
        Is buffer allocated? 0

  After destroy buffer call:
        Num elements in view: 10
        Is view allocated? 0
        Num buffers in datastore: 0

The point of this final example is to emphasize the point that most Sidre 
operations are atomic and do precisely what their name describes with few side 
effects. If this were not the case, we believe that the Sidre API would become 
bloated and confusing to be able to support the range of usage patterns 
employed by Axom users.

