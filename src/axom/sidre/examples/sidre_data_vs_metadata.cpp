// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**************************************************************************
 *************************************************************************/

#include "axom/sidre.hpp"

using namespace axom;
using namespace sidre;

/**************************************************************************
 * Subroutine:  main
 * Purpose   : Provide examples showing when deleting a sidre group or view
 *             also deletes the associated data and when it does not.
 *************************************************************************/

/* 
 * This example code contains snippets used in the Sidre Sphinx documentation.
 * They begin and end with comments, such as
 *
 * _datastore_initial_start
 * _datastore_initial_end
 */

int main(int argc, char* argv[])
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  //
  // Create a datastore with a simple group hierarchy.
  //
  // _ex_datastore_initial_start
  DataStore* ds = new DataStore();
  Group* root_grp = ds->getRoot();

  Group* A_grp = root_grp->createGroup("A");
  Group* B_grp = root_grp->createGroup("B");

  //
  // Initially, datastore contains no buffer objects since no data has
  // been allocated yet. Also, the groups have no views.
  //
  std::cout << "Datastore start state...." << std::endl;
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  std::cout << "\tNum views in group B: " << B_grp->getNumViews() << std::endl;
  // _ex_datastore_initial_end

  std::cout << std::endl;

  // -----------------------------------------------------------------------
  // Example 1: One-to-one Buffer to View relationship
  // -----------------------------------------------------------------------
  //
  // Create a view with an integer array of length dat_size. Verify that the
  // group has one view and the datastore has one buffer (holding the
  // view's data).
  //
  // Initialize the elements of the integer array.
  //
  // Then, deallocate the view. The view and its description of the data
  // remains, but the buffer holding the array is deallocated since there
  // was only one view referencing its data.
  //
  // When the view is allocated again, the associated buffer is re-allocated
  // using the same data description as before
  //
  // Lastly, we destroy the view and data with a single method call.
  //

  std::cout << "Example 1: One-to-one Buffer to View relationship\n\n";

  // _ex1_oneview_onebuffer_create_start
  const int dat_size = 10;

  View* aview = A_grp->createViewAndAllocate("aview", INT_ID, dat_size);
  std::cout << "After A_grp->createViewAndAllocate() call\n";
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  axom::IndexType nelems = aview->getNumElements();
  std::cout << "\tNum elements in view: " << nelems << std::endl;

  Buffer* buf1 = aview->getBuffer();
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views attached to buffer: " << buf1->getNumViews() << std::endl;
  std::cout << "\tNum elements in buffer array: " << buf1->getNumElements() << std::endl;

  int* a_array = aview->getArray();
  for(axom::IndexType i = 0; i < nelems; ++i)
  {
    a_array[i] = i + 2;
  }

  std::cout << std::endl;

  std::cout << "After initialization of view array\n";
  int* buf1_ptr = buf1->getData();
  std::cout << "\tValue of elt 5 in buffer array (expect 7): " << buf1_ptr[5] << std::endl;
  // _ex1_oneview_onebuffer_create_end

  std::cout << std::endl;

  // _ex1_oneview_onebuffer_deallocalloc_start
  aview->deallocate();
  std::cout << "After view deallocate call, the data no longer exists,\n"
            << "but the view description remains." << std::endl;
  std::cout << "\tNum elements in view: " << aview->getNumElements() << std::endl;
  std::cout << "\tView has buffer? " << aview->hasBuffer() << std::endl;
  std::cout << "\tIs view allocated? " << aview->isAllocated() << std::endl;
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tIs buffer allocated? " << buf1->isAllocated() << std::endl;
  std::cout << "\tNum views attached to buffer: " << buf1->getNumViews() << std::endl;
  std::cout << "\tNum elements in buffer array: " << buf1->getNumElements() << std::endl;

  std::cout << std::endl;

  aview->allocate();
  std::cout << "After allocating view again...\n";
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tIs buffer allocated? " << buf1->isAllocated() << std::endl;
  std::cout << "\tIs view allocated? " << aview->isAllocated() << std::endl;
  std::cout << "\tNum elements in view: " << aview->getNumElements() << std::endl;
  // _ex1_oneview_onebuffer_deallocalloc_end

  std::cout << std::endl;

  // _ex1_oneview_onebuffer_destroy_start
  A_grp->destroyViewAndData("aview");
  std::cout << "After destroyViewAndData() call\n";
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  // _ex1_oneview_onebuffer_destroy_end

  // -----------------------------------------------------------------------
  // Example 2: One-to-many Buffer to View relationships
  // -----------------------------------------------------------------------
  //
  // Create a buffer with a double array of length dat_size and initialize its
  // data. Attach the buffer to two views, with data descriptions that vary in
  // offset. Verify that the data associated with each view is correct.
  //
  // Destroy one of the views. Verify that the data is still accessible via
  // the other view and that all data is still accessible via buffer.
  //
  std::cout << "\nExample 2: One-to-one Buffer to View relationships\n";

  // _ex2_twoviews_onebuffer_start
  std::cout << "\nDatastore start state\n";
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  std::cout << "\tNum views in group B: " << B_grp->getNumViews() << std::endl;

  Buffer* buf2 = ds->createBuffer(DOUBLE_ID, dat_size)->allocate();

  double* dat2 = buf2->getData();
  for(axom::IndexType i = 0; i < dat_size; ++i)
  {
    dat2[i] = i;
  }

  View* aview1 = A_grp->createView("aview1", buf2);
  View* aview2 = A_grp->createView("aview2", buf2);

  //
  // aview1 data gets even values, aview2 data gets odd values.
  //
  axom::IndexType view_nelem = dat_size / 2;

  aview1->apply(DOUBLE_ID, view_nelem, 0 /*offset*/, 2 /*stride*/);
  aview2->apply(DOUBLE_ID, view_nelem, 1 /*offset*/, 2 /*stride*/);

  std::cout << "\nAfter buffer allocation and attaching to views\n";
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tBuffer num elements: " << buf2->getNumElements() << std::endl;

  std::cout << "\n\taview1 data has even values:\n";
  std::cout << "\taview1 num elements: " << aview1->getNumElements() << std::endl;
  std::cout << "\taview1 offset: " << aview1->getOffset() << std::endl;
  std::cout << "\taview1 stride: " << aview1->getStride() << std::endl;
  double* arr1 = aview1->getArray();
  axom::IndexType vlen = aview1->getNumElements();
  axom::IndexType vstr = aview1->getStride();
  std::cout << "\taview1 data:\t";
  for(axom::IndexType i = 0; i < vlen * vstr; i += vstr)
  {
    std::cout << arr1[i] << "   ";
  }
  std::cout << std::endl;

  std::cout << "\n\taview2 data has odd values:\n";
  std::cout << "\taview2 num elements: " << aview2->getNumElements() << std::endl;
  std::cout << "\taview2 offset: " << aview2->getOffset() << std::endl;
  std::cout << "\taview2 stride: " << aview2->getStride() << std::endl;
  double* arr2 = aview2->getArray();
  vlen = aview2->getNumElements();
  vstr = aview2->getStride();
  std::cout << "\taview2 data:\t";
  for(axom::IndexType i = 0; i < vlen * vstr; i += vstr)
  {
    std::cout << arr2[i] << "   ";
  }
  std::cout << std::endl;

  A_grp->destroyViewAndData("aview1");

  std::cout << "\nAfter destroyViewAndData(aview1) call\n";
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tBuffer num elements: " << buf2->getNumElements() << std::endl;
  std::cout << "\taview2 data still has its odd values:\t";
  arr2 = aview2->getArray();
  vlen = aview2->getNumElements();
  vstr = aview2->getStride();
  for(axom::IndexType i = 0; i < vlen * vstr; i += vstr)
  {
    std::cout << arr2[i] << "   ";
  }
  // _ex2_twoviews_onebuffer_end
  std::cout << std::endl;

  // -----------------------------------------------------------------------
  // Example 3: One-to-many Buffer to View relationships (view copy)
  // -----------------------------------------------------------------------
  //
  // Create a copy of a group's view in a different group. Verify that the
  // views share the same data (i.e., they have the same base address, length,
  // offset, and stride).
  //
  // Then, destroy one of the groups and verify that the data is still
  // accessible via the view in the other group.
  //
  // Lastly, destroy the remaining group without explicitly destroying its
  // data. We see that the data buffer remains intact and allocated.
  //
  std::cout << "\nExample 3: One-to-many Buffer to View (view copy)\n";

  // _ex3_twoviews_onebuffer_copy_start
  std::cout << "\nDatastore start state\n";
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  std::cout << "\tNum views in group B: " << B_grp->getNumViews() << std::endl;

  B_grp->copyView(aview2);

  std::cout << "\nAfter copying aview2 to group B\n";
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  std::cout << "\tNum views in group B: " << B_grp->getNumViews() << std::endl;

  View* aview2_in_Agrp = A_grp->getView("aview2");
  View* aview2_in_Bgrp = B_grp->getView("aview2");

  std::cout << "\taview2 in A group has values:\t";
  double* arr_A = aview2_in_Agrp->getArray();
  vlen = aview2_in_Agrp->getNumElements();
  vstr = aview2_in_Agrp->getStride();
  for(axom::IndexType i = 0; i < vlen * vstr; i += vstr)
  {
    std::cout << arr_A[i] << "   ";
  }
  std::cout << std::endl;

  std::cout << "\taview2 in B group has values:\t";
  double* arr_B = aview2_in_Bgrp->getArray();
  vlen = aview2_in_Bgrp->getNumElements();
  vstr = aview2_in_Bgrp->getStride();
  for(axom::IndexType i = 0; i < vlen * vstr; i += vstr)
  {
    std::cout << arr_B[i] << "   ";
  }
  std::cout << std::endl;

  std::cout << std::endl;

  std::cout << "\tBase address of array in A group: " << arr_A << std::endl;
  std::cout << "\tBase address of array in B group: " << arr_B << std::endl;

  root_grp->destroyGroup("A");

  Buffer* buf_aview2 = aview2_in_Bgrp->getBuffer();

  std::cout << "\nAfter destroyGroup(A) call:\n";
  std::cout << "\tNum views in group B: " << B_grp->getNumViews() << std::endl;

  std::cout << "\taview2 in B group has values:\t";
  aview2_in_Bgrp = B_grp->getView("aview2");
  arr_B = aview2_in_Bgrp->getArray();
  vlen = aview2_in_Bgrp->getNumElements();
  vstr = aview2_in_Bgrp->getStride();
  for(axom::IndexType i = 0; i < vlen * vstr; i += vstr)
  {
    std::cout << arr_B[i] << "   ";
  }

  std::cout << std::endl;

  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tIs buffer allocated? " << buf_aview2->isAllocated() << std::endl;
  std::cout << "\tNum views attached to buffer: " << buf_aview2->getNumViews() << std::endl;
  std::cout << std::endl;

  root_grp->destroyGroup("B");

  std::cout << "\nAfter destroyGroup(B) call:\n";
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tIs buffer allocated? " << buf_aview2->isAllocated() << std::endl;
  std::cout << "\tNum views attached to buffer: " << buf_aview2->getNumViews() << std::endl;
  // _ex3_twoviews_onebuffer_copy_end
  std::cout << std::endl;

  //
  // Destroy datastore and all contents to start fresh for last example...
  //
  delete ds;

  // -----------------------------------------------------------------------
  // Example 4: More on basic mechanics
  // -----------------------------------------------------------------------
  //
  // We first create a new datastore and one group "A" in the root group.
  // Then, we create and allocate a view in that group, and grab a pointer to
  // the associated buffer.
  //
  // When we destroy the view (but not its data as we did in the first example),
  // we verify that the buffer is still in the datastore, described and
  // allocated.
  //
  // We recreate the view, attach the buffer to it, and describe the view
  // data as before. We verify that the view is allocated since its buffer is.
  //
  // When we deallocate the buffer, we see that it remains in the datastore
  // and that the view and buffer are deallocated, but their descriptions
  // remain.
  //
  // Lastly, we destroy the buffer and see that the datastore has zero buffers.
  // However, the view remains in the group described and unallocated.
  //
  std::cout << "\nExample 4: More on basic mechanics\n";

  // _ex4_more_mechanics_start
  ds = new DataStore();
  root_grp = ds->getRoot();

  A_grp = root_grp->createGroup("A");

  aview = A_grp->createViewAndAllocate("aview", INT_ID, dat_size);
  std::cout << "\nAfter A_grp->createViewAndAllocate() call:\n";
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  nelems = aview->getNumElements();
  std::cout << "\tNum elements in view: " << nelems << std::endl;
  std::cout << "\tIs view allocated? " << aview->isAllocated() << std::endl;

  Buffer* bufa = aview->getBuffer();
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views attached to buffer: " << bufa->getNumViews() << std::endl;
  std::cout << "\tNum elements in buffer array: " << bufa->getNumElements() << std::endl;
  std::cout << "\tIs buffer allocated? " << bufa->isAllocated() << std::endl;

  std::cout << std::endl;

  A_grp->destroyView("aview");
  std::cout << "After A_grp->destroyView() call:\n";
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views attached to buffer: " << bufa->getNumViews() << std::endl;
  std::cout << "\tNum elements in buffer array: " << bufa->getNumElements() << std::endl;
  std::cout << "\tIs buffer allocated? " << bufa->isAllocated() << std::endl;

  std::cout << std::endl;

  aview = A_grp->createView("aview1", bufa);
  aview->apply(INT_ID, dat_size);
  std::cout << "After recreating view, attaching buffer, and describing data:\n";
  std::cout << "\tNum views in group A: " << A_grp->getNumViews() << std::endl;
  nelems = aview->getNumElements();
  std::cout << "\tNum elements in view: " << nelems << std::endl;
  std::cout << "\tIs view allocated? " << aview->isAllocated() << std::endl;

  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views attached to buffer: " << bufa->getNumViews() << std::endl;
  std::cout << "\tNum elements in buffer array: " << bufa->getNumElements() << std::endl;
  std::cout << "\tIs buffer allocated? " << bufa->isAllocated() << std::endl;

  std::cout << std::endl;

  bufa->deallocate();
  std::cout << "After buffer deallocate call:\n";
  nelems = aview->getNumElements();
  std::cout << "\tNum elements in view: " << nelems << std::endl;
  std::cout << "\tIs view allocated? " << aview->isAllocated() << std::endl;

  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  std::cout << "\tNum views attached to buffer: " << bufa->getNumViews() << std::endl;
  std::cout << "\tNum elements in buffer array: " << bufa->getNumElements() << std::endl;
  std::cout << "\tIs buffer allocated? " << bufa->isAllocated() << std::endl;

  std::cout << std::endl;

  ds->destroyBuffer(bufa);
  std::cout << "After destroy buffer call:\n";
  nelems = aview->getNumElements();
  std::cout << "\tNum elements in view: " << nelems << std::endl;
  std::cout << "\tIs view allocated? " << aview->isAllocated() << std::endl;

  std::cout << "\tNum buffers in datastore: " << ds->getNumBuffers() << std::endl;
  // _ex4_more_mechanics_end

  //
  // Clean up...
  //
  delete ds;

  return 0;
}
