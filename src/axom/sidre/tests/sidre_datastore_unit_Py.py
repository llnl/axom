# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import random

def verify_empty_group_named(dg, name):
	assert dg.getName() == name

	assert dg.getNumGroups() == 0
	assert not dg.hasGroup(-1)
	assert not dg.hasGroup(0)
	assert not dg.hasGroup(1)
	assert not dg.hasGroup("some_name")
	assert dg.getGroupIndex("some_other_name") == pysidre.InvalidIndex
	assert dg.getFirstValidGroupIndex() == pysidre.InvalidIndex
	assert dg.getNextValidGroupIndex(0) == pysidre.InvalidIndex
	assert dg.getNextValidGroupIndex(4) == pysidre.InvalidIndex

	assert dg.getNumViews() == 0
	assert not dg.hasView(-1)
	assert not dg.hasView(0)
	assert not dg.hasView(1)
	assert not dg.hasView("some_name")
	assert dg.getViewIndex("some_other_name") == pysidre.InvalidIndex
	assert dg.getFirstValidViewIndex() == pysidre.InvalidIndex
	assert dg.getNextValidViewIndex(0) == pysidre.InvalidIndex
	assert dg.getNextValidViewIndex(4) == pysidre.InvalidIndex

def verify_buffer_identity(ds, bs):
	bufcount = len(bs)

	# Does ds contain the number of buffers we expect?
	assert ds.getNumBuffers() == bufcount

	# Does ds contain the buffer IDs and pointers we expect?
	iterated_count = 0
	idx = ds.getFirstValidBufferIndex()
	while idx != pysidre.InvalidIndex and iterated_count < bufcount:
		assert idx in bs
		if idx in bs:
			assert bs[idx] == ds.getBuffer(idx)
		idx = ds.getNextValidBufferIndex(idx)
		iterated_count += 1

	# Have we iterated over exactly the number of buffers we expect, finishing on InvalidIndex?
	assert iterated_count == bufcount
	assert idx == pysidre.InvalidIndex


def test_default_ctor():
	ds = pysidre.DataStore()

	# After construction, the DataStore should contain no buffers.
	assert ds.getNumBuffers() == 0
	assert not ds.hasBuffer(-15)
	assert not ds.hasBuffer(-1)
	assert not ds.hasBuffer(0)
	assert not ds.hasBuffer(1)
	assert not ds.hasBuffer(8)

	assert ds.getFirstValidBufferIndex() == pysidre.InvalidIndex
	assert ds.getNextValidBufferIndex(0) == pysidre.InvalidIndex
	assert ds.getNextValidBufferIndex(4) == pysidre.InvalidIndex

	# The new DataStore should contain exactly one group, the root group.
	# The root group should be named "" and should contain no views and no groups.
	dg = ds.getRoot()

	assert dg is not None
	assert dg == dg.getParent()
	assert dg.getDataStore() == ds

	verify_empty_group_named(dg, "")


# The dtor destroys all buffers and deletes the root group.
# An outside tool should be used to check for proper memory cleanup.
def test_create_destroy_buffers_basic():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    # Basic tests
    dbuff = ds.createBuffer()
    assert ds.getNumBuffers() == 1

    buffer_index = ds.getFirstValidBufferIndex()
    assert dbuff.getIndex() == 0
    assert dbuff.getIndex() == buffer_index
    assert ds.getNextValidBufferIndex(buffer_index) == pysidre.InvalidIndex

    # Do we get the buffer we expect?
    assert dbuff == ds.getBuffer(buffer_index)
    bad_buffer_index = 9999
    assert ds.getBuffer(bad_buffer_index) is None

    ds.destroyBuffer(buffer_index)
    # should be no buffers
    assert ds.getNumBuffers() == 0
    assert ds.getFirstValidBufferIndex() == pysidre.InvalidIndex
    assert not ds.hasBuffer(buffer_index)
    assert ds.getBuffer(buffer_index) is None
    assert ds.getBuffer(bad_buffer_index) is None


def test_create_destroy_buffers_order():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    dbuff = ds.createBuffer()
    assert ds.getNumBuffers() == 1

    buffer_index = ds.getFirstValidBufferIndex()
    ds.destroyBuffer(dbuff)

    # After destroy, test that buffer index should be available again for reuse.
    dbuff2 = ds.createBuffer(pysidre.TypeID.FLOAT32_ID, 16)
    d2_index = dbuff2.getIndex()
    assert ds.getFirstValidBufferIndex() == buffer_index
    assert d2_index == buffer_index
    assert ds.hasBuffer(d2_index)
    assert ds.getBuffer(buffer_index) == dbuff2

    dbuff3 = ds.createBuffer()
    d3_index = dbuff3.getIndex()
    assert ds.getNumBuffers() == 2
    assert ds.hasBuffer(d3_index)

    # Try destroying the first valid buffer; see if we have the correct count and indices
    ds.destroyBuffer(buffer_index)
    assert ds.getNumBuffers() == 1
    assert not ds.hasBuffer(d2_index)
    assert ds.hasBuffer(d3_index)

    # Add some more buffers, then try destroying the second one; see if we have the correct count and indices
    dbuff4 = ds.createBuffer()
    dbuff5 = ds.createBuffer()
    d4_index = dbuff4.getIndex()
    d5_index = dbuff5.getIndex()
    assert ds.getNumBuffers() == 3
    assert d4_index == buffer_index  # dbuff4 should have recycled index 0
    assert ds.hasBuffer(d3_index)
    assert ds.hasBuffer(d4_index)
    assert ds.hasBuffer(d5_index)

    # Destroy dbuff3 (not dbuff4) because we already tested destroying index 0
    ds.destroyBuffer(d3_index)
    assert ds.getNumBuffers() == 2
    assert not ds.hasBuffer(d3_index)
    assert ds.hasBuffer(d4_index)
    assert ds.hasBuffer(d5_index)

    # Can we destroy all buffers?
    ds.destroyAllBuffers()
    assert ds.getNumBuffers() == 0
    assert not ds.hasBuffer(d2_index)
    assert not ds.hasBuffer(d4_index)
    assert not ds.hasBuffer(d5_index)


def test_create_destroy_buffers_views():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    dbuff3 = ds.createBuffer()
    d3_index = dbuff3.getIndex()
    dbuff4 = ds.createBuffer()
    d4_index = dbuff4.getIndex()
    dbuff5 = ds.createBuffer()
    d5_index = dbuff5.getIndex()
    dbuff6 = ds.createBuffer()
    d6_index = dbuff6.getIndex()

    assert ds.getBuffer(d3_index) == dbuff3
    assert ds.getBuffer(d4_index) == dbuff4
    assert ds.getBuffer(d5_index) == dbuff5
    assert ds.getBuffer(d6_index) == dbuff6

    # Create and verify views referencing buffers
    vA = ds.getRoot().createView("vA", dbuff3)
    vB = ds.getRoot().createView("vB", dbuff3)
    vC = ds.getRoot().createView("vC", dbuff4)
    vD = ds.getRoot().createView("vD", dbuff6)
    vE = ds.getRoot().createView("vE", dbuff6)
    assert vA.getBuffer() == dbuff3
    assert vB.getBuffer() == dbuff3
    assert dbuff3.getNumViews() == 2
    assert vC.getBuffer() == dbuff4
    assert vD.getBuffer() == dbuff6
    assert vE.getBuffer() == dbuff6
    assert dbuff6.getNumViews() == 2

    # Destroying a buffer should detach it from the view
    ds.destroyBuffer(d4_index)
    assert ds.getNumBuffers() == 3
    assert ds.hasBuffer(d3_index)
    assert not ds.hasBuffer(d4_index)
    assert ds.hasBuffer(d5_index)
    assert ds.hasBuffer(d6_index)
    assert vA.getBuffer() == dbuff3
    assert vB.getBuffer() == dbuff3
    assert not vC.hasBuffer()
    assert vD.getBuffer() == dbuff6
    assert vE.getBuffer() == dbuff6

    # Destroying a buffer should detach it from all of its views
    ds.destroyBuffer(d3_index)
    assert ds.getNumBuffers() == 2
    assert not ds.hasBuffer(d3_index)
    assert not ds.hasBuffer(d4_index)
    assert ds.hasBuffer(d5_index)
    assert ds.hasBuffer(d6_index)
    assert not vA.hasBuffer()
    assert not vB.hasBuffer()
    assert not vC.hasBuffer()
    assert vD.getBuffer() == dbuff6
    assert vE.getBuffer() == dbuff6

    # Destroying all buffers should detach them from all of their views
    dbuff3 = ds.createBuffer()
    dbuff4 = ds.createBuffer()
    vA.attachBuffer(dbuff3)
    vB.attachBuffer(dbuff3)
    vC.attachBuffer(dbuff4)
    ds.destroyAllBuffers()
    assert ds.getNumBuffers() == 0
    assert not vA.hasBuffer()
    assert not vB.hasBuffer()
    assert not vC.hasBuffer()
    assert not vD.hasBuffer()
    assert not vE.hasBuffer()


def psrand(min_val, max_val):
    # Returns a pseudorandom int in [min_val, max_val] (closed interval)
    return random.randint(min_val, max_val)


def irhall(n):
    # Approximates the Irwin-Hall distribution, a sum of uniformly-distributed
    # samples that approaches Gaussian distribution
    retval = 0
    for i in range(n):
        retval += psrand(-5, 5)
    return retval


# Test iteration through buffers, as well as proper index and buffer behavior
# while buffers are created and deleted
def test_iterate_buffers_basic():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    bad_buffer_index = 9999
    # Do we get sidre::InvalidIndex for several queries with no buffers?
    assert ds.getFirstValidBufferIndex() == pysidre.InvalidIndex
    assert ds.getNextValidBufferIndex(0) == pysidre.InvalidIndex
    assert ds.getNextValidBufferIndex(bad_buffer_index) == pysidre.InvalidIndex
    assert ds.getNextValidBufferIndex(pysidre.InvalidIndex) == pysidre.InvalidIndex

    # Create one data buffer, verify its index is zero, and that iterators behave as expected
    initial = ds.createBuffer()
    assert initial.getIndex() == 0
    assert ds.getNumBuffers() == 1
    assert ds.getFirstValidBufferIndex() == 0
    assert ds.getNextValidBufferIndex(0) == pysidre.InvalidIndex

    # Destroy the data buffer, verify that iterators behave as expected
    ds.destroyBuffer(initial)
    assert ds.getNumBuffers() == 0
    assert ds.getFirstValidBufferIndex() == pysidre.InvalidIndex
    assert ds.getNextValidBufferIndex(0) == pysidre.InvalidIndex


def test_iterate_buffers_simple():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    bs = {}
    bufcount = 20

    for i in range(bufcount):
        b = ds.createBuffer(pysidre.TypeID.FLOAT64_ID, 400 * i)
        idx = b.getIndex()
        bs[idx] = b

    verify_buffer_identity(ds, bs)


def test_iterate_buffers_iterators():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    bs = {}
    bufcount = 20

    for i in range(bufcount):
        b = ds.createBuffer(pysidre.TypeID.FLOAT64_ID, 400 * i)
        idx = b.getIndex()
        bs[idx] = b

    found_buffers = 0
    for buff in ds.buffers():
        idx = buff.getIndex()
        found_buffers += 1
        assert pysidre.indexIsValid(idx)
        assert ds.getBuffer(idx) == buff
        assert bs[idx] == buff
    assert found_buffers == bufcount


# Test creating and allocating buffers, then destroying several of them
def test_create_delete_buffers_iterate():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    bs = {}
    bufcount = 50  # Arbitrary number of buffers

    # Initially, create some buffers of varying size
    for i in range(bufcount):
        b = ds.createBuffer(pysidre.TypeID.FLOAT64_ID, (400 * i) % 10000)
        b.allocate()
        idx = b.getIndex()
        bs[idx] = b

    nbs = {}
    for i, (idx, buf) in enumerate(bs.items()):
        # Eliminate some buffers (arbitrarily chosen)
        if i % 5 and i % 7:
            nbs[idx] = buf
        else:
            ds.destroyBuffer(idx)

    verify_buffer_identity(ds, nbs)


def test_iterate_buffers_with_delete_iterators():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    init_buff_count = 22
    for i in range(init_buff_count):
        ds.createBuffer(pysidre.TypeID.FLOAT64_ID, 400 * i)
    assert ds.getNumBuffers() == init_buff_count

    # Remove a few buffers by index
    for idx in [5, 10, 15, 20, 25]:
        if ds.hasBuffer(idx):
            ds.destroyBuffer(idx)
    exp_buff_count = 18
    assert ds.getNumBuffers() == exp_buff_count

    # Add a buffer, expect it to reuse a lower index
    assert not ds.hasBuffer(5)
    buff = ds.createBuffer(pysidre.TypeID.FLOAT64_ID, 10)
    idx = buff.getIndex()
    assert idx < init_buff_count
    assert ds.hasBuffer(idx)
    exp_buff_count = 19
    assert ds.getNumBuffers() == exp_buff_count

    # Remove a few more buffers, some by object
    for idx in [3, 6, 9, 12, 15, 18, 21]:
        if ds.hasBuffer(idx):
            buff = ds.getBuffer(idx)
            ds.destroyBuffer(buff)
    exp_buff_count = 13
    assert ds.getNumBuffers() == exp_buff_count

    # Iterate using standard Python for-loop
    found_buffers = 0
    for buff in ds.buffers():
        idx = buff.getIndex()
        found_buffers += 1
        assert pysidre.indexIsValid(idx)
        assert ds.getBuffer(idx) == buff
    assert found_buffers == exp_buff_count


# Test creating+allocating buffers, then destroying several of them, repeatedly
def test_loop_create_delete_buffers_iterate():
    ds = pysidre.DataStore()
    assert ds.getNumBuffers() == 0

    bs = {}
    idxlist = []
    initbufcount = 50 # Arbitrary number of buffers

    # Initially, create some buffers of varying size
    for i in range(initbufcount):
        b = ds.createBuffer(pysidre.TypeID.FLOAT64_ID, (400 * i) % 10000)
        b.allocate()
        idx = b.getIndex()
        bs[idx] = b
        idxlist.append(idx)

    totalrounds = 100 # Arbitrary number of rounds
    for round in range(totalrounds):
        # In each round, choose a random number of buffers to delete or create
        delta = irhall(5)
        if delta < 0:
            rmvcount = abs(delta)
            rmvcount = min(rmvcount, len(bs))
            for i in range(rmvcount):
                if not idxlist:
                    break
                rmvidx = psrand(0, len(idxlist) - 1)
                rmvid = idxlist[rmvidx]
                assert ds.hasBuffer(rmvid)
                assert rmvid in bs
                ds.destroyBuffer(rmvid)
                bs.pop(rmvid)
                idxlist.pop(rmvidx)
                assert not ds.hasBuffer(rmvid)
                assert rmvid not in bs
        elif delta > 0:
            addcount = delta
            for _ in range(addcount):
                buf = ds.createBuffer(pysidre.TypeID.FLOAT64_ID, 400)
                buf.allocate()
                addid = buf.getIndex()
                assert ds.hasBuffer(addid)
                assert addid not in bs
                bs[addid] = buf
                idxlist.append(addid)

        verify_buffer_identity(ds, bs)