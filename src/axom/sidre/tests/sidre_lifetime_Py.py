# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Lifetime-soundness regression tests for the sidre python bindings.

Each test obtains a sidre-owned object (child proxy, ancestor proxy, harvested
iterator element, or zero-copy numpy array), drops every Python owner, forces a
garbage collection, and then *uses* the object. Before the lifetime audit these
patterns dereferenced freed memory and segfaulted; with reference_internal on
owner-chain accessors, keep_alive on iterator elements, and self-as-owner on
returned arrays, the keep_alive graph keeps the backing DataStore alive and the
accesses are safe.

These tests therefore only "pass" against the audited bindings; against the
prior bindings they crash the interpreter (the failure mode the audit fixes).
"""

import gc
import weakref

import numpy as np
import pytest

import pysidre


def _force_gc():
    # A couple of passes to be robust to reference cycles in the proxies.
    for _ in range(3):
        gc.collect()


# ---------------------------------------------------------------------------
# Child proxies must pin their owner chain (parent -> owned child)
# ---------------------------------------------------------------------------
def test_root_outlives_datastore():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    del ds
    _force_gc()
    # root must keep its DataStore alive
    assert root.getNumViews() == 0
    root.createView("v")
    assert root.hasView("v")


def test_child_group_outlives_datastore():
    ds = pysidre.DataStore()
    grp = ds.getRoot().createGroup("a/b/c")
    del ds
    _force_gc()
    # grp pins parent group pins root pins DataStore
    assert grp.getName() == "c"
    assert grp.getPathName() == "a/b/c"


def test_view_outlives_datastore():
    ds = pysidre.DataStore()
    view = ds.getRoot().createViewAndAllocate("field", pysidre.TypeID.INT_ID, 4)
    del ds
    _force_gc()
    assert view.getNumElements() == 4
    assert view.getName() == "field"


def test_buffer_outlives_datastore():
    ds = pysidre.DataStore()
    buff = ds.createBuffer(pysidre.TypeID.INT_ID, 8)
    buff.allocate()
    del ds
    _force_gc()
    assert buff.getNumElements() == 8


# ---------------------------------------------------------------------------
# Ancestor proxies (child -> ancestor) must pin the object they were minted from
# ---------------------------------------------------------------------------
def test_owning_group_outlives_datastore():
    ds = pysidre.DataStore()
    view = ds.getRoot().createView("v")
    owner = view.getOwningGroup()
    del ds
    del view
    _force_gc()
    assert owner.hasView("v")


def test_get_datastore_back_reference():
    ds = pysidre.DataStore()
    grp = ds.getRoot().createGroup("child")
    back = grp.getDataStore()
    del ds
    del grp
    _force_gc()
    # back-reference keeps the store alive; reach a fresh proxy through it
    assert back.getRoot().hasGroup("child")


def test_view_buffer_back_reference():
    ds = pysidre.DataStore()
    view = ds.getRoot().createViewAndAllocate("field", pysidre.TypeID.INT_ID, 4)
    buff = view.getBuffer()
    del ds
    del view
    _force_gc()
    assert buff.getNumElements() == 4


# ---------------------------------------------------------------------------
# Iterator elements harvested into a list must outlive the collection + store
# ---------------------------------------------------------------------------
def test_harvested_views_outlive_datastore():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    for i in range(4):
        root.createView(f"v{i}")
    harvested = list(root.views())
    del ds
    del root
    _force_gc()
    names = sorted(v.getName() for v in harvested)
    assert names == ["v0", "v1", "v2", "v3"]


def test_harvested_groups_outlive_datastore():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    for i in range(3):
        root.createGroup(f"g{i}")
    harvested = list(root.groups())
    del ds
    del root
    _force_gc()
    names = sorted(g.getName() for g in harvested)
    assert names == ["g0", "g1", "g2"]


def test_harvested_buffers_outlive_datastore():
    ds = pysidre.DataStore()
    for _ in range(3):
        ds.createBuffer(pysidre.TypeID.INT_ID, 2).allocate()
    harvested = list(ds.buffers())
    del ds
    _force_gc()
    assert sorted(b.getNumElements() for b in harvested) == [2, 2, 2]


def test_harvested_attributes_outlive_datastore():
    ds = pysidre.DataStore()
    ds.createAttributeString("a0", "x")
    ds.createAttributeScalar("a1", 42)
    harvested = list(ds.attributes())
    del ds
    _force_gc()
    assert sorted(a.getName() for a in harvested) == ["a0", "a1"]


# ---------------------------------------------------------------------------
# Iterator adaptors must outlive the owning Group/DataStore
# ---------------------------------------------------------------------------
def test_views_adaptor_outlives_group_and_datastore():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    for i in range(3):
        root.createView(f"v{i}")
    adaptor = root.views()
    del root
    del ds
    _force_gc()
    assert sorted(v.getName() for v in adaptor) == ["v0", "v1", "v2"]


def test_groups_adaptor_outlives_group_and_datastore():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    for i in range(3):
        root.createGroup(f"g{i}")
    adaptor = root.groups()
    del root
    del ds
    _force_gc()
    assert sorted(g.getName() for g in adaptor) == ["g0", "g1", "g2"]


def test_buffers_adaptor_outlives_datastore():
    ds = pysidre.DataStore()
    for _ in range(3):
        ds.createBuffer(pysidre.TypeID.INT_ID, 2).allocate()
    adaptor = ds.buffers()
    del ds
    _force_gc()
    assert sorted(b.getNumElements() for b in adaptor) == [2, 2, 2]


def test_attributes_adaptor_outlives_datastore():
    ds = pysidre.DataStore()
    ds.createAttributeString("a0", "x")
    ds.createAttributeScalar("a1", 42)
    adaptor = ds.attributes()
    del ds
    _force_gc()
    assert sorted(a.getName() for a in adaptor) == ["a0", "a1"]


# ---------------------------------------------------------------------------
# Lookup accessors should return proxies that pin their owner chain
# ---------------------------------------------------------------------------
def test_get_view_outlives_group_and_datastore():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    root.createView("v")
    view = root.getView("v")
    del root
    del ds
    _force_gc()
    assert view.getName() == "v"


def test_get_group_outlives_group_and_datastore():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    root.createGroup("g")
    grp = root.getGroup("g")
    del root
    del ds
    _force_gc()
    assert grp.getName() == "g"


def test_get_buffer_outlives_datastore():
    ds = pysidre.DataStore()
    ds.createBuffer(pysidre.TypeID.INT_ID, 7).allocate()
    buff = ds.getBuffer(0)
    del ds
    _force_gc()
    assert buff.getNumElements() == 7


def test_get_attribute_outlives_datastore():
    ds = pysidre.DataStore()
    ds.createAttributeString("a0", "x")
    attr = ds.getAttribute("a0")
    del ds
    _force_gc()
    assert attr.getName() == "a0"


def test_parent_group_outlives_datastore():
    ds = pysidre.DataStore()
    child = ds.getRoot().createGroup("a/b")
    parent = child.getParent()
    del ds
    del child
    _force_gc()
    assert parent.getName() == "a"
    assert parent.hasGroup("b")


def test_moved_group_outlives_owner_chain():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    src = root.createGroup("src")
    dst = root.createGroup("dst")
    child = src.createGroup("child")

    moved = dst.moveGroup(child)
    assert moved.getPathName() == "dst/child"

    del ds
    del root
    del src
    del dst
    del child
    _force_gc()

    assert moved.getName() == "child"
    assert moved.getPathName() == "dst/child"


# ---------------------------------------------------------------------------
# Zero-copy numpy arrays must pin their backing View / Buffer (and DataStore)
# ---------------------------------------------------------------------------
def test_view_array_outlives_datastore():
    ds = pysidre.DataStore()
    view = ds.getRoot().createViewAndAllocate("field", pysidre.TypeID.INT_ID, 4)
    arr = view.getDataArray()
    arr[:] = [10, 20, 30, 40]
    del ds
    del view
    _force_gc()
    # array still references valid memory: read and write
    assert list(arr) == [10, 20, 30, 40]
    arr[0] = 99
    assert arr[0] == 99


def test_buffer_array_outlives_datastore():
    ds = pysidre.DataStore()
    buff = ds.createBuffer(pysidre.TypeID.INT_ID, 4)
    buff.allocate()
    arr = buff.getDataArray()
    arr[:] = [1, 2, 3, 4]
    del ds
    del buff
    _force_gc()
    assert list(arr) == [1, 2, 3, 4]


def test_view_array_survives_owner_chain_collection():
    # Keep only the array; let the entire DataStore/Group/View chain be dropped.
    def make_array():
        ds = pysidre.DataStore()
        view = ds.getRoot().createViewAndAllocate("field", pysidre.TypeID.FLOAT64_ID, 5)
        a = view.getDataArray()
        a[:] = np.arange(5, dtype=np.float64)
        return a

    arr = make_array()
    _force_gc()
    np.testing.assert_array_equal(arr, np.arange(5, dtype=np.float64))


# ---------------------------------------------------------------------------
# External numpy storage borrowed by Sidre must stay alive with the C++ View
# ---------------------------------------------------------------------------
def test_create_view_external_array_owner_survives_discarded_proxy():
    ds = pysidre.DataStore()
    root = ds.getRoot()

    def create_external_view():
        external = np.arange(6, dtype=np.int64)
        ref = weakref.ref(external)
        root.createView("external", external).apply(pysidre.TypeID.INT64_ID, 6)
        return ref

    ref = create_external_view()
    _force_gc()
    assert ref() is not None
    np.testing.assert_array_equal(root.getView("external").getDataArray(), np.arange(6))


def test_create_view_with_shape_external_array_owner_survives_discarded_proxy():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    shape = np.array([2, 3])

    def create_external_view():
        external = np.arange(6, dtype=np.int64)
        ref = weakref.ref(external)
        root.createViewWithShape("shaped", pysidre.TypeID.INT64_ID, 2, shape, external)
        return ref

    ref = create_external_view()
    _force_gc()
    assert ref() is not None
    np.testing.assert_array_equal(root.getView("shaped").getDataArray(), np.arange(6).reshape(2, 3))


def test_set_external_data_array_owner_survives_discarded_proxy():
    ds = pysidre.DataStore()
    root = ds.getRoot()

    def set_external_data():
        view = root.createView("external")
        external = np.arange(6, dtype=np.int64)
        ref = weakref.ref(external)
        view.setExternalData(pysidre.TypeID.INT64_ID, 6, external)
        return ref

    ref = set_external_data()
    _force_gc()
    assert ref() is not None
    np.testing.assert_array_equal(root.getView("external").getDataArray(), np.arange(6))


def test_set_external_data_with_shape_array_owner_survives_discarded_proxy():
    ds = pysidre.DataStore()
    root = ds.getRoot()
    shape = np.array([2, 3])

    def set_external_data():
        view = root.createView("shaped")
        external = np.arange(6, dtype=np.int64)
        ref = weakref.ref(external)
        view.setExternalData(pysidre.TypeID.INT64_ID, 2, shape, external)
        return ref

    ref = set_external_data()
    _force_gc()
    assert ref() is not None
    np.testing.assert_array_equal(root.getView("shaped").getDataArray(), np.arange(6).reshape(2, 3))


def test_clear_releases_external_array_owner():
    ds = pysidre.DataStore()
    view = ds.getRoot().createView("external")
    external = np.arange(6, dtype=np.int64)
    ref = weakref.ref(external)
    view.setExternalData(pysidre.TypeID.INT64_ID, 6, external)
    del external
    _force_gc()
    assert ref() is not None

    view.clear()
    _force_gc()
    assert ref() is None


def test_copy_view_with_external_data_preserves_pin():
    """copyView on an external View should copy the pin to prevent premature collection."""
    ds = pysidre.DataStore()
    root = ds.getRoot()
    src_group = root.createGroup("src")
    dst_group = root.createGroup("dst")

    # Create view with external data
    external = np.arange(10, dtype=np.int64)
    ref = weakref.ref(external)
    src_view = src_group.createView("original", external)
    src_view.apply(pysidre.TypeID.INT64_ID, 10)
    del external
    _force_gc()
    assert ref() is not None  # Pin keeps it alive

    # Copy the view - should copy the pin
    copied_view = dst_group.copyView(src_view)
    assert copied_view.isExternal()

    # Delete source view and DS - copied view should keep pin alive
    del src_view
    del src_group
    del ds
    del root
    del dst_group
    _force_gc()

    # Pin should still be valid
    assert ref() is not None
    np.testing.assert_array_equal(copied_view.getDataArray(), np.arange(10))


def test_copy_group_with_external_data_preserves_pins():
    """copyGroup should recursively copy pins for all external Views in hierarchy."""
    ds = pysidre.DataStore()
    root = ds.getRoot()
    src = root.createGroup("source")

    # Create nested groups with external data
    external1 = np.arange(5, dtype=np.int32)
    external2 = np.arange(8, dtype=np.int64)
    ref1 = weakref.ref(external1)
    ref2 = weakref.ref(external2)

    src.createView("view1", external1).apply(pysidre.TypeID.INT32_ID, 5)
    child = src.createGroup("child")
    child.createView("view2", external2).apply(pysidre.TypeID.INT64_ID, 8)

    del external1
    del external2
    _force_gc()
    assert ref1() is not None
    assert ref2() is not None

    # Copy the entire group hierarchy - copyGroup creates a group with the source name
    dst_parent = root.createGroup("copies")
    dst = dst_parent.copyGroup(src)
    assert dst is not None
    assert dst.getName() == "source"
    assert dst.hasView("view1")
    assert dst.hasGroup("child")
    assert dst.getGroup("child").hasView("view2")

    # Delete source hierarchy - copied pins should keep arrays alive
    del src
    del child
    del ds
    del root
    del dst_parent
    _force_gc()

    # Both pins should still be valid
    assert ref1() is not None
    assert ref2() is not None
    np.testing.assert_array_equal(dst.getView("view1").getDataArray(), np.arange(5, dtype=np.int32))
    np.testing.assert_array_equal(
        dst.getGroup("child").getView("view2").getDataArray(), np.arange(8, dtype=np.int64))


def test_move_view_with_external_data_preserves_pin():
    """moveView should preserve the pin since the View* pointer doesn't change."""
    ds = pysidre.DataStore()
    root = ds.getRoot()
    src = root.createGroup("src")
    dst = root.createGroup("dst")

    # Create view with external data
    external = np.arange(12, dtype=np.int64)
    ref = weakref.ref(external)
    view = src.createView("moveable", external)
    view.apply(pysidre.TypeID.INT64_ID, 12)
    del external
    _force_gc()
    assert ref() is not None  # Pin keeps it alive

    # Move the view to dst - View* stays the same, pin should remain valid
    moved = dst.moveView(view)
    assert moved.getPathName() == "dst/moveable"
    assert moved.isExternal()
    assert not src.hasView("moveable")
    assert dst.hasView("moveable")

    # Delete original references - moved view should keep pin alive
    del view
    del src
    del ds
    del root
    del dst
    _force_gc()

    # Pin should still be valid
    assert ref() is not None
    np.testing.assert_array_equal(moved.getDataArray(), np.arange(12))


def test_destroy_view_by_index_releases_external_pin():
    """destroyView(IndexType) should release the external data pin."""
    ds = pysidre.DataStore()
    root = ds.getRoot()

    external = np.arange(7, dtype=np.int32)
    ref = weakref.ref(external)
    view = root.createView("indexed", external)
    view.apply(pysidre.TypeID.INT32_ID, 7)
    view_idx = view.getIndex()
    del external
    del view
    _force_gc()
    assert ref() is not None  # Pin keeps it alive

    # Destroy by index - should release the pin
    root.destroyView(view_idx)
    _force_gc()
    assert ref() is None  # Pin released, array collected


def test_pin_overwrite_warning():
    """Setting external data twice on the same View correctly replaces the pin."""
    ds = pysidre.DataStore()
    view = ds.getRoot().createView("test")

    # First external data
    external1 = np.arange(5, dtype=np.int32)
    view.setExternalData(pysidre.TypeID.INT32_ID, 5, external1)

    # Second external data on same view - old pin released, new pin created
    external2 = np.arange(10, dtype=np.int64)
    view.setExternalData(pysidre.TypeID.INT64_ID, 10, external2)

    # The second pin should be active; first pin was automatically released
    np.testing.assert_array_equal(view.getDataArray(), external2)


def test_registry_cleanup_on_explicit_destroy():
    """Pins are released when Views are explicitly destroyed, preventing registry bloat."""
    ds = pysidre.DataStore()
    root = ds.getRoot()

    weak_refs = []
    for i in range(5):
        external = np.arange(10, dtype=np.int32)
        weak_refs.append(weakref.ref(external))
        view = root.createView(f"view_{i}")
        view.setExternalData(pysidre.TypeID.INT32_ID, 10, external)
        del external  # Pin keeps it alive
        del view  # Don't hold view reference

    # Explicitly destroy all views - this releases pins
    for i in range(5):
        root.destroyView(f"view_{i}")

    _force_gc()

    # Verify all arrays were collected after explicit destruction
    collected_count = sum(1 for ref in weak_refs if ref() is None)
    assert collected_count == len(weak_refs), \
        f"Only {collected_count}/{len(weak_refs)} arrays collected after explicit destroy"


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
