# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Lifetime-soundness regression tests for the sidre python bindings.

Each test obtains a sidre-owned object (child or ancestor proxy, iterator, zero-copy numpy array),
drops every Python owner, forces a garbage collection, and then uses the object.
"""

import gc
import weakref

import numpy as np
import pytest

import axom.sidre as pysidre


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


def test_external_pins_released_when_datastore_destroyed():
    """Pins are released when the DataStore is destroyed without explicit destroy*().

    This is the implicit counterpart to test_registry_cleanup_on_explicit_destroy:
    the Views are torn down by the DataStore destructor (the C++ path), not by a
    bound destroyView()/destroyGroup(). A weak reference on the DataStore must
    still clear its pins, so the external arrays are collected and the registry
    does not accumulate dangling entries.
    """
    weak_refs = []

    def build_and_drop():
        ds = pysidre.DataStore()
        root = ds.getRoot()
        for i in range(5):
            external = np.arange(10, dtype=np.int32)
            weak_refs.append(weakref.ref(external))
            # Mix createView(external) and setExternalData() entry points.
            if i % 2 == 0:
                root.createView(f"view_{i}", external).apply(pysidre.TypeID.INT32_ID, 10)
            else:
                root.createView(f"view_{i}").setExternalData(pysidre.TypeID.INT32_ID, 10, external)
        # Pins keep the arrays alive while ds is alive...
        gc.collect()
        assert all(ref() is not None for ref in weak_refs)
        # ...and ds goes out of scope here without any explicit destroy call.

    build_and_drop()
    _force_gc()

    collected = sum(1 for ref in weak_refs if ref() is None)
    assert collected == len(weak_refs), \
        f"Only {collected}/{len(weak_refs)} external arrays collected after DataStore destruction"


def test_external_pins_released_for_nested_groups_on_datastore_destruction():
    """DataStore destruction releases pins for Views nested in child Groups too."""
    weak_refs = []

    def build_and_drop():
        ds = pysidre.DataStore()
        root = ds.getRoot()
        grp = root.createGroup("a/b/c")
        for i in range(3):
            external = np.arange(8, dtype=np.int64)
            weak_refs.append(weakref.ref(external))
            grp.createView(f"deep_{i}", external).apply(pysidre.TypeID.INT64_ID, 8)
        gc.collect()
        assert all(ref() is not None for ref in weak_refs)

    build_and_drop()
    _force_gc()

    assert all(ref() is None for ref in weak_refs), \
        "Nested-Group external arrays were not released on DataStore destruction"


def test_external_pins_isolated_between_datastores():
    """A View* address reused across DataStores must not cross-associate pins.

    Each DataStore owns a private pin scope. Destroying one DataStore releases
    only its own pins; a concurrently live DataStore is unaffected, even though
    the allocator may hand out overlapping View* addresses across them.
    """
    keep_alive = []
    surviving_refs = []

    # Build and drop several DataStores in sequence, encouraging View* reuse.
    for _ in range(4):
        ds = pysidre.DataStore()
        a = np.arange(6, dtype=np.int64)
        r = weakref.ref(a)
        ds.getRoot().createView("v", a).apply(pysidre.TypeID.INT64_ID, 6)
        del a, ds
        _force_gc()
        # Each dropped DataStore must release its own array.
        assert r() is None

    # A long-lived DataStore created afterwards (possibly at a reused address)
    # must hold its own pin independently.
    survivor = pysidre.DataStore()
    b = np.arange(6, dtype=np.int64)
    surviving_refs.append(weakref.ref(b))
    survivor.getRoot().createView("v", b).apply(pysidre.TypeID.INT64_ID, 6)
    keep_alive.append(survivor)
    del b
    _force_gc()
    assert surviving_refs[0]() is not None, \
        "Survivor DataStore's pin was wrongly released (cross-datastore misattribution)"

    # Cleanup releases the survivor's pin.
    del survivor
    keep_alive.clear()
    _force_gc()
    assert surviving_refs[0]() is None


def test_multiple_concurrent_datastores():
    """Multiple active DataStores with external data should not interfere with each other."""
    # Create multiple DataStores simultaneously
    datastores = []
    views = []
    arrays = []

    for ds_idx in range(3):
        ds = pysidre.DataStore()
        root = ds.getRoot()
        datastores.append(ds)

        # Each DataStore has multiple Views with external data
        ds_views = []
        ds_arrays = []
        for view_idx in range(3):
            external = np.arange(view_idx * 10, (view_idx + 1) * 10, dtype=np.int32)
            ds_arrays.append(external)
            view = root.createView(f"view_{view_idx}")
            view.setExternalData(pysidre.TypeID.INT32_ID, 10, external)
            ds_views.append(view)

        views.append(ds_views)
        arrays.append(ds_arrays)

    # Verify all views can access their data correctly
    for ds_idx in range(3):
        for view_idx in range(3):
            view = views[ds_idx][view_idx]
            expected = arrays[ds_idx][view_idx]
            retrieved = view.getDataArray()
            np.testing.assert_array_equal(retrieved,
                                          expected,
                                          err_msg=f"DS{ds_idx} view{view_idx} data mismatch")

    # Destroy one DataStore while others remain active
    root0 = datastores[0].getRoot()
    for view_idx in range(3):
        root0.destroyView(f"view_{view_idx}")

    # Keep references to DS1 and DS2 before deleting DS0
    ds1_views = views[1]
    ds2_views = views[2]
    ds1_arrays = arrays[1]
    ds2_arrays = arrays[2]

    del datastores[0], views[0], arrays[0]
    _force_gc()

    # Verify remaining DataStores still work correctly
    for local_idx, (ds_views, ds_arrays, ds_label) in enumerate([(ds1_views, ds1_arrays, "DS1"),
                                                                 (ds2_views, ds2_arrays, "DS2")]):
        for view_idx in range(3):
            view = ds_views[view_idx]
            expected = ds_arrays[view_idx]
            retrieved = view.getDataArray()
            np.testing.assert_array_equal(
                retrieved,
                expected,
                err_msg=f"After DS0 destroy: {ds_label} view{view_idx} data mismatch")

    # Create a new DataStore and verify it doesn't conflict
    ds_new = pysidre.DataStore()
    root_new = ds_new.getRoot()
    external_new = np.arange(100, 110, dtype=np.int32)
    view_new = root_new.createView("new_view")
    view_new.setExternalData(pysidre.TypeID.INT32_ID, 10, external_new)

    # Verify new DataStore works
    np.testing.assert_array_equal(view_new.getDataArray(), external_new)

    # Verify old DataStores still work
    for ds_views, ds_arrays, ds_label in [(ds1_views, ds1_arrays, "DS1"),
                                          (ds2_views, ds2_arrays, "DS2")]:
        for view_idx in range(3):
            view = ds_views[view_idx]
            expected = ds_arrays[view_idx]
            retrieved = view.getDataArray()
            np.testing.assert_array_equal(
                retrieved,
                expected,
                err_msg=f"After new DS: {ds_label} view{view_idx} data mismatch")


def test_concurrent_datastores_with_copy_move():
    """Copy/move operations should work correctly with multiple concurrent DataStores."""
    ds1 = pysidre.DataStore()
    ds2 = pysidre.DataStore()

    root1 = ds1.getRoot()
    root2 = ds2.getRoot()

    # Create view with external data in DS1
    external1 = np.arange(10, dtype=np.int32)
    view1 = root1.createView("src")
    view1.setExternalData(pysidre.TypeID.INT32_ID, 10, external1)

    # Copy view to a group in DS1
    grp1 = root1.createGroup("grp1")
    copied = grp1.copyView(view1)
    np.testing.assert_array_equal(copied.getDataArray(), external1)

    # Create view with different external data in DS2
    external2 = np.arange(20, 30, dtype=np.int64)
    view2 = root2.createView("other")
    view2.setExternalData(pysidre.TypeID.INT64_ID, 10, external2)

    # Verify both DataStores maintain correct data
    np.testing.assert_array_equal(view1.getDataArray(), external1)
    np.testing.assert_array_equal(copied.getDataArray(), external1)
    np.testing.assert_array_equal(view2.getDataArray(), external2)

    # Move view within DS1
    grp2 = root1.createGroup("grp2")
    moved = grp2.moveView(copied)
    np.testing.assert_array_equal(moved.getDataArray(), external1)

    # Verify DS2 unaffected
    np.testing.assert_array_equal(view2.getDataArray(), external2)


def test_concurrent_datastores_registry_isolation():
    """Registry should correctly isolate pins between different DataStores."""
    ds1 = pysidre.DataStore()
    ds2 = pysidre.DataStore()

    external1 = np.arange(5, dtype=np.int32)
    external2 = np.arange(5, dtype=np.int64)

    ref1 = weakref.ref(external1)
    ref2 = weakref.ref(external2)

    # Both DataStores use external data
    view1 = ds1.getRoot().createView("v1")
    view1.setExternalData(pysidre.TypeID.INT32_ID, 5, external1)

    view2 = ds2.getRoot().createView("v2")
    view2.setExternalData(pysidre.TypeID.INT64_ID, 5, external2)

    del external1, external2  # Only pins keep them alive
    _force_gc()

    # Both should still be alive
    assert ref1() is not None, "DS1 external data collected prematurely"
    assert ref2() is not None, "DS2 external data collected prematurely"

    # Destroy view in DS1
    ds1.getRoot().destroyView("v1")
    _force_gc()

    # DS1's external data should be collected, DS2's should not
    assert ref1() is None, "DS1 external data not collected after destroyView"
    assert ref2() is not None, "DS2 external data incorrectly collected"

    # Destroy view in DS2
    ds2.getRoot().destroyView("v2")
    _force_gc()

    # Now DS2's external data should be collected too
    assert ref2() is None, "DS2 external data not collected after destroyView"


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
