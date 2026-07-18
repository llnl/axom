# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""IOManager communicator-construction tests.

These exercise the IOManager constructor's optional mpi4py communicator argument:
the default (MPI_COMM_WORLD) path and an explicit split-communicator path.
They are skipped unless the bindings were built with MPI and mpi4py is importable.

Run multi-rank, e.g.:
    mpirun -n 2 python -m pytest sidre_spio_Py.py
"""

import os

import pytest

import axom.sidre as sidre

if not sidre.AXOM_ENABLE_MPI:
    pytest.skip("sidre built without MPI", allow_module_level=True)

mpi4py = pytest.importorskip("mpi4py")
from mpi4py import MPI  # noqa: E402


def _shared_base(tmp_path, name):
    world = MPI.COMM_WORLD
    base_dir = world.bcast(str(tmp_path) if world.Get_rank() == 0 else None, root=0)
    return os.path.join(base_dir, name)


def _fill_datastore():
    ds = sidre.DataStore()
    root = ds.getRoot()
    view = root.createViewAndAllocate("field", sidre.TypeID.INT_ID, 4)
    rank = MPI.COMM_WORLD.Get_rank()
    view.getDataArray()[:] = [rank, rank + 1, rank + 2, rank + 3]
    return ds


def test_iomanager_legacy_use_scr_constructor():
    # Preserve the previous IOManager(use_scr=False) positional API.
    sidre.IOManager(False)


def test_iomanager_rejects_non_communicator():
    with pytest.raises(AttributeError):
        sidre.IOManager(object())


def test_iomanager_default_communicator(tmp_path):
    # No communicator argument -> MPI_COMM_WORLD (preserves the prior behavior).
    world = MPI.COMM_WORLD
    ds = _fill_datastore()
    iom = sidre.IOManager()
    base = _shared_base(tmp_path, "default_comm")
    iom.write(ds.getRoot(), 1, base, sidre.Group.getDefaultIOProtocol())
    world.Barrier()
    assert iom.getNumFilesFromRoot(base + ".root") == 1
    assert iom.getNumGroupsFromRoot(base + ".root") == world.Get_size()


def test_iomanager_explicit_world_communicator(tmp_path):
    # Passing COMM_WORLD explicitly must match the default path.
    world = MPI.COMM_WORLD
    ds = _fill_datastore()
    iom = sidre.IOManager(MPI.COMM_WORLD)
    base = _shared_base(tmp_path, "explicit_world")
    iom.write(ds.getRoot(), 1, base, sidre.Group.getDefaultIOProtocol())
    world.Barrier()

    ds_in = sidre.DataStore()
    iom.read(ds_in.getRoot(), base + ".root")
    arr = ds_in.getRoot().getView("field").getDataArray()
    rank = world.Get_rank()
    assert list(arr) == [rank, rank + 1, rank + 2, rank + 3]


def test_iomanager_owned_duplicate_survives_comm_free(tmp_path):
    # IOManager duplicates the input communicator, so callers may free their
    # mpi4py communicator after construction.
    comm = MPI.COMM_SELF.Dup()
    iom = sidre.IOManager(comm)
    comm.Free()

    ds = _fill_datastore()
    base = _shared_base(tmp_path, f"freed_comm_rank{MPI.COMM_WORLD.Get_rank()}")
    iom.write(ds.getRoot(), 1, base, sidre.Group.getDefaultIOProtocol())
    assert iom.getNumFilesFromRoot(base + ".root") == 1
    assert iom.getNumGroupsFromRoot(base + ".root") == 1
    MPI.COMM_WORLD.Barrier()


def test_iomanager_split_communicator(tmp_path):
    # Construct from a split communicator; each writes its own dataset.
    world = MPI.COMM_WORLD
    if world.Get_size() < 2:
        pytest.skip("split-communicator test needs >= 2 ranks")

    color = world.Get_rank() % 2
    sub = world.Split(color=color, key=world.Get_rank())
    sub_freed = False
    try:
        sub_size = sub.Get_size()
        ds = _fill_datastore()
        iom = sidre.IOManager(sub)
        sub.Free()
        sub_freed = True
        # tmp_path differs per rank; rendezvous on a shared, rank-0-broadcast dir
        base = _shared_base(tmp_path, f"split_color{color}")
        iom.write(ds.getRoot(), 1, base, sidre.Group.getDefaultIOProtocol())
        assert iom.getNumFilesFromRoot(base + ".root") == 1
        assert iom.getNumGroupsFromRoot(base + ".root") == sub_size
    finally:
        if not sub_freed:
            sub.Free()


def test_distributed_generate_blueprint_index(tmp_path):
    # The distributed generateBlueprintIndex overload is built only under
    # nanobind >= 2.10; skip cleanly if this build omitted it.
    if not sidre.AXOM_HAS_DISTRIBUTED_BLUEPRINT_INDEX_BINDING:
        pytest.skip("generateBlueprintIndex not bound")

    ds = sidre.DataStore()
    root = ds.getRoot()
    mesh = root.createGroup("mesh")
    coords = mesh.createGroup("coordsets/coords")
    coords.createViewString("type", "explicit")
    topo = mesh.createGroup("topologies/topo")
    topo.createViewString("type", "unstructured")
    topo.createViewString("coordset", "coords")

    base_dir = MPI.COMM_WORLD.bcast(str(tmp_path) if MPI.COMM_WORLD.Get_rank() == 0 else None,
                                    root=0)
    out = os.path.join(base_dir, "bp_index")

    # Distributed signature: (comm, domain_path, mesh_name, index_path).
    # We only assert the call is reachable and returns a bool; full Blueprint
    # validity is covered by the C++ spio tests.
    result = ds.generateBlueprintIndex(MPI.COMM_WORLD, "mesh", "mesh", out)
    assert isinstance(result, bool)


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
