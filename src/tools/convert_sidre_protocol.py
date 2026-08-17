# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""
Convert a Sidre datastore or Conduit Blueprint mesh to another protocol.

Users must supply a path to an input rootfile and base name for the output.
Sidre datastore conversion remains the default. Use ``--input-type blueprint``
to convert Conduit Blueprint meshes, including HDF5 Blueprint root files.
Optional command line arguments include a ``--protocol`` option (the default is ``json``)
and a ``--strip`` option to truncate the array data to at most N elements.
The strip option also prepends each array with its original size, the new
size and a filler entry of 0 for integer arrays or nan for floating point
arrays. E.g. if the array had 6 entries [1.01, 2.02, 3.03, 4.04, 5.05, 6.06]
and the user passed in ``--strip 3``, the array would be converted to
[6, 3, nan, 1.01, 2.02, 3.03].

For Blueprint meshes, ``--strip`` applies to numeric array leaves with more than
one element and leaves scalar/string metadata intact.
A ``state/Note`` node is added to each output domain to record that the mesh was stripped.

The resulting stripped output is intended for debugging, not for use as a valid mesh.
In the future, the conversion and truncation/display functionality may be
separated into distinct utilities.
"""

from __future__ import annotations

import sys
import argparse
import math
from pathlib import Path

import numpy as np
import axom.sidre as sidre

SIDRE_PROTOCOLS = (
    "json",
    "sidre_hdf5",
    "sidre_conduit_json",
    "sidre_json",
    "conduit_hdf5",
    "conduit_bin",
    "conduit_json",
)

BLUEPRINT_PROTOCOLS = (
    "hdf5",
    "json",
    "yaml",
    "conduit_bin",
    "conduit_hdf5",
    "conduit_json",
)

VALID_PROTOCOLS = tuple(dict.fromkeys(SIDRE_PROTOCOLS + BLUEPRINT_PROTOCOLS))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sidre/Blueprint protocol converter")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Filename of input sidre-hdf5 datastore or Blueprint mesh root file",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Filename of output datastore/mesh (without extension)",
    )
    parser.add_argument(
        "--input-type",
        choices=("sidre", "blueprint"),
        default="sidre",
        help="Type of input file to convert",
    )
    parser.add_argument(
        "-p",
        "--protocol",
        default="json",
        choices=VALID_PROTOCOLS,
        help="Desired output protocol; valid protocols depend on --input-type",
    )
    parser.add_argument(
        "-s",
        "--strip",
        type=int,
        default=None,
        help="If provided, output arrays will be stripped to first N entries",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Sets output to verbose",
    )
    return parser.parse_args()


def initialize_mpi() -> tuple[object | None, bool, int, int]:
    if not sidre.AXOM_ENABLE_MPI:
        return None, False, 1, 0

    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError(
            "convert_sidre_protocol.py requires mpi4py when Axom is built with MPI support",
        ) from exc

    initialized_mpi = False
    if not MPI.Is_initialized():
        MPI.Init()
        initialized_mpi = True

    return MPI, initialized_mpi, MPI.COMM_WORLD.Get_size(), MPI.COMM_WORLD.Get_rank()


#
# Allocate storage for external data of the input datastore.
#
# Iterates recursively through the views and groups of the provided group to
# find the external data views and allocates the required storage within the
# holders list.
#
# Also initializes the data in each allocated array to zeros.
#
def allocate_external_data(group: sidre.Group, holders: list[np.ndarray], verbose: bool) -> None:
    # for each view
    for view in group.views():
        if view.isExternal():
            if verbose:
                print(
                    f"Allocating external storage for {view.getPathName()} "
                    f"({view.getNumElements()} elements, {view.getTotalBytes()} bytes)", )
            storage = np.zeros(view.getTotalBytes(), dtype=np.uint8)
            view.setExternalData(storage)
            holders.append(storage)

    # for each group
    for child in group.groups():
        allocate_external_data(child, holders, verbose)


#
# Shift the data to the right by three elements.
#
# The new first value will be the size of the original array.
# The next value will be the number of retained elements and the
# third value will be 0 for integer data and Nan for float data.
# This is followed by the values in the truncated original dataset.
#
# This function creates a copy of the data since there could be
# several views in the original dataset pointing to the same memory.
#
def modify_final_values(
    view: sidre.View,
    original_size: int,
    retained_size: int | None = None,
) -> None:
    flattened = np.asarray(view.getDataArray()).reshape(-1)
    if retained_size is None:
        retained = flattened
    else:
        retained = flattened[:retained_size]

    retained_size = int(retained.size)

    # Create a new buffer for copied data.
    datastore = view.getOwningGroup().getDataStore()
    type_id = view.getTypeID()
    new_size = retained_size + 3
    buffer = datastore.createBuffer(type_id, new_size)
    buffer.allocate()

    # Explicitly set the first two elements and copy elements over.
    new_array = np.asarray(buffer.getDataArray()).reshape(-1)
    np.copyto(new_array[:2], np.asarray([original_size, retained_size]), casting="unsafe")
    new_array[2] = math.nan if np.issubdtype(new_array.dtype, np.floating) else 0
    if retained_size > 0:
        np.copyto(new_array[3:], retained, casting="unsafe")

    # Update view's buffer to the new data.
    # The C++ utility uses detachBuffer() here because this path is valid for
    # buffer-backed views. In Python we also need to support external views, so
    # we make the state transition explicit before attaching the new buffer.
    if view.hasBuffer():
        view.attachBuffer(None)
    elif view.isExternal():
        view.setExternalData(None)

    view.attachBuffer(type_id, new_size, buffer)


#
# Recursively traverse views and groups in group and truncate views to
# have at most max_size elements.
#
# Within the truncated arrays, the first element will be the size of the
# original array, the second will be the number of retained elements and the
# third will be 0 for integers or nan for floating points.
# This will be followed by at most the first max_size elements of the
# original array.
#
def truncate_bulk_data(group: sidre.Group, max_size: int, verbose: bool) -> None:
    # for each view
    for view in group.views():
        is_array = view.hasBuffer() or view.isExternal()

        if is_array:
            original_size = view.getNumElements()
            retained_size = min(max_size, original_size)

            if view.hasBuffer() and original_size > retained_size:
                view.apply(retained_size, view.getOffset(), view.getStride())
            elif view.isExternal() and original_size > retained_size:
                data = np.asarray(view.getDataArray()).reshape(-1)
                view.setExternalData(view.getTypeID(), retained_size, data)

            if verbose:
                print(
                    f"Truncating view {view.getPathName()} from {original_size} to {retained_size}",
                )
            modify_final_values(view, original_size, retained_size)

    # for each group
    for child in group.groups():
        truncate_bulk_data(child, max_size, verbose)


def blueprint_protocol(protocol: str) -> str:
    if protocol == "conduit_hdf5":
        return "hdf5"
    if protocol in BLUEPRINT_PROTOCOLS:
        return protocol

    valid = ", ".join(BLUEPRINT_PROTOCOLS)
    raise RuntimeError(
        f"Protocol '{protocol}' is not valid for Blueprint mesh output. "
        f"Use one of: {valid}", )


def sidre_protocol(protocol: str) -> str:
    if protocol in SIDRE_PROTOCOLS:
        return protocol

    valid = ", ".join(SIDRE_PROTOCOLS)
    raise RuntimeError(
        f"Protocol '{protocol}' is not valid for Sidre datastore output. "
        f"Use one of: {valid}", )


def blueprint_domain_count(mesh) -> int:
    if mesh.has_path("coordsets") and mesh.has_path("topologies"):
        return 1
    return mesh.number_of_children()


def strip_note(data_kind: str, max_size: int) -> str:
    return (f"This {data_kind} was created by axom's 'convert_sidre_protocol' utility "
            f"with option '--strip {max_size}'. To simplify debugging, the bulk "
            f"data in this {data_kind} has been truncated to have at most {max_size} "
            "original values per array. Three values have been prepended to each array: "
            "the size of the original array, the number of retained elements and a zero/Nan.")


def add_blueprint_strip_note(mesh, note: str) -> None:
    if mesh.has_path("coordsets") and mesh.has_path("topologies"):
        mesh["state/Note"] = note
        return

    for child_idx in range(mesh.number_of_children()):
        domain = mesh.child(child_idx)
        if domain.has_path("coordsets") and domain.has_path("topologies"):
            domain["state/Note"] = note


def truncate_conduit_numeric_array(node, max_size: int, verbose: bool) -> int:
    dtype = node.dtype()
    original_size = dtype.number_of_elements()
    if not dtype.is_number() or original_size <= 1:
        return 0

    retained_size = min(max_size, original_size)
    values = np.asarray(node.value()).reshape(-1)
    retained = values[:retained_size].copy()

    new_values = np.empty(retained_size + 3, dtype=values.dtype)
    np.copyto(new_values[:2], np.asarray([original_size, retained_size]), casting="unsafe")
    new_values[2] = math.nan if np.issubdtype(new_values.dtype, np.floating) else 0
    if retained_size > 0:
        np.copyto(new_values[3:], retained, casting="unsafe")

    if verbose:
        print(f"Truncating node {node.path()} from {original_size} to {retained_size}")

    node.reset()
    node.set(new_values)
    return 1


def truncate_conduit_bulk_data(node, max_size: int, verbose: bool) -> int:
    if node.number_of_children() == 0:
        return truncate_conduit_numeric_array(node, max_size, verbose)

    truncated_count = 0
    for child_idx in range(node.number_of_children()):
        truncated_count += truncate_conduit_bulk_data(node.child(child_idx), max_size, verbose)
    return truncated_count


def convert_blueprint_mesh(args: argparse.Namespace, MPI: object | None, comm_size: int,
                           rank: int) -> int:
    import conduit
    import conduit.blueprint
    import conduit.relay.io.blueprint

    input_path = Path(args.input)
    protocol = blueprint_protocol(args.protocol)
    mesh = conduit.Node()

    if MPI is not None and comm_size > 1:
        import conduit.blueprint.mpi.mesh
        import conduit.relay.mpi.io.blueprint

        comm = MPI.COMM_WORLD.py2f()
        if rank == 0:
            print(f"Loading Blueprint mesh from {input_path} on {comm_size} MPI rank(s)", )
        conduit.relay.mpi.io.blueprint.load_mesh(mesh, str(input_path), comm)

        info = conduit.Node()
        valid = conduit.blueprint.mpi.mesh.verify(mesh, info, comm)
        local_domains = blueprint_domain_count(mesh)
        total_domains = MPI.COMM_WORLD.allreduce(local_domains)
        if rank == 0:
            print(f"Input Blueprint mesh layout: {total_domains} domain(s)")
        if not valid:
            raise RuntimeError(f"Input Blueprint mesh failed verification:\n{info.to_yaml()}")

        if args.strip is not None:
            if rank == 0:
                print(f"Truncating numeric Blueprint arrays to at most {args.strip} elements.")
            local_truncated = truncate_conduit_bulk_data(mesh, args.strip, args.verbose)
            total_truncated = MPI.COMM_WORLD.allreduce(local_truncated)
            if rank == 0:
                print(f"Truncated {total_truncated} numeric Blueprint array(s).")
            add_blueprint_strip_note(mesh, strip_note("Blueprint mesh", args.strip))

        if rank == 0:
            print(
                f"Writing out Blueprint mesh in '{protocol}' protocol to file(s) "
                f"with base name {args.output}", )
        conduit.relay.mpi.io.blueprint.save_mesh(mesh, args.output, comm, protocol)
    else:
        print(f"Loading Blueprint mesh from {input_path}")
        conduit.relay.io.blueprint.load_mesh(mesh, str(input_path))

        info = conduit.Node()
        valid = conduit.blueprint.mesh.verify(mesh, info)
        domains = blueprint_domain_count(mesh)
        print(f"Input Blueprint mesh layout: {domains} domain(s)")
        if not valid:
            raise RuntimeError(f"Input Blueprint mesh failed verification:\n{info.to_yaml()}")

        if args.strip is not None:
            print(f"Truncating numeric Blueprint arrays to at most {args.strip} elements.")
            truncated_count = truncate_conduit_bulk_data(mesh, args.strip, args.verbose)
            print(f"Truncated {truncated_count} numeric Blueprint array(s).")
            add_blueprint_strip_note(mesh, strip_note("Blueprint mesh", args.strip))

        print(
            f"Writing out Blueprint mesh in '{protocol}' protocol to file(s) "
            f"with base name {args.output}", )
        conduit.relay.io.blueprint.save_mesh(mesh, args.output, protocol)

    return 0


def convert_sidre_datastore(args: argparse.Namespace, comm_size: int) -> int:
    if not sidre.AXOM_ENABLE_MPI:
        raise RuntimeError("sidre.IOManager bindings require an MPI-enabled Axom build")

    protocol = sidre_protocol(args.protocol)
    input_path = Path(args.input)
    manager = sidre.IOManager()
    datastore = sidre.DataStore()
    root = datastore.getRoot()

    num_files = manager.getNumFilesFromRoot(str(input_path))
    num_groups = manager.getNumGroupsFromRoot(str(input_path))

    print(
        "Input datastore layout: "
        f"{num_groups} rank group(s) across {num_files} file(s); "
        f"running on {comm_size} MPI rank(s)", )

    if comm_size != num_groups:
        print(
            "Warning: current MPI size does not match the input datastore rank count. "
            "sidre_hdf5 supports some rank-count mismatches, but not every layout.", )
        if comm_size < num_groups and num_files != num_groups:
            print(
                "Warning: this run has fewer MPI ranks than the input datastore. "
                "Sidre only supports that case for file-per-rank sidre_hdf5 datasets "
                f"(number_of_files == number_of_trees, but here {num_files} != {num_groups}), "
                "so IOManager.read() is expected to fail.", )
        elif comm_size < num_groups:
            print(
                "Warning: this run has fewer MPI ranks than the input datastore, but the "
                "dataset is file-per-rank "
                f"({num_files} file(s) for {num_groups} rank group(s)), so Sidre can load it. "
                "In that reduced-rank case, one output rank may absorb data from multiple input "
                "ranks, and the loaded hierarchy may be reshaped under "
                "'rank_%07d/sidre_input' groups.", )

    print(f"Loading datastore from {input_path}")
    manager.read(root, str(input_path))

    print("Loading external data from datastore")
    external_holders: list[np.ndarray] = []
    allocate_external_data(root, external_holders, args.verbose)
    manager.loadExternalData(root, str(input_path))

    if args.strip is not None:
        print(f"Truncating views to at most {args.strip} elements.")
        truncate_bulk_data(root, args.strip, args.verbose)
        root.createViewString("Note", strip_note("datastore", args.strip))

    print(
        f"Writing out datastore in '{protocol}' protocol to file(s) with base name {args.output}", )
    manager.write(root, num_files, args.output, protocol)

    return 0


def main() -> int:
    args = parse_args()
    if args.strip is not None and args.strip < 0:
        raise RuntimeError("--strip must be nonnegative")

    MPI, initialized_mpi, comm_size, rank = initialize_mpi()
    try:
        if args.input_type == "blueprint":
            return convert_blueprint_mesh(args, MPI, comm_size, rank)
        return convert_sidre_datastore(args, comm_size)
    finally:
        if MPI is not None and initialized_mpi and not MPI.Is_finalized():
            MPI.Finalize()


if __name__ == "__main__":
    sys.exit(main())
