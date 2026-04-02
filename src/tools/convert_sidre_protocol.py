# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Convert a Sidre datastore to another supported protocol."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from mpi4py import MPI
import pysidre

VALID_PROTOCOLS = (
    "json",
    "sidre_hdf5",
    "sidre_conduit_json",
    "sidre_json",
    "conduit_hdf5",
    "conduit_bin",
    "conduit_json",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sidre protocol converter")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Filename of input sidre-hdf5 datastore",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Filename of output datastore (without extension)",
    )
    parser.add_argument(
        "-p",
        "--protocol",
        default="json",
        choices=VALID_PROTOCOLS,
        help="Desired protocol for output datastore",
    )
    parser.add_argument(
        "-s",
        "--strip",
        type=positive_int,
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


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("strip value must be positive")
    return parsed


def log(enabled: bool, message: str) -> None:
    if enabled:
        print(message)


def iter_views(group: pysidre.Group):
    idx = group.getFirstValidViewIndex()
    while pysidre.indexIsValid(idx):
        yield group.getView(idx)
        idx = group.getNextValidViewIndex(idx)


def iter_groups(group: pysidre.Group):
    idx = group.getFirstValidGroupIndex()
    while pysidre.indexIsValid(idx):
        yield group.getGroup(idx)
        idx = group.getNextValidGroupIndex(idx)


def allocate_external_data(group: pysidre.Group, holders: list[np.ndarray], verbose: bool) -> None:
    for view in iter_views(group):
        if view.isExternal():
            log(
                verbose,
                f"Allocating external storage for {view.getPathName()} "
                f"({view.getNumElements()} elements, {view.getTotalBytes()} bytes)",
            )
            storage = np.zeros(view.getTotalBytes(), dtype=np.uint8)
            view.setExternalData(storage)
            holders.append(storage)

    for child in iter_groups(group):
        allocate_external_data(child, holders, verbose)


def filler_value(dtype: np.dtype):
    if np.issubdtype(dtype, np.floating):
        return math.nan
    return 0


def modify_final_values(
    view: pysidre.View,
    original_size: int,
    retained_size: int | None = None,
) -> None:
    flattened = np.asarray(view.getDataArray()).reshape(-1)
    if retained_size is None:
        retained = flattened
    else:
        retained = flattened[:retained_size]

    retained_size = int(retained.size)
    datastore = view.getOwningGroup().getDataStore()
    type_id = view.getTypeID()
    new_size = retained_size + 3
    buffer = datastore.createBuffer(type_id, new_size)
    buffer.allocate()

    new_array = np.asarray(buffer.getDataArray()).reshape(-1)
    np.copyto(new_array[:2], np.asarray([original_size, retained_size]), casting="unsafe")
    new_array[2] = filler_value(new_array.dtype)
    if retained_size > 0:
        np.copyto(new_array[3:], retained, casting="unsafe")

    view.replaceDataWithBuffer(type_id, new_size, buffer)


def truncate_bulk_data(group: pysidre.Group, max_size: int, verbose: bool) -> None:
    for view in iter_views(group):
        if not (view.hasBuffer() or view.isExternal()):
            continue

        original_size = view.getNumElements()
        retained_size = min(max_size, original_size)

        if view.hasBuffer() and original_size > retained_size:
            view.apply(retained_size, view.getOffset(), view.getStride())

        log(verbose,
            f"Truncating view {view.getPathName()} from {original_size} to {retained_size}")
        modify_final_values(view, original_size, retained_size)

    for child in iter_groups(group):
        truncate_bulk_data(child, max_size, verbose)


def add_strip_note(root: pysidre.Group, num_elements: int) -> None:
    note = ("This datastore was created by axom's 'convert_sidre_protocol' utility "
            f"with option '--strip {num_elements}'. To simplify debugging, the bulk "
            f"data in this datastore has been truncated to have at most {num_elements} "
            "original values per array. Three values have been prepended to each "
            "array: the size of the original array, the number of retained elements "
            "and a zero/Nan.")
    root.createViewString("Note", note)


def main() -> int:
    args = parse_args()

    if not pysidre.AXOM_ENABLE_MPI:
        raise RuntimeError("pysidre.IOManager bindings require an MPI-enabled Axom build")

    initialized_mpi = False
    if not MPI.Is_initialized():
        MPI.Init()
        initialized_mpi = True

    input_path = Path(args.input)
    manager = pysidre.IOManager()
    datastore = pysidre.DataStore()
    root = datastore.getRoot()

    log(args.verbose, f"Loading datastore from {input_path}")
    manager.read(root, str(input_path))
    num_files = manager.getNumFilesFromRoot(str(input_path))

    log(args.verbose, "Loading external data from datastore")
    external_holders: list[np.ndarray] = []
    allocate_external_data(root, external_holders, args.verbose)
    manager.loadExternalData(root, str(input_path))

    if args.strip is not None:
        log(args.verbose, f"Truncating views to at most {args.strip} elements")
        truncate_bulk_data(root, args.strip, args.verbose)
        add_strip_note(root, args.strip)

    log(
        args.verbose,
        f"Writing out datastore in '{args.protocol}' protocol to file(s) with base name {args.output}",
    )
    manager.write(root, num_files, args.output, args.protocol)

    if initialized_mpi and not MPI.Is_finalized():
        MPI.Finalize()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
