# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""
Convert a Sidre datastore from the sidre_hdf5 protocol to another protocol.

Users must supply a path to a sidre_hdf5 rootfile and base name for
the output datastores. Optional command line arguments include
a ``--protocol`` option (the default is ``json``)
and a ``--strip`` option to truncate the array data to at most N elements.
The strip option also prepends each array with its original size, the new
size and a filler entry of 0 for integer arrays or nan for floating point
arrays. E.g. if the array had 6 entries [1.01, 2.02, 3.03, 4.04, 5.05, 6.06]
and the user passed in ``--strip 3``, the array would be converted to
[6, 3, nan, 1.01, 2.02, 3.03].

The strip option is intended as a temporary solution to truncating
a dataset to allow easier debugging. In the future, the conversion and
truncation/display functionality may be separated into distinct utilities.
"""

from __future__ import annotations

import sys
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


#
# Allocate storage for external data of the input datastore.
#
# Iterates recursively through the views and groups of the provided group to
# find the external data views and allocates the required storage within the
# holders list.
#
# Also initializes the data in each allocated array to zeros.
#
def allocate_external_data(group: pysidre.Group, holders: list[np.ndarray], verbose: bool) -> None:
    # for each view
    for view in iter_views(group):
        if view.isExternal():
            if verbose:
                print(
                    f"Allocating external storage for {view.getPathName()} "
                    f"({view.getNumElements()} elements, {view.getTotalBytes()} bytes)", )
            storage = np.zeros(view.getTotalBytes(), dtype=np.uint8)
            view.setExternalData(storage)
            holders.append(storage)

    # for each group
    for child in iter_groups(group):
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
def truncate_bulk_data(group: pysidre.Group, max_size: int, verbose: bool) -> None:
    # for each view
    for view in iter_views(group):
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
    for child in iter_groups(group):
        truncate_bulk_data(child, max_size, verbose)


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

    print(f"Loading datastore from {input_path}")
    manager.read(root, str(input_path))
    num_files = manager.getNumFilesFromRoot(str(input_path))

    print("Loading external data from datastore")
    external_holders: list[np.ndarray] = []
    allocate_external_data(root, external_holders, args.verbose)
    manager.loadExternalData(root, str(input_path))

    if args.strip is not None:
        print(f"Truncating views to at most {args.strip} elements.")
        truncate_bulk_data(root, args.strip, args.verbose)
        note = ("This datastore was created by axom's 'convert_sidre_protocol' utility "
                f"with option '--strip {args.strip}'. To simplify debugging, the bulk "
                f"data in this datastore has been truncated to have at most {args.strip} "
                "original values per array. Three values have been prepended to each "
                "array: the size of the original array, the number of retained elements "
                "and a zero/Nan.")
        root.createViewString("Note", note)

    print(
        f"Writing out datastore in '{args.protocol}' protocol to file(s) with base name {args.output}",
    )
    manager.write(root, num_files, args.output, args.protocol)

    if initialized_mpi and not MPI.Is_finalized():
        MPI.Finalize()

    return 0


if __name__ == "__main__":
    sys.exit(main())
