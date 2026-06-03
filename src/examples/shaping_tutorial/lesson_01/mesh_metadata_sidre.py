#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

from argparse import ArgumentParser

# An example that uses Sidre to describe a 2D Cartesian mesh.
#
# This example demonstrates creating a Sidre DataStore and using Groups and Views
# to represent mesh metadata including bounding box coordinates and resolution.
# It also optionally generates the mesh following the conduit mesh blueprint convention.
#
# Example run:
#     ./mesh_metadata_sidre --min_x 0.0 --min_y 0.0 --max_x 2.0 --max_y 3.0 --res_x 20 --res_y 30 [-o]

try:
    import conduit
    import conduit.blueprint
    import conduit.relay
except ModuleNotFoundError as exc:
    raise RuntimeError("This lesson requires the Conduit Python modules. Run it through "
                       "Axom's build/bin wrapper so PYTHONPATH includes Conduit.") from exc

import pysidre


class Input:
    """Class representing user input parameters for mesh bounding box and resolution.

    This struct holds the minimum and maximum x and y coordinates defining the bounding box,
    as well as the resolution (number of cells) in both x and y directions.
    """

    def __init__(self):
        self.min_x = 0.0
        self.min_y = 0.0
        self.max_x = 1.0
        self.max_y = 1.0
        self.res_x = 10
        self.res_y = 20


def setup_mesh_metadata(datastore, input_data):
    """Set up the Sidre hierarchy for mesh metadata.

    This creates a hierarchical structure within the Sidre DataStore to represent
    mesh metadata, including bounding box coordinates (min and max) and resolution in both
    x and y directions. It stores these parameters in Groups and Views for easy access and
    manipulation.

    Args:
        datastore: Reference to the Sidre DataStore to populate.
        input_data: Bounding box coordinates and resolution data.
    """

    # Create a root group for the mesh metadata
    mesh_group = datastore.getRoot().createGroup("mesh")

    # Create bounding box groups and views
    min_group = mesh_group.createGroup("bounding_box/min")
    min_group.createViewScalar("x", input_data.min_x)
    min_group.createViewScalar("y", input_data.min_y)

    max_group = mesh_group.createGroup("bounding_box/max")
    max_group.createViewScalar("x", input_data.max_x)
    max_group.createViewScalar("y", input_data.max_y)

    # Create resolution group and views
    res_group = mesh_group.createGroup("resolution")
    res_group.createViewScalar("x", input_data.res_x)
    res_group.createViewScalar("y", input_data.res_y)


def verify_mesh_metadata(mesh_group):
    """Verify that the Sidre mesh metadata hierarchy is complete.

    This function accesses the mesh group and verifies that all required metadata
    (bounding box coordinates and resolution) exists.

    Args:
        mesh_group: Sidre Group representing the mesh metadata.

    Returns:
        True if the expected metadata is present, False otherwise.
    """

    if not mesh_group:
        print("Missing mesh group")
        return False

    valid = True

    # check bounding_box group
    if mesh_group.hasGroup("bounding_box"):
        bb_group = mesh_group.getGroup("bounding_box")
        for path in ["min/x", "min/y", "max/x", "max/y"]:
            if not bb_group.hasView(path):
                print(f"Missing '{path}' view in 'bounding_box' group")
                valid = False
    else:
        print("Missing 'bounding_box' group in mesh metadata")

    # check resolution group
    if mesh_group.hasGroup("resolution"):
        res_group = mesh_group.getGroup("resolution")
        for path in ["x", "y"]:
            if not res_group.hasView(path):
                print(f"Missing '{path}' view in 'resolution' group")
                valid = False
    else:
        print("Missing 'resolution' group in mesh metadata")
        return False

    return valid


def print_mesh_metadata(mesh_group):
    """Prints mesh metadata stored in Sidre DataStore

    This function accesses the mesh group and prints bounding box coordinates and resolution.

    Args:
        mesh_group: Sidre Group representing the mesh metadata.
    """

    if not mesh_group:
        return

    bb_group = mesh_group.getGroup("bounding_box")
    if not bb_group:
        return

    res_group = mesh_group.getGroup("resolution")
    if not res_group:
        return

    print(f"Bounding Box Min: ({bb_group.getView('min/x').getDataFloat()}, "
          f"{bb_group.getView('min/y').getDataFloat()})")
    print(f"Bounding Box Max: ({bb_group.getView('max/x').getDataFloat()}, "
          f"{bb_group.getView('max/y').getDataFloat()})")
    print(f"Resolution: ({res_group.getView('x').getDataInt()}, "
          f"{res_group.getView('y').getDataInt()})")


def create_mesh_blueprint(mesh_group):
    """Convert Sidre mesh metadata to a Conduit mesh blueprint.

    This function extracts mesh metadata from Sidre groups and views and creates
    a Conduit Node that represents a uniform 2D Cartesian mesh following
    the Conduit mesh blueprint conventions.

    Args:
        mesh_group: Sidre Group containing the mesh metadata.

    Returns:
        A Conduit Node containing the mesh blueprint representation.
    """

    blueprint = conduit.Node()

    if not mesh_group:
        print("Invalid mesh group provided")
        return blueprint

    # Get bounding box information
    bb_group = mesh_group.getGroup("bounding_box")
    if not bb_group:
        print("Missing bounding_box group in mesh metadata")
        return blueprint

    x_min = bb_group.getView("min/x").getDataFloat()
    y_min = bb_group.getView("min/y").getDataFloat()
    x_max = bb_group.getView("max/x").getDataFloat()
    y_max = bb_group.getView("max/y").getDataFloat()

    # Get resolution information
    res_group = mesh_group.getGroup("resolution")
    if not res_group:
        print("Missing resolution group in mesh metadata")
        return blueprint

    res_x = res_group.getView("x").getDataInt()
    res_y = res_group.getView("y").getDataInt()
    if res_x < 1:
        raise RuntimeError(f"Resolution in x-coordinate ({res_x}) must be positive")
    if res_y < 1:
        raise RuntimeError(f"Resolution in y-coordinate ({res_y}) must be positive")

    # Create coordset
    blueprint["coordsets/coords/type"] = "uniform"
    blueprint["coordsets/coords/dims/i"] = res_x + 1
    blueprint["coordsets/coords/dims/j"] = res_y + 1
    blueprint["coordsets/coords/origin/x"] = x_min
    blueprint["coordsets/coords/origin/y"] = y_min
    blueprint["coordsets/coords/spacing/dx"] = (x_max - x_min) / res_x
    blueprint["coordsets/coords/spacing/dy"] = (y_max - y_min) / res_y

    # Create topology
    blueprint["topologies/mesh/type"] = "uniform"
    blueprint["topologies/mesh/coordset"] = "coords"

    return blueprint


def save_blueprint(blueprint, file_name, file_format="yaml"):
    """Save a Conduit mesh blueprint to a file.

    This function verifies that the provided Node conforms to the Conduit mesh blueprint
    specification and saves it to the specified file format if valid.

    Args:
        blueprint: Conduit Node containing the mesh blueprint to save.
        file_name: Base name for the output file, without extension.
        file_format: Output format to save, such as ``yaml`` or ``json``.

    Returns:
        True if the save succeeded, False otherwise.
    """

    if blueprint.number_of_children() > 0:
        print(f"Saving mesh blueprint to '{file_name}.{file_format}'")
        info = conduit.Node()
        if conduit.blueprint.mesh.verify(blueprint, info):
            conduit.relay.io.blueprint.save_mesh(blueprint, file_name, file_format)
            return True

        print("Blueprint verification failed")
        print(info.to_string())
        return False

    print("Unable to create mesh blueprint, cannot save file")
    return False


def main():
    input_data = Input()
    output_blueprint = False

    # parse input using argparse
    app = ArgumentParser(description="Mesh Metadata Setup")
    app.add_argument("--min_x", type=float, default=input_data.min_x)
    app.add_argument("--min_y", type=float, default=input_data.min_y)
    app.add_argument("--max_x", type=float, default=input_data.max_x)
    app.add_argument("--max_y", type=float, default=input_data.max_y)
    app.add_argument("--res_x", type=int, default=input_data.res_x)
    app.add_argument("--res_y", type=int, default=input_data.res_y)
    app.add_argument("-o", "--output_blueprint", action="store_true")
    args = app.parse_args()

    input_data.min_x = args.min_x
    input_data.min_y = args.min_y
    input_data.max_x = args.max_x
    input_data.max_y = args.max_y
    input_data.res_x = args.res_x
    input_data.res_y = args.res_y
    output_blueprint = args.output_blueprint

    # load parsed data into sidre datastore
    datastore = pysidre.DataStore()
    setup_mesh_metadata(datastore, input_data)

    # validate and print results
    root = datastore.getRoot()
    if not root.hasGroup("mesh"):
        raise RuntimeError("Missing expected 'mesh' group")

    mesh_group = root.getGroup("mesh")
    is_valid = verify_mesh_metadata(mesh_group)
    if is_valid:
        print("Sidre hierarchy was properly set up")
        print_mesh_metadata(mesh_group)
    else:
        print("Sidre hierarchy was not properly set up")
        print("Sidre hierarchy:")
        root.print()

    # Optionally, create mesh blueprint and save as a yaml file
    if output_blueprint:
        blueprint = create_mesh_blueprint(mesh_group)
        save_blueprint(blueprint, "uniform_bp")

    return 0 if is_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
