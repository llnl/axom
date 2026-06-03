#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path

try:
    import conduit
    import conduit.blueprint
    import conduit.relay
except ModuleNotFoundError as exc:
    raise RuntimeError("This lesson requires the Conduit Python modules. Run it through "
                       "Axom's build/bin wrapper so PYTHONPATH includes Conduit.") from exc

import pysidre


def setup_mesh_metadata(root, args):
    mesh_group = root.createGroup("mesh")
    min_group = mesh_group.createGroup("bounding_box/min")
    min_group.createViewScalar("x", args.min_x)
    min_group.createViewScalar("y", args.min_y)

    max_group = mesh_group.createGroup("bounding_box/max")
    max_group.createViewScalar("x", args.max_x)
    max_group.createViewScalar("y", args.max_y)

    res_group = mesh_group.createGroup("resolution")
    res_group.createViewScalar("x", args.res_x)
    res_group.createViewScalar("y", args.res_y)
    return mesh_group


def verify_mesh_metadata(mesh_group) -> bool:
    required_views = ["min/x", "min/y", "max/x", "max/y"]
    if not mesh_group.hasGroup("bounding_box"):
        return False
    if not mesh_group.hasGroup("resolution"):
        return False

    bb_group = mesh_group.getGroup("bounding_box")
    res_group = mesh_group.getGroup("resolution")
    return all(bb_group.hasView(path)
               for path in required_views) and all(res_group.hasView(path) for path in ["x", "y"])


def create_blueprint(mesh_group) -> conduit.Node:
    bb_group = mesh_group.getGroup("bounding_box")
    res_group = mesh_group.getGroup("resolution")
    x_min = bb_group.getView("min/x").getDataFloat()
    y_min = bb_group.getView("min/y").getDataFloat()
    x_max = bb_group.getView("max/x").getDataFloat()
    y_max = bb_group.getView("max/y").getDataFloat()
    res_x = res_group.getView("x").getDataInt()
    res_y = res_group.getView("y").getDataInt()

    blueprint = conduit.Node()
    blueprint["coordsets/coords/type"] = "uniform"
    blueprint["coordsets/coords/dims/i"] = res_x + 1
    blueprint["coordsets/coords/dims/j"] = res_y + 1
    blueprint["coordsets/coords/origin/x"] = x_min
    blueprint["coordsets/coords/origin/y"] = y_min
    blueprint["coordsets/coords/spacing/dx"] = (x_max - x_min) / res_x
    blueprint["coordsets/coords/spacing/dy"] = (y_max - y_min) / res_y
    blueprint["topologies/mesh/type"] = "uniform"
    blueprint["topologies/mesh/coordset"] = "coords"
    return blueprint


def main() -> int:
    parser = ArgumentParser(description="Set up and inspect mesh metadata in Sidre.")
    parser.add_argument("--min_x", type=float, default=0.0)
    parser.add_argument("--min_y", type=float, default=0.0)
    parser.add_argument("--max_x", type=float, default=1.0)
    parser.add_argument("--max_y", type=float, default=1.0)
    parser.add_argument("--res_x", type=int, default=10)
    parser.add_argument("--res_y", type=int, default=20)
    parser.add_argument("--output_blueprint", action="store_true")
    parser.add_argument(
        "--output_file",
        default="uniform_bp",
        help="Base filename for the blueprint rootfile output",
    )
    args = parser.parse_args()

    datastore = pysidre.DataStore()
    mesh_group = setup_mesh_metadata(datastore.getRoot(), args)
    if not verify_mesh_metadata(mesh_group):
        raise RuntimeError("Sidre hierarchy validation failed.")

    print("Bounding Box Min:", (args.min_x, args.min_y))
    print("Bounding Box Max:", (args.max_x, args.max_y))
    print("Resolution:", (args.res_x, args.res_y))

    if args.output_blueprint:
        output_base = Path(args.output_file).with_suffix("")
        blueprint = create_blueprint(mesh_group)
        info = conduit.Node()
        if not conduit.blueprint.mesh.verify(blueprint, info):
            raise RuntimeError(f"Blueprint verification failed: {info.to_string()}")

        conduit.relay.io.blueprint.save_mesh(blueprint, str(output_base), "yaml")
        print("Wrote blueprint:", output_base.with_suffix(".root"))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
