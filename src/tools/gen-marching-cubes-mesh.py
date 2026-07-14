#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

# Write a single-domain structured or unstructured Blueprint mesh for
# quest::MarchingCubes testing.
#
# The generated Blueprint hierarchy is:
#   <bp_root>
#    ├── state
#    │   └─• domain_id              == 0
#    ├── topologies
#    │   └── <topologyName>
#    │        ├─• coordset          == <coordsetName>
#    │        ├─• type              == "structured" or "unstructured"
#    │        └── elements
#    │             ├─• dims/{i,j,[k]}      (structured cell dimensions)
#    │             ├─• shape               (unstructured "quad" or "hex")
#    │             └─• connectivity        (unstructured flat int64 connectivity)
#    ├── coordsets
#    │   └── <coordsetName>
#    │        ├─• type              == "explicit"
#    │        └── values            (i-fastest node ordering)
#    │             ├─• x            (float64, node_count)
#    │             ├─• y            (float64, node_count)
#    │             └─• [z]          (float64, node_count, present in 3D)
#    └── fields
#        └── <fieldName>
#             ├─• topology          == <topologyName>
#             ├─• association       == "vertex"
#             └─• values            (float64 signed distance samples)

try:
    import conduit
    import conduit.blueprint
    import conduit.relay
except ModuleNotFoundError as e:
    print(f'{e}\nMake sure your PYTHONPATH includes /path/to/conduit/install/python-modules\n'
          'Conduit must be configured with python and hdf5.\n'
          'Alternatively, use the build directory convenience script:\n'
          '/path/to/axom_build_dir/bin/run_python_with_axom.sh')
    exit(-1)

import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def parse_component_list(values, cast):
    '''Convert comma-separated and/or space-separated values to a list.'''
    result = []
    for value in values:
        for component in value.split(','):
            component = component.strip()
            if component:
                result.append(cast(component))
    return result


def parse_args():
    ps = ArgumentParser(description='Write a single-domain MarchingCubes Blueprint mesh.',
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ps.add_argument('-ml',
                    '--min',
                    dest='ml',
                    nargs='+',
                    default=('0.', '0.'),
                    help='Mesh lower coordinates, space- or comma-separated')
    ps.add_argument('-mu',
                    '--max',
                    dest='mu',
                    nargs='+',
                    default=('1.', '1.'),
                    help='Mesh upper coordinates, space- or comma-separated')
    ps.add_argument('-ms',
                    '--res',
                    dest='ms',
                    nargs='+',
                    default=('20', '20'),
                    help='Logical size of mesh (cells), space- or comma-separated')
    ps.add_argument('-o', '--output', type=str, default='mcmesh', help='Output file base name')
    ps.add_argument('--topology',
                    choices=('structured', 'unstructured'),
                    default='structured',
                    help='Topology representation to write')
    ps.add_argument('--field',
                    choices=('sphere', 'plane'),
                    default='sphere',
                    help='Vertex field to sample')
    ps.add_argument(
        '--center',
        nargs='+',
        default=None,
        help='Sphere center or point on plane, space- or comma-separated. Defaults to mesh center')
    ps.add_argument(
        '--radius',
        type=float,
        default=None,
        help='Sphere/circle radius. Defaults to one quarter of the shortest mesh extent')
    ps.add_argument(
        '--normal',
        nargs='+',
        default=None,
        help='Plane normal, space- or comma-separated. Defaults to +z in 3D or +y in 2D')
    ps.add_argument('--fieldName', type=str, default='fcn', help='Output vertex field name')
    ps.add_argument('--topologyName', type=str, default='mesh', help='Output topology name')
    ps.add_argument('--coordsetName', type=str, default='coords', help='Output coordset name')
    ps.add_argument('--verbose', action='store_true', help='Print additional info')

    opts, unkn = ps.parse_known_args()
    if opts.verbose:
        print(opts, unkn)
    if unkn:
        print("Unrecognized arguments:", *unkn)
        quit(1)
    opts.ml = parse_component_list(opts.ml, float)
    opts.mu = parse_component_list(opts.mu, float)
    opts.ms = parse_component_list(opts.ms, int)
    if opts.center is not None:
        opts.center = parse_component_list(opts.center, float)
    if opts.normal is not None:
        opts.normal = parse_component_list(opts.normal, float)
    return opts


def validated_mesh_options(opts):
    dim = len(opts.ms)
    if dim not in (2, 3) or len(opts.ml) != dim or len(opts.mu) != dim:
        raise RuntimeError('ms, ml and mu options must have the same dimensions (2 or 3)')
    if any(s <= 0 for s in opts.ms):
        raise RuntimeError(f'ms ({opts.ms}) entries must be positive')

    mesh_size = np.array(opts.ms, dtype=np.int64)
    mesh_lower = np.array(opts.ml, dtype=np.float64)
    mesh_upper = np.array(opts.mu, dtype=np.float64)
    mesh_extent = mesh_upper - mesh_lower
    if np.any(mesh_extent <= 0.0):
        raise RuntimeError(f'mu ({opts.mu}) must be greater than ml ({opts.ml}) in all directions')

    center = np.array(opts.center if opts.center is not None else 0.5 * (mesh_lower + mesh_upper),
                      dtype=np.float64)
    if len(center) != dim:
        raise RuntimeError(f'center must have {dim} components')

    radius = opts.radius if opts.radius is not None else 0.25 * np.min(mesh_extent)

    normal = np.array(opts.normal if opts.normal is not None else ((0., 1.) if dim == 2 else
                                                                   (0., 0., 1.)),
                      dtype=np.float64)
    if len(normal) != dim:
        raise RuntimeError(f'normal must have {dim} components')
    normal_norm = np.linalg.norm(normal)
    if normal_norm == 0.0:
        raise RuntimeError('normal must be nonzero')
    normal = normal / normal_norm

    return {
        'dim': dim,
        'mesh_size': mesh_size,
        'mesh_lower': mesh_lower,
        'mesh_extent': mesh_extent,
        'center': center,
        'radius': radius,
        'normal': normal,
    }


def sample_field(pt, field_kind, center, radius, normal):
    if field_kind == 'sphere':
        return np.linalg.norm(pt - center) - radius
    return np.dot(pt - center, normal)


def node_index(mesh_size, i, j, k=0):
    ni = mesh_size[0] + 1
    nj = mesh_size[1] + 1
    return i + j * ni + k * ni * nj


def generate_coordset(mesh, opts, context):
    dim = context['dim']
    mesh_size = context['mesh_size']
    mesh_lower = context['mesh_lower']
    mesh_extent = context['mesh_extent']

    node_counts = mesh_size + 1
    num_nodes = int(np.prod(node_counts))
    coords = np.empty((num_nodes, dim), dtype=np.float64)

    idx = 0
    if dim == 2:
        for j in range(node_counts[1]):
            for i in range(node_counts[0]):
                logical = np.array((i, j), dtype=np.float64) / mesh_size
                coords[idx, :] = mesh_lower + logical * mesh_extent
                idx += 1
    else:
        for k in range(node_counts[2]):
            for j in range(node_counts[1]):
                for i in range(node_counts[0]):
                    logical = np.array((i, j, k), dtype=np.float64) / mesh_size
                    coords[idx, :] = mesh_lower + logical * mesh_extent
                    idx += 1

    coordset = mesh[f'coordsets/{opts.coordsetName}']
    coordset['type'] = 'explicit'
    coordset['values/x'].set(coords[:, 0])
    coordset['values/y'].set(coords[:, 1])
    if dim == 3:
        coordset['values/z'].set(coords[:, 2])

    return coords


def generate_topology(mesh, opts, context):
    dim = context['dim']
    mesh_size = context['mesh_size']

    topo = mesh[f'topologies/{opts.topologyName}']
    topo['coordset'] = opts.coordsetName
    if opts.topology == 'structured':
        topo['type'] = 'structured'
        topo['elements/dims/i'] = int(mesh_size[0])
        topo['elements/dims/j'] = int(mesh_size[1])
        if dim == 3:
            topo['elements/dims/k'] = int(mesh_size[2])
        return

    topo['type'] = 'unstructured'
    topo['elements/shape'] = 'quad' if dim == 2 else 'hex'
    connectivity = []
    if dim == 2:
        for j in range(mesh_size[1]):
            for i in range(mesh_size[0]):
                connectivity.extend([
                    node_index(mesh_size, i, j),
                    node_index(mesh_size, i + 1, j),
                    node_index(mesh_size, i + 1, j + 1),
                    node_index(mesh_size, i, j + 1),
                ])
    else:
        for k in range(mesh_size[2]):
            for j in range(mesh_size[1]):
                for i in range(mesh_size[0]):
                    connectivity.extend([
                        node_index(mesh_size, i, j, k),
                        node_index(mesh_size, i + 1, j, k),
                        node_index(mesh_size, i + 1, j + 1, k),
                        node_index(mesh_size, i, j + 1, k),
                        node_index(mesh_size, i, j, k + 1),
                        node_index(mesh_size, i + 1, j, k + 1),
                        node_index(mesh_size, i + 1, j + 1, k + 1),
                        node_index(mesh_size, i, j + 1, k + 1),
                    ])
    topo['elements/connectivity'].set(np.array(connectivity, dtype=np.int64))


def generate_fields(mesh, opts, context, coords):
    values = np.empty(coords.shape[0], dtype=np.float64)
    for idx, pt in enumerate(coords):
        values[idx] = sample_field(pt, opts.field, context['center'], context['radius'],
                                   context['normal'])

    field = mesh[f'fields/{opts.fieldName}']
    field['topology'] = opts.topologyName
    field['association'] = 'vertex'
    field['values'].set(values)


def generate_mesh(opts):
    context = validated_mesh_options(opts)
    mesh = conduit.Node()

    coords = generate_coordset(mesh, opts, context)
    generate_topology(mesh, opts, context)
    generate_fields(mesh, opts, context, coords)
    mesh['state/domain_id'] = 0

    return mesh


def main():
    opts = parse_args()
    mesh = generate_mesh(opts)

    info = conduit.Node()
    if not conduit.blueprint.mesh.verify(mesh, info):
        print("Mesh failed blueprint verification. Info:")
        print(info)
        return 2

    if opts.verbose:
        print(mesh)

    conduit.relay.io.blueprint.save_mesh(mesh, opts.output, "hdf5")
    print(f'Wrote mesh {opts.output}')
    return 0


if __name__ == '__main__':
    exit(main())
