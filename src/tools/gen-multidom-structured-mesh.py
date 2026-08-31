#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

# Write a multidomain Blueprint mesh for testing.
#
# You need Conduit's Python module available on PYTHONPATH.
# Axom's build tree usually provides `bin/run_python_with_axom.sh` for that.
#
# The generated Blueprint hierarchy is:
#   <bp_root>
#    ├── <domain_i_j[_k]>          Map entry, or list child with --useList
#    │    ├── topologies
#    │    │   └── mesh
#    │    │        ├─• type        == "structured" or "unstructured"
#    │    │        ├─• coordset    == "coords"
#    │    │        └── elements
#    │    │             ├── dims              Structured only
#    │    │             │    ├─• i
#    │    │             │    ├─• j
#    │    │             │    └─• [k]
#    │    │             ├─• shape             Unstructured only, quad or hex
#    │    │             └─• connectivity      Unstructured only, flat int64 connectivity
#    │    ├── coordsets
#    │    │   └── coords
#    │    │        ├─• type        == "explicit"
#    │    │        └── values      i-fastest node ordering, ghost padded for --strided
#    │    │             ├─• x
#    │    │             ├─• y
#    │    │             └─• [z]
#    │    └── fields
#    │        ├── field                        Conduit's example field
#    │        │    ├─• association == "element"
#    │        │    ├─• topology    == "mesh"
#    │        │    └─• values
#    │        └── <fieldName>                  Only with --field. Nodal for MarchingCubes
#    │             ├─• association == "vertex"
#    │             ├─• topology    == "mesh"
#    │             └─• values
#    └── ...

try:
    import conduit
    import conduit.blueprint
    import conduit.relay
except ModuleNotFoundError as e:
    print(f'{e}\n'
          'Add Conduit to PYTHONPATH, for example:\n'
          '  export PYTHONPATH=/path/to/conduit/install/python-modules:$PYTHONPATH\n'
          'If you have an Axom build directory, you can also run:\n'
          '  /path/to/axom_build_dir/bin/run_python_with_axom.sh <this_script.py> ...\n'
          'Note: HDF5 support is only needed when using `--protocol hdf5`.')
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
    ps = ArgumentParser(description='Write a multidomain Blueprint mesh.',
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ps.add_argument('--useList', action='store_true', help='Put domains in a list instead of a map')
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
                    default=('3', '3'),
                    help='Logical size of mesh (cells), space- or comma-separated')
    ps.add_argument('-dc',
                    '--domains',
                    dest='dc',
                    nargs='+',
                    default=('1', '1'),
                    help='Domain counts in each index direction, space- or comma-separated')
    ps.add_argument('-o', '--output', type=str, default='mdmesh', help='Output file base name')
    ps.add_argument('--strided', action='store_true', help='Use strided_structured (has ghosts)')
    ps.add_argument(
        '--topology',
        choices=('structured', 'unstructured'),
        default='structured',
        help='Topology type. "unstructured" emits single-shape quad or hex connectivity '
        'over the same nodes. Incompatible with --strided.')
    ps.add_argument(
        '--field',
        choices=('none', 'sphere', 'plane'),
        default='none',
        help='Add an analytic nodal field, which MarchingCubes needs. The Conduit example '
        'field is element-associated.')
    ps.add_argument('--fieldName', type=str, default='fcn', help='Name of the analytic nodal field')
    ps.add_argument('--center',
                    nargs='+',
                    default=None,
                    help='Center for --field sphere, or point on plane for --field plane. '
                    'Defaults to the mesh center.')
    ps.add_argument(
        '--radius',
        type=float,
        default=None,
        help='Radius for --field sphere. Defaults to one quarter of the shortest mesh extent.')
    ps.add_argument('--normal',
                    nargs='+',
                    default=None,
                    help='Normal direction for --field plane. Defaults to +y in 2D and +z in 3D.')
    ps.add_argument('--protocol',
                    choices=('hdf5', 'json', 'yaml'),
                    default='hdf5',
                    help='Conduit relay output protocol. json/yaml let readers run without HDF5.')
    ps.add_argument('-v', '--verbose', action='store_true', help='Print additional info')

    opts, unkn = ps.parse_known_args()
    opts.ml = parse_component_list(opts.ml, float)
    opts.mu = parse_component_list(opts.mu, float)
    opts.ms = parse_component_list(opts.ms, int)
    opts.dc = parse_component_list(opts.dc, int)

    if opts.verbose:
        print(opts, unkn)
    if unkn:
        print("Unrecognized arguments:", *unkn)
        quit(1)
    return opts


def validated_mesh_options(opts):
    dim = len(opts.dc)

    if dim not in (2, 3) or len(opts.ms) != dim or len(opts.ml) != dim or len(opts.mu) != dim:
        raise RuntimeError('dc, ms, ml and mu options must have the same dimensions (2 or 3)')

    if any(s <= 0 for s in opts.ms):
        raise RuntimeError(f'ms ({opts.ms}) entries must be positive')
    if any(d <= 0 for d in opts.dc):
        raise RuntimeError(f'dc ({opts.dc}) entries must be positive')

    # Must have enough cells for requested partitioning.
    if any(opts.ms[i] < opts.dc[i] for i in range(dim)):
        raise RuntimeError(f'ms ({opts.ms}) must be >= dc ({opts.dc}) in all directions.')

    mesh_size = np.array(opts.ms, dtype=int)
    mesh_lower = np.array(opts.ml, dtype=float)
    mesh_upper = np.array(opts.mu, dtype=float)
    mesh_extent = mesh_upper - mesh_lower
    if np.any(mesh_extent <= 0.0):
        raise RuntimeError(f'mu ({opts.mu}) must be greater than ml ({opts.ml}) in all directions')

    domain_counts = opts.dc if dim == 3 else (*opts.dc, 1)
    domain_counts = np.array(domain_counts, dtype=int)

    domain_size = mesh_size // domain_counts[:dim]
    domain_size_remainder = mesh_size % domain_counts[:dim]

    if opts.topology == 'unstructured' and opts.strided:
        raise RuntimeError(
            '--topology unstructured is incompatible with --strided. The ghost padded '
            'coordset does not have compact node numbering to build connectivity over.')

    mesh_center = 0.5 * (mesh_lower + mesh_upper)
    if opts.center is None:
        center = mesh_center
    else:
        center = np.array(parse_component_list(opts.center, float), dtype=float)
        if len(center) < dim:
            raise RuntimeError(f'--center ({opts.center}) needs at least {dim} components')

    if opts.radius is None:
        radius = 0.25 * float(np.min(mesh_extent))
    else:
        radius = float(opts.radius)
    if radius <= 0.0:
        raise RuntimeError(f'--radius must be positive (got {radius})')

    if opts.normal is None:
        normal = np.array((0.0, 1.0) if dim == 2 else (0.0, 0.0, 1.0), dtype=float)
    else:
        normal = np.array(parse_component_list(opts.normal, float), dtype=float)
    if opts.field != 'none':
        if opts.field == 'plane' and len(normal) < dim:
            raise RuntimeError(f'--normal ({opts.normal}) needs at least {dim} components')
        if opts.field == 'plane':
            normal_norm = float(np.linalg.norm(normal[:dim]))
            if normal_norm == 0.0:
                raise RuntimeError('--normal must be nonzero for --field plane')
            normal = normal / normal_norm

    return {
        'center': center,
        'normal': normal,
        'radius': radius,
        'dim': dim,
        'domain_counts': domain_counts,
        'mesh_size': mesh_size,
        'mesh_lower': mesh_lower,
        'mesh_upper': mesh_upper,
        'cell_physical_size': mesh_extent / mesh_size,
        'domain_size': domain_size,
        'domain_size_remainder': domain_size_remainder,
        'num_phony_nodes_left': 2 if opts.strided else 0,
        'num_phony_nodes_right': 1 if opts.strided else 0,
    }


def domain_index_begin(context, di, dj, dk=None):
    '''Compute first cell index of the domain with multi-dimensional index (di, dj, dk).'''
    dim = context['dim']
    ds = (di, dj) if dim == 2 else (di, dj, dk)
    idx = np.array(ds)
    std = context['domain_size'] * ds
    extra = np.where(idx < context['domain_size_remainder'][:dim], idx,
                     context['domain_size_remainder'][:dim])
    return std + extra


def domain_node(md_mesh, opts, di, dj, dk):
    if opts.useList:
        return md_mesh.append()

    dom_name = f'domain_{di:1d}_{dj:1d}'
    if len(opts.dc) == 3:
        dom_name += f'_{dk:1d}'
    return md_mesh[dom_name]


def generate_topology(dom, opts, context, point_counts, cell_start, cell_end):
    '''Generate a structured Blueprint topology and seed matching example data.'''
    dim = context['dim']
    npnl = context['num_phony_nodes_left']
    npnr = context['num_phony_nodes_right']

    point_counts_3 = point_counts if len(point_counts) == 3 else (*point_counts, 0)
    if opts.strided:
        elem_extents = (cell_end - cell_start) + (npnl + npnr + 1)
        vert_extents = np.array(point_counts) + (npnl + npnr)
        elem_offset = np.full(dim, npnl)
        vert_offset = np.full(dim, npnl)

        desc = conduit.Node()
        desc['vertex_data/shape'].set(vert_extents)
        desc['vertex_data/origin'].set(vert_offset)
        desc['element_data/shape'].set(elem_extents)
        desc['element_data/origin'].set(elem_offset)
        conduit.blueprint.mesh.examples.strided_structured(desc, *point_counts_3, dom)
        if dom.has_child("state"):
            dom.remove_child("state")
    else:
        conduit.blueprint.mesh.examples.basic('structured', *point_counts_3, dom)


def generate_coordset(dom, context, start_coord, end_coord):
    '''Scale and shift the generated explicit coordset to the requested domain bounds.'''
    npnl = context['num_phony_nodes_left']
    npnr = context['num_phony_nodes_right']

    ndim = dom['coordsets/coords/values'].number_of_children()

    dom_lens_node = dom['topologies/mesh/elements/dims']
    dirs = 'ij' if ndim == 2 else 'ijk'
    dom_lens = np.array([dom_lens_node[d] for d in dirs])
    domain_physical_size = np.array(end_coord) - np.array(start_coord)

    assert (dom['topologies/mesh/type'] == 'structured')
    assert (len(start_coord) >= ndim)
    assert (len(domain_physical_size) >= ndim)

    coord_array_lens = dom_lens + 1 + npnl + npnr

    xyz = 'xyz'
    for d in range(ndim):
        coords = dom['coordsets/coords/values'][d]
        coords = np.reshape(coords, np.flip(coord_array_lens))

        # real_coords excludes the ghost layers.
        if ndim == 2:
            real_coords = coords[npnl:, npnl:] if npnr == 0 else coords[npnl:-npnr, npnl:-npnr]
        else:
            real_coords = coords[npnl:, npnl:,
                                 npnl:] if npnr == 0 else coords[npnl:-npnr, npnl:-npnr, npnl:-npnr]

        min_coord, max_coord = np.amin(real_coords), np.amax(real_coords)
        cur_range = max_coord - min_coord
        coords = (coords - min_coord) * domain_physical_size[d] / cur_range + start_coord[d]
        dom['coordsets/coords/values'][xyz[d]] = coords


def add_analytic_nodal_field(dom, opts, context):
    '''Add a nodal scalar field sampled at every coordset node. The implementation is vectorized'''

    if opts.field == 'none':
        return

    dim = context['dim']
    vals = dom['coordsets/coords/values']
    comps = [np.asarray(vals[c], dtype=np.float64) for c in 'xyz'[:dim]]
    pts = np.stack(comps, axis=1)

    center = np.array(context['center'][:dim], dtype=np.float64)
    if opts.field == 'sphere':
        values = np.linalg.norm(pts - center, axis=1) - context['radius']
    else:
        normal = np.array(context['normal'][:dim], dtype=np.float64)
        values = (pts - center) @ normal

    field = dom[f'fields/{opts.fieldName}']
    field['topology'] = 'mesh'
    field['association'] = 'vertex'
    field['values'] = values


def structured_to_unstructured(dom, context):
    '''Rewrite the structured topology as single-shape quad/hex connectivity.

    This uses numpy broadcasting instead of a per-cell Python loop. Node order
    within a cell matches Blueprint's quad/hex convention. Node ids follow the
    coordset's i-fastest numbering, so the coordset stays as-is.
    '''
    dim = context['dim']
    topo = dom['topologies/mesh']
    dims = topo['elements/dims']
    cells = np.array([dims['i'], dims['j']] + ([dims['k']] if dim == 3 else []), dtype=np.int64)
    pts = cells + 1

    if dim == 2:
        j, i = np.meshgrid(np.arange(cells[1]), np.arange(cells[0]), indexing='ij')
        base = (i + j * pts[0]).ravel()
        offs = [0, 1, 1 + pts[0], pts[0]]
    else:
        k, j, i = np.meshgrid(np.arange(cells[2]),
                              np.arange(cells[1]),
                              np.arange(cells[0]),
                              indexing='ij')
        base = (i + j * pts[0] + k * pts[0] * pts[1]).ravel()
        pij = pts[0] * pts[1]
        offs = [0, 1, 1 + pts[0], pts[0], pij, pij + 1, pij + 1 + pts[0], pij + pts[0]]

    conn = (base[:, None] + np.array(offs, dtype=np.int64)[None, :]).ravel()

    topo.remove_child('elements')
    topo['type'] = 'unstructured'
    topo['elements/shape'] = 'quad' if dim == 2 else 'hex'
    topo['elements/connectivity'] = conn


def generate_fields(dom, opts, context):
    '''Keep Conduit's example element field and ensure it references the generated topology.'''
    del opts, context
    if dom.has_path('fields/field'):
        dom['fields/field/topology'] = 'mesh'
        dom['fields/field/association'] = 'element'


def generate_domain(md_mesh, opts, context, di, dj, dk):
    dim = context['dim']
    mesh_lower = context['mesh_lower']
    cell_physical_size = context['cell_physical_size']

    dom = domain_node(md_mesh, opts, di, dj, dk)

    cell_start = domain_index_begin(context, di, dj, dk)
    cell_end = domain_index_begin(context, di + 1, dj + 1, dk + 1 if dim == 3 else 0)
    point_counts = cell_end - cell_start + 1

    generate_topology(dom, opts, context, point_counts, cell_start, cell_end)

    dom_lower = mesh_lower[:dim] + cell_start * cell_physical_size[:dim]
    dom_upper = mesh_lower[:dim] + cell_end * cell_physical_size[:dim]
    generate_coordset(dom, context, dom_lower, dom_upper)
    generate_fields(dom, opts, context)
    # Sample the field before rewriting the topology. The unstructured rewrite drops
    # elements/dims, and we still need those to size the coordset.
    add_analytic_nodal_field(dom, opts, context)
    if opts.topology == 'unstructured':
        structured_to_unstructured(dom, context)


def generate_mesh(opts):
    context = validated_mesh_options(opts)
    domain_counts = context['domain_counts']

    if opts.verbose:
        print(f"meshSize={context['mesh_size']} cells, domCounts={domain_counts[0:context['dim']]}"
              f" domSize={context['domain_size']}"
              f" domSizeRem={context['domain_size_remainder']}")

    md_mesh = conduit.Node()
    for dk in range(domain_counts[2]):
        for dj in range(domain_counts[1]):
            for di in range(domain_counts[0]):
                generate_domain(md_mesh, opts, context, di, dj, dk)

    return md_mesh


def main():
    opts = parse_args()
    md_mesh = generate_mesh(opts)

    if opts.verbose:
        print('mdMesh:')
        print(md_mesh)

    info = conduit.Node()
    if not conduit.blueprint.mesh.verify(md_mesh, info):
        print("Mesh failed blueprint verification. Info:")
        print(info)
        return 2

    conduit.relay.io.blueprint.save_mesh(md_mesh, opts.output, opts.protocol)
    print(f'Wrote mesh {opts.output}')
    return 0


if __name__ == '__main__':
    exit(main())
