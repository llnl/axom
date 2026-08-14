#!/usr/bin/env python3

# gen-multidom-point-mesh.py
# Write a point mesh following the Conduit mesh blueprint.
#
# The generated Blueprint hierarchy is:
#   <bp_root>                         (single-domain output)
#   <bp_root>/<domainName>            (multidomain output)
#    |-- state
#    |   `-- domain_id                == <global domain id>
#    |-- topologies
#    |   `-- <topologyName>
#    |        |-- coordset            == <coordsetName>
#    |        |-- type                == "points" or "unstructured"
#    |        `-- elements            (only for "unstructured")
#    |             |-- shape          == "point"
#    |             |-- connectivity   (int32, point_count)
#    |             |-- sizes          (int32, point_count, all 1)
#    |             `-- offsets        (int32, point_count)
#    |-- coordsets
#    |   `-- <coordsetName>
#    |        |-- type                == "explicit"
#    |        `-- values
#    |             |-- x              (float64, point_count)
#    |             |-- y              (float64, point_count)
#    |             `-- [z]            (float64, point_count, present in 3D)
#    `-- fields                       (optional; omitted when there are no fields)
#        `-- global_id                (only with --id-field)
#             |-- topology            == <topologyName>
#             |-- association         == "vertex"
#             `-- values              (int64, point_count)

# This script requires a conduit installation configured with python3 and hdf5.
# Make sure PYTHONPATH includes /path/to/conduit/install/python-modules,
# or use Axom's convenience script /path/to/axom_build_dir/bin/run_python_with_axom.sh
# that includes Conduit in PYTHONPATH.

import itertools
import math
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, ArgumentTypeError

try:
    import conduit
    import conduit.blueprint
    import conduit.relay
except ModuleNotFoundError as e:
    print(
        f'{e}\nMake sure your PYTHONPATH includes /path/to/conduit/install/python-modules\n'
        'Conduit must be configured with python and hdf5.\n'
        'Alternatively, you can use the convenience script\n'
        '/path/to/axom_build_dir/bin/run_python_with_axom.sh\n'
        'that includes Conduit in PYTHONPATH.'
    )
    sys.exit(-1)

try:
    import numpy as np
except ModuleNotFoundError as e:
    print(f'{e}\nThis script requires numpy.')
    sys.exit(-1)


AXES = 'xyz'


def get_mpi_context():
    '''Return MPI context when mpi4py and Conduit MPI relay are available.'''
    try:
        from mpi4py import MPI
        import conduit.blueprint.mpi.mesh
        import conduit.relay.mpi
        import conduit.relay.mpi.io
        import conduit.relay.mpi.io.blueprint
    except ModuleNotFoundError:
        return {'enabled': False, 'comm': None, 'rank': 0, 'size': 1}

    comm = MPI.COMM_WORLD.py2f()
    return {
        'enabled': True,
        'comm': comm,
        'rank': conduit.relay.mpi.rank(comm),
        'size': conduit.relay.mpi.size(comm),
    }


def csv_values(s, converter, value_name):
    '''Convert a comma-separated string to a list.'''
    try:
        values = [converter(token) for token in s.split(',')]
    except ValueError as e:
        raise ArgumentTypeError(f'{value_name} must be a comma-separated list') from e

    if len(values) == 0:
        raise ArgumentTypeError(f'{value_name} must not be empty')

    return values


def i_c(s):
    '''Convert comma-separated string to list of integers.'''
    return csv_values(s, int, 'integer values')


def f_c(s):
    '''Convert comma-separated string to list of floating point numbers.'''
    return csv_values(s, float, 'floating point values')


def positive_int(s):
    '''Convert a string to a positive integer.'''
    try:
        value = int(s)
    except ValueError as e:
        raise ArgumentTypeError('value must be an integer') from e

    if value < 1:
        raise ArgumentTypeError('value must be positive')

    return value


def parse_args():
    parser = ArgumentParser(
        description='Write a single-domain or multidomain blueprint point mesh.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'shape',
        choices=('circle', 'sphere', 'torus', 'grid', 'gaussian', 'uniform'),
        help='Point distribution to generate')
    parser.add_argument(
        '-d',
        '--dimension',
        type=int,
        choices=(2, 3),
        help='Spatial dimension for dimension-flexible shapes')
    parser.add_argument(
        '-n',
        '--point-count',
        type=positive_int,
        default=1024,
        help='Number of points, or target number for grid when --grid-size is omitted')
    parser.add_argument(
        '--grid-size',
        type=i_c,
        help='Point counts in each grid direction, e.g. 101,101 or 65,65,65')
    parser.add_argument(
        '-dc',
        '--domain-counts',
        type=i_c,
        help='Domain counts in one chunk direction or each coordinate direction')
    parser.add_argument(
        '--domains',
        type=positive_int,
        help='Total domain count for a one-dimensional chunk partition')
    parser.add_argument(
        '--single-domain',
        action='store_true',
        help='Write a single-domain point mesh instead of multidomain output')
    parser.add_argument('-ml', type=f_c, help='Lower coordinates for grid and uniform shapes')
    parser.add_argument('-mu', type=f_c, help='Upper coordinates for grid and uniform shapes')
    parser.add_argument('--center', type=f_c, help='Center for circle, sphere, torus, and gaussian')
    parser.add_argument('--radius', type=float, default=1.0, help='Radius for circle or sphere')
    parser.add_argument('--major-radius', type=float, default=1.0, help='Major radius for torus')
    parser.add_argument('--minor-radius', type=float, default=0.25, help='Minor radius for torus')
    parser.add_argument('--stddev',
                        type=f_c,
                        help='Gaussian standard deviation; scalar or one value per dimension')
    parser.add_argument('--seed', type=int, default=0, help='Random seed for stochastic shapes')
    parser.add_argument('--use-list',
                        action='store_true',
                        help='Put domains in a list instead of a map')
    parser.add_argument('--topology-name', default='mesh', help='Blueprint topology name')
    parser.add_argument('--coordset-name', default='coords', help='Blueprint coordset name')
    parser.add_argument(
        '--topology-type',
        choices=('unstructured', 'points'),
        default='unstructured',
        help='Blueprint point topology representation')
    parser.add_argument(
        '--id-field',
        action='store_true',
        help='Add a vertex-associated int64 global_id field')
    parser.add_argument('--protocol', default='hdf5', help='Conduit relay protocol for save_mesh')
    parser.add_argument('-o', '--output', default='mdpoint', help='Output file base name')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print additional info')

    opts, unkn = parser.parse_known_args()
    if opts.verbose:
        print(opts, unkn)
    if unkn:
        print('Unrecognized arguments:', *unkn)
        sys.exit(1)

    return opts


def infer_dimension(opts):
    '''Infer and validate the spatial dimension.'''
    fixed_dims = {'circle': 2, 'sphere': 3, 'torus': 3}
    if opts.shape in fixed_dims:
        dim = fixed_dims[opts.shape]
        if opts.dimension is not None and opts.dimension != dim:
            raise RuntimeError(f'{opts.shape} requires dimension {dim}')
        return dim

    dim = opts.dimension
    inferred = []
    for values in (opts.grid_size, opts.ml, opts.mu, opts.center, opts.domain_counts):
        if values is not None:
            inferred.append(len(values))

    if dim is None and inferred:
        dim = inferred[0]
    if dim is None:
        dim = 2

    if dim not in (2, 3):
        raise RuntimeError('dimension must be 2 or 3')

    for value_dim in inferred:
        if value_dim not in (1, dim):
            raise RuntimeError(
                'dimensioned options must have one value or one value per spatial dimension')

    return dim


def vector_option(values, dim, default, name):
    '''Return a dimension-sized numpy vector from an optional scalar/list argument.'''
    if values is None:
        return np.full(dim, default, dtype=float)
    if len(values) == 1:
        return np.full(dim, values[0], dtype=float)
    if len(values) != dim:
        raise RuntimeError(f'{name} must have one value or {dim} values')
    return np.array(values, dtype=float)


def domain_counts_from_options(opts, dim, mpi_size):
    '''Return domain counts as an integer vector.'''
    if opts.single_domain:
        return np.array([1], dtype=int)

    if opts.domain_counts is not None and opts.domains is not None:
        raise RuntimeError('Use --domain-counts or --domains, not both')

    if opts.domain_counts is not None:
        counts = np.array(opts.domain_counts, dtype=int)
        if len(counts) not in (1, dim):
            raise RuntimeError(
                f'--domain-counts must have one value or {dim} values; use --domains for chunks')
    elif opts.domains is not None:
        counts = np.array([opts.domains], dtype=int)
    else:
        counts = np.array([mpi_size], dtype=int)

    if np.any(counts < 1):
        raise RuntimeError('domain counts must be positive')

    return counts


def grid_counts_from_options(opts, dim):
    '''Return point counts for the regular grid shape.'''
    if opts.grid_size is not None:
        counts = np.array(opts.grid_size, dtype=int)
        if len(counts) != dim:
            raise RuntimeError(f'--grid-size must have {dim} values')
    else:
        base = max(1, int(round(opts.point_count**(1.0 / dim))))
        counts = np.full(dim, base, dtype=int)
        while int(np.prod(counts)) < opts.point_count:
            counts[np.argmin(counts)] += 1

    if np.any(counts < 1):
        raise RuntimeError('grid sizes must be positive')

    return counts


def split_range(total, parts, part):
    '''Return the half-open item range for one chunk in an even partition.'''
    base = total // parts
    rem = total % parts
    begin = part * base + min(part, rem)
    end = begin + base + (1 if part < rem else 0)
    return begin, end


def spatially_sort_points(points):
    '''Return points ordered so contiguous chunks have reasonable spatial coherence.'''
    if points.shape[0] == 0:
        return points

    spans = np.ptp(points, axis=0)
    primary = int(np.argmax(spans))
    secondary = [d for d in range(points.shape[1]) if d != primary]
    keys = [points[:, d] for d in reversed(secondary)]
    keys.append(points[:, primary])
    return points[np.lexsort(tuple(keys))]


def split_points(points, domain_counts):
    '''Split points into nearly equal contiguous chunks.'''
    points = spatially_sort_points(points)
    point_count = points.shape[0]
    domain_count = int(np.prod(domain_counts))
    if point_count < domain_count:
        raise RuntimeError(
            f'point count ({point_count}) must be at least the domain count ({domain_count})')

    chunks = []
    domain_indices = list(itertools.product(*[range(c) for c in domain_counts]))
    for domain_id, domain_index in enumerate(domain_indices):
        begin, end = split_range(point_count, domain_count, domain_id)
        chunks.append((domain_id, domain_index, points[begin:end], begin))
    return chunks


def domain_name(domain_index):
    '''Create a stable map key for a domain.'''
    if len(domain_index) == 1:
        return f'domain_{domain_index[0]:06d}'
    return 'domain_' + '_'.join(f'{idx:03d}' for idx in domain_index)


def generate_circle(point_count, center, radius):
    theta = np.linspace(0.0, 2.0 * math.pi, point_count, endpoint=False)
    points = np.empty((point_count, 2), dtype=np.float64)
    points[:, 0] = center[0] + radius * np.cos(theta)
    points[:, 1] = center[1] + radius * np.sin(theta)
    return points


def generate_sphere(point_count, center, radius):
    # Fibonacci sphere points give a deterministic, nearly uniform surface sampling.
    idx = np.arange(point_count, dtype=np.float64) + 0.5
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    z = 1.0 - 2.0 * idx / point_count
    rxy = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    theta = golden_angle * idx

    points = np.empty((point_count, 3), dtype=np.float64)
    points[:, 0] = center[0] + radius * rxy * np.cos(theta)
    points[:, 1] = center[1] + radius * rxy * np.sin(theta)
    points[:, 2] = center[2] + radius * z
    return points


def generate_torus(point_count, center, major_radius, minor_radius):
    # Use a deterministic irrational sequence to avoid requiring two surface counts.
    idx = np.arange(point_count, dtype=np.float64)
    golden_ratio_conj = (math.sqrt(5.0) - 1.0) / 2.0
    theta = 2.0 * math.pi * idx / point_count
    phi = 2.0 * math.pi * np.mod(idx * golden_ratio_conj, 1.0)
    ring_radius = major_radius + minor_radius * np.cos(phi)

    points = np.empty((point_count, 3), dtype=np.float64)
    points[:, 0] = center[0] + ring_radius * np.cos(theta)
    points[:, 1] = center[1] + ring_radius * np.sin(theta)
    points[:, 2] = center[2] + minor_radius * np.sin(phi)
    return points


def generate_gaussian(point_count, center, stddev, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(loc=center, scale=stddev, size=(point_count, len(center)))


def generate_uniform(point_count, lower, upper, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(low=lower, high=upper, size=(point_count, len(lower)))


def generate_grid_points(grid_counts, lower, upper):
    axes = [np.linspace(lower[d], upper[d], grid_counts[d]) for d in range(len(grid_counts))]
    coords = np.meshgrid(*axes, indexing='ij')
    return np.column_stack([coord.ravel() for coord in coords])


def generate_grid_domains(grid_counts, lower, upper, domain_counts):
    '''Generate grid points one spatial domain at a time.'''
    for d in range(len(grid_counts)):
        if domain_counts[d] > grid_counts[d]:
            raise RuntimeError(
                f'domain count {domain_counts[d]} exceeds grid size {grid_counts[d]} '
                f'in {AXES[d]} direction')

    axes = [np.linspace(lower[d], upper[d], grid_counts[d]) for d in range(len(grid_counts))]
    chunks = []
    global_start = 0
    for domain_id, domain_index in enumerate(itertools.product(*[range(c) for c in domain_counts])):
        local_axes = []
        for d, idx in enumerate(domain_index):
            begin, end = split_range(grid_counts[d], domain_counts[d], idx)
            local_axes.append(axes[d][begin:end])
        coords = np.meshgrid(*local_axes, indexing='ij')
        points = np.column_stack([coord.ravel() for coord in coords])
        chunks.append((domain_id, domain_index, points, global_start))
        global_start += points.shape[0]
    return chunks


def generated_domains(opts, dim, lower, upper, center, stddev, domain_counts):
    '''Return generated domains as tuples: id, index, points, global_id_start.'''
    if opts.shape == 'grid':
        grid_counts = grid_counts_from_options(opts, dim)
        if opts.verbose:
            print(f'grid_size={grid_counts.tolist()} point_count={int(np.prod(grid_counts))}')
        if len(domain_counts) == dim:
            return generate_grid_domains(grid_counts, lower, upper, domain_counts)

        points = generate_grid_points(grid_counts, lower, upper)
        return split_points(points, domain_counts)

    if opts.shape == 'circle':
        points = generate_circle(opts.point_count, center, opts.radius)
    elif opts.shape == 'sphere':
        points = generate_sphere(opts.point_count, center, opts.radius)
    elif opts.shape == 'torus':
        points = generate_torus(opts.point_count, center, opts.major_radius, opts.minor_radius)
    elif opts.shape == 'gaussian':
        points = generate_gaussian(opts.point_count, center, stddev, opts.seed)
    elif opts.shape == 'uniform':
        points = generate_uniform(opts.point_count, lower, upper, opts.seed)
    else:
        raise RuntimeError(f'Unknown shape: {opts.shape}')

    return split_points(points, domain_counts)


def fill_domain(dom, opts, domain_id, points, global_id_start):
    '''Fill one blueprint domain with point mesh data.'''
    point_count = points.shape[0]
    dim = points.shape[1]

    dom['state/domain_id'] = int(domain_id)
    dom['coordsets/' + opts.coordset_name + '/type'] = 'explicit'
    for d in range(dim):
        values_path = f'coordsets/{opts.coordset_name}/values/{AXES[d]}'
        dom[values_path].set(np.ascontiguousarray(points[:, d], dtype=np.float64))

    topo = dom['topologies/' + opts.topology_name]
    topo['type'] = opts.topology_type
    topo['coordset'] = opts.coordset_name
    if opts.topology_type == 'unstructured':
        topo['elements/shape'] = 'point'
        topo['elements/connectivity'].set(np.arange(point_count, dtype=np.int32))
        topo['elements/sizes'].set(np.ones(point_count, dtype=np.int32))
        topo['elements/offsets'].set(np.arange(point_count, dtype=np.int32))

    if opts.id_field:
        field = dom['fields/global_id']
        field['association'] = 'vertex'
        field['topology'] = opts.topology_name
        field['values'].set(
            np.arange(global_id_start, global_id_start + point_count, dtype=np.int64))


def create_domain(md_mesh, opts, domain_id, domain_index, points, global_id_start):
    '''Append one domain to the multidomain mesh.'''
    if opts.use_list:
        dom = md_mesh.append()
    else:
        dom = md_mesh[domain_name(domain_index)]

    fill_domain(dom, opts, domain_id, points, global_id_start)


def local_chunks_for_rank(chunks, rank, size, single_domain):
    '''Return the chunks that should be written by this rank.'''
    if size <= 1:
        return chunks
    if single_domain:
        return chunks if rank == 0 else []

    begin, end = split_range(len(chunks), size, rank)
    return chunks[begin:end]


def verify_mesh(mesh, mpi):
    '''Verify a mesh with the serial or MPI Blueprint verifier.'''
    info = conduit.Node()
    if mpi['enabled'] and mpi['size'] > 1:
        valid = conduit.blueprint.mpi.mesh.verify(mesh, info, mpi['comm'])
    else:
        valid = conduit.blueprint.mesh.verify(mesh, info)

    return valid, info


def save_mesh(mesh, opts, mpi):
    '''Save a mesh with the serial or MPI Blueprint relay.'''
    if mpi['enabled'] and mpi['size'] > 1:
        try:
            conduit.relay.mpi.io.blueprint.save_mesh(mesh, opts.output, opts.protocol, mpi['comm'])
        except TypeError:
            conduit.relay.mpi.io.blueprint.save_mesh(mesh, opts.output, mpi['comm'], opts.protocol)
    else:
        conduit.relay.io.blueprint.save_mesh(mesh, opts.output, opts.protocol)


def main():
    mpi = get_mpi_context()
    opts = parse_args()
    dim = infer_dimension(opts)
    lower = vector_option(opts.ml, dim, -1.0, '-ml')
    upper = vector_option(opts.mu, dim, 1.0, '-mu')
    center = vector_option(opts.center, dim, 0.0, '--center')
    stddev = vector_option(opts.stddev, dim, 1.0, '--stddev')
    domain_counts = domain_counts_from_options(opts, dim, mpi['size'])

    if np.any(upper <= lower):
        raise RuntimeError('-mu values must be greater than -ml values')
    if opts.radius <= 0.0:
        raise RuntimeError('--radius must be positive')
    if opts.major_radius <= 0.0 or opts.minor_radius <= 0.0:
        raise RuntimeError('--major-radius and --minor-radius must be positive')
    if np.any(stddev <= 0.0):
        raise RuntimeError('--stddev values must be positive')

    chunks = generated_domains(opts, dim, lower, upper, center, stddev, domain_counts)
    local_chunks = local_chunks_for_rank(chunks, mpi['rank'], mpi['size'], opts.single_domain)

    mesh = conduit.Node()
    if opts.single_domain:
        if local_chunks:
            domain_id, _, points, global_id_start = local_chunks[0]
            fill_domain(mesh, opts, domain_id, points, global_id_start)
    else:
        for domain_id, domain_index, points, global_id_start in local_chunks:
            create_domain(mesh, opts, domain_id, domain_index, points, global_id_start)

    if opts.verbose:
        print(f'rank {mpi["rank"]} mesh:')
        print(mesh)

    valid, info = verify_mesh(mesh, mpi)
    if not valid:
        print('Mesh failed blueprint verification.  Info:')
        print(info)
        return 2

    save_mesh(mesh, opts, mpi)
    if mpi['rank'] == 0:
        total_points = sum(points.shape[0] for _, _, points, _ in chunks)
        mesh_kind = 'single-domain' if opts.single_domain else 'multidomain'
        print(
            f'Wrote {total_points} {dim}D {opts.shape} points as a {mesh_kind} mesh '
            f'with {len(chunks)} domain(s) to {opts.output} using {opts.protocol}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
