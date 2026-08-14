#!/usr/bin/env python3

# gen-multidom-point-mesh.py
# Write a mesh following the Conduit mesh blueprint for distributed closest point (DCP) testing.
#
# The generated Blueprint hierarchy is:
#   <bp_root>                         (single-domain output)
#   <bp_root>/<domainName>            (multidomain output)
#    |-- state
#    |   `-- domain_id                == <global domain id>
#    |   `-- description              (analytic/input description string; stored in state
#    |                                  so Blueprint relay preserves it)
#    |   `-- verification             (optional analytic verification metadata; stored in state)
#    |        |-- shape               == "circle", "sphere", "torus", "plane", or "annulus"
#    |        |-- dimension           == 2 or 3
#    |        |-- center              (float64, dimension)
#    |        |-- normal              (float64, dimension; plane only)
#    |        |-- radius              (circle/sphere)
#    |        |-- {major,minor}_radius (torus)
#    |        |-- {inner,outer}_radius (annulus)
#    |        |-- surface_tolerance   (float64 residual tolerance for shape membership)
#    |        `-- distance_tolerance  (float64 sampling slack for DCP checks)
#    |-- topologies
#    |   `-- <topologyName>
#    |        |-- coordset            == <coordsetName>
#    |        |-- type                == "points", "structured", or "unstructured"
#    |        `-- elements            (present for structured/unstructured)
#    |             |-- dims/{i,j,[k]}  (structured cell dimensions)
#    |             |-- shape          (unstructured "point", "quad", or "hex")
#    |             |-- connectivity   (unstructured int32 connectivity)
#    |             |-- sizes          (unstructured point meshes only)
#    |             `-- offsets        (unstructured point meshes only)
#    |-- coordsets
#    |   `-- <coordsetName>
#    |        |-- type                == "explicit"
#    |        `-- values              (i-fastest node ordering for grids)
#    |             |-- x              (float64, point_count)
#    |             |-- y              (float64, point_count)
#    |             `-- [z]            (float64, point_count, present in 3D)
#    `-- fields                       (optional; omitted when there are no fields)
#        `-- global_id                (only with --id-field)
#             |-- topology            == <topologyName>
#             |-- association         == "vertex"
#             `-- values              (int64, point_count)
#
# Usage examples:
#   # 2D circle object mesh split over N MPI ranks, with some empty domains:
#   srun -n N build-axom/bin/run_python_with_axom.sh \\
#     src/tools/gen-multidom-point-mesh.py circle --point-count 64 --domains 4 \\
#     --center 0.7,0.9 --radius 0.9 --output dcp_object_circle
#
#   # 3D sphere object mesh split over N MPI ranks using latitude/longitude samples:
#   srun -n N build-axom/bin/run_python_with_axom.sh \\
#     src/tools/gen-multidom-point-mesh.py sphere --long-point-count 60 \\
#     --lat-point-count 30 --center 0.7,0.9,0.5 --radius 0.9 \\
#     --domains 6 --output dcp_object_sphere
#
#   # 3D torus object mesh:
#   build-axom/bin/run_python_with_axom.sh \\
#     src/tools/gen-multidom-point-mesh.py torus --point-count 1024 \\
#     --center 0,0,0 --major-radius 0.75 --minor-radius 0.2 \\
#     --domains 4 --output dcp_object_torus
#
#   # 3D codimension-1 plane object mesh:
#   build-axom/bin/run_python_with_axom.sh \\
#     src/tools/gen-multidom-point-mesh.py plane --dimension 3 --grid-size 32,32 \\
#     --center 0,0,0 --normal 0,0,1 --extent 2,2 --domains 4 \\
#     --output dcp_object_plane
#
#   # 2D structured quad query mesh; DCP queries are evaluated over the vertices:
#   build-axom/bin/run_python_with_axom.sh \\
#     src/tools/gen-multidom-point-mesh.py grid --grid-size 100,100 \\
#     --min 0,0 --max 2,2 --domain-counts 2,1 --grid-topology structured \\
#     --output dcp_query_structured_quads
#
#   # 3D unstructured hex query mesh:
#   build-axom/bin/run_python_with_axom.sh \\
#     src/tools/gen-multidom-point-mesh.py grid --grid-size 20,20,15 \\
#     --min 0,0,0 --max 2,2,2 --domain-counts 2,2,1 \\
#     --grid-topology unstructured --output dcp_query_unstructured_hexes

# This script requires a conduit installation configured with python3 and hdf5.
# Make sure PYTHONPATH includes /path/to/conduit/install/python-modules,
# or use Axom's convenience script /path/to/axom_build_dir/bin/run_python_with_axom.sh
# that includes Conduit in PYTHONPATH.

import itertools
import math
import os
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, ArgumentTypeError

try:
    import conduit
    import conduit.blueprint
    import conduit.relay
except ModuleNotFoundError as e:
    print(f'{e}\nMake sure your PYTHONPATH includes /path/to/conduit/install/python-modules\n'
          'Conduit must be configured with python and hdf5.\n'
          'Alternatively, you can use the convenience script\n'
          '/path/to/axom_build_dir/bin/run_python_with_axom.sh\n'
          'that includes Conduit in PYTHONPATH.')
    sys.exit(-1)

try:
    import numpy as np
except ModuleNotFoundError as e:
    print(f'{e}\nThis script requires numpy.')
    sys.exit(-1)

AXES = 'xyz'
LOGICAL_AXES = 'ijk'


def get_mpi_context():
    '''Return MPI context when mpi4py and Conduit MPI relay are available.'''
    mpi_size_env = 1
    for var_name in ('OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_SIZE', 'SLURM_STEP_NUM_TASKS'):
        try:
            mpi_size_env = max(mpi_size_env, int(os.environ.get(var_name, '1')))
        except ValueError:
            pass

    if mpi_size_env <= 1:
        return {'enabled': False, 'comm': None, 'rank': 0, 'size': 1}

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


def add_common_options(parser):
    '''Add output and distribution options shared by all mesh generators.'''
    dist = parser.add_argument_group('distribution')
    dist.add_argument('-dc',
                      '--domain-counts',
                      type=i_c,
                      help='Domain counts in one chunk direction or each coordinate direction')
    dist.add_argument('--domains',
                      type=positive_int,
                      help='Total domain count for a one-dimensional chunk partition')
    dist.add_argument('--single-domain',
                      action='store_true',
                      help='Write a single-domain mesh instead of multidomain output')
    dist.add_argument('--use-list',
                      action='store_true',
                      help='Put multidomain domains in a list instead of a map')

    bp = parser.add_argument_group('blueprint names and output')
    bp.add_argument('--topology-name', default='mesh', help='Blueprint topology name')
    bp.add_argument('--coordset-name', default='coords', help='Blueprint coordset name')
    bp.add_argument('--id-field',
                    action='store_true',
                    help='Add a vertex-associated int64 global_id field')
    bp.add_argument('--no-verification',
                    action='store_true',
                    help='Omit analytic verification metadata')
    bp.add_argument('--protocol', default='hdf5', help='Conduit relay protocol for save_mesh')
    bp.add_argument('-o', '--output', default='mdmesh', help='Output file base name')
    bp.add_argument('-v', '--verbose', action='store_true', help='Print additional info')


def add_point_topology_options(parser):
    parser.add_argument('--topology-type',
                        choices=('unstructured', 'points'),
                        default='unstructured',
                        help='Point topology representation')


def add_center_radius_options(parser, fixed_dim):
    parser.add_argument('--center',
                        type=f_c,
                        help=f'Center coordinates; scalar or {fixed_dim} comma-separated values')
    parser.add_argument('--radius', type=float, default=1.0, help='Circle/sphere radius')


def add_min_max_options(parser):
    parser.add_argument('--min',
                        dest='lower',
                        type=f_c,
                        help='Lower coordinate bounds; scalar or one value per dimension')
    parser.add_argument('--max',
                        dest='upper',
                        type=f_c,
                        help='Upper coordinate bounds; scalar or one value per dimension')


def parse_args():
    parser = ArgumentParser(description='Write blueprint meshes for DCP tests.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(dest='shape', required=True)

    circle = subparsers.add_parser('circle',
                                   formatter_class=ArgumentDefaultsHelpFormatter,
                                   help='2D analytic circle point mesh')
    add_common_options(circle)
    add_point_topology_options(circle)
    add_center_radius_options(circle, 2)
    circle.add_argument('-n',
                        '--point-count',
                        '--long-point-count',
                        dest='point_count',
                        type=positive_int,
                        default=1024,
                        help='Number of points around the circle')
    circle.add_argument('--random-spacing',
                        action='store_true',
                        help='Use random angular spacing instead of uniform spacing')
    circle.add_argument('--seed', type=int, default=0, help='Random seed for random spacing')

    sphere = subparsers.add_parser('sphere',
                                   formatter_class=ArgumentDefaultsHelpFormatter,
                                   help='3D analytic sphere point mesh')
    add_common_options(sphere)
    add_point_topology_options(sphere)
    add_center_radius_options(sphere, 3)
    sphere.add_argument('-n',
                        '--point-count',
                        type=positive_int,
                        default=1024,
                        help='Number of Fibonacci sphere points when not using latitude sampling')
    sphere.add_argument('--long-point-count',
                        type=positive_int,
                        help='Number of longitudinal samples for latitude sampling')
    sphere.add_argument('--lat-point-count',
                        type=positive_int,
                        help='Number of latitudinal samples for latitude sampling')
    sphere.add_argument('--lat-range',
                        type=f_c,
                        default=(-90.0, 90.0),
                        help='Latitude range in degrees for latitude sampling')
    sphere.add_argument('--random-spacing',
                        action='store_true',
                        help='Use random longitudinal spacing with latitude sampling')
    sphere.add_argument('--seed', type=int, default=0, help='Random seed for random spacing')

    torus = subparsers.add_parser('torus',
                                  formatter_class=ArgumentDefaultsHelpFormatter,
                                  help='3D analytic torus point mesh')
    add_common_options(torus)
    add_point_topology_options(torus)
    torus.add_argument('-n',
                       '--point-count',
                       type=positive_int,
                       default=1024,
                       help='Number of torus surface points')
    torus.add_argument('--center', type=f_c, help='Center coordinates; scalar or 3 values')
    torus.add_argument('--major-radius', type=float, default=1.0, help='Torus major radius')
    torus.add_argument('--minor-radius', type=float, default=0.25, help='Torus minor radius')

    annulus = subparsers.add_parser('annulus',
                                    formatter_class=ArgumentDefaultsHelpFormatter,
                                    help='2D analytic annulus boundary point mesh')
    add_common_options(annulus)
    add_point_topology_options(annulus)
    annulus.add_argument('-n',
                         '--point-count',
                         type=positive_int,
                         default=1024,
                         help='Total number of points on the inner and outer circles')
    annulus.add_argument('--center', type=f_c, help='Center coordinates; scalar or 2 values')
    annulus.add_argument('--inner-radius', type=float, default=0.5, help='Inner radius')
    annulus.add_argument('--outer-radius', type=float, default=1.0, help='Outer radius')
    annulus.add_argument('--random-spacing',
                         action='store_true',
                         help='Use random angular spacing instead of uniform spacing')
    annulus.add_argument('--seed', type=int, default=0, help='Random seed for random spacing')

    plane = subparsers.add_parser('plane',
                                  formatter_class=ArgumentDefaultsHelpFormatter,
                                  help='2D line or 3D plane point mesh')
    add_common_options(plane)
    add_point_topology_options(plane)
    plane.add_argument('-d', '--dimension', type=int, choices=(2, 3), help='Spatial dimension')
    plane.add_argument('-n',
                       '--point-count',
                       type=positive_int,
                       default=1024,
                       help='Target point count when --grid-size is omitted')
    plane.add_argument('--grid-size',
                       type=i_c,
                       help='Point count along the line, or two counts on the plane')
    plane.add_argument('--center', type=f_c, help='Plane center; scalar or per dimension')
    plane.add_argument('--normal', type=f_c, help='Plane normal; scalar or per dimension')
    plane.add_argument('--extent', type=f_c, help='Line length in 2D, or two side lengths in 3D')

    grid = subparsers.add_parser('grid',
                                 formatter_class=ArgumentDefaultsHelpFormatter,
                                 help='Regular grid mesh over --min/--max')
    add_common_options(grid)
    add_min_max_options(grid)
    grid.add_argument('-d', '--dimension', type=int, choices=(2, 3), help='Spatial dimension')
    grid.add_argument('-n',
                      '--point-count',
                      type=positive_int,
                      default=1024,
                      help='Target point count when --grid-size is omitted')
    grid.add_argument('--grid-size',
                      type=i_c,
                      help='Point counts in each grid direction, e.g. 101,101')
    grid.add_argument('--grid-topology',
                      choices=('points', 'unstructured-points', 'structured', 'unstructured'),
                      default='unstructured-points',
                      help='Topology to generate over the regular grid coordinates')

    gaussian = subparsers.add_parser('gaussian',
                                     formatter_class=ArgumentDefaultsHelpFormatter,
                                     help='Random Gaussian point mesh')
    add_common_options(gaussian)
    add_point_topology_options(gaussian)
    gaussian.add_argument('-d', '--dimension', type=int, choices=(2, 3), help='Spatial dimension')
    gaussian.add_argument('-n',
                          '--point-count',
                          type=positive_int,
                          default=1024,
                          help='Number of points')
    gaussian.add_argument('--center', type=f_c, help='Distribution center; scalar or per dimension')
    gaussian.add_argument('--stddev',
                          type=f_c,
                          help='Gaussian standard deviation; scalar or per dimension')
    gaussian.add_argument('--seed', type=int, default=0, help='Random seed')

    uniform = subparsers.add_parser('uniform',
                                    formatter_class=ArgumentDefaultsHelpFormatter,
                                    help='Random uniform point mesh over --min/--max')
    add_common_options(uniform)
    add_point_topology_options(uniform)
    add_min_max_options(uniform)
    uniform.add_argument('-d', '--dimension', type=int, choices=(2, 3), help='Spatial dimension')
    uniform.add_argument('-n',
                         '--point-count',
                         type=positive_int,
                         default=1024,
                         help='Number of points')
    uniform.add_argument('--seed', type=int, default=0, help='Random seed')

    opts, unkn = parser.parse_known_args()
    if opts.verbose:
        print(opts, unkn)
    if unkn:
        print('Unrecognized arguments:', *unkn)
        sys.exit(1)

    return opts


def infer_dimension(opts):
    '''Infer and validate the spatial dimension.'''
    fixed_dims = {'circle': 2, 'sphere': 3, 'torus': 3, 'annulus': 2}
    if opts.shape in fixed_dims:
        return fixed_dims[opts.shape]

    dim = opts.dimension
    inferred = []
    names = ('grid_size', 'lower', 'upper', 'center', 'stddev', 'domain_counts')
    if opts.shape == 'plane':
        # Plane grid_size/extent are intrinsic coordinates, not spatial dimension.
        names = ('center', 'normal', 'domain_counts')
    for name in names:
        values = getattr(opts, name, None)
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
        if hasattr(default, '__len__'):
            default = np.array(default, dtype=float)
            if len(default) != dim:
                raise RuntimeError(f'default for {name} must have {dim} values')
            return default
        return np.full(dim, default, dtype=float)
    if len(values) == 1:
        return np.full(dim, values[0], dtype=float)
    if len(values) != dim:
        raise RuntimeError(f'{name} must have one value or {dim} values')
    return np.array(values, dtype=float)


def normalized(vec, name):
    '''Return a normalized copy of a vector.'''
    norm = np.linalg.norm(vec)
    if norm <= 0.0:
        raise RuntimeError(f'{name} must be nonzero')
    return vec / norm


def default_normal(dim):
    normal = np.zeros(dim, dtype=float)
    normal[-1] = 1.0
    return normal


def domain_counts_from_options(opts, dim, mpi_size):
    '''Return domain counts as an integer vector.'''
    if opts.single_domain:
        return np.ones(dim, dtype=int)

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
        base = max(1, int(round(opts.point_count ** (1.0 / dim))))
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

    chunks = []
    domain_indices = list(itertools.product(*[range(c) for c in domain_counts]))
    for domain_id, domain_index in enumerate(domain_indices):
        begin, end = split_range(point_count, domain_count, domain_id)
        chunks.append({
            'domain_id': domain_id,
            'domain_index': domain_index,
            'points': points[begin:end],
            'global_id_start': begin,
            'grid_counts': None,
        })
    return chunks


def domain_name(domain_index):
    '''Create a stable map key for a domain.'''
    if len(domain_index) == 1:
        return f'domain_{domain_index[0]:06d}'
    return 'domain_' + '_'.join(f'{idx:03d}' for idx in domain_index)


def generate_circle(point_count, center, radius, random_spacing, seed):
    if random_spacing:
        rng = np.random.default_rng(seed)
        theta = np.sort(rng.uniform(0.0, 2.0 * math.pi, point_count))
    else:
        theta = np.linspace(0.0, 2.0 * math.pi, point_count, endpoint=False)

    points = np.empty((point_count, 2), dtype=np.float64)
    points[:, 0] = center[0] + radius * np.cos(theta)
    points[:, 1] = center[1] + radius * np.sin(theta)
    return points


def generate_sphere(opts, center, radius):
    if opts.long_point_count is None and opts.lat_point_count is None:
        point_count = opts.point_count
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

    long_count = opts.long_point_count if opts.long_point_count is not None else opts.point_count
    lat_count = opts.lat_point_count if opts.lat_point_count is not None else 1
    if len(opts.lat_range) != 2:
        raise RuntimeError('--lat-range must have two values')

    min_lat = math.radians(opts.lat_range[0])
    max_lat = math.radians(opts.lat_range[1])
    lat_spacing = 0.0 if lat_count == 1 else (max_lat - min_lat) / (lat_count - 1)
    long_spacing = 2.0 * math.pi / long_count
    rng = np.random.default_rng(opts.seed)

    points = np.empty((long_count * lat_count, 3), dtype=np.float64)
    idx = 0
    for li in range(lat_count):
        lat = min_lat + li * lat_spacing
        xy_radius = radius * math.cos(lat)
        z = center[2] + radius * math.sin(lat)
        if opts.random_spacing:
            theta_values = np.sort(rng.uniform(0.0, 2.0 * math.pi, long_count))
        else:
            theta_values = np.arange(long_count, dtype=np.float64) * long_spacing
        for theta in theta_values:
            points[idx, 0] = center[0] + xy_radius * math.cos(theta)
            points[idx, 1] = center[1] + xy_radius * math.sin(theta)
            points[idx, 2] = z
            idx += 1

    return points


def generate_torus(point_count, center, major_radius, minor_radius):
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


def generate_annulus(point_count, center, inner_radius, outer_radius, random_spacing, seed):
    if point_count < 2:
        raise RuntimeError('annulus requires at least two points')

    inner_count = point_count // 2
    outer_count = point_count - inner_count
    rng = np.random.default_rng(seed)

    def ring_points(count, radius):
        if random_spacing:
            theta = np.sort(rng.uniform(0.0, 2.0 * math.pi, count))
        else:
            theta = np.linspace(0.0, 2.0 * math.pi, count, endpoint=False)
        ring = np.empty((count, 2), dtype=np.float64)
        ring[:, 0] = center[0] + radius * np.cos(theta)
        ring[:, 1] = center[1] + radius * np.sin(theta)
        return ring

    return np.vstack((ring_points(inner_count,
                                  inner_radius), ring_points(outer_count, outer_radius)))


def plane_basis(normal):
    '''Return one tangent in 2D or two orthonormal tangents in 3D.'''
    if len(normal) == 2:
        return (np.array([-normal[1], normal[0]], dtype=np.float64), )

    ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if abs(np.dot(ref, normal)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    tangent0 = normalized(np.cross(normal, ref), '--normal')
    tangent1 = normalized(np.cross(normal, tangent0), '--normal')
    return tangent0, tangent1


def plane_grid_counts(opts, dim):
    if opts.grid_size is not None:
        counts = np.array(opts.grid_size, dtype=int)
        expected = 1 if dim == 2 else 2
        if len(counts) != expected:
            raise RuntimeError(f'--grid-size for a {dim}D plane must have {expected} value(s)')
        if np.any(counts < 1):
            raise RuntimeError('plane grid sizes must be positive')
        return counts

    if dim == 2:
        return np.array([opts.point_count], dtype=int)

    base = max(1, int(round(math.sqrt(opts.point_count))))
    counts = np.array([base, base], dtype=int)
    while int(np.prod(counts)) < opts.point_count:
        counts[np.argmin(counts)] += 1
    return counts


def generate_plane(opts, dim, center, normal):
    counts = plane_grid_counts(opts, dim)
    extent = vector_option(getattr(opts, 'extent', None), len(counts), 2.0, '--extent')
    if np.any(extent <= 0.0):
        raise RuntimeError('--extent values must be positive')

    basis = plane_basis(normal)
    if dim == 2:
        t = np.linspace(-0.5 * extent[0], 0.5 * extent[0], counts[0])
        return center + np.outer(t, basis[0])

    axes = [np.linspace(-0.5 * extent[d], 0.5 * extent[d], counts[d]) for d in range(2)]
    coords = np.meshgrid(*axes, indexing='ij')
    weights = np.column_stack([coord.ravel(order='F') for coord in coords])
    points = np.empty((weights.shape[0], 3), dtype=np.float64)
    for i, weight in enumerate(weights):
        points[i] = center + weight[0] * basis[0] + weight[1] * basis[1]
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
    return np.column_stack([coord.ravel(order='F') for coord in coords])


def generate_grid_domains(grid_counts, lower, upper, domain_counts):
    '''Generate regular grid points one spatial domain at a time.'''
    axes = [np.linspace(lower[d], upper[d], grid_counts[d]) for d in range(len(grid_counts))]
    chunks = []
    global_start = 0
    for domain_id, domain_index in enumerate(itertools.product(*[range(c) for c in domain_counts])):
        local_axes = []
        for d, idx in enumerate(domain_index):
            begin, end = split_range(grid_counts[d], domain_counts[d], idx)
            local_axes.append(axes[d][begin:end])

        local_counts = np.array([len(axis) for axis in local_axes], dtype=int)
        if np.any(local_counts == 0):
            points = np.empty((0, len(grid_counts)), dtype=np.float64)
        else:
            coords = np.meshgrid(*local_axes, indexing='ij')
            points = np.column_stack([coord.ravel(order='F') for coord in coords])

        chunks.append({
            'domain_id': domain_id,
            'domain_index': domain_index,
            'points': points,
            'global_id_start': global_start,
            'grid_counts': local_counts,
        })
        global_start += points.shape[0]

    return chunks


def generated_domains(opts, dim, lower, upper, center, normal, stddev, domain_counts):
    '''Return generated domain chunk dictionaries.'''
    if opts.shape == 'grid':
        grid_counts = grid_counts_from_options(opts, dim)
        if opts.verbose:
            print(f'grid_size={grid_counts.tolist()} point_count={int(np.prod(grid_counts))}')

        if opts.grid_topology in ('structured', 'unstructured') and len(domain_counts) != dim:
            raise RuntimeError(
                'structured and unstructured grid topologies require --domain-counts with '
                f'{dim} values')

        if len(domain_counts) == dim:
            return generate_grid_domains(grid_counts, lower, upper, domain_counts)

        points = generate_grid_points(grid_counts, lower, upper)
        return split_points(points, domain_counts)

    if opts.shape == 'circle':
        points = generate_circle(opts.point_count, center, opts.radius, opts.random_spacing,
                                 opts.seed)
    elif opts.shape == 'sphere':
        points = generate_sphere(opts, center, opts.radius)
    elif opts.shape == 'torus':
        points = generate_torus(opts.point_count, center, opts.major_radius, opts.minor_radius)
    elif opts.shape == 'annulus':
        points = generate_annulus(opts.point_count, center, opts.inner_radius, opts.outer_radius,
                                  opts.random_spacing, opts.seed)
    elif opts.shape == 'plane':
        points = generate_plane(opts, dim, center, normal)
    elif opts.shape == 'gaussian':
        points = generate_gaussian(opts.point_count, center, stddev, opts.seed)
    elif opts.shape == 'uniform':
        points = generate_uniform(opts.point_count, lower, upper, opts.seed)
    else:
        raise RuntimeError(f'Unknown shape: {opts.shape}')

    return split_points(points, domain_counts)


def add_unstructured_point_topology(topo, point_count):
    topo['type'] = 'unstructured'
    topo['elements/shape'] = 'point'
    topo['elements/connectivity'].set(np.arange(point_count, dtype=np.int32))
    topo['elements/sizes'].set(np.ones(point_count, dtype=np.int32))
    topo['elements/offsets'].set(np.arange(point_count, dtype=np.int32))


def point_index(i, j, k, counts):
    return i + counts[0] * (j + counts[1] * k)


def grid_connectivity(counts):
    '''Return unstructured quad/hex connectivity for a regular grid domain.'''
    dim = len(counts)
    if np.any(counts < 2):
        return np.empty(0, dtype=np.int32), 'quad' if dim == 2 else 'hex'

    conn = []
    if dim == 2:
        for j in range(counts[1] - 1):
            for i in range(counts[0] - 1):
                conn.extend([
                    point_index(i, j, 0, counts),
                    point_index(i + 1, j, 0, counts),
                    point_index(i + 1, j + 1, 0, counts),
                    point_index(i, j + 1, 0, counts),
                ])
        return np.array(conn, dtype=np.int32), 'quad'

    for k in range(counts[2] - 1):
        for j in range(counts[1] - 1):
            for i in range(counts[0] - 1):
                conn.extend([
                    point_index(i, j, k, counts),
                    point_index(i + 1, j, k, counts),
                    point_index(i + 1, j + 1, k, counts),
                    point_index(i, j + 1, k, counts),
                    point_index(i, j, k + 1, counts),
                    point_index(i + 1, j, k + 1, counts),
                    point_index(i + 1, j + 1, k + 1, counts),
                    point_index(i, j + 1, k + 1, counts),
                ])
    return np.array(conn, dtype=np.int32), 'hex'


def format_vector(values):
    return '(' + ', '.join(f'{value:g}' for value in values) + ')'


def shape_description(opts, dim, lower, upper, center, normal, stddev):
    if opts.shape == 'circle':
        return f'2D circle centered at {format_vector(center)} with radius {opts.radius:g}'
    if opts.shape == 'sphere':
        return f'3D sphere centered at {format_vector(center)} with radius {opts.radius:g}'
    if opts.shape == 'torus':
        return (f'3D torus centered at {format_vector(center)} with major radius '
                f'{opts.major_radius:g} and minor radius {opts.minor_radius:g}')
    if opts.shape == 'annulus':
        return (f'2D annulus centered at {format_vector(center)} with inner radius '
                f'{opts.inner_radius:g} and outer radius {opts.outer_radius:g}')
    if opts.shape == 'plane':
        return f'{dim}D codimension-1 plane centered at {format_vector(center)}'
    if opts.shape == 'grid':
        return f'{dim}D regular grid over [{format_vector(lower)}, {format_vector(upper)}]'
    if opts.shape == 'gaussian':
        return (f'{dim}D Gaussian point cloud centered at {format_vector(center)} with stddev '
                f'{format_vector(stddev)}')
    if opts.shape == 'uniform':
        return f'{dim}D uniform random point cloud over [{format_vector(lower)}, {format_vector(upper)}]'
    return f'{dim}D {opts.shape} point mesh'


def circle_spacing(radius, point_count):
    return 2.0 * math.pi * radius / max(point_count, 1)


def sphere_spacing(opts, radius):
    if opts.long_point_count is not None or opts.lat_point_count is not None:
        long_count = opts.long_point_count if opts.long_point_count is not None else opts.point_count
        lat_count = opts.lat_point_count if opts.lat_point_count is not None else 1
        long_spacing = 2.0 * math.pi * radius / max(long_count, 1)
        lat_spacing = 0.0
        if lat_count > 1:
            lat_spacing = radius * math.radians(opts.lat_range[1] -
                                                opts.lat_range[0]) / (lat_count - 1)
        return math.sqrt(long_spacing * long_spacing + lat_spacing * lat_spacing)

    return radius * math.sqrt(4.0 * math.pi / max(opts.point_count, 1))


def torus_spacing(major_radius, minor_radius, point_count):
    area = 4.0 * math.pi * math.pi * major_radius * minor_radius
    return 2.0 * math.sqrt(area / max(point_count, 1))


def annulus_spacing(opts):
    inner_count = max(opts.point_count // 2, 1)
    outer_count = max(opts.point_count - inner_count, 1)
    return max(circle_spacing(opts.inner_radius, inner_count),
               circle_spacing(opts.outer_radius, outer_count))


def plane_spacing(opts, dim):
    counts = plane_grid_counts(opts, dim)
    extent = vector_option(getattr(opts, 'extent', None), len(counts), 2.0, '--extent')
    spacing = []
    for d, count in enumerate(counts):
        spacing.append(extent[d] / max(count - 1, 1))
    return max(spacing)


def verification_metadata(opts, dim, center, normal):
    if opts.no_verification:
        return None

    metadata = {'shape': opts.shape, 'dimension': dim, 'surface_tolerance': 1.0e-8}
    if opts.shape == 'circle':
        metadata.update({
            'center': center,
            'radius': opts.radius,
            'distance_tolerance': circle_spacing(opts.radius, opts.point_count),
        })
    elif opts.shape == 'sphere':
        metadata.update({
            'center': center,
            'radius': opts.radius,
            'distance_tolerance': sphere_spacing(opts, opts.radius),
        })
    elif opts.shape == 'torus':
        metadata.update({
            'center':
            center,
            'major_radius':
            opts.major_radius,
            'minor_radius':
            opts.minor_radius,
            'distance_tolerance':
            torus_spacing(opts.major_radius, opts.minor_radius, opts.point_count),
        })
    elif opts.shape == 'annulus':
        metadata.update({
            'center': center,
            'inner_radius': opts.inner_radius,
            'outer_radius': opts.outer_radius,
            'distance_tolerance': annulus_spacing(opts),
        })
    elif opts.shape == 'plane':
        metadata.update({
            'center': center,
            'normal': normal,
            'distance_tolerance': plane_spacing(opts, dim),
        })
    else:
        return None

    return metadata


def add_metadata(dom, metadata):
    dom['state/description'] = metadata['description']
    verification = metadata['verification']
    if verification is None:
        return

    ver_node = dom['state/verification']
    for key, value in verification.items():
        if isinstance(value, np.ndarray):
            ver_node[key].set(np.ascontiguousarray(value, dtype=np.float64))
        else:
            ver_node[key] = value


def fill_domain(dom, opts, chunk, metadata):
    '''Fill one blueprint domain with mesh data.'''
    points = chunk['points']
    point_count = points.shape[0]
    dim = points.shape[1]

    dom['state/domain_id'] = int(chunk['domain_id'])
    add_metadata(dom, metadata)
    dom[f'coordsets/{opts.coordset_name}/type'] = 'explicit'
    for d in range(dim):
        values_path = f'coordsets/{opts.coordset_name}/values/{AXES[d]}'
        dom[values_path].set(np.ascontiguousarray(points[:, d], dtype=np.float64))

    topo = dom[f'topologies/{opts.topology_name}']
    topo['coordset'] = opts.coordset_name

    if opts.shape == 'grid':
        grid_topology = opts.grid_topology
        if grid_topology == 'points':
            topo['type'] = 'points'
        elif grid_topology == 'unstructured-points':
            add_unstructured_point_topology(topo, point_count)
        elif grid_topology == 'structured':
            topo['type'] = 'structured'
            counts = chunk['grid_counts']
            for d in range(dim):
                topo[f'elements/dims/{LOGICAL_AXES[d]}'] = max(int(counts[d]) - 1, 0)
        elif grid_topology == 'unstructured':
            topo['type'] = 'unstructured'
            counts = chunk['grid_counts']
            conn, shape = grid_connectivity(counts)
            topo['elements/shape'] = shape
            topo['elements/connectivity'].set(conn)
        else:
            raise RuntimeError(f'Unsupported grid topology: {grid_topology}')
    elif opts.topology_type == 'points':
        topo['type'] = 'points'
    else:
        add_unstructured_point_topology(topo, point_count)

    if opts.id_field:
        field = dom['fields/global_id']
        field['association'] = 'vertex'
        field['topology'] = opts.topology_name
        field['values'].set(
            np.arange(chunk['global_id_start'],
                      chunk['global_id_start'] + point_count,
                      dtype=np.int64))


def create_domain(md_mesh, opts, chunk, metadata):
    '''Append one domain to the multidomain mesh.'''
    if opts.use_list:
        dom = md_mesh.append()
    else:
        dom = md_mesh[domain_name(chunk['domain_index'])]

    fill_domain(dom, opts, chunk, metadata)


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
    lower = vector_option(getattr(opts, 'lower', None), dim, -1.0, '--min')
    upper = vector_option(getattr(opts, 'upper', None), dim, 1.0, '--max')
    center = vector_option(getattr(opts, 'center', None), dim, 0.0, '--center')
    normal = vector_option(getattr(opts, 'normal', None), dim, default_normal(dim), '--normal')
    if opts.shape == 'plane':
        normal = normalized(normal, '--normal')
    stddev = vector_option(getattr(opts, 'stddev', None), dim, 1.0, '--stddev')
    domain_counts = domain_counts_from_options(opts, dim, mpi['size'])

    if np.any(upper <= lower):
        raise RuntimeError('--max values must be greater than --min values')
    if hasattr(opts, 'radius') and opts.radius <= 0.0:
        raise RuntimeError('--radius must be positive')
    if hasattr(opts, 'major_radius') and (opts.major_radius <= 0.0 or opts.minor_radius <= 0.0):
        raise RuntimeError('--major-radius and --minor-radius must be positive')
    if hasattr(opts, 'inner_radius') and (opts.inner_radius <= 0.0
                                          or opts.outer_radius <= opts.inner_radius):
        raise RuntimeError('--outer-radius must be greater than positive --inner-radius')
    if np.any(stddev <= 0.0):
        raise RuntimeError('--stddev values must be positive')

    metadata = {
        'description': shape_description(opts, dim, lower, upper, center, normal, stddev),
        'verification': verification_metadata(opts, dim, center, normal),
    }

    chunks = generated_domains(opts, dim, lower, upper, center, normal, stddev, domain_counts)
    local_chunks = local_chunks_for_rank(chunks, mpi['rank'], mpi['size'], opts.single_domain)

    mesh = conduit.Node()
    if opts.single_domain:
        if local_chunks:
            fill_domain(mesh, opts, local_chunks[0], metadata)
    else:
        for chunk in local_chunks:
            create_domain(mesh, opts, chunk, metadata)

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
        total_points = sum(chunk['points'].shape[0] for chunk in chunks)
        mesh_kind = 'single-domain' if opts.single_domain else 'multidomain'
        print(f'Wrote {total_points} {dim}D {opts.shape} points as a {mesh_kind} mesh '
              f'with {len(chunks)} domain(s) to {opts.output} using {opts.protocol}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
