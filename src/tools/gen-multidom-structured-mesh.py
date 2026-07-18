#!/usr/bin/env python3

# gen-multidom-structured-mesh.py
# Write a simple multidomain structured blueprint mesh for testing.

# This script requires a conduit installation configured with python3 and hdf5.
# Make sure PYTHONPATH includes /path/to/conduit/install/python-modules,
# or use Axom's convenience script /path/to/axom_build_dir/bin/run_python_with_axom.sh
# that includes Conduit in PYTHONPATH.

try:
    import conduit
    import conduit.blueprint
    import conduit.relay
except ModuleNotFoundError as e:
    print(
        f'{e}\nMake sure your PYTHONPATH includes /path/to/conduit/install/python-modules\nConduit must be configured with python and hdf5.\nAlternatively, you can use the convenience script\n/path/to/axom_build_dir/bin/run_python_with_axom.sh\nthat includes Conduit in PYTHONPATH.'
    )
    exit(-1)

import numpy as np


def i_c(s):
    '''Convert comma-separated string to list of integers.'''
    return list(map(int, s.split(',')))


def f_c(s):
    '''Convert comma-separated string to list of floating point numbers.'''
    return list(map(float, s.split(',')))


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

ps = ArgumentParser(description='Write a blueprint strided-unstructured mesh.',
                    formatter_class=ArgumentDefaultsHelpFormatter)
ps.add_argument('--useList', action='store_true', help='Put domains in a list instead of a map')
ps.add_argument('-ml', type=f_c, default=(0., 0.), help='Mesh lower coordinates')
ps.add_argument('-mu', type=f_c, default=(1., 1.), help='Mesh upper coordinates')
ps.add_argument('-ms', type=i_c, default=(3, 3), help='Logical size of mesh (cells)')
ps.add_argument('-dc', type=i_c, default=(1, 1), help='Domain counts in each index direction')
ps.add_argument('-o', '--output', type=str, default='mdmesh', help='Output file base name')
ps.add_argument('--strided', action='store_true', help='Use strided_structured (has ghosts)')
ps.add_argument('-v', '--verbose', action='store_true', help='Print additional info')
opts, unkn = ps.parse_known_args()
if (opts.verbose): print(opts, unkn)
if (unkn):
    print("Unrecognized arguments:", *unkn)
    quit(1)

dim = len(opts.dc)

if dim not in (2,3)	\
   or len(opts.ms) != dim	\
   or len(opts.ml) != dim	\
   or len(opts.mu) != dim:
    raise RuntimeError('dc, ms, ml and mu options must have the same dimensions (2 or 3)')

# Must have enough cells for requested partitioning.
goodDc = [opts.ms[i] >= opts.dc[i] for i in range(dim)]
if sum(goodDc) < dim:
    raise RuntimeError(f'ms ({opts.ms}) must be >= dc ({opts.dc}) in all directions.')

# Number of phony nodes on left and right sides, for strided option
if opts.strided:
    npnl, npnr = 2, 1
else:
    npnl, npnr = 0, 0


def scale_structured_domain(n, startCoord, endCoord):
    '''This function scales and shifts a blueprint structured domain after
  it has been created.  There's no way to specify the physical extent
  of a domain using conduit.blueprint.mesh.examples methods, as far as
  I can tell.
  '''
    #print(f'Rescaling to {startCoord} -> {endCoord}')

    ndim = n['coordsets/coords/values'].number_of_children()

    domLens = n['topologies/mesh/elements/dims']
    dirs = 'ij' if ndim == 2 else 'ijk'
    domLens = [domLens[d] for d in dirs]
    domLens = np.array(domLens)
    domPhysicalSize = np.array(endCoord) - np.array(startCoord)
    #print(f'domLens={domLens}   domPhysicalSize={domPhysicalSize}')
    assert (n['topologies/mesh/type'] == 'structured')
    assert (len(startCoord) >= ndim)
    assert (len(domPhysicalSize) >= ndim)

    coordArrayLens = domLens + 1 + npnl + npnr
    #print(f'coordArrayLens={coordArrayLens}')

    xyz = 'xyz'
    for d in range(ndim):
        coords = n['coordsets/coords/values'][d]
        coords = np.reshape(coords, np.flip(coordArrayLens))

        # realCoords excludes the ghost layers.
        if ndim == 2:
            if npnr == 0:
                realCoords = coords[npnl:, npnl:]
            else:
                realCoords = coords[npnl:-npnr, npnl:-npnr]
        else:
            if npnr == 0:
                realCoords = coords[npnl:, npnl:, npnl:]
            else:
                realCoords = coords[npnl:-npnr, npnl:-npnr, npnl:-npnr]

        minC, maxC = np.amin(realCoords), np.amax(realCoords)
        curRange = maxC - minC
        shift = startCoord[d] - minC
        scale = domPhysicalSize[d] / curRange
        coords = (coords - minC) * domPhysicalSize[d] / curRange + startCoord[d]
        n['coordsets/coords/values'][xyz[d]] = coords


domType = 'structured'

domCounts = opts.dc if dim == 3 else (*opts.dc, 1)  # domCounts must be length 3, even for 2D.
meshSize = opts.ms
meshLower = opts.ml
meshUpper = opts.mu

# Convert to np.array to use element-wise arithmetic.
domCounts = np.array(domCounts, dtype=int)
meshSize = np.array(meshSize, dtype=int)
meshLower = np.array(meshLower)
meshUpper = np.array(meshUpper)

domPhysicalSize = (meshUpper - meshLower) / domCounts[:dim]
cellPhysicalSize = (meshUpper - meshLower) / meshSize

domSize = meshSize // domCounts[:dim]
domSizeRem = meshSize % domCounts[:dim]
if opts.verbose:
    print(f'meshSize={meshSize} cells, domCounts={domCounts[0:dim]}'
          f' domSize={domSize} domSizeRem={domSizeRem}')


def domain_index_begin(di, dj, dk=None):
    '''Compute first cell index of the domain with multi-dimensional index (di, dj, dk).'''
    ds = (di, dj) if dim == 2 else (di, dj, dk)
    idx = np.array(ds)
    std = domSize * ds
    extra = np.where(idx < domSizeRem[:dim], idx, domSizeRem[:dim])
    begin = std + extra
    return begin


mdMesh = conduit.Node()
for dk in range(domCounts[2]):
    for dj in range(domCounts[1]):
        for di in range(domCounts[0]):
            if opts.useList:
                dom = mdMesh.append()
            else:
                domName = f'domain_{di:1d}_{dj:1d}'
                if len(opts.dc) == 3: domName += f'_{dk:1d}'
                dom = mdMesh[domName]

            cellStart = domain_index_begin(di, dj, dk)
            cellEnd = domain_index_begin(di + 1, dj + 1, dk + 1 if dim == 3 else 0)
            pointCounts = cellEnd - cellStart + 1
            #print(f'cellStart={cellStart}  cellEnd={cellEnd}  pointCounts={pointCounts}')

            elemExtents = (cellEnd - cellStart) + (npnl + npnr + 1)
            vertExtents = np.array(pointCounts) + (npnl + npnr)
            elemOffset = np.full(dim, npnl)
            vertOffset = np.full(dim, npnl)
            #print(f'\n{domName}: {cellStart} -> {cellEnd}')

            pointCounts3 = pointCounts if len(pointCounts) == 3 else (*pointCounts, 0)
            if opts.strided:
                desc = conduit.Node()
                desc['vertex_data/shape'].set(vertExtents)
                desc['vertex_data/origin'].set(vertOffset)
                desc['element_data/shape'].set(elemExtents)
                desc['element_data/origin'].set(elemOffset)
                #print(f'\ndesc({di},{dj},{dk}):', end=''); print(desc)
                conduit.blueprint.mesh.examples.strided_structured(desc, *pointCounts3, dom)
                if dom.has_child("state"): dom.remove_child("state")
            else:
                conduit.blueprint.mesh.examples.basic(domType, *pointCounts3, dom)

            domLower = meshLower[:dim] + cellStart * cellPhysicalSize[:dim]
            domUpper = meshLower[:dim] + cellEnd * cellPhysicalSize[:dim]
            scale_structured_domain(dom, domLower, domUpper)
            # if opts.verbose: print(f'Domain [{di},{dj},{dk}]: {dom}')

if opts.verbose:
    print('mdMesh:')
    print(mdMesh)

info = conduit.Node()
if not conduit.blueprint.mesh.verify(mdMesh, info):
    print("Mesh failed blueprint verification.  Info:")
    print(info)

conduit.relay.io.blueprint.save_mesh(mdMesh, opts.output, "hdf5")
print(f'Wrote mesh {opts.output}')
