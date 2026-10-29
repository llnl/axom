#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Generate Hilbert finite state machine (FSM) tables for Axom's Skilling convention.

Each transition reads one bit plane from the input coordinates.
The state holds an axis permutation, axis flips, and a Gray correction parity bit.
Prefix XOR of the oriented input bits gives the Gray digit; XOR with
the saved parity gives the Hilbert digit. Skilling's swaps and flips
update the orientation of lower planes. The final Gray bit, before correction,
updates the saved parity.

The reachable machines have 8, 48, and 384 states in 2D, 3D, and 4D.
Partition refinement merges states whose outputs match and whose next states
belong to the same partition for every input digit. This leaves 4, 24, and 192 states.
Grouping composes one-plane transitions and keeps the same state numbering.

The 2D tables use the minimized machine. The 3D/4D tables store only axis
permutation rows, with reflections held in a separate source-axis XOR mask.
Each entry returns output digits and an XOR update to the packed state.
This needs 6 rows in 3D and 24 in 4D. The emitted tables process groups of
1/2/3 planes in 2D, 1/2 in 3D, and 1 in 4D.

Input digits have x in the least-significant bit, matching Axom MortonIndex.
Output digits have axis zero most significant, matching the scalar reference.
State zero is identity at every depth. Do not pad a short leading group with zero planes
since a zero plane changes the orientation.

The transition rules follow Skilling's algorithm appendix:
https://doi.org/10.1063/1.1751381
"""
import argparse
import itertools
from pathlib import Path
import random

PRODUCTION_LEVELS = {2: (1, 2, 3), 3: (1, 2), 4: (1, )}
RAW_STATE_COUNTS = {2: 8, 3: 48, 4: 384}
MINIMIZED_STATE_COUNTS = {2: 4, 3: 24, 4: 192}
VERIFICATION_SEED = 42
RANDOM_SAMPLES_PER_DEPTH = 1000


def require(condition, message):
    """Keep generation and verification checks active with python -O."""
    if not condition:
        raise ValueError(message)


def scalar_reference(point, bits):
    """Skilling's axes/transpose transform with explicit bit packing.

    This is separate from the state machine, but uses the same curve algorithm
    as Axom's scalar implementation. Convention fixtures remain in the C++ tests.
    """
    axes = list(point)
    dim = len(axes)
    for bit in range(bits - 1, 0, -1):
        q = 1 << bit
        p = q - 1
        for axis in range(dim):
            if axes[axis] & q:
                axes[0] ^= p
            else:
                exchange = (axes[0] ^ axes[axis]) & p
                axes[0] ^= exchange
                axes[axis] ^= exchange
    for axis in range(1, dim):
        axes[axis] ^= axes[axis - 1]
    correction = 0
    for bit in range(bits - 1, 0, -1):
        if axes[-1] & (1 << bit):
            correction ^= (1 << bit) - 1
    axes = [axis ^ correction for axis in axes]
    key = 0
    for bit in range(bits - 1, -1, -1):
        for axis in axes:
            key = (key << 1) | ((axis >> bit) & 1)
    return key


def transition(state, digit, dim):
    """Emit one Hilbert digit and update the orientation of lower input planes."""
    permutation, flips, parity = state
    plane = [((digit >> permutation[a]) ^ (flips >> a)) & 1 for a in range(dim)]
    gray = output = 0
    for value in plane:
        gray ^= value
        output = (output << 1) | (gray ^ parity)
    next_permutation = list(permutation)
    next_flips = [(flips >> a) & 1 for a in range(dim)]
    for axis in range(dim):
        if plane[axis]:
            next_flips[0] ^= 1
        else:
            next_permutation[0], next_permutation[axis] = (next_permutation[axis],
                                                           next_permutation[0])
            next_flips[0], next_flips[axis] = next_flips[axis], next_flips[0]
    next_state = (tuple(next_permutation), sum(value << a for a, value in enumerate(next_flips)),
                  parity ^ gray)
    return output, next_state


def machine(dim):
    initial = (tuple(range(dim)), 0, 0)
    states = [initial]
    ids = {initial: 0}
    rows = []
    for state in states:  # breadth-first closure, deterministic digit order
        row = []
        for digit in range(1 << dim):
            output, next_state = transition(state, digit, dim)
            if next_state not in ids:
                ids[next_state] = len(states)
                states.append(next_state)
            row.append((output, ids[next_state]))
        rows.append(row)
    require(len(states) == RAW_STATE_COUNTS[dim], f'{dim}D raw state count')
    partitions = [0] * len(states)
    while True:
        signatures = [tuple((out, partitions[dest]) for out, dest in row) for row in rows]
        groups = {signature: i for i, signature in enumerate(dict.fromkeys(signatures))}
        refined = [groups[signature] for signature in signatures]
        if refined == partitions:
            break
        partitions = refined
    minimized = [None] * len(set(partitions))
    for state_id, row in enumerate(rows):
        mapped = [(out, partitions[dest]) for out, dest in row]
        block = partitions[state_id]
        require(minimized[block] is None or minimized[block] == mapped,
                f'{dim}D partition has inconsistent transitions')
        minimized[block] = mapped
    require(partitions[0] == 0, f'{dim}D identity must be state zero')
    require(len(minimized) == MINIMIZED_STATE_COUNTS[dim], f'{dim}D minimized state count')
    return minimized


def grouped(rows, dim, levels):
    """Pack each output group below its absolute next-state index."""
    result = []
    for state in range(len(rows)):
        row = []
        for digits in range(1 << (dim * levels)):
            current = state
            output = 0
            for level in range(levels - 1, -1, -1):
                out, current = rows[current][(digits >> (level * dim)) & ((1 << dim) - 1)]
                output = (output << dim) | out
            row.append((current << (dim * levels)) | output)
        result.extend(row)
    return result


def inverse(table, dim, levels):
    """Invert each output row, preserving its absolute state or XOR-update bits."""
    group_bits = dim * levels
    row_size = 1 << group_bits
    mask = row_size - 1
    result = [None] * len(table)
    for base in range(0, len(table), row_size):
        for digit, entry in enumerate(table[base:base + row_size]):
            output, state_update = entry & mask, entry >> group_bits
            require(result[base + output] is None, 'Output row must be bijective')
            result[base + output] = (state_update << group_bits) | digit
    require(all(entry is not None for entry in result), 'Inverse row is incomplete')
    return result


def factored(dim, levels):
    """Keep source-axis reflections in a mask, outside the permutation table.

    XOR the incoming mask into every input plane before lookup. The table needs
    only DIM! rows, one per permutation with mask zero, and returns an XOR state
    update. Fold the Gray correction parity into the next orientation's axis-zero
    flip, which complements every prefix-XOR output bit. Store the resulting
    flips in source-axis order for the next lookup.
    """
    permutations = list(itertools.permutations(range(dim)))
    ids = {permutation: i for i, permutation in enumerate(permutations)}
    width = dim * levels
    repeat = ((1 << width) - 1) // ((1 << dim) - 1)
    table = []
    for permutation in permutations:
        for digits in range(1 << width):
            state = (permutation, 0, 0)
            output = 0
            for level in range(levels - 1, -1, -1):
                out, state = transition(state, (digits >> (dim * level)) & ((1 << dim) - 1), dim)
                output = (output << dim) | out
            next_permutation, flips, parity = state
            delta = sum((((flips >> a) & 1) ^ (parity if a == 0 else 0)) << next_permutation[a]
                        for a in range(dim))
            state_delta = ((ids[permutation] ^ ids[next_permutation]) << width) | delta * repeat
            table.append(output | (state_delta << width))
    return table


def verify_factored(dim):
    """Check every permutation, reflection mask, and input group against transition().

    Matching outputs and next orientations establish equivalence for streams
    of any length, including partial leading groups.
    """
    permutations = list(itertools.permutations(range(dim)))
    ids = {permutation: i for i, permutation in enumerate(permutations)}
    count = 0
    for levels in PRODUCTION_LEVELS[dim]:
        width = dim * levels
        mask = (1 << width) - 1
        repeat = mask // ((1 << dim) - 1)
        forward = factored(dim, levels)
        backward = inverse(forward, dim, levels)
        for p, permutation in enumerate(permutations):
            for flips in range(1 << dim):
                for digits in range(1 << width):
                    state = (p << width) | flips * repeat
                    entry = forward[state ^ digits]
                    next_state = state ^ (entry >> width)
                    output = entry & mask
                    inv = backward[(state & ~mask) | output]
                    context = f'{dim}D/{levels} planes, permutation {p}, flips {flips}, digit {digits}'
                    require(((inv ^ state) & mask) == digits, f'Inverse output: {context}')
                    require(state ^ (inv >> width) == next_state, f'Inverse state: {context}')
                    # Apply the one-plane rule with the incoming reflections.
                    axes_flips = sum(((flips >> permutation[a]) & 1) << a for a in range(dim))
                    raw = (permutation, axes_flips, 0)
                    expected = 0
                    for level in range(levels - 1, -1, -1):
                        out, raw = transition(raw, (digits >> (dim * level)) & ((1 << dim) - 1),
                                              dim)
                        expected = (expected << dim) | out
                    perm, axis_flips, parity = raw
                    source_flips = sum(
                        (((axis_flips >> a) & 1) ^ (parity if a == 0 else 0)) << perm[a]
                        for a in range(dim))
                    require(output == expected, f'Factored output: {context}')
                    require(next_state == (ids[perm] << width) | source_flips * repeat,
                            f'Factored state: {context}')
                    count += 1
    return count


def convert(key, bits, dim, levels, tables):
    """MSB-first grouped conversion with an exact leading partial group."""
    state = result = 0
    remaining = bits
    while remaining:
        group = remaining % levels or levels
        remaining -= group
        width = dim * group
        mask = (1 << width) - 1
        entry = tables[group][(state << width) | ((key >> (remaining * dim)) & mask)]
        result = (result << width) | (entry & mask)
        state = entry >> width
    return result


def convert_factored(key, bits, dim, levels, tables, inverse_direction=False):
    """Exercise factored tables, repacking reflections for each group's width."""
    permutation = flips = result = 0
    remaining = bits
    plane_mask = (1 << dim) - 1
    while remaining:
        group = remaining % levels or levels
        remaining -= group
        width = dim * group
        mask = (1 << width) - 1
        state = (permutation << width) | flips * (mask // plane_mask)
        digits = (key >> (remaining * dim)) & mask
        offset = ((state & ~mask) | digits) if inverse_direction else (state ^ digits)
        entry = tables[group][offset]
        output = ((entry ^ state) if inverse_direction else entry) & mask
        result = (result << width) | output
        state ^= entry >> width
        permutation, flips = state >> width, state & plane_mask
    return result


def verify(dim):
    """Check the emitted table representations against the scalar transform."""
    count = 0
    groups = PRODUCTION_LEVELS[dim]
    rows = machine(dim) if dim == 2 else None
    forward = {
        levels: grouped(rows, dim, levels) if dim == 2 else factored(dim, levels)
        for levels in groups
    }
    backward = {levels: inverse(table, dim, levels) for levels, table in forward.items()}

    def check(point, bits):
        expected = scalar_reference(point, bits)
        morton = sum(((value >> bit) & 1) << (dim * bit + axis) for axis, value in enumerate(point)
                     for bit in range(bits))
        context = f'{dim}D, depth {bits}, point {point}'
        for levels in groups:
            if dim == 2:
                encoded = convert(morton, bits, dim, levels, forward)
                decoded = convert(expected, bits, dim, levels, backward)
            else:
                encoded = convert_factored(morton, bits, dim, levels, forward)
                decoded = convert_factored(expected, bits, dim, levels, backward, True)
            require(encoded == expected, f'{levels}-plane encode: {context}')
            require(decoded == morton, f'{levels}-plane decode: {context}')

    # Exhaust 4D through depth 4, which covers 65,536 points at that depth.
    for bits in range(1, (5 if dim < 4 else 4) + 1):
        for point in itertools.product(range(1 << bits), repeat=dim):
            check(point, bits)
            count += 1
    rng = random.Random(VERIFICATION_SEED + dim)
    for bits in range(1, 64 // dim + 1):
        for point in itertools.product((0, (1 << bits) - 1), repeat=dim):
            check(point, bits)
            count += 1
        for axis in range(dim):
            for bit in range(bits):
                point = [0] * dim
                point[axis] = 1 << bit
                check(point, bits)
                count += 1
        for _ in range(RANDOM_SAMPLES_PER_DEPTH):
            point = [rng.getrandbits(bits) for _ in range(dim)]
            check(point, bits)
            count += 1
    return count


def pack_2d(table):
    """Pack the four state bytes for each three-plane digit into one word."""
    require(len(table) == 256 and max(table) < 256, '2D packing requires 256 byte entries')
    packed = [
        sum(table[state * 64 + digit] << (state * 8) for state in range(4)) for digit in range(64)
    ]
    for state in range(4):
        for digit in range(64):
            require((packed[digit] >> (state * 8)) & 255 == table[state * 64 + digit],
                    '2D state byte packing')
    return packed


def append_table(lines, name, table, ctype):
    lines.append(f'AXOM_HILBERT_TABLE {ctype} {name}[{len(table)}] = {{')
    line = ' '
    for value in table:
        item = f' {value},'
        if len(line) + len(item) > 100:
            lines.append(line)
            line = ' '
        line += item
    lines += [line, '};']


def render():
    lines = [
        '// Copyright (c) Lawrence Livermore National Security, LLC and other',
        '// Axom Project Contributors. See top-level LICENSE and COPYRIGHT',
        '// files for dates and other details.', '//', '// SPDX-License-Identifier: (BSD-3-Clause)',
        '', '#pragma once', '', '/*!', ' * \\file HilbertTables.hpp',
        ' * \\brief Generated Hilbert transition tables.', ' */', '',
        '// Generated by internal/generate_hilbert_tables.py. Do not edit.',
        '// Verify: python src/axom/spin/internal/generate_hilbert_tables.py --check --verify',
        '//', "// Lookup tables for Hilbert encoding and decoding using Skilling's algorithm.",
        '// The generator follows the algorithm one bit per coordinate at a time, recording',
        '// the output bits and the axis orientation needed for the next step. It combines',
        '// consecutive steps into a single lookup. To reduce table size, it merges',
        '// equivalent states in 2D and handles axis reflections separately in 3D/4D.', '//',
        '// table* converts groups of Morton index bits to Hilbert index bits;',
        '// inverseTable* converts them back. Each entry also records how to update the',
        '// axis order and reflections before converting the next group.', '//',
        '// Names identify the dimension and bits per coordinate handled by a lookup.',
        '// For example, table3D2 processes two bits from each of three coordinates.',
        '// Conversion starts at state zero and processes the most significant bits first.',
        '// Use a smaller table for a short leading group; padding with zeros changes the result.',
        '//', '// See the generator for table construction, entry packing, and verification.', '',
        '#include "axom/core/Macros.hpp"', '', '#include <cstdint>', '',
        '// Dynamic device lookups need device storage; the host copy stays constexpr.',
        '#if defined(AXOM_DEVICE_CODE)', '  #define AXOM_HILBERT_TABLE inline __device__ const',
        '#else', '  #define AXOM_HILBERT_TABLE inline constexpr', '#endif', '',
        'namespace axom::spin::internal::hilbert', '{', '// clang-format off'
    ]
    for dim, groups in PRODUCTION_LEVELS.items():
        rows = machine(dim) if dim == 2 else None
        for levels in groups:
            forward = grouped(rows, dim, levels) if dim == 2 else factored(dim, levels)
            for prefix, table in [('table', forward),
                                  ('inverseTable', inverse(forward, dim, levels))]:
                ctype = 'std::uint8_t' if max(table) < 256 else 'std::uint16_t'
                require(max(table) < 65536, f'{dim}D/{levels} planes exceed 16-bit entries')
                if dim == 2 and levels == 3:
                    table = pack_2d(table)
                    ctype = 'std::uint32_t'
                    lines.append(
                        '// One byte per state, packed into a 32-bit word; select with state * 8.')
                append_table(lines, f'{prefix}{dim}D{levels}', table, ctype)
    lines += [
        '// clang-format on', '}  // namespace axom::spin::internal::hilbert', '',
        '#undef AXOM_HILBERT_TABLE', ''
    ]
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--check',
        action='store_true',
        help='compare the header with generated content, allowing CRLF line endings')
    parser.add_argument('--verify',
                        action='store_true',
                        help='check transitions, inverses, and scalar equivalence')
    parser.add_argument('--output',
                        type=Path,
                        default=Path(__file__).with_name('HilbertTables.hpp'))
    args = parser.parse_args()
    try:
        generated = render().encode('utf-8')
        if args.check:
            require(args.output.read_bytes().replace(b'\r\n', b'\n') == generated,
                    'Generated header differs; rerun generator')
            print('Generated tables match the checked-in header (allowing CRLF).', flush=True)
        if args.verify:
            for dim in PRODUCTION_LEVELS:
                count = verify(dim)
                print(f'{dim}D: {count} exhaustive/boundary/random scalar comparisons passed.',
                      flush=True)
                if dim >= 3:
                    count = verify_factored(dim)
                    print(f'{dim}D: {count} factored transitions and inverses verified.',
                          flush=True)
        if not args.check:
            args.output.write_bytes(generated)
            print(f'Wrote {args.output}', flush=True)
    except (OSError, ValueError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
