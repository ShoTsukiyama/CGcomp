"""Native oxRNA reader and RNA-specific atomic backmapping constraints.

Positions remain in native simulation units until reconstruction. The atomistic
proxy convention follows the official analysis/oxDNA_PDB RNA branch, not DNA.
"""
import itertools
from pathlib import Path
import numpy as np
from .generic import Nucleotide

ANGSTROM_PER_UNIT = 8.518
COM_SHIFT = 0.5  # Å; official backmapper's documented empirical origin shift.


def fragment_proxies(base, atoms):
    """Official RNA template COM, pairing-edge axis, and signed ring normal."""
    origin = np.mean(list(atoms.values()), axis=0)
    edge = ('C6', 'N1', 'C2') if base in 'AG' else ('C4', 'N3', 'C2')
    a1 = np.mean([atoms[n] for n in edge], axis=0) - origin
    a1 /= np.linalg.norm(a1)
    base_atoms = [xyz for n, xyz in atoms.items() if "'" not in n and 'P' not in n]
    positive = atoms["O4'"] - np.mean(base_atoms, axis=0)
    normal = np.zeros(3)
    for p, q, r in itertools.permutations(('C2', 'C4', 'C5', 'C6', 'N1', 'N3'), 3):
        direction = np.cross(atoms[p] - atoms[q], atoms[p] - atoms[r])
        length = np.linalg.norm(direction)
        if length < 1e-10:
            continue
        direction /= length
        normal += direction if np.dot(direction, positive) >= 0 else -direction
    if np.linalg.norm(normal) < 1e-10:
        raise ValueError('Cannot determine oxRNA fragment base normal')
    a3 = normal / np.linalg.norm(normal)
    a2 = np.cross(a3, a1)
    # The official RNA proxies are C5' and the pairing edge N1(A/G)/N3(C/U).
    proxy = np.array([a1, a3, a2, atoms["C5'"] - origin,
                      atoms['N1' if base in 'AG' else 'N3'] - origin])
    return origin, proxy


def target_proxies(residue):
    a1, a2, a3 = residue.orientation
    return np.array([a1, a3, a2,
                     ANGSTROM_PER_UNIT * (-0.4 * a1 + 0.2 * a3),
                     ANGSTROM_PER_UNIT * (0.4 * a1)])


def read_topology(path):
    lines = [l.strip() for l in Path(path).read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')]
    header = lines[0].split()
    count, strands = map(int, header[:2])
    result = []
    if len(header) == 3 and header[2] == '5->3':
        if len(lines) != strands + 1:
            raise ValueError('oxRNA topology strand count mismatch')
        for strand, line in enumerate(lines[1:], 1):
            fields = line.split()
            sequence = fields[0]
            settings = dict(f.split('=', 1) for f in fields[1:] if '=' in f)
            if settings.get('type', 'RNA') != 'RNA':
                raise ValueError('Topology explicitly selects a non-RNA model')
            start = len(result)
            circular = settings.get('circular', 'false').lower() in ('true', '1')
            for i, base in enumerate(sequence):
                n5 = start + i - 1 if i else (start + len(sequence) - 1 if circular else -1)
                n3 = start + i + 1 if i < len(sequence) - 1 else (start if circular else -1)
                result.append((strand, base, n3, n5))
    elif len(header) == 2:
        for line in lines[1:]:
            fields = line.split()
            if len(fields) != 4:
                raise ValueError('Expected strand base n3 n5 in legacy topology')
            strand, base, n3, n5 = fields
            result.append((int(strand), base, int(n3), int(n5)))
    else:
        raise ValueError('Unsupported oxRNA topology header')
    if len(result) != count or len({r[0] for r in result}) != strands:
        raise ValueError('oxRNA topology nucleotide/strand count mismatch')
    for i, (strand, base, n3, n5) in enumerate(result):
        if base not in ('A', 'U', 'G', 'C'):
            raise ValueError(f'Unsupported oxRNA base {base}')
        for neighbor, reciprocal in ((n3, 3), (n5, 2)):
            if neighbor != -1 and (not 0 <= neighbor < count or neighbor == i or result[neighbor][0] != strand or result[neighbor][reciprocal] != i):
                raise ValueError('Invalid or nonreciprocal oxRNA strand connectivity')
    return result


def parse(path, topology):
    if topology is None:
        raise ValueError('oxRNA requires its matching native topology; DAT contains no sequence')
    records = read_topology(topology)
    lines = Path(path).read_text().splitlines()
    if len(lines) != len(records) + 3 or not all(lines[i].strip().startswith(key) for i, key in enumerate(('t =', 'b =', 'E ='))):
        raise ValueError('Expected one oxRNA configuration with t/b/E headers and one row per nucleotide')
    if not np.isfinite([float(v) for l in lines[:3] for v in l.split('=', 1)[1].split()]).all():
        raise ValueError('Non-finite oxRNA header')
    residues, counters = [], {}
    for record, line in zip(records, lines[3:]):
        values = np.asarray([float(x) for x in line.split()])
        if len(values) != 15 or not np.isfinite(values).all():
            raise ValueError('oxRNA rows must have 15 finite position/orientation/momentum values')
        a1, a3 = values[3:6], values[6:9]
        if abs(np.linalg.norm(a1) - 1) > 1e-4 or abs(np.linalg.norm(a3) - 1) > 1e-4 or abs(np.dot(a1, a3)) > 1e-4:
            raise ValueError('oxRNA a1/a3 must be orthonormal')
        strand, base, n3, n5 = record
        counters[strand] = counters.get(strand, 0) + 1
        residues.append(Nucleotide(base, str(strand), counters[strand], sites={'COM': values[:3].copy()},
                                   orientation=np.array([a1, np.cross(a3, a1), a3]), neighbors=(n3, n5)))
    return residues
