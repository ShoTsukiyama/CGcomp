"""SPQR native whitespace PDB, including its nonstandard coordinate fields."""
from pathlib import Path
import numpy as np
from .generic import Nucleotide


def frame_sites(base, atoms):
    center = np.mean([atoms[n] for n in ('C2', 'C4', 'C6')], axis=0)
    x = atoms['C2'] - center
    x /= np.linalg.norm(x)
    y = atoms['C6' if base in 'AG' else 'C4'] - center
    y -= np.dot(y, x) * x
    y /= np.linalg.norm(y)
    return np.array([center, center + x, center + y, center + np.cross(x, y)])


def reference_sites(base, atoms):
    """Native simulation A3 forward mapping (MC_map_sugar / mc_integrate.h).

    SUGR is an ideal anti-state position in the base frame, even when sugar
    atoms are available. PHOS is P when present (the converter synthesizes a
    state-dependent default only when P is absent).
    """
    frame = frame_sites(base, atoms)
    center = frame[0]
    axes = frame[1:] - center
    # pdb2spqr.py uses -4.8467 for purines, but the native simulation
    # recomputes SUGR with -4.48467. These are distinct source conventions.
    sugar = (-1.56584, -4.48467, .709821) if base in 'AG' else (1.05136, -3.42329, .652454)
    return dict(BASE=center, XVEC=frame[1], YVEC=frame[2],
                SUGR=center + np.asarray(sugar) @ axes, PHOS=atoms['P'])


def target_frame(r):
    center = r.sites['BASE']
    x, y = r.sites['XVEC'] - center, r.sites['YVEC'] - center
    if abs(np.linalg.norm(x) - 1) > .02 or abs(np.linalg.norm(y) - 1) > .02 or abs(np.dot(x, y)) > .02:
        raise ValueError(f'Invalid SPQR orientation at {r.chain}:{r.number}')
    x = x / np.linalg.norm(x)
    y = y - np.dot(x, y) * x
    y /= np.linalg.norm(y)
    return np.array([center, center + x, center + y, center + np.cross(x, y)])


def parse(path):
    if Path(path).suffix.lower() != '.pdb':
        raise ValueError('SPQR requires its native PDB representation')
    residues = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith('ATOM'):
            continue
        fields = line.split()
        if len(fields) < 9:
            raise ValueError('Malformed native SPQR record')
        _, serial, name, base, chain, number, x, y, z, *state = fields
        if base not in ('A', 'G', 'C', 'U'):
            raise ValueError(f'Unsupported SPQR base {base}')
        key = (chain, int(number))
        r = residues.setdefault(key, Nucleotide(base, chain, int(number)))
        if r.base != base or name in r.sites:
            raise ValueError(f'Conflicting SPQR residue {key}')
        r.sites[name] = np.array([float(x), float(y), float(z)])
        if name == 'BASE':
            r.state = state[0][:2] if state else ''
    if not residues:
        raise ValueError('No SPQR nucleotides found')
    for r in residues.values():
        if set(r.sites) != {'BASE', 'XVEC', 'YVEC', 'SUGR', 'PHOS'}:
            raise ValueError(f'Incomplete SPQR residue {r.chain}:{r.number}')
        if not all(np.isfinite(v).all() for v in r.sites.values()):
            raise ValueError('Non-finite SPQR coordinates')
        if r.state != 'A3':
            raise ValueError(f'SPQR state {r.state!r} unsupported by the supplied anti/C3-endo fragments')
        target_frame(r)
    return list(residues.values())
