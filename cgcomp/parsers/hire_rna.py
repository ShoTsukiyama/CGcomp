"""HiRE: five real atoms and mass-weighted aromatic-ring centers."""
import numpy as np
from .generic import read_structure

BACKBONE = ('P', "O5'", "C5'", "C4'", "C1'")
SIX_RING = ('N1', 'C2', 'N3', 'C4', 'C5', 'C6')
FIVE_RING = ('N9', 'C8', 'N7', 'C5', 'C4')


def ring_center(atoms, names):
    weights = [14.007 if n.startswith('N') else 12.011 for n in names]
    return np.average([atoms[n] for n in names], axis=0, weights=weights)


def virtual_sites(base, atoms):
    sites = {n: atoms[n] for n in BACKBONE}
    sites['B1'] = ring_center(atoms, FIVE_RING if base in 'AG' else SIX_RING)
    if base in 'AG':
        sites['B2'] = ring_center(atoms, SIX_RING)
    return sites


def parse(path):
    residues = read_structure(path)
    for r in residues:
        expected = set(BACKBONE) | ({'B1', 'B2'} if r.base in 'AG' else {'B1'})
        if set(r.sites) != expected:
            raise ValueError(f'HiRE-RNA {r.chain}:{r.number} requires {sorted(expected)}')
    return residues
