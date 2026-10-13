"""Archive bead definitions, shared by PDB and mmCIF readers.

COM means mass-weighted heavy-atom center. Geometric averages are explicitly
separate. iFoldRNA's supplied export groups differ from the original native
three-bead definition; the audited archive profile is the RNAgrainy-compatible
default. Native publication profiles remain distinct and are not implemented.
"""
from functools import partial
import numpy as np
from .generic import read_structure
from .hire_rna import virtual_sites as hire_sites

MASSES = {'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.974}
SUGAR_RING = ("C1'", "C2'", "C3'", "C4'", "O4'")
PURINE_RING = ('N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4')
PYRIMIDINE_RING = ('N1', 'C2', 'N3', 'C4', 'C5', 'C6')
COMMON_R = ('N9', 'C8', 'N7', 'C5', 'C4', 'N3')
COMMON_Y = ('N1', 'C5', 'C6')
GROUPS = {
    'isrna1': {
        'A': {'R_C': COMMON_R, 'A_C': ('C6', 'N6'), 'A_N': ('C2', 'N1')},
        'G': {'R_C': COMMON_R, 'G_O': ('C6', 'O6'), 'G_N': ('N1', 'C2', 'N2')},
        'C': {'Y_C': COMMON_Y, 'C_N': ('C2', 'O2', 'N3', 'C4', 'N4')},
        'U': {'Y_C': COMMON_Y, 'U_O': ('C2', 'O2', 'N3', 'C4', 'O4')},
    },
    'isrna2': {
        'A': {'R_1': COMMON_R, 'A_1': ('C6', 'N6'), 'A_2': ('C2', 'N1')},
        'G': {'R_1': COMMON_R, 'G_1': ('C6', 'O6', 'N1'), 'G_2': ('C2', 'N2')},
        'C': {'Y_1': COMMON_Y, 'Y_2': ('C2', 'O2'), 'C_1': ('N3', 'C4', 'N4')},
        'U': {'Y_1': COMMON_Y, 'Y_2': ('C2', 'O2'), 'U_1': ('N3', 'C4', 'O4')},
    },
}


def center(atoms, names, mass=False):
    weights = [MASSES[n[0]] for n in names] if mass else None
    return np.average([atoms[n] for n in names], axis=0, weights=weights)


def atomic_names(mode, base):
    if mode in ('vfold', 'febrna'):
        return ('P', "C4'", 'N9' if base in 'AG' else 'N1')
    if mode in ('simrna', 'rna-jp'):
        return ('P', "C4'", 'N9', 'C6', 'C2') if base in 'AG' else ('P', "C4'", 'N1', 'C2', 'C4')
    return ()


def reference_sites(mode, base, atoms):
    if mode == 'generic':
        return atoms
    if mode == 'hire-rna':
        return hire_sites(base, atoms)
    if mode in ('isrna1', 'isrna2'):
        return {'P': atoms['P'], 'S': atoms["C4'"],
                **{n: center(atoms, group, True) for n, group in GROUPS[mode][base].items()}}
    if mode == 'toprna':
        return {'P': atoms['P'], 'S': center(atoms, (*SUGAR_RING, "C5'")),
                'B': center(atoms, PURINE_RING if base in 'AG' else
                            (*PYRIMIDINE_RING, 'N4') if base == 'C' else PYRIMIDINE_RING)}
    if mode == 'ifoldrna':
        return {'P': center(atoms, ('P', 'OP1', 'OP2'), True),
                'S': center(atoms, (*SUGAR_RING, "C5'", "O3'", "O2'"), True),
                'B': center(atoms, PURINE_RING if base in 'AG' else (*PYRIMIDINE_RING, 'O2'), True)}
    return {n: atoms[n] for n in atomic_names(mode, base)}


def real_atom_sites(mode, residue):
    """Atom name -> original coordinate; bead name need not equal atom name."""
    if mode in ('isrna1', 'isrna2'):
        return {'P': residue.sites['P'], "C4'": residue.sites['S']}
    if mode == 'toprna':
        return {'P': residue.sites['P']}
    if mode == 'ifoldrna':
        return {}  # Even the P-labeled bead is NOT a phosphorus atom here.
    if mode == 'hire-rna':
        return {n: p for n, p in residue.sites.items() if n not in ('B1', 'B2')}
    return {n: p for n, p in residue.sites.items() if not n.lstrip('0123456789').startswith('H')}


def parse(path, mode):
    residues = read_structure(path)
    for r in residues:
        if mode in ('isrna1', 'isrna2'):
            expected = {'P', 'S', *GROUPS[mode][r.base]}
        elif mode in ('ifoldrna', 'toprna'):
            expected = {'P', 'S', 'B'}
        else:
            expected = set(atomic_names(mode, r.base))
        if set(r.sites) != expected:
            raise ValueError(f'{mode} {r.chain}:{r.number} requires sites {sorted(expected)}')
    return residues


PARSERS = {mode: partial(parse, mode=mode) for mode in
           ('vfold', 'febrna', 'rna-jp', 'toprna', 'ifoldrna', 'isrna1', 'isrna2')}
