from importlib.resources import files
from pathlib import Path
from .parsers.generic import read_structure

BACKBONE = ('P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'")


def load_fragments(directory=None):
    directory = Path(directory) if directory else Path(str(files('cgcomp').joinpath('data/fragments')))
    result = {}
    for base in 'AUGC':
        residues = read_structure(directory / (base + '.pdb'))
        if len(residues) != 1 or residues[0].base != base:
            raise ValueError(f'{base}.pdb must contain exactly one {base} nucleotide')
        atoms = {n: c.copy() for n, c in residues[0].sites.items() if not n.lstrip('0123456789').startswith('H')}
        if not set(BACKBONE).issubset(atoms):
            raise ValueError(f'{base}.pdb is missing backbone heavy atoms')
        required = {'N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'} if base in 'AG' else {'N1', 'C2', 'N3', 'C4', 'C5', 'C6'}
        required |= {'A': {'N6'}, 'G': {'O6', 'N2'}, 'C': {'O2', 'N4'}, 'U': {'O2', 'O4'}}[base]
        if not required.issubset(atoms):
            raise ValueError(f'{base}.pdb is missing base heavy atoms')
        result[base] = atoms
    return result
