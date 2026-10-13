from .generic import read_structure


def parse(path):
    residues = read_structure(path)
    for r in residues:
        expected = {'P', "C4'", 'N9', 'C6', 'C2'} if r.base in 'AG' else {'P', "C4'", 'N1', 'C2', 'C4'}
        if set(r.sites) != expected:
            raise ValueError(f'SimRNA {r.chain}:{r.number} requires {sorted(expected)}; use generic for partial atomic input')
    return residues
