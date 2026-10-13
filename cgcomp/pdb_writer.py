from pathlib import Path
import numpy as np
from .fragments import BACKBONE


def pdb_text(residues, backbone=False):
    chains = list(dict.fromkeys(r.chain for r in residues))
    if any(len(c) != 1 for c in chains):
        raise ValueError('PDB requires one-character chain identifiers; mmCIF chains cannot be truncated safely')
    lines, serial, previous_chain = [], 1, None
    for r in residues:
        if not -999 <= r.number <= 9999 or len(r.insertion) != 1:
            raise ValueError('Residue identifier exceeds PDB field limits')
        if previous_chain is not None and previous_chain != r.chain:
            lines.append('TER')
        previous_chain = r.chain
        for name, xyz in r.sites.items():
            if backbone and name not in BACKBONE:
                continue
            if len(name) > 4 or serial > 99999 or not np.isfinite(xyz).all():
                raise ValueError('Atom cannot be represented in PDB')
            coords = [f'{x:8.3f}' for x in xyz]
            if any(len(c) != 8 for c in coords):
                raise ValueError('Coordinate exceeds PDB field limits')
            atom = f' {name:<3}' if len(name) < 4 else name
            lines.append(f'ATOM  {serial:5d} {atom} {r.base:>3} {r.chain}{r.number:4d}{r.insertion}   '
                         + ''.join(coords) + f'  1.00  0.00          {name[0]:>2}  ')
            serial += 1
    lines.extend(['TER', 'END'])
    return '\n'.join(lines) + '\n'


def write_outputs(residues, directory):
    # Validate both strings before creating files.
    texts = {'backbone.pdb': pdb_text(residues, True), 'all_atom.pdb': pdb_text(residues)}
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, content in texts.items():
        (directory / name).write_text(content)
