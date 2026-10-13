"""Read atomic PDB/mmCIF without assigning pseudo-atoms atomic meanings."""
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def atom_name(name):
    name = name.strip().replace('*', "'")
    return {"O1P": "OP1", "O2P": "OP2"}.get(name, name)


@dataclass
class Nucleotide:
    base: str
    chain: str
    number: int
    insertion: str = " "
    sites: dict = field(default_factory=dict)
    state: str = ""
    orientation: object = None
    neighbors: object = None


def read_structure(path):
    path = Path(path)
    records = []
    if path.suffix.lower() in ('.cif', '.mmcif'):
        d = MMCIF2Dict(str(path))
        if '_atom_site.Cartn_x' not in d:
            raise ValueError('mmCIF has no atom_site coordinates')
        def column(auth, label=None, default=''):
            return d.get('_atom_site.' + auth, d.get('_atom_site.' + str(label), [default] * len(d['_atom_site.Cartn_x'])))
        cols = [column('auth_comp_id', 'label_comp_id'), column('auth_asym_id', 'label_asym_id'),
                column('auth_seq_id', 'label_seq_id'), column('pdbx_PDB_ins_code', default=' '),
                column('auth_atom_id', 'label_atom_id'), column('Cartn_x'), column('Cartn_y'),
                column('Cartn_z'), column('label_alt_id', default='.'), column('pdbx_PDB_model_num', default='1')]
        if not cols[-1] or len({len(c) for c in cols}) != 1:
            raise ValueError('mmCIF atom_site columns are empty or have different lengths')
        first_model = cols[-1][0]
        for base, chain, num, ins, name, x, y, z, alt, model in zip(*cols):
            if model == first_model and alt in ('.', '?', '', 'A'):
                records.append((base, chain, int(num), ' ' if ins in ('.', '?') else ins, name, [float(x), float(y), float(z)]))
    elif path.suffix.lower() == '.pdb':
        for line in path.read_text().splitlines():
            if line.startswith('ENDMDL'):
                break
            if line.startswith(('ATOM  ', 'HETATM')) and line[16:17] in (' ', 'A'):
                records.append((line[17:20].strip(), line[21:22], int(line[22:26]), line[26:27],
                                line[12:16], [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    else:
        raise ValueError('This mode requires PDB or mmCIF input')
    residues = {}
    for base, chain, num, ins, name, xyz in records:
        base = {'RA': 'A', 'RU': 'U', 'RG': 'G', 'RC': 'C'}.get(base, base)
        if base not in 'AUGC' or len(base) != 1:
            raise ValueError(f'Unsupported residue {base!r} at {chain}:{num}')
        key = (chain, num, ins)
        res = residues.setdefault(key, Nucleotide(base, chain, num, ins))
        if res.base != base:
            raise ValueError(f'Conflicting identities at {key}')
        name = atom_name(name)
        if name in res.sites:
            raise ValueError(f'Duplicate atom {name} at {key}')
        xyz = np.asarray(xyz, dtype=float)
        if not np.isfinite(xyz).all():
            raise ValueError(f'Non-finite coordinates at {key}')
        res.sites[name] = xyz
    if not residues:
        raise ValueError('No RNA residues found')
    return list(residues.values())


def read_fasta(path):
    entries = []
    current = ''
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if current:
                entries.append(current)
            current = ''
        else:
            current += ''.join(line.split()).upper()
    if current:
        entries.append(current)
    if not entries or any(b not in 'AUGC' for seq in entries for b in seq):
        raise ValueError('FASTA must contain standard RNA A/U/G/C sequences')
    return entries
