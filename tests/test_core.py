import zipfile
import numpy as np
import pytest
from Bio.PDB import PDBParser
from cgcomp.fragments import load_fragments, BACKBONE
from cgcomp.parsers.generic import Nucleotide, read_structure, atom_name
from cgcomp.parsers.hire_rna import virtual_sites
from cgcomp.reconstruction import kabsch, reconstruct
from cgcomp.pdb_writer import pdb_text
from conftest import ROOT


def test_original_data_and_bundled_fragments_unchanged():
    for archive in ('CG_models.zip', 'single_nucleotide_fragments.zip'):
        with zipfile.ZipFile(ROOT / archive) as z:
            for item in z.infolist():
                if item.is_dir() or item.filename.startswith('__MACOSX/'):
                    continue
                assert (ROOT / item.filename).read_bytes() == z.read(item)
    for base in 'AUGC':
        assert (ROOT / 'single_nucleotide_fragments' / (base + '.pdb')).read_bytes() == (ROOT / 'cgcomp/data/fragments' / (base + '.pdb')).read_bytes()


def test_fragment_loading(fragments):
    assert {b: len(a) for b, a in fragments.items()} == {'A': 22, 'U': 20, 'G': 23, 'C': 20}
    bundled = load_fragments()
    for base in fragments:
        assert set(BACKBONE) <= fragments[base].keys()
        for atom, xyz in fragments[base].items():
            assert np.array_equal(xyz, bundled[base][atom])


@pytest.mark.parametrize('base', list('AUGC'))
def test_kabsch_known_transform(base, fragments):
    source = np.array(list(fragments[base].values()))
    rng = np.random.default_rng(412)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    translation = np.array([4., -12., 8.])
    target = source @ rotation + translation
    fitted, shift, rmsd = kabsch(source, target)
    assert np.allclose(fitted, rotation, atol=1e-12)
    assert np.allclose(shift, translation, atol=1e-12)
    assert rmsd < 1e-12
    assert np.linalg.det(fitted) == pytest.approx(1)


def test_kabsch_rejects_underdetermined():
    with pytest.raises(ValueError, match='three'):
        kabsch([[0, 0, 0]] * 2, [[0, 0, 0]] * 2)
    with pytest.raises(ValueError, match='collinear'):
        kabsch([[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]])


def test_kabsch_no_reflection(fragments):
    source = np.array(list(fragments['A'].values()))
    target = source * [-1, 1, 1]
    rotation, _, rmsd = kabsch(source, target)
    assert np.linalg.det(rotation) == pytest.approx(1)
    assert rmsd > .1


@pytest.mark.parametrize('stem', ['exp', 'exp_SimRNA', 'exp_HiRE-RNA', 'exp_VFold', 'exp_YUP', 'exp_Nares-2P', 'exp_NAST', 'exp_iFoldRNA', 'exp_FebRNA', 'exp_TopRNA', 'exp_isRNA2', 'exp_RNA-JP', 'exp_isRNA1'])
def test_actual_pdb_and_mmcif(stem):
    pdb = read_structure(ROOT / 'CG_models' / (stem + '.pdb'))
    assert len(pdb) == 102
    if stem == 'exp':
        return  # The archive contains no exp.cif.
    cif = read_structure(ROOT / 'CG_models' / (stem + '.cif'))
    for a, b in zip(pdb, cif):
        assert (a.base, a.chain, a.number, a.insertion) == (b.base, b.chain, b.number, b.insertion)
        assert a.sites.keys() == b.sites.keys()
        for name in a.sites:
            assert np.allclose(a.sites[name], b.sites[name], atol=.001)


def test_atom_name_aliases():
    assert atom_name(' C4* ') == "C4'"
    assert atom_name('O1P') == 'OP1'
    assert atom_name('O2P') == 'OP2'
    assert atom_name('B1') == 'B1'


@pytest.mark.parametrize('base', list('AUGC'))
def test_hire_ring_mapping_against_supplied_reference(base):
    atomic = read_structure(ROOT / 'CG_models/exp.pdb')
    cg = read_structure(ROOT / 'CG_models/exp_HiRE-RNA.pdb')
    for a, c in zip(atomic, cg):
        if a.base == base:
            sites = virtual_sites(base, a.sites)
            for name in sites:
                assert np.linalg.norm(sites[name] - c.sites[name]) < .0015


@pytest.mark.parametrize('base', list('AUGC'))
@pytest.mark.parametrize('group', ['backbone1', 'backbone2', 'base', 'all'])
def test_generic_recovers_supplied_fragment(base, group, fragments, tmp_path):
    atoms = fragments[base]
    names = {'backbone1': ['P', "C4'", "C3'"], 'backbone2': ['P', "O5'", "C4'"],
             'base': ['N9', 'C4', 'C8'] if base in 'AG' else ['N1', 'C2', 'C6'], 'all': list(atoms)}[group]
    angle = .731
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    original = {n: xyz @ rotation + [2, 5, -7] for n, xyz in atoms.items()}
    partial = Nucleotide(base, 'B', 7, 'A', {n: original[n] for n in names})
    path = tmp_path / 'partial.pdb'
    path.write_text(pdb_text([partial]))
    parsed = read_structure(path)
    output, info = reconstruct(parsed, fragments, 'generic')
    assert info[0]['constraints'] == names
    assert max(np.linalg.norm(output[0].sites[n] - xyz) for n, xyz in original.items()) < .008
    for name in names:
        assert np.array_equal(output[0].sites[name], parsed[0].sites[name])
    independent = PDBParser(QUIET=True).get_structure('partial', path)
    r = next(independent.get_residues())
    assert r.id == (' ', 7, 'A') and r.parent.id == 'B'


def test_writer_independent_parser(fragments, tmp_path):
    residues = [Nucleotide('A', 'A', 1, sites=fragments['A']), Nucleotide('C', 'B', 9, sites=fragments['C'])]
    for backbone in (False, True):
        path = tmp_path / 'output.pdb'
        path.write_text(pdb_text(residues, backbone))
        parsed = PDBParser(QUIET=True).get_structure('test', path)
        assert len(list(parsed.get_chains())) == 2
        for source, residue in zip(residues, parsed.get_residues()):
            expected = set(BACKBONE) if backbone else set(source.sites)
            assert {a.name for a in residue} == expected
            assert all(a.element in ('P', 'O', 'C', 'N') for a in residue)
            for atom in residue:
                assert np.allclose(atom.coord, source.sites[atom.name], atol=.00051)


def test_generic_rejects_pseudoatoms_and_insufficient_atoms(fragments):
    with pytest.raises(ValueError, match='Unmapped'):
        reconstruct(read_structure(ROOT / 'CG_models/exp_Nares-2P.pdb'), fragments, 'generic')
    with pytest.raises(ValueError, match='three'):
        reconstruct(read_structure(ROOT / 'CG_models/exp_NAST.pdb'), fragments, 'generic')
