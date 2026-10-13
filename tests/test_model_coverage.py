import json
import shlex
from pathlib import Path
import numpy as np
import pytest
from cgcomp.models import MODELS
from cgcomp.parsers import mapped, oxrna
from cgcomp.parsers.generic import read_structure
from cgcomp.reconstruction import kabsch, parse_input
from conftest import ROOT


def test_complete_archive_model_catalog():
    # Inputs are top-level files; *_full.pdb directories hold generated outputs.
    stems = {path.stem for path in (ROOT / 'CG_models').iterdir()
             if path.is_file() and path.suffix in ('.pdb', '.cif', '.dat')}
    assert stems == {model.stem for model in MODELS.values()}
    inventory = json.loads((ROOT / 'docs/model_inventory.json').read_text())
    assert len(inventory['models']) == len(MODELS) == 15
    assert sum(not m.blocker for m in MODELS.values()) == 12  # 11 CG modes + generic.


def test_inventory_cli_readme_examples_outputs_and_coverage_agree():
    inventory = json.loads((ROOT / 'docs/model_inventory.json').read_text())['models']
    commands = json.loads((ROOT / 'docs/example_commands.json').read_text())
    report = json.loads((ROOT / 'docs/validation_report.json').read_text())
    readme = (ROOT / 'README.md').read_text()
    assert {entry['mode'] for entry in inventory.values()} == set(MODELS)
    assert {entry['mode'] for entry in report['modes']} == set(MODELS)
    assert {entry['mode'] for entry in commands} == {m for m, model in MODELS.items() if not model.blocker}
    from Bio.PDB import PDBParser
    for mode, model in MODELS.items():
        assert f'`{mode}`' in readme
        if model.blocker:
            row = next(r for r in report['modes'] if r['mode'] == mode)
            assert row['outputs'] == [] and row['blocker'] == model.blocker
            continue
        examples = [entry for entry in commands if entry['mode'] == mode]
        documented = [shlex.split(line) for line in readme.splitlines() if line.startswith('cgcomp -i ')]
        def options(command):
            result = dict(zip(command[1::2], command[2::2]))
            if mode in ('ifoldrna', 'toprna'):
                result.pop('--bead-profile', None)  # Historical manifest uses the now-default archive mapping.
            return result
        assert any(options(entry['command']) == options(command) for entry in examples for command in documented)
        for entry in examples:
            for name, atoms in [('backbone.pdb', 1224), ('all_atom.pdb', 2181)]:
                structure = PDBParser(QUIET=True).get_structure('example', ROOT / entry['output'] / name)
                assert len(list(structure.get_residues())) == 102
                assert len(list(structure.get_atoms())) == atoms


@pytest.mark.parametrize('mode', ['simrna', 'rna-jp', 'vfold', 'febrna', 'hire-rna', 'toprna', 'ifoldrna', 'isrna1', 'isrna2'])
def test_mapping_against_every_supplied_atomic_residue(mode):
    reference = read_structure(ROOT / 'CG_models/exp.pdb')
    cg = parse_input(mode, ROOT / 'CG_models' / (MODELS[mode].stem + '.pdb'))
    for atomic, bead in zip(reference, cg):
        assert atomic.base == bead.base
        sites = mapped.reference_sites(mode, atomic.base, atomic.sites)
        assert sites.keys() == bead.sites.keys()
        for name, xyz in sites.items():
            assert np.linalg.norm(xyz - bead.sites[name]) < .0015


def test_nares_midpoints_and_orientation_ambiguity(fragments):
    residues = parse_input('nares-2p', ROOT / 'CG_models/exp_Nares-2P.pdb')
    for previous, current in zip(residues, residues[1:]):
        assert np.linalg.norm(current.sites['P'] - (previous.sites['S'] + current.sites['S']) / 2) < .001
    # S-B rotation preserves both local sites, hence every inter-residue midpoint
    # as well. It changes real atom positions, proving the inverse is not unique.
    atoms = fragments['A']
    sites = mapped.reference_sites('nares-2p', 'A', atoms)
    axis = sites['B'] - sites['S']
    axis /= np.linalg.norm(axis)
    rotation = 2 * np.outer(axis, axis) - np.eye(3)  # Proper half-turn about S-B.
    rotated = {n: (xyz - sites['S']) @ rotation + sites['S'] for n, xyz in atoms.items()}
    new_sites = mapped.reference_sites('nares-2p', 'A', rotated)
    assert all(np.allclose(new_sites[n], xyz, atol=1e-12) for n, xyz in sites.items())
    assert np.linalg.norm(rotated['P'] - atoms['P']) > 1


def test_ifold_archive_profile_is_not_native_ding2008():
    atomic = read_structure(ROOT / 'CG_models/exp.pdb')
    cg = parse_input('ifoldrna', ROOT / 'CG_models/exp_iFoldRNA.pdb')
    sugar_differences = [np.linalg.norm(mapped.center(a.sites, mapped.SUGAR_RING, True) - c.sites['S']) for a, c in zip(atomic, cg)]
    purine_differences = [np.linalg.norm(mapped.center(a.sites, mapped.PYRIMIDINE_RING) - c.sites['B'])
                          for a, c in zip(atomic, cg) if a.base in 'AG']
    assert max(sugar_differences) > .5
    assert max(purine_differences) > .5


@pytest.mark.parametrize('base', list('AUGC'))
def test_oxrna_reference_features_against_official_source_oracle(base, fragments):
    oracle = json.loads((ROOT / 'docs/oxrna_reference_frames.json').read_text())
    origin, proxies = oxrna.fragment_proxies(base, fragments[base])
    assert np.allclose(origin, oracle['bases'][base]['origin'], atol=1e-12)
    assert np.allclose(proxies, oracle['bases'][base]['proxies'], atol=1e-12)
    a1, a2, a3 = np.eye(3)
    from cgcomp.parsers.generic import Nucleotide
    target = oxrna.target_proxies(Nucleotide(base, '1', 1, orientation=np.eye(3)))
    rotation, translation, residual = kabsch(proxies, target, center=False)
    assert np.allclose(rotation, oracle['bases'][base]['canonical_rotation'], atol=1e-12)
    assert residual == pytest.approx(oracle['bases'][base]['proxy_rmsd'], abs=1e-12)
    assert np.array_equal(translation, np.zeros(3))
