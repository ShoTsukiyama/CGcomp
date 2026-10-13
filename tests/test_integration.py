import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from cgcomp.fragments import BACKBONE
from cgcomp.reconstruction import PARSERS, parse_input
from cgcomp.models import MODELS
from cgcomp.parsers import mapped, spqr, oxrna
from conftest import ROOT

CLI = Path(sys.executable).parent / 'cgcomp'
CASES = [(mode, model.stem, ext) for mode, model in MODELS.items() if mode != 'generic'
         for ext in model.formats if (ROOT / 'CG_models' / f'{model.stem}.{ext}').is_file()]
CASES += [('generic', stem, ext) for stem in ('exp_VFold', 'exp_FebRNA', 'exp_RNA-JP', 'exp_SimRNA') for ext in ('pdb', 'cif')]
CASES += [('generic', 'exp', 'pdb')]


def run_cli(args, cwd):
    assert CLI.is_file(), 'Install CGcomp into the test Python environment before running tests'
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    return subprocess.run([str(CLI), *map(str, args)], cwd=cwd, env=env, capture_output=True, text=True)


@pytest.mark.parametrize('mode,stem,ext', CASES)
def test_installed_cli_actual_archive_models(mode, stem, ext, tmp_path, fragments):
    input_file = ROOT / 'CG_models' / f'{stem}.{ext}'
    output = tmp_path / 'output'
    options = ['--topology', ROOT / 'examples/auxiliary/exp_oxRNA.top'] if mode == 'oxrna' else []
    completed = run_cli(['--input', input_file, '--mode', mode, '--output', output, *options], tmp_path)
    assert completed.returncode == 0, completed.stderr
    cg = parse_input(mode, input_file, ROOT / 'examples/auxiliary/exp_oxRNA.top' if mode == 'oxrna' else None)
    report = json.loads((output / 'fit_report.json').read_text())['residues']
    for filename, backbone in [('backbone.pdb', True), ('all_atom.pdb', False)]:
        # This parser is independent of CGcomp's PDB reader/writer.
        structure = PDBParser(QUIET=True).get_structure('output', output / filename)
        residues = list(structure.get_residues())
        assert len(residues) == len(cg) == 102
        assert ''.join(r.resname for r in residues) == ''.join(r.base for r in cg)
        for r, original in zip(residues, cg):
            assert r.parent.id == original.chain and r.id[1:] == (original.number, original.insertion)
            expected = set(BACKBONE) if backbone else set(fragments[original.base])
            assert {a.name for a in r} == expected
            xyz = np.array([a.coord for a in r])
            assert np.isfinite(xyz).all()
            assert np.linalg.matrix_rank(xyz - xyz.mean(0), tol=.001) == 3
            distances = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=2)
            assert distances[np.triu_indices(len(xyz), 1)].min() > .1
            if mode not in ('spqr', 'oxrna'):
                for name, xyz in mapped.real_atom_sites(mode, original).items():
                    if name in expected:
                        assert np.allclose(r[name].coord, xyz, atol=.0006, rtol=0)
            elif mode == 'spqr' and 'P' in expected:
                # Independent PHOS cannot also be imposed on a rigid fragment
                # whose position/orientation are fixed by the base frame.
                from cgcomp.reconstruction import kabsch
                rotation, translation, _ = kabsch(spqr.frame_sites(original.base, fragments[original.base]), spqr.target_frame(original))
                for name in expected:
                    assert np.allclose(r[name].coord, fragments[original.base][name] @ rotation + translation, atol=.0006, rtol=0)
    for original, residue, fit in zip(cg, residues, report):
        fragment = fragments[original.base]
        if mode == 'oxrna':
            # Oracle computed by the official RNA backmapper helpers, separately
            # from CGcomp, using these very same supplied fragment coordinates.
            oracle = json.loads((ROOT / 'docs/oxrna_reference_frames.json').read_text())['bases'][original.base]
            reference_rotation = np.asarray(oracle['canonical_rotation']) @ original.orientation
            origin = np.asarray(oracle['origin'])
            proxies = np.asarray(oracle['proxies'])
            pos = original.sites['COM'] * 8.518
            for name, xyz in fragment.items():
                expected = (xyz - origin - .5 * proxies[0]) @ reference_rotation + pos
                assert np.allclose(residue[name].coord, expected, atol=.0006, rtol=0)
            assert fit['proxy_fit_rmsd'] == pytest.approx(oracle['proxy_rmsd'], abs=1e-12)
            output_atoms = {a.name: a.coord.astype(float) for a in residue}
            output_origin, output_proxies = oxrna.fragment_proxies(original.base, output_atoms)
            assert np.allclose(output_origin + .5 * output_proxies[0], pos, atol=.0015, rtol=0)
            differences = np.linalg.norm(output_proxies[:3] - oxrna.target_proxies(original)[:3], axis=1)
            assert np.allclose(differences, fit['orientation_vector_residuals'], atol=.0015, rtol=0)
            assert fit['proxy_fit_rmsd'] < .73
        elif mode == 'spqr':
            actual_sites = spqr.frame_sites(original.base, {a.name: a.coord.astype(float) for a in residue})
            assert np.allclose(actual_sites, spqr.target_frame(original), atol=.002)
            assert fit['fit_rmsd_angstrom'] < 1e-10
            assert fit['sugar_center_deviation_angstrom'] < .18
            assert fit['phosphate_site_deviation_angstrom'] < 1.5
            sites = spqr.reference_sites(original.base, {a.name: a.coord.astype(float) for a in residue})
            for name, xyz in original.sites.items():
                # Recomputing SUGR amplifies PDB atom rounding through its
                # ~5 Å frame offset; direct sites retain the usual tolerance.
                tolerance = .006 if name == 'SUGR' else .0015
                assert np.linalg.norm(sites[name] - xyz) == pytest.approx(fit['output_site_residuals_angstrom'][name], abs=tolerance)
        else:
            mapping = mapped.reference_sites(mode, original.base, fragment)
            names = [n for n in original.sites if n in mapping]
            source = np.array([mapping[n] for n in names])
            target = np.array([original.sites[n] for n in names])
            independent = SVDSuperimposer()
            independent.set(target, source)
            independent.run()
            assert fit['fit_rmsd_angstrom'] == pytest.approx(independent.get_rms(), abs=1e-9)
            rotation, translation = independent.get_rotran()
            for name in set(fragment) - set(mapped.real_atom_sites(mode, original)):
                assert np.allclose(residue[name].coord, fragment[name] @ rotation + translation, atol=.0006)
            final_mapping = mapped.reference_sites(mode, original.base, {a.name: a.coord.astype(float) for a in residue})
            for name in names:
                residual = np.linalg.norm(final_mapping[name] - original.sites[name])
                assert residual == pytest.approx(fit['output_site_residuals_angstrom'][name], abs=.0015)
            # A fully resolved experimental nucleotide may differ substantially
            # from the template; its coordinates must still be preserved.
            if stem != 'exp':
                assert fit['fit_rmsd_angstrom'] < 1.6


def test_cli_help_and_all_mode_listing(tmp_path):
    help_result = run_cli(['--help'], tmp_path)
    assert help_result.returncode == 0
    assert '--list-modes' in help_result.stdout
    listed = run_cli(['--list-modes'], tmp_path)
    assert listed.returncode == 0
    assert {line.split()[0] for line in listed.stdout.splitlines()} == set(MODELS)
    for mode, model in MODELS.items():
        assert mode in listed.stdout
        assert mode in help_result.stdout
        assert model.auxiliary in listed.stdout


@pytest.mark.parametrize('mode,stem', [('nast', 'exp_NAST'), ('yup', 'exp_YUP'), ('nares-2p', 'exp_Nares-2P')])
@pytest.mark.parametrize('ext', ['pdb', 'cif'])
def test_removed_modes_are_not_cli_choices(mode, stem, ext, tmp_path):
    result = run_cli(['-i', ROOT / 'CG_models' / f'{stem}.{ext}', '-m', mode, '-o', tmp_path / 'out'], tmp_path)
    assert result.returncode == 2 and 'invalid choice' in result.stderr
    assert mode not in MODELS and mode not in PARSERS
    assert not (tmp_path / 'out').exists()
    listed = run_cli(['--list-modes'], tmp_path)
    assert listed.returncode == 0
    assert mode not in listed.stdout


@pytest.mark.parametrize('mode', [mode for mode, model in MODELS.items()
                                 if 'cif' in model.formats and mode != 'generic'])
def test_equivalent_pdb_cif_and_mmcif(mode, tmp_path):
    model = MODELS[mode]
    outputs = []
    mmcif = tmp_path / 'native.mmcif'
    mmcif.write_bytes((ROOT / 'CG_models' / f'{model.stem}.cif').read_bytes())
    for i, source in enumerate((ROOT / 'CG_models' / f'{model.stem}.pdb', ROOT / 'CG_models' / f'{model.stem}.cif', mmcif)):
        directory = tmp_path / str(i)
        result = run_cli(['-i', source, '-m', mode, '-o', directory], tmp_path)
        assert result.returncode == 0, result.stderr
        s = PDBParser(QUIET=True).get_structure('output', directory / 'all_atom.pdb')
        outputs.append(np.array([a.coord for a in s.get_atoms()]))
    assert np.allclose(outputs[0], outputs[1], atol=.002, rtol=0)
    assert np.array_equal(outputs[1], outputs[2])


def test_oxrna_cli_requires_topology(tmp_path):
    result = run_cli(['-i', ROOT / 'CG_models/exp_oxRNA.dat', '-m', 'oxrna', '-o', tmp_path / 'out'], tmp_path)
    assert result.returncode == 2 and 'matching native topology' in result.stderr
    assert not (tmp_path / 'out').exists()


@pytest.mark.parametrize('mode', ['ifoldrna', 'toprna'])
def test_default_matches_previous_archive_profile(mode, tmp_path):
    args = ['-i', ROOT / 'CG_models' / (MODELS[mode].stem + '.pdb'), '-m', mode]
    default, explicit = tmp_path / 'default', tmp_path / 'explicit'
    for directory, options in [(default, []), (explicit, ['--bead-profile', 'archive'])]:
        result = run_cli([*args, '-o', directory, *options], tmp_path)
        assert result.returncode == 0, result.stderr
    for name in ('backbone.pdb', 'all_atom.pdb', 'fit_report.json'):
        assert (default / name).read_bytes() == (explicit / name).read_bytes()
    # Preserved results were generated before the default changed.
    for name in ('backbone.pdb', 'all_atom.pdb'):
        previous = ROOT / 'examples/reconstructed' / mode / 'pdb' / name
        assert (default / name).read_bytes() == previous.read_bytes()


def test_dat_is_not_silently_parsed_as_another_mode(tmp_path):
    result = run_cli(['-i', ROOT / 'CG_models/exp_oxRNA.dat', '-m', 'simrna', '-o', tmp_path / 'out'], tmp_path)
    assert result.returncode == 2 and 'does not interpret this file format' in result.stderr
    assert not (tmp_path / 'out').exists()


def test_custom_fragments_and_sequence(tmp_path):
    sequence = ''.join(r.base for r in PARSERS['simrna'](ROOT / 'CG_models/exp_SimRNA.pdb'))
    fasta = tmp_path / 'sequence.fasta'
    fasta.write_text('>archive_sequence\n' + sequence + '\n')
    args = ['-i', ROOT / 'CG_models/exp_SimRNA.pdb', '-m', 'simrna', '-o', tmp_path / 'out',
            '--fragments', ROOT / 'single_nucleotide_fragments', '--sequence', fasta]
    assert run_cli(args, tmp_path).returncode == 0
    fasta.write_text('>wrong\nA\n')
    assert 'does not match' in run_cli(args, tmp_path).stderr


def test_cli_does_not_overwrite_input_through_symlink(tmp_path):
    original = ROOT / 'CG_models/exp_SimRNA.pdb'
    local_input = tmp_path / 'input.pdb'
    local_input.write_bytes(original.read_bytes())
    output = tmp_path / 'out'
    output.mkdir()
    (output / 'all_atom.pdb').symlink_to(local_input)
    result = run_cli(['-i', local_input, '-m', 'simrna', '-o', output], tmp_path)
    assert result.returncode == 2 and 'overwrite an input' in result.stderr
    assert local_input.read_bytes() == original.read_bytes()
