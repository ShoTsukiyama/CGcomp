"""Structural regressions: distinguish rigid fitting from exact-atom restoration."""
import importlib.util
import json
import numpy as np
import pytest
from cgcomp.models import MODELS
from cgcomp.parsers import mapped, spqr, oxrna
from cgcomp.parsers.generic import Nucleotide
from cgcomp.reconstruction import parse_input, reconstruct
from conftest import ROOT

spec = importlib.util.spec_from_file_location('geometry_audit', ROOT/'docs/validate_reconstruction.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
SUPPORTED = [m for m, model in MODELS.items() if not model.blocker]


def actual_input(mode):
    model = MODELS[mode]
    return parse_input(mode, ROOT/'CG_models'/f'{model.stem}.{model.formats[0]}',
                       ROOT/'examples/auxiliary/exp_oxRNA.top' if mode == 'oxrna' else None)


@pytest.mark.parametrize('mode', SUPPORTED)
def test_actual_intranucleotide_geometry_and_preservation_causality(mode, fragments, record_property):
    cg = actual_input(mode)
    rebuilt, _ = reconstruct(cg, fragments, mode)
    result = audit.geometry(mode, cg, rebuilt, fragments)
    record_property('bond_deviation_summary', json.dumps(result['absolute_bond_deviation']))
    record_property('bonds_deviating_over_0.1_A', result['bonds_exceeding_tolerance'])
    for r, out, detail in zip(cg, rebuilt, result['residues']):
        rigid = audit.rigid_atoms(mode, r, fragments[r.base])
        preserved = {} if mode in ('spqr', 'oxrna') else mapped.real_atom_sites(mode, r)
        # Keep real atomic constraints in models that require them. Virtual beads
        # must never be assigned directly to atoms just because names coincide.
        for name, xyz in out.sites.items():
            assert np.allclose(xyz, preserved.get(name, rigid[name]), atol=1e-10), (mode,r.number,name)
        expected_atoms = {n: preserved.get(n, xyz) for n, xyz in rigid.items()}
        expected_lengths = audit.bond_lengths(r.base, expected_atoms)
        for bond, actual in detail['lengths_angstrom'].items():
            assert actual == pytest.approx(expected_lengths[bond], abs=1e-10), (mode,r.number,bond)
            assert detail['before_preservation_lengths_angstrom'][bond] == pytest.approx(
                detail['template_lengths_angstrom'][bond], abs=1e-10)
        if mode in ('spqr', 'oxrna', 'ifoldrna'):
            assert max(abs(v) for v in detail['deviations_angstrom'].values()) < 1e-10, (mode,r.number,detail)


@pytest.mark.parametrize('mode', SUPPORTED)
def test_actual_internucleotide_diagnostics(mode, fragments, record_property):
    cg = actual_input(mode)
    rebuilt, _ = reconstruct(cg, fragments, mode)
    result = audit.geometry(mode, cg, rebuilt, fragments)['internucleotide']
    assert result['count'] == 101
    assert np.isfinite([result[n] for n in ('mean','median','maximum')]).all()
    assert 0 <= result['exceeding_4_angstrom'] <= result['exceeding_2_angstrom'] <= 101
    record_property('O3_P_diagnostic', json.dumps({k:v for k,v in result.items() if k!='pairs'}))
    # Distances are diagnostic, never an acceptance criterion for bond closure.


@pytest.mark.parametrize('mode', SUPPORTED)
def test_actual_forward_cg_mapping(mode, fragments, record_property):
    cg = actual_input(mode)
    rebuilt, fits = reconstruct(cg, fragments, mode)
    report = audit.mapping(mode, cg, rebuilt, fragments)
    record_property('CG_mapping', json.dumps(report['primary_error_angstrom']))
    for index, (r, out, fit) in enumerate(zip(cg, rebuilt, fits)):
        if mode == 'oxrna':
            origin, proxies = oxrna.fragment_proxies(r.base, out.sites)
            assert np.allclose(origin + .5*proxies[0], r.sites['COM']*8.518, atol=1e-10)
            # Exact inversion of template registration, distinct from approximate
            # physical atomic proxies. Detects order, handedness and translation.
            info = report['residues'][index]['errors']
            assert max(info['registered_axis_errors']) < 1e-10
            assert np.allclose(np.linalg.norm(proxies[:3]-oxrna.target_proxies(r)[:3],axis=1),
                               fit['orientation_vector_residuals'],atol=1e-10)
        elif mode == 'spqr':
            sites = spqr.reference_sites(r.base, out.sites)
            for n in ('BASE','XVEC','YVEC'):
                # Native coordinates are rounded; target_frame orthonormalizes.
                assert np.linalg.norm(sites[n]-r.sites[n]) < .02
            assert np.allclose(spqr.frame_sites(r.base, out.sites), spqr.target_frame(r), atol=1e-10)
            # Native MC_map_sugar uses raw endpoint vectors, including the
            # small nonorthonormality in this supplied simulation. Verify that
            # definition independently, then measure normalized-frame residual.
            center = r.sites['BASE']
            x, y = r.sites['XVEC']-center, r.sites['YVEC']-center
            offset = np.array([-1.56584,-4.48467,.709821] if r.base in 'AG' else [1.05136,-3.42329,.652454])
            native_sugar = center + offset @ np.array([x,y,np.cross(x,y)])
            # Native output uses variable decimal precision (down to 0.001 Å
            # at large coordinates). Subtracting endpoints and multiplying by
            # a ~5 Å lever arm amplifies that rounding; allow 0.01 Å here.
            assert np.linalg.norm(native_sugar-r.sites['SUGR']) < .01
            ideal_sugar = spqr.target_frame(r)[0] + offset @ (spqr.target_frame(r)[1:]-center)
            assert np.allclose(sites['SUGR'],ideal_sugar,atol=1e-10)
            for n, xyz in r.sites.items():
                assert np.linalg.norm(sites[n]-xyz) == pytest.approx(fit['output_site_residuals_angstrom'][n],abs=1e-10)
        else:
            sites = mapped.reference_sites(mode, r.base, out.sites)
            for n, xyz in r.sites.items():
                assert np.linalg.norm(sites[n]-xyz) == pytest.approx(fit['output_site_residuals_angstrom'][n],abs=1e-10)
                if n in mapped.real_atom_sites(mode,r) and n in out.sites:
                    assert np.array_equal(out.sites[n],xyz)


@pytest.mark.parametrize('base', list('AUGC'))
def test_spqr_independent_phosphate_does_not_break_fragment(base, fragments):
    atoms = fragments[base]
    sites = spqr.reference_sites(base, atoms)
    displaced = {n:xyz.copy() for n,xyz in sites.items()}
    displaced['PHOS'] += np.array([1.,-.7,.4])
    r = Nucleotide(base,'A',1,sites=displaced,state='A3')
    out, fits = reconstruct([r],fragments,'spqr')
    names = list(atoms)
    source, target = np.array([atoms[n] for n in names]), np.array([out[0].sites[n] for n in names])
    # Test ALL atom-pair distances, including P–OP1/OP2, not just atom counts.
    assert np.allclose(np.linalg.norm(source[:,None]-source[None,:],axis=2),
                       np.linalg.norm(target[:,None]-target[None,:],axis=2),atol=1e-10)
    assert fits[0]['phosphate_site_deviation_angstrom'] == pytest.approx(np.linalg.norm([1.,-.7,.4]))
    assert fits[0]['output_site_residuals_angstrom']['SUGR'] < 1e-10
    assert np.allclose(spqr.frame_sites(base,out[0].sites),spqr.target_frame(r),atol=1e-10)


@pytest.mark.parametrize('mode', ['nares-2p','nast','yup'])
def test_missing_orientation_stays_blocked_for_geometry(mode, fragments):
    with pytest.raises(ValueError, match='unconstrained|undetermined'):
        reconstruct(actual_input(mode),fragments,mode)


@pytest.mark.parametrize('mode', SUPPORTED)
def test_experimental_residue_mapping_and_rmsd_diagnostics(mode, fragments, record_property):
    from cgcomp.parsers.generic import read_structure
    cg = actual_input(mode)
    reference = read_structure(ROOT/'CG_models/exp.pdb')
    rebuilt, _ = reconstruct(cg,fragments,mode)
    result = audit.accuracy(cg,rebuilt,reference)
    assert result['matched_heavy_atoms'] == 2181
    assert result['matched_backbone_atoms'] == 1224
    assert len(result['per_residue']) == 102
    assert all(a.base==b.base for a,b in zip(cg,reference))
    if mode == 'spqr':
        assert [(r.chain,r.number) for r in cg] == [('0',i) for i in range(102)]
    elif mode == 'oxrna':
        assert [(r.chain,r.number) for r in cg] == [('1',i) for i in range(1,103)]
    else:
        assert [(r.chain,r.number) for r in cg] == [(r.chain,r.number) for r in reference]
        # The same-frame forward mapping from exp.pdb has independent coverage
        # tests; here ensure RMSD uses atom correspondences without superposition.
        if mode == 'generic':
            assert result['all_heavy_rmsd_angstrom'] == 0
    record_property('comparison_applicable', mode not in ('spqr','oxrna'))
    record_property('unsuperimposed_all_heavy_RMSD_A',result['all_heavy_rmsd_angstrom'])
    record_property('unsuperimposed_backbone_RMSD_A',result['backbone_rmsd_angstrom'])


@pytest.mark.parametrize('base', list('AUGC'))
def test_supplied_spqr_fragments_are_anti_c3_endo(base, fragments):
    atoms = fragments[base]
    def dihedral(names):
        a,b,c,d = [atoms[n] for n in names]
        u,v,w = b-a,c-b,d-c
        p,q = np.cross(u,v),np.cross(v,w)
        return np.degrees(np.arctan2(np.dot(v,np.cross(p,q))/np.linalg.norm(v),np.dot(p,q)))
    chi = dihedral(["O4'","C1'",'N9' if base in 'AG' else 'N1','C4' if base in 'AG' else 'C2'])
    delta = dihedral(["C5'","C4'","C3'","O3'"])
    # Official pdb2spqr.py class_glp thresholds; these fragments are A3.
    assert -180 <= chi <= -120 or 155 < chi <= 180
    assert 0 < delta < 100
    cg = actual_input('spqr')
    assert all(r.state=='A3' for r in cg)
    out, _ = reconstruct(cg,fragments,'spqr')
    assert all(r.state=='A3' for r in out)
