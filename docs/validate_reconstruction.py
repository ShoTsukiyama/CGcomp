"""Reproducible installed-CLI geometry audit; never overwrite prior outputs.

Run: .venv/bin/python docs/validate_reconstruction.py before|after
No whole-structure superposition is performed.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cgcomp.fragments import load_fragments, BACKBONE
from cgcomp.models import MODELS
from cgcomp.parsers import mapped, spqr, oxrna
from cgcomp.parsers.generic import read_structure
from cgcomp.reconstruction import parse_input, kabsch

BONDS = [('P', "O5'"), ("O5'", "C5'"), ("C5'", "C4'"),
         ("C4'", "C3'"), ("C3'", "O3'"), ("C1'", 'N')]
BOND_TOLERANCE = .10  # Å, significant against unchanged fragment bonds (~1.4–1.6 Å).


def summary(values):
    v = np.asarray(values)
    return dict(count=len(v), mean=float(v.mean()), median=float(np.median(v)),
                maximum=float(v.max()))


def bond_lengths(base, atoms):
    return {f'{a}-{b}': float(np.linalg.norm(atoms[a] - atoms[
        ('N9' if base in 'AG' else 'N1') if b == 'N' else b])) for a, b in BONDS}


def rigid_atoms(mode, r, fragment):
    if mode == 'oxrna':
        origin, source = oxrna.fragment_proxies(r.base, fragment)
        rot, _, _ = kabsch(source, oxrna.target_proxies(r), center=False)
        shift = r.sites['COM'] * 8.518 - (origin + .5 * source[0]) @ rot
    else:
        if mode == 'spqr':
            source, target = spqr.frame_sites(r.base, fragment), spqr.target_frame(r)
        else:
            refs = mapped.reference_sites(mode, r.base, fragment)
            names = [n for n in r.sites if n in refs]
            source, target = [refs[n] for n in names], [r.sites[n] for n in names]
        rot, shift, _ = kabsch(source, target)
    return {n: xyz @ rot + shift for n, xyz in fragment.items()}


def geometry(mode, cg, rebuilt, fragments):
    details, errors = [], []
    for r, out in zip(cg, rebuilt):
        ideal = bond_lengths(r.base, fragments[r.base])
        actual = bond_lengths(r.base, out.sites)
        rigid = bond_lengths(r.base, rigid_atoms(mode, r, fragments[r.base]))
        deviations = {n: actual[n] - ideal[n] for n in ideal}
        details.append(dict(chain=r.chain, residue=r.number, base=r.base,
                            lengths_angstrom=actual, template_lengths_angstrom=ideal,
                            before_preservation_lengths_angstrom=rigid,
                            deviations_angstrom=deviations,
                            preservation_changes_angstrom={n: actual[n]-rigid[n] for n in ideal}))
        errors.extend(abs(x) for x in deviations.values())
    pairs = []
    if mode == 'oxrna':
        links = [(i, r.neighbors[0]) for i, r in enumerate(cg) if r.neighbors[0] != -1]
    else:
        links = [(i, i+1) for i, r in enumerate(cg[:-1])
                 if r.chain == cg[i+1].chain and cg[i+1].number == r.number+1]
    for i, j in links:
        pairs.append(dict(chain=cg[i].chain, residue=cg[i].number, next_residue=cg[j].number,
                          distance_angstrom=float(np.linalg.norm(rebuilt[i].sites["O3'"]-rebuilt[j].sites['P']))))
    distances = [p['distance_angstrom'] for p in pairs]
    return dict(bond_tolerance_angstrom=BOND_TOLERANCE, residues=details,
                absolute_bond_deviation=summary(errors),
                bonds_exceeding_tolerance=sum(e > BOND_TOLERANCE for e in errors),
                internucleotide={**summary(distances), 'exceeding_2_angstrom':sum(d>2 for d in distances),
                                 'exceeding_4_angstrom':sum(d>4 for d in distances), 'pairs':pairs})


def mapping(mode, cg, rebuilt, fragments):
    details, errors = [], []
    for r, out in zip(cg, rebuilt):
        if mode == 'oxrna':
            origin, proxies = oxrna.fragment_proxies(r.base, out.sites)
            target = oxrna.target_proxies(r)
            # Invert the official template registration; raw atomic axes/sites
            # are proxies, not an exact inverse of its mixed least-squares fit.
            _, template = oxrna.fragment_proxies(r.base, fragments[r.base])
            canonical, _, _ = kabsch(template, oxrna.target_proxies(type(r)(r.base, '1', 1, orientation=np.eye(3))), center=False)
            names = list(fragments[r.base])
            rot, _, _ = kabsch([fragments[r.base][n] for n in names], [out.sites[n] for n in names])
            axes = canonical.T @ rot
            delta = origin + .5*proxies[0] - r.sites['COM']*8.518
            physical_back = out.sites["C5'"] - r.sites['COM']*8.518
            physical_base = out.sites['N1' if r.base in 'AG' else 'N3'] - r.sites['COM']*8.518
            d = dict(com_registration_error_angstrom=float(np.linalg.norm(delta)),
                     registered_axis_errors=np.linalg.norm(axes-r.orientation, axis=1).tolist(),
                     raw_atomic_axis_errors=np.linalg.norm(proxies[:3]-target[:3], axis=1).tolist(),
                     interaction_site_errors_angstrom=np.linalg.norm(np.array([physical_back, physical_base])-target[3:], axis=1).tolist())
            errors.append(d['com_registration_error_angstrom'])
        elif mode == 'spqr':
            frame = spqr.frame_sites(r.base, out.sites)
            target = spqr.target_frame(r)
            d = {n: float(np.linalg.norm(a-b)) for n,a,b in zip(('BASE','XVEC','YVEC','ZVEC'), frame, target)}
            # SUGR is the native state-dependent ideal position, NOT a measured
            # sugar-ring centroid. Report that centroid only as a proxy residual.
            d['sugar_ring_proxy'] = float(np.linalg.norm(mapped.center(out.sites, mapped.SUGAR_RING)-r.sites['SUGR']))
            d['PHOS'] = float(np.linalg.norm(out.sites['P']-r.sites['PHOS']))
            if hasattr(spqr, 'reference_sites'):
                d['SUGR'] = float(np.linalg.norm(spqr.reference_sites(r.base, out.sites)['SUGR']-r.sites['SUGR']))
            errors.extend(d[n] for n in ('BASE','XVEC','YVEC','ZVEC'))
        else:
            refs = mapped.reference_sites(mode, r.base, out.sites)
            d = {n:float(np.linalg.norm(refs[n]-xyz)) for n,xyz in r.sites.items()}
            errors.extend(d.values())
        details.append(dict(chain=r.chain, residue=r.number, errors=d))
    return dict(primary_error_angstrom=summary(errors), residues=details)


def accuracy(cg, rebuilt, reference):
    assert len(cg) == len(reference)
    assert ''.join(r.base for r in cg) == ''.join(r.base for r in reference)
    all_errors, back_errors, details, residue_map = [], [], [], []
    for r, out, exp in zip(cg, rebuilt, reference):
        names = sorted(set(out.sites) & set(exp.sites))
        values = [float(np.sum((out.sites[n]-exp.sites[n])**2)) for n in names]
        backbone = [v for n,v in zip(names,values) if n in BACKBONE]
        all_errors.extend(values)
        back_errors.extend(backbone)
        details.append(dict(chain=r.chain, residue=r.number, exp_chain=exp.chain, exp_residue=exp.number,
                            matched_atoms=len(names), all_heavy_rmsd_angstrom=float(np.sqrt(np.mean(values))),
                            backbone_rmsd_angstrom=float(np.sqrt(np.mean(backbone)))))
        residue_map.append([r.chain,r.number,exp.chain,exp.number,r.base])
    return dict(all_heavy_rmsd_angstrom=float(np.sqrt(np.mean(all_errors))),
                backbone_rmsd_angstrom=float(np.sqrt(np.mean(back_errors))), matched_heavy_atoms=len(all_errors),
                matched_backbone_atoms=len(back_errors), residue_mapping=residue_map, per_residue=details,
                residues_above_3_angstrom=[d for d in details if d['all_heavy_rmsd_angstrom']>3],
                largest_errors=sorted(details,key=lambda d:d['all_heavy_rmsd_angstrom'],reverse=True)[:10])


def run(label):
    dest = ROOT/'examples/reconstruction_validation'/label
    dest.mkdir(parents=True, exist_ok=False)
    fragments, reference = load_fragments(), read_structure(ROOT/'CG_models/exp.pdb')
    results = {}
    for mode, model in MODELS.items():
        ext = model.formats[0]
        source = ROOT/'CG_models'/f'{model.stem}.{ext}'
        output = dest/mode
        options = ['--topology', ROOT/'examples/auxiliary/exp_oxRNA.top'] if mode=='oxrna' else []
        if mode in ('ifoldrna','toprna'):
            options += ['--bead-profile','archive']
        cmd = [str(Path(sys.executable).parent/'cgcomp'),'-i',str(source),'-m',mode,'-o',str(output),*map(str,options)]
        p = subprocess.run(cmd,cwd='/tmp',capture_output=True,text=True)
        assert p.returncode==0, p.stderr
        cg = parse_input(mode,source,ROOT/'examples/auxiliary/exp_oxRNA.top' if mode=='oxrna' else None)
        rebuilt = read_structure(output/'all_atom.pdb')
        backbone = read_structure(output/'backbone.pdb')
        assert [(r.base,r.chain,r.number) for r in cg] == [(r.base,r.chain,r.number) for r in rebuilt]
        assert all(set(r.sites)==set(fragments[r.base]) for r in rebuilt)
        assert all(set(r.sites)==set(BACKBONE) for r in backbone)
        assert all(np.array_equal(b.sites[n],a.sites[n]) for a,b in zip(rebuilt,backbone) for n in BACKBONE)
        fits = json.loads((output/'fit_report.json').read_text())['residues']
        fit_key = 'proxy_fit_rmsd' if mode=='oxrna' else 'fit_rmsd_angstrom'
        rms = accuracy(cg,rebuilt,reference)
        rms['applicable_reconstruction_accuracy'] = mode not in ('spqr','oxrna')
        rms['frame_note'] = ('Sequence-only mapping; independently initialized/simulated CG example, not derived from exp.pdb. Raw unsuperimposed displacement is NOT reconstruction accuracy.' if mode in ('spqr','oxrna') else 'Atomic forward mappings verified against exp.pdb; same original coordinate frame, no superposition.')
        results[mode] = dict(command=cmd, stdout=p.stdout, input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                             output_sha256={n:hashlib.sha256((output/n).read_bytes()).hexdigest() for n in ('all_atom.pdb','backbone.pdb','fit_report.json')},
                             fit_metric=fit_key, fit_residual=summary([f[fit_key] for f in fits]),
                             geometry=geometry(mode,cg,rebuilt,fragments), cg_mapping=mapping(mode,cg,rebuilt,fragments),
                             experimental_comparison=rms)
    report = ROOT/'docs/reconstruction_validation'/f'{label}.json'
    report.write_text(json.dumps(results,indent=2)+'\n')
    for mode,r in results.items():
        a,g = r['experimental_comparison'],r['geometry']
        print(mode, 'fit',round(r['fit_residual']['mean'],4),'backbone',round(a['backbone_rmsd_angstrom'],3),
              'all',round(a['all_heavy_rmsd_angstrom'],3),'O3-P',round(g['internucleotide']['median'],3),
              'bond max deviation',round(g['absolute_bond_deviation']['maximum'],3),
              'bad bonds',g['bonds_exceeding_tolerance'])

if __name__=='__main__':
    run(sys.argv[1])
