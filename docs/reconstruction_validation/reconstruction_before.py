"""Common row-vector Kabsch engine and model-specific fragment constraints."""
import numpy as np
from .parsers import generic, simrna, spqr, oxrna, mapped
from .models import MODELS

from .parsers.hire_rna import parse as parse_hire

PARSERS = {'generic': generic.read_structure, 'simrna': simrna.parse,
           'hire-rna': parse_hire, 'spqr': spqr.parse, **mapped.PARSERS}


def parse_input(mode, path, topology=None):
    if mode == 'oxrna':
        return oxrna.parse(path, topology)
    if topology:
        raise ValueError('--topology is only applicable to oxrna')
    return PARSERS[mode](path)


def kabsch(source, target, center=True):
    source, target = np.asarray(source, dtype=float), np.asarray(target, dtype=float)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3 or len(source) < 3:
        raise ValueError('At least three matching points are required')
    if not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ValueError('Fit coordinates must be finite')
    a, b = (source.mean(0), target.mean(0)) if center else (np.zeros(3), np.zeros(3))
    if np.linalg.matrix_rank(source - a, tol=1e-8) < 2 or np.linalg.matrix_rank(target - b, tol=1e-8) < 2:
        raise ValueError('Orientation cannot be determined: matching points are collinear')
    u, singular, vt = np.linalg.svd((source - a).T @ (target - b))
    rotation = u @ np.diag([1., 1., np.linalg.det(u @ vt)]) @ vt
    translation = b - a @ rotation
    rmsd = float(np.sqrt(np.mean(np.sum((source @ rotation + translation - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def reconstruct(residues, fragments, mode):
    if MODELS[mode].blocker:
        raise ValueError(MODELS[mode].blocker)
    output, diagnostics = [], []
    for r in residues:
        atoms = fragments[r.base]
        centered = True
        if mode == 'oxrna':
            origin, source = oxrna.fragment_proxies(r.base, atoms)
            target = oxrna.target_proxies(r)
            preserve = {}
            names = ['a1', 'a3', 'a2', 'RNA_BACK_PROXY', 'RNA_BASE_PROXY']
            centered = False
        elif mode == 'spqr':
            source, target = spqr.frame_sites(r.base, atoms), spqr.target_frame(r)
            preserve = {'P': r.sites['PHOS']}
            names = ['BASE', 'XVEC', 'YVEC', 'ZVEC']
        else:
            reference = mapped.reference_sites(mode, r.base, atoms)
            names = [n for n in r.sites if n in reference]
            unknown = [n for n in r.sites if n not in reference and not n.lstrip('0123456789').startswith('H')]
            if unknown:
                raise ValueError(f'Unmapped heavy atoms/beads {unknown} at {r.chain}:{r.number}')
            source = np.array([reference[n] for n in names])
            target = np.array([r.sites[n] for n in names])
            preserve = {n: xyz for n, xyz in mapped.real_atom_sites(mode, r).items() if n in atoms}
        try:
            rotation, translation, rmsd = kabsch(source, target, center=centered)
        except ValueError as e:
            raise ValueError(f'{r.chain}:{r.number} ({r.base}): {e}') from e
        if mode == 'oxrna':
            translation = r.sites['COM'] * oxrna.ANGSTROM_PER_UNIT - (origin + oxrna.COM_SHIFT * source[0]) @ rotation
        transformed = {n: xyz @ rotation + translation for n, xyz in atoms.items()}
        # Input real atoms are exact constraints in the final PDB. The reported
        # rigid-fit residual is measured BEFORE this required preservation.
        transformed.update({n: xyz.copy() for n, xyz in preserve.items()})
        out = generic.Nucleotide(r.base, r.chain, r.number, r.insertion, transformed)
        output.append(out)
        info = {'chain': r.chain, 'residue': r.number, 'base': r.base,
                'constraints': names, 'fit_rmsd_angstrom': rmsd,
                'mapping_profile': MODELS[mode].profile}
        if mode == 'oxrna':
            # Direction proxies have length one; the site proxies use Å, as in
            # the official converter. The combined residual is a mixed metric.
            del info['fit_rmsd_angstrom']
            fitted = source @ rotation
            info.update(proxy_fit_rmsd=rmsd,
                orientation_vector_residuals=np.linalg.norm(fitted[:3] - target[:3], axis=1).tolist(),
                interaction_site_residuals_angstrom=np.linalg.norm(fitted[3:] - target[3:], axis=1).tolist(),
                origin_convention='atomic geometric centroid + 0.5 Å * template a1',
                input_com_native=r.sites['COM'].tolist(), input_axes=r.orientation.tolist())
        elif mode != 'spqr':
            final_sites = mapped.reference_sites(mode, r.base, transformed)
            info['output_site_residuals_angstrom'] = {n: float(np.linalg.norm(final_sites[n] - r.sites[n])) for n in names}
        if mode == 'spqr':
            sugar = np.mean([transformed[n] for n in ("C1'", "C2'", "C3'", "C4'", "O4'")], axis=0)
            fitted_p = atoms['P'] @ rotation + translation
            info.update(sugar_center_deviation_angstrom=float(np.linalg.norm(sugar - r.sites['SUGR'])),
                        phosphate_preservation_shift_angstrom=float(np.linalg.norm(fitted_p - r.sites['PHOS'])))
        output[-1].state = r.state
        diagnostics.append(info)
    return output, diagnostics
