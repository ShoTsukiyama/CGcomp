"""Complete archive model catalog; recognition is distinct from reconstruction."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Model:
    name: str
    stem: str
    formats: tuple = ('pdb', 'cif', 'mmcif')
    auxiliary: str = 'none'
    profile: str = 'named atomic sites'
    blocker: str = ''


MODELS = {
    'febrna': Model('FebRNA', 'exp_FebRNA'),
    'hire-rna': Model('HiRE-RNA', 'exp_HiRE-RNA', profile='atomic backbone and ring COMs'),
    'ifoldrna': Model('iFoldRNA', 'exp_iFoldRNA', auxiliary='--bead-profile archive', profile='archive export COM groups; differs from Ding2008 native groups'),
    'isrna1': Model('isRNA1', 'exp_isRNA1', profile='P/C4\' and grouped base COMs'),
    'isrna2': Model('isRNA2', 'exp_isRNA2', profile='P/C4\' and three grouped base COMs'),
    'nares-2p': Model('Nares-2P', 'exp_Nares-2P', profile='sugar/base centers and inter-residue phosphate midpoints',
        blocker='NARES-2P P is a midpoint between adjacent S sites, not an intra-nucleotide atomic landmark. S and B supply only two independent local points; rotation about S-B is unconstrained. A base normal or a third intra-nucleotide landmark is required.'),
    'nast': Model('NAST', 'exp_NAST', profile='one C3\' site per nucleotide',
        blocker='NAST supplies only C3\' per nucleotide and no orientation. Three rotational degrees of freedom remain undetermined; two additional non-collinear intra-nucleotide landmarks or an authoritative oriented frame are required.'),
    'oxrna': Model('oxRNA', 'exp_oxRNA', ('dat',), 'matching native topology', 'RNA-specific official atomic backmapping proxies'),
    'rna-jp': Model('RNA-JP', 'exp_RNA-JP'),
    'simrna': Model('SimRNA', 'exp_SimRNA'),
    'spqr': Model('SPQR', 'exp_SPQR', ('pdb',), profile='BASE/XVEC/YVEC frame, A3 rigid fragment; independent PHOS residual'),
    'toprna': Model('TopRNA', 'exp_TopRNA', auxiliary='--bead-profile archive', profile='archive geometric centers; cytosine B includes exocyclic N4'),
    'vfold': Model('VFold', 'exp_VFold'),
    'yup': Model('YUP', 'exp_YUP', profile='one P site per nucleotide',
        blocker='The supplied YUP export contains only one P site per nucleotide and no orientation. Three rotational degrees of freedom remain undetermined; two additional non-collinear intra-nucleotide landmarks or an authoritative oriented frame are required.'),
    'generic': Model('Generic atomic RNA', 'exp', profile='all matched atomic sites; three non-collinear matches required'),
}
