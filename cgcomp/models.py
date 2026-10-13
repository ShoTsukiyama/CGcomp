"""Supported reconstruction model catalog."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Model:
    name: str
    stem: str
    formats: tuple = ('pdb', 'cif', 'mmcif')
    auxiliary: str = 'none'
    profile: str = 'named atomic sites'


MODELS = {
    'febrna': Model('FebRNA', 'exp_FebRNA'),
    'hire-rna': Model('HiRE-RNA', 'exp_HiRE-RNA', profile='atomic backbone and ring COMs'),
    'ifoldrna': Model('iFoldRNA', 'exp_iFoldRNA', profile='archive export COM groups; differs from Ding2008 native groups'),
    'isrna1': Model('isRNA1', 'exp_isRNA1', profile='P/C4\' and grouped base COMs'),
    'isrna2': Model('isRNA2', 'exp_isRNA2', profile='P/C4\' and three grouped base COMs'),
    'oxrna': Model('oxRNA', 'exp_oxRNA', ('dat',), 'matching native topology', 'RNA-specific official atomic backmapping proxies'),
    'rna-jp': Model('RNA-JP', 'exp_RNA-JP'),
    'simrna': Model('SimRNA', 'exp_SimRNA'),
    'spqr': Model('SPQR', 'exp_SPQR', ('pdb',), profile='BASE/XVEC/YVEC frame, A3 rigid fragment; independent PHOS residual'),
    'toprna': Model('TopRNA', 'exp_TopRNA', profile='archive geometric centers; cytosine B includes exocyclic N4'),
    'vfold': Model('VFold', 'exp_VFold'),
    'generic': Model('Generic atomic RNA', 'exp', profile='all matched atomic sites; three non-collinear matches required'),
}
