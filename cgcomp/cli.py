import argparse
import json
from pathlib import Path
import sys
from .fragments import load_fragments
from .reconstruction import parse_input, reconstruct
from .models import MODELS
from .parsers.generic import read_fasta
from .pdb_writer import write_outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description='Reconstruct RNA using supplied single-nucleotide fragments. Supports 11 CG models and generic atomic fitting; --list-modes reports supported representations.')
    parser.add_argument('-i', '--input', help='Input PDB/mmCIF or native oxRNA DAT')
    parser.add_argument('-m', '--mode', choices=sorted(MODELS))
    parser.add_argument('-o', '--output', help='Output directory')
    parser.add_argument('--list-modes', action='store_true', help='List formats, auxiliary inputs, reconstruction status, and profiles')
    parser.add_argument('--topology', help='Matching native topology, required by oxrna')
    parser.add_argument('--bead-profile', choices=['archive'], help='Optional for ifoldrna/toprna: archive is the default RNAgrainy-compatible mapping; native publication definitions differ')
    parser.add_argument('--fragments', help='Directory containing A.pdb/U.pdb/G.pdb/C.pdb')
    parser.add_argument('--sequence', help='Optional FASTA used to verify the input sequence')
    args = parser.parse_args(argv)
    if args.list_modes:
        for mode, model in sorted(MODELS.items()):
            status = 'SUPPORTED'
            print(f'{mode:12} {status:11} formats={"/".join(model.formats):14} auxiliary={model.auxiliary}; {model.profile}')
        return 0
    if not all((args.input, args.mode, args.output)):
        parser.error('--input, --mode, and --output are required unless --list-modes is used')
    try:
        if args.bead_profile and args.mode not in ('ifoldrna', 'toprna'):
            raise ValueError('--bead-profile is only applicable to ifoldrna/toprna')
        if Path(args.input).suffix.lower().lstrip('.') not in MODELS[args.mode].formats:
            raise ValueError(f'{args.mode} accepts {"/".join(MODELS[args.mode].formats)}; it does not interpret this file format')
        residues = parse_input(args.mode, args.input, args.topology)
        if args.sequence and ''.join(read_fasta(args.sequence)) != ''.join(r.base for r in residues):
            raise ValueError('FASTA sequence does not match the input residue order')
        out = Path(args.output).resolve()
        protected = {Path(args.input).resolve()}
        if args.sequence:
            protected.add(Path(args.sequence).resolve())
        if args.topology:
            protected.add(Path(args.topology).resolve())
        if args.fragments:
            protected.update(p.resolve() for p in Path(args.fragments).glob('*.pdb'))
        if any((out / name).resolve() in protected for name in ('all_atom.pdb', 'backbone.pdb', 'fit_report.json')):
            raise ValueError('Output would overwrite an input file')
        rebuilt, diagnostics = reconstruct(residues, load_fragments(args.fragments), args.mode)
        write_outputs(rebuilt, out)
        (out / 'fit_report.json').write_text(json.dumps({'mode': args.mode, 'residues': diagnostics}, indent=2) + '\n')
        if args.bead_profile:
            print('Bead profile: ' + MODELS[args.mode].profile)
        print(f'Reconstructed {len(rebuilt)} residues: {out / "backbone.pdb"}, {out / "all_atom.pdb"}')
        return 0
    except (ValueError, OSError, KeyError) as e:
        print(f'cgcomp: {e}', file=sys.stderr)
        return 2
