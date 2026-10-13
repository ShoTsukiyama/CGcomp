# CGcomp

**CGcomp** reconstructs RNA backbone and all-atom heavy-atom structures from coarse-grained (CG) RNA models. It supports **11 CG models** and a `generic` mode for structures containing known atomic coordinates.

CGcomp places nucleotide templates to complete missing atoms. It does **not** perform structural refinement or energy minimization.

## Installation

From the CGcomp source directory:

```bash
conda env create -f environment.yml
conda activate cgcomp
pip install .
```

Check the installation:

```bash
cgcomp --help
cgcomp --list-modes
```

## Usage

```bash
cgcomp -i INPUT -m MODEL -o OUTPUT_DIR
```

- `-i`, `--input`: input CG structure
- `-m`, `--mode`: CG model (see the table below)
- `-o`, `--output`: output directory

For example, to reconstruct a SimRNA structure:

```bash
cgcomp -i CG_models/exp_SimRNA.pdb -m simrna -o results/simrna
```

**oxRNA** requires the topology file corresponding to its DAT structure:

```bash
cgcomp -i CG_models/exp_oxRNA.dat -m oxrna \
  --topology examples/auxiliary/exp_oxRNA.top -o results/oxrna
```

For all available options, run `cgcomp --help`.

## Supported CG models

| CG model | Mode | Input format |
|---|---|---|
| FebRNA | `febrna` | PDB, CIF/mmCIF |
| HiRE-RNA | `hire-rna` | PDB, CIF/mmCIF |
| iFoldRNA | `ifoldrna` | PDB, CIF/mmCIF |
| isRNA1 | `isrna1` | PDB, CIF/mmCIF |
| isRNA2 | `isrna2` | PDB, CIF/mmCIF |
| oxRNA | `oxrna` | DAT + topology file |
| RNA-JP | `rna-jp` | PDB, CIF/mmCIF |
| SimRNA | `simrna` | PDB, CIF/mmCIF |
| SPQR | `spqr` | PDB (native SPQR format) |
| TopRNA | `toprna` | PDB, CIF/mmCIF |
| VFold | `vfold` | PDB, CIF/mmCIF |
| Generic atomic coordinates | `generic` | PDB, CIF/mmCIF |

The `generic` mode requires at least three non-collinear matching heavy-atom positions per nucleotide. SPQR reconstruction supports the **A3 state only**. The iFoldRNA and TopRNA modes use bead definitions compatible with the supplied RNAgrainy-generated examples; other bead definitions are not necessarily interchangeable.

## Output files

CGcomp writes the following files to the specified output directory:

```text
OUTPUT_DIR/
├── backbone.pdb
├── all_atom.pdb
└── fit_report.json
```

- **`backbone.pdb`**: reconstructed RNA backbone heavy atoms
- **`all_atom.pdb`**: reconstructed RNA backbone, ribose, and nucleobase heavy atoms
- **`fit_report.json`**: information about template fitting and reconstruction

Hydrogen atoms are not included.

## Template attribution

The nucleotide fragments used by CGcomp originate from [Arena](https://github.com/pylelab/Arena). Please check the original source and its redistribution terms before sharing the fragment files.

## Citation
If you use CGcomp in your research, please cite the following paper:
Tsukiyama, S., Li, Y., Sato, K., Kurata, H., & Zhang, Y. (2026). RNArefine: AI-guided Atomic-Level Refinement of RNA Structures. bioRxiv. https://doi.org/10.64898/2026.06.26.734804
