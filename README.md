# CGcomp

CGcomp reconstructs RNA backbone and all-atom structures from coarse-grained (CG) models. It adds missing heavy atoms by placing nucleotide fragments according to the input CG coordinates. CGcomp does not perform energy minimization or structural refinement.

## Installation

Clone the repository and install CGcomp in a Conda environment:

```bash
git clone https://github.com/USERNAME/CGcomp.git
cd CGcomp
conda env create -f environment.yml
conda activate cgcomp
pip install .
```

Replace `USERNAME` with the repository owner's GitHub username.

Check the installation:

```bash
cgcomp --help
```

## Usage

```bash
cgcomp -i INPUT -m MODEL -o OUTPUT_DIR
```

- `-i`, `--input`: Input CG structure file.
- `-m`, `--mode`: CG model name (see below).
- `-o`, `--output`: Output directory.

For example:

```bash
cgcomp -i input.pdb -m simrna -o results
```

For oxRNA, provide the corresponding topology file, which contains nucleotide identities and strand connectivity:

```bash
cgcomp -i input.dat -m oxrna --topology input.top -o results
```

To display available model names and their input requirements:

```bash
cgcomp --list-modes
```

### Supported CG models

| CG model | Mode | Input format |
|---|---|---|
| FebRNA | `febrna` | PDB, CIF, mmCIF |
| HiRE-RNA | `hire-rna` | PDB, CIF, mmCIF |
| iFoldRNA | `ifoldrna` | PDB, CIF, mmCIF |
| isRNA1 | `isrna1` | PDB, CIF, mmCIF |
| isRNA2 | `isrna2` | PDB, CIF, mmCIF |
| oxRNA | `oxrna` | DAT + topology file |
| RNA-JP | `rna-jp` | PDB, CIF, mmCIF |
| SimRNA | `simrna` | PDB, CIF, mmCIF |
| SPQR | `spqr` | PDB (native SPQR format) |
| TopRNA | `toprna` | PDB, CIF, mmCIF |
| VFold | `vfold` | PDB, CIF, mmCIF |
| Generic atomic input | `generic` | PDB, CIF, mmCIF |

The `generic` mode requires at least three non-collinear matched heavy-atom coordinates per nucleotide. SPQR reconstruction currently supports the A3 state only.

NAST, NARES-2P, and YUP are not supported for atomic reconstruction because their CG representations do not provide sufficient information to determine nucleotide orientations.

## Output

CGcomp writes three files to the selected output directory:

```text
results/
├── backbone.pdb
├── all_atom.pdb
└── fit_report.json
```

- **`backbone.pdb`**: Reconstructed RNA backbone heavy-atom coordinates.
- **`all_atom.pdb`**: Reconstructed RNA backbone, ribose, and nucleobase heavy-atom coordinates.
- **`fit_report.json`**: Reconstruction information and fitting diagnostics.

## Template attribution

Nucleotide fragments used by CGcomp originate from [Arena](https://github.com/pylelab/Arena). Check the source's redistribution terms before redistributing these fragment files.
