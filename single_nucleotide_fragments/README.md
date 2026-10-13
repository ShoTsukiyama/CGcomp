# Arena single-nucleotide fragments

These files contain the original single-nucleotide heavy-atom fragments used by Arena. They were extracted from `../IdealRNA/ACGU.rna.pdb` according to the residue-name and chain-selection logic in `../IdealRNA/createIdealRNA2.py`, as confirmed against `../IdealRNA.hpp` and `../MissingRNAatom.hpp`.

| Output | Source residue | Heavy atoms |
| --- | --- | ---: |
| `A.pdb` | chain A, residue 1 (A) | 22 |
| `C.pdb` | chain A, residue 2 (C) | 20 |
| `G.pdb` | chain B, residue 7 (G) | 23 |
| `U.pdb` | chain B, residue 8 (U) | 20 |

## Extraction procedure

For each source residue, every `ATOM` record whose element was not hydrogen was selected. No records from neighboring residues were included. Atom names, element symbols, occupancies, temperature factors, and the x, y, and z coordinate fields were preserved exactly. Only atom serial numbers, chain IDs, and residue numbers were normalized for standalone files: each fragment is chain A, residue 1, with atom serials beginning at 1. `TER` and `END` records were added.

Validation compared each output atom, in order, with its source residue. All atom-name lists and coordinate strings match the source exactly, each file contains one residue, and no hydrogen atoms are present.

## Arena source and publication

Original repository: [pylelab/Arena](https://github.com/pylelab/Arena)

Zion R. Perry, Anna Marie Pyle, and Chengxin Zhang (2023), “Arena: rapid and accurate reconstruction of full atomic RNA structures from coarse-grained models,” *Journal of Molecular Biology*. DOI: [10.1016/j.jmb.2023.168210](https://doi.org/10.1016/j.jmb.2023.168210).

## License and redistribution

No Arena-wide license file or redistribution notice was found in this source distribution. The bundled `INF.py` and `pstream.h` files contain their own GPL and Boost Software License notices, respectively; those notices apply to those files and do not state a license for Arena or `IdealRNA/ACGU.rna.pdb`.
