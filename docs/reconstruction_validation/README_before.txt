# CGcomp

CGcomp reconstructs RNA backbone, ribose, and nucleobase heavy atoms by rigidly
superimposing the supplied **single-nucleotide Arena fragments**. All supported
models use the same four templates and Kabsch engine. This independent Python
package requires no RNArefine, other reconstruction software, dynamics, or refinement.

Version 1.1 accounts for every model in `CG_models.zip`: eleven CG models pass
installed-CLI reconstruction tests; three are recognized and parsed but lack
the geometric information needed to orient individual fragments. Generic atomic
fitting is also supported. Unsupported entries never produce fabricated output.

## Complete inventory and coverage

The recursively inspected archive contains **14 distinct CG models**, **27
structural files**, and a full-atom reference (`exp.pdb`). All examples have
102 nucleotides. There are no topology or FASTA files inside the ZIP. Twelve
CG models have PDB/CIF pairs; SPQR supplies native PDB, and oxRNA native DAT.

| CG model / mode | Bead representation | Input format | Additional inputs | Reconstruction method | Test status |
|---|---|---|---|---|---|
| FebRNA / `febrna` | P/C4′/N9(A/G) or N1(C/U), atomic | PDB/CIF/MMCIF | None | Three-atom Kabsch; preserve atoms | Passed |
| HiRE-RNA / `hire-rna` | Five atomic backbone sites plus B1/B2 ring COMs | PDB/CIF/MMCIF | None | Atomic/virtual-site Kabsch | Passed |
| iFoldRNA / `ifoldrna` | Archive P/S/B COM groups, detailed below | PDB/CIF/MMCIF | `--bead-profile archive` | COM-site Kabsch | Passed, explicit archive profile |
| isRNA1 / `isrna1` | P and S=C4′; 3 purine or 2 pyrimidine COMs | PDB/CIF/MMCIF | None | Atomic/virtual-site Kabsch | Passed |
| isRNA2 / `isrna2` | P and S=C4′; 3 base COMs for every nucleotide | PDB/CIF/MMCIF | None | Atomic/virtual-site Kabsch | Passed |
| Nares-2P / `nares-2p` | S/B centers; P is an adjacent-sugar midpoint | PDB/CIF/MMCIF | Missing base normal or third local landmark | Local S/B axis leaves axial rotation free | Unsupported; rejection tested |
| NAST / `nast` | One C3′ point per nucleotide | PDB/CIF/MMCIF | Missing nucleotide orientation | One point leaves 3 rotations free | Unsupported; rejection tested |
| oxRNA / `oxrna` | COM, a1/a3 axes, strand connectivity | DAT | Matching native `--topology` | Official RNA atomic proxies; shared Kabsch | Passed |
| RNA-JP / `rna-jp` | P/C4′ + N9/C6/C2(A/G) or N1/C2/C4(C/U), atomic | PDB/CIF/MMCIF | None | Five-atom Kabsch; preserve atoms | Passed |
| SimRNA / `simrna` | Five named atomic sites listed above | PDB/CIF/MMCIF | None | Five-atom Kabsch; preserve atoms | Passed |
| SPQR / `spqr` | BASE/XVEC/YVEC, virtual SUGR, PHOS | Native PDB | A3 state | Oriented base-frame Kabsch; preserve P | Passed, A3 only |
| TopRNA / `toprna` | P atom, geometric S/B centers; cytosine B includes N4 | PDB/CIF/MMCIF | `--bead-profile archive` | Geometric-site Kabsch | Passed, explicit archive profile |
| VFold / `vfold` | P/C4′/N9(A/G) or N1(C/U), atomic | PDB/CIF/MMCIF | None | Three-atom Kabsch; preserve atoms | Passed |
| YUP / `yup` | One P point per nucleotide | PDB/CIF/MMCIF | Missing nucleotide orientation | One point leaves 3 rotations free | Unsupported; rejection tested |
| Atomic reference / `generic` | Any matched template heavy atoms | PDB/CIF/MMCIF | ≥3 non-collinear matches | All matched atoms enter Kabsch | Passed |

`--mode` recognizes every discovered model. The three unsupported entries give
precise technical errors and create no atomic output. `cgcomp --list-modes`
displays formats, required auxiliary information, profiles, and blockers.
Recognition does not claim successful reconstruction; no placeholder reconstruction
or fallback to another model exists.

`docs/archive_inventory.json` lists every ZIP member, including macOS metadata.
`docs/input_inventory.json` records all extracted files, checksums, sequences,
identifiers, and bead names. `docs/model_inventory.json` groups their files and
records atomic/virtual sites, orientations, auxiliary requirements, and mapping
checks. Original archives and datasets are unchanged.

## Installation

Python 3.10+, NumPy, and Biopython are required; setuptools builds the package.
The wheel includes unchanged fragment files, so execution does not depend on
the current directory or extracted library.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install .
cgcomp --help
cgcomp --list-modes
```

Install tests with `pip install '.[test]'`, then run `python -m pytest -q`.
On the offline development host, existing Anaconda dependencies were exposed
through a workspace venv created with `--system-site-packages`. Installation was
tested with `pip install --no-build-isolation --no-deps .`. System environments
and RNArefine were not modified. The release wheel is
`dist/cgcomp-1.1.0-py3-none-any.whl`.

## CLI and tested examples

```text
cgcomp -i INPUT -m MODEL -o OUTPUT_DIR
cgcomp --input INPUT --mode MODEL --output OUTPUT_DIR
cgcomp --list-modes

--topology FILE           Matching native topology, required by oxrna
--bead-profile archive    Required by ifoldrna/toprna; supplied export definitions
--fragments DIR           Custom A.pdb/U.pdb/G.pdb/C.pdb library
--sequence FILE           Optional FASTA verifying identity and input order
--list-modes              Formats, auxiliary inputs, profiles, reconstruction status
-h, --help               Usage
```

Both long and short input/mode/output options are accepted. DAT is interpreted
only by the selected native mode, never automatically as oxRNA. PDB/mmCIF parsers
preserve bead names and coordinates internally. Atomic aliases O1P/O2P and sugar
`*` suffixes normalize to OP1/OP2 and `'`; RA/RU/RG/RC normalize to A/U/G/C.
Residue numbers, insertion codes, one-character chain identifiers, sequence, and
input order are retained. The first coordinate model and blank/A alternate
locations are read; duplicate ambiguous sites are errors. FASTA verifies the
concatenated sequence and does not invent residues or override identities.

Run these commands from the project root after installation. Every command has
been executed through the installed entry point. Both archived PDB/CIF versions
of every applicable supported model were additionally run; the full command
manifest is `docs/example_commands.json`.

```bash
cgcomp -i CG_models/exp_FebRNA.pdb -m febrna -o examples/reconstructed/febrna/pdb
cgcomp -i CG_models/exp_HiRE-RNA.cif -m hire-rna -o examples/reconstructed/hire-rna/cif
cgcomp -i CG_models/exp_iFoldRNA.pdb -m ifoldrna --bead-profile archive -o examples/reconstructed/ifoldrna/pdb
cgcomp -i CG_models/exp_isRNA1.pdb -m isrna1 -o examples/reconstructed/isrna1/pdb
cgcomp -i CG_models/exp_isRNA2.pdb -m isrna2 -o examples/reconstructed/isrna2/pdb
cgcomp -i CG_models/exp_oxRNA.dat -m oxrna --topology examples/auxiliary/exp_oxRNA.top -o examples/reconstructed/oxrna/dat
cgcomp -i CG_models/exp_RNA-JP.pdb -m rna-jp -o examples/reconstructed/rna-jp/pdb
cgcomp -i CG_models/exp_SimRNA.pdb -m simrna -o examples/reconstructed/simrna/pdb
cgcomp -i CG_models/exp_SPQR.pdb -m spqr -o examples/reconstructed/spqr/pdb
cgcomp -i CG_models/exp_TopRNA.pdb -m toprna --bead-profile archive -o examples/reconstructed/toprna/pdb
cgcomp -i CG_models/exp_VFold.pdb -m vfold -o examples/reconstructed/vfold/pdb
cgcomp -i CG_models/exp.pdb -m generic -o examples/reconstructed/generic/pdb
```

Generic mode also passes on PDB/CIF versions of VFold, FebRNA, RNA-JP, and
SimRNA. It reads each atom combination from the input rather than fixing it.

## Shared reconstruction method

Each residue selects its supplied A/U/G/C fragment. Atomic sites, virtual
centers, or orientation proxies enter the common proper-rotation Kabsch SVD
engine. Every template heavy atom receives the resulting rigid transformation.
Reflections and underdetermined fits are rejected. Sites explicitly representing
real atoms replace their fitted positions exactly; virtual sites never replace
individual atoms. Templates remain unchanged; no new ideal templates, base-pair
fragments, or stacking fragments are generated.

Fit RMSD is measured **before** real-atom preservation. Deformed CG geometry
can differ from the one template available for each base, so fitting residuals
are expected and atom preservation can change local geometry. Final forward-mapped
site residuals are reported too. Independent nucleotide placement performs no
inter-residue optimization or refinement.

### Atomic models and grouped COMs

VFold/FebRNA use all three named atomic sites; SimRNA/RNA-JP use all five.
Generic uses every matched known heavy atom and requires three non-collinear
matches. Unknown heavy-atom names cause an error instead of an inferred mapping.
The definitions are supported by the [VFold scaffold description](https://rnajournal.cshlp.org/content/18/4/610.full.html),
[FebRNA paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9515226/), and
[RNA-JP author manual](https://rna.physics.missouri.edu/RNAJP/RNAJP_User_Manual.pdf),
and were checked against supplied `exp.pdb`.

COMs use C=12.011, N=14.007, O=15.999, P=30.974. HiRE-RNA preserves
P/O5′/C5′/C4′/C1′. B1 is the COM of N9/C8/N7/C5/C4 for A/G, or
N1/C2/N3/C4/C5/C6 for C/U; purine B2 is the six-member-ring COM.
This follows the [HiRE-RNA definition](https://pmc.ncbi.nlm.nih.gov/articles/PMC5529179/).

isRNA1/isRNA2 preserve P and C4′, despite the latter bead's name S. Base
beads are mass-weighted centers of these groups:

| Model | Bead(s) | Heavy-atom groups |
|---|---|---|
| Both | R_C / R_1 | N9/C8/N7/C5/C4/N3 |
| isRNA1 | A_C; A_N | C6/N6; C2/N1 |
| isRNA1 | G_O; G_N | C6/O6; N1/C2/N2 |
| isRNA1 | Y_C | N1/C5/C6 |
| isRNA1 | C_N; U_O | C2/O2/N3/C4/N4; C2/O2/N3/C4/O4 |
| isRNA2 | A_1; A_2 | C6/N6; C2/N1 |
| isRNA2 | G_1; G_2 | C6/O6/N1; C2/N2 |
| isRNA2 | Y_1; Y_2 | N1/C5/C6; C2/O2 |
| isRNA2 | C_1; U_1 | N3/C4/N4; N3/C4/O4 |

[isRNA1 methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8064582/) establish its
atomic backbone and grouped-COM bases; [isRNA2 Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC9731381/)
gives its groups. Supplied isRNA1 groups were identified against the full-atom
reference, consistently with its grouped-COM representation. Every forward-mapped
HiRE/isRNA site of every archived residue agrees within 0.0015 Å. Force-field
mass assignments do not change a site's explicitly defined coordinate location.

### Explicit archive profiles: iFoldRNA and TOPRNA

These exports differ from native publications. Their exact groups were identified
by forward-mapping all 102 supplied reference residues, reproducing every bead
within 0.0015 Å. `--bead-profile archive` explicitly selects those definitions
and records them in the report; it is required because P/S/B names alone do not
distinguish profiles.

* iFoldRNA archive P is the COM of P/OP1/OP2. S uses
  C1′/C2′/C3′/C4′/O4′/C5′/O3′/O2′. Purine B uses the nine ring atoms
  N9/C8/N7/C5/C6/N1/C2/N3/C4; pyrimidine B uses
  N1/C2/N3/C4/C5/C6/**O2**. All three are virtual centers.
* TOPRNA archive P is atomic P. S geometrically averages
  C1′/C2′/C3′/C4′/O4′/C5′. B geometrically averages the purine ring atoms
  above or the six pyrimidine ring atoms, adding **N4 for cytosine**.

The [original iFoldRNA/Ding2008 model](https://pmc.ncbi.nlm.nih.gov/articles/PMC2390798/)
uses a five-member sugar-ring COM and six-member base-ring center. The
[original TOPRNA methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC3983386/)
specify cyclic N/C averages, excluding exocyclic N4. Native outputs using those
different definitions are **not** validated by these archive profiles and must
not be passed as `archive`. No different CG model is silently substituted.

### SPQR

A dedicated parser handles SPQR's native nonstandard PDB coordinate fields.
BASE averages C2/C4/C6. X points toward C2; Y is orthogonalized toward C6 for
A/G or C4 for C/U; Z=X×Y. XVEC/YVEC are absolute endpoints one Å from BASE.
These form four frame constraints for Kabsch, reproducing position/orientation
to numerical precision before rounding. PHOS is preserved as P; virtual SUGR
is never assigned to a sugar atom. Its deviation from the fitted ring center
and the displacement needed to preserve phosphate are reported.

All supplied states are A3 (anti/C3′-endo). Other states are rejected because
matching templates are absent. A single fragment cannot generally fit frame,
sugar, and phosphate exactly together. The maximum supplied sugar deviation is
0.166451 Å and phosphate preservation shift 1.460181 Å. Mapping follows the
[SPQR paper](https://arxiv.org/html/1711.00741) and official `pdb2spqr.py`.

### oxRNA

Native DAT stores position, a1, a3, velocity, and angular velocity (15 values/body).
The parser retains native coordinates and axes, with a2=a3×a1. One configuration
is accepted. Topology provides sequence/connectivity; both `N S 5->3` strand
format (`type=RNA`) and legacy `strand base n3 n5` are parsed, with reciprocal
connectivity checks. Storage order is retained.

The ZIP omits topology. `examples/auxiliary/exp_oxRNA.top` is an unchanged authentic
companion recovered from the neighboring development run. That run's DAT is
byte-identical to the archive (SHA-256
`f38347e68809333c4d96b496c80dac86a74109dede7cfe73c45d0d4b7e71306c`).
The topology declares 102 nucleotides, one 5′→3′ strand, and `type=RNA`.
Sources/checksums are recorded in `docs/auxiliary_provenance.json`; the external
topology is not described as archive content.

Atomic registration follows the RNA branch of the
[official backmapper](https://github.com/lorenzo-rovigatti/oxDNA/blob/master/analysis/src/oxDNA_analysis_tools/oxDNA_PDB.py).
Its template origin is the geometric atom centroid. a1 points toward mean pairing-edge
atoms (C6/N1/C2 for A/G; C4/N3/C2 for C/U); a3 averages ring normals signed
toward O4′. Proxies comprise a1/a3/a2, relative C5′, and N1(A/G) or N3(C/U).
Targets include native axes and RNA-specific sites 8.518×(−0.4a1+0.2a3) and
8.518×0.4a1 Å. Kabsch fits without centering; translation scales native position
by 8.518 Å/unit and retains the source-defined empirical 0.5 Å template-a1 shift.
No DNA backbone convention or invented offset is used. Direction proxies and
Å site vectors give a mixed `proxy_fit_rmsd`; separate vector/site residuals are
reported. COM registration follows this backmapping origin convention, not a
claim that the atomic centroid physically equals COM. The native center is never
assigned to P/C4′. Official-source reference-frame data, independently computed
from the supplied fragments, verify the mapping in `docs/oxrna_reference_frames.json`.
External atomic templates are not used.

## Specific technical blockers

* **YUP:** only P is supplied per nucleotide; no orientation or other local
  coordinates exist. All rotations about that point remain possible.
* **NAST:** only C3′ is supplied; three rotational degrees of freedom remain.
  Neighbor-derived tangents would impose an unsupported base-orientation rule.
* **NARES-2P:** P for residues 2–102 agrees with the midpoint of previous/current
  S within 0.001 Å. Local S/B fix an axis but allow rotation around it, preserving
  those centers and every adjacent midpoint. The terminal first P cannot orient
  the other 101 nucleotides. The supplied S averages nine sugar-region atoms
  (C1′/C2′/C3′/C4′/O4′/C5′/O5′/O3′/O2′); B geometrically averages purine ring
  atoms or all eight pyrimidine base atoms. This sugar export also differs from
  the native ring-only definition. The [NARES-2P paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4608961/)
  confirms that P is an adjacent-sugar midpoint, not atomic phosphorus.

YUP/NAST need an authoritative oriented frame or two additional non-collinear
local landmarks. NARES needs a base normal or third independent local landmark.
Sequence/topology cannot provide missing geometry. Reconstruction never borrows
coordinates from `exp.pdb`; that reference audits mappings and tests only.
All three entries have actual-input CLI rejection tests, and a rotation-invariance
test proves the NARES ambiguity. None is marked as successful reconstruction.

## Outputs and validation

```text
examples/reconstructed/oxrna/dat/
├── backbone.pdb
├── all_atom.pdb
└── fit_report.json
```

Backbone atoms are P, OP1, OP2, O5′, C5′, C4′, O4′, C3′, O3′, C2′, O2′, C1′.
Full-atom output contains every template heavy atom: 22/A, 20/U, 23/G, 20/C.
Every successful archived example has **102 residues, 1,224 backbone atoms, and
2,181 full-atom atoms**. Standard PDB names, chain/residue identifiers, 0.001 Å
precision, TER/END records, and finite non-degenerate coordinates are independently
validated with Biopython. No hydrogens or virtual beads are written.

Inspect PDBs in a molecular viewer and `fit_report.json` for per-residue identities,
constraints, profiles, fitting residuals, and model-specific diagnostics.
`docs/validation_report.json` records all models, files, outputs/counts/checksums,
or blockers; `docs/coverage.tsv` is the final coverage table.

```bash
python -m pytest -q
```

Tests execute the installed CLI outside the source directory on every supported
model's actual archived files. Applicable PDB/CIF pairs and byte-identical MMCIF
aliases are compared. Both outputs are checked for sequence/order, identifiers,
atom sets, coordinates, preserved atoms, virtual constraints, and orientation
registration. Biopython's separate SVD fitter checks centered fits; the official
source oracle checks oxRNA. Inventory, CLI catalog, mappings, tests, README example
manifest, and generated outputs are compared for complete coverage.

Additional tests cover unchanged archives/fragments, native topology parsing,
known transformations of supplied fragments, degeneracy/reflection rejection,
aliases, custom fragments, FASTA, overwrite protection, and independent PDB
parsing. Generic recovery uses P/C4′/C3′, P/O5′/C4′, N9/C4/C8, and N1/C2/C6
triplets. Synthetic CG files do not replace authentic integration examples.

## Limits, provenance, and licensing

Independent placement does not enforce inter-residue bond geometry, remove
clashes, or establish force-field suitability. One template per base limits
conformational diversity. Known-atom/phosphate preservation can distort local
geometry. Modified residues, DNA, missing residues, hydrogens, trajectories,
different bead profiles, and alternative SPQR states are outside this release.
PDB field limits apply; long mmCIF chain IDs and oversized fields cause explicit
errors instead of truncation. Output-directory reuse replaces generated files;
input collisions, including symlink aliases, are rejected.

The supplied fragments originate from [pylelab/Arena](https://github.com/pylelab/Arena),
extracted from `IdealRNA/ACGU.rna.pdb`: A=A/1, C=A/2, G=B/7, U=B/8. Atom names
and coordinate bytes remain unchanged. The supplied [fragment README](single_nucleotide_fragments/README.md)
reports no Arena-wide license or redistribution notice covering these coordinates.
No fragment/dataset license is invented or inferred from unrelated source notices;
resolve redistribution permission with the authors before public redistribution.
Mapping provenance includes the cited primary definitions, supplied reference,
and hashed official SPQR/oxRNA helpers. Those external programs and their atomic
template libraries are not required for installation, reconstruction, or tests.
