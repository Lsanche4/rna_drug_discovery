# rbp-pseudobind — druggable-pocket extractor for RNA-binding proteins (ProFSA-style)

`rbp_pseudobind.py` finds the **pockets on an RNA-binding protein** — the patches of
the protein surface that grip its RNA. **The protein is the target.** The goal is to
locate the protein pockets that a small molecule could later occupy to block the
protein–RNA interaction. Each RNA interface fragment is used only as a **probe that
marks where a functional pocket sits on the protein** — the pseudo-ligand, exactly
as ProFSA uses its peptide fragment.

## Target vs. probe (read this first)

- **Target = the protein.** Every pocket the tool outputs is a set of *protein*
  residues — the site on the protein where the RNA binds, i.e. the site you would
  drug. This is what the tool extracts, characterizes, and hands downstream.
- **Probe = the RNA fragment.** The short stretch of RNA is *not* the target. It is
  a marker: wherever the protein contacts RNA, there is a functional pocket worth
  drugging, and the fragment shows us where. We keep the fragment only to locate and
  define the pocket.

We are targeting the protein that binds the RNA — not the RNA.

## Why (the ProFSA idea, applied to RNA-binding proteins)

Experimentally determined pocket–ligand complexes are scarce, which limits
large-scale pocket pretraining. **ProFSA** (Gao et al., ICLR 2024,
arXiv:2310.07229) gets around this by mining abundant protein-only structures:
it cuts a short **peptide fragment** out of a protein, treats that fragment as a
"pseudo-ligand" probe, and defines the surrounding protein residues as the pocket —
manufacturing millions of (pocket, fragment) pairs without needing a bound drug.

`rbp-pseudobind` applies the same idea to **RNA-binding proteins**, with the RNA
interface fragment playing the role of ProFSA's peptide fragment. If anything the
mapping is cleaner here: the probe (RNA) and the target (protein) are genuinely
*different* molecules, so each fragment marks a real, functional binding interface
on the protein we want to drug. From one protein–RNA complex we obtain many protein
pockets — the raw material for characterizing, comparing, and later pretraining on
these druggable sites.

**This is a pocket-extraction tool.** It does *not* generate drug-like small
molecules — there is no molecule library, no BRICS recombination, and no
Lipinski / Rule-of-Five drug-likeness filter. Those belong to a later phase; here we
only extract the protein pockets and the RNA probe that marks each one, matching
ProFSA directly.

## What a "pocket" is here

For each RNA-contact site the tool records:

- **the pocket (the target)** — every *protein* residue with at least one heavy atom
  within **6 Å** of the RNA fragment (the 6 Å heavy-atom cutoff follows
  ProFSA / Uni-Mol). This residue set is the druggable site on the protein;
- **a pharmacophore profile** of that pocket — counts of H-bond donors/acceptors,
  cationic (Arg/Lys), anionic (Asp/Glu), aromatic, and hydrophobic residues, so
  pockets can be compared and clustered;
- **the RNA fragment (the probe)** — the contiguous run of 1–4 nucleotides that
  marks the pocket, recorded so the site can be relocated in the structure.

Short fragments are slid along each RNA chain to enumerate every contact site, so
one complex yields many protein pockets.

### Cleanup and de-duplication

- **Non-RNA entities are frozen out.** Only standard amino acids and standard RNA
  nucleotides are used. Metal ions, cofactors, bound small molecules,
  crystallization additives, and water are excluded (and recorded in the manifest)
  so non-RNA functional pockets can't skew the set.
- **Redundant pockets are collapsed.** Overlapping sliding-window pockets are merged
  by a Jaccard residue-set filter (`--jaccard`, default 0.7).
- **Symmetry copies are collapsed (default on).** A crystal often contains several
  copies of the same protein–RNA complex (e.g. chains A and B), so an identical
  binding site is enumerated once per copy. By default the tool keeps one pocket per
  unique site (same RNA fragment at the same position), retaining the best-resolved
  copy. Pass `--keep-symmetry-copies` to keep every copy. *Example: 1M8Y has two
  protein copies, so its 42 raw pockets collapse to **21 unique sites**; 1FXL has a
  single copy and stays at **12**.*

## The two example complexes

The tool is validated on two well-characterized single-stranded-RNA-binding proteins.
Both are X-ray structures with the RNA resolved in the binding groove:

| PDB | Protein | RNA-binding module | RNA in the structure | What it is |
|-----|---------|--------------------|----------------------|------------|
| **1M8Y** | **Pumilio-1** (human PUM1) | PUF / Pumilio-homology domain | `5'-AUUGUACAUA-3'` (Nanos response element, NRE) | The textbook sequence-specific ssRNA reader: eight PUF repeats each recognize one base. Gives large, well-defined pockets. |
| **1FXL** | **HuD** (ELAV-like neuronal antigen ELAVL4) | two tandem RRM domains | `5'-UUUUAUUUU-3'` (AU-rich element from c-fos mRNA) | Canonical RRM recognition of an AU-rich element; the classic model for AU-rich-element regulation. |

Between them they cover the two dominant ssRNA-recognition folds (PUF repeats and
RRMs), which is why they make a good pilot: different pocket shapes, different base
preferences, both with unambiguous interfaces.

## Install

```bash
conda create -n rbp -c conda-forge python=3.11 biotite numpy pandas requests
conda activate rbp
```

## Run

Extraction is the whole tool — just point it at one or more complexes:

```bash
# Extract protein pockets from the two example complexes
python rbp_pseudobind.py --pdb-ids 1M8Y,1FXL --out ./pockets

# Your own structures (downloaded from RCSB by ID)
python rbp_pseudobind.py --pdb-list mypdbs.txt --out ./pockets

# Local .cif / .mmcif files (no network needed)
python rbp_pseudobind.py --cif-dir ./my_structures --out ./pockets
```

Structures are fetched from RCSB by PDB ID, so ID-based runs need internet access;
`--cif-dir` runs entirely offline on local files.

### Useful options

| Flag | Meaning | Default |
|------|---------|---------|
| `--cutoff` | interface distance cutoff in Å (pocket = protein residues with a heavy atom within this of the RNA fragment) | 6.0 |
| `--jaccard` | pocket redundancy threshold; drop pockets whose residue set overlaps a kept pocket above this Jaccard (1.0 disables) | 0.7 |
| `--max-sites` | cap the number of pockets per complex | none |
| `--out` | output directory | `./pockets` |
| `--seed` | random seed (reproducibility) | 42 |

## Output

Written to `--out`:

- **`pockets_nonredundant.csv`** — one row per protein pocket:
  `pocket_id, pdb_id, rna_fragment, frag_start, n_pocket_res, pocket_residues`
  (the pocket's protein residues as `chain:resid:name`, semicolon-separated), plus
  the pharmacophore feature counts (`pocket_hbd`, `pocket_hba`, `pocket_cationic`, …).
- **`manifest.json`** — run parameters, the ProFSA reference, and the list of
  entities frozen out of each structure, for provenance.

Each row is a self-contained **protein pocket**: the listed protein residues are the
druggable site, and the RNA fragment tells you where on the structure it is. Pull the
atoms (by PDB ID + the listed residues) to render, cluster, or feed the pockets into
downstream analysis or docking.

## Relationship to ProFSA — what matches and what differs

| | ProFSA (the paper) | rbp-pseudobind |
|---|---|---|
| Target (what the pocket is on) | the protein | **the protein** ✅ same |
| Probe (pseudo-ligand marking the pocket) | peptide fragment, used directly | **RNA interface fragment, used directly** ✅ same idea |
| Pocket definition | protein residues with a heavy atom within 6 Å of the fragment | 6 Å heavy-atom cutoff ✅ same |
| Fragment enumeration | many fragments slid along the chain → many pockets | slide 1–4-nt RNA fragments → many pockets ✅ same |
| Fragment length | 1–8 residues | 1–4 nucleotides (RNA adaptation) |
| Drug-molecule generation | none | **none** ✅ |

ProFSA additionally pairs this data with a contrastive-pretraining objective and
applies peptide-specific steps (long-range residue exclusion, N-/C-terminal capping)
that don't apply to a separate RNA chain. `rbp-pseudobind` implements the
**pocket / fragment extraction** half; model training is out of scope for this script.

## Reference

Gao, Jia, Mo, Ni, Ma, Ma, Lan. *ProFSA: Self-Supervised Pocket Pretraining via
Protein Fragment-Surroundings Alignment.* ICLR 2024. arXiv:2310.07229.
