# rbp-pseudobind

**PROFSA-style pseudo-ligand generator for RNA-binding proteins.**

From real protein–RNA complexes, this tool manufactures training-ready test sets of
drug-like, Lipinski-compliant small molecules that target the protein's RNA-binding
sites — a data-augmentation engine for building models that drug RNA-binding proteins
(i.e. disrupt protein–RNA interactions with small molecules).

## Why

Experimentally-determined pocket–ligand complexes are scarce, which limits models
that need interaction data. ProFSA (Gao et al., ICLR 2024, arXiv:2310.07229) solved
this for generic pockets by **synthesizing** pseudo-ligand/pocket complexes from
protein-only structures. `rbp-pseudobind` adapts that idea to **RNA-binding
proteins**: it treats each protein↔RNA interface as a druggable pocket and generates
50+ small-molecule pseudo-ligands per site, so a scarce set of real complexes becomes
a large labelled dataset.

## How it works

1. **Fetch** protein–RNA complexes (RCSB by PDB ID, or local `.cif`).
2. **Freeze non-RNA entities** — only standard amino acids and standard RNA
   nucleotides (A/U/G/C) are used. Metal ions, cofactors, bound small-molecule
   ligands, crystallization additives, and water are excluded from pocket and
   ligand construction (recorded in `manifest.json → frozen_entities`) so a
   protein's non-RNA functional pockets cannot skew the dataset.
3. **Detect interfaces** — slide 1–4-nt fragments along each RNA chain; the
   surrounding protein residues (heavy atom within a cutoff, default 6 Å) are the
   pocket. One complex yields many sites. Overlapping sliding-window pockets are
   collapsed by a **Jaccard residue-set redundancy filter** (`--jaccard`, default
   0.7): a pocket is dropped if >70% of its residues match an already-kept pocket.
4. **Fingerprint** each site's pharmacophore (H-bond donors/acceptors, cationic,
   anionic, aromatic, hydrophobic residues; the RNA phosphate backbone is anionic).
5. **Generate** small molecules by BRICS recombination of a curated drug-like
   fragment pool (every molecule is novel), then select per site with two engines:
   - **Engine A — pocket-complementary**: molecules whose features *complement* the
     pocket (pocket donor→ligand acceptor, cation-rich pocket→anionic/acceptor
     ligand, aromatic stacking, hydrophobic packing). It satisfies the same contacts
     the RNA does **without** copying the RNA's own negative charge.
   - **Engine B — RNA-fragment surrogate**: molecules that mimic the pharmacophore of
     the interface nucleotide fragment plus its anionic backbone (the direct analog
     of ProFSA's peptide fragment).
6. **Filter** to Lipinski Ro5-compliant molecules; **package** per-site and combined
   test sets (CSV + SDF) with full provenance.

**Pockets-only mode** (`--pockets-only`): run steps 1–4 to enumerate the
non-redundant RNA-binding pockets and their pharmacophore profiles, writing
`pockets_nonredundant.csv`, **without** generating any decoys.

The two engines produce chemically distinct decoy sets (see `engine_chemspace.png`),
which enriches the augmented dataset.

## Install

```bash
conda create -n rbp -c conda-forge python=3.11 rdkit biotite freesasa numpy pandas requests
conda activate rbp
# or: pip install -r requirements.txt   (rdkit via pip wheels)
```

## Usage

```bash
# From PDB IDs (downloads from RCSB)
python rbp_pseudobind.py --pdb-ids 1M8Y,1FXL --out ./testsets

# From a list file (one PDB ID per line)
python rbp_pseudobind.py --pdb-list rbp_pdbs.txt --out ./testsets

# From local structures
python rbp_pseudobind.py --cif-dir ./my_structures --out ./testsets

# Pockets only — enumerate non-redundant RNA-binding pockets, no decoy generation
python rbp_pseudobind.py --pdb-ids 1M8Y,1FXL --pockets-only --out ./pockets

# Tune: 50 molecules per engine (100/site), cap sites per complex, wider interface
python rbp_pseudobind.py --pdb-ids 1M8Y --per-engine 50 --max-sites 10 --cutoff 6.5 --out ./testsets
```

Key options: `--per-engine N` (molecules per engine per site; N×2 total),
`--max-sites`, `--cutoff` (Å), `--jaccard` (pocket redundancy threshold; 1.0
disables), `--pockets-only` (pockets, no decoys), `--lib-size`, `--seed`
(reproducible), `--no-sdf` (faster).

## Output (in `--out`)

| file | contents |
|---|---|
| `combined_testset.csv` | all pseudo-ligands, one row per molecule (SMILES + Ro5 properties + site provenance) |
| `combined_testset.sdf` | unique molecules as 3D-embedded SDF with properties |
| `per_site/<site_id>.csv` | one test set per RNA-binding site |
| `sites_summary.csv` | one row per site (pocket size, residue list, pharmacophore features) |
| `pockets_nonredundant.csv` | *(pockets-only mode)* one row per non-redundant pocket, with residue list + pharmacophore features |
| `manifest.json` | run parameters + provenance, including `frozen_entities` per structure |

Each row carries `site_id, pdb_id, rna_fragment, engine, rank, score, smiles`, plus
`MW, logP, HBD, HBA, TPSA, RotB, ro5_violations, ro5_pass`.

## Next phase (docking validation)

The generated candidates are pharmacophore-matched, not energy-placed. To add
geometric fidelity, dock each site's molecules into the pocket (e.g. DiffDock) and
keep only confidently-scored poses. This is a documented later phase, not in the
current script.

## Caveats

- Pseudo-ligands are **decoys for augmentation**, not validated binders. As in ProFSA,
  their value is teaching a model interaction patterns at scale; downstream validation
  (docking, assays) is required before any experimental claim.
- Pharmacophore features are assigned per residue/base at neutral pH; explicit
  protonation/tautomer handling is not modelled.
- Engine B surrogates approximate a nucleotide's base pharmacophore with drug-like
  chemistry; they intentionally drift from the (large, charged) real nucleotide.

Inspired by ProFSA (Gao, Jia, Mo, Ni, Ma, Ma, Lan; ICLR 2024).
