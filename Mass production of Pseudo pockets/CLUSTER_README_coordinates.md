# Building the coordinate dataset on the cluster (thorshammer)

**Why this exists.** The shared `rbp_pseudo_pockets.csv` stores each pocket as a
list of **residues** (`chain:resid:resname`), not 3D atom coordinates. A
DrugCLIP-style encoder needs per-atom elements + XYZ for both the protein pocket
and the RNA fragment. Those coordinates live only in the source structures, so we
re-extract them from RCSB **on the cluster** — this keeps the multi-GB coordinate
file next to the GPUs and off github (only code is shared).

The re-extraction uses the **same validated engine** (`rbp_pseudobind.py`,
Tanimoto 0.98) that built the CSV, so the pocket set is identical — it just now
carries coordinates. Pilot check: 1FXL → 19 pockets, 1M8Y → 20 pockets (same as
the CSV).

## Files
| file | role |
|---|---|
| `rbp_pseudobind.py` | extraction engine (now with `emit_coords=True`) |
| `rbp_entry_ids.json` | frozen list of 4,599 PDB ids (the dataset definition) |
| `build_drugclip_lmdb.py` | builds the coordinate LMDB in DrugCLIP schema |
| `fetch_uniprot_map.py` | pdb_id → UniProt, for a leak-free protein-level split |

## Steps (on the cluster, in the repo folder)
```bash
# 0. deps (extraction needs no GPU)
pip install lmdb msgpack biotite numpy requests

# 1. protein->UniProt map for the split (fast, ~a few min)
python fetch_uniprot_map.py --ids rbp_entry_ids.json --out pdb2uniprot.json

# 2. build the coordinate LMDB (network-bound; resumable via progress jsonl)
python build_drugclip_lmdb.py \
    --ids rbp_entry_ids.json \
    --out rbp_drugclip.lmdb \
    --workers 6 --redundancy tanimoto --redundancy-max 0.98
# re-run the same command if interrupted; it skips entries already in
# drugclip_build_progress.jsonl
```

## LMDB record schema (one value per pocket)
```
pocket             str    pocket_id "1FXL_pk00"  (grouping / split key)
pdb_id             str    source PDB entry
pocket_atoms       [str]  protein pocket heavy-atom elements
pocket_coordinates [[f]]  (N,3) XYZ, angstrom
atoms              [str]  RNA-fragment (pseudo-ligand) heavy-atom elements
coordinates        [[f]]  (M,3) XYZ, angstrom
smi                str    RNA fragment sequence, e.g. "AUUU" (ligand id)
label              int    1  (all rows are true pairs; DrugCLIP makes negatives
                              by in-batch shuffling — no negative rows stored)
```

## Split (answering the "split key" decision)
Group pockets by **UniProt** (from `pdb2uniprot.json`) so all pockets of one
protein stay on one side — splitting on `pdb_id` alone leaks the same protein
across train/test. Entries with no protein UniProt (empty list) fall back to
`pdb_id` grouping.

Fully ProFSA-faithful alternative: ProFSA clusters at 70% sequence identity;
RCSB serves those cluster memberships via the same API (see the note in
`fetch_uniprot_map.py`) if you want cluster-level grouping instead.

## The three decisions the training agent asked about
1. **RNA-fragment-as-pseudo-ligand** — yes; this is the ProFSA design, and the
   only framing under which the 222,255 pockets support contrastive training.
2. **Split key** — UniProt (protein-level), not `pdb_id`; see above.
3. **Re-extract coordinates** — yes, via `build_drugclip_lmdb.py`; the CSV has no
   coordinates and nothing was cached, so the structures must be re-read. Doing
   it on the cluster avoids moving a multi-GB file through github.
