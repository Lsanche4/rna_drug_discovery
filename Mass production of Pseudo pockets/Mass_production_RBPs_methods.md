# Mass Production of RBP Pseudo-Pockets — Methods & Provenance

**Purpose.** Scale the validated two-protein pilot (1M8Y, 1FXL) to every
RNA-binding protein structure in the PDB that carries bound RNA, producing a
curated pseudo-pocket dataset for downstream protein-encoder / contrastive
learning (the DrugCLIP-style setup discussed for this project).

Everything here is reproducible from public, allowlisted APIs. No value in the
output tables is typed by hand; every number comes from an API call or a
deterministic computation on downloaded coordinates. This document records the
exact endpoints, query bodies, and parameters so a coworker can regenerate the
dataset from scratch and audit each step.

---

## 1. What counts as an "RBP" here

The pilot method defines a pocket **from the RNA it grips**: the pocket is the
set of protein residues whose heavy atoms fall within 6.0 Å of a bound RNA
fragment. That definition *requires RNA to be present in the structure*. So the
reproducible, method-consistent definition of the RBP set is:

> **Every PDB entry containing at least one protein chain AND at least one RNA
> chain, excluding megacomplexes (> 10 protein entities).**

The megacomplex cutoff removes ribosomes and spliceosomes — assemblies with
tens of protein chains where a single RNA-contact footprint is not a meaningful
"pocket." The cutoff is a documented parameter (`polymer_entity_count_protein
<= 10`) and easy to change.

**Honest scope caveat (state this to anyone who uses the dataset):** this is the
set of proteins seen *in contact with RNA in a solved structure*. It is NOT the
full universe of RNA-binding proteins — many true RBPs were crystallized in
their apo (RNA-free) form, and those cannot be footprinted by this method (no
RNA to define the pocket). Capturing apo RBPs would need a different pocket
detector (e.g. fpocket on the bare protein) and a functional RBP definition
(Pfam RRM/KH/dsRBD families or GO "RNA binding"); that is a separate pipeline.

---

## 2. Entry list — RCSB Search API v2

**Endpoint:** `POST https://search.rcsb.org/rcsbsearch/v2/query`

**Verbatim query body** (also frozen inside `rbp_entry_ids.json`):

```json
{
  "query": {
    "type": "group", "logical_operator": "and",
    "nodes": [
      {"type":"terminal","service":"text","parameters":{
        "attribute":"rcsb_entry_info.polymer_entity_count_protein",
        "operator":"greater_or_equal","value":1}},
      {"type":"terminal","service":"text","parameters":{
        "attribute":"rcsb_entry_info.polymer_entity_count_RNA",
        "operator":"greater_or_equal","value":1}},
      {"type":"terminal","service":"text","parameters":{
        "attribute":"rcsb_entry_info.polymer_entity_count_protein",
        "operator":"less_or_equal","value":10}}
    ]
  },
  "return_type": "entry",
  "request_options": {
    "paginate": {"start": 0, "rows": 10000},
    "results_content_type": ["experimental"]
  }
}
```

**Result (frozen `RETRIEVED_DATE` in `rbp_entry_ids.json`):**
`TOTAL_ENTRIES` entries. Method breakdown: X-ray `N_XRAY`, cryo-EM `N_EM`,
NMR `N_NMR` (+ a handful of hybrid/fiber/neutron). The ID list is frozen so the
dataset is reproducible even as the PDB grows.

---

## 3. Per-entry metadata — RCSB Data API (GraphQL)

**Endpoint:** `POST https://data.rcsb.org/graphql` (batched 50 IDs/request).
Per-entry REST equivalent (one entry, browser-auditable):
`https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}` — this URL is stored in the
`data_api_url` column of `rbp_master_table.csv` for every row, so any row can be
verified by pasting its URL into a browser.

Fields pulled: experimental method, resolution, structure title, source
organism(s), protein/RNA/DNA entity counts, initial release date.

---

## 4. Pocket extraction — the pilot engine, applied unchanged in structure

The extraction engine is `rbp_pseudobind.py`, the *same* code validated on the
pilot. Parameters for mass production:

| Parameter | Value | Notes |
|---|---|---|
| Interface cutoff | **6.0 Å** | protein residue kept if any heavy atom ≤ 6.0 Å of an RNA fragment atom |
| RNA fragment lengths | **1–4 nt** | sliding window along each RNA chain; each fragment is a pseudo-ligand probe |
| Min pocket size | **4 residues** | pockets smaller than this are dropped |
| Symmetry-copy dedup | **ON** | one pocket per unique (RNA-fragment-sequence, start) across crystal copies |
| Redundancy metric | **pharmacophore Tanimoto** | see §5 |
| Redundancy threshold | **0.98** | pockets more similar than this are collapsed |

Non-standard entities (metals, ligands, cofactors, water) are frozen out of
pocket construction automatically and logged for provenance.

---

## 5. Redundancy filter: Tanimoto vs. the pilot's Jaccard, and vs. ProFSA

This is the one deliberate methods change from the pilot, made at the
researcher's request ("Tanimoto not Jaccard"). It matters enough to spell out.

**What the pilot did.** The pilot collapsed near-duplicate pockets with a
**Jaccard index on the pocket residue set** (|A∩B|/|A∪B|) at threshold 0.7.

**Why a rename alone does nothing.** On a *binary set* (residue present / absent)
the Tanimoto coefficient and the Jaccard index are the *same number*. To make
Tanimoto a genuinely different metric it must be computed over a *count /
continuous vector*. We use each pocket's **pharmacophore feature-count vector** —
counts of (cationic, anionic, H-bond donor, H-bond acceptor, aromatic,
hydrophobic) residues — with the cheminformatics Tanimoto form:

> T(a, b) = (a·b) / (|a|² + |b|² − a·b)

So the filter now collapses pockets with near-identical *binding chemistry*
rather than overlapping *residue identity*. Threshold **0.98** was chosen from a
threshold sweep on the pilot proteins so the retained-site counts stay close to
the validated pilot (1M8Y ≈ 20 vs 21 Jaccard; 1FXL ≈ 19 vs 12 Jaccard) while
using the requested Tanimoto metric. The counts differ because the two metrics
collapse *different* groups of pockets — there is no threshold that reproduces
Jaccard exactly, and that is expected.

**What ProFSA actually does (checked against the paper, arXiv:2310.07229).**
ProFSA does **not** use Tanimoto or Jaccard to deduplicate pockets. Its
redundancy control is:
  1. **Upstream, at the sequence level:** it starts from non-redundant PDB
     entries *clustered at 70% sequence identity* ("Our data construction
     process starts with non-redundant PDB entries clustered at a 70 percent
     sequence identity threshold").
  2. **Downstream, for the task:** pocket similarity is measured as **cosine
     similarity between learned pocket representations** (following DeeplyTough),
     not a set-overlap metric on residues.

So our per-pocket Tanimoto filter is *our own* redundancy step, not a copy of
ProFSA's. A future version that wants to be faithful to ProFSA's redundancy
strategy would add a 70%-sequence-identity clustering pass over the entry list
(RCSB provides sequence-cluster memberships via the Search API), and reserve
representation-space cosine similarity for the trained encoder — this is noted
as future work, not implemented here.

---

## 6. Outputs

| File | Contents |
|---|---|
| `rbp_entry_ids.json` | frozen entry ID list + verbatim Search query + timestamp |
| `rbp_master_table.csv` | one row per PDB entry: method, resolution, chain counts, title, organism, release date, `data_api_url` |
| `rbp_dataset_overview.png` | entries-per-method bar + resolution histogram |
| `rbp_pseudo_pockets.csv` | one row per pseudo-pocket: `pdb_id, pocket_id, rna_fragment, n_pocket_res, pocket_residues`, pharmacophore feature counts |
| `mass_production_RBPs.py` | the consolidated, resumable pipeline (query → fetch → footprint) |
| `rbp_footprint_failures.csv` | entries that failed extraction, with the error (nothing is silently dropped) |

---

## 7. Reproduce from scratch

```bash
# 1. freeze the entry list + fetch metadata + build master table
python mass_production_RBPs.py entries   --out .
python mass_production_RBPs.py metadata  --out .
# 2. extract pockets (resumable; re-run to continue after interruption)
python mass_production_RBPs.py footprint --out . --workers 12
```

All parameters above are defaults in the script and are echoed into a
`manifest.json` on every run for provenance.
