#!/usr/bin/env python
"""
rbp_pseudobind.py — druggable-pocket extractor for RNA-binding proteins (ProFSA-style).

TARGET = THE PROTEIN. This tool finds the pockets on an RNA-binding protein -- the
patches of protein surface that grip its RNA and that a small molecule could later
occupy to block the protein-RNA interaction. Every pocket it outputs is a set of
PROTEIN residues.

The RNA interface fragment is only a PROBE: it marks where a functional pocket sits
on the protein. It is the pseudo-ligand, used directly, exactly as ProFSA (Gao et
al., ICLR 2024, arXiv:2310.07229) uses its peptide fragment -- we keep it only to
locate and define the protein pocket. We target the protein that binds the RNA, not
the RNA.

This is a POCKET-EXTRACTION tool. It does not generate drug-like small molecules:
there is no molecule library, no BRICS recombination, and no Lipinski / Rule-of-Five
filter. The RNA fragment taken from the structure is the pseudo-ligand, matching
ProFSA directly.

From one protein-RNA complex we obtain many (pocket, RNA-fragment) pairs:
  1. slide short fragments (1-4 nt) along each RNA chain,
  2. define the pocket as the protein residues within 6 A of the fragment,
  3. profile each pocket's pharmacophore (donor/acceptor/charge/aromatic/hydrophobic),
  4. collapse near-duplicate pockets by a Jaccard residue-set filter.

USAGE
  python rbp_pseudobind.py --pdb-ids 1M8Y,1FXL --out ./pockets
  python rbp_pseudobind.py --pdb-list mypdbs.txt --out ./pockets
  python rbp_pseudobind.py --cif-dir ./my_structures --out ./pockets   # local .cif files

OUTPUT (in --out)
  pockets_nonredundant.csv   one row per RNA-binding pocket (fragment, residues, features)
  manifest.json              run parameters + provenance (incl. frozen-out entities)

Requires: biotite, numpy, pandas, requests  (Python >= 3.9)
Reproducible given a fixed --seed.
"""

import os
import sys
import io
import json
import argparse
import random
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

# ----------------------------------------------------------------------------- #
#  Constants: residue pharmacophore features
# ----------------------------------------------------------------------------- #
RNA_RES = {"A", "U", "G", "C"}

# Non-RNA / non-amino-acid entities that are FROZEN OUT of pocket construction, so a
# pilot is not skewed by non-RNA functional pockets. Anything that is neither a
# standard amino acid nor a standard RNA nucleotide (metal ions, cofactors, bound
# small-molecule ligands, crystallization additives, water) is excluded
# automatically; these sets are used only to LABEL what was frozen.
METAL_IONS = {"NA", "K", "MG", "CA", "MN", "ZN", "FE", "CU", "NI",
              "CO", "CD", "HG", "BA", "SR"}
WATER_RES = {"HOH", "WAT"}

# Pharmacophore features presented by each amino-acid side chain. Used to profile a
# pocket so pockets can be compared, clustered, or fed into downstream analysis.
AA_FEATURES = {
    "ARG": {"cationic", "hbd"}, "LYS": {"cationic", "hbd"},
    "HIS": {"cationic", "aromatic", "hbd", "hba"},
    "ASP": {"anionic", "hba"}, "GLU": {"anionic", "hba"},
    "SER": {"hbd", "hba"}, "THR": {"hbd", "hba"}, "TYR": {"aromatic", "hbd", "hba"},
    "ASN": {"hbd", "hba"}, "GLN": {"hbd", "hba"},
    "PHE": {"aromatic", "hydrophobic"}, "TRP": {"aromatic", "hbd", "hydrophobic"},
    "MET": {"hydrophobic"}, "LEU": {"hydrophobic"}, "ILE": {"hydrophobic"},
    "VAL": {"hydrophobic"}, "ALA": {"hydrophobic"}, "PRO": {"hydrophobic"},
    "CYS": {"hbd"}, "GLY": set(),
}

# ----------------------------------------------------------------------------- #
#  Structure fetch + entity freezing
# ----------------------------------------------------------------------------- #
def load_structure(pdb_id=None, cif_path=None, timeout=60):
    """Load a single-model structure from a local .cif or by download from RCSB."""
    if cif_path:
        f = pdbx.CIFFile.read(cif_path)
    else:
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        txt = requests.get(url, timeout=timeout).text
        f = pdbx.CIFFile.read(io.StringIO(txt))
    return pdbx.get_structure(f, model=1)


def report_frozen(arr):
    """Identify entities excluded from pocket construction (everything that is
    neither a standard amino acid nor a standard RNA nucleotide). Returns a dict
    labelling frozen metals, other ligands/cofactors, and water -- for provenance."""
    heavy = arr[arr.element != "H"]
    is_aa = struc.filter_amino_acids(heavy)
    is_rna = np.isin(heavy.res_name, list(RNA_RES))
    other = sorted(set(heavy.res_name[~(is_aa | is_rna)].tolist()))
    return {"frozen_metals": [r for r in other if r in METAL_IONS],
            "frozen_water": [r for r in other if r in WATER_RES],
            "frozen_ligands": [r for r in other
                               if r not in METAL_IONS and r not in WATER_RES]}

# ----------------------------------------------------------------------------- #
#  Pocket enumeration
# ----------------------------------------------------------------------------- #
def _nonredundant(sites, jaccard_max=0.7):
    """Greedy redundancy filter: keep a pocket only if its residue set differs from
    every already-kept pocket by Jaccard <= jaccard_max. Because sliding windows
    generate heavily overlapping pockets, this collapses near-duplicate sites.
    jaccard_max >= 1.0 disables the filter."""
    if jaccard_max >= 1.0:
        return sites
    kept = []
    for s in sorted(sites, key=lambda x: -len(x["pocket_res"])):
        ps = {(c, r) for c, r, _ in s["pocket_res"]}
        redundant = False
        for q in kept:
            qs = {(c, r) for c, r, _ in q["pocket_res"]}
            if len(ps & qs) / len(ps | qs) > jaccard_max:
                redundant = True
                break
        if not redundant:
            kept.append(s)
    return kept


def _dedup_symmetry(sites):
    """Collapse symmetry-copy pockets. A crystal asymmetric unit often contains
    several copies of the same protein-RNA complex (e.g. chains A and B), so the
    identical binding site is enumerated once per copy. Two pockets are the same
    binding site when they share the RNA fragment sequence AND its start position;
    we keep the best-resolved copy (largest pocket). Returns one pocket per unique
    binding site."""
    best = {}
    for s in sites:
        key = (s["seq"], s["frag_start"])
        if key not in best or len(s["pocket_res"]) > len(best[key]["pocket_res"]):
            best[key] = s
    # preserve original ordering by first appearance
    seen, out = set(), []
    for s in sites:
        key = (s["seq"], s["frag_start"])
        if key not in seen:
            seen.add(key)
            out.append(best[key])
    return out


def segment_sites(arr, frag_len=(1, 2, 3, 4), cutoff=6.0, min_pocket=4, jaccard_max=0.7):
    """Slide short RNA fragments along each RNA chain; the surrounding PROTEIN
    residues (heavy atom within `cutoff` A) are the pocket -- the druggable site on
    the protein (the target). One complex yields many pockets (as in ProFSA's pocket
    construction). The RNA fragment is only the probe/pseudo-ligand that marks each
    pocket."""
    heavy = arr[arr.element != "H"]
    rna = heavy[np.isin(heavy.res_name, list(RNA_RES))]
    prot = heavy[struc.filter_amino_acids(heavy)]
    if rna.array_length() == 0 or prot.array_length() == 0:
        return []
    rna_res = sorted(set(zip(rna.chain_id.tolist(), rna.res_id.tolist(), rna.res_name.tolist())),
                     key=lambda x: (x[0], x[1]))
    by_chain = defaultdict(list)
    for c, rid, rn in rna_res:
        by_chain[c].append((rid, rn))
    cell = struc.CellList(prot, cell_size=cutoff)
    sites = []
    for c, reslist in by_chain.items():
        for L in frag_len:
            for i in range(len(reslist) - L + 1):
                frag = reslist[i:i + L]
                if [r[0] for r in frag] != list(range(frag[0][0], frag[0][0] + L)):
                    continue  # require contiguous numbering
                frag_rids = [rid for rid, _ in frag]
                fatoms = rna[np.isin(rna.res_id, frag_rids) & (rna.chain_id == c)]
                pocket = set()
                for k in range(fatoms.array_length()):
                    idx = cell.get_atoms(fatoms.coord[k], radius=cutoff)
                    idx = idx[idx != -1]
                    for a in idx:
                        pocket.add((prot.chain_id[a], int(prot.res_id[a]), prot.res_name[a]))
                if len(pocket) >= min_pocket:
                    pf = defaultdict(int)
                    for _, _, rn in pocket:
                        for feat in AA_FEATURES.get(rn, set()):
                            pf[feat] += 1
                    sites.append({"chain": c, "seq": "".join(rn for _, rn in frag),
                                  "frag_start": frag[0][0],
                                  "rna_frag": [(c, rid, rn) for rid, rn in frag],
                                  "pocket_res": pocket, "pocket_feats": dict(pf)})
    # dedupe identical (seq, pocket-residue-set) ...
    uniq = {}
    for s in sites:
        key = (s["seq"], frozenset((c, r) for c, r, _ in s["pocket_res"]))
        uniq.setdefault(key, s)
    # ... then collapse near-duplicate pockets by Jaccard residue-set overlap
    return _nonredundant(list(uniq.values()), jaccard_max=jaccard_max)

# ----------------------------------------------------------------------------- #
#  Per-complex driver
# ----------------------------------------------------------------------------- #
def process_complex(pdb_id=None, cif_path=None, cutoff=6.0, jaccard_max=0.7,
                    max_sites=None, dedup_symmetry=True):
    """Extract non-redundant RNA-binding pockets from one complex. Returns
    (pocket_rows, n_pockets, frozen_entities). When dedup_symmetry is True (default),
    binding sites duplicated across symmetry copies (chains) are collapsed to one."""
    arr = load_structure(pdb_id=pdb_id, cif_path=cif_path)
    label = pdb_id or os.path.splitext(os.path.basename(cif_path))[0]
    frozen = report_frozen(arr)
    sites = segment_sites(arr, cutoff=cutoff, jaccard_max=jaccard_max)
    if dedup_symmetry:
        sites = _dedup_symmetry(sites)
    if max_sites:
        sites = sites[:max_sites]
    rows = []
    for si, site in enumerate(sites):
        pocket_id = f"{label}_pk{si:02d}"
        rows.append({"pocket_id": pocket_id, "pdb_id": label,
                     "rna_fragment": site["seq"], "frag_start": site["frag_start"],
                     "n_pocket_res": len(site["pocket_res"]),
                     "pocket_residues": ";".join(
                         f"{c}:{r}:{n}" for c, r, n in
                         sorted(site["pocket_res"], key=lambda x: (x[0], x[1]))),
                     **{f"pocket_{k}": v for k, v in site["pocket_feats"].items()}})
    return rows, len(sites), frozen

# ----------------------------------------------------------------------------- #
#  CLI
# ----------------------------------------------------------------------------- #
def parse_pdb_inputs(args):
    ids, cifs = [], []
    if args.pdb_ids:
        ids += [x.strip().upper() for x in args.pdb_ids.split(",") if x.strip()]
    if args.pdb_list:
        with open(args.pdb_list) as fh:
            ids += [ln.strip().upper() for ln in fh if ln.strip() and not ln.startswith("#")]
    if args.cif_dir:
        for fn in sorted(os.listdir(args.cif_dir)):
            if fn.lower().endswith((".cif", ".mmcif")):
                cifs.append(os.path.join(args.cif_dir, fn))
    return ids, cifs


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="RNA-binding-pocket extractor (ProFSA-style). The RNA interface "
                    "fragment is the pseudo-ligand; no molecule generation.")
    src = ap.add_argument_group("input (choose one or more)")
    src.add_argument("--pdb-ids", help="comma-separated PDB IDs, e.g. 1M8Y,1FXL")
    src.add_argument("--pdb-list", help="text file, one PDB ID per line")
    src.add_argument("--cif-dir", help="directory of local .cif/.mmcif files")
    ap.add_argument("--out", default="./pockets", help="output directory")
    ap.add_argument("--cutoff", type=float, default=6.0,
                    help="interface distance cutoff (A): pocket = protein residues "
                         "with a heavy atom within this of the RNA fragment. default 6.0")
    ap.add_argument("--jaccard", type=float, default=0.7,
                    help="pocket redundancy threshold; drop pockets with residue-set "
                         "Jaccard above this vs a kept pocket (1.0 disables). default 0.7")
    ap.add_argument("--max-sites", type=int, default=None, help="cap pockets per complex")
    ap.add_argument("--keep-symmetry-copies", action="store_true",
                    help="keep binding sites duplicated across crystal symmetry "
                         "copies (chains); by default these are collapsed to one "
                         "pocket per unique site")
    ap.add_argument("--seed", type=int, default=42, help="random seed (reproducibility)")
    args = ap.parse_args(argv)

    random.seed(args.seed)
    np.random.seed(args.seed)

    ids, cifs = parse_pdb_inputs(args)
    if not ids and not cifs:
        ap.error("provide --pdb-ids, --pdb-list, or --cif-dir")
    os.makedirs(args.out, exist_ok=True)

    all_rows, frozen_all = [], {}
    targets = [("id", x) for x in ids] + [("cif", x) for x in cifs]
    for kind, tgt in targets:
        try:
            rows, npk, frozen = process_complex(
                pdb_id=tgt if kind == "id" else None,
                cif_path=tgt if kind == "cif" else None,
                cutoff=args.cutoff, jaccard_max=args.jaccard, max_sites=args.max_sites,
                dedup_symmetry=not args.keep_symmetry_copies)
        except Exception as e:
            print(f"      ! {tgt}: FAILED ({e})"); continue
        all_rows += rows
        label = tgt if kind == "id" else os.path.basename(tgt)
        frozen_all[label] = frozen
        froz = ", ".join(f"{k.split('_')[1]}={len(v)}" for k, v in frozen.items() if v) or "none"
        print(f"      {label}: {npk} non-redundant RNA-binding pockets  [frozen: {froz}]")

    if not all_rows:
        print("no RNA-binding pockets found; check inputs."); sys.exit(1)

    pk = pd.DataFrame(all_rows)
    pk.to_csv(os.path.join(args.out, "pockets_nonredundant.csv"), index=False)
    manifest = {
        "tool": "rbp_pseudobind", "version": "2.0", "mode": "pocket_extraction",
        "run_utc": datetime.now(timezone.utc).isoformat(), "params": vars(args),
        "n_complexes": len(targets), "n_pockets": len(all_rows),
        "frozen_entities": frozen_all,
        "inspired_by": "ProFSA (Gao et al., ICLR 2024, arXiv:2310.07229)",
        "note": "The RNA interface fragment is the pseudo-ligand, used directly "
                "(no molecule generation), matching ProFSA.",
    }
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nDONE. {len(all_rows)} non-redundant RNA-binding pockets across "
          f"{len(targets)} complex(es).")
    print(f"  pockets_nonredundant.csv, manifest.json  -> {args.out}")


if __name__ == "__main__":
    main()
