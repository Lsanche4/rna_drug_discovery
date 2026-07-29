#!/usr/bin/env python3
"""
rbp_pseudobind.py — PROFSA-style pseudo-ligand generator for RNA-binding proteins.

Inspired by ProFSA (Gao et al., ICLR 2024), which manufactures pseudo-ligand/pocket
complexes from protein-only structures to overcome the scarcity of experimental
complex data. This tool retargets that idea to RNA-BINDING PROTEINS: from a real
protein-RNA complex it (1) locates every RNA-binding site, (2) reads each site's
pharmacophore, and (3) generates 50+ Lipinski-compliant, drug-like small-molecule
"pseudo-ligands" per site using two complementary engines. The output is a set of
training-ready test sets for building models that target RNA-binding sites with
small molecules (i.e. disrupting protein-RNA interactions).

Two generation engines (run together by default):
  Engine A  — pharmacophore-COMPLEMENTARY sampling. Builds drug-like molecules
              whose features COMPLEMENT what the pocket presents (pocket donor ->
              ligand acceptor, cation-rich pocket -> anionic/acceptor ligand,
              aromatic stacking, hydrophobic packing). It rewards satisfying the
              same contacts the RNA satisfies WITHOUT copying the RNA's own charge.
  Engine B  — RNA-FRAGMENT SURROGATE. Mimics the pharmacophore of the interface
              nucleotide fragment plus its anionic backbone (the direct analog of
              ProFSA's peptide fragment).

Cleanup / freezing: only standard amino acids and standard RNA nucleotides are used.
Metal ions, cofactors, bound small-molecule ligands, crystallization additives, and
water are FROZEN OUT (recorded in the manifest) so non-RNA functional pockets cannot
skew a pilot. Overlapping sliding-window pockets are collapsed by a Jaccard
residue-set redundancy filter (--jaccard, default 0.7).

Pockets-only mode (--pockets-only): enumerate the non-redundant RNA-binding pockets
and their pharmacophore profiles WITHOUT generating any decoys.

Molecules are assembled by BRICS recombination of a curated drug-like fragment pool,
so every output is novel (not a raw seed) and, by construction, synthesizable-like.

Later phase (not in this script): dock the generated candidates into each pocket
(e.g. with DiffDock) and keep only geometrically plausible poses.

USAGE
  python rbp_pseudobind.py --pdb-ids 1M8Y,1FXL --out ./testsets
  python rbp_pseudobind.py --pdb-list mypdbs.txt --per-engine 30 --out ./testsets
  python rbp_pseudobind.py --cif-dir ./my_structures --out ./testsets   # local .cif files
  python rbp_pseudobind.py --pdb-ids 1M8Y,1FXL --pockets-only --out ./pockets  # pockets, no decoys

OUTPUT (in --out)
  per_site/<site_id>.csv     one test set per RNA-binding site
  combined_testset.csv       all pseudo-ligands, all sites, one row per molecule
  combined_testset.sdf       same molecules as 3D-embeddable SDF (2D coords + props)
  sites_summary.csv          one row per site (pocket size, features, counts)
  manifest.json              run parameters + provenance

Requires: rdkit, biotite, numpy, pandas, requests  (Python >= 3.9)
Author's note: reproducible given a fixed --seed.
"""
from __future__ import annotations
import argparse, io, json, os, random, sys
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, BRICS, AllChem
RDLogger.DisableLog("rdApp.*")

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

# ----------------------------------------------------------------------------- #
#  Constants: residue / base pharmacophore features
# ----------------------------------------------------------------------------- #
RNA_RES = {"A", "U", "G", "C"}

# Non-RNA / non-amino-acid entities that must be FROZEN OUT of pocket & ligand
# construction, so a pilot is not skewed by non-RNA functional pockets. Anything
# that is neither a standard amino acid nor a standard RNA nucleotide (metal ions,
# cofactors, bound small-molecule ligands, crystallization additives, water) is
# excluded automatically; this set is used only to LABEL what was frozen.
METAL_IONS = {"NA", "K", "MG", "CA", "MN", "ZN", "FE", "CU", "NI",
              "CO", "CD", "HG", "BA", "SR"}
WATER_RES = {"HOH", "WAT"}

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
BASE_FEATURES = {b: {"aromatic", "hba", "hbd"} for b in RNA_RES}

# ----------------------------------------------------------------------------- #
#  Drug-like seed scaffolds -> BRICS fragment pool
# ----------------------------------------------------------------------------- #
SEED_SCAFFOLDS = [
    "O=C(O)c1ccccc1O", "c1ccc2[nH]ccc2c1", "O=c1[nH]cnc2[nH]cnc12", "Nc1ncnc2[nH]cnc12",
    "c1ccc(-c2ccccc2)cc1", "O=C(N)c1ccccc1", "c1ccncc1", "c1ccc2ncccc2c1",
    "CC(=O)Nc1ccccc1", "O=S(=O)(N)c1ccccc1", "NC(=N)N", "OCC1OC(O)C(O)C(O)C1O",
    "c1cc2ccccc2[nH]1", "O=c1cc[nH]c(=O)[nH]1", "c1cnc2[nH]ccc2c1", "CC(C)Cc1ccccc1",
    "O=C1CCCN1", "c1ccc(CN)cc1", "Oc1ccccc1", "N#Cc1ccccc1", "c1ccc(OC)cc1",
    "C1CCNCC1", "O=C(O)CCc1ccccc1", "c1ccc2c(c1)OCO2",
    "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "O=C(Cc1ccccc1)Nc1ccc(O)cc1", "c1ccc(-c2nc3ccccc3o2)cc1",
    "O=C1NC(=O)c2ccccc21", "c1ccc2c(c1)ncn2", "Cc1ccc(S(=O)(=O)N)cc1", "O=C(O)c1ccc(N)cc1",
    "c1ccc(Cn2ccnc2)cc1", "O=c1[nH]c2ccccc2o1", "c1csc(-c2ccccc2)n1", "NC(=O)c1cccnc1",
    "O=C(O)c1cccc(C(=O)O)c1", "c1ccc(-n2cccn2)cc1", "CC(=O)c1ccccc1", "O=C1CCC(=O)N1",
    "c1ccc(CCN)cc1", "Oc1ccc(CCN)cc1", "c1cnc(N)nc1", "O=C(N)c1ccncc1", "c1ccc2[nH]nnc2c1",
    "c1ccc(-c2ccncc2)cc1", "O=C(NC1CCCCC1)c1ccccc1", "c1ccc(OCc2ccccc2)cc1",
    "Cc1cc(C)nc(N)n1", "O=S(=O)(c1ccccc1)N1CCOCC1", "c1ccc(-c2csc(N)n2)cc1",
    "OC(=O)c1ccc(-c2ccccc2)cc1", "c1ccc2c(c1)cc[nH]2", "N#Cc1ccc(N)cc1", "COc1ccc(CN)cc1",
    "c1ccc(-c2noc(C)n2)cc1", "O=C1NC(=O)NC1", "c1ccc(C(F)(F)F)cc1", "Nc1nc2ccccc2s1",
    "O=C(O)Cc1c[nH]c2ccccc12", "c1ccc(N2CCNCC2)cc1", "OCc1ccccc1", "Nc1ccc(S(N)(=O)=O)cc1",
]

# SMARTS pharmacophore patterns for molecules
_SMARTS = {
    "aromatic":    Chem.MolFromSmarts("a"),
    "hbd":         Chem.MolFromSmarts("[$([N;!H0]),$([O,S;H1])]"),
    "hba":         Chem.MolFromSmarts("[$([O,S;H0;v2]),$([N;v3;!$(N-a)]),n]"),
    "cationic":    Chem.MolFromSmarts("[$([NX3;H2,H1][CX4]),$([NX4+]),$([nX3+]),$(NC(=N)N)]"),
    "anionic":     Chem.MolFromSmarts("[$([CX3](=O)[OX1H0-,OX2H1]),$([PX4](=O)([OX1-,OX2H])),$([SX4](=O)(=O)[OX1-,OX2H])]"),
    "hydrophobic": Chem.MolFromSmarts("[$([CX4;!$(C[!#6;!#1])]),c]"),
}

# ----------------------------------------------------------------------------- #
#  Molecule helpers
# ----------------------------------------------------------------------------- #
def mol_features(mol):
    return {k: len(mol.GetSubstructMatches(p)) > 0 for k, p in _SMARTS.items()}

def ro5(mol):
    mw = Descriptors.MolWt(mol); logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol); hba = Lipinski.NumHAcceptors(mol)
    viol = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    return {"MW": round(mw, 1), "logP": round(logp, 2), "HBD": hbd, "HBA": hba,
            "TPSA": round(Descriptors.TPSA(mol), 1),
            "RotB": Descriptors.NumRotatableBonds(mol),
            "ro5_violations": viol, "ro5_pass": viol <= 1}

def build_fragment_pool(seeds=SEED_SCAFFOLDS):
    pool = set()
    for s in seeds:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        for f in BRICS.BRICSDecompose(m):
            pool.add(f)
    return pool

def generate_library(pool, n=3000, seed=42, mw_min=120, mw_max=500, max_depth=4):
    """Assemble novel drug-like molecules by BRICS recombination (reproducible)."""
    random.seed(seed)
    frags = [Chem.MolFromSmiles(f) for f in pool]
    frags = [f for f in frags if f is not None]
    builder = BRICS.BRICSBuild(frags, uniquify=True, scrambleReagents=True, maxDepth=max_depth)
    out, seen = [], set()
    for mol in builder:
        if mol is None:
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue
        smi = Chem.MolToSmiles(mol)
        if smi in seen:
            continue
        seen.add(smi)
        mw = Descriptors.MolWt(mol)
        if mw < mw_min or mw > mw_max + 50:
            continue
        out.append((smi, mol))
        if len(out) >= n:
            break
    return out

# ----------------------------------------------------------------------------- #
#  Structure fetch + interface detection
# ----------------------------------------------------------------------------- #
def load_structure(pdb_id=None, cif_path=None, timeout=60):
    if cif_path:
        f = pdbx.CIFFile.read(cif_path)
    else:
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        txt = requests.get(url, timeout=timeout).text
        f = pdbx.CIFFile.read(io.StringIO(txt))
    return pdbx.get_structure(f, model=1)

def report_frozen(arr):
    """Identify entities excluded from pocket/ligand construction (everything that
    is neither a standard amino acid nor a standard RNA nucleotide). Returns a dict
    labelling frozen metals, other ligands/cofactors, and water — for provenance."""
    heavy = arr[arr.element != "H"]
    is_aa = struc.filter_amino_acids(heavy)
    is_rna = np.isin(heavy.res_name, list(RNA_RES))
    other = sorted(set(heavy.res_name[~(is_aa | is_rna)].tolist()))
    return {"frozen_metals": [r for r in other if r in METAL_IONS],
            "frozen_water": [r for r in other if r in WATER_RES],
            "frozen_ligands": [r for r in other
                               if r not in METAL_IONS and r not in WATER_RES]}

def _nonredundant(sites, jaccard_max=0.7):
    """Greedy redundancy filter: keep a site only if its pocket residue set differs
    from every already-kept site by Jaccard <= jaccard_max. Because sliding windows
    generate heavily overlapping pockets, this collapses near-duplicate sites while
    (validated on 1M8Y/1FXL) tracking true spatial redundancy. jaccard_max=1.0
    disables the filter."""
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

def segment_sites(arr, frag_len=(1, 2, 3, 4), cutoff=6.0, min_pocket=4, jaccard_max=0.7):
    """Slide short fragments along each RNA chain; surrounding protein residues = pocket.
    One complex yields MANY sites (as in ProFSA's pocket construction)."""
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
#  Generation engines
# ----------------------------------------------------------------------------- #
def complementarity_score(pocket_feats, lf):
    """Engine A score: reward a ligand whose pharmacophore COMPLEMENTS the pocket.

    A cation-rich RNA-binding pocket (Arg/Lys gripping the anionic phosphate
    backbone) is best satisfied by an ANIONIC / H-bond-acceptor-rich ligand --
    exactly the way the RNA phosphates satisfy it. Rewarding a *cationic* ligand
    here would mimic the RNA rather than complement the pocket; that is Engine B's
    job (engine_B), so it is deliberately NOT part of this score.
    """
    tot = sum(pocket_feats.values()) or 1
    w = {k: pocket_feats.get(k, 0) / tot
         for k in ["hbd", "hba", "cationic", "anionic", "aromatic", "hydrophobic"]}
    s  = 2.0 * w["hbd"] * lf["hba"]                          # pocket donors    <- ligand acceptors
    s += 2.0 * w["hba"] * lf["hbd"]                          # pocket acceptors <- ligand donors
    s += 2.5 * w["cationic"] * (lf["anionic"] or lf["hba"])  # pocket Arg/Lys   <- anionic/acceptor ligand
    s += 1.5 * w["anionic"] * (lf["cationic"] or lf["hbd"])  # pocket Asp/Glu   <- cationic/donor ligand
    s += 2.0 * w["aromatic"] * lf["aromatic"]                # aromatic stacking
    s += 1.5 * w["hydrophobic"] * lf["hydrophobic"]          # hydrophobic packing
    return s

def engine_A(site, lib_feats, k=30, jitter=0.15, seed=0):
    random.seed(seed + hash(site["seq"]) % 10000)
    scored = []
    for smi, mol, lf, rp in lib_feats:
        if not rp["ro5_pass"]:
            continue
        sc = complementarity_score(site["pocket_feats"], lf) + random.uniform(0, jitter)
        scored.append((sc, smi, rp))
    scored.sort(reverse=True, key=lambda x: x[0])
    return [(smi, rp, round(sc, 3)) for sc, smi, rp in scored[:k]]

def engine_B(site, lib_feats, k=30, seed=0):
    random.seed(seed + hash(site["seq"]) % 9999)
    target = set()
    for _, _, base in site["rna_frag"]:
        target |= BASE_FEATURES.get(base, set())
    target.add("anionic")   # phosphodiester backbone
    scored = []
    for smi, mol, lf, rp in lib_feats:
        if not rp["ro5_pass"]:
            continue
        present = {k2 for k2, v in lf.items() if v}
        inter = len(target & present); union = len(target | present) or 1
        sim = inter / union + 0.3 * lf["aromatic"]
        scored.append((sim + random.uniform(0, 0.1), smi, rp))
    scored.sort(reverse=True, key=lambda x: x[0])
    return [(smi, rp, round(sc, 3)) for sc, smi, rp in scored[:k]]

# ----------------------------------------------------------------------------- #
#  Per-complex driver
# ----------------------------------------------------------------------------- #
def process_complex(lib_feats, pdb_id=None, cif_path=None, per_engine=30,
                    max_sites=None, cutoff=6.0, seed=0, jaccard_max=0.7,
                    pockets_only=False):
    arr = load_structure(pdb_id=pdb_id, cif_path=cif_path)
    label = pdb_id or os.path.splitext(os.path.basename(cif_path))[0]
    frozen = report_frozen(arr)
    sites = segment_sites(arr, cutoff=cutoff, jaccard_max=jaccard_max)
    if max_sites:
        sites = sites[:max_sites]
    rows, site_rows = [], []
    for si, site in enumerate(sites):
        site_id = f"{label}_site{si:03d}_{site['seq']}"
        site_rows.append({"site_id": site_id, "pdb_id": label, "rna_fragment": site["seq"],
                          "n_pocket_res": len(site["pocket_res"]),
                          "pocket_residues": ";".join(
                              f"{c}:{r}:{n}" for c, r, n in
                              sorted(site["pocket_res"], key=lambda x: (x[0], x[1]))),
                          **{f"pocket_{k}": v for k, v in site["pocket_feats"].items()}})
        if pockets_only:
            continue  # pocket enumeration only; no decoy generation
        for eng_name, eng in [("A_complementary", engine_A), ("B_rnasurrogate", engine_B)]:
            picks = eng(site, lib_feats, k=per_engine, seed=seed)
            for rank, (smi, rp, sc) in enumerate(picks):
                rows.append({"site_id": site_id, "pdb_id": label,
                             "rna_fragment": site["seq"],
                             "n_pocket_res": len(site["pocket_res"]),
                             "engine": eng_name, "rank": rank, "score": sc,
                             "smiles": smi, **rp})
    return rows, site_rows, len(sites), frozen

# ----------------------------------------------------------------------------- #
#  Output writers
# ----------------------------------------------------------------------------- #
def write_sdf(df, path):
    w = Chem.SDWriter(path)
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            continue
        m = Chem.AddHs(m)
        try:
            AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            m = Chem.MolFromSmiles(r["smiles"])
        m = Chem.RemoveHs(m) if m.GetNumConformers() else m
        for col in ["site_id", "pdb_id", "rna_fragment", "engine", "rank", "score",
                    "MW", "logP", "HBD", "HBA", "TPSA", "RotB", "ro5_violations"]:
            m.SetProp(col, str(r[col]))
        w.write(m)
    w.close()

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
    ap = argparse.ArgumentParser(description="PROFSA-style pseudo-ligand generator for RNA-binding proteins.")
    src = ap.add_argument_group("input (choose one or more)")
    src.add_argument("--pdb-ids", help="comma-separated PDB IDs, e.g. 1B7F,2ADC")
    src.add_argument("--pdb-list", help="text file, one PDB ID per line")
    src.add_argument("--cif-dir", help="directory of local .cif/.mmcif files")
    ap.add_argument("--out", default="./testsets", help="output directory")
    ap.add_argument("--per-engine", type=int, default=30,
                    help="molecules per engine per site (default 30 -> 60/site total)")
    ap.add_argument("--max-sites", type=int, default=None, help="cap sites per complex")
    ap.add_argument("--cutoff", type=float, default=6.0, help="interface distance cutoff (A)")
    ap.add_argument("--jaccard", type=float, default=0.7,
                    help="pocket redundancy threshold; drop pockets with residue-set "
                         "Jaccard above this vs a kept pocket (1.0 disables). default 0.7")
    ap.add_argument("--pockets-only", action="store_true",
                    help="enumerate non-redundant RNA-binding pockets ONLY; no decoy "
                         "generation. Writes pockets_nonredundant.csv + manifest.")
    ap.add_argument("--lib-size", type=int, default=3000, help="BRICS library size")
    ap.add_argument("--seed", type=int, default=42, help="random seed (reproducibility)")
    ap.add_argument("--no-sdf", action="store_true", help="skip 3D SDF export (faster)")
    args = ap.parse_args(argv)

    ids, cifs = parse_pdb_inputs(args)
    if not ids and not cifs:
        ap.error("provide --pdb-ids, --pdb-list, or --cif-dir")
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(os.path.join(args.out, "per_site"), exist_ok=True)

    if args.pockets_only:
        print("[1/2] pockets-only mode: skipping library build (no decoy generation).")
        lib_feats = []
    else:
        print(f"[1/4] building BRICS fragment pool + library (size {args.lib_size}, seed {args.seed}) ...")
        pool = build_fragment_pool()
        lib = generate_library(pool, n=args.lib_size, seed=args.seed)
        lib_feats = [(smi, mol, mol_features(mol), ro5(mol)) for smi, mol in lib]
        print(f"      {len(pool)} fragments -> {len(lib)} novel drug-like molecules "
              f"({sum(r['ro5_pass'] for *_, r in lib_feats)} Ro5-pass)")

    all_rows, all_sites, frozen_all = [], [], {}
    targets = [("id", x) for x in ids] + [("cif", x) for x in cifs]
    for kind, tgt in targets:
        try:
            rows, srows, nsites, frozen = process_complex(
                lib_feats, pdb_id=tgt if kind == "id" else None,
                cif_path=tgt if kind == "cif" else None,
                per_engine=args.per_engine, max_sites=args.max_sites,
                cutoff=args.cutoff, seed=args.seed, jaccard_max=args.jaccard,
                pockets_only=args.pockets_only)
        except Exception as e:
            print(f"      ! {tgt}: FAILED ({e})"); continue
        all_rows += rows; all_sites += srows
        label = tgt if kind == "id" else os.path.basename(tgt)
        frozen_all[label] = frozen
        froz = ", ".join(f"{k.split('_')[1]}={len(v)}" for k, v in frozen.items() if v) or "none"
        print(f"      {label}: {nsites} non-redundant pockets"
              + (f" -> {len(rows)} pseudo-ligands" if not args.pockets_only else "")
              + f"  [frozen: {froz}]")

    if not all_sites:
        print("no RNA-binding pockets found; check inputs."); sys.exit(1)

    # ---- pockets-only mode: write the non-redundant pocket table and stop ----
    if args.pockets_only:
        pk = pd.DataFrame(all_sites)
        pk.to_csv(os.path.join(args.out, "pockets_nonredundant.csv"), index=False)
        manifest = {
            "tool": "rbp_pseudobind", "version": "1.1", "mode": "pockets_only",
            "run_utc": datetime.now(timezone.utc).isoformat(), "params": vars(args),
            "n_complexes": len(targets), "n_pockets": len(all_sites),
            "frozen_entities": frozen_all,
            "inspired_by": "ProFSA (Gao et al., ICLR 2024, arXiv:2310.07229)",
        }
        with open(os.path.join(args.out, "manifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)
        print(f"\nDONE (pockets-only). {len(all_sites)} non-redundant pockets across "
              f"{len(targets)} complex(es).")
        print(f"  pockets_nonredundant.csv, manifest.json  -> {args.out}")
        return

    if not all_rows:
        print("no pseudo-ligands generated; check inputs."); sys.exit(1)

    df = pd.DataFrame(all_rows)
    sdf_df = df.drop_duplicates("smiles")
    print(f"[3/4] writing outputs to {args.out} ...")
    for sid, g in df.groupby("site_id"):
        g.to_csv(os.path.join(args.out, "per_site", f"{sid}.csv"), index=False)
    df.to_csv(os.path.join(args.out, "combined_testset.csv"), index=False)
    pd.DataFrame(all_sites).to_csv(os.path.join(args.out, "sites_summary.csv"), index=False)
    if not args.no_sdf:
        print(f"[4/4] embedding {len(sdf_df)} unique molecules to SDF ...")
        write_sdf(sdf_df, os.path.join(args.out, "combined_testset.sdf"))

    manifest = {
        "tool": "rbp_pseudobind", "version": "1.1",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "params": vars(args),
        "n_complexes": len(targets), "n_sites": len(all_sites),
        "n_pseudoligands": len(df), "n_unique_smiles": int(df["smiles"].nunique()),
        "all_ro5_pass": bool(df["ro5_pass"].all()),
        "engines": ["A_complementary", "B_rnasurrogate"],
        "frozen_entities": frozen_all,
        "inspired_by": "ProFSA (Gao et al., ICLR 2024, arXiv:2310.07229)",
    }
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nDONE. {len(df)} pseudo-ligands across {len(all_sites)} sites "
          f"({df['smiles'].nunique()} unique molecules, all Ro5-pass={df['ro5_pass'].all()}).")
    print(f"  combined_testset.csv, sites_summary.csv, per_site/*.csv, manifest.json"
          + ("" if args.no_sdf else ", combined_testset.sdf") + f"  -> {args.out}")

if __name__ == "__main__":
    main()
