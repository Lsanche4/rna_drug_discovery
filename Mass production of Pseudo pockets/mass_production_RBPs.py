#!/usr/bin/env python
"""
mass_production_RBPs.py
=======================
Mass-produce RNA-binding-protein (RBP) pseudo-pockets from every PDB structure
that contains a protein in contact with bound RNA. This is the scale-up of the
validated two-protein pilot (1M8Y, 1FXL) to the whole database, built for a
downstream protein-encoder / contrastive-learning dataset.

REPRODUCIBILITY IS THE POINT. Every entry ID, every metadata field, and every
pocket comes from a public, allowlisted API or a deterministic computation on
downloaded coordinates -- nothing is typed by hand. Each stage writes its inputs
and parameters to disk so a coworker can audit or rerun any step independently.
Full method rationale (incl. the Tanimoto-vs-Jaccard-vs-ProFSA discussion) is in
Mass_production_RBPs_methods.md.

PIPELINE (three subcommands, run in order)
------------------------------------------
  entries    Query the RCSB Search API v2 for all protein+RNA entries with
             <=10 protein chains (excludes ribosomes/spliceosomes). Freeze the
             ID list + the verbatim query body + a UTC timestamp to
             rbp_entry_ids.json.
  metadata   Fetch per-entry metadata (method, resolution, organism, chain
             counts, title, release date) from the RCSB Data API (GraphQL,
             batched). Build rbp_master_table.csv, one row per entry, each row
             carrying the literal Data-API URL that reproduces it.
  footprint  Apply the pilot pocket-extraction engine (rbp_pseudobind.py) to
             every entry: slide 1-4 nt RNA fragments, take protein residues
             within 6.0 A as the pocket, collapse redundant pockets by
             pharmacophore-Tanimoto (threshold 0.98). Resumable: writes one
             JSONL line per entry and skips entries already done. Emits
             rbp_pseudo_pockets.csv and rbp_footprint_failures.csv.

USAGE
-----
  python mass_production_RBPs.py entries   --out .
  python mass_production_RBPs.py metadata  --out .
  python mass_production_RBPs.py footprint --out . --workers 12

PARAMETERS (all documented, all overridable)
--------------------------------------------
  --max-protein-chains 10     megacomplex exclusion cutoff
  --cutoff 6.0                interface distance (A) for the pocket definition
  --redundancy tanimoto       pocket redundancy metric (tanimoto|jaccard|none)
  --redundancy-max 0.98       collapse pockets more similar than this
  --workers 12                parallel download/parse threads (footprint stage)

Requires: requests, pandas, numpy, biotite, and rbp_pseudobind.py in the same
directory (the pilot engine, imported unchanged).
"""
import os, sys, json, time, argparse
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import pandas as pd

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
ENTRY_REST = "https://data.rcsb.org/rest/v1/core/entry/{}"

IDS_FILE = "rbp_entry_ids.json"
META_RAW = "meta_raw_entries.json"
MASTER_CSV = "rbp_master_table.csv"
POCKETS_CSV = "rbp_pseudo_pockets.csv"
FAILS_CSV = "rbp_footprint_failures.csv"
JSONL = "footprint_results.jsonl"


# --------------------------------------------------------------------------- #
#  Stage 1: entries  (RCSB Search API v2)
# --------------------------------------------------------------------------- #
def _term(attr, op, val):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attr, "operator": op, "value": val}}


def build_search_body(max_protein_chains, start=0, rows=10000):
    """The verbatim Search-API query body. Frozen into rbp_entry_ids.json so the
    entry set is reproducible even as the PDB grows."""
    nodes = [
        _term("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
        _term("rcsb_entry_info.polymer_entity_count_RNA", "greater_or_equal", 1),
        _term("rcsb_entry_info.polymer_entity_count_protein", "less_or_equal", max_protein_chains),
    ]
    return {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry",
            "request_options": {"paginate": {"start": start, "rows": rows},
                                "results_content_type": ["experimental"]}}


def stage_entries(args):
    body = build_search_body(args.max_protein_chains, start=0, rows=10000)
    all_ids, start = [], 0
    while True:
        b = build_search_body(args.max_protein_chains, start=start, rows=10000)
        r = requests.post(SEARCH_URL, json=b, timeout=120)
        if r.status_code == 204:      # no (more) results
            break
        r.raise_for_status()
        data = r.json()
        ids = [x["identifier"] for x in data["result_set"]]
        all_ids += ids
        if not ids or len(all_ids) >= data["total_count"]:
            break
        start += 10000
        time.sleep(0.3)
    all_ids = sorted(set(all_ids))
    out = {
        "description": "RBP mass-production entry list: PDB entries with >=1 protein "
                       "entity, >=1 RNA entity, and <=%d protein entities "
                       "(megacomplex exclusion)." % args.max_protein_chains,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "search_api": SEARCH_URL,
        "query_body": body,
        "n_ids": len(all_ids),
        "pdb_ids": all_ids,
    }
    path = os.path.join(args.out, IDS_FILE)
    json.dump(out, open(path, "w"), indent=2)
    print(f"entries: {len(all_ids)} frozen -> {path}")
    return all_ids


# --------------------------------------------------------------------------- #
#  Stage 2: metadata  (RCSB Data API, GraphQL batched)
# --------------------------------------------------------------------------- #
GQL = """query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    struct{title}
    exptl{method}
    rcsb_entry_info{resolution_combined
      polymer_entity_count_protein polymer_entity_count_RNA polymer_entity_count_DNA}
    rcsb_accession_info{initial_release_date}
    polymer_entities{entity_poly{rcsb_entity_polymer_type}
      rcsb_entity_source_organism{ncbi_scientific_name}}
  }
}"""


def _fetch_batch(ids):
    for attempt in range(3):
        try:
            r = requests.post(GRAPHQL_URL, json={"query": GQL, "variables": {"ids": ids}}, timeout=120)
            d = r.json()
            if d.get("errors"):
                raise RuntimeError(str(d["errors"])[:200])
            return d["data"]["entries"]
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2)


def _flatten(e):
    info = e["rcsb_entry_info"]
    res = info.get("resolution_combined")
    res = res[0] if res else None
    methods = "; ".join(sorted(set(m["method"] for m in (e.get("exptl") or []))))
    orgs = []
    for pe in (e.get("polymer_entities") or []):
        for o in (pe.get("rcsb_entity_source_organism") or []):
            if o.get("ncbi_scientific_name"):
                orgs.append(o["ncbi_scientific_name"])
    pid = e["rcsb_id"]
    return {"pdb_id": pid, "experimental_method": methods, "resolution_A": res,
            "n_protein_chains": info.get("polymer_entity_count_protein"),
            "n_RNA_chains": info.get("polymer_entity_count_RNA"),
            "n_DNA_chains": info.get("polymer_entity_count_DNA"),
            "title": (e.get("struct") or {}).get("title", ""),
            "organism": "; ".join(sorted(set(orgs))) if orgs else "",
            "release_date": (e.get("rcsb_accession_info") or {}).get("initial_release_date", ""),
            "data_api_url": ENTRY_REST.format(pid)}


def stage_metadata(args):
    ids = json.load(open(os.path.join(args.out, IDS_FILE)))["pdb_ids"]
    entries = []
    for i in range(0, len(ids), 50):
        entries += _fetch_batch(ids[i:i + 50])
        if i % 500 == 0:
            print(f"  metadata {len(entries)}/{len(ids)}")
        time.sleep(0.15)
    json.dump(entries, open(os.path.join(args.out, META_RAW), "w"))
    df = pd.DataFrame([_flatten(e) for e in entries]).sort_values("pdb_id").reset_index(drop=True)
    path = os.path.join(args.out, MASTER_CSV)
    df.to_csv(path, index=False)
    print(f"metadata: {len(df)} rows -> {path}")
    return df


# --------------------------------------------------------------------------- #
#  Stage 3: footprint  (pilot engine, parallel + resumable)
# --------------------------------------------------------------------------- #
def stage_footprint(args):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import rbp_pseudobind as pb

    ids = json.load(open(os.path.join(args.out, IDS_FILE)))["pdb_ids"]
    jsonl = os.path.join(args.out, JSONL)
    done = set()
    if os.path.exists(jsonl):
        for ln in open(jsonl):
            try:
                done.add(json.loads(ln)["pdb_id"])
            except Exception:
                pass
    todo = [i for i in ids if i not in done]
    print(f"footprint: total={len(ids)} done={len(done)} todo={len(todo)}")

    def work(pid):
        # Retry with exponential backoff on transient network errors: RCSB
        # throttles bulk coordinate downloads, so a connection failure is
        # retried rather than mis-recorded as a real structural failure.
        last = ""
        for attempt in range(5):
            try:
                rows, npk, frozen = pb.process_complex(
                    pdb_id=pid, cutoff=args.cutoff,
                    redundancy=args.redundancy, redundancy_max=args.redundancy_max)
                return {"pdb_id": pid, "status": "ok", "n_pockets": npk, "rows": rows}
            except Exception as e:
                last = str(e)[:200]
                if any(s in last for s in ("Connection", "Max retries", "timed out", "HTTPSConnectionPool")):
                    time.sleep(2 ** attempt)   # 1,2,4,8,16 s
                    continue
                break   # a real structural failure -> record it, don't retry
        return {"pdb_id": pid, "status": "fail", "error": last, "rows": []}

    t0, n = time.time(), 0
    with open(jsonl, "a") as fh, ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(work, pid): pid for pid in todo}
        for fut in as_completed(futs):
            res = fut.result()
            fh.write(json.dumps(res) + "\n")
            fh.flush()
            n += 1
            if n % 200 == 0:
                print(f"  {n}/{len(todo)}  {time.time()-t0:.0f}s")
    consolidate(args)


def consolidate(args):
    """Fold the JSONL log into the two final CSVs + a manifest."""
    jsonl = os.path.join(args.out, JSONL)
    pocket_rows, fails = [], []
    for ln in open(jsonl):
        rec = json.loads(ln)
        if rec["status"] == "ok":
            pocket_rows += rec["rows"]
        else:
            fails.append({"pdb_id": rec["pdb_id"], "error": rec.get("error", "")})
    pk = pd.DataFrame(pocket_rows)
    pk.to_csv(os.path.join(args.out, POCKETS_CSV), index=False)
    pd.DataFrame(fails).to_csv(os.path.join(args.out, FAILS_CSV), index=False)
    manifest = {"tool": "mass_production_RBPs", "run_utc": datetime.now(timezone.utc).isoformat(),
                "params": {k: getattr(args, k) for k in
                           ("cutoff", "redundancy", "redundancy_max", "max_protein_chains")},
                "n_pockets": len(pk), "n_failures": len(fails),
                "n_entries_with_pockets": pk["pdb_id"].nunique() if len(pk) else 0}
    json.dump(manifest, open(os.path.join(args.out, "manifest.json"), "w"), indent=2)
    print(f"consolidate: {len(pk)} pockets across "
          f"{manifest['n_entries_with_pockets']} entries; {len(fails)} failures")


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("PIPELINE")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["entries", "metadata", "footprint", "consolidate"])
    ap.add_argument("--out", default=".", help="output directory")
    ap.add_argument("--max-protein-chains", type=int, default=10,
                    help="megacomplex exclusion: drop entries with more protein entities. default 10")
    ap.add_argument("--cutoff", type=float, default=6.0,
                    help="interface distance cutoff (A). default 6.0")
    ap.add_argument("--redundancy", choices=["tanimoto", "jaccard", "none"], default="tanimoto",
                    help="pocket redundancy metric. default tanimoto")
    ap.add_argument("--redundancy-max", type=float, default=0.98,
                    help="collapse pockets more similar than this. default 0.98")
    ap.add_argument("--workers", type=int, default=6,
                    help="footprint parallel threads. default 6 (RCSB throttles higher)")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)

    {"entries": stage_entries, "metadata": stage_metadata,
     "footprint": stage_footprint, "consolidate": consolidate}[args.stage](args)


if __name__ == "__main__":
    main()
