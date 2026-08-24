#!/usr/bin/env python
"""
build_drugclip_lmdb.py
======================
Build a DrugCLIP-style LMDB of RNA-binding pseudo-pockets, ON THE CLUSTER, so no
multi-GB coordinate file has to move through github. Run this AFTER cloning the
"Mass production Pseudo pocket" repo; it re-extracts coordinates directly from
RCSB (the residue-level rbp_pseudo_pockets.csv has no coordinates) and writes the
per-atom record every DrugCLIP pocket/ligand encoder expects.

WHY RE-EXTRACT: the shared CSV stores pocket RESIDUES (chain:resid:resname), not
3D atom coordinates. Coordinates only exist in the source structures. This script
pulls each structure once and emits atoms+XYZ using the SAME validated engine
(rbp_pseudobind.process_complex, Tanimoto 0.98) that built the CSV, so the pocket
set is identical -- only now it carries coordinates.

RECORD SCHEMA (one LMDB value per pocket, msgpack/pickle dict):
    pocket            str    pocket_id, e.g. "1FXL_pk00"   (grouping / split key)
    pdb_id            str    source PDB entry
    pocket_atoms      [str]  heavy-atom element symbols of the PROTEIN pocket
    pocket_coordinates [[f]] (N,3) XYZ of those atoms, A
    atoms             [str]  heavy-atom elements of the RNA fragment (pseudo-ligand)
    coordinates       [[[f]]] conformer list [(M,3)] -- DrugCLIP samples one
                             conformer, so this is a 1-element list holding the
                             single (M,3) fragment coordinate array, NOT a bare (M,3)
    smi               str    RNA fragment sequence as the ligand id (e.g. "AUUU")
    label             int    1 (all rows are true pocket-fragment pairs; DrugCLIP
                             draws negatives by in-batch shuffling, so no negative
                             rows are stored here)
Keys are the pocket_id encoded utf-8. A companion split file (see --split-col)
lists which keys go to train/valid/test.

USAGE (on the cluster)
----------------------
  python build_drugclip_lmdb.py \
      --ids rbp_entry_ids.json \
      --out rbp_drugclip.lmdb \
      --workers 6 \
      --split-col uniprot          # protein-level split; needs pdb2uniprot.json
                                    # (produced by fetch_uniprot_map.py)

Requires: rbp_pseudobind.py (same dir), lmdb, msgpack, requests, biotite, numpy.
"""
import os, sys, json, time, argparse, pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rbp_pseudobind as pb

try:
    import lmdb
except ImportError:
    sys.exit("pip install lmdb  (not found)")
try:
    import msgpack
    def _dumps(o): return msgpack.packb(o, use_bin_type=True)
except ImportError:
    def _dumps(o): return pickle.dumps(o)   # fallback; set --serializer pickle


def one(pid, cutoff, redundancy, redundancy_max):
    """Extract coordinate-bearing pockets for one entry, with retry/backoff on
    RCSB throttling (same policy as the CSV pipeline)."""
    last = ""
    for attempt in range(5):
        try:
            rows, npk, frozen = pb.process_complex(
                pdb_id=pid, cutoff=cutoff, redundancy=redundancy,
                redundancy_max=redundancy_max, emit_coords=True)
            recs = []
            for r in rows:
                recs.append({
                    "pocket": r["pocket_id"], "pdb_id": r["pdb_id"],
                    "pocket_atoms": r["pocket_atoms"],
                    # pocket side is NEVER conformer-sampled -> flat (N,3)
                    "pocket_coordinates": r["pocket_coordinates"],
                    "atoms": r["frag_atoms"],
                    # mol/ligand side: DrugCLIP AffinityDataset samples a conformer
                    # (len == n_conformers, then [idx] -> one (M,3) array), so the
                    # RNA-fragment coords MUST be wrapped as a single-conformer list.
                    "coordinates": [r["frag_coordinates"]],
                    "smi": r["rna_fragment"], "label": 1,
                })
            return pid, recs, None
        except Exception as e:
            last = str(e)[:200]
            if any(s in last for s in ("Connection", "Max retries", "timed out", "HTTPSConnectionPool")):
                time.sleep(2 ** attempt); continue
            break
    return pid, [], last


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("RECORD SCHEMA")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", default="rbp_entry_ids.json", help="frozen entry-id json")
    ap.add_argument("--out", default="rbp_drugclip.lmdb", help="output LMDB path")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--cutoff", type=float, default=6.0)
    ap.add_argument("--redundancy", choices=["tanimoto", "jaccard", "none"], default="tanimoto")
    ap.add_argument("--redundancy-max", type=float, default=0.98)
    ap.add_argument("--map-size-gb", type=float, default=50.0, help="LMDB max map size")
    ap.add_argument("--progress-jsonl", default="drugclip_build_progress.jsonl",
                    help="resumable per-entry log (which entries are already written)")
    args = ap.parse_args()

    ids = json.load(open(args.ids))["pdb_ids"]
    done = set()
    if os.path.exists(args.progress_jsonl):
        for ln in open(args.progress_jsonl):
            try: done.add(json.loads(ln)["pdb_id"])
            except Exception: pass
    todo = [i for i in ids if i not in done]
    print(f"entries total={len(ids)} done={len(done)} todo={len(todo)}", flush=True)

    env = lmdb.open(args.out, map_size=int(args.map_size_gb * (1 << 30)))
    n_pockets, n_fail, t0 = 0, 0, time.time()
    with open(args.progress_jsonl, "a") as plog, \
         ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(one, pid, args.cutoff, args.redundancy, args.redundancy_max): pid
                for pid in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            pid, recs, err = fut.result()
            if err:
                n_fail += 1
            else:
                with env.begin(write=True) as txn:
                    for r in recs:
                        txn.put(r["pocket"].encode("utf-8"), _dumps(r))
                n_pockets += len(recs)
            plog.write(json.dumps({"pdb_id": pid, "n_pockets": len(recs),
                                   "error": err}) + "\n"); plog.flush()
            if i % 200 == 0:
                print(f"  {i}/{len(todo)}  pockets={n_pockets}  fails={n_fail}  {time.time()-t0:.0f}s", flush=True)
    with env.begin() as txn:
        total_keys = txn.stat()["entries"]
    env.close()
    print(f"DONE: {total_keys} pocket records in {args.out}  "
          f"({n_fail} entry failures)  {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
