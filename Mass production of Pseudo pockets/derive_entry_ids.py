#!/usr/bin/env python
"""
derive_entry_ids.py
===================
Regenerate rbp_entry_ids.json DETERMINISTICALLY from the already-committed
footprint_results.jsonl -- NOT from a live RCSB search.

Why: the original entry list came from a live RCSB search
(mass_production_RBPs.py `entries` stage). Re-running that search today could
return a different / larger set than the 4,599 the dataset was actually built
on, silently changing the dataset. The committed footprint_results.jsonl is the
ground truth of exactly which entries were processed, so we derive the id list
from it. This keeps the whole pipeline reproducible from committed files alone.

USAGE
-----
  python derive_entry_ids.py \
      --results footprint_results.jsonl \
      --out rbp_entry_ids.json
"""
import json, argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("USAGE")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="footprint_results.jsonl")
    ap.add_argument("--out", default="rbp_entry_ids.json")
    args = ap.parse_args()

    ids = sorted(set(json.loads(l)["pdb_id"] for l in open(args.results) if l.strip()))
    json.dump({"description": "Frozen RBP entry ids, derived deterministically "
                              "from footprint_results.jsonl (the entries actually "
                              "processed), not from a live RCSB search.",
               "source": args.results, "n_ids": len(ids), "pdb_ids": ids},
              open(args.out, "w"), indent=2)
    print(f"DONE: {len(ids)} entry ids -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
