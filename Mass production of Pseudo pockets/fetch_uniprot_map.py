#!/usr/bin/env python
"""
fetch_uniprot_map.py
====================
Map every RBP entry to its protein UniProt accession(s), for a LEAK-FREE split
of the DrugCLIP dataset. Splitting the pocket set by pdb_id alone lets the SAME
protein (crystallized as several PDB entries) land in both train and test, which
inflates measured performance. Grouping by UniProt (or, better, by sequence
cluster) puts all pockets of one protein on one side of the split.

Uses the RCSB Data API GraphQL endpoint, batched, one protein-entity lookup per
entry. RNA/DNA entities are ignored (they have no UniProt). An entry with no
protein UniProt (e.g. engineered/synthetic construct) is recorded with an empty
list and should be grouped by pdb_id as a fallback.

OUTPUT: pdb2uniprot.json  ->  {pdb_id: [uniprot_acc, ...], ...}

USAGE
-----
  python fetch_uniprot_map.py --ids rbp_entry_ids.json --out pdb2uniprot.json

NOTE on the fully ProFSA-faithful split: ProFSA clusters PDB entries at 70%
sequence identity. RCSB also serves those cluster memberships
(rcsb_cluster_membership, identity 70) via the same API; swap the query below to
pull cluster ids if you want cluster-level rather than UniProt-level grouping.
"""
import sys, json, time, argparse
import requests

GRAPHQL = "https://data.rcsb.org/graphql"
Q = """query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    polymer_entities{
      entity_poly{rcsb_entity_polymer_type}
      rcsb_polymer_entity_container_identifiers{uniprot_ids}
    }
  }
}"""


def fetch_batch(ids):
    for attempt in range(4):
        try:
            d = requests.post(GRAPHQL, json={"query": Q, "variables": {"ids": ids}}, timeout=120).json()
            if d.get("errors"):
                raise RuntimeError(str(d["errors"])[:200])
            return d["data"]["entries"]
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("USAGE")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", default="rbp_entry_ids.json")
    ap.add_argument("--out", default="pdb2uniprot.json")
    args = ap.parse_args()

    ids = json.load(open(args.ids))["pdb_ids"]
    mapping, n_empty = {}, 0
    for i in range(0, len(ids), 50):
        for e in fetch_batch(ids[i:i + 50]):
            accs = []
            for pe in (e.get("polymer_entities") or []):
                if (pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein":
                    accs += (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []
            accs = sorted(set(accs))
            mapping[e["rcsb_id"]] = accs
            if not accs:
                n_empty += 1
        if i % 500 == 0:
            print(f"  {len(mapping)}/{len(ids)}", flush=True)
        time.sleep(0.15)
    json.dump(mapping, open(args.out, "w"), indent=2)
    n_uni = len(set(a for v in mapping.values() for a in v))
    print(f"DONE: {len(mapping)} entries -> {n_uni} unique UniProt accessions "
          f"({n_empty} entries with no protein UniProt) -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
