"""Resumable, parallel driver that applies the pilot footprint engine
(rbp_pseudobind.process_complex, 6.0 A cutoff, 1-4 nt fragments, Jaccard 0.7,
symmetry dedup) to every frozen RBP entry. Writes one JSONL line per entry so
the run can resume after interruption without re-processing done entries."""
import os, sys, json, time
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.getcwd())
import rbp_pseudobind as pb

ids = json.load(open("rbp_entry_ids.json"))["pdb_ids"]
OUT = "footprint_results.jsonl"

done = set()
if os.path.exists(OUT):
    for ln in open(OUT):
        try: done.add(json.loads(ln)["pdb_id"])
        except Exception: pass
todo = [i for i in ids if i not in done]
print(f"total={len(ids)} done={len(done)} todo={len(todo)}", flush=True)

# Locked mass-production parameters (see Mass_production_RBPs_methods.md):
#   cutoff 6.0 A, RNA fragments 1-4 nt, symmetry dedup ON,
#   redundancy = continuous pharmacophore-Tanimoto, threshold 0.98
CUTOFF = 6.0
REDUNDANCY = "tanimoto"
REDUNDANCY_MAX = 0.98

def work(pid):
    # Retry with exponential backoff: RCSB throttles bulk coordinate downloads,
    # so a transient connection failure is retried rather than recorded as a
    # real (structural) failure.
    last = ""
    for attempt in range(5):
        try:
            rows, npk, frozen = pb.process_complex(
                pdb_id=pid, cutoff=CUTOFF,
                redundancy=REDUNDANCY, redundancy_max=REDUNDANCY_MAX)
            return {"pdb_id": pid, "status": "ok", "n_pockets": npk, "rows": rows}
        except Exception as e:
            last = str(e)[:200]
            if "Connection" in last or "Max retries" in last or "timed out" in last or "HTTPSConnectionPool" in last:
                time.sleep(2 ** attempt)   # 1,2,4,8,16 s
                continue
            break   # a real structural failure, don't retry
    return {"pdb_id": pid, "status": "fail", "error": last, "rows": []}

t0 = time.time(); n = 0
with open(OUT, "a") as fh, ThreadPoolExecutor(max_workers=6) as ex:
    futs = {ex.submit(work, pid): pid for pid in todo}
    for fut in as_completed(futs):
        res = fut.result()
        fh.write(json.dumps(res) + "\n"); fh.flush()
        n += 1
        if n % 200 == 0:
            print(f"{n}/{len(todo)}  {time.time()-t0:.0f}s  last={res['pdb_id']}:{res['status']}", flush=True)
print(f"DONE {n} entries in {time.time()-t0:.0f}s", flush=True)
