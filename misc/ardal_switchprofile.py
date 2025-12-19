#!/usr/bin/env python3
"""
Ardal: machine-local backend switchpoint profiler.

Generates synthetic CSR matrices over a grid of column counts (m) and densities (p),
times bit-packed vs roaring backends for each hot operation, and writes a JSON
profile that Ardal can load at runtime to choose a backend per operation.

Examples
--------
# Typical sweep (single-thread) for pairwise + radii + topK:
ardal-switchprofile \
  --nrows 5000 \
  --m 5000 10000 25000 50000 100000 250000 \
  --p 0.0005 0.001 0.0025 0.005 0.01 0.02 0.05 0.1 \
  --threads 1 \
  --ops pairwise radius topk \
  --r 25 50 100 \
  --k 10 50

# Repeat at your common thread count (e.g. 8 cores):
OMP_NUM_THREADS=8 taskset -c 0-7 ardal-switchprofile --threads 8 --same-grid
"""

import os, sys, time, math, argparse, json, random, platform, shutil, subprocess
from statistics import median
from pathlib import Path
from typing import List, Tuple, Dict, Any
import numpy as np
from ardal import Ardal, ardal




# ------------------------------- RNG ------------------------------------------
def _rng(seed: int):
    py = random.Random(seed)
    np.random.seed(seed & 0x7FFFFFFF)
    return py




# ---------------------------- Synthetic matrix -----------------------------------
def gen_iid(rows: int, cols: int, density: float, seed: int):
    py = _rng(seed)
    k = int(round(density*cols))
    mat = np.zeros((rows, cols), dtype=int)
    if k <= 0:
        return mat
    if k < cols/4:
        base = list(range(cols)) 
        for i in range(rows):
            indx = sorted(py.sample(base, k))
            for j in indx:
                mat[i, j] = 1
    else:
        for i in range(rows):
            for j in range(cols):
                if py.random() < density:
                    mat[i, j] = 1
    return mat




# ------------------------------ Pair sampler ----------------------------------
def sample_pairs(n: int, n_pairs: int, seed: int) -> List[Tuple[int,int]]:
    py = _rng(seed)
    out, seen = [], set()
    limit = n*(n-1)//2
    target = min(n_pairs, limit)
    while len(out) < target:
        i = py.randrange(n); j = py.randrange(n)
        if i < j and (i,j) not in seen:
            seen.add((i,j)); out.append((i,j))
    return out




# -------------------------------- Timing --------------------------------------
def time_median(fn, repeats=7, warmup=1) -> float:
    for _ in range(warmup): fn()
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter_ns(); fn(); t1 = time.perf_counter_ns()
        times.append((t1-t0)/1e9)
    return median(times)




# --------------------------- Adapters (ONE place to wire) ---------------------
class BackendAdapter:
    """
    Thin shim so benchmark code is backend-agnostic.

    The underlying object must implement:
      pairwise_hamming_pairs(pairs, threads)
      neighbourhood_radius(r, threads)
      neighbourhood_topk(k, threads)
    """
    def __init__(self, name: str, obj: Any, n_rows: int, n_cols: int):
        self.name = name; self.obj = obj; self.n_rows = n_rows; self.n_cols = n_cols

    def bench_pairwise(self, pairs, threads, adapter):
        def _run(): self.obj.hamming(True, threads, adapter)
        t = time_median(_run)
        return {"func":"hamming","time_s":t,"throughput":len(pairs)/t,"pairs":len(pairs), "adapter":adapter}

    def bench_neighbourhood(self, r, threads, adapter):
        def _run(): self.obj.neighbourhood(int(self.n_rows/2), r, True, threads, adapter)
        t = time_median(_run)
        return {"func":f"nh_radius_{r}","time_s":t,"throughput":self.n_rows/t,"pairs":self.n_rows, "adapter":adapter}

    def bench_topk(self, k, threads):
        def _run(): self.obj.neighbourhood_topk(k, threads=threads)
        t = time_median(_run)
        return {"func":f"nh_topk_{k}","time_s":t,"throughput":self.n_rows/t,"pairs":self.n_rows}




# ---- FILL THESE TWO STUBS WITH YOUR REAL CONSTRUCTORS (BitMatrix/Roaring etc.) ----
def _build_backend_from_mat(mat, n_rows, n_cols: int) -> BackendAdapter:
    hm_obj = ardal.HybridMatrix(mat, is_bitpacked=False, n_cols_bits=n_cols, use_roaring_if_sparse=True, density_threshold=1)
    return BackendAdapter("hybrid_matrix", hm_obj, n_rows, n_cols)


# ------------------------------ System stamp -----------------------------------
def system_fingerprint() -> Dict[str,str]:
    stamp = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
    }
    if shutil.which("lscpu"):
        try:
            out = subprocess.check_output(["lscpu"], text=True)
            for line in out.splitlines():
                if "Model name:" in line: stamp["cpu_model"] = line.split(":",1)[1].strip()
                if "Architecture:" in line: stamp["arch"] = line.split(":",1)[1].strip()
                if "CPU(s):" in line: stamp["cpus"] = line.split(":",1)[1].strip()
        except Exception:
            pass
    return stamp




# -------------------------- Switchpoint calculation ----------------------------
def _interp_zero(p0, s0, p1, s1):
    return p0 + (0.0 - s0) * (p1 - p0) / (s1 - s0)


def compute_switchpoints(rows: List[Dict[str,Any]]) -> Dict[Tuple[str,int], List[Dict[str,Any]]]:
    """
    rows: list with keys func, backend, n, m, p, threads, throughput

    Returns: {(func,threads): [ {m, p_star, note}, ... ]}, sorted by m.
    """
    # Index by (func,threads,m,p,backend)
    idx = {}
    for r in rows:
        idx[(r["func"], r["threads"], r["m"], r["p"], r["backend"])] = r

    # For each (func,threads,m), compute s = log2(bit/roaring) across p and find where s crosses 0
    comb = {}
    for (func,thr,m,p,be) in idx:
        comb.setdefault((func,thr,m), set()).add(p)

    out: Dict[Tuple[str,int], List[Dict[str,Any]]] = {}
    for (func,thr,m), pset in comb.items():
        ps = sorted(pset)
        series = []
        ok = True
        for p in ps:
            try:
                b = idx[(func,thr,m,p,"bit")]["throughput"]
                r = idx[(func,thr,m,p,"roaring")]["throughput"]
            except KeyError:
                ok = False; break
            series.append((p, math.log2(b/r)))
        if not ok: continue

        p_star, note = None, ""
        for (p0,s0),(p1,s1) in zip(series[:-1], series[1:]):
            if s0 == 0.0: p_star, note = p0, "exact"; break
            if (s0 < 0.0 and s1 > 0.0) or (s0 > 0.0 and s1 < 0.0):
                p_star, note = _interp_zero(p0,s0,p1,s1), "interpolated"; break
        if p_star is None:
            p_star = min(series, key=lambda t: abs(t[1]))[0]
            dom = "bit>roaring" if series[0][1] > 0 else "roaring>bit"
            note = f"no-crossing; nearest ({dom})"
        out.setdefault((func,thr), []).append({"m": m, "p_star": float(p_star), "note": note})

    for k in out:
        out[k] = sorted(out[k], key=lambda d: d["m"])
    return out




# --------------------------------- CLI ----------------------------------------
def print_csv_header():
    print("func,backend,n,m,p,threads,time_s,throughput,pairs", flush=True)


def print_csv_row(r):
    print("{func},{backend},{n},{m},{p:.6g},{threads},{time_s:.6f},{throughput:.3f},{pairs}".format(**r), flush=True)


def main():
    ap = argparse.ArgumentParser(prog="ardal-switchprofile",
        description="Probe machine-local backend switchpoints for Ardal hot operations.")
    ap.add_argument("--ops", nargs="+", choices=["pairwise","neighbourhood","topk"], required=True,
                    help="Which ops to sweep")
    ap.add_argument("--nrows", type=int, default=5000)
    ap.add_argument("--m", type=int, nargs="+", required=True)
    ap.add_argument("--p", type=float, nargs="+", required=True)
    ap.add_argument("--threads", type=int, default=max(1, int(os.environ.get("OMP_NUM_THREADS","1"))))
    ap.add_argument("--pairs", type=int, default=200000, help="#pairs for pairwise")
    ap.add_argument("--r", type=int, nargs="+", default=[25,50,100], help="radii for neighbourhood")
    ap.add_argument("--k", type=int, nargs="+", default=[10,50], help="K for top-K")
    ap.add_argument("--seed", type=int, default=1337)
    ap.add_argument("--repeats", type=int, default=7)
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--out", type=Path, help="Optional NDJSON of all measurements")
    ap.add_argument("--switchout", type=Path, default=Path("switch_thresholds_local.json"))
    args = ap.parse_args()

    # Stamp
    stamp = system_fingerprint()
    stamp["env_threads"] = args.threads

    out_f = open(args.out, "w") if args.out else None
    if out_f: out_f.write("")
    print_csv_header()

    all_rows: List[Dict[str,Any]] = []
    for m in sorted(set(args.m)):
        for p in sorted(set(args.p)):
            mat = gen_iid(args.nrows, m, p, seed=args.seed)
            backend = _build_backend_from_mat(mat, args.nrows, m)
            
            meta = {"n": args.nrows, "m": m, "p": p, "threads": args.threads}

            if "pairwise" in args.ops:
                pairs = sample_pairs(args.nrows, args.pairs, seed=args.seed)
                for adapter in ("bit", "roaring"):
                    res = backend.bench_pairwise(pairs, args.threads, adapter)
                    row = dict(meta, func=res["func"], backend=adapter,
                               time_s=res["time_s"], throughput=res["throughput"], pairs=res["pairs"])
                    print_csv_row(row); all_rows.append(row)
                    if out_f: out_f.write(json.dumps(row)+"\n")

            if "neighbourhood" in args.ops:
                for epsilon in args.r:
                    for adapter in ("bit", "roaring"):
                        res = backend.bench_neighbourhood(epsilon, args.threads, adapter)
                        row = dict(meta, func=res["func"], backend=adapter,
                                   time_s=res["time_s"], throughput=res["throughput"], pairs=res["pairs"])
                        print_csv_row(row); all_rows.append(row)
                        if out_f: out_f.write(json.dumps(row)+"\n")

            # if "topk" in args.ops:
            #     for k in args.k:
            #         for adapter in (bit, roa):
            #             res = adapter.bench_topk(k, args.threads)
            #             row = dict(meta, func=res["func"], backend=adapter.name,
            #                        time_s=res["time_s"], throughput=res["throughput"], pairs=res["pairs"])
            #             print_csv_row(row); all_rows.append(row)
            #             if out_f: out_f.write(json.dumps(row)+"\n")

    if out_f: out_f.close()

    switches = compute_switchpoints(all_rows)
    payload = {
        "stamp": stamp,
        "grid": {"n": args.nrows, "m": sorted(set(args.m)), "p": sorted(set(args.p)),
                 "repeats": args.repeats, "warmup": args.warmup},
        "switchpoints": {f"{func}|T{thr}": vals for (func,thr), vals in switches.items()}
    }
    with open(args.switchout, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"# wrote switchpoints => {args.switchout}", file=sys.stderr)


if __name__ == "__main__":
    main()
