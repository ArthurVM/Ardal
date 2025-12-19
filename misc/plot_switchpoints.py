#!/usr/bin/env python3
"""
Plot Ardal backend switchpoints.

Subcommands
-----------
1) decision: heatmap of log2(bit/roaring) across (m, p) for one operation+thread count,
   with zero-level contour (the switch line) and measured grid points.
   Input: NDJSON measurements produced by ardal-switchprofile (or ardal-switchprobe).

2) boundary: p*(m) switchpoint curve from switch_thresholds_local.json.
   Input: the JSON file emitted by your profiler.

Usage
-----
# Decision field for pairwise (T=10):
python plot_switchpoints.py decision \
  --measurements results.ndjson \
  --op pairwise_hamming --threads 10 \
  --out decision_pairwise_T10.png

# Decision field with overlay of precomputed p*(m) from JSON:
python plot_switchpoints.py decision \
  --measurements results.ndjson \
  --op nh_radius_50 --threads 10 \
  --profile switch_thresholds_local.json \
  --out decision_radius50_T10.png

# Switchpoint curve for pairwise (exact key seen in your JSON):
python plot_switchpoints.py boundary \
  --profile switch_thresholds_local.json \
  --key "hamming|T10" \
  --out pstar_hamming_T10.png

# Or list available keys in the JSON:
python plot_switchpoints.py boundary --profile switch_thresholds_local.json --list-keys
"""
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import CenteredNorm

# ------------- helpers -------------
def _read_ndjson(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            rows.append(json.loads(line))
    df = pd.DataFrame(rows)
    # normalize column names if needed
    need = {"func","backend","n","m","p","threads","throughput"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"NDJSON missing columns: {missing}")
    return df

def _avg_cell(df):
    # average duplicates within a cell (func, backend, n, m, p, threads)
    grp = df.groupby(["func","backend","n","m","p","threads"], as_index=False)["throughput"].mean()
    return grp

def _pivot_log2_ratio(df):
    # df is a single func+threads subset, both backends present
    bit = df[df["backend"]=="bit"].rename(columns={"throughput":"bit"})
    roa = df[df["backend"]=="roaring"].rename(columns={"throughput":"roaring"})
    merged = pd.merge(bit[["m","p","bit"]], roa[["m","p","roaring"]], on=["m","p"], how="inner")
    merged["log2_ratio"] = np.log2(merged["bit"] / merged["roaring"])
    return merged

def _grid_for_pcolormesh(tab: pd.DataFrame):
    # We want a full (p×m) grid for pcolormesh; missing cells become NaN.
    ms = np.unique(tab["m"].to_numpy())
    ps = np.unique(tab["p"].to_numpy())
    M, P = np.meshgrid(ms, ps)
    Z = np.full_like(M, np.nan, dtype=float)
    # map (m,p) -> value
    key = {(int(m), float(p)): v for m,p,v in zip(tab["m"], tab["p"], tab["log2_ratio"])}
    for i, p in enumerate(ps):
        for j, m in enumerate(ms):
            v = key.get((int(m), float(p)))
            if v is not None:
                Z[i,j] = v
    return ms, ps, Z

def _interp_zero(p0,s0,p1,s1):
    return p0 + (0.0 - s0) * (p1 - p0) / (s1 - s0)

def _zero_contour_from_tab(tab: pd.DataFrame):
    # Produce p*(m) from the same merged table (per m, scan p in order for sign change)
    out = []
    for m, sub in tab.groupby("m"):
        sub = sub.sort_values("p")
        ps = sub["p"].to_numpy()
        ss = sub["log2_ratio"].to_numpy()
        pstar = None
        for (p0,s0),(p1,s1) in zip(ps[:-1], ps[1:]):
            if s0 == 0.0:
                pstar = p0; break
            if (s0<0 and s1>0) or (s0>0 and s1<0):
                pstar = _interp_zero(p0,s0,p1,s1); break
        if pstar is None:
            # pick nearest-to-zero
            pstar = ps[np.argmin(np.abs(ss))]
        out.append((m, float(pstar)))
    out = sorted(out, key=lambda t:t[0])
    return np.array([m for m,_ in out]), np.array([p for _,p in out])

# ------------- plotting -------------
def plot_decision_field(meas_path: Path, op: str, threads: int, out_path: Path, profile_path: Path|None):
    df = _read_ndjson(meas_path)
    # filter to operation + threads
    sub = df[(df["func"]==op) & (df["threads"]==threads)]
    if sub.empty:
        raise SystemExit(f"No rows for func='{op}', threads={threads} in {meas_path}")
    sub = _avg_cell(sub)
    merged = _pivot_log2_ratio(sub)
    ms, ps, Z = _grid_for_pcolormesh(merged)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    # pcolormesh expects bin edges; build simple edges in log space
    m_edges = np.geomspace(ms.min(), ms.max(), num=len(ms)+1)
    p_edges = np.geomspace(ps.min(), ps.max(), num=len(ps)+1)

    # heatmap centered at 0 (bit==roaring)
    pcm = ax.pcolormesh(m_edges, p_edges, Z, norm=CenteredNorm(0.0), shading="auto")
    cb = fig.colorbar(pcm, ax=ax, label="log2(throughput_bit / throughput_roaring)")

    # scatter actual measured cells
    ax.scatter(merged["m"], merged["p"], s=16, marker="o", alpha=0.7)

    # zero contour estimated from same table
    m_star, p_star = _zero_contour_from_tab(merged)
    ax.plot(m_star, p_star, linestyle="-", linewidth=2)

    # optional overlay from precomputed profile JSON (if provided)
    if profile_path is not None:
        prof = json.loads(Path(profile_path).read_text())
        key = f"{op}|T{threads}"
        if key in prof.get("switchpoints", {}):
            rows = prof["switchpoints"][key]
            m_prof = np.array([r["m"] for r in rows])
            p_prof = np.array([r["p_star"] for r in rows])
            ax.plot(m_prof, p_prof, linestyle="--", linewidth=1.5)  # overlay
        else:
            print(f"[warn] key '{key}' not found in profile; skipping overlay")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Number of columns (m)")
    ax.set_ylabel("Density p")
    ax.set_title(f"Decision field: {op}  (threads={threads})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(f"[ok] wrote {out_path}")

def plot_boundary(profile_path: Path, key: str, out_path: Path, list_keys: bool):
    prof = json.loads(Path(profile_path).read_text())
    if list_keys:
        print("Available keys in 'switchpoints':")
        for k in sorted(prof.get("switchpoints", {}).keys()):
            print("  ", k)
        return
    rows = prof.get("switchpoints", {}).get(key)
    if not rows:
        raise SystemExit(f"No switchpoints for key '{key}' in {profile_path}")

    ms = np.array([r["m"] for r in rows], dtype=float)
    ps = np.array([r["p_star"] for r in rows], dtype=float)
    notes = [r.get("note","") for r in rows]

    # choose marker style based on note
    markers = []
    for note in notes:
        if note.startswith("no-crossing"):
            markers.append(("o", False))  # open circle
        elif note == "exact":
            markers.append(("s", True))   # filled square
        else:
            markers.append(("o", True))   # filled circle

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    # draw segments in order
    ax.plot(ms, ps, linewidth=1.5)
    # draw points with per-point fill
    for m,p,(mk,filled) in zip(ms, ps, markers):
        ax.plot([m],[p], marker=mk, markersize=6,
                markerfacecolor="none" if not filled else None)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Number of columns (m)")
    ax.set_ylabel("Switch density p*")
    ax.set_title(f"Switchpoint curve: {key}")
    # small legend hint
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0],[0], marker='o', color='k', label='interpolated p*', markersize=6, linewidth=0),
        Line2D([0],[0], marker='s', color='k', label='exact tie', markersize=6, linewidth=0),
        Line2D([0],[0], marker='o', markerfacecolor='none', color='k',
               label='no crossing (nearest)', markersize=6, linewidth=0),
    ]
    ax.legend(handles=legend_elems, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(f"[ok] wrote {out_path}")

# ------------- CLI -------------
def main():
    ap = argparse.ArgumentParser(prog="plot_switchpoints.py")
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("decision", help="Decision field heatmap with zero contour")
    d.add_argument("--measurements", type=Path, required=True, help="NDJSON from profiler")
    d.add_argument("--op", required=True, help="func name (e.g. pairwise_hamming, nh_radius_50, nh_topk_10)")
    d.add_argument("--threads", type=int, required=True)
    d.add_argument("--out", type=Path, required=True)
    d.add_argument("--profile", type=Path, help="Optional switch_thresholds_local.json to overlay p*(m)")

    b = sub.add_parser("boundary", help="Switchpoint p*(m) from JSON")
    b.add_argument("--profile", type=Path, required=True)
    b.add_argument("--key", help="Exact key in JSON (e.g. 'hamming|T10' or 'nh_radius_50|T10')")
    b.add_argument("--out", type=Path)
    b.add_argument("--listkeys", action="store_true")

    args = ap.parse_args()
    if args.cmd == "decision":
        plot_decision_field(args.measurements, args.op, args.threads, args.out, args.profile)
    else:
        if args.listkeys:
            plot_boundary(args.profile, key="", out_path=Path("noop.png"), list_keys=True)
        else:
            if not args.key or not args.out:
                raise SystemExit("--key and --out are required unless --list-keys is provided")
            plot_boundary(args.profile, key=args.key, out_path=args.out, list_keys=False)

if __name__ == "__main__":
    main()
