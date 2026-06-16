#!/usr/bin/env python3
import argparse
import json
import random
import sys
from pathlib import Path


DNA = ("A", "C", "G", "T")


def build_alleles(n_alleles: int, n_chroms: int, seed: int):
    rng = random.Random(seed)
    alleles = []
    site_keys = []
    per_chrom = max(1, n_alleles // n_chroms)
    extra = n_alleles - (per_chrom * n_chroms)
    for chrom_idx in range(1, n_chroms + 1):
        chrom = f"chr{chrom_idx}"
        count = per_chrom + (1 if chrom_idx <= extra else 0)
        for i in range(count):
            pos = 100 + (i * 10)
            ref = rng.choice(DNA)
            alt = rng.choice([b for b in DNA if b != ref])
            allele_id = f"{chrom}.{pos}.{ref}.{alt}"
            alleles.append(allele_id)
            site_keys.append(f"{chrom}.{pos}")
    return alleles, site_keys


def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate simulated per-sample allele JSON files with missing positions."
    )
    ap.add_argument("--outdir", type=Path, default=Path("sim_jsons"),
                    help="Output directory for JSON files")
    ap.add_argument("--n-samples", type=int, default=10_000,
                    help="Number of sample JSON files to generate (default: 10000)")
    ap.add_argument("--n-alleles", type=int, default=20_000,
                    help="Number of allele IDs to generate (default: 20000)")
    ap.add_argument("--n-chroms", type=int, default=10,
                    help="Number of chromosomes to distribute alleles across (default: 10)")
    ap.add_argument("--seed", type=int, default=1337,
                    help="Random seed for reproducibility (default: 1337)")
    ap.add_argument("--min-alleles", type=int, default=50,
                    help="Minimum alleles per sample (default: 50)")
    ap.add_argument("--max-alleles", type=int, default=500,
                    help="Maximum alleles per sample (default: 500)")
    ap.add_argument("--min-missing", type=int, default=0,
                    help="Minimum missing positions per sample (default: 0)")
    ap.add_argument("--max-missing", type=int, default=2000,
                    help="Maximum missing positions per sample (default: 1000)")
    ap.add_argument("--sample-prefix", type=str, default="sample",
                    help="Sample ID/file prefix (default: sample)")
    return ap.parse_args()


def main():
    args = parse_args()
    if args.n_samples < 1 or args.n_alleles < 1:
        print("n-samples and n-alleles must be >= 1", file=sys.stderr)
        sys.exit(1)
    if args.min_alleles > args.max_alleles:
        print("min-alleles must be <= max-alleles", file=sys.stderr)
        sys.exit(1)
    if args.min_missing > args.max_missing:
        print("min-missing must be <= max-missing", file=sys.stderr)
        sys.exit(1)

    rng = random.Random(args.seed)
    alleles, site_keys = build_alleles(args.n_alleles, args.n_chroms, args.seed + 1)
    if len(alleles) != args.n_alleles:
        alleles = alleles[: args.n_alleles]
        site_keys = site_keys[: args.n_alleles]

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    for i in range(args.n_samples):
        sample_id = f"{args.sample_prefix}_{i:05d}"
        allele_count = rng.randint(args.min_alleles, min(args.max_alleles, len(alleles)))
        missing_count = rng.randint(args.min_missing, min(args.max_missing, len(site_keys)))

        allele_ids = rng.sample(alleles, allele_count)
        missing = rng.sample(site_keys, missing_count)

        payload = {
            "sample_id": sample_id,
            "alleles": allele_ids,
            "missing": missing,
        }
        path = outdir / f"{sample_id}.json"
        path.write_text(json.dumps(payload))

        if (i + 1) % 500 == 0 or i == args.n_samples - 1:
            print(f"\r[Write] {i+1}/{args.n_samples}", end="", flush=True)

    print(f"\nDone. Wrote {args.n_samples} JSONs to {outdir}")


if __name__ == "__main__":
    main()
