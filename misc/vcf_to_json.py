#!/usr/bin/env python3
"""
vcf_snps_to_json.py

Read one or more VCF/VCF.GZ files, filter by quality, keep SNPs only,
and write a JSON mapping: {sample_id: [ "CHROM.POS.REF.ALT", ... ], ... }.

Defaults:
  - QUAL >= 30
  - FILTER must be PASS or '.'
  - SNPs only: REF and each ALT are single bases in {A,C,G,T}

You can relax/tweak via CLI flags below.
"""

from __future__ import annotations
import argparse
import gzip
import json
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

DNA4 = {"A", "C", "G", "T"}

def open_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")

def is_simple_snp(ref: str, alt: str) -> bool:
    ref = ref.upper()
    alt = alt.upper()
    return len(ref) == 1 and len(alt) == 1 and ref in DNA4 and alt in DNA4

def parse_info_dp(info: str) -> int | None:
    # INFO like "DP=42;AC=1" → 42
    for part in info.split(";"):
        if part.startswith("DP="):
            val = part[3:]
            try:
                return int(val)
            except ValueError:
                return None
    return None

def passes_filters(
    qual_field: str,
    filter_field: str,
    min_qual: float,
    allow_filtered: bool,
    accept_missing_qual: bool,
    min_dp: int | None,
    info_field: str
) -> bool:
    # QUAL
    if qual_field == ".":
        if not accept_missing_qual:
            return False
    else:
        try:
            if float(qual_field) < min_qual:
                return False
        except ValueError:
            return False

    # FILTER
    if not allow_filtered:
        if filter_field not in ("PASS", "."):
            return False

    # INFO DP (optional)
    if min_dp is not None:
        dp = parse_info_dp(info_field)
        if dp is None or dp < min_dp:
            return False

    return True

def snp_ids_from_vcf(
    vcf_path: Path,
    min_qual: float,
    allow_filtered: bool,
    accept_missing_qual: bool,
    min_dp: int | None,
) -> List[str]:
    snp_ids: List[str] = []
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            # VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT ...]
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                # malformed
                continue
            chrom, pos, _id, ref, alt_field, qual, filt, info = parts[:8]

            # Quick reject: if ref has length != 1 or non-ACGT, we can still keep
            # multi-allelic lines if any ALT is a length-1 A/C/G/T AND ref is length-1 A/C/G/T.
            # Otherwise skip early.
            if not (len(ref) == 1 and ref.upper() in DNA4):
                continue

            if not passes_filters(qual, filt, min_qual, allow_filtered, accept_missing_qual, min_dp, info):
                continue

            # ALT may be comma-separated; generate one ID per simple SNP ALT.
            for alt in alt_field.split(","):
                if is_simple_snp(ref, alt):
                    snp_ids.append(f"{chrom}.{pos}.{ref.upper()}.{alt.upper()}")

    return snp_ids

def derive_sample_id(path: Path, sample_regex: str | None, sample_regex_group: int) -> str:
    base = path.name
    # Strip .vcf or .vcf.gz
    if base.endswith(".vcf.gz"):
        base = base[:-7]
    elif base.endswith(".vcf"):
        base = base[:-4]
    if sample_regex:
        m = re.search(sample_regex, base)
        if m:
            try:
                return m.group(sample_regex_group)
            except IndexError:
                pass
    return base

def iter_vcf_files(input_path: Path) -> Iterable[Path]:
    if input_path.is_file():
        if input_path.suffix in (".vcf",) or input_path.suffixes[-2:] == [".vcf", ".gz"]:
            yield input_path
        else:
            raise ValueError(f"Unsupported file extension for {input_path}")
    else:
        for p in sorted(input_path.rglob("*")):
            if p.is_file() and (p.suffix == ".vcf" or p.suffixes[-2:] == [".vcf", ".gz"]):
                yield p

def main():
    ap = argparse.ArgumentParser(
        description="Collect SNP IDs from VCFs into JSON: {sample_id: [chr.pos.ref.alt, ...], ...}"
    )
    ap.add_argument("input", type=Path, help="VCF file or directory containing VCF/VCF.GZ files")
    ap.add_argument("-o", "--outdir", type=Path, default="./", help="Output directory")
    ap.add_argument("--min-qual", type=float, default=30.0, help="Minimum QUAL (default: 30)")
    ap.add_argument("--allow-filtered", action="store_true",
                    help="Allow FILTER values other than PASS/. (default: require PASS or .)")
    ap.add_argument("--accept-missing-qual", action="store_true",
                    help="Treat missing QUAL '.' as passing (default: fail)")
    ap.add_argument("--min-dp", type=int, default=None,
                    help="Require INFO/DP >= this value if present; if absent, variant fails this check")
    ap.add_argument("--sample-regex", type=str, default=None,
                    help="Optional regex to extract sample_id from filename (applied to basename without .vcf[.gz])")
    ap.add_argument("--sample-regex-group", type=int, default=1,
                    help="Capture group index to take from --sample-regex (default: 1)")
    args = ap.parse_args()

    inputs = list(iter_vcf_files(args.input))
    if not inputs:
        raise SystemExit(f"No VCF files found in {args.input}")

    il = len(inputs)
    for i, vcf in enumerate(inputs):
        
        result: Dict[str, List[str]] = {}
        
        sample_id = derive_sample_id(vcf, args.sample_regex, args.sample_regex_group)
        print(f"\r {i}/{il} : {sample_id}", end="             ", flush=True)
        
        out = f"{args.outdir}/{sample_id}_snps.json"
        
        if os.path.exists(out):
            continue
            
        snps = snp_ids_from_vcf(
            vcf_path=vcf,
            min_qual=args.min_qual,
            allow_filtered=args.allow_filtered,
            accept_missing_qual=args.accept_missing_qual,
            min_dp=args.min_dp,
        )
        # Deduplicate and sort for stability
        result[sample_id] = sorted(set(snps))

        with open(out, "w") as outfh:
            json.dump(result, outfh, indent=2)
        # print(f"Wrote {out} with {len(result[sample_id])} SNPs.")

if __name__ == "__main__":
    main()
