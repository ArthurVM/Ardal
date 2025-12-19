#!/usr/bin/env python3
"""
vcf_snps_to_json.py

Read one or more VCF/VCF.GZ files, filter by quality, keep SNPs only,
and write JSON payloads containing both allele calls and non-callable (missing)
positions inferred from `samtools depth -aa` BED/BED.GZ files. Sample IDs are
taken from the VCF header; merged VCFs emit one JSON per sample.

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
import re
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

DNA4 = {"A", "C", "G", "T"}
VCF_SUFFIXES: Tuple[str, ...] = (".raw.vcf", ".raw.vcf.gz", ".vcf.gz", ".vcf")
DEPTH_SUFFIXES: Tuple[str, ...] = (".depth.bed.gz", ".depth.bed", ".bed.gz", ".bed")

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

def apply_sample_regex(sample_name: str, sample_regex: str | None, sample_regex_group: int) -> str:
    if sample_regex:
        m = re.search(sample_regex, sample_name)
        if m:
            try:
                return m.group(sample_regex_group)
            except IndexError:
                pass
    return sample_name

def sample_alleles_from_vcf(
    vcf_path: Path,
    min_qual: float,
    allow_filtered: bool,
    accept_missing_qual: bool,
    min_dp: int | None,
    sample_regex: str | None,
    sample_regex_group: int,
) -> Dict[str, Set[str]]:
    sample_ids: List[str] | None = None
    allele_map: Dict[str, Set[str]] = {}
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == "#":
                if line.startswith("#CHROM"):
                    header_parts = line.rstrip("\n").split("\t")
                    if len(header_parts) <= 9:
                        raise SystemExit(f"No sample columns found in {vcf_path}")
                    header_samples = header_parts[9:]
                    sample_ids = [
                        apply_sample_regex(sample, sample_regex, sample_regex_group)
                        for sample in header_samples
                    ]
                    allele_map = {sample: set() for sample in sample_ids}
                continue

            if sample_ids is None:
                raise SystemExit(f"Missing #CHROM header in {vcf_path}")

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            # VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLES...]
            chrom, pos, _id, ref, alt_field, qual, filt, info, fmt = parts[:9]

            # Quick reject: if ref has length != 1 or non-ACGT, we can still keep
            # multi-allelic lines if any ALT is a length-1 A/C/G/T AND ref is length-1 A/C/G/T.
            # Otherwise skip early.
            if not (len(ref) == 1 and ref.upper() in DNA4):
                continue

            if not passes_filters(qual, filt, min_qual, allow_filtered, accept_missing_qual, min_dp, info):
                continue

            format_fields = fmt.split(":")
            try:
                gt_index = format_fields.index("GT")
            except ValueError:
                # No genotype; cannot attribute alleles to samples.
                continue

            alt_alleles = alt_field.split(",")
            if not alt_alleles:
                continue

            for sample_offset, sample_id in enumerate(sample_ids):
                sample_col_idx = 9 + sample_offset
                if sample_col_idx >= len(parts):
                    continue
                sample_fields = parts[sample_col_idx].split(":")
                if len(sample_fields) <= gt_index:
                    continue
                gt_value = sample_fields[gt_index]
                if not gt_value or gt_value == ".":
                    continue

                for allele_token in re.split(r"[|/]", gt_value):
                    if not allele_token or allele_token == ".":
                        continue
                    try:
                        allele_idx = int(allele_token)
                    except ValueError:
                        continue
                    if allele_idx <= 0:
                        continue
                    alt_idx = allele_idx - 1
                    if alt_idx < 0 or alt_idx >= len(alt_alleles):
                        continue
                    alt = alt_alleles[alt_idx]
                    if not is_simple_snp(ref, alt):
                        continue
                    allele_map[sample_id].add(f"{chrom}.{pos}.{ref.upper()}.{alt.upper()}")

    return allele_map

def strip_suffix(name: str, suffixes: Tuple[str, ...]) -> str:
    for suffix in sorted(suffixes, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name

def collect_files(path: Path, suffixes: Tuple[str, ...]) -> List[Path]:
    if path.is_file():
        if any(path.name.endswith(sfx) for sfx in suffixes):
            return [path]
        raise ValueError(f"{path} does not match expected suffixes: {suffixes}")

    files: List[Path] = []
    for candidate in path.rglob("*"):
        if candidate.is_file() and any(candidate.name.endswith(sfx) for sfx in suffixes):
            files.append(candidate)
    return sorted(files)

def pair_files(vcf_path: Path, depth_path: Path | None) -> List[Tuple[Path, Path | None]]:
    if vcf_path is None:
        raise ValueError("A VCF input path must be provided.")

    vcf_files = collect_files(vcf_path, VCF_SUFFIXES)
    if not vcf_files:
        raise SystemExit(f"No VCF files found under {vcf_path}")

    if depth_path is None:
        depth_path = vcf_path if vcf_path.is_dir() else vcf_path.parent
    depth_files = collect_files(depth_path, DEPTH_SUFFIXES)
    if not depth_files:
        print(f"Warning: no BED/BED.GZ depth files found under {depth_path}; missing lists will be empty.",
              file=sys.stderr)

    depth_lookup: Dict[str, Path] = {}
    for depth in depth_files:
        key = strip_suffix(depth.name, DEPTH_SUFFIXES)
        if key in depth_lookup:
            raise SystemExit(f"Duplicate depth files with base '{key}': {depth_lookup[key]} vs {depth}")
        depth_lookup[key] = depth

    pairs: List[Tuple[Path, Path | None]] = []
    missing_depth: List[str] = []
    for vcf in vcf_files:
        key = strip_suffix(vcf.name, VCF_SUFFIXES)
        depth = depth_lookup.get(key)
        if depth is None:
            missing_depth.append(vcf.name)
        pairs.append((vcf, depth))

    if missing_depth:
        print(
            "Warning: missing depth BED files for the following VCFs (matched by basename): "
            + ", ".join(missing_depth),
            file=sys.stderr,
        )

    return sorted(pairs, key=lambda pair: str(pair[0]))

def missing_positions_from_bed(depth_path: Path, min_depth: int) -> List[str]:
    missing: List[str] = []
    seen: Set[str] = set()
    with open_maybe_gzip(depth_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, pos, depth_str = parts[:3]
            try:
                depth_val = int(float(depth_str))
            except ValueError:
                continue
            if depth_val >= min_depth:
                continue
            key = f"{chrom}.{pos}"
            if key not in seen:
                seen.add(key)
                missing.append(key)
    return missing

def main():
    ap = argparse.ArgumentParser(
        description="Collect SNP IDs and missing coverage positions into JSON payloads."
    )
    ap.add_argument("-v", "--vcf_in", type=Path, required=True,
                    help="VCF file or directory containing VCF/VCF.GZ files")
    ap.add_argument("-d", "--depth_in", type=Path, default=None,
                    help="Depth BED/BED.GZ file(s) from `samtools depth -aa`; defaults to scanning the VCF directory")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("./"), help="Output directory")
    ap.add_argument("--min-qual", type=float, default=30.0, help="Minimum QUAL (default: 30)")
    ap.add_argument("--allow-filtered", action="store_true",
                    help="Allow FILTER values other than PASS/. (default: require PASS or .)")
    ap.add_argument("--accept-missing-qual", action="store_true",
                    help="Treat missing QUAL '.' as passing (default: fail)")
    ap.add_argument("--min-dp", type=int, default=None,
                    help="Require INFO/DP >= this value if present; if absent, variant fails this check")
    ap.add_argument("--sample-regex", type=str, default=None,
                    help="Optional regex to transform sample IDs from the VCF header")
    ap.add_argument("--sample-regex-group", type=int, default=1,
                    help="Capture group index to take from --sample-regex (default: 1)")
    ap.add_argument("--min-missing-depth", type=int, default=10,
                    help="Positions with coverage < this value are marked as missing (default: 10)")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    inputs = pair_files(args.vcf_in, args.depth_in)
    total = len(inputs)
    for idx, (vcf, depth) in enumerate(inputs, start=1):
        allele_map = sample_alleles_from_vcf(
            vcf_path=vcf,
            min_qual=args.min_qual,
            allow_filtered=args.allow_filtered,
            accept_missing_qual=args.accept_missing_qual,
            min_dp=args.min_dp,
            sample_regex=args.sample_regex,
            sample_regex_group=args.sample_regex_group,
        )
        print(f"\r {idx}/{total} : {vcf.name} ({len(allele_map)} samples)", end="             ", flush=True)

        if depth is None:
            missing: List[str] = []
        else:
            missing = missing_positions_from_bed(depth, args.min_missing_depth)

        for sample_id, alleles in allele_map.items():
            out_path = args.outdir / f"{sample_id}_snps.json"
            if out_path.exists():
                continue

            payload = {
                "alleles": sorted(alleles),
                "missing": missing,
            }

            with open(out_path, "w") as outfh:
                json.dump(payload, outfh, indent=2)

if __name__ == "__main__":
    main()
