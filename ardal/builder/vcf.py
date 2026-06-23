"""
VCF parsing and filtering helpers for Ardal sample JSON generation.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple

DNA4 = {"A", "C", "G", "T"}
VCF_SUFFIXES: Tuple[str, ...] = (".raw.vcf", ".raw.vcf.gz", ".vcf.gz", ".vcf", ".filt.vcf", ".filt.vcf.gz")
DEPTH_SUFFIXES: Tuple[str, ...] = (".depth.bed.gz", ".depth.bed", ".bed.gz", ".bed", ".depth")
READ_EXCEPTIONS = (gzip.BadGzipFile, EOFError, OSError)
AlleleRecord = Tuple[str, int, int, str, str]
NonsynonymousRecord = Tuple[AlleleRecord, str, int, str, str]

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

def sample_ids_from_vcf_header(
    vcf_path: Path,
    sample_regex: str | None,
    sample_regex_group: int,
) -> List[str]:
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if not line:
                continue
            if not line.startswith("#"):
                break
            if line.startswith("#CHROM"):
                header_parts = line.rstrip("\n").split("\t")
                if len(header_parts) <= 9:
                    raise SystemExit(f"No sample columns found in {vcf_path}")
                sample_ids = [
                    apply_sample_regex(sample, sample_regex, sample_regex_group)
                    for sample in header_parts[9:]
                ]
                return sorted(set(sample_ids))
    raise SystemExit(f"Missing #CHROM header in {vcf_path}")

def sample_alleles_from_vcf(
    vcf_path: Path,
    min_qual: float,
    allow_filtered: bool,
    accept_missing_qual: bool,
    min_dp: int | None,
    sample_regex: str | None,
    sample_regex_group: int,
    annotator: VariantAnnotator | None = None,
) -> Dict[str, Dict[str, Set]]:
    sample_ids: List[str] | None = None
    allele_map: Dict[str, Dict[str, Set]] = {}
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
                    allele_map = {
                        sample: {
                            "alleles": set(),
                            "genic": set(),
                            "nonsynonymous": set(),
                        }
                        for sample in sample_ids
                    }
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

            try:
                pos_int = int(pos)
            except ValueError:
                continue
            annotations_by_alt: Dict[int, Tuple[AlleleRecord, bool, Set[NonsynonymousRecord]]] = {}
            for alt_idx, alt in enumerate(alt_alleles):
                if not is_simple_snp(ref, alt):
                    continue
                allele_record: AlleleRecord = (chrom, pos_int, pos_int, ref.upper(), alt.upper())
                is_genic = False
                nonsynonymous_records: Set[NonsynonymousRecord] = set()
                if annotator is not None and annotator.is_genic(chrom, pos_int):
                    is_genic = True
                    for gene_id, aa_pos, ref_aa, alt_aa in annotator.nonsynonymous_aa_changes(chrom, pos_int, ref, alt):
                        nonsynonymous_records.add((allele_record, gene_id, aa_pos, ref_aa, alt_aa))
                annotations_by_alt[alt_idx] = (allele_record, is_genic, nonsynonymous_records)

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
                    annotations = annotations_by_alt.get(alt_idx)
                    if annotations is None:
                        continue
                    allele_record, is_genic, nonsynonymous_records = annotations
                    allele_map[sample_id]["alleles"].add(allele_record)
                    if is_genic:
                        allele_map[sample_id]["genic"].add(allele_record)
                    allele_map[sample_id]["nonsynonymous"].update(nonsynonymous_records)

    return allele_map

def samples_with_candidate_alleles_from_vcf(
    vcf_path: Path,
    min_qual: float,
    allow_filtered: bool,
    accept_missing_qual: bool,
    min_dp: int | None,
    sample_regex: str | None,
    sample_regex_group: int,
) -> Set[str]:
    """Fast preflight: return samples with at least one passing SNP genotype."""
    sample_ids: List[str] | None = None
    samples_with_alleles: Set[str] = set()
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == "#":
                if line.startswith("#CHROM"):
                    header_parts = line.rstrip("\n").split("\t")
                    if len(header_parts) <= 9:
                        raise SystemExit(f"No sample columns found in {vcf_path}")
                    sample_ids = [
                        apply_sample_regex(sample, sample_regex, sample_regex_group)
                        for sample in header_parts[9:]
                    ]
                continue

            if sample_ids is None:
                raise SystemExit(f"Missing #CHROM header in {vcf_path}")

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            _chrom, _pos, _id, ref, alt_field, qual, filt, info, fmt = parts[:9]
            if not (len(ref) == 1 and ref.upper() in DNA4):
                continue
            if not passes_filters(qual, filt, min_qual, allow_filtered, accept_missing_qual, min_dp, info):
                continue

            callable_alt_indices = {
                alt_idx
                for alt_idx, alt in enumerate(alt_field.split(","))
                if is_simple_snp(ref, alt)
            }
            if not callable_alt_indices:
                continue

            format_fields = fmt.split(":")
            try:
                gt_index = format_fields.index("GT")
            except ValueError:
                continue

            for sample_offset, sample_id in enumerate(sample_ids):
                if sample_id in samples_with_alleles:
                    continue
                sample_col_idx = 9 + sample_offset
                if sample_col_idx >= len(parts):
                    continue
                sample_fields = parts[sample_col_idx].split(":")
                if len(sample_fields) <= gt_index:
                    continue
                for allele_token in re.split(r"[|/]", sample_fields[gt_index]):
                    if not allele_token or allele_token == ".":
                        continue
                    try:
                        allele_idx = int(allele_token)
                    except ValueError:
                        continue
                    if allele_idx > 0 and allele_idx - 1 in callable_alt_indices:
                        samples_with_alleles.add(sample_id)
                        break

            if len(samples_with_alleles) == len(set(sample_ids)):
                break

    return samples_with_alleles

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
    # print(pairs)
    return sorted(pairs, key=lambda pair: str(pair[0]))
