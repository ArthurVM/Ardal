#!/usr/bin/env python3
"""
vcf_snps_to_json.py

Read one or more VCF/VCF.GZ files, filter by quality, keep SNPs only,
and write JSON payloads containing both allele calls and non-callable (missing)
positions inferred from `samtools depth -aa` BED/BED.GZ files. Sample IDs are
taken from the VCF header; merged VCFs emit one JSON per sample.

If a matching reference FASTA and GFF are supplied, each sample also gets
separate JSON outputs for genic SNP calls and non-synonymous amino-acid changes.

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
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

DNA4 = {"A", "C", "G", "T"}
VCF_SUFFIXES: Tuple[str, ...] = (".raw.vcf", ".raw.vcf.gz", ".vcf.gz", ".vcf", ".filt.vcf", ".filt.vcf.gz")
DEPTH_SUFFIXES: Tuple[str, ...] = (".depth.bed.gz", ".depth.bed", ".bed.gz", ".bed", ".depth")
READ_EXCEPTIONS = (gzip.BadGzipFile, EOFError, OSError)
ALLELE_OUTPUTS: Tuple[Tuple[str, str], ...] = (
    ("alleles", "snps"),
    ("genic_snps", "genic_snps"),
    ("nonsynonymous", "nonsynonymous"),
)
GENIC_FEATURES = {
    "gene",
    "pseudogene",
    "mrna",
    "transcript",
    "cds",
    "exon",
    "trna",
    "rrna",
    "ncrna",
}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class CdsFeature:
    chrom: str
    start: int
    end: int
    strand: str
    parent: str
    gene_id: str


@dataclass
class TranscriptCds:
    gene_id: str
    chrom: str
    strand: str
    sequence: str
    genomic_to_cds_index: Dict[int, int]
    cds_index_to_genomic: Dict[int, int]

def open_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")

def is_simple_snp(ref: str, alt: str) -> bool:
    ref = ref.upper()
    alt = alt.upper()
    return len(ref) == 1 and len(alt) == 1 and ref in DNA4 and alt in DNA4

def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1].upper()

def read_fasta(path: Path) -> Dict[str, str]:
    sequences: Dict[str, List[str]] = {}
    current_id: str | None = None
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                sequences.setdefault(current_id, [])
                continue
            if current_id is None:
                raise ValueError(f"FASTA sequence encountered before header in {path}")
            sequences[current_id].append(line.upper())
    return {seq_id: "".join(parts) for seq_id, parts in sequences.items()}

def parse_gff_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_field.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs

def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[Tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged

def interval_contains(intervals: List[Tuple[int, int]], pos: int) -> bool:
    # GFF inputs here are 1-based inclusive, same as VCF POS.
    for start, end in intervals:
        if pos < start:
            return False
        if start <= pos <= end:
            return True
    return False

def translate_codon(codon: str) -> str | None:
    codon = codon.upper()
    if len(codon) != 3 or any(base not in DNA4 for base in codon):
        return None
    return CODON_TABLE.get(codon)

def preferred_gff_label(attrs: Dict[str, str], fallback: str, gene_id_attribute: str) -> str:
    for key in (gene_id_attribute, "gene", "locus_tag", "Name", "ID", "Parent"):
        value = attrs.get(key)
        if value:
            return value.split(",")[0]
    return fallback

def resolve_gene_id(
    feature_id: str,
    parent_by_id: Dict[str, str],
    gene_id_by_id: Dict[str, str],
) -> str | None:
    seen: Set[str] = set()
    current: str | None = feature_id
    while current and current not in seen:
        seen.add(current)
        gene_id = gene_id_by_id.get(current)
        if gene_id:
            return gene_id
        current = parent_by_id.get(current)
    return None

def load_gff_annotations(
    path: Path,
    gene_id_attribute: str,
) -> Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, List[CdsFeature]]]:
    genic_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    cds_by_parent: Dict[str, List[CdsFeature]] = defaultdict(list)
    parent_by_id: Dict[str, str] = {}
    gene_id_by_id: Dict[str, str] = {}
    cds_records: List[Tuple[str, int, int, str, Dict[str, str]]] = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature_type, start_s, end_s, _score, strand, _phase, attrs_s = parts[:9]
            feature_type_l = feature_type.lower()
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            if start > end:
                start, end = end, start

            attrs = parse_gff_attributes(attrs_s)
            feature_id = attrs.get("ID")
            parents = attrs.get("Parent")
            if feature_id and parents:
                parent_by_id[feature_id] = parents.split(",")[0]
            if feature_id:
                if attrs.get(gene_id_attribute):
                    gene_id_by_id[feature_id] = attrs[gene_id_attribute].split(",")[0]
                elif feature_type_l == "gene":
                    gene_id_by_id[feature_id] = feature_id

            if feature_type_l in GENIC_FEATURES:
                genic_intervals[chrom].append((start, end))

            if feature_type_l == "cds":
                cds_records.append((chrom, start, end, strand, attrs))

    for chrom, start, end, strand, attrs in cds_records:
        parent = attrs.get("Parent") or attrs.get("ID") or attrs.get("gene") or f"{chrom}:{start}-{end}:{strand}"
        fallback_gene_id = preferred_gff_label(attrs, parent, gene_id_attribute)
        # A CDS may list multiple transcript parents.
        for parent_id in parent.split(","):
            gene_id = (
                attrs.get(gene_id_attribute)
                or attrs.get("gene")
                or attrs.get("locus_tag")
                or resolve_gene_id(parent_id, parent_by_id, gene_id_by_id)
                or fallback_gene_id
            ).split(",")[0]
            cds_by_parent[parent_id].append(
                CdsFeature(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    parent=parent_id,
                    gene_id=gene_id,
                )
            )

    return (
        {chrom: merge_intervals(intervals) for chrom, intervals in genic_intervals.items()},
        dict(cds_by_parent),
    )

def build_transcript_cds(reference: Dict[str, str], cds_by_parent: Dict[str, List[CdsFeature]]) -> Dict[str, List[TranscriptCds]]:
    transcripts_by_chrom: Dict[str, List[TranscriptCds]] = defaultdict(list)
    for _parent, features in cds_by_parent.items():
        features = [feature for feature in features if feature.chrom in reference]
        if not features:
            continue
        chroms = {feature.chrom for feature in features}
        strands = {feature.strand for feature in features}
        if len(chroms) != 1 or len(strands) != 1:
            continue

        chrom = features[0].chrom
        strand = features[0].strand
        gene_id = features[0].gene_id
        seq_parts: List[str] = []
        genomic_to_cds_index: Dict[int, int] = {}
        cds_index_to_genomic: Dict[int, int] = {}
        ref_seq = reference[chrom]

        cds_offset = 0
        if strand == "-":
            ordered = sorted(features, key=lambda feature: feature.start, reverse=True)
            for feature in ordered:
                segment = ref_seq[feature.start - 1:feature.end]
                oriented = reverse_complement(segment)
                for offset, _base in enumerate(oriented):
                    genomic_pos = feature.end - offset
                    cds_idx = cds_offset + offset
                    genomic_to_cds_index[genomic_pos] = cds_idx
                    cds_index_to_genomic[cds_idx] = genomic_pos
                seq_parts.append(oriented)
                cds_offset += len(oriented)
        else:
            ordered = sorted(features, key=lambda feature: feature.start)
            for feature in ordered:
                segment = ref_seq[feature.start - 1:feature.end].upper()
                for offset, _base in enumerate(segment):
                    genomic_pos = feature.start + offset
                    cds_idx = cds_offset + offset
                    genomic_to_cds_index[genomic_pos] = cds_idx
                    cds_index_to_genomic[cds_idx] = genomic_pos
                seq_parts.append(segment)
                cds_offset += len(segment)

        sequence = "".join(seq_parts)
        if len(sequence) >= 3:
            transcripts_by_chrom[chrom].append(
                TranscriptCds(
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                    sequence=sequence,
                    genomic_to_cds_index=genomic_to_cds_index,
                    cds_index_to_genomic=cds_index_to_genomic,
                )
            )

    return dict(transcripts_by_chrom)

class VariantAnnotator:
    def __init__(self, fasta_path: Path, gff_path: Path, gene_id_attribute: str = "gene_id"):
        self.reference = read_fasta(fasta_path)
        self.genic_intervals, cds_by_parent = load_gff_annotations(gff_path, gene_id_attribute)
        self.transcripts_by_chrom = build_transcript_cds(self.reference, cds_by_parent)

    def is_genic(self, chrom: str, pos: int) -> bool:
        return interval_contains(self.genic_intervals.get(chrom, []), pos)

    def nonsynonymous_aa_changes(self, chrom: str, pos: int, ref: str, alt: str) -> Set[Tuple[str, int, str, str]]:
        ref = ref.upper()
        alt = alt.upper()
        if chrom not in self.reference or pos < 1 or pos > len(self.reference[chrom]):
            return set()
        if self.reference[chrom][pos - 1].upper() != ref:
            return set()

        changes: Set[Tuple[str, int, str, str]] = set()
        for transcript in self.transcripts_by_chrom.get(chrom, []):
            cds_idx = transcript.genomic_to_cds_index.get(pos)
            if cds_idx is None:
                continue
            codon_start = (cds_idx // 3) * 3
            codon = transcript.sequence[codon_start:codon_start + 3]
            ref_aa = translate_codon(codon)
            if ref_aa is None:
                continue
            alt_coding_base = reverse_complement(alt) if transcript.strand == "-" else alt
            alt_codon = codon[:cds_idx - codon_start] + alt_coding_base + codon[cds_idx - codon_start + 1:]
            alt_aa = translate_codon(alt_codon)
            if alt_aa is not None and alt_aa != ref_aa:
                aa_pos = (codon_start // 3) + 1
                changes.add((transcript.gene_id, aa_pos, ref_aa, alt_aa))
        return changes

    def aa_missing_positions(self, missing: List[List], threshold: int) -> List[List]:
        missing_by_chrom: Dict[str, Set[int]] = defaultdict(set)
        for chrom, pos in missing:
            try:
                missing_by_chrom[str(chrom)].add(int(pos))
            except ValueError:
                continue

        aa_missing: Set[Tuple[str, str, int]] = set()
        for chrom, transcripts in self.transcripts_by_chrom.items():
            chrom_missing = missing_by_chrom.get(chrom)
            if not chrom_missing:
                continue
            for transcript in transcripts:
                complete_cds_len = (len(transcript.sequence) // 3) * 3
                for codon_start in range(0, complete_cds_len, 3):
                    codon_positions = [
                        transcript.cds_index_to_genomic.get(codon_start + offset)
                        for offset in range(3)
                    ]
                    missing_count = sum(1 for pos in codon_positions if pos in chrom_missing)
                    if missing_count >= threshold:
                        aa_pos = (codon_start // 3) + 1
                        aa_missing.add((chrom, transcript.gene_id, aa_pos))

        return [[chrom, f"{gene_id}.{aa_pos}"] for chrom, gene_id, aa_pos in sorted(aa_missing)]

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
) -> Dict[str, Dict[str, Set[str]]]:
    sample_ids: List[str] | None = None
    allele_map: Dict[str, Dict[str, Set[str]]] = {}
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
                            "genic_snps": set(),
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
            annotations_by_alt: Dict[int, Tuple[str, Set[str], Set[str]]] = {}
            for alt_idx, alt in enumerate(alt_alleles):
                if not is_simple_snp(ref, alt):
                    continue
                snp_id = f"{chrom}.{pos}.{ref.upper()}.{alt.upper()}"
                genic_ids: Set[str] = set()
                nonsynonymous_ids: Set[str] = set()
                if annotator is not None and annotator.is_genic(chrom, pos_int):
                    genic_ids.add(snp_id)
                    for gene_id, aa_pos, ref_aa, alt_aa in annotator.nonsynonymous_aa_changes(chrom, pos_int, ref, alt):
                        nonsynonymous_ids.add(f"{chrom}.{gene_id}.{aa_pos}.{ref_aa}.{alt_aa}")
                annotations_by_alt[alt_idx] = (snp_id, genic_ids, nonsynonymous_ids)

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
                    snp_id, genic_ids, nonsynonymous_ids = annotations
                    allele_map[sample_id]["alleles"].add(snp_id)
                    allele_map[sample_id]["genic_snps"].update(genic_ids)
                    allele_map[sample_id]["nonsynonymous"].update(nonsynonymous_ids)

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

def missing_positions_from_bed(depth_path: Path, min_depth: int) -> List[List]:
    missing: List[List] = []
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
                missing.append([chrom, pos])
    return missing

def genic_missing_positions(missing: List[List], annotator: VariantAnnotator | None) -> List[List]:
    if annotator is None:
        return []
    genic_missing: List[List] = []
    for chrom, pos in missing:
        try:
            pos_int = int(pos)
        except ValueError:
            continue
        if annotator.is_genic(str(chrom), pos_int):
            genic_missing.append([chrom, pos])
    return genic_missing

def aa_missing_positions(
    missing: List[List],
    annotator: VariantAnnotator | None,
    threshold: int,
) -> List[List]:
    if annotator is None:
        return []
    return annotator.aa_missing_positions(missing, threshold)

def active_allele_outputs(annotator: VariantAnnotator | None) -> Tuple[Tuple[str, str], ...]:
    if annotator is None:
        return ALLELE_OUTPUTS[:1]
    return ALLELE_OUTPUTS

def sample_output_path(outdir: Path, sample_id: str, suffix: str) -> Path:
    return outdir / f"{sample_id}_{suffix}.json"

def missing_sample_outputs(
    outdir: Path,
    sample_ids: List[str],
    allele_outputs: Tuple[Tuple[str, str], ...],
) -> Dict[str, Set[str]]:
    missing: Dict[str, Set[str]] = {}
    for sample_id in sample_ids:
        for allele_key, suffix in allele_outputs:
            if not sample_output_path(outdir, sample_id, suffix).exists():
                missing.setdefault(sample_id, set()).add(allele_key)
    return missing

def write_sample_payload(out_path: Path, alleles: Set[str], missing: List[List]) -> None:
    payload = {
        "alleles": sorted(alleles),
        "missing": missing,
    }
    try:
        with open(out_path, "x") as outfh:
            json.dump(payload, outfh, indent=2)
    except FileExistsError as exc:
        raise RuntimeError(f"Refusing to overwrite existing output file: {out_path}") from exc
    except Exception as exc:
        raise RuntimeError(
            f"Failed to write {out_path} "
            f"({len(payload['alleles'])} alleles, {len(payload['missing'])} missing positions): "
            f"{type(exc).__name__}: {exc}"
        ) from exc

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
    ap.add_argument("--aa-missing-codon-threshold", type=int, choices=(1, 2, 3), default=1,
                    help="For non-synonymous AA outputs, mark an AA site missing when at least this many codon nucleotides are missing (default: 1)")
    ap.add_argument("--ref-fasta", type=Path, default=None,
                    help="Reference FASTA used with --gff to annotate genic and non-synonymous SNPs")
    ap.add_argument("--gff", type=Path, default=None,
                    help="Reference GFF/GFF3 used with --ref-fasta to annotate genic and non-synonymous SNPs")
    ap.add_argument("--gene-id-attribute", type=str, default="gene_id",
                    help="GFF attribute to use as gene ID in non-synonymous allele IDs (default: gene_id)")
    args = ap.parse_args()

    if (args.ref_fasta is None) != (args.gff is None):
        raise SystemExit("--ref-fasta and --gff must be provided together")

    args.outdir.mkdir(parents=True, exist_ok=True)
    annotator = None
    if args.ref_fasta is not None and args.gff is not None:
        annotator = VariantAnnotator(args.ref_fasta, args.gff, args.gene_id_attribute)

    inputs = pair_files(args.vcf_in, args.depth_in)
    total = len(inputs)
    allele_outputs = active_allele_outputs(annotator)
    for idx, (vcf, depth) in enumerate(inputs, start=1):
        try:
            sample_ids = sample_ids_from_vcf_header(
                vcf,
                sample_regex=args.sample_regex,
                sample_regex_group=args.sample_regex_group,
            )
        except READ_EXCEPTIONS as exc:
            print(f"Warning: skipping unreadable VCF {vcf}: {exc}", file=sys.stderr)
            continue

        outputs_to_write = missing_sample_outputs(args.outdir, sample_ids, allele_outputs)
        missing_file_count = sum(len(output_keys) for output_keys in outputs_to_write.values())
        expected_file_count = len(sample_ids) * len(allele_outputs)
        if not outputs_to_write:
            print(f"\r {idx}/{total} : {vcf.name} ({len(sample_ids)} samples; outputs exist, skipping)", end="             ", flush=True)
            continue

        try:
            allele_map = sample_alleles_from_vcf(
                vcf_path=vcf,
                min_qual=args.min_qual,
                allow_filtered=args.allow_filtered,
                accept_missing_qual=args.accept_missing_qual,
                min_dp=args.min_dp,
                sample_regex=args.sample_regex,
                sample_regex_group=args.sample_regex_group,
                annotator=annotator,
            )
        except READ_EXCEPTIONS as exc:
            print(f"Warning: skipping unreadable VCF {vcf}: {exc}", file=sys.stderr)
            continue
        print(
            f"\r {idx}/{total} : {vcf.name} "
            f"({len(allele_map)} samples; writing {missing_file_count}/{expected_file_count} outputs)",
            end="             ",
            flush=True,
        )

        if depth is None:
            print(f"Warning: no depth BED file found for {vcf.name}; missing positions will be empty", file=sys.stderr)
            missing: List[List] = []
        else:
            try:
                missing = missing_positions_from_bed(depth, args.min_missing_depth)
            except READ_EXCEPTIONS as exc:
                print(f"Warning: skipping unreadable depth BED {depth}: {exc}", file=sys.stderr)
                missing = []

        missing_by_output = {
            "alleles": missing,
            "genic_snps": genic_missing_positions(missing, annotator),
            "nonsynonymous": aa_missing_positions(missing, annotator, args.aa_missing_codon_threshold),
        }

        for sample_id, allele_sets in allele_map.items():
            sample_missing_outputs = outputs_to_write.get(sample_id)
            if not sample_missing_outputs:
                continue
            for allele_key, suffix in allele_outputs:
                if allele_key not in sample_missing_outputs:
                    continue
                write_sample_payload(
                    sample_output_path(args.outdir, sample_id, suffix),
                    allele_sets[allele_key],
                    missing_by_output[allele_key],
                )

if __name__ == "__main__":
    main()
