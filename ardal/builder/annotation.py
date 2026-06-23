"""
Reference annotation helpers for Ardal matrix building.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from .vcf import DNA4, open_maybe_gzip

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

class CdsFeature:
    chrom: str
    start: int
    end: int
    strand: str
    parent: str
    gene_id: str

class TranscriptCds:
    gene_id: str
    chrom: str
    strand: str
    sequence: str
    genomic_to_cds_index: Dict[int, int]
    cds_index_to_genomic: Dict[int, int]

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
