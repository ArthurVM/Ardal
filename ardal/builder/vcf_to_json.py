"""
Compatibility entrypoint for VCF-to-sample-JSON building.
"""

import argparse
from pathlib import Path

from .api import vcf_to_sample_json


def parse_args():
    """Parse command-line arguments for the VCF-to-JSON workflow."""
    ap = argparse.ArgumentParser(
        description="Collect structured SNP calls and missing coverage blocks into JSON payloads."
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
    ap.add_argument("--overwrite", action="store_true",
                    help="Overwrite existing per-sample JSON outputs instead of skipping them")
    ap.add_argument("--sample-regex", type=str, default=None,
                    help="Optional regex to transform sample IDs from the VCF header")
    ap.add_argument("--sample-regex-group", type=int, default=1,
                    help="Capture group index to take from --sample-regex (default: 1)")
    ap.add_argument("--min-missing-depth", type=int, default=10,
                    help="Contiguous positions with coverage < this value are written as missing blocks (default: 10)")
    ap.add_argument("--aa-missing-codon-threshold", type=int, choices=(1, 2, 3), default=1,
                    help="Deprecated; retained for CLI compatibility and ignored by schema v2 output")
    ap.add_argument("--ref-fasta", type=Path, default=None,
                    help="Reference FASTA used with --gff to annotate genic and non-synonymous SNPs")
    ap.add_argument("--gff", type=Path, default=None,
                    help="Reference GFF/GFF3 used with --ref-fasta to annotate genic and non-synonymous SNPs")
    ap.add_argument("--gene-id-attribute", type=str, default="gene_id",
                    help="GFF attribute to use as gene ID in non-synonymous allele IDs (default: gene_id)")
    
    return ap.parse_args()


def main():
    """Entry point for invoking the VCF-to-JSON workflow."""
    args = parse_args()
    try:
        vcf_to_sample_json(
            vcf_in=args.vcf_in,
            depth_in=args.depth_in,
            outdir=args.outdir,
            min_qual=args.min_qual,
            allow_filtered=args.allow_filtered,
            accept_missing_qual=args.accept_missing_qual,
            min_dp=args.min_dp,
            overwrite=args.overwrite,
            sample_regex=args.sample_regex,
            sample_regex_group=args.sample_regex_group,
            min_missing_depth=args.min_missing_depth,
            ref_fasta=args.ref_fasta,
            gff=args.gff,
            gene_id_attribute=args.gene_id_attribute,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc


__all__ = [
    "main",
    "parse_args",
    "vcf_to_sample_json",
]
