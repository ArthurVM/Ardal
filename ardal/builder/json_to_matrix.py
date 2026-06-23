"""
Compatibility entrypoint for JSON-to-matrix building.
"""

import argparse
import sys
from pathlib import Path

from .api import create_all_outputs, json_directory_to_matrices
from .schemas import DATA_TYPES, DEFAULT_ALLELE_ID_FORMAT, DEFAULT_NONSYN_ALLELE_ID_FORMAT


def parse_args():
    """Parse command-line arguments for the JSON-to-matrix workflow."""
    ap = argparse.ArgumentParser(
        description="Convert per-sample allele JSON files into Ardal dense/bitpacked matrices."
    )
    ap.add_argument("json_directory", type=Path, help="Directory containing sample JSON payloads")
    ap.add_argument("output_prefix", type=str, help="Prefix/path for emitted artifacts")
    ap.add_argument("--sha256", action="store_true", dest="compute_sha256",
                    help="Compute sha256 for .bin matrices")
    ap.add_argument("--npy", action="store_true", dest="emit_npy",
                    help="Keep dense .npy matrices (memmap files) as outputs")
    ap.add_argument("--npz", action="store_true", dest="emit_npz",
                    help="Write compressed dense .npz archives")
    ap.add_argument("--csv", action="store_true", dest="emit_csv",
                    help="Write dense matrices as wide CSV files")
    ap.add_argument("--agnostic", action="store_true",
                    help="Also generate reference-agnostic/minor allele matrices")
    ap.add_argument("--no-agnostic", action="store_true",
                    help="Deprecated compatibility flag; agnostic matrices are skipped unless --agnostic is set")
    ap.add_argument("--allele-id-format", type=str, default=None,
                    help=(
                        f"Allele ID format string (default: {DEFAULT_ALLELE_ID_FORMAT}; "
                        "{chr} optional for single-chrom references; {gene} supported for gene-relative IDs)"
                    ))
    ap.add_argument("--nonsyn-allele-id-format", type=str, default=DEFAULT_NONSYN_ALLELE_ID_FORMAT,
                    help=(
                        f"Allele ID format for nonsyn matrices "
                        f"(default: {DEFAULT_NONSYN_ALLELE_ID_FORMAT})"
                    ))
    ap.add_argument("--data-types", nargs="+", choices=DATA_TYPES, default=["snps"],
                    help="Data projections to build from schema-v2 JSONs (default: snps)")
    ap.add_argument("--indent", type=int, dest="json_indent",
                    help="Indentation level for JSON artifacts (default uses compact JSON)")
    
    return ap.parse_args()


def main():
    """Entry point for invoking the JSON-to-matrix conversion workflow."""
    args = parse_args()
    try:
        json_directory_to_matrices(
            args.json_directory,
            args.output_prefix,
            compute_sha256=args.compute_sha256,
            emit_npy=args.emit_npy,
            emit_npz=args.emit_npz,
            emit_csv=args.emit_csv,
            emit_agnostic=args.agnostic and not args.no_agnostic,
            allele_id_format=args.allele_id_format,
            nonsyn_allele_id_format=args.nonsyn_allele_id_format,
            data_types=args.data_types,
            json_indent=args.json_indent,
        )
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)


__all__ = [
    "create_all_outputs",
    "json_directory_to_matrices",
    "main",
    "parse_args",
]
