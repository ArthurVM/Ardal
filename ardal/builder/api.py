"""
High-level Python API for Ardal matrix-building workflows.
"""

import sys
from pathlib import Path
from typing import Sequence

from .alleles import configure_allele_id_format, reorder_alleles_by_position
from .annotation import VariantAnnotator
from .matrix_plan import allele_genomic_sites_from_json, derive_site_statistics, index_samples_and_alleles
from .missingness import build_columns_by_genomic_site, column_masks_from_missing_blocks, expand_index_ranges, missing_blocks_from_bed
from .sample_json import count_v2_annotation_entries, missing_sample_outputs, sample_output_path, write_sample_payload
from .schemas import DATA_TYPES, DEFAULT_NONSYN_ALLELE_ID_FORMAT, build_allele_model, column_masks_payload
from .vcf import READ_EXCEPTIONS, pair_files, sample_alleles_from_vcf, sample_ids_from_vcf_header, samples_with_candidate_alleles_from_vcf
from .writers import apply_reference_minor_columns, build_dense_matrices, dense_work_path, emit_bitpack_artifacts, emit_direct_bitpack_artifacts, write_dense_artifacts, write_site_metadata


def create_all_outputs(
    json_files,
    output_prefix,
    compute_sha256 = False,
    emit_npy = False,
    emit_npz = False,
    emit_csv = False,
    emit_agnostic = True,
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
    data_type : str = "snps",
    json_indent : int | None = None,
):
    """
    Description:
        Orchestrate the full pipeline from sample JSON payloads to dense/bitpacked matrices and metadata outputs.

    Inputs:
        json_files: Iterable of JSON file paths.
        output_prefix: Prefix used for every emitted artifact group.
        compute_sha256 / emit_npy / emit_npz / emit_csv / emit_agnostic: Flags controlling serialization steps.
        json_indent: Optional indentation for JSON outputs.

    Outputs:
        None; writes files to disk.

    Exceptions:
        Propagates exceptions from subordinate helper routines.
    """
    if data_type == "nonsyns":
        effective_allele_id_format = nonsyn_allele_id_format or DEFAULT_NONSYN_ALLELE_ID_FORMAT
    else:
        effective_allele_id_format = allele_id_format
    configure_allele_id_format(effective_allele_id_format)
    (
        sample_to_idx,
        allele_to_idx,
        site_reference,
        site_coordinates,
        site_alt_counts,
        site_alt_samples,
        sample_missing_positions,
    ) = index_samples_and_alleles(
        json_files,
        data_type=data_type,
        allele_id_format=allele_id_format,
        nonsyn_allele_id_format=nonsyn_allele_id_format,
        collect_site_alt_samples=emit_agnostic,
    )

    allele_to_idx, ordered_alleles = reorder_alleles_by_position(allele_to_idx)

    n_rows, n_ref_cols = len(sample_to_idx), len(allele_to_idx)
    site_metadata, minor_columns, minor_alt_to_idx, minor_ref_to_idx = derive_site_statistics(
        sample_to_idx,
        allele_to_idx,
        site_reference,
        site_coordinates,
        site_alt_counts,
        site_alt_samples,
    )
    n_minor_cols = len(minor_columns)

    guids = [None] * n_rows
    
    for s, i in sample_to_idx.items():
        guids[i] = s
    alleles = ordered_alleles[:]
    ref_allele_model = build_allele_model(
        data_type=data_type,
        allele_id_format=allele_id_format,
        nonsyn_allele_id_format=nonsyn_allele_id_format,
        matrix_projection="allele_presence",
    )
    agnostic_allele_model = build_allele_model(
        data_type=data_type,
        allele_id_format=allele_id_format,
        nonsyn_allele_id_format=nonsyn_allele_id_format,
        matrix_projection="reference_agnostic_minor",
    )

    ref_prefix = f"{output_prefix}.ref"
    minor_prefix = f"{output_prefix}.agnostic"

    emit_dense_artifacts = emit_npy or emit_npz or emit_csv
    if ( not emit_dense_artifacts and not emit_agnostic ):
        allele_genomic_sites = allele_genomic_sites_from_json(
            json_files,
            data_type,
            allele_id_format,
            nonsyn_allele_id_format,
        )
        columns_by_site = build_columns_by_genomic_site(
            ordered_alleles,
            allele_to_idx,
            allele_genomic_sites,
        )
        emit_direct_bitpack_artifacts(
            json_files,
            sample_to_idx,
            allele_to_idx,
            sample_missing_positions,
            columns_by_site,
            ref_prefix,
            guids,
            alleles,
            compute_sha256,
            json_indent,
            allele_model=ref_allele_model,
            data_type=data_type,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
        )
        write_site_metadata(
            output_prefix,
            n_rows,
            minor_columns,
            site_metadata,
            ref_prefix,
            None,
            json_indent=json_indent,
        )
        
        return

    allele_genomic_sites = allele_genomic_sites_from_json(
        json_files,
        data_type,
        allele_id_format,
        nonsyn_allele_id_format,
    )
    columns_by_site = build_columns_by_genomic_site(
        ordered_alleles,
        allele_to_idx,
        allele_genomic_sites,
    )
    column_masks = column_masks_from_missing_blocks(
        sample_to_idx,
        sample_missing_positions,
        columns_by_site,
    )
    
    missing_sites_payload = column_masks_payload(
        data_type=data_type,
        column_masks=column_masks,
    )

    sample_masked_columns_by_idx = [
        expand_index_ranges(column_masks.get(guid, [])) for guid in guids
    ]

    ref_npy_path = dense_work_path(ref_prefix, emit_npy)
    minor_npy_path = dense_work_path(minor_prefix, emit_npy)

    dense_ref, dense_minor = build_dense_matrices(
        json_files,
        sample_to_idx,
        allele_to_idx,
        ref_npy_path,
        minor_npy_path,
        minor_columns,
        minor_alt_to_idx,
        sample_masked_columns_by_idx,
        emit_agnostic,
        data_type=data_type,
        allele_id_format=allele_id_format,
        nonsyn_allele_id_format=nonsyn_allele_id_format,
    )

    if ( dense_minor is not None ):
        apply_reference_minor_columns(dense_minor, minor_ref_to_idx, site_alt_samples)

    if ( emit_dense_artifacts ):
        write_dense_artifacts(
            dense_ref,
            guids,
            alleles,
            ref_prefix,
            emit_npy,
            emit_npz,
            emit_csv,
            missing_sites=missing_sites_payload,
            allele_model=ref_allele_model,
            json_indent=json_indent,
        )
        if ( dense_minor is not None ):
            write_dense_artifacts(
                dense_minor,
                guids,
                minor_columns,
                minor_prefix,
                emit_npy,
                emit_npz,
                emit_csv,
                missing_sites=missing_sites_payload,
                allele_model=agnostic_allele_model,
                json_indent=json_indent,
            )

    ## bitpack from dense (reference)
    emit_bitpack_artifacts(
        dense_ref,
        n_ref_cols,
        ref_prefix,
        guids,
        alleles,
        n_rows,
        missing_sites_payload,
        compute_sha256,
        json_indent,
        keep_dense_file=emit_npy,
        npy_path=ref_npy_path,
        allele_model=ref_allele_model,
    )

    ## bitpack from dense (reference-agnostic)
    emit_bitpack_artifacts(
        dense_minor,
        n_minor_cols,
        minor_prefix,
        guids,
        minor_columns,
        n_rows,
        missing_sites_payload,
        compute_sha256,
        json_indent,
        keep_dense_file=emit_npy,
        npy_path=minor_npy_path,
        allele_model=agnostic_allele_model,
    )

    ## site metadata summarising major/minor alleles
    write_site_metadata(
        output_prefix,
        n_rows,
        minor_columns,
        site_metadata,
        ref_prefix,
        minor_prefix if emit_agnostic else None,
        json_indent=json_indent,
    )

def json_directory_to_matrices(
    json_directory : Path | str,
    output_prefix : str,
    *,
    compute_sha256 : bool = False,
    emit_npy : bool = False,
    emit_npz : bool = False,
    emit_csv : bool = False,
    emit_agnostic : bool = False,
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = DEFAULT_NONSYN_ALLELE_ID_FORMAT,
    data_types : Sequence[str] = ("snps",),
    json_indent : int | None = None,
):
    """
    Description:
        Convert a directory of per-sample JSON payloads into Ardal matrix artifacts.

    Inputs:
        json_directory: Directory containing sample JSON payloads.
        output_prefix: Prefix/path used for emitted matrix artifacts.
        compute_sha256 / emit_npy / emit_npz / emit_csv / emit_agnostic:
            Serialization options matching the compatibility script.
        allele_id_format / nonsyn_allele_id_format:
            Identifier formats used for nucleotide and amino-acid projections.
        data_types: One or more projections to build.
        json_indent: Optional indentation for emitted JSON metadata.

    Outputs:
        None; matrix artifacts are written to disk.

    Exceptions:
        Raises ValueError when inputs are invalid and propagates downstream IO/parsing errors.
    """
    json_dir = Path(json_directory)
    if ( not json_dir.is_dir() ):
        raise ValueError(f"{json_dir} is not a directory")

    json_files = sorted(p for p in json_dir.iterdir() if p.is_file() and p.suffix == ".json")
    if ( not json_files ):
        raise ValueError(f"No JSON files found in {json_dir}")

    for data_type in data_types:
        if ( data_type not in DATA_TYPES ):
            raise ValueError(
                f"Unsupported data type: {data_type}. Expected one of: {', '.join(DATA_TYPES)}"
            )

    if ( emit_npy or emit_npz ):
        requested = []
        if ( emit_npy ):
            requested.append(".npy")
        if ( emit_npz ):
            requested.append(".npz")
        print(
            "[Warning] "
            + " and ".join(requested)
            + " dense outputs are redundant with the default .bin v2 output for Ardal loading. "
            + "They are still written because you requested them, but .bin v2 should be preferred.",
            file=sys.stderr,
        )

    requested_types = set(data_types)
    if ( requested_types.intersection({"genic_snps", "nonsyns"}) ):
        payload_count, genic_count, nonsyn_count = count_v2_annotation_entries(json_files)
        if ( payload_count ):
            if ( "genic_snps" in requested_types and genic_count == 0 ):
                raise ValueError(
                    "Requested genic_snps, but no genic annotations were found in the schema-v2 JSON files. "
                    "These JSONs were likely created without a reference FASTA and GFF pair; refusing to build an empty genic matrix."
                )
            if ( "nonsyns" in requested_types and nonsyn_count == 0 ):
                raise ValueError(
                    "Requested nonsyns, but no nonsynonymous annotations were found in the schema-v2 JSON files. "
                    "These JSONs were likely created without a reference FASTA and GFF pair; refusing to build an empty nonsynonymous matrix."
                )

    for data_type in data_types:
        typed_prefix = f"{output_prefix}.{data_type}"
        print(f"[Projection] {data_type} -> {typed_prefix}")
        create_all_outputs(
            json_files,
            typed_prefix,
            compute_sha256=compute_sha256,
            emit_npy=emit_npy,
            emit_npz=emit_npz,
            emit_csv=emit_csv,
            emit_agnostic=emit_agnostic,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
            data_type=data_type,
            json_indent=json_indent,
        )

def vcf_to_sample_json(
    *,
    vcf_in : Path | str,
    depth_in : Path | str | None = None,
    outdir : Path | str = Path("./"),
    min_qual : float = 30.0,
    allow_filtered : bool = False,
    accept_missing_qual : bool = False,
    min_dp : int | None = None,
    overwrite : bool = False,
    sample_regex : str | None = None,
    sample_regex_group : int = 1,
    min_missing_depth : int = 10,
    ref_fasta : Path | str | None = None,
    gff : Path | str | None = None,
    gene_id_attribute : str = "gene_id",
) -> None:
    """
    Description:
        Convert VCF and optional depth BED inputs into per-sample Ardal JSON payloads.

    Inputs:
        vcf_in: VCF file or directory containing VCF/VCF.GZ files.
        depth_in: Optional depth BED/BED.GZ file or directory.
        outdir: Destination directory for per-sample JSON payloads.
        min_qual / allow_filtered / accept_missing_qual / min_dp:
            Variant-level filtering controls.
        overwrite: Whether existing JSON outputs may be replaced.
        sample_regex / sample_regex_group: Optional sample ID rewrite rule.
        min_missing_depth: Coverage threshold below which depth positions are missing.
        ref_fasta / gff / gene_id_attribute:
            Optional reference annotation inputs for genic and nonsynonymous projections.

    Outputs:
        None; per-sample JSON files are written to `outdir`.

    Exceptions:
        Raises ValueError when reference annotation inputs are incomplete.
        Propagates RuntimeError from failed output writes.
    """
    vcf_in = Path(vcf_in)
    depth_in = Path(depth_in) if depth_in is not None else None
    outdir = Path(outdir)
    ref_fasta = Path(ref_fasta) if ref_fasta is not None else None
    gff = Path(gff) if gff is not None else None

    if ( ref_fasta is None ) != ( gff is None ):
        raise ValueError("--ref-fasta and --gff must be provided together")

    outdir.mkdir(parents=True, exist_ok=True)
    annotator = None
    if ref_fasta is not None and gff is not None:
        annotator = VariantAnnotator(ref_fasta, gff, gene_id_attribute)

    inputs = pair_files(vcf_in, depth_in)
    total = len(inputs)
    for idx, (vcf, depth) in enumerate(inputs, start=1):
        try:
            sample_ids = sample_ids_from_vcf_header(
                vcf,
                sample_regex=sample_regex,
                sample_regex_group=sample_regex_group,
            )
        except READ_EXCEPTIONS as exc:
            print(f"Warning: skipping unreadable VCF {vcf}: {exc}", file=sys.stderr)
            continue

        outputs_to_write = missing_sample_outputs(outdir, sample_ids, overwrite)
        expected_file_count = len(sample_ids)
        if not outputs_to_write:
            print(f"\r {idx}/{total} : {vcf.name} ({len(sample_ids)} samples; outputs exist, skipping)", end="             ", flush=True)
            continue

        try:
            candidate_samples = samples_with_candidate_alleles_from_vcf(
                vcf_path=vcf,
                min_qual=min_qual,
                allow_filtered=allow_filtered,
                accept_missing_qual=accept_missing_qual,
                min_dp=min_dp,
                sample_regex=sample_regex,
                sample_regex_group=sample_regex_group,
            )
        except READ_EXCEPTIONS as exc:
            print(f"Warning: skipping unreadable VCF {vcf}: {exc}", file=sys.stderr)
            continue

        outputs_to_write = {
            sample_id
            for sample_id in outputs_to_write
            if sample_id in candidate_samples
        }
        if not outputs_to_write:
            print(
                f"\r {idx}/{total} : {vcf.name} "
                f"({len(sample_ids)} samples; no passing SNP alleles, skipping)",
                end="             ",
                flush=True,
            )
            continue

        try:
            allele_map = sample_alleles_from_vcf(
                vcf_path=vcf,
                min_qual=min_qual,
                allow_filtered=allow_filtered,
                accept_missing_qual=accept_missing_qual,
                min_dp=min_dp,
                sample_regex=sample_regex,
                sample_regex_group=sample_regex_group,
                annotator=annotator,
            )
        except READ_EXCEPTIONS as exc:
            print(f"Warning: skipping unreadable VCF {vcf}: {exc}", file=sys.stderr)
            continue

        outputs_to_write = {
            sample_id
            for sample_id in outputs_to_write
            if allele_map.get(sample_id, {}).get("alleles")
        }
        missing_file_count = len(outputs_to_write)
        if not outputs_to_write:
            print(
                f"\r {idx}/{total} : {vcf.name} "
                f"({len(allele_map)} samples; no alleles to write, skipping)",
                end="             ",
                flush=True,
            )
            continue

        print(
            f"\r {idx}/{total} : {vcf.name} "
            f"({len(allele_map)} samples; writing {missing_file_count}/{expected_file_count} outputs)",
            end="             ",
            flush=True,
        )

        if depth is None:
            print(f"Warning: no depth BED file found for {vcf.name}; missing blocks will be empty", file=sys.stderr, end="")
            missing: List[List] = []
        else:
            try:
                missing = missing_blocks_from_bed(depth, min_missing_depth)
            except READ_EXCEPTIONS as exc:
                print(f"Warning: skipping unreadable depth BED {depth}: {exc}", file=sys.stderr)
                missing = []

        for sample_id, allele_sets in allele_map.items():
            if sample_id not in outputs_to_write:
                continue
            write_sample_payload(
                sample_output_path(outdir, sample_id),
                sample_id,
                allele_sets,
                missing,
                overwrite,
            )
