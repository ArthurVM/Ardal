"""
Shared schemas for Ardal matrix-building artifacts.
"""

from __future__ import annotations

from typing import Dict


DEFAULT_ALLELE_ID_FORMAT = "{chr}.{pos}.{ref}.{alt}"
DEFAULT_NONSYN_ALLELE_ID_FORMAT = "{chr}.{gene}.{pos}.{ref}.{alt}"
SAMPLE_SCHEMA_VERSION = "ardal.sample_variants.v2"
V2_SAMPLE_SCHEMA_VERSION = SAMPLE_SCHEMA_VERSION
DATA_TYPES = ("snps", "genic_snps", "nonsyns")

ALLELE_RECORD_FIELDS = ["chrom", "start", "end", "ref", "alt"]
NONSYNONYMOUS_RECORD_FIELDS = ["allele_index", "gene_id", "aa_pos", "ref_aa", "alt_aa"]
MISSING_BLOCK_FIELDS = ["chrom", "start", "end"]


def build_matrix_metadata(
    *,
    format_name : str,
    dtype : str,
    n_rows : int,
    n_cols : int,
    matrix_file : str,
    row_major : bool = True,
    data_nbytes : int | None = None,
    data_sha256 : str | None = None,
    words_per_row : int | None = None,
    bits_per_word : int | None = None,
    row_stride_bytes : int | None = None,
    endianness : str | None = None,
    matrix_kind : str | None = None,
    allele_model : Dict[str, object] | None = None,
    sections : Dict[str, object] | None = None,
) -> Dict[str, object]:
    """
    Description:
        Build the shared `.meta` sidecar metadata block for serialized matrices.

    Inputs:
        format_name / dtype / n_rows / n_cols / matrix_file:
            Core matrix identity and shape fields.
        row_major / data_nbytes / data_sha256 / words_per_row / bits_per_word /
        row_stride_bytes / endianness / matrix_kind / allele_model / sections:
            Optional serialization details used by matrix loaders.

    Outputs:
        Returns the metadata dictionary ready to place under the `meta` key.

    Exceptions:
        None.
    """
    def _as_int( value ):
        
        return None if value is None else int(value)

    meta = {
        "format": format_name,
        "dtype": dtype,
        "endianness": endianness,
        "row_major": bool(row_major),
        "n_rows": int(n_rows),
        "n_cols": int(n_cols),
        "matrix_file": matrix_file,
        "data_nbytes": _as_int(data_nbytes),
        "data_sha256": data_sha256,
        "words_per_row": _as_int(words_per_row),
        "bits_per_word": _as_int(bits_per_word),
        "row_stride_bytes": _as_int(row_stride_bytes),
        "generated_by": "ardal.builder",
    }
    if ( matrix_kind is not None ):
        meta["matrix_kind"] = matrix_kind
    if ( allele_model is not None ):
        meta["allele_model"] = allele_model
    if ( sections is not None ):
        meta["sections"] = sections

    return meta


def build_allele_model(
    *,
    data_type : str,
    allele_id_format : str | None,
    nonsyn_allele_id_format : str | None,
    matrix_projection : str | None = None,
) -> Dict[str, object]:
    """
    Description:
        Describe the biological alphabet and coordinate system encoded by matrix alleles.

    Inputs:
        data_type: Matrix projection being generated from sample JSON payloads.
        allele_id_format / nonsyn_allele_id_format:
            Identifier formats used for nucleotide and amino-acid projections.
        matrix_projection: Optional matrix-level projection label.

    Outputs:
        Returns an allele-model dictionary suitable for matrix sidecars.

    Exceptions:
        None.
    """
    if ( data_type == "nonsyns" ):
        model = {
            "alphabet": "amino_acid",
            "coordinate_system": "protein",
            "variant_type": "residue_substitution",
            "allele_id_format": nonsyn_allele_id_format or DEFAULT_NONSYN_ALLELE_ID_FORMAT,
            "data_projection": data_type,
            "compatible_metrics": ["hamming", "jaccard", "cosine"],
        }
    else:
        model = {
            "alphabet": "nucleotide",
            "coordinate_system": "genomic",
            "variant_type": "snv",
            "allele_id_format": allele_id_format or DEFAULT_ALLELE_ID_FORMAT,
            "data_projection": data_type,
            "compatible_metrics": ["hamming", "jaccard", "cosine", "snv"],
        }

    if ( matrix_projection is not None ):
        model["matrix_projection"] = matrix_projection

    return model


def column_masks_payload(
    *,
    data_type : str,
    column_masks : Dict[str, object],
) -> Dict[str, object]:
    """
    Description:
        Build the legacy JSON column-mask missingness payload.

    Inputs:
        data_type: Matrix projection the masks apply to.
        column_masks: Per-sample compact column ranges.

    Outputs:
        Returns the sidecar missingness payload.

    Exceptions:
        None.
    """
    return {
        "missing": {
            "schema": "ardal.column_masks.v1",
            "data_type": data_type,
            "encoding": "column_index_ranges",
            "coordinate_system": "matrix_columns_0_based_inclusive",
            "column_masks": column_masks,
        }
    }


def binary_missing_ranges_payload() -> Dict[str, object]:
    """
    Description:
        Build the `.bin v2` binary-section missingness payload.

    Inputs:
        None.

    Outputs:
        Returns the sidecar missingness payload pointing at binary sections.

    Exceptions:
        None.
    """
    return {
        "missing": {
            "schema": "ardal.missing_column_ranges.v1",
            "encoding": "binary_sections",
            "coordinate_system": "matrix_columns_0_based",
            "range_semantics": "start_inclusive_end_exclusive",
            "offsets_section": "missing_offsets",
            "ranges_section": "missing_ranges",
        }
    }
