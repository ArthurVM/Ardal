"""
Ardal sample JSON v2 reading, writing, and projection helpers.
"""

import json
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Set, Tuple

from .alleles import allele_id_from_fields, parse_allele_key
from .missingness import missing_blocks_from_entries, normalize_missing_entries
from .schemas import (
    ALLELE_RECORD_FIELDS,
    DEFAULT_NONSYN_ALLELE_ID_FORMAT,
    MISSING_BLOCK_FIELDS,
    NONSYNONYMOUS_RECORD_FIELDS,
    SAMPLE_SCHEMA_VERSION,
    V2_SAMPLE_SCHEMA_VERSION,
)

def derive_sample_id_from_path(
    json_path : Path,
) -> str:
    """Infer the sample identifier from a JSON filename."""
    stem = json_path.stem
    if ( stem.endswith("_snps") ):
        stem = stem[:-5]

    return stem

def format_v2_alleles(
    payload : Dict,
    data_type : str,
    allele_id_format : str | None,
    nonsyn_allele_id_format : str | None,
) -> List[str]:
    """Project a schema-v2 sample payload to the requested matrix allele IDs."""
    raw_alleles = payload.get("alleles") or []
    formatted_snps: List[str] = []
    for record in raw_alleles:
        if not isinstance(record, (list, tuple)) or len(record) != 5:
            continue
        chrom, start, end, ref, alt = record
        formatted_snps.append(
            allele_id_from_fields(
                chrom=chrom,
                start=start,
                end=end,
                ref=ref,
                alt=alt,
                allele_id_format=allele_id_format,
            )
        )

    if ( data_type == "snps" ):
        return formatted_snps

    if ( data_type == "genic_snps" ):
        alleles: List[str] = []
        for idx in payload.get("genic") or []:
            try:
                alleles.append(formatted_snps[int(idx)])
            except (TypeError, ValueError, IndexError):
                continue
        return alleles

    if ( data_type == "nonsyns" ):
        fmt = nonsyn_allele_id_format or DEFAULT_NONSYN_ALLELE_ID_FORMAT
        alleles: List[str] = []
        for record in payload.get("nonsynonymous") or []:
            if not isinstance(record, (list, tuple)) or len(record) != 5:
                continue
            allele_idx, gene_id, aa_pos, ref_aa, alt_aa = record
            try:
                allele_record = raw_alleles[int(allele_idx)]
            except (TypeError, ValueError, IndexError):
                continue
            if not isinstance(allele_record, (list, tuple)) or len(allele_record) != 5:
                continue
            chrom, _start, _end, _ref, _alt = allele_record
            alleles.append(
                allele_id_from_fields(
                    chrom=chrom,
                    start=aa_pos,
                    end=aa_pos,
                    ref=ref_aa,
                    alt=alt_aa,
                    gene=gene_id,
                    allele_id_format=fmt,
                )
            )
        return alleles

    raise ValueError(f"Unsupported data type: {data_type}")

def iter_json_payloads(json_path: Path) -> Iterator[Dict]:
    """Yield one or more sample payload dictionaries from a JSON file."""
    with open(json_path) as fh:
        payload = json.load(fh)

    if isinstance(payload, list):
        for item in payload:
            if isinstance(item, dict):
                yield item
        return

    if isinstance(payload, dict):
        yield payload
        return

    raise ValueError(f"Unsupported JSON structure in {json_path}")

def count_v2_annotation_entries(
    json_files : Iterable[Path],
) -> Tuple[int, int, int]:
    """
    Description:
        Count schema-v2 payloads and their genic / nonsynonymous annotation entries.

    Inputs:
        json_files: Iterable of sample JSON paths to inspect.

    Outputs:
        Returns `(payload_count, genic_entry_count, nonsynonymous_entry_count)`.

    Exceptions:
        Propagates JSON parsing errors from malformed payload files.
    """
    payload_count = 0
    genic_count = 0
    nonsyn_count = 0
    for json_path in json_files:
        for payload in iter_json_payloads(Path(json_path)):
            if ( payload.get("schema_version") != V2_SAMPLE_SCHEMA_VERSION ):
                continue
            payload_count += 1
            genic_count += len(payload.get("genic") or [])
            nonsyn_count += len(payload.get("nonsynonymous") or [])

    return payload_count, genic_count, nonsyn_count

def project_v2_entries(
    payload : Dict,
    data_type : str,
    allele_id_format : str | None,
    nonsyn_allele_id_format : str | None,
) -> List[Tuple[str, str, int]]:
    """Return (matrix_allele_id, genomic_chrom, genomic_pos) records for v2 payloads."""
    raw_alleles = payload.get("alleles") or []
    snp_entries: List[Tuple[str, str, int]] = []
    for record in raw_alleles:
        if not isinstance(record, (list, tuple)) or len(record) != 5:
            continue
        chrom, start, end, ref, alt = record
        try:
            genomic_pos = int(start)
        except (TypeError, ValueError):
            continue
        snp_entries.append(
            (
                allele_id_from_fields(
                    chrom=chrom,
                    start=start,
                    end=end,
                    ref=ref,
                    alt=alt,
                    allele_id_format=allele_id_format,
                ),
                str(chrom),
                genomic_pos,
            )
        )

    if ( data_type == "snps" ):
        return snp_entries

    if ( data_type == "genic_snps" ):
        entries: List[Tuple[str, str, int]] = []
        for idx in payload.get("genic") or []:
            try:
                entries.append(snp_entries[int(idx)])
            except (TypeError, ValueError, IndexError):
                continue
        return entries

    if ( data_type == "nonsyns" ):
        fmt = nonsyn_allele_id_format or DEFAULT_NONSYN_ALLELE_ID_FORMAT
        entries: List[Tuple[str, str, int]] = []
        for record in payload.get("nonsynonymous") or []:
            if not isinstance(record, (list, tuple)) or len(record) != 5:
                continue
            allele_idx, gene_id, aa_pos, ref_aa, alt_aa = record
            try:
                allele_record = raw_alleles[int(allele_idx)]
            except (TypeError, ValueError, IndexError):
                continue
            if not isinstance(allele_record, (list, tuple)) or len(allele_record) != 5:
                continue
            chrom, start, _end, _ref, _alt = allele_record
            try:
                genomic_pos = int(start)
            except (TypeError, ValueError):
                continue
            entries.append(
                (
                    allele_id_from_fields(
                        chrom=chrom,
                        start=aa_pos,
                        end=aa_pos,
                        ref=ref_aa,
                        alt=alt_aa,
                        gene=gene_id,
                        allele_id_format=fmt,
                    ),
                    str(chrom),
                    genomic_pos,
                )
            )
        return entries

    raise ValueError(f"Unsupported data type: {data_type}")

def iter_sample_entries(
    json_path : Path,
    data_type : str = "snps",
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
) -> Iterator[Tuple[str, List[str], List[str]]]:
    """
    Description:
        Iterate over sample payloads, yielding allele and missing-site entries per sample object.

    Inputs:
        json_path: Path to the sample JSON file on disk.

    Outputs:
        Generates `(sample_id, alleles, missing_positions)` tuples for consumers.

    Exceptions:
        Raises ValueError when the JSON layout is not supported.
    """
    for payload in iter_json_payloads(json_path):
        if ( payload.get("schema_version") == V2_SAMPLE_SCHEMA_VERSION ):
            sample_id = payload.get("sample_id") or derive_sample_id_from_path(json_path)
            alleles = format_v2_alleles(payload, data_type, allele_id_format, nonsyn_allele_id_format)
            missing = missing_blocks_from_entries(payload.get("missing") or [])
            yield sample_id, alleles, missing
            continue

        if ( "alleles" in payload or "missing" in payload ):
            sample_id = payload.get("sample_id") or derive_sample_id_from_path(json_path)
            alleles = list(payload.get("alleles") or [])
            missing = normalize_missing_entries(payload.get("missing") or [])
            yield sample_id, alleles, missing
            continue

        for sample_id, alleles in payload.items():
            yield sample_id, list(alleles or []), []

def iter_sample_projected_entries(
    json_path : Path,
    data_type : str = "snps",
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
) -> Iterator[Tuple[str, List[Tuple[str, str, int]], List]]:
    """Yield allele IDs with genomic sites used for missing-mask projection."""
    for payload in iter_json_payloads(json_path):
        if ( payload.get("schema_version") == V2_SAMPLE_SCHEMA_VERSION ):
            sample_id = payload.get("sample_id") or derive_sample_id_from_path(json_path)
            entries = project_v2_entries(payload, data_type, allele_id_format, nonsyn_allele_id_format)
            missing = missing_blocks_from_entries(payload.get("missing") or [])
            yield sample_id, entries, missing
            continue

        if ( "alleles" in payload or "missing" in payload ):
            sample_id = payload.get("sample_id") or derive_sample_id_from_path(json_path)
            allele_ids = list(payload.get("alleles") or [])
            missing = normalize_missing_entries(payload.get("missing") or [])
            payload_items = [(sample_id, allele_ids, missing)]
        else:
            payload_items = [
                (sample_id, list(alleles or []), [])
                for sample_id, alleles in payload.items()
            ]

        for sample_id, allele_ids, missing in payload_items:
            entries = []
            for allele_id in allele_ids:
                chrom, pos, _ref, _alt = parse_allele_key(allele_id)
                try:
                    entries.append((allele_id, str(chrom), int(pos)))
                except (TypeError, ValueError):
                    entries.append((allele_id, str(chrom), pos))
            yield sample_id, entries, missing

def sample_output_path(outdir: Path, sample_id: str) -> Path:
    return outdir / f"{sample_id}.json"

def missing_sample_outputs(outdir: Path, sample_ids: List[str], overwrite: bool) -> Set[str]:
    if overwrite:
        return set(sample_ids)
    return {
        sample_id
        for sample_id in sample_ids
        if not sample_output_path(outdir, sample_id).exists()
    }

def sample_payload(sample_id: str, allele_sets: Dict[str, Set], missing: List[List]) -> Dict:
    alleles = sorted(allele_sets["alleles"])
    allele_to_idx = {allele: idx for idx, allele in enumerate(alleles)}
    genic = sorted(allele_to_idx[allele] for allele in allele_sets["genic"] if allele in allele_to_idx)
    nonsynonymous = sorted(
        [
            allele_to_idx[allele],
            gene_id,
            aa_pos,
            ref_aa,
            alt_aa,
        ]
        for allele, gene_id, aa_pos, ref_aa, alt_aa in allele_sets["nonsynonymous"]
        if allele in allele_to_idx
    )
    payload = {
        "schema_version": SAMPLE_SCHEMA_VERSION,
        "sample_id": sample_id,
        "coordinate_system": "1-based inclusive",
        "allele_record_fields": ALLELE_RECORD_FIELDS,
        "nonsynonymous_record_fields": NONSYNONYMOUS_RECORD_FIELDS,
        "missing_block_fields": MISSING_BLOCK_FIELDS,
        "alleles": [list(allele) for allele in alleles],
        "genic": genic,
        "nonsynonymous": nonsynonymous,
        "missing": missing,
    }
    return payload

def write_sample_payload(
    out_path: Path,
    sample_id: str,
    allele_sets: Dict[str, Set],
    missing: List[List],
    overwrite: bool,
) -> None:
    payload = sample_payload(sample_id, allele_sets, missing)
    try:
        open_mode = "w" if overwrite else "x"
        with open(out_path, open_mode) as outfh:
            json.dump(payload, outfh, indent=2)
    except FileExistsError as exc:
        raise RuntimeError(f"Refusing to overwrite existing output file: {out_path}") from exc
    except Exception as exc:
        raise RuntimeError(
            f"Failed to write {out_path} "
            f"({len(payload['alleles'])} alleles, {len(payload['missing'])} missing blocks): "
            f"{type(exc).__name__}: {exc}"
        ) from exc
