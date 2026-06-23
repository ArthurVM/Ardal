#!/usr/bin/env python3
import argparse
import gc
import re
import sys, json, hashlib, csv
import bisect
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple
import numpy as np
from numpy.lib.format import open_memmap
from ardal.builder.schemas import (
    DATA_TYPES,
    DEFAULT_ALLELE_ID_FORMAT,
    DEFAULT_NONSYN_ALLELE_ID_FORMAT,
    V2_SAMPLE_SCHEMA_VERSION,
    binary_missing_ranges_payload,
    build_allele_model,
    build_matrix_metadata,
    column_masks_payload,
)

_ALLELE_ID_FORMAT: str | None = None
_ALLELE_ID_PATTERN: re.Pattern | None = None
_ALLELE_ID_POS_KEY: str | None = None
_ALLELE_ID_HAS_CHR: bool = True
_ALLELE_ID_HAS_GENE: bool = False


def dump_json(
    path : Path | str,
    payload,
    indent : int | None = None,
):
    """Persist a JSON payload to disk with the configured indentation."""
    Path(path).write_text(json.dumps(payload, indent=indent))


def write_matrix_meta(
    *,
    matrix_path : str,
    format_name : str,
    dtype : str,
    n_rows : int,
    n_cols : int,
    headers : Dict[str, List[str]],
    missing_sites : Dict | None,
    json_indent : int | None,
    row_major : bool = True,
    data_sha256 : str | None = None,
    words_per_row : int | None = None,
    bits_per_word : int | None = None,
    row_stride_bytes : int | None = None,
    endianness : str | None = None,
    matrix_kind : str | None = None,
    allele_model : Dict[str, object] | None = None,
    sections : Dict[str, object] | None = None,
) -> str:
    """
    Description:
        Emit a `.meta` sidecar describing a saved matrix artifact.

    Inputs:
        matrix_path: Absolute or relative path to the serialized matrix file.
        format_name: Identifier describing the encoding (e.g. ardal.bin.v1).
        dtype: Data type string for entries in the matrix.
        n_rows / n_cols: Matrix dimensions for validation downstream.
        headers: Mapping that exposes GUIDs and allele identifiers.
        missing_sites: Optional structure tracking per-guid missing positions.
        json_indent: Indentation level for the JSON payload.
        row_major / data_sha256 / words_per_row / bits_per_word / row_stride_bytes / endianness:
            Additional serialization metadata used by consumers.

    Outputs:
        Returns the generated `.meta` path for convenience.

    Exceptions:
        Propagates any IOErrors raised while reading the matrix size or writing JSON.
    """
    matrix_file = Path(matrix_path)
    data_nbytes = matrix_file.stat().st_size if matrix_file.exists() else None

    meta = build_matrix_metadata(
        format_name=format_name,
        dtype=dtype,
        n_rows=n_rows,
        n_cols=n_cols,
        matrix_file=matrix_file.name,
        row_major=row_major,
        data_nbytes=data_nbytes,
        data_sha256=data_sha256,
        words_per_row=words_per_row,
        bits_per_word=bits_per_word,
        row_stride_bytes=row_stride_bytes,
        endianness=endianness,
        matrix_kind=matrix_kind,
        allele_model=allele_model,
        sections=sections,
    )

    payload = {"meta": meta, "headers": headers}

    if ( missing_sites is not None ) :
        payload.update(missing_sites)
    else:
        payload["column_masks"] = None

    meta_path = f"{matrix_path}.meta"
    dump_json(meta_path, payload, indent=json_indent)

    print(f"[Meta] {meta_path}")

    return meta_path


def close_memmap(
    mm : np.memmap,
):
    """Flush and close a numpy memmap safely."""
    if ( mm is None ):
        
        return
    mm.flush()
    mmap_obj = getattr(mm, "_mmap", None)
    if ( mmap_obj is not None ):
        mmap_obj.close()


def remove_file(
    path : Path,
):
    """Delete a file path while ignoring missing files."""
    try:
        Path(path).unlink()
    except FileNotFoundError:
        pass


def dense_work_path(
    prefix: str,
    keep_dense_output: bool,
) -> str:
    """Return the working dense path, using a hidden staging suffix unless .npy output is requested."""
    if keep_dense_output:
        return f"{prefix}.npy"
    return f"{prefix}.working.npy"


def derive_sample_id_from_path(
    json_path : Path,
) -> str:
    """Infer the sample identifier from a JSON filename."""
    stem = json_path.stem
    if ( stem.endswith("_snps") ):
        stem = stem[:-5]

    return stem


def normalize_missing_entries(
    entries : Iterable[object],
) -> List[str]:
    """Normalize missing site markers to strings.

    Supports string site keys as well as [chrom, pos] and (chrom, pos) tuples.
    """
    normalized = []
    for entry in entries or []:
        if isinstance(entry, (list, tuple)):
            if len(entry) == 2:
                chrom, pos = entry
                entry_str = position_key(str(chrom), pos)
            elif len(entry) == 1:
                entry_str = str(entry[0])
            else:
                entry_str = str(entry)
        else:
            entry_str = str(entry)

        if ( not _ALLELE_ID_HAS_CHR ):
            if ( "." in entry_str ):
                _, tail = entry_str.rsplit(".", 1)
                if ( tail.isdigit() ):
                    entry_str = tail
        normalized.append(entry_str)

    return normalized

def allele_id_from_fields(
    *,
    chrom,
    start,
    end,
    ref,
    alt,
    gene = None,
    allele_id_format : str | None = None,
) -> str:
    """Render an allele record using the requested matrix-time identifier format."""
    fmt = allele_id_format or DEFAULT_ALLELE_ID_FORMAT
    pos = start
    values = {
        "chr": chrom,
        "start": start,
        "end": end,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gene": gene,
    }
    try:
        return fmt.format(**values)
    except KeyError as exc:
        raise ValueError(f"Unsupported placeholder in allele_id_format: {exc}") from exc


def missing_blocks_from_entries(
    entries : Iterable[object],
) -> List[Tuple[str, int, int]]:
    """Parse v2 missing blocks as (chrom, start, end), preserving only valid numeric intervals."""
    blocks: List[Tuple[str, int, int]] = []
    for entry in entries or []:
        if not isinstance(entry, (list, tuple)) or len(entry) != 3:
            continue
        chrom, start, end = entry
        try:
            start_i = int(start)
            end_i = int(end)
        except (TypeError, ValueError):
            continue
        if end_i < start_i:
            start_i, end_i = end_i, start_i
        blocks.append((str(chrom), start_i, end_i))

    return blocks


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


def position_key(
    chrom : str,
    pos : str | int,
) -> str:
    """Compose a string key for positional lookups."""
    if ( chrom is None or chrom == "" ):
        return f"{pos}"
    
    return f"{chrom}.{pos}"

def missing_entry_contains_position(entry, chrom: str, pos: str | int) -> bool:
    """Return whether a legacy missing key or v2 missing block contains a site."""
    try:
        pos_i = int(pos)
    except (TypeError, ValueError):
        pos_i = None

    if isinstance(entry, (list, tuple)) and len(entry) == 3:
        block_chrom, start, end = entry
        if str(block_chrom) != str(chrom) or pos_i is None:
            return False
        try:
            return int(start) <= pos_i <= int(end)
        except (TypeError, ValueError):
            return False

    return str(entry) == position_key(str(chrom), pos)


def any_missing_entry_contains_position(entries, chrom: str, pos: str | int) -> bool:
    return any(missing_entry_contains_position(entry, chrom, pos) for entry in entries or [])


def sort_position_key(
    loc : str,
) -> Tuple[str, int, str]:
    """Produce a tuple that allows consistent ordering of positional keys."""
    if ( "." in loc ):
        chrom, pos = loc.rsplit(".", 1)
    else:
        chrom, pos = "", loc
    try:
        pos_val = int(pos)
        
        return chrom, 0, pos_val
    except ValueError:
        
        return chrom, 1, pos


def parse_allele_key(
    allele_id : str,
):
    """Split an allele identifier into chrom, position, reference, and alternate entries."""
    if ( _ALLELE_ID_PATTERN is None ):
        try:
            chrom, pos, ref, alt = allele_id.rsplit(".", 3)
        except ValueError as exc:
            raise ValueError(f"Unexpected allele identifier format: {allele_id}") from exc
        
        return chrom, pos, ref, alt

    match = _ALLELE_ID_PATTERN.match(allele_id)
    if ( not match ):
        raise ValueError(
            f"Allele ID '{allele_id}' does not match the format '{_ALLELE_ID_FORMAT}'."
        )

    parts = match.groupdict()
    chrom = parts.get("chr") if _ALLELE_ID_HAS_CHR else None
    pos = parts.get(_ALLELE_ID_POS_KEY)
    gene = parts.get("gene") if _ALLELE_ID_HAS_GENE else None
    ref = parts.get("ref")
    alt = parts.get("alt")

    if ( pos is None or ref is None or alt is None or ( _ALLELE_ID_HAS_CHR and chrom is None ) ):
        raise ValueError(
            "allele_id_format must include {pos}/{start}, {ref}, and {alt} placeholders."
        )

    if ( gene is not None ):
        pos = f"{gene}.{pos}"

    return chrom, pos, ref, alt


def compile_allele_id_format(
    allele_id_format : str,
) -> Tuple[re.Pattern, str, bool, bool]:
    """Compile an allele_id_format string into a regex pattern and return the position placeholder key."""
    if ( not allele_id_format ):
        raise ValueError("allele_id_format cannot be empty.")

    placeholders = re.findall(r"\{([^}]+)\}", allele_id_format)
    if ( not placeholders ):
        raise ValueError(
            "allele_id_format must include placeholders like {pos}, {ref}, {alt}."
        )

    allowed = {"chr", "gene", "pos", "start", "end", "ref", "alt"}
    invalid = sorted({p for p in placeholders if p not in allowed})
    if ( invalid ):
        raise ValueError(
            "Unsupported placeholders in allele_id_format: "
            + ", ".join(invalid)
            + f". Allowed: {', '.join(sorted(allowed))}."
        )

    seen = set()
    dupes = set()
    for p in placeholders:
        if p in seen:
            dupes.add(p)
        seen.add(p)
    if ( dupes ):
        raise ValueError(
            "allele_id_format contains duplicate placeholders: "
            + ", ".join(sorted(dupes))
        )

    if ( "pos" in placeholders and "start" in placeholders ):
        raise ValueError("allele_id_format must use either {pos} or {start}, not both.")

    pos_key = "pos" if "pos" in placeholders else "start"
    if ( pos_key not in placeholders ):
        raise ValueError("allele_id_format must include {pos} or {start}.")

    missing_required = []
    if ( "ref" not in placeholders ):
        missing_required.append("{ref}")
    if ( "alt" not in placeholders ):
        missing_required.append("{alt}")
    if ( missing_required ):
        raise ValueError(
            "allele_id_format missing required placeholders: "
            + ", ".join(missing_required)
        )

    pattern = re.escape(allele_id_format)
    replacements = {
        "ref": r"(?P<ref>.+)",
        "alt": r"(?P<alt>.+)",
        "chr": r"(?P<chr>.+)",
        "gene": r"(?P<gene>.+)",
        "start": r"(?P<start>\d+)",
        "pos": r"(?P<pos>\d+)",
        "end": r"(?P<end>\d+)",
    }
    for key, regex_pattern in replacements.items():
        escaped_placeholder = re.escape(f"{{{key}}}")
        if escaped_placeholder in pattern:
            pattern = pattern.replace(escaped_placeholder, regex_pattern)

    pattern = f"^{pattern}$"

    return re.compile(pattern), pos_key, ("chr" in placeholders), ("gene" in placeholders)


def configure_allele_id_format(
    allele_id_format : str | None,
):
    """Configure how allele identifiers are parsed within this module."""
    global _ALLELE_ID_FORMAT, _ALLELE_ID_PATTERN, _ALLELE_ID_POS_KEY, _ALLELE_ID_HAS_CHR, _ALLELE_ID_HAS_GENE

    if ( allele_id_format is None ):
        _ALLELE_ID_FORMAT = None
        _ALLELE_ID_PATTERN = None
        _ALLELE_ID_POS_KEY = None
        _ALLELE_ID_HAS_CHR = True
        _ALLELE_ID_HAS_GENE = False
        
        return

    pattern, pos_key, has_chr, has_gene = compile_allele_id_format(allele_id_format)
    _ALLELE_ID_FORMAT = allele_id_format
    _ALLELE_ID_PATTERN = pattern
    _ALLELE_ID_POS_KEY = pos_key
    _ALLELE_ID_HAS_CHR = has_chr
    _ALLELE_ID_HAS_GENE = has_gene


def allele_sort_key(
    allele_id : str,
) -> Tuple[str, int, object, str, str]:
    """Provide a deterministic ordering tuple for allele identifiers."""
    chrom, pos, ref, alt = parse_allele_key(allele_id)
    try:
        pos_val = int(pos)
        
        return chrom, 0, pos_val, ref, alt
    except ValueError:
        
        return chrom, 1, pos, ref, alt


def reorder_alleles_by_position(
    allele_to_idx : Dict[str, int],
) -> Tuple[Dict[str, int], List[str]]:
    """Provide a remapped allele index dict plus the ordered allele list."""
    ordered = sorted(allele_to_idx.keys(), key=allele_sort_key)
    remapped = {allele_id: idx for idx, allele_id in enumerate(ordered)}
    
    return remapped, ordered

def compact_index_ranges(indices: Iterable[int]) -> List[List[int]]:
    """Compact sorted integer indices into inclusive [start, end] ranges."""
    ranges: List[List[int]] = []
    sorted_indices = sorted(set(int(i) for i in indices))
    if not sorted_indices:
        return ranges
    start = prev = sorted_indices[0]
    for idx in sorted_indices[1:]:
        if idx == prev + 1:
            prev = idx
            continue
        ranges.append([start, prev])
        start = prev = idx
    ranges.append([start, prev])
    return ranges

def expand_index_ranges(ranges: Iterable[Sequence[int]]) -> set:
    expanded = set()
    for item in ranges or []:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            continue
        start, end = item
        expanded.update(range(int(start), int(end) + 1))
    return expanded


def allele_genomic_sites_from_json(
    json_files : Sequence[Path],
    data_type : str,
    allele_id_format : str | None,
    nonsyn_allele_id_format : str | None,
) -> Dict[str, Tuple[str, int]]:
    """Map matrix allele IDs to the genomic positions used for missing-block masks."""
    allele_sites: Dict[str, Tuple[str, int]] = {}
    total = len(json_files)
    for i, jf in enumerate(json_files, 1):
        if ( i % 25 == 0 or i == total ):
            print(f"\r[Allele genomic sites] {i}/{total}: {jf}\x1b[K", end="", flush=True)
        for _sample_id, allele_entries, _missing in iter_sample_projected_entries(
            jf,
            data_type=data_type,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
        ):
            for allele_id, chrom, pos in allele_entries:
                if allele_id not in allele_sites:
                    try:
                        allele_sites[allele_id] = (str(chrom), int(pos))
                    except (TypeError, ValueError):
                        pass
    print(flush=True)

    return allele_sites

def build_columns_by_genomic_site(
    ordered_alleles : Sequence[str],
    allele_to_idx : Dict[str, int],
    allele_genomic_sites : Dict[str, Tuple[str, int]],
) -> Dict[str, Tuple[List[int], List[List[int]]]]:
    """Index matrix columns by genomic coordinate for fast missing-block projection."""
    print("[Build missing-site column index]", flush=True)
    grouped: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(list))
    for allele_id in ordered_alleles:
        genomic_site = allele_genomic_sites.get(allele_id)
        if genomic_site is None:
            continue
        chrom, pos = genomic_site
        grouped[str(chrom)][int(pos)].append(allele_to_idx[allele_id])

    indexed: Dict[str, Tuple[List[int], List[List[int]]]] = {}
    for chrom, pos_to_cols in grouped.items():
        positions = sorted(pos_to_cols)
        columns = [pos_to_cols[pos] for pos in positions]
        indexed[chrom] = (positions, columns)
    print(f"[Build missing-site column index] {len(indexed)} chromosomes", flush=True)

    return indexed

def column_masks_from_missing_blocks(
    sample_to_idx,
    sample_missing_positions : Dict[str, list],
    columns_by_site : Dict[str, Tuple[List[int], List[List[int]]]],
) -> Dict[str, List[List[int]]]:
    """Project missing intervals onto existing matrix allele columns."""
    column_masks = {}
    total = len(sample_to_idx)
    for i, sample_id in enumerate(sample_to_idx.keys(), 1):
        if ( i == 1 or i % 25 == 0 or i == total ):
            print(f"\r[Missing masks] {i}/{total}: {sample_id}\x1b[K", end="", flush=True)
        masked_cols = []
        for entry in sample_missing_positions.get(sample_id, []):
            if isinstance(entry, (list, tuple)) and len(entry) == 3:
                chrom, start, end = entry
                site_index = columns_by_site.get(str(chrom))
                if site_index is None:
                    continue
                positions, columns = site_index
                left = bisect.bisect_left(positions, int(start))
                right = bisect.bisect_right(positions, int(end))
                for cols_at_pos in columns[left:right]:
                    masked_cols.extend(cols_at_pos)
            else:
                if ( "." not in str(entry) ):
                    continue
                chrom, pos = str(entry).rsplit(".", 1)
                site_index = columns_by_site.get(chrom)
                if site_index is None:
                    continue
                positions, columns = site_index
                try:
                    pos_i = int(pos)
                except ValueError:
                    continue
                idx = bisect.bisect_left(positions, pos_i)
                if idx < len(positions) and positions[idx] == pos_i:
                    masked_cols.extend(columns[idx])
        compacted = compact_index_ranges(masked_cols)
        if ( compacted ):
            column_masks[sample_id] = compacted
    print(flush=True)

    return column_masks


def append_column_range(
    ranges : List[List[int]],
    start : int,
    end : int,
):
    """Append a [start, end) column range, merging with the previous range when adjacent."""
    if ( end <= start ):
        return
    if ( ranges and start <= ranges[-1][1] ):
        if ( end > ranges[-1][1] ):
            ranges[-1][1] = end
        return
    ranges.append([int(start), int(end)])


def missing_column_ranges_from_entries(
    entries : Iterable[object],
    columns_by_site : Dict[str, Tuple[List[int], List[List[int]]]],
) -> List[List[int]]:
    """Project one sample's missing genomic entries to compact matrix-column [start, end) ranges."""
    ranges: List[List[int]] = []
    for entry in entries or []:
        if isinstance(entry, (list, tuple)) and len(entry) == 3:
            chrom, start, end = entry
            site_index = columns_by_site.get(str(chrom))
            if site_index is None:
                continue
            positions, columns = site_index
            left = bisect.bisect_left(positions, int(start))
            right = bisect.bisect_right(positions, int(end))
            for cols_at_pos in columns[left:right]:
                if ( not cols_at_pos ):
                    continue
                append_column_range(ranges, min(cols_at_pos), max(cols_at_pos) + 1)
            continue

        if ( "." not in str(entry) ):
            continue
        chrom, pos = str(entry).rsplit(".", 1)
        site_index = columns_by_site.get(chrom)
        if site_index is None:
            continue
        positions, columns = site_index
        try:
            pos_i = int(pos)
        except ValueError:
            continue
        idx = bisect.bisect_left(positions, pos_i)
        if idx < len(positions) and positions[idx] == pos_i:
            cols_at_pos = columns[idx]
            if ( cols_at_pos ):
                append_column_range(ranges, min(cols_at_pos), max(cols_at_pos) + 1)

    return ranges


def file_sha256(
    path : Path | str,
) -> str:
    """Compute sha256 for a binary artifact."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b""):
            h.update(chunk)

    return h.hexdigest()


def append_missing_range_sections(
    *,
    binpath : str,
    guids : Sequence[str],
    sample_missing_positions : Dict[str, list],
    columns_by_site : Dict[str, Tuple[List[int], List[List[int]]]],
    n_cols : int,
) -> Dict[str, object]:
    """Append per-sample missing column ranges to the allele matrix .bin file."""
    range_dtype = np.dtype("<u4") if n_cols <= np.iinfo(np.uint32).max else np.dtype("<u8")
    offset_dtype = np.dtype("<u8")
    n_samples = len(guids)
    offsets_offset = Path(binpath).stat().st_size
    offsets_nbytes = (n_samples + 1) * offset_dtype.itemsize

    print(
        f"[Append missing ranges] offsets={n_samples + 1:,} dtype={offset_dtype.str} "
        f"ranges dtype={range_dtype.str}",
        flush=True,
    )
    with open(binpath, "ab") as fh:
        fh.write(b"\0" * offsets_nbytes)
    offsets = np.memmap(
        binpath,
        mode="r+",
        dtype=offset_dtype,
        offset=offsets_offset,
        shape=(n_samples + 1,),
    )
    offsets[0] = 0

    ranges_offset = offsets_offset + offsets_nbytes
    total_ranges = 0
    with open(binpath, "ab") as ranges_fh:
        for i, sample_id in enumerate(guids, 1):
            if ( i == 1 or i % 25 == 0 or i == n_samples ):
                print(f"\r[Append missing ranges] {i}/{n_samples}: {sample_id}\x1b[K", end="", flush=True)
            sample_ranges = missing_column_ranges_from_entries(
                sample_missing_positions.get(sample_id, []),
                columns_by_site,
            )
            if ( sample_ranges ):
                arr = np.asarray(sample_ranges, dtype=range_dtype)
                ranges_fh.write(arr.tobytes(order="C"))
                total_ranges += int(arr.shape[0])
            offsets[i] = total_ranges
    print(flush=True)
    offsets.flush()
    close_memmap(offsets)

    ranges_nbytes = total_ranges * 2 * range_dtype.itemsize
    return {
        "missing_offsets": {
            "offset": int(offsets_offset),
            "nbytes": int(offsets_nbytes),
            "dtype": offset_dtype.str,
            "length": int(n_samples + 1),
            "units": "range_rows",
        },
        "missing_ranges": {
            "offset": int(ranges_offset),
            "nbytes": int(ranges_nbytes),
            "dtype": range_dtype.str,
            "shape": [int(total_ranges), 2],
            "coordinate_system": "matrix_columns_0_based",
            "range_semantics": "start_inclusive_end_exclusive",
        },
    }


def site_sort_key(
    site_id : str,
) -> Tuple[str, int, object, str]:
    """Return a deterministic ordering tuple for site identifiers."""
    parts = site_id.rsplit(".", 2)
    if ( len(parts) == 3 ):
        chrom, pos, ref = parts
    elif ( len(parts) == 2 ):
        chrom, pos, ref = "", parts[0], parts[1]
    else:
        raise ValueError(f"Unexpected site identifier format: {site_id}")
    try:
        pos_val = int(pos)
        
        return chrom, 0, pos_val, ref
    except ValueError:
        
        return chrom, 1, pos, ref


def site_identifier(
    chrom : str,
    pos : int | str,
    ref : str,
):
    """Construct the canonical identifier for a site."""
    if ( chrom is None or chrom == "" ):
        return f"{pos}.{ref}"
    return f"{chrom}.{pos}.{ref}"


def index_samples_and_alleles(
    json_files : Sequence[Path],
    data_type : str = "snps",
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
    collect_site_alt_samples : bool = False,
):
    """
    Description:
        Traverse all sample JSON artifacts to build lookup tables for samples, alleles, and site-level stats.

    Inputs:
        json_files: Sequence of JSON paths emitted by upstream allele callers.

    Outputs:
        Returns mapping dictionaries plus site/sample metadata required for downstream matrix construction.

    Exceptions:
        Propagates parsing errors raised while decoding malformed JSON.
    """
    sample_to_idx, allele_to_idx = {}, {}
    site_alt_counts = defaultdict(lambda: defaultdict(int))
    site_alt_samples = defaultdict(set) if collect_site_alt_samples else None
    site_alt_sample_counts = defaultdict(int)
    site_reference = {}
    site_coordinates = {}
    sample_missing_positions: Dict[str, list] = defaultdict(list)
    
    for ( i, jf ) in enumerate(json_files, 1):
        print(f"\r[Indexing] {i}/{len(json_files)}: {jf}\x1b[K", end="", flush=True)
        for sample_id, allele_entries, missing in iter_sample_projected_entries(
            jf,
            data_type=data_type,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
        ):
            if ( sample_id not in sample_to_idx ):
                sample_to_idx[sample_id] = len(sample_to_idx)
            sample_idx = sample_to_idx[sample_id]
            if ( missing ):
                sample_missing_positions[sample_id].extend(missing)
            if ( not allele_entries ):
                continue
            seen = set()
            per_site_in_sample = set()
            for allele_id, _genomic_chrom, _genomic_pos in allele_entries:
                if ( allele_id in seen ):
                    continue
                seen.add(allele_id)
                chrom, pos, ref, _ = parse_allele_key(allele_id)
                site_id = site_identifier(chrom, pos, ref)
                site_reference[site_id] = ref
                site_coordinates[site_id] = {"chrom": chrom, "pos": pos}
                if ( allele_id not in allele_to_idx ):
                    allele_to_idx[allele_id] = len(allele_to_idx)
                site_alt_counts[site_id][allele_id] += 1
                per_site_in_sample.add(site_id)
            for site_id in per_site_in_sample:
                if ( site_alt_samples is not None ):
                    site_alt_samples[site_id].add(sample_idx)
                site_alt_sample_counts[site_id] += 1
    print()
    
    return (
        sample_to_idx,
        allele_to_idx,
        site_reference,
        site_coordinates,
        site_alt_counts,
        site_alt_samples if site_alt_samples is not None else site_alt_sample_counts,
        sample_missing_positions,
    )


def derive_site_statistics(
    sample_to_idx,
    allele_to_idx,
    site_reference,
    site_coordinates,
    site_alt_counts,
    site_alt_samples,
):
    """
    Description:
        Summarize the major/minor allele makeup for every site observed in the dataset.

    Inputs:
        sample_to_idx / allele_to_idx: Index lookups constructed earlier.
        site_reference / site_coordinates: Context for reference alleles and positions.
        site_alt_counts / site_alt_samples: Aggregated counts for alt alleles and carriers.

    Outputs:
        Tuple of (site_metadata dict, minor column labels, alt-to-idx map, ref-to-idx map).

    Exceptions:
        Propagates parsing errors arising from malformed allele identifiers.
    """
    n_samples = len(sample_to_idx)
    site_metadata = OrderedDict()
    minor_columns = []
    minor_alt_to_idx = {}
    minor_ref_to_idx = {}

    ## iterate sites in sorted order for deterministic output
    for site_id in sorted(site_reference.keys(), key=site_sort_key):
        ref_base = site_reference[site_id]
        coords = site_coordinates[site_id]
        alt_counts_map = site_alt_counts.get(site_id, {})
        ## convert to regular dict for stable iteration later
        alt_counts = dict(alt_counts_map)
        alt_sample_entry = site_alt_samples.get(site_id, set())
        alt_sample_count = (
            len(alt_sample_entry)
            if isinstance(alt_sample_entry, set)
            else int(alt_sample_entry)
        )
        ref_count = n_samples - alt_sample_count

        major_label = "REF"
        major_count = ref_count

        for allele_id, count in alt_counts.items():
            if ( count > major_count ):
                major_label = allele_id
                major_count = count
            elif ( count == major_count ):
                if ( major_label == "REF" ):
                    continue
                if ( major_label != "REF" and allele_id < major_label ):
                    major_label = allele_id

        if ( major_label == "REF" ):
            major_type = "REF"
            major_base = ref_base
            major_allele_id = None
        else:
            major_type = "ALT"
            _, _, _, alt_base = parse_allele_key(major_label)
            major_base = alt_base
            major_allele_id = major_label

        minor_alt_labels = []
        alt_items = sorted(alt_counts.items(), key=lambda kv: allele_to_idx.get(kv[0], sys.maxsize))
        for allele_id, _count in alt_items:
            if ( allele_id == major_label ):
                continue
            idx = len(minor_columns)
            minor_columns.append(allele_id)
            minor_alt_to_idx[allele_id] = idx
            minor_alt_labels.append(allele_id)

        ref_minor_label = None
        if ( major_label != "REF" ):
            ref_minor_label = f"{site_id}.REF"
            idx = len(minor_columns)
            minor_columns.append(ref_minor_label)
            minor_ref_to_idx[site_id] = idx

        alt_allele_details = []
        for allele_id, count in alt_items:
            _, _, _, alt_base = parse_allele_key(allele_id)
            alt_allele_details.append(
                {
                    "allele_id": allele_id,
                    "allele": alt_base,
                    "count": count,
                    "is_minor": allele_id in minor_alt_to_idx,
                }
            )

        site_metadata[site_id] = {
            "chrom": coords["chrom"],
            "pos": coords["pos"],
            "ref": ref_base,
            "total_samples": n_samples,
            "alt_sample_count": alt_sample_count,
            "major": {
                "type": major_type,
                "allele": major_base,
                "allele_id": major_allele_id,
                "count": major_count,
            },
            "alt_alleles": alt_allele_details,
            "minor_columns": {
                "ref": ref_minor_label,
                "alt": minor_alt_labels,
            },
        }
    
    return site_metadata, minor_columns, minor_alt_to_idx, minor_ref_to_idx


def build_dense_matrices(
    json_files : Sequence[Path],
    sample_to_idx,
    allele_to_idx,
    ref_npy_path,
    minor_npy_path,
    minor_columns,
    minor_alt_to_idx,
    sample_masked_columns : Sequence[set],
    emit_minor_matrix : bool,
    data_type : str = "snps",
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
):
    """
    Description:
        Populate dense reference and minor-allele matrices directly from sample JSON payloads.

    Inputs:
        json_files: Sequence of per-sample allele JSON files.
        sample_to_idx / allele_to_idx: Lookup maps derived earlier.
        ref_npy_path / minor_npy_path: Memmap destinations for each dense matrix.
        minor_columns / minor_alt_to_idx: Metadata describing the agnostic/minor layout.
        sample_masked_columns: Indexed sets describing matrix columns masked per sample.
        emit_minor_matrix: Flag toggling whether the agnostic matrix should be produced.

    Outputs:
        Returns the reference and minor memmap handles for downstream serialization.

    Exceptions:
        Propagates parsing errors from malformed allele identifiers.
    """
    n_rows, n_ref_cols = len(sample_to_idx), len(allele_to_idx)
    n_minor_cols = len(minor_columns)
    print(f"[Plan] Dense (ref): rows={n_rows:,} cols={n_ref_cols:,} bytes≈{n_rows*n_ref_cols/2**30:.2f} GiB (uint8)")
    if ( n_minor_cols ):
        print(f"[Plan] Dense (minor): rows={n_rows:,} cols={n_minor_cols:,} bytes≈{n_rows*n_minor_cols/2**30:.2f} GiB (uint8)")
    else:
        print("[Plan] Dense (minor): no columns (all sites reference-major)")

    dense_ref = open_memmap(ref_npy_path, mode="w+", dtype=np.uint8, shape=(n_rows, n_ref_cols))
    dense_ref[:] = 0

    dense_minor = None
    if ( emit_minor_matrix and n_minor_cols ):
        dense_minor = open_memmap(minor_npy_path, mode="w+", dtype=np.uint8, shape=(n_rows, n_minor_cols))
        dense_minor[:] = 0

    for ( i, jf ) in enumerate(json_files, 1):
        print(f"\r[Populate dense] {i}/{len(json_files)}: {jf}\x1b[K", end="", flush=True)
        for sample_id, allele_entries, _missing in iter_sample_projected_entries(
            jf,
            data_type=data_type,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
        ):
            r = sample_to_idx.get(sample_id)
            if ( r is None or not allele_entries ):
                continue
            seen = set()
            filtered_alleles = []
            masked_columns = sample_masked_columns[r] if sample_masked_columns else set()
            for allele_id, genomic_chrom, genomic_pos in allele_entries:
                if ( allele_id in seen ):
                    continue
                allele_idx = allele_to_idx.get(allele_id)
                if ( allele_idx is None ):
                    continue
                chrom, pos, ref, alt = parse_allele_key(allele_id)
                if ( masked_columns and allele_idx in masked_columns ):
                    continue
                seen.add(allele_id)
                filtered_alleles.append(allele_id)
            if ( not filtered_alleles ):
                continue

            ref_cols = [allele_to_idx[a] for a in filtered_alleles]
            dense_ref[r, ref_cols] = 1

            if ( dense_minor is None ):
                continue
            for allele_id in filtered_alleles:
                minor_col = minor_alt_to_idx.get(allele_id)
                if ( minor_col is not None ):
                    dense_minor[r, minor_col] = 1

    print("\n[Done] Dense matrices built.")
    dense_ref.flush()
    if ( dense_minor is not None ):
        dense_minor.flush()

    return dense_ref, dense_minor


def apply_reference_minor_columns(
    dense_minor,
    minor_ref_to_idx,
    site_alt_samples,
):
    """Populate columns representing reference alleles within the minor matrix."""
    if ( dense_minor is None or not minor_ref_to_idx ):
        
        return
    for site_id, col_idx in minor_ref_to_idx.items():
        dense_minor[:, col_idx] = 1
        alt_samples = site_alt_samples.get(site_id, set())
        if ( alt_samples ):
            dense_minor[list(alt_samples), col_idx] = 0
    dense_minor.flush()


def write_dense_artifacts(
    dense_mm,
    guids,
    alleles,
    output_prefix,
    emit_npy: bool,
    emit_npz: bool,
    emit_csv: bool,
    missing_sites=None,
    allele_model=None,
    json_indent: int | None = None,
):
    """
    Description:
        Emit optional dense `.npy` / `.npz` matrices (and CSVs) plus aligned metadata files.

    Inputs:
        dense_mm: Memmap view of the dense matrix being serialized.
        guids / alleles: Header vectors applied to every output.
        output_prefix: Prefix used for file naming.
        emit_npy / emit_npz / emit_csv: Flags controlling which artifacts are written.
        missing_sites: Optional structure capturing per-guid missing site indices.
        json_indent: Indentation level for metadata outputs.

    Outputs:
        None; side effects include files on disk.

    Exceptions:
        Propagates IOErrors raised while writing disk artifacts.
    """
    if ( dense_mm is None ):
        
        return

    dense_mm.flush()
    headers = {"guids": guids, "alleles": alleles}
    dtype_str = str(dense_mm.dtype)
    n_rows = len(guids)
    n_cols = len(alleles)

    if ( emit_npz ):
        npz_path = f"{output_prefix}.npz"
        print("[Write] .npz (dense)…")
        np.savez_compressed(
            npz_path,
            matrix=dense_mm,
            guids=np.array(guids, dtype=object),
            alleles=np.array(alleles, dtype=object),
        )
        write_matrix_meta(
            matrix_path=npz_path,
            format_name="ardal.npz.v1",
            dtype=dtype_str,
            n_rows=n_rows,
            n_cols=n_cols,
            headers=headers,
            missing_sites=missing_sites,
            allele_model=allele_model,
            json_indent=json_indent,
        )

    if ( emit_csv ):
        csv_path = f"{output_prefix}.csv"
        print("[Write] .csv (dense wide)…")
        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_id", *alleles])
            for r, sid in enumerate(guids):
                writer.writerow([sid, *dense_mm[r, :]])

    if ( emit_npy ):
        npy_path = f"{output_prefix}.npy"
        write_matrix_meta(
            matrix_path=npy_path,
            format_name="ardal.npy.v1",
            dtype=dtype_str,
            n_rows=n_rows,
            n_cols=n_cols,
            headers=headers,
            missing_sites=missing_sites,
            allele_model=allele_model,
            json_indent=json_indent,
        )


def bitpack_from_dense(
    dense_mm,
    n_cols,
    output_prefix,
    compute_sha256 = False,
):
    """
    Description:
        Convert a dense boolean matrix into a bitpacked `.bin` file for efficient storage.

    Inputs:
        dense_mm: Dense matrix memmap view.
        n_cols: Number of allele columns encoded inside the matrix.
        output_prefix: Prefix used when writing the `.bin` artifact.
        compute_sha256: Toggles digest computation for integrity verification.

    Outputs:
        Returns `(binpath, sha256, words_per_row)` for downstream metadata.

    Exceptions:
        Propagates IOErrors encountered while writing the bitpacked file or hashing it.
    """
    words = (n_cols + 63) // 64
    dtype = np.dtype("<u8")  ## little-endian uint64
    binpath = f"{output_prefix}.bin"

    print(f"[Plan] Bitpack: words/row={words}, file≈{dense_mm.shape[0]*words*8/2**30:.2f} GiB")
    if ( words == 0 ):
        ## No columns → nothing to bitpack; emit empty file (or reuse existing)
        Path(binpath).write_bytes(b"")
        sha256 = hashlib.sha256(b"").hexdigest() if compute_sha256 else None
        
        return binpath, sha256, words

    mm_bin = np.memmap(binpath, mode="w+", dtype=dtype, shape=(dense_mm.shape[0], words), order="C")
    mm_bin[:] = 0

    ## Precompute powers of two up to 64 bits
    pow2 = (np.uint64(1) << np.arange(64, dtype=np.uint64))

    for r in range(dense_mm.shape[0]):
        row = dense_mm[r, :]
        ## pack 64-bit chunks
        for w in range(words):
            start = w * 64
            end   = min(start + 64, n_cols)
            chunk = row[start:end]
            if ( chunk.size == 0 ):
                continue
            ## build mask: sum of 2^bit for bits where chunk == 1
            mask = (chunk.astype(np.uint64) * pow2[:chunk.size]).sum(dtype=np.uint64)
            mm_bin[r, w] = mask

        if ( (r + 1) % 1000 == 0 or r == dense_mm.shape[0] - 1 ):
            print(f"\r[Bitpack] row {r+1}/{dense_mm.shape[0]}\x1b[K", end="", flush=True)
    print("\n[Done] Bitpack written.")
    mm_bin.flush()

    ## Optional: sha256 of .bin
    sha256 = None
    if ( compute_sha256 ):
        print("[Hash] sha256 of .bin …")
        h = hashlib.sha256()
        with open(binpath, "rb") as f:
            for chunk in iter(lambda: f.read(8*1024*1024), b""):
                h.update(chunk)
        sha256 = h.hexdigest()

    ## Return paths and hash
    
    return binpath, sha256, words


def write_bitpacked_from_json(
    json_files : Sequence[Path],
    sample_to_idx,
    allele_to_idx,
    output_prefix,
    guids,
    alleles,
    compute_sha256 = False,
    data_type : str = "snps",
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
):
    """
    Write the bitpacked matrix directly from JSON payloads, avoiding dense staging.

    This is the scalable path for large cohorts when dense .npy/.npz/.csv
    artifacts are not requested.
    """
    n_rows = len(sample_to_idx)
    n_cols = len(allele_to_idx)
    words = (n_cols + 63) // 64
    dtype = np.dtype("<u8")
    binpath = f"{output_prefix}.bin"

    print(f"[Plan] Direct bitpack: rows={n_rows:,} cols={n_cols:,} words/row={words}, file≈{n_rows*words*8/2**30:.2f} GiB")
    if ( words == 0 ):
        Path(binpath).write_bytes(b"")
        allele_nbytes = 0
        return binpath, None, words, allele_nbytes

    mm_bin = np.memmap(binpath, mode="w+", dtype=dtype, shape=(n_rows, words), order="C")

    total = len(json_files)
    for ( i, jf ) in enumerate(json_files, 1):
        print(f"\r[Direct bitpack] {i}/{total}: {jf}\x1b[K", end="", flush=True)
        for sample_id, allele_entries, missing in iter_sample_projected_entries(
            jf,
            data_type=data_type,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
        ):
            r = sample_to_idx.get(sample_id)
            if ( r is None or not allele_entries ):
                continue
            row = mm_bin[r, :]
            seen = set()
            for allele_id, genomic_chrom, genomic_pos in allele_entries:
                if ( allele_id in seen ):
                    continue
                allele_idx = allele_to_idx.get(allele_id)
                if ( allele_idx is None ):
                    continue
                if ( missing and any_missing_entry_contains_position(missing, genomic_chrom, genomic_pos) ):
                    continue
                seen.add(allele_id)
                word_idx = allele_idx // 64
                bit_idx = allele_idx % 64
                row[word_idx] |= np.uint64(1) << np.uint64(bit_idx)

    print("\n[Done] Direct bitpack written.")
    mm_bin.flush()
    close_memmap(mm_bin)

    allele_nbytes = n_rows * words * dtype.itemsize

    return binpath, None, words, allele_nbytes


def emit_direct_bitpack_artifacts(
    json_files,
    sample_to_idx,
    allele_to_idx,
    sample_missing_positions,
    columns_by_site,
    prefix,
    guids,
    alleles,
    compute_sha256,
    json_indent : int | None,
    allele_model : Dict[str, object] | None = None,
    data_type : str = "snps",
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
):
    """Write direct bitpacked output and matching metadata."""
    binpath, _sha256, words, allele_nbytes = write_bitpacked_from_json(
        json_files,
        sample_to_idx,
        allele_to_idx,
        prefix,
        guids,
        alleles,
        compute_sha256=compute_sha256,
        data_type=data_type,
        allele_id_format=allele_id_format,
        nonsyn_allele_id_format=nonsyn_allele_id_format,
    )
    missing_sections = append_missing_range_sections(
        binpath=binpath,
        guids=guids,
        sample_missing_positions=sample_missing_positions,
        columns_by_site=columns_by_site,
        n_cols=len(alleles),
    )
    sections = {
        "allele_matrix": {
            "offset": 0,
            "nbytes": int(allele_nbytes),
            "dtype": "<u8",
            "shape": [int(len(sample_to_idx)), int(words)],
            "words_per_row": int(words),
            "bits_per_word": 64,
            "row_stride_bytes": int(words * 8),
        },
        **missing_sections,
    }
    sha256 = file_sha256(binpath) if compute_sha256 else None
    missing_sites = binary_missing_ranges_payload()
    write_matrix_meta(
        matrix_path=binpath,
        format_name="ardal.bin.v2",
        dtype="<u8",
        n_rows=len(sample_to_idx),
        n_cols=len(alleles),
        headers={"guids": guids, "alleles": alleles},
        missing_sites=missing_sites,
        json_indent=json_indent,
        data_sha256=sha256,
        words_per_row=words,
        bits_per_word=64,
        row_stride_bytes=words * 8,
        endianness="little",
        matrix_kind="allele_presence",
        allele_model=allele_model,
        sections=sections,
    )

    return None


def emit_bitpack_artifacts(
    dense_mm,
    n_cols,
    prefix,
    guids,
    alleles,
    n_rows,
    missing_sites,
    compute_sha256,
    json_indent : int | None,
    keep_dense_file : bool,
    npy_path : str,
    allele_model : Dict[str, object] | None = None,
):
    """
    Description:
        Bitpack a dense matrix, write the `.bin` artifact, emit its `.meta`, and clean up the memmap.

    Inputs:
        dense_mm: Memmap reference matrix ready for packing.
        n_cols / n_rows: Dimensions for metadata generation.
        prefix: Output prefix for the `.bin` file.
        guids / alleles: Headers stored alongside the metadata.
        missing_sites: Optional missing-site payload.
        compute_sha256: Whether to compute a digest for integrity verification.
        json_indent / keep_dense_file / npy_path: Serialization knobs controlling header indentation
            and whether the dense memmap is removed afterwards.

    Outputs:
        None; operates via side effects.

    Exceptions:
        Propagates IOErrors during serialization or hashing.
    """
    if ( dense_mm is None ):
        
        return None
    binpath, sha256, words = bitpack_from_dense(dense_mm, n_cols, prefix, compute_sha256)
    write_matrix_meta(
        matrix_path=binpath,
        format_name="ardal.bin.v1",
        dtype="<u8",
        n_rows=n_rows,
        n_cols=n_cols,
        headers={"guids": guids, "alleles": alleles},
        missing_sites=missing_sites,
        json_indent=json_indent,
        data_sha256=sha256,
        words_per_row=words,
        bits_per_word=64,
        row_stride_bytes=words * 8,
        endianness="little",
        matrix_kind="allele_presence",
        allele_model=allele_model,
    )
    close_memmap(dense_mm)
    gc.collect()
    if ( not keep_dense_file ):
        remove_file(Path(npy_path))

    return None


def write_site_metadata(
    output_prefix,
    n_samples,
    minor_columns,
    site_metadata,
    ref_prefix,
    minor_prefix,
    json_indent: int | None = None,
):
    """
    Description:
        Emit the metadata JSON describing which alleles are major/minor across sites.

    Inputs:
        output_prefix: Base output path for the metadata artifact.
        n_samples: Number of samples represented in the dataset.
        minor_columns / site_metadata: Structures captured earlier summarizing per-site behavior.
        ref_prefix / minor_prefix: Prefixes of the dense matrices referenced by the metadata.
        json_indent: Optional indentation for readability.

    Outputs:
        Returns the written metadata path.

    Exceptions:
        Propagates IOErrors while writing the JSON file.
    """
    meta_path = f"{output_prefix}.site_major_minor.json"
    payload = {
        "meta": {
            "format": "ardal.site_major_minor.v1",
            "n_samples": n_samples,
            "n_sites": len(site_metadata),
            "minor_matrix_columns": len(minor_columns),
            "reference_matrix_prefix": Path(ref_prefix).name,
            "minor_matrix_prefix": (
                Path(minor_prefix).name if ( minor_columns and minor_prefix ) else None
            ),
        },
        "minor_columns": minor_columns,
        "sites": site_metadata,
    }
    dump_json(meta_path, payload, indent=json_indent)
    print(f"[Site metadata] {meta_path}")

    return meta_path


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

def parse_args():
    """Parse command-line arguments for the conversion script."""
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

if ( __name__ == "__main__" ) :
    main()
