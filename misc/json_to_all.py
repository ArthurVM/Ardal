#!/usr/bin/env python3
import argparse
import gc
import re
import sys, json, hashlib, csv
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple
import numpy as np
from numpy.lib.format import open_memmap

DEFAULT_ALLELE_ID_FORMAT = "{chr}.{pos}.{ref}.{alt}"
_ALLELE_ID_FORMAT: str | None = None
_ALLELE_ID_PATTERN: re.Pattern | None = None
_ALLELE_ID_POS_KEY: str | None = None
_ALLELE_ID_HAS_CHR: bool = True
_ALLELE_ID_HAS_GENE: bool = False


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
) -> Dict[str, object]:
    """Create a compact metadata dict for any serialized matrix artifact."""
    def _as_int( value ):
        
        return None if value is None else int(value)

    return {
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
        "generated_by": "ardal::json_to_all.py",
    }


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


def iter_sample_entries(
    json_path : Path,
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
    with open(json_path) as fh:
        payload = json.load(fh)

    if ( isinstance(payload, dict) and ( "alleles" in payload or "missing" in payload ) ):
        sample_id = payload.get("sample_id") or derive_sample_id_from_path(json_path)
        alleles = list(payload.get("alleles") or [])
        missing = normalize_missing_entries(payload.get("missing") or [])
        yield sample_id, alleles, missing
        
        return

    if ( isinstance(payload, dict) ):
        for sample_id, alleles in payload.items():
            yield sample_id, list(alleles or []), []
        
        return

    raise ValueError(f"Unsupported JSON structure in {json_path}")


def position_key(
    chrom : str,
    pos : str | int,
) -> str:
    """Compose a string key for positional lookups."""
    if ( chrom is None or chrom == "" ):
        return f"{pos}"
    
    return f"{chrom}.{pos}"


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
    site_alt_samples = defaultdict(set)
    site_reference = {}
    site_coordinates = {}
    sample_missing_positions: Dict[str, set] = defaultdict(set)
    
    for ( i, jf ) in enumerate(json_files, 1):
        print(f"\r[Indexing] {i}/{len(json_files)}: {jf}\x1b[K", end="", flush=True)
        for sample_id, alleles, missing in iter_sample_entries(jf):
            if ( sample_id not in sample_to_idx ):
                sample_to_idx[sample_id] = len(sample_to_idx)
            sample_idx = sample_to_idx[sample_id]
            if ( missing ):
                sample_missing_positions[sample_id].update(missing)
            if ( not alleles ):
                continue
            seen = set()
            per_site_in_sample = set()
            for allele_id in alleles:
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
                site_alt_samples[site_id].add(sample_idx)
    print()
    
    return (
        sample_to_idx,
        allele_to_idx,
        site_reference,
        site_coordinates,
        site_alt_counts,
        site_alt_samples,
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
        alt_sample_set = site_alt_samples.get(site_id, set())
        ref_count = n_samples - len(alt_sample_set)

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
            "alt_sample_count": len(alt_sample_set),
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
    sample_missing_sets : Sequence[set],
    emit_minor_matrix : bool,
):
    """
    Description:
        Populate dense reference and minor-allele matrices directly from sample JSON payloads.

    Inputs:
        json_files: Sequence of per-sample allele JSON files.
        sample_to_idx / allele_to_idx: Lookup maps derived earlier.
        ref_npy_path / minor_npy_path: Memmap destinations for each dense matrix.
        minor_columns / minor_alt_to_idx: Metadata describing the agnostic/minor layout.
        sample_missing_sets: Indexed list describing per-sample missing positions.
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
        for sample_id, allele_ids, _missing in iter_sample_entries(jf):
            r = sample_to_idx.get(sample_id)
            if ( r is None or not allele_ids ):
                continue
            seen = set()
            filtered_alleles = []
            missing_positions = sample_missing_sets[r] if sample_missing_sets else set()
            for allele_id in allele_ids:
                if ( allele_id in seen ):
                    continue
                allele_idx = allele_to_idx.get(allele_id)
                if ( allele_idx is None ):
                    continue
                chrom, pos, ref, alt = parse_allele_key(allele_id)
                if ( missing_positions and position_key(chrom, pos) in missing_positions ):
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
    configure_allele_id_format(allele_id_format)
    (
        sample_to_idx,
        allele_to_idx,
        site_reference,
        site_coordinates,
        site_alt_counts,
        site_alt_samples,
        sample_missing_positions,
    ) = index_samples_and_alleles(json_files)

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

    variant_positions = {
        position_key(coords["chrom"], coords["pos"]) for coords in site_coordinates.values()
    }
    
    filtered_missing_by_sample: Dict[str, List[str]] = {}
    all_missing_positions = set()
    
    for sample_id in sample_to_idx.keys():
        raw_positions = sample_missing_positions.get(sample_id, set())
        filtered = sorted(
            (pos for pos in raw_positions if pos in variant_positions),
            key=sort_position_key,
        )
        filtered_missing_by_sample[sample_id] = filtered
        all_missing_positions.update(filtered)

    sorted_missing_positions = sorted(all_missing_positions, key=sort_position_key)
    # missing_site_key = {pos: idx for idx, pos in enumerate(sorted_missing_positions)}
    
    ## create a dictionary of sites and alleles which land in those sites
    missing_site_key = defaultdict(list)
    for a in ordered_alleles:
        chrom, pos, _ref, _alt = parse_allele_key(a)
        site = position_key(chrom, pos)
        missing_site_key[site].append(a)

    # print(missing_site_key)
    
    column_masks = defaultdict(list)
    for sample_id, positions in sample_missing_positions.items():
        for position in positions:
            for a in missing_site_key[position]:
                # print(sample_id, a)
                col = allele_to_idx[a]
                column_masks[sample_id].append(col)
    
    missing_sites_payload = {
        "column_masks" : column_masks
    }
    

    guids = [None] * n_rows
    
    for s, i in sample_to_idx.items():
        guids[i] = s
    alleles = ordered_alleles[:]

    sample_missing_sets_by_idx = [
        set(filtered_missing_by_sample.get(guid, [])) for guid in guids
    ]

    ref_prefix = f"{output_prefix}.ref"
    minor_prefix = f"{output_prefix}.agnostic"

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
        sample_missing_sets_by_idx,
        emit_agnostic,
    )

    if ( dense_minor is not None ):
        apply_reference_minor_columns(dense_minor, minor_ref_to_idx, site_alt_samples)

    emit_dense_artifacts = emit_npy or emit_npz or emit_csv
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
    ap.add_argument("--no-agnostic", action="store_true",
                    help="Skip generating reference-agnostic/minor allele matrices")
    ap.add_argument("--allele-id-format", type=str, default=None,
                    help=(
                        f"Allele ID format string (default: {DEFAULT_ALLELE_ID_FORMAT}; "
                        "{chr} optional for single-chrom references; {gene} supported for gene-relative IDs)"
                    ))
    ap.add_argument("--indent", type=int, dest="json_indent",
                    help="Indentation level for JSON artifacts (default uses compact JSON)")
    
    return ap.parse_args()

def main():
    """Entry point for invoking the JSON-to-matrix conversion workflow."""
    args = parse_args()
    json_dir = args.json_directory
    if ( not json_dir.is_dir() ):
        print(f"{json_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    json_files = sorted(p for p in json_dir.iterdir() if p.is_file() and p.suffix == ".json")
    if ( not json_files ):
        print(f"No JSON files found in {json_dir}")
        sys.exit(1)

    create_all_outputs(
        json_files,
        args.output_prefix,
        compute_sha256=args.compute_sha256,
        emit_npy=args.emit_npy,
        emit_npz=args.emit_npz,
        emit_csv=args.emit_csv,
        emit_agnostic=not args.no_agnostic,
        allele_id_format=args.allele_id_format,
        json_indent=args.json_indent,
    )

if ( __name__ == "__main__" ) :
    main()
