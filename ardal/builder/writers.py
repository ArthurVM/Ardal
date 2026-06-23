"""
Matrix artifact writers for dense, bitpacked, and metadata outputs.
"""

import csv
import gc
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from numpy.lib.format import open_memmap

from .alleles import parse_allele_key
from .missingness import any_missing_entry_contains_position, missing_column_ranges_from_entries
from .sample_json import iter_sample_projected_entries
from .schemas import build_matrix_metadata, binary_missing_ranges_payload

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
