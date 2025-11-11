#!/usr/bin/env python3
import os, sys, json, hashlib, csv
from collections import defaultdict, OrderedDict
from pathlib import Path
import numpy as np
from numpy.lib.format import open_memmap

def parse_allele_key(allele_id):
    try:
        chrom, pos, ref, alt = allele_id.rsplit(".", 3)
    except ValueError as exc:
        raise ValueError(f"Unexpected allele identifier format: {allele_id}") from exc
    return chrom, pos, ref, alt

def site_identifier(chrom, pos, ref):
    return f"{chrom}.{pos}.{ref}"

def index_samples_and_alleles(json_files):
    sample_to_idx, allele_to_idx = {}, {}
    site_alt_counts = defaultdict(lambda: defaultdict(int))
    site_alt_samples = defaultdict(set)
    site_reference = {}
    site_coordinates = {}
    for i, jf in enumerate(json_files, 1):
        print(f"\r[Indexing] {i}/{len(json_files)}: {jf}", end="", flush=True)
        with open(jf) as f:
            data = json.load(f)
        for s, alleles in data.items():
            if s not in sample_to_idx:
                sample_to_idx[s] = len(sample_to_idx)
            sample_idx = sample_to_idx[s]
            if not alleles:
                continue
            seen = set()
            per_site_in_sample = set()
            for a in alleles:
                if a in seen:
                    continue
                seen.add(a)
                chrom, pos, ref, alt = parse_allele_key(a)
                site_id = site_identifier(chrom, pos, ref)
                site_reference[site_id] = ref
                site_coordinates[site_id] = {"chrom": chrom, "pos": pos}
                if a not in allele_to_idx:
                    allele_to_idx[a] = len(allele_to_idx)
                site_alt_counts[site_id][a] += 1
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
    )

def derive_site_statistics(
    sample_to_idx,
    allele_to_idx,
    site_reference,
    site_coordinates,
    site_alt_counts,
    site_alt_samples,
):
    n_samples = len(sample_to_idx)
    site_metadata = OrderedDict()
    minor_columns = []
    minor_alt_to_idx = {}
    minor_ref_to_idx = {}

    # Iterate sites in sorted order for deterministic output
    for site_id in sorted(site_reference.keys()):
        ref_base = site_reference[site_id]
        coords = site_coordinates[site_id]
        alt_counts_map = site_alt_counts.get(site_id, {})
        # Convert to regular dict for stable iteration later
        alt_counts = dict(alt_counts_map)
        alt_sample_set = site_alt_samples.get(site_id, set())
        ref_count = n_samples - len(alt_sample_set)

        major_label = "REF"
        major_count = ref_count

        for allele_id, count in alt_counts.items():
            if count > major_count:
                major_label = allele_id
                major_count = count
            elif count == major_count:
                if major_label == "REF":
                    continue
                if major_label != "REF" and allele_id < major_label:
                    major_label = allele_id

        if major_label == "REF":
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
            if allele_id == major_label:
                continue
            idx = len(minor_columns)
            minor_columns.append(allele_id)
            minor_alt_to_idx[allele_id] = idx
            minor_alt_labels.append(allele_id)

        ref_minor_label = None
        if major_label != "REF":
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
    json_files,
    sample_to_idx,
    allele_to_idx,
    ref_npy_path,
    minor_npy_path,
    minor_columns,
    minor_alt_to_idx,
):
    n_rows, n_ref_cols = len(sample_to_idx), len(allele_to_idx)
    n_minor_cols = len(minor_columns)
    print(f"[Plan] Dense (ref): rows={n_rows:,} cols={n_ref_cols:,} bytes≈{n_rows*n_ref_cols/2**30:.2f} GiB (uint8)")
    if n_minor_cols:
        print(f"[Plan] Dense (minor): rows={n_rows:,} cols={n_minor_cols:,} bytes≈{n_rows*n_minor_cols/2**30:.2f} GiB (uint8)")
    else:
        print("[Plan] Dense (minor): no columns (all sites reference-major)")

    dense_ref = open_memmap(ref_npy_path, mode="w+", dtype=np.uint8, shape=(n_rows, n_ref_cols))
    dense_ref[:] = 0

    dense_minor = None
    if n_minor_cols:
        dense_minor = open_memmap(minor_npy_path, mode="w+", dtype=np.uint8, shape=(n_rows, n_minor_cols))
        dense_minor[:] = 0

    for i, jf in enumerate(json_files, 1):
        print(f"\r[Populate dense] {i}/{len(json_files)}: {jf}", end="", flush=True)
        with open(jf) as f:
            data = json.load(f)
        for s, allele_ids in data.items():
            r = sample_to_idx.get(s)
            if r is None or not allele_ids:
                continue
            seen = set()
            unique_alleles = []
            for allele_id in allele_ids:
                if allele_id in seen:
                    continue
                if allele_id not in allele_to_idx:
                    continue
                seen.add(allele_id)
                unique_alleles.append(allele_id)
            if not unique_alleles:
                continue

            ref_cols = [allele_to_idx[a] for a in unique_alleles]
            dense_ref[r, ref_cols] = 1

            if dense_minor is None:
                continue
            for allele_id in unique_alleles:
                minor_col = minor_alt_to_idx.get(allele_id)
                if minor_col is not None:
                    dense_minor[r, minor_col] = 1

    print("\n[Done] Dense matrices built.")
    dense_ref.flush()
    if dense_minor is not None:
        dense_minor.flush()
    return dense_ref, dense_minor

def apply_reference_minor_columns(dense_minor, minor_ref_to_idx, site_alt_samples):
    if dense_minor is None or not minor_ref_to_idx:
        return
    for site_id, col_idx in minor_ref_to_idx.items():
        dense_minor[:, col_idx] = 1
        alt_samples = site_alt_samples.get(site_id, set())
        if alt_samples:
            dense_minor[list(alt_samples), col_idx] = 0
    dense_minor.flush()

def write_dense_artifacts(dense_mm, guids, alleles, output_prefix):
    npy_path = f"{output_prefix}.npy"
    npz_path = f"{output_prefix}.npz"
    csv_path = f"{output_prefix}.csv"
    dense_head = f"{output_prefix}.dense.headers.json"

    # NPZ (compressed): matrix + headers
    print("[Write] .npz (dense)…")
    # Passing the memmap is fine; savez_compressed will stream from it.
    np.savez_compressed(
        npz_path,
        matrix=dense_mm,
        guids=np.array(guids, dtype=object),
        alleles=np.array(alleles, dtype=object),
    )

    # CSV (wide): header + 0/1 rows
    # print("[Write] .csv (dense wide)…")
    # with open(csv_path, "w", newline="") as fh:
    #     w = csv.writer(fh)
    #     w.writerow(["sample_id", *alleles])
    #     for r, sid in enumerate(guids):
    #         # csv.writer can take numpy scalars; this will stream row-wise
    #         w.writerow([sid, *dense_mm[r, :]])

    # Dense headers JSON
    meta = {
        "meta": {
            "format": "ardal.dense.v1",
            "dtype": "uint8",
            "row_major": True,
            "n_rows": len(guids),
            "n_cols": len(alleles),
            "matrix_file": Path(npy_path).name,
            "npz_file": Path(npz_path).name,
            "csv_file": Path(csv_path).name,
        },
        "headers": {
            "guids": guids,
            "alleles": alleles,
        },
    }
    Path(dense_head).write_text(json.dumps(meta))
    print(f"[Headers] {dense_head}")
    return npy_path, npz_path, csv_path, dense_head

def bitpack_from_dense(dense_mm, n_cols, output_prefix, compute_sha256=False):
    words = (n_cols + 63) // 64
    dtype = np.dtype("<u8")  # little-endian uint64
    binpath = f"{output_prefix}.bin"
    headpath = f"{output_prefix}.bitpack.headers.json"

    print(f"[Plan] Bitpack: words/row={words}, file≈{dense_mm.shape[0]*words*8/2**30:.2f} GiB")
    if words == 0:
        # No columns → nothing to bitpack; emit empty file (or reuse existing)
        Path(binpath).write_bytes(b"")
        sha256 = hashlib.sha256(b"").hexdigest() if compute_sha256 else None
        return binpath, headpath, sha256, words

    mm_bin = np.memmap(binpath, mode="w+", dtype=dtype, shape=(dense_mm.shape[0], words), order="C")
    mm_bin[:] = 0

    # Precompute powers of two up to 64 bits
    pow2 = (np.uint64(1) << np.arange(64, dtype=np.uint64))

    for r in range(dense_mm.shape[0]):
        row = dense_mm[r, :]
        # pack 64-bit chunks
        for w in range(words):
            start = w * 64
            end   = min(start + 64, n_cols)
            chunk = row[start:end]
            if chunk.size == 0:
                continue
            # build mask: sum of 2^bit for bits where chunk == 1
            mask = (chunk.astype(np.uint64) * pow2[:chunk.size]).sum(dtype=np.uint64)
            mm_bin[r, w] = mask

        if (r + 1) % 1000 == 0 or r == dense_mm.shape[0] - 1:
            print(f"\r[Bitpack] row {r+1}/{dense_mm.shape[0]}", end="", flush=True)
    print("\n[Done] Bitpack written.")
    mm_bin.flush()

    # Optional: sha256 of .bin
    sha256 = None
    if compute_sha256:
        print("[Hash] sha256 of .bin …")
        h = hashlib.sha256()
        with open(binpath, "rb") as f:
            for chunk in iter(lambda: f.read(8*1024*1024), b""):
                h.update(chunk)
        sha256 = h.hexdigest()

    # Return paths and hash
    return binpath, headpath, sha256, words

def write_bitpack_headers(headpath, binpath, sha256, n_rows, n_cols, words, guids, alleles, output_prefix):
    est_bytes = n_rows * words * 8
    meta = {
        "meta": {
            "format": "ardal.bitpack.v1",
            "dtype": "<u8",
            "endianness": "little",
            "row_major": True,
            "n_rows": n_rows,
            "n_cols": n_cols,
            "words_per_row": words,
            "bits_per_word": 64,
            "row_stride_bytes": words * 8,
            "data_file": Path(binpath).name,
            "data_nbytes": est_bytes,
            "data_sha256": sha256,
        },
        "headers": {
            "guids": guids,
            "alleles": alleles,
        },
    }
    Path(headpath).write_text(json.dumps(meta))
    Path(f"{output_prefix}.headers.json").write_text(json.dumps(meta["headers"]))
    print(f"[Headers] {headpath}")

def write_site_metadata(output_prefix, n_samples, minor_columns, site_metadata, ref_prefix, minor_prefix):
    meta_path = f"{output_prefix}.site_major_minor.json"
    payload = {
        "meta": {
            "format": "ardal.site_major_minor.v1",
            "n_samples": n_samples,
            "n_sites": len(site_metadata),
            "minor_matrix_columns": len(minor_columns),
            "reference_matrix_prefix": Path(ref_prefix).name,
            "minor_matrix_prefix": Path(minor_prefix).name if minor_columns else None,
        },
        "minor_columns": minor_columns,
        "sites": site_metadata,
    }
    Path(meta_path).write_text(json.dumps(payload))
    print(f"[Site metadata] {meta_path}")
    return meta_path

def create_all_outputs(json_files, output_prefix, compute_sha256=False):
    (
        sample_to_idx,
        allele_to_idx,
        site_reference,
        site_coordinates,
        site_alt_counts,
        site_alt_samples,
    ) = index_samples_and_alleles(json_files)

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
    alleles = [None] * n_ref_cols
    for a, j in allele_to_idx.items():
        alleles[j] = a

    ref_prefix = f"{output_prefix}.ref"
    minor_prefix = f"{output_prefix}.agnostic"

    ref_npy_path = f"{ref_prefix}.npy"
    minor_npy_path = f"{minor_prefix}.npy"

    dense_ref, dense_minor = build_dense_matrices(
        json_files,
        sample_to_idx,
        allele_to_idx,
        ref_npy_path,
        minor_npy_path,
        minor_columns,
        minor_alt_to_idx,
    )

    apply_reference_minor_columns(dense_minor, minor_ref_to_idx, site_alt_samples)

    # ---- dense artifacts (reference) ----
    write_dense_artifacts(dense_ref, guids, alleles, ref_prefix)

    # ---- dense artifacts (reference-agnostic/minor alleles) ----
    if dense_minor is not None:
        write_dense_artifacts(dense_minor, guids, minor_columns, minor_prefix)

    # ---- bitpack from dense (reference) ----
    ref_binpath, ref_bithead_path, ref_sha256, ref_words = bitpack_from_dense(
        dense_ref, n_ref_cols, ref_prefix, compute_sha256
    )
    write_bitpack_headers(
        ref_bithead_path, ref_binpath, ref_sha256, n_rows, n_ref_cols, ref_words, guids, alleles, ref_prefix
    )

    # ---- bitpack from dense (reference-agnostic) ----
    if dense_minor is not None:
        minor_binpath, minor_bithead_path, minor_sha256, minor_words = bitpack_from_dense(
            dense_minor, n_minor_cols, minor_prefix, compute_sha256
        )
        write_bitpack_headers(
            minor_bithead_path,
            minor_binpath,
            minor_sha256,
            n_rows,
            n_minor_cols,
            minor_words,
            guids,
            minor_columns,
            minor_prefix,
        )

    # ---- site metadata summarising major/minor alleles ----
    write_site_metadata(
        output_prefix,
        n_rows,
        minor_columns,
        site_metadata,
        ref_prefix,
        minor_prefix,
    )

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <json_directory> <output_prefix> [--sha256]")
        print("Outputs: <output_prefix>.ref.* (reference-based) and <output_prefix>.agnostic.* (minor allele) artifacts.")
        sys.exit(1)
    json_dir = sys.argv[1]
    output_prefix = sys.argv[2]
    compute_sha256 = ("--sha256" in sys.argv[3:])

    json_files = sorted(
        os.path.join(json_dir, f) for f in os.listdir(json_dir) if f.endswith(".json")
    )
    if not json_files:
        print(f"No JSON files found in {json_dir}")
        sys.exit(1)

    create_all_outputs(json_files, output_prefix, compute_sha256)

if __name__ == "__main__":
    main()
