import pytest
import numpy as np
import json
import gzip
import os
from tempfile import NamedTemporaryFile, TemporaryDirectory

from ardal.core.ArdalParser import ArdalParser
from ardal.utils.exceptions import MalformedInputError, LoadMatrixError, UnsupportedFormatError


def valid_headers(n_rows, n_cols):
    return {
        "guids": [f"s{i}" for i in range(n_rows)],
        "alleles": [f"a{j}" for j in range(n_cols)]
    }


## GENERAL TEST

def test_rejects_invalid_list_length():
    with pytest.raises(MalformedInputError):
        ArdalParser([np.ones((3, 5))])

def test_rejects_non_matching_shapes():
    matrix = np.ones((3, 5), dtype=np.uint8)
    headers = valid_headers(2, 5)  # wrong number of guids
    with pytest.raises(LoadMatrixError):
        ArdalParser([matrix, headers])

def test_rejects_nonunique_guids():
    matrix = np.ones((2, 5), dtype=np.uint8)
    headers = {"guids": ["s1", "s1"], "alleles": [f"a{i}" for i in range(5)]}
    with pytest.raises(LoadMatrixError):
        ArdalParser([matrix, headers])

def test_rejects_nonunique_alleles():
    matrix = np.ones((2, 5), dtype=np.uint8)
    headers = {"guids": ["s1", "s2"], "alleles": ["a", "a", "b", "c", "d"]}
    with pytest.raises(LoadMatrixError):
        ArdalParser([matrix, headers])

## TEST IN MEMORY
def test_accepts_valid_memory_input():
    matrix = np.ones((3, 5), dtype=np.uint8)
    headers = valid_headers(3, 5)
    parser = ArdalParser([matrix, headers])
    assert parser.matrix.shape == (3, 5)
    assert parser.headers == headers

def test_accepts_dict_matrix_order():
    matrix = np.ones((2, 4), dtype=np.uint8)
    headers = valid_headers(2, 4)
    parser = ArdalParser([headers, matrix])
    assert parser.matrix.shape == (2, 4)

def test_rejects_none_input():
    with pytest.raises(MalformedInputError):
        ArdalParser(None)

## TEST ON DISK
def test_rejects_missing_file():
    with pytest.raises(FileNotFoundError):
        ArdalParser(["nonexistent.json", "nonexistent.npy"])

def test_rejects_unsupported_format():
    with NamedTemporaryFile(suffix=".foo") as f:
        f.write(b"bad format")
        f.flush()
        with pytest.raises(UnsupportedFormatError):
            ArdalParser(f.name)

def test_valid_npy_db(test_data_matrix_npy):
    npy_file, header_file = test_data_matrix_npy
    npy_array = np.load(npy_file)
    with open(header_file, 'r') as fin:
        headers = json.load(fin)
    parser = ArdalParser([npy_file, header_file])
    assert (parser.matrix == npy_array).all()
    assert parser.headers == headers

def test_valid_npz_db(test_data_matrix_npz):
    npz_file, header_file = test_data_matrix_npz
    npz_array = np.load(npz_file)['matrix']
    with open(header_file, 'r') as fin:
        headers = json.load(fin)
    parser = ArdalParser([npz_file, header_file])
    assert (parser.matrix == npz_array).all()
    assert parser.headers == headers

def test_valid_csv_db(test_data_matrix_csv):
    import pandas as pd
    csv_file = test_data_matrix_csv
    df = pd.read_csv(csv_file, index_col=0)
    matrix = df.values
    headers = {"guids": list(df.index), "alleles": list(df.columns)}
    parser = ArdalParser(csv_file)
    assert (parser.matrix == matrix).all()
    assert parser.headers == headers


def test_valid_npy_with_gzipped_meta(tmp_path, test_data_matrix_npy):
    npy_file, header_file = test_data_matrix_npy
    npy_array = np.load(npy_file)
    with open(header_file, "r") as fin:
        headers = json.load(fin)

    meta_path = tmp_path / "sim_matrix.npy.meta.gz"
    payload = {
        "meta": {
            "format": "ardal.npy.v1",
            "dtype": str(npy_array.dtype),
            "endianness": "little",
            "row_major": True,
            "n_rows": npy_array.shape[0],
            "n_cols": npy_array.shape[1],
            "matrix_file": os.path.basename(npy_file),
            "data_nbytes": int(npy_array.nbytes),
            "data_sha256": None,
            "generated_by": "test",
        },
        "headers": headers,
    }
    with gzip.open(meta_path, "wt") as fout:
        json.dump(payload, fout)

    parser = ArdalParser([npy_file, str(meta_path)])
    assert (parser.matrix == npy_array).all()
    assert parser.headers == headers


def test_valid_bin_zst_with_gzipped_meta(tmp_path, test_data_mem):
    zstandard = pytest.importorskip("zstandard")

    matrix, headers = test_data_mem
    packed_bytes = np.packbits(matrix, axis=1, bitorder="little")
    words_per_row = (matrix.shape[1] + 63) // 64
    packed = np.zeros((matrix.shape[0], words_per_row), dtype="<u8")
    packed.view(np.uint8).reshape(matrix.shape[0], words_per_row * 8)[:, :packed_bytes.shape[1]] = packed_bytes

    bin_zst_path = tmp_path / "sim_matrix.bin.zst"
    compressed = zstandard.ZstdCompressor().compress(packed.tobytes(order="C"))
    bin_zst_path.write_bytes(compressed)

    meta_path = tmp_path / "sim_matrix.bin.meta.gz"
    payload = {
        "meta": {
            "format": "ardal.bitpack.v1",
            "dtype": "<u8",
            "endianness": "little",
            "row_major": True,
            "n_rows": matrix.shape[0],
            "n_cols": matrix.shape[1],
            "matrix_file": "sim_matrix.bin",
            "data_nbytes": int(packed.nbytes),
            "data_sha256": None,
            "words_per_row": words_per_row,
            "bits_per_word": 64,
            "row_stride_bytes": words_per_row * 8,
            "generated_by": "test",
        },
        "headers": headers,
    }
    with gzip.open(meta_path, "wt") as fout:
        json.dump(payload, fout)

    parser = ArdalParser([str(bin_zst_path), str(meta_path)])
    np.testing.assert_array_equal(parser.matrix, packed)
    assert parser.headers == headers
    assert parser.is_bitpacked


def test_valid_bin_v2_with_binary_missing_sections(tmp_path):
    headers = {
        "guids": ["s0", "s1"],
        "alleles": [f"chr1.{i}.A.T" for i in range(10)],
    }
    packed = np.array([[0b10101], [0b01010]], dtype="<u8")
    bin_path = tmp_path / "matrix.bin"
    bin_path.write_bytes(packed.tobytes(order="C"))

    offsets_offset = bin_path.stat().st_size
    offsets = np.array([0, 2, 2], dtype="<u8")
    ranges_offset = offsets_offset + offsets.nbytes
    ranges = np.array([[1, 3], [7, 8]], dtype="<u4")
    with bin_path.open("ab") as fout:
        fout.write(offsets.tobytes(order="C"))
        fout.write(ranges.tobytes(order="C"))

    meta_path = tmp_path / "matrix.bin.meta"
    payload = {
        "meta": {
            "format": "ardal.bin.v2",
            "dtype": "<u8",
            "endianness": "little",
            "row_major": True,
            "n_rows": 2,
            "n_cols": 10,
            "matrix_file": "matrix.bin",
            "data_nbytes": int(bin_path.stat().st_size),
            "data_sha256": None,
            "words_per_row": 1,
            "bits_per_word": 64,
            "row_stride_bytes": 8,
            "sections": {
                "allele_matrix": {
                    "offset": 0,
                    "nbytes": int(packed.nbytes),
                    "dtype": "<u8",
                    "shape": [2, 1],
                    "words_per_row": 1,
                    "bits_per_word": 64,
                    "row_stride_bytes": 8,
                },
                "missing_offsets": {
                    "offset": int(offsets_offset),
                    "nbytes": int(offsets.nbytes),
                    "dtype": "<u8",
                    "length": 3,
                    "units": "range_rows",
                },
                "missing_ranges": {
                    "offset": int(ranges_offset),
                    "nbytes": int(ranges.nbytes),
                    "dtype": "<u4",
                    "shape": [2, 2],
                    "coordinate_system": "matrix_columns_0_based",
                    "range_semantics": "start_inclusive_end_exclusive",
                },
            },
        },
        "headers": headers,
        "missing": {
            "schema": "ardal.missing_column_ranges.v1",
            "encoding": "binary_sections",
            "coordinate_system": "matrix_columns_0_based",
            "range_semantics": "start_inclusive_end_exclusive",
            "offsets_section": "missing_offsets",
            "ranges_section": "missing_ranges",
        },
    }
    meta_path.write_text(json.dumps(payload))

    parser = ArdalParser([str(bin_path), str(meta_path)])

    np.testing.assert_array_equal(parser.matrix, packed)
    assert parser.headers == headers
    assert parser.is_bitpacked
    assert parser.missing_masks.has_missing_mask()
    assert parser.missing_masks.get_guid_missing_mask("s0") == [1, 2, 7]
    assert parser.missing_masks.get_guid_missing_mask("s1") == []
    assert parser.missing_masks.get_missing_mask_rows() == [[1, 2, 7], []]
    backend_payload = parser.missing_masks.to_backend_payload()
    assert backend_payload["encoding"] == "range_sections"
    np.testing.assert_array_equal(backend_payload["offsets"], offsets)
    np.testing.assert_array_equal(backend_payload["ranges"], ranges)


def test_valid_bin_v3_with_compressed_allele_section_and_uncompressed_missing_sections(tmp_path):
    zstandard = pytest.importorskip("zstandard")

    headers = {
        "guids": ["s0", "s1"],
        "alleles": [f"chr1.{i}.A.T" for i in range(10)],
    }
    packed = np.array([[0b10101], [0b01010]], dtype="<u8")
    compressed_alleles = zstandard.ZstdCompressor().compress(packed.tobytes(order="C"))

    bin_path = tmp_path / "matrix.bin"
    bin_path.write_bytes(compressed_alleles)

    offsets_offset = bin_path.stat().st_size
    offsets = np.array([0, 1, 2], dtype="<u8")
    ranges_offset = offsets_offset + offsets.nbytes
    ranges = np.array([[1, 3], [7, 8]], dtype="<u4")
    with bin_path.open("ab") as fout:
        fout.write(offsets.tobytes(order="C"))
        fout.write(ranges.tobytes(order="C"))

    meta_path = tmp_path / "matrix.bin.meta"
    payload = {
        "meta": {
            "format": "ardal.bin.v3",
            "dtype": "<u8",
            "endianness": "little",
            "row_major": True,
            "n_rows": 2,
            "n_cols": 10,
            "matrix_file": "matrix.bin",
            "data_nbytes": int(bin_path.stat().st_size),
            "data_sha256": None,
            "words_per_row": 1,
            "bits_per_word": 64,
            "row_stride_bytes": 8,
            "sections": {
                "allele_matrix": {
                    "offset": 0,
                    "nbytes": int(len(compressed_alleles)),
                    "compression": "zstd",
                    "uncompressed_nbytes": int(packed.nbytes),
                    "dtype": "<u8",
                    "shape": [2, 1],
                    "words_per_row": 1,
                    "bits_per_word": 64,
                    "row_stride_bytes": 8,
                },
                "missing_offsets": {
                    "offset": int(offsets_offset),
                    "nbytes": int(offsets.nbytes),
                    "compression": None,
                    "dtype": "<u8",
                    "length": 3,
                    "units": "range_rows",
                },
                "missing_ranges": {
                    "offset": int(ranges_offset),
                    "nbytes": int(ranges.nbytes),
                    "compression": None,
                    "dtype": "<u4",
                    "shape": [2, 2],
                    "coordinate_system": "matrix_columns_0_based",
                    "range_semantics": "start_inclusive_end_exclusive",
                },
            },
        },
        "headers": headers,
        "missing": {
            "schema": "ardal.missing_column_ranges.v1",
            "encoding": "binary_sections",
            "coordinate_system": "matrix_columns_0_based",
            "range_semantics": "start_inclusive_end_exclusive",
            "offsets_section": "missing_offsets",
            "ranges_section": "missing_ranges",
        },
    }
    meta_path.write_text(json.dumps(payload))

    parser = ArdalParser([str(bin_path), str(meta_path)])

    np.testing.assert_array_equal(parser.matrix, packed)
    assert parser.headers == headers
    assert parser.is_bitpacked
    assert parser.missing_masks.has_missing_mask()
    assert parser.missing_masks.get_guid_missing_mask("s0") == [1, 2]
    assert parser.missing_masks.get_guid_missing_mask("s1") == [7]
