import pytest
import pandas as pd
import numpy as np
import json
import os
from Bio import SeqIO

from ardal import Ardal


def test_to_dataframe(io_component, test_data_mem):
    """
    Tests the toDataFrame method of the ArdalIO class.
    """
    df = io_component.to_dataframe()
    matrix, headers = test_data_mem

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (10, 10)
    assert df.index.tolist() == headers["guids"]
    assert df.columns.tolist() == headers["alleles"]
    np.testing.assert_array_equal(df.values, matrix)


def test_to_dict(io_component):
    """
    Tests the toDict method of the ArdalIO class.
    """
    data_dict = io_component.to_dict()
    assert isinstance(data_dict, dict)
    ## based on test_data in conftest.py
    ## GUID1 has SNP1, SNP5, SNP6, SNP7, SNP8, SNP10
    assert set(data_dict["GUID1"]) == {"chr1.1.A.T", "chr1.5.A.T", "chr1.6.A.T", "chr1.7.A.T", "chr1.8.A.T", "chr1.10.A.T"}
    ## GUID10 has SNP4, SNP5, SNP7, SNP9
    assert set(data_dict["GUID10"]) == {"chr1.4.A.T", "chr1.7.A.T", "chr1.9.A.T", "chr1.5.A.T"}


def test_write_npy(io_component, test_data_mem, tmp_path):
    """
    Tests the write method with npy output.
    """
    matrix, headers = test_data_mem
    output_dir = tmp_path / "test_output_npy"
    output_dir.mkdir()
    io_component.write("test_data", str(output_dir), format="npy")

    ## verify files were created
    json_file = output_dir / "test_data.json"
    npy_file = output_dir / "test_data.npy"
    assert json_file.exists()
    assert npy_file.exists()

    ## verify content
    with open(json_file, 'r') as f:
        loaded_headers = json.load(f)
    assert loaded_headers == headers

    loaded_matrix = np.load(npy_file)
    np.testing.assert_array_equal(loaded_matrix, matrix)


def test_write_npz(io_component, test_data_mem, tmp_path):
    """
    Tests the write method with npz output.
    """
    matrix, _ = test_data_mem
    output_dir = tmp_path / "test_output_npz"
    output_dir.mkdir()
    io_component.write("test_data", str(output_dir), format="npz")

    npz_file = output_dir / "test_data.npz"
    assert npz_file.exists()

    loaded_data = np.load(npz_file)
    assert "matrix" in loaded_data
    np.testing.assert_array_equal(loaded_data["matrix"], matrix)
    

def test_write_bin(io_component, test_data_mem, tmp_path):
    """
    Tests the write method with npz output.
    """
    matrix, _ = test_data_mem
    output_dir = tmp_path / "test_output_bin"
    output_dir.mkdir()
    io_component.write("test_data", str(output_dir), format="bin")

    bin_file = output_dir / "test_data.bin"
    assert bin_file.exists()
    
    json_file = output_dir / "test_data.json"
    assert json_file.exists()

    loaded_json = json.load(open(json_file, 'r'))
    assert "meta" in loaded_json
    assert "headers" in loaded_json


def test_write_metadata_shape(io_component, tmp_path):
    output_dir = tmp_path / "test_output_meta_shape"
    output_dir.mkdir()

    io_component.write("test_data", str(output_dir), format="bin")
    meta_file = output_dir / "test_data.bin.meta"
    assert meta_file.exists()

    with open(meta_file, "r") as f:
        payload = json.load(f)

    assert set(payload.keys()) == {"meta", "headers", "column_masks"}
    assert set(payload["headers"].keys()) == {"guids", "alleles"}
    assert isinstance(payload["column_masks"], dict)

    expected_meta_keys = {
        "format",
        "dtype",
        "endianness",
        "row_major",
        "n_rows",
        "n_cols",
        "matrix_file",
        "data_nbytes",
        "data_sha256",
        "words_per_row",
        "bits_per_word",
        "row_stride_bytes",
        "generated_by",
    }
    assert set(payload["meta"].keys()) == expected_meta_keys


def test_write_metadata_column_masks_format(hybrid_matrix, headers, meta, tmp_path):
    from ardal.core.ArdalHeaderUtils import ArdalHeaderUtils
    from ardal.core.ArdalIO import ArdalIO

    masks = {"column_masks": {"GUID1": [0, 7], "GUID2": [0], "GUID10": []}}
    header_utils = ArdalHeaderUtils(
        headers=headers,
        meta=meta,
        allele_id_format="{chr}.{start}.{ref}.{alt}",
        missing_masks=masks,
    )
    io_component = ArdalIO(header_utils, hybrid_matrix, hybrid_matrix.roaringEnabled())

    output_dir = tmp_path / "test_output_meta_masks"
    output_dir.mkdir()
    io_component.write("test_data", str(output_dir), format="npy")

    meta_file = output_dir / "test_data.npy.meta"
    with open(meta_file, "r") as f:
        payload = json.load(f)

    assert payload["column_masks"]["GUID1"] == [0, 7]
    assert payload["column_masks"]["GUID2"] == [0]
    assert payload["column_masks"]["GUID10"] == []


def test_make_alignment_non_snp_uses_reference_sequence(tmp_path):
    matrix = np.asarray(
        [
            [1, 0],
            [0, 1],
        ],
        dtype=np.uint8,
    )
    headers = {
        "guids": ["GUID1", "GUID2"],
        "alleles": ["chr1.2.A.T", "chr1.5.A.G"],
    }

    ardal = Ardal(
        data_source=[matrix, headers],
        allele_id_format="{chr}.{start}.{ref}.{alt}",
        quiet_init=True,
    )

    ref_path = tmp_path / "ref.fa"
    ref_path.write_text(">chr1\nAAAAAA\n")

    out_path = ardal.io.make_alignment(
        output_prefix="alignment",
        ref=str(ref_path),
        allele_id_format="{chr}.{start}.{ref}.{alt}",
        out_directory=str(tmp_path),
        snp_only=False,
    )

    records = list(SeqIO.parse(out_path, "fasta"))
    assert [record.id for record in records] == ["GUID1", "GUID2"]
    assert str(records[0].seq) == "ATAAAA"
    assert str(records[1].seq) == "AAAAGA"
