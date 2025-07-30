import pytest
import pandas as pd
import numpy as np
import json
import os


def test_to_dataframe(io_component, test_data):
    """
    Tests the toDataFrame method of the ArdalIO class.
    """
    df = io_component.toDataFrame()
    matrix, headers = test_data

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (10, 10)
    assert df.index.tolist() == headers["guids"]
    assert df.columns.tolist() == headers["alleles"]
    np.testing.assert_array_equal(df.values, matrix)


def test_to_dict(io_component):
    """
    Tests the toDict method of the ArdalIO class.
    """
    data_dict = io_component.toDict()
    assert isinstance(data_dict, dict)
    # Based on test_data in conftest.py
    # GUID1 has SNP1, SNP5, SNP6, SNP7, SNP8, SNP10
    assert set(data_dict["GUID1"]) == {"SNP1", "SNP5", "SNP6", "SNP7", "SNP8", "SNP10"}
    # GUID10 has SNP4, SNP5, SNP7, SNP9
    assert set(data_dict["GUID10"]) == {"SNP4", "SNP5", "SNP7", "SNP9"}


def test_write_npy(io_component, headers, test_data, tmp_path):
    """
    Tests the write method with npy output.
    """
    matrix, _ = test_data
    output_dir = tmp_path / "test_output_npy"
    output_dir.mkdir()
    io_component.write(str(output_dir), "test_data", npz=False)

    # Verify files were created
    json_file = output_dir / "test_data_headers.json"
    npy_file = output_dir / "test_data_matrix.npy"
    assert json_file.exists()
    assert npy_file.exists()

    # Verify content
    with open(json_file, 'r') as f:
        loaded_headers = json.load(f)
    assert loaded_headers == headers

    loaded_matrix = np.load(npy_file)
    np.testing.assert_array_equal(loaded_matrix, matrix)


def test_write_npz(io_component, test_data, tmp_path):
    """
    Tests the write method with npz output.
    """
    matrix, _ = test_data
    output_dir = tmp_path / "test_output_npz"
    output_dir.mkdir()
    io_component.write(str(output_dir), "test_data", npz=True)

    npz_file = output_dir / "test_data_matrix.npz"
    assert npz_file.exists()

    loaded_data = np.load(npz_file)
    assert "matrix" in loaded_data
    np.testing.assert_array_equal(loaded_data["matrix"], matrix)