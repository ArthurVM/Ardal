import pytest
import numpy as np
import json
import os
from tempfile import NamedTemporaryFile, TemporaryDirectory

from ardal.core.ArdalParser import ArdalParser
from ardal.exceptions.exceptions import MalformedInputError, LoadMatrixError, UnsupportedFormatError


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

def test_valid_npy_db():
    npy_file = "../data/sim_matrix.npy"
    header_file = "../data/sim_headers.json"
    npy_array = np.load(npy_file)
    with open(header_file, 'r') as fin:
        headers = json.load(fin)
    parser = ArdalParser([npy_file, header_file])
    assert (parser.matrix == npy_array).all()
    assert parser.headers == headers

def test_valid_npy_db():
    npz_file = "../data/sim_matrix.npz"
    header_file = "../data/sim_headers.json"
    npz_array = np.load(npz_file)['matrix']
    with open(header_file, 'r') as fin:
        headers = json.load(fin)
    parser = ArdalParser([npz_file, header_file])
    assert (parser.matrix == npz_array).all()
    assert parser.headers == headers

def test_valid_csv_db():
    import pandas as pd
    csv_file = "../data/sim_matrix.csv"
    df = pd.read_csv(csv_file, index_col=0)
    matrix = df.values
    headers = {"guids": list(df.index), "alleles": list(df.columns)}
    parser = ArdalParser(csv_file)
    assert (parser.matrix == matrix).all()
    assert parser.headers == headers