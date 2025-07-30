import sys
import pytest
import numpy as np
from ardal import ArdalParser


def test_parse_npy_json(test_data):
    """
    Tests parsing of NPY/JSON file pairs.
    """
    matrix, headers = test_data
    parser = ArdalParser(input_file_structure=[headers, matrix])
    assert parser.matrix is not None
    assert parser.headers is not None


def test_invalid_format():
    """
    Tests handling of invalid file formats.
    """
    with pytest.raises(ValueError):
        ArdalParser(input_file_structure="test.invalid", file_format="invalid")


def test_invalid_list():
    """
    Tests handling of invalid file lists.
    """
    with pytest.raises(ValueError):
        ArdalParser(input_file_structure=["test.json"])


def test_invalid_list_types():
    """
    Tests handling of invalid file list types.
    """
    with pytest.raises(ValueError):
        ArdalParser(input_file_structure=[1, 2])
