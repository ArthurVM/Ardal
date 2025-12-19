""" make_meta.py
Centralised helpers for constructing matrix metadata for the Ardal framework.
"""
from typing import Dict, Any
import numpy as np

from .exceptions import LoadMatrixError


def make_meta( matrix: np.ndarray,
               headers: Dict[str, Any],
               *,
               generated_by: str,
               format_name: str,
               matrix_file: str | None = None ) -> Dict[str, Any]:
    """
    Construct metadata for either a dense or bitpacked matrix.

    Rules:
      - Matrix must be 2D.
      - Headers must contain 'guids' (list[str]) and 'alleles' (list[str]).
      - If matrix dtype is uint64 and shape[1] matches words_per_row derived from alleles,
        treat as bitpacked; otherwise treat as dense.
    """
    if not isinstance(headers, dict):
        raise LoadMatrixError("Headers must be a Dictionary.")
    if "guids" not in headers or "alleles" not in headers:
        raise LoadMatrixError("Headers must contain 'guids' and 'alleles' keys.")
    guids = headers["guids"]
    alleles = headers["alleles"]
    if not isinstance(guids, list) or not all(isinstance(g, str) for g in guids):
        raise LoadMatrixError("GUIDs must be a list of strings.")
    if not isinstance(alleles, list) or not all(isinstance(a, str) for a in alleles):
        raise LoadMatrixError("Alleles must be a list of strings.")
    if matrix.ndim != 2:
        raise LoadMatrixError("Matrix must be 2-dimensional.")

    n_rows = int(matrix.shape[0])
    n_cols = len(alleles)

    ## detect packed
    is_uint64 = matrix.dtype == np.dtype("<u8") or matrix.dtype == np.dtype("uint64")
    words_per_row = (n_cols + 63) // 64 if n_cols else 0
    looks_packed = is_uint64 and matrix.shape[1] == words_per_row

    if looks_packed:
        ## verify words_per_row matches header derived expectation
        if matrix.shape[1] != words_per_row:
            raise LoadMatrixError(
                f"Bitpack matrix words_per_row mismatch: matrix has {matrix.shape[1]}, expected {words_per_row} for {n_cols} columns."
            )
        return {
            "format": "ardal.bitpack.v1",
            "dtype": "<u8",
            "endianness": "little",
            "row_major": True,
            "n_rows": n_rows,
            "n_cols": n_cols,
            "words_per_row": int(words_per_row),
            "bits_per_word": 64,
            "row_stride_bytes": int(words_per_row * 8),
            "data_file": matrix_file,
            "data_file_resolved": None,
            "data_nbytes": int(matrix.nbytes),
            "data_sha256": None,
            "generated_by": generated_by,
        }

    ## dense path
    endianness = "little" if matrix.dtype.byteorder in ("<", "=", "|") else "big"
    stride = int(matrix.strides[0]) if matrix.ndim >= 2 else None
    return {
        "format": "ardal.dense.v1",
        "dtype": str(matrix.dtype),
        "endianness": endianness,
        "row_major": True,
        "n_rows": n_rows,
        "n_cols": n_cols,
        "words_per_row": None,
        "bits_per_word": None,
        "row_stride_bytes": stride,
        "data_file": matrix_file,
        "data_file_resolved": None,
        "data_nbytes": int(matrix.nbytes),
        "data_sha256": None,
        "generated_by": generated_by,
    }