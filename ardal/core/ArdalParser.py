# core/ArdalParser.py
import os
import csv
import gzip
import io
import json
import hashlib
import tempfile
from sys import byteorder
from pathlib import Path
from typing import Union, Tuple, Dict
from importlib.metadata import version, PackageNotFoundError

import pandas as pd
import numpy as np

from ..utils.misc import require_package
from ..utils.exceptions import MalformedInputError, UnsupportedFormatError, LoadMatrixError
from ..utils.logger import get_logger
from ..utils.make_meta import make_meta

try:
    from _ardal_pack import pack_dense_to_words as _pack_words_cpp
except Exception:
    _pack_words_cpp = None
    
    
log = get_logger()


class BinaryMissingMaskStore:
    """Lazy view over v2 binary missing-column range sections."""

    def __init__(self, guids, offsets, ranges):
        self.guids = list(guids)
        self._guid_to_idx = {guid: idx for idx, guid in enumerate(self.guids)}
        self.offsets = offsets
        self.ranges = ranges
        self._cache = {}

    def has_missing_mask(self) -> bool:
        return bool(self.ranges.shape[0])

    def get_guid_missing_mask(self, guid: str) -> list:
        if guid in self._cache:
            return self._cache[guid]
        row_idx = self._guid_to_idx.get(guid)
        if row_idx is None:
            return []
        start = int(self.offsets[row_idx])
        end = int(self.offsets[row_idx + 1])
        cols = []
        for range_start, range_end in self.ranges[start:end]:
            cols.extend(range(int(range_start), int(range_end)))
        self._cache[guid] = cols
        return cols

    def get_missing_mask_rows(self) -> list:
        return [self.get_guid_missing_mask(guid) for guid in self.guids]

    def to_backend_payload(self) -> dict:
        return {
            "encoding": "range_sections",
            "offsets": self.offsets,
            "ranges": self.ranges,
        }


class ArdalParser:
    """
    Parses Ardal inputs:
      - Dense: CSV, Parquet, NPY/JSON, NPZ/JSON
      - Bitpacked: BIN/JSON, BIN.ZST/JSON

    Returns:
      self.matrix  -> np.ndarray
        - dense: shape (n_rows, n_cols), dtype arbitrary
        - bitpack: memmap '<u8' or in-memory '<u8', shape (n_rows, words_per_row)
      self.headers -> Dict with "guids", "alleles"
      self.meta    -> Dict (bitpack only), else {}
      self.is_bitpacked -> bool
    """

    def __init__( self,
                  input_data_structure: Union[list, str],
                  verify_hash: bool = False,
                  is_packed_mem: bool = False ):
        self.input_data = input_data_structure
        self.file_format = None
        self.matrix: np.ndarray | None = None
        self.headers: Dict = {}
        self.meta: Dict = {}
        self.missing_masks: Dict = {"missing_masks": {}}
        self.is_bitpacked: bool = False
        self.verify_hash = verify_hash
        self.is_packed_mem = is_packed_mem
        self._temporary_matrix_files: list[str] = []
        self._parse()

    def __del__( self ):
        for path in getattr(self, "_temporary_matrix_files", []):
            try:
                Path(path).unlink(missing_ok=True)
            except Exception:
                pass



    ## ------------- parsing -------------

    def _parse( self ) -> Union[int, None]:
        
        if self.input_data is None:
            raise MalformedInputError("Input data structure cannot be None.")

        ## two objects in a list
        ## e.g. [np.ndarray, Dict] or [Path, Path]
        if isinstance(self.input_data, list):
            
            if len(self.input_data) != 2:
                raise MalformedInputError("Input list must contain two elements: matrix and headers (or two file paths).")

            a, b = self.input_data


            ## in-memory [array, headers] or [headers, array]
            if (isinstance(a, np.ndarray) and isinstance(b, Dict)) or (isinstance(b, np.ndarray) and isinstance(a, Dict)):
                matrix, headers_meta_raw = self._order_input(a, b)
                log.info(f"Parsing matrix data from list in memory of form: {[type(i) for i in self.input_data]}")
                headers_clean, missing_obj, meta_obj = self._load_headers_dict(headers_meta_raw)
                self.headers = headers_clean
                self.missing_masks = self._normalise_missing_masks(missing_obj, self.headers)

                want_bitpack = self.is_packed_mem or self._is_bitpacked_candidate(matrix)
                # print(want_bitpack)
                self.matrix = np.ascontiguousarray(matrix)

                if want_bitpack:
                    if not self._is_bitpacked_candidate(self.matrix):
                        raise LoadMatrixError("Bitpack matrix must be a 2D uint64 array.")
                    if self.matrix.dtype != np.dtype("<u8"):
                        if self.matrix.dtype.kind == "u" and self.matrix.dtype.itemsize == 8:
                            self.matrix = self.matrix.astype(np.dtype("<u8"), copy=False)
                        else:
                            raise LoadMatrixError(f"Bitpack matrix has incompatible dtype: {self.matrix.dtype}")
                    self.file_format = "bitpack"
                    self.is_bitpacked = True
                    self.meta = meta_obj if meta_obj else self._ensure_meta(self.matrix, self.headers, None, is_bitpacked=True)
                    # self._validate_bitpack(self.matrix, headers_clean, meta_obj)
                    return 0

                self.file_format = "memory_npy"
                self.is_bitpacked = False
                self._validate_dense()
                self.meta = meta_obj if meta_obj else self._ensure_meta(self.matrix, self.headers, None, is_bitpacked=False)
                return 0


            ## two file paths
            if isinstance(a, str) and isinstance(b, str):
                log.info(f"Parsing matrix data from two file paths: {self.input_data}")
                
                a, b = Path(a), Path(b)
                if not a.exists() or not b.exists():
                    raise FileNotFoundError(f"One or more file paths do not exist: {a}, {b}")
                meta_path = a if self._is_meta_path(a) else b
                mat_path  = b if meta_path == a else a
                matrix_format = self._matrix_path_format(mat_path)

                ## bitpack pair
                ## 1. load pair
                ## 2. validate data structures
                ## 3. validate missing structures
                ## 4. define file format
                ## 5. construct variables
                if matrix_format in (".bin", ".bin.zst"):
                    log.info(f"Detected matrix format: {matrix_format}")
                    matrix, headers, meta, missing_raw = self._load_bitpack_pair(meta_path, mat_path)
                    self._validate_bitpack(matrix, headers, meta)
                    self.missing_masks = self._normalise_missing_masks(missing_raw, headers)
                    self.matrix = matrix
                    self.headers = headers
                    self.meta = self._ensure_meta(matrix, headers, meta, is_bitpacked=True)
                    self.file_format = "bitpack"
                    self.is_bitpacked = True
                    return 0

                ## dense pairs
                if matrix_format == ".npy":
                    log.info(f"Detected matrix format: .npy")
                    matrix, headers, meta, missing_raw = self._load_npy_pair(meta_path, mat_path)
                    self._validate_dense_pair(matrix, headers)
                    self.missing_masks = self._normalise_missing_masks(missing_raw, headers)
                    self.matrix = matrix
                    self.headers = headers
                    self.meta = self._ensure_meta(matrix, headers, meta, is_bitpacked=False)
                    self.file_format = "npy"
                    return 0

                if matrix_format == ".npz":
                    log.info(f"Detected matrix format: .npz")
                    matrix, headers, meta, missing_raw = self._load_npz_pair(meta_path, mat_path)
                    self._validate_dense_pair(matrix, headers)
                    self.missing_masks = self._normalise_missing_masks(missing_raw, headers)
                    self.matrix = matrix
                    self.headers = headers
                    self.meta = self._ensure_meta(matrix, headers, meta, is_bitpacked=False)
                    self.file_format = "npz"
                    return 0

                raise UnsupportedFormatError(f"Unrecognized file pair: {a.name}, {b.name}")

            raise MalformedInputError(
                "If list input, it must be [headers::Dict, matrix::np.ndarray] or "
                "[headers.json, matrix.bin|matrix.bin.zst] or dense pairs."
            )

        ## single object
        if isinstance(self.input_data, str):
            log.info(f"Parsing matrix data from one file path: {self.input_data}")
            
            p = Path(self.input_data)
            suffix = p.suffix.lower()
            supported_suffixes = {".csv"}
            if suffix not in supported_suffixes:
                raise UnsupportedFormatError(f"Unsupported file extension: {p.suffix}")
            if not p.exists():
                raise FileNotFoundError(f"File does not exist: {p}")

            ## csv
            if suffix == ".csv":
                log.info(f"Detected matrix format: .csv")
                matrix, headers = self._load_csv(str(p))
                headers_clean, missing_obj, meta_obj = self._load_headers_dict(headers)
                self.matrix, self.headers = matrix, headers_clean
                self.meta = meta_obj if meta_obj else self._ensure_meta(self.matrix, self.headers, None, is_bitpacked=False)
                self.file_format = "csv"
                self.is_bitpacked = False
                self.missing_masks = self._normalise_missing_masks(missing_obj, self.headers)
                return 0

            if suffix in (".npy", ".npz", ".bin", ".zst"):
                raise MalformedInputError("Binary file provided without matching headers.json; provide both.")

            raise UnsupportedFormatError(f"Unsupported file extension: {p.suffix}")

        raise MalformedInputError("Input must be list or string.")


    def _load_headers_dict( self,
                            headers_meta_raw: Dict ) -> Tuple[Dict, Union[Dict, None], Union[Dict, None]]:
        """
        Validate an in-memory headers dictionary and extract missing-mask and metadata payloads.
        Expected shapes:
          1) {\"headers\": {...}, \"column_masks\": ..., \"meta\": ...}
          2) {\"guids\": [...], \"alleles\": [...], \"column_masks\": ..., \"meta\": ...}
        """
        if not isinstance(headers_meta_raw, Dict):
            raise LoadMatrixError("Headers metadata must be a Dictionary.")

        # detect payload layout
        if "headers" in headers_meta_raw and isinstance(headers_meta_raw["headers"], Dict):
            headers_raw = headers_meta_raw["headers"]
            if "column_masks" in headers_meta_raw:
                missing_raw = {"column_masks": headers_meta_raw.get("column_masks")}
            else:
                missing_raw = headers_meta_raw.get("missing")
            meta_raw = headers_meta_raw.get("meta")
        else:
            headers_raw = headers_meta_raw
            if "column_masks" in headers_meta_raw:
                missing_raw = {"column_masks": headers_meta_raw.get("column_masks")}
            else:
                missing_raw = headers_meta_raw.get("missing")
            meta_raw = headers_meta_raw.get("meta")

        if "guids" not in headers_raw or "alleles" not in headers_raw:
            raise LoadMatrixError("Headers must contain 'guids' and 'alleles' keys.")

        guids = headers_raw["guids"]
        alleles = headers_raw["alleles"]

        if not isinstance(guids, list) or not all(isinstance(g, str) for g in guids):
            raise LoadMatrixError("GUIDs must be a list of strings.")
        if not isinstance(alleles, list) or not all(isinstance(a, str) for a in alleles):
            raise LoadMatrixError("Alleles must be a list of strings.")
        if len(set(guids)) != len(guids):
            raise LoadMatrixError("GUIDs must be unique.")
        if len(set(alleles)) != len(alleles):
            raise LoadMatrixError("Alleles must be unique.")

        headers_clean = {"guids": guids, "alleles": alleles}
        return headers_clean, missing_raw, meta_raw


    ## ------------- validators -------------

    def _validate_dense( self ) -> None:
        if not isinstance(self.matrix, np.ndarray):
            raise LoadMatrixError("Matrix must be a NumPy array.")
        if self.matrix.ndim != 2:
            raise LoadMatrixError("Matrix must be 2-dimensional.")
        if not self.matrix.flags['C_CONTIGUOUS']:
            try:
                self.matrix = np.ascontiguousarray(self.matrix)
            except Exception as e:
                raise LoadMatrixError(f"Failed to make array contiguous in memory : {e}")

        n_rows, n_cols_bits = self.matrix.shape

        if not isinstance(self.headers, Dict):
            raise LoadMatrixError("Headers must be a Dictionary.")
        if "guids" not in self.headers or "alleles" not in self.headers:
            raise LoadMatrixError("Headers must contain 'guids' and 'alleles' keys.")

        guids = self.headers["guids"]
        alleles = self.headers["alleles"]

        if not isinstance(guids, list) or not all(isinstance(g, str) for g in guids):
            raise LoadMatrixError("GUIDs must be a list of strings.")
        if not isinstance(alleles, list) or not all(isinstance(a, str) for a in alleles):
            raise LoadMatrixError("Alleles must be a list of strings.")

        if len(guids) != n_rows:
            raise LoadMatrixError(f"Mismatch: {n_rows} matrix rows vs {len(guids)} GUIDs.")
        if len(alleles) != n_cols_bits:
            raise LoadMatrixError(f"Mismatch: {n_cols_bits} matrix columns vs {len(alleles)} alleles.")

        if len(set(guids)) != len(guids):
            raise LoadMatrixError("GUIDs must be unique.")
        if len(set(alleles)) != len(alleles):
            raise LoadMatrixError("Alleles must be unique.")


    def _validate_bitpack( self,
                           mmap_arr : np.memmap,
                           headers_raw : Dict,
                           meta_raw : Dict ) -> None:
        supported_formats = ["ardal.dense.v1", "ardal.bitpack.v1", "ardal.bin.v1", "ardal.bin.v2", "ardal.bin.v3", "ardal.npy.v1", "ardal.npz.v1"]
        
        if not isinstance(mmap_arr, np.ndarray):
            raise LoadMatrixError("Bitpack matrix must be a NumPy array.")
        if not self._is_le_uint64(mmap_arr):
            raise LoadMatrixError(f"Bitpack dtype must be little-endian uint64 ('<u8'), got {mmap_arr.dtype}.")
        if mmap_arr.ndim != 2:
            raise LoadMatrixError("Bitpack matrix must be 2-dimensional.")
        if not mmap_arr.flags['C_CONTIGUOUS']:
            raise LoadMatrixError("Bitpack matrix must be C-contiguous.")

        n_rows, words = mmap_arr.shape

        if not isinstance(meta_raw, Dict):
            raise LoadMatrixError("Bitpack JSON missing 'meta' Dictionary.")
        if meta_raw.get("format") not in  supported_formats:
            raise LoadMatrixError(f"Unsupported bitpack format: {meta_raw.get('format')}")

        if meta_raw.get("dtype") != "<u8" or meta_raw.get("endianness") != "little":
            raise LoadMatrixError("dtype/endianness must be '<u8' and 'little'.")
        if not bool(meta_raw.get("row_major", True)):
            raise LoadMatrixError("Only row-major bitpack is supported.")
        if int(meta_raw.get("bits_per_word", 64)) != 64:
            raise LoadMatrixError("bits_per_word must be 64 for uint64 packing.")

        n_cols_bits = int(meta_raw["n_cols"])
        expected_words = (n_cols_bits + 63) // 64
        if words != expected_words:
            raise LoadMatrixError(f"words_per_row mismatch: header {expected_words}, file {words}")

        ## file size check if loaded from an uncompressed .bin. v2 files can
        ## contain appended binary sections after the allele matrix.
        bin_resolved = meta_raw.get("data_file_resolved")
        if bin_resolved:
            bin_path = Path(bin_resolved)
            expected_bytes = n_rows * expected_words * 8
            size = bin_path.stat().st_size
            if meta_raw.get("format") in ("ardal.bin.v2", "ardal.bin.v3"):
                self._validate_binary_sections(meta_raw, size, expected_bytes)
            elif size != expected_bytes:
                raise LoadMatrixError(f"Binary size mismatch: expected {expected_bytes}, got {size}")
        else:
            expected_bytes = n_rows * expected_words * 8

        ## headers check
        if not isinstance(headers_raw, Dict):
            raise LoadMatrixError("Headers must be a Dictionary.")
        if "guids" not in headers_raw or "alleles" not in headers_raw:
            raise LoadMatrixError("Headers must contain 'guids' and 'alleles' keys.")

        guids = headers_raw["guids"]
        alleles = headers_raw["alleles"]
        if len(guids) != n_rows:
            raise LoadMatrixError(f"Mismatch: {n_rows} rows vs {len(guids)} GUIDs.")
        if len(alleles) != n_cols_bits:
            raise LoadMatrixError(f"Mismatch: {n_cols_bits} allele bits vs {len(alleles)} labels.")
        if len(set(guids)) != len(guids):
            raise LoadMatrixError("GUIDs must be unique.")
        if len(set(alleles)) != len(alleles):
            raise LoadMatrixError("Alleles must be unique.")

        ## optional integrity check
        want_hash = self.verify_hash and meta_raw.get("data_sha256")
        if want_hash:
            hash_target = meta_raw.get("data_file_resolved") or meta_raw.get("compressed_data_file_resolved")
            if not hash_target:
                raise LoadMatrixError("Hash verification requested but no matrix file path was recorded.")
            digest = self._sha256_file(Path(hash_target))
            if digest != meta_raw["data_sha256"]:
                raise LoadMatrixError(
                    f"SHA256 mismatch for {hash_target}: expected {self.meta['data_sha256']} got {digest}"
                )
            
    
    def _validate_dense_pair( self,
                              dense: np.ndarray,
                              headers: Dict ) -> None:
        """
        Validate a dense matrix + headers without mutating them.

        Requirements:
        - dense: 2D, bool or integer with values in {0,1}
        - headers: has 'guids' and 'alleles' (lists of str), lengths match, unique
        """
        ## --- matrix checks ---
        if not isinstance(dense, np.ndarray):
            raise LoadMatrixError("Dense matrix must be a NumPy array.")
        if dense.ndim != 2:
            raise LoadMatrixError("Dense matrix must be 2-dimensional.")

        n_rows, n_cols_bits = dense.shape
        kind = dense.dtype.kind

        if dense.dtype == np.bool_:
            pass  ## {0,1} by construction
        elif kind in ("i", "u"):
            ## allocation-free bounds check
            minv = int(dense.min())
            maxv = int(dense.max())
            if minv < 0 or maxv > 1:
                raise LoadMatrixError(f"Dense matrix contains values outside {{0,1}} "
                                    f"(min={minv}, max={maxv}).")
        else:
            raise LoadMatrixError(f"Dense dtype {dense.dtype} not supported; use bool or integer {{0,1}}.")

        ## --- headers checks ---
        if not isinstance(headers, Dict):
            raise LoadMatrixError("Headers must be a Dictionary.")
        if "guids" not in headers or "alleles" not in headers:
            raise LoadMatrixError("Headers must contain 'guids' and 'alleles' keys.")

        guids = headers["guids"]
        alleles = headers["alleles"]

        if not isinstance(guids, list) or not all(isinstance(g, str) for g in guids):
            raise LoadMatrixError("GUIDs must be a list of strings.")
        if not isinstance(alleles, list) or not all(isinstance(a, str) for a in alleles):
            raise LoadMatrixError("Alleles must be a list of strings.")

        if len(guids) != n_rows:
            raise LoadMatrixError(f"Mismatch: {n_rows} matrix rows vs {len(guids)} GUIDs.")
        if len(alleles) != n_cols_bits:
            raise LoadMatrixError(f"Mismatch: {n_cols_bits} matrix columns vs {len(alleles)} alleles.")

        if len(set(guids)) != len(guids):
            raise LoadMatrixError("GUIDs must be unique.")
        if len(set(alleles)) != len(alleles):
            raise LoadMatrixError("Alleles must be unique.")


    ## ------------- loaders -------------
    def _load_bitpack_pair( self,
                            header_meta_path: Path,
                            bin_path: Path ) -> Tuple[np.ndarray, Dict, Dict, Union[Dict, None]]:
        """ Load bitpacked matrix/headers meta pair
        """
        ## parse and load header metadata
        headers_raw, meta_raw, missing_raw = self._load_header_meta(header_meta_path)            
        
        ## check metadata and the input binary are congruous
        hdr_bin = meta_raw.get("matrix_file") or meta_raw.get("data_file")
        matrix_format = self._matrix_path_format(bin_path)
        hdr_name = Path(hdr_bin).name if hdr_bin else None
        compressed_match = (
            matrix_format == ".bin.zst"
            and hdr_name == bin_path.with_suffix("").name
        )
        if hdr_name and hdr_name != bin_path.name and not compressed_match:
            log.warning(f"Header data_file='{hdr_bin}' != provided bin '{bin_path.name}'. Using provided bin.")

        n_rows = int(meta_raw["n_rows"])
        words  = int(meta_raw["words_per_row"])
        allele_section = self._get_allele_matrix_section(meta_raw, n_rows, words)

        if matrix_format == ".bin":
            meta_raw["data_file_resolved"] = str(bin_path.resolve())
            compression = allele_section.get("compression")
            if compression in (None, "", "none"):
                ## memmap the array
                try:
                    mmap_arr = np.memmap(
                        bin_path,
                        mode="r",
                        dtype=np.dtype(allele_section["dtype"]),
                        offset=int(allele_section["offset"]),
                        shape=tuple(allele_section["shape"]),
                        order="C",
                    )
                except Exception as e:
                    raise LoadMatrixError(f"Failed to load packed matrix: {e}")
            elif compression == "zstd":
                mmap_arr = self._load_compressed_allele_section(bin_path, allele_section)
            else:
                raise UnsupportedFormatError(f"Unsupported allele_matrix compression: {compression}")
            missing_raw = self._load_binary_missing_masks(meta_raw, headers_raw, bin_path, missing_raw)
            return mmap_arr, headers_raw, meta_raw, missing_raw

        if matrix_format == ".bin.zst":
            if meta_raw.get("format") in ("ardal.bin.v2", "ardal.bin.v3"):
                raise UnsupportedFormatError("Sectioned Ardal binary formats are not supported for whole-file .bin.zst inputs.")
            meta_raw["compressed_data_file_resolved"] = str(bin_path.resolve())
            meta_raw["data_file_resolved"] = None
            return self._load_bitpack_zstd(bin_path, n_rows, words), headers_raw, meta_raw, missing_raw
        
        raise UnsupportedFormatError(f"Unsupported bitpacked matrix format: {bin_path.name}")
    
    
    def _load_npy_pair( self,
                        header_meta_path: Path,
                        npy_path: Path ) -> Tuple[np.ndarray, Dict, Dict, Union[Dict, None]]:
        """ Load npy matrix/headers meta pair
        """
        ## parse and load header metadata
        headers_raw, meta_raw, missing_raw = self._load_header_meta(header_meta_path) 
        
        ## check metadata and the input binary are congruous
        hdr_bin = meta_raw.get("matrix_file") or meta_raw.get("data_file")
        if hdr_bin and Path(hdr_bin).name != npy_path.name:
            log.warning(f"Header data_file='{hdr_bin}' != provided bin '{npy_path.name}'. Using provided bin.")
            
        meta_raw["data_file_resolved"] = str(npy_path.resolve())
        
        try:
            matrix = np.ascontiguousarray(np.load(npy_path))
        except Exception as e:
            raise LoadMatrixError(f"Failed to load npy matrix: {e}")
                
        return matrix, headers_raw, meta_raw, missing_raw


    def _load_npz_pair( self,
                        header_meta_path: Path,
                        npz_path: Path ) -> Tuple[np.ndarray, Dict, Dict, Union[Dict, None]]:
        sp_sparse = require_package("scipy", attr="sparse")
        
        ## parse and load header metadata
        headers_raw, meta_raw, missing_raw = self._load_header_meta(header_meta_path) 
        
        ## check metadata and the input binary are congruous
        hdr_bin = meta_raw.get("matrix_file") or meta_raw.get("data_file")
        if hdr_bin and Path(hdr_bin).name != npz_path.name:
            log.warning(f"Header data_file='{hdr_bin}' != provided bin '{npz_path.name}'. Using provided bin.")
            
        meta_raw["data_file_resolved"] = str(npz_path.resolve())
        
        try:
            sp_mat = sp_sparse.load_npz(npz_path)
            matrix = np.ascontiguousarray(sp_mat.toarray())
        except Exception as sp_error:
            try:
                data = np.load(npz_path)
                if 'matrix' not in data:
                    raise ValueError("Key 'matrix' not found in .npz file.")
                matrix = np.ascontiguousarray(data['matrix'])
            except Exception as np_error:
                raise LoadMatrixError(
                    f"Failed to load matrix from npz: "
                    f"scipy.sparse.load_npz failed with {sp_error}; "
                    f"np.load failed with {np_error}"
                )
                
        return matrix, headers_raw, meta_raw, missing_raw
    
    
    def _load_csv( self,
                   csv_path: str ) -> Tuple[np.ndarray, Dict]:
        """
        Load a dense CSV where first column is GUID and remaining columns are {0,1}.
        Returns (matrix: np.ndarray[C-contig, uint8], headers: Dict).
        """
        path = Path(csv_path)
        if not path.exists():
            raise LoadMatrixError(f"CSV file does not exist: {path}")

        with path.open("r", newline="") as f:
            reader = csv.reader(f)
            try:
                header = next(reader)
            except StopIteration:
                raise LoadMatrixError("CSV is empty.")
        if len(header) < 2:
            raise LoadMatrixError("CSV must have an index column + at least one allele column.")

        alleles = [h.strip() for h in header[1:]]
        n_cols_bits = len(alleles)
        if len(set(alleles)) != n_cols_bits:
            raise LoadMatrixError("Allele headers must be unique.")

        guids: list[str] = []
        seen_guids: set[str] = set()
        rows: list[np.ndarray] = []

        with path.open("r", newline="") as f:
            reader = csv.reader(f)
            _ = next(reader, None)  # skip header
            r_global = 0
            for row in reader:
                if not row:
                    continue
                if len(row) != n_cols_bits + 1:
                    raise LoadMatrixError(f"Row {r_global} has {len(row)-1} allele columns; expected {n_cols_bits}.")

                guid = row[0].strip()
                if guid in seen_guids:
                    raise LoadMatrixError(f"Duplicate GUID encountered: {guid}")
                seen_guids.add(guid)
                guids.append(guid)

                vals = row[1:]
                arr = np.zeros(n_cols_bits, dtype=np.uint8)
                for j, tok in enumerate(vals):
                    t = tok.strip()
                    if t == "1" or t == "1.0":
                        arr[j] = 1
                    elif t == "" or t == "0" or t == "0.0":
                        continue
                    else:
                        raise LoadMatrixError(f"Non-binary token at row {r_global}, col {j+1}: '{tok}'")
                rows.append(arr)
                r_global += 1

        if rows:
            matrix = np.vstack(rows).astype(np.uint8, copy=False)
        else:
            matrix = np.zeros((0, n_cols_bits), dtype=np.uint8)

        matrix = np.ascontiguousarray(matrix)
        headers = {"guids": guids, "alleles": alleles}
        return matrix, headers


    def _get_allele_matrix_section( self,
                                    meta_raw: Dict,
                                    n_rows: int,
                                    words: int ) -> Dict:
        sections = meta_raw.get("sections")
        if not isinstance(sections, Dict) or "allele_matrix" not in sections:
            return {
                "offset": 0,
                "dtype": "<u8",
                "shape": [n_rows, words],
                "nbytes": n_rows * words * 8,
            }

        section = sections["allele_matrix"]
        if not isinstance(section, Dict):
            raise LoadMatrixError("ardal.bin.v2 metadata has malformed allele_matrix section.")

        shape = section.get("shape")
        if not isinstance(shape, list) or len(shape) != 2:
            raise LoadMatrixError("ardal.bin.v2 allele_matrix section must define shape [n_rows, words_per_row].")
        if int(shape[0]) != n_rows or int(shape[1]) != words:
            raise LoadMatrixError(
                f"ardal.bin.v2 allele_matrix shape mismatch: section {shape}, meta {[n_rows, words]}"
            )

        dtype = np.dtype(section.get("dtype", "<u8"))
        if dtype != np.dtype("<u8"):
            raise LoadMatrixError(f"ardal.bin.v2 allele_matrix dtype must be '<u8', got {dtype}.")

        return {
            "offset": int(section.get("offset", 0)),
            "dtype": dtype,
            "shape": [n_rows, words],
            "nbytes": int(section.get("nbytes", n_rows * words * dtype.itemsize)),
            "compression": section.get("compression"),
            "uncompressed_nbytes": int(section.get("uncompressed_nbytes", n_rows * words * dtype.itemsize)),
        }


    def _validate_binary_sections( self,
                                   meta_raw: Dict,
                                   file_size: int,
                                   expected_matrix_bytes: int ) -> None:
        sections = meta_raw.get("sections")
        if not isinstance(sections, Dict):
            raise LoadMatrixError("ardal.bin.v2 metadata missing 'sections'.")

        required = ("allele_matrix",)
        for name in required:
            if name not in sections:
                raise LoadMatrixError(f"ardal.bin.v2 metadata missing section '{name}'.")

        for name, section in sections.items():
            if not isinstance(section, Dict):
                raise LoadMatrixError(f"ardal.bin.v2 section '{name}' is malformed.")
            try:
                offset = int(section.get("offset", 0))
                nbytes = int(section.get("nbytes", 0))
            except (TypeError, ValueError) as exc:
                raise LoadMatrixError(f"ardal.bin.v2 section '{name}' has invalid offset/nbytes.") from exc
            if offset < 0 or nbytes < 0 or offset + nbytes > file_size:
                raise LoadMatrixError(
                    f"ardal.bin.v2 section '{name}' exceeds file bounds: offset={offset}, nbytes={nbytes}, file={file_size}."
                )

        allele_section = sections["allele_matrix"]
        compression = allele_section.get("compression")
        if compression in (None, "", "none"):
            allele_nbytes = int(allele_section.get("nbytes", expected_matrix_bytes))
            if allele_nbytes != expected_matrix_bytes:
                raise LoadMatrixError(
                    f"ardal.bin.v2 allele_matrix byte size mismatch: expected {expected_matrix_bytes}, got {allele_nbytes}."
                )
        elif compression == "zstd":
            uncompressed_nbytes = int(allele_section.get("uncompressed_nbytes", 0))
            if uncompressed_nbytes != expected_matrix_bytes:
                raise LoadMatrixError(
                    f"Compressed allele_matrix uncompressed byte size mismatch: expected {expected_matrix_bytes}, got {uncompressed_nbytes}."
                )
        else:
            raise UnsupportedFormatError(f"Unsupported allele_matrix compression: {compression}")


    def _load_binary_missing_masks( self,
                                    meta_raw: Dict,
                                    headers_raw: Dict,
                                    bin_path: Path,
                                    missing_raw: Union[Dict, None] ) -> Union[Dict, None]:
        if not isinstance(missing_raw, Dict):
            return missing_raw

        descriptor = missing_raw.get("missing") if "missing" in missing_raw else missing_raw
        if not isinstance(descriptor, Dict):
            return missing_raw
        if descriptor.get("encoding") != "binary_sections":
            return missing_raw

        sections = meta_raw.get("sections")
        if not isinstance(sections, Dict):
            raise LoadMatrixError("Binary missing masks require metadata sections.")

        offsets_name = descriptor.get("offsets_section", "missing_offsets")
        ranges_name = descriptor.get("ranges_section", "missing_ranges")
        offsets_section = sections.get(offsets_name)
        ranges_section = sections.get(ranges_name)
        
        if not isinstance(offsets_section, Dict) or not isinstance(ranges_section, Dict):
            raise LoadMatrixError("Binary missing masks require missing_offsets and missing_ranges sections.")

        guids = headers_raw.get("guids", []) if isinstance(headers_raw, Dict) else []
        n_rows = len(guids)
        
        offsets = np.memmap(
            bin_path,
            mode="r",
            dtype=np.dtype(offsets_section.get("dtype", "<u8")),
            offset=int(offsets_section["offset"]),
            shape=(int(offsets_section.get("length", n_rows + 1)),),
        )
        
        if offsets.shape[0] != n_rows + 1:
            raise LoadMatrixError(
                f"Binary missing offsets length mismatch: expected {n_rows + 1}, got {offsets.shape[0]}."
            )

        ranges_shape = ranges_section.get("shape", [0, 2])
        if not isinstance(ranges_shape, list) or len(ranges_shape) != 2 or int(ranges_shape[1]) != 2:
            raise LoadMatrixError("Binary missing ranges section must define shape [n_ranges, 2].")
        
        ranges = np.memmap(
            bin_path,
            mode="r",
            dtype=np.dtype(ranges_section.get("dtype", "<u4")),
            offset=int(ranges_section["offset"]),
            shape=(int(ranges_shape[0]), 2),
        )
        log.debug(
            "Binary missing interval sections mapped for lazy inflation: "
            f"offsets_len={offsets.shape[0]}, ranges_shape={ranges.shape}, "
            f"ranges_dtype={ranges.dtype}."
        )

        for row_idx, guid in enumerate(guids):
            start = int(offsets[row_idx])
            end = int(offsets[row_idx + 1])
            
            if start < 0 or end < start or end > ranges.shape[0]:
                raise LoadMatrixError(f"Binary missing range offsets malformed for GUID '{guid}'.")
            
            for range_start, range_end in ranges[start:end]:
                left = int(range_start)
                right = int(range_end)
                
                if right < left:
                    raise LoadMatrixError(f"Binary missing range malformed for GUID '{guid}': [{left}, {right}).")
        log.debug(
            "Prepared lazy binary missing interval store during matrix load: "
            f"rows={n_rows}, ranges={ranges.shape[0]}."
        )

        return BinaryMissingMaskStore(guids, offsets, ranges)
        

    def _load_parquet( self,
                       filepath: str ) -> Tuple[np.ndarray, Dict]:
        df = pd.read_parquet(filepath, engine="fastparquet")
        matrix = df.values
        headers = {"guids": list(df.index), "alleles": list(df.columns)}
        return matrix, headers


    def _load_header_meta( self,
                           header_meta_path: Path ) -> Tuple[Dict, Dict, Union[Dict, None]]:
        """ Load metadata headers JSON
        """
        with self._open_text_metadata(header_meta_path) as f:
            header_meta_data = json.load(f)
        
        headers_raw, missing_raw, meta_raw = self._load_headers_dict(header_meta_data)
        if meta_raw is None:
            meta_raw = {}
        
        ## handle missing sites data
        ## this supports legacy formats where the headers-meta JSON does not contain missing sites data
        if missing_raw is None:
            log.debug("Legacy headers detected. No missing sites provided.")
        
        return headers_raw, meta_raw, missing_raw


    ## ----- helpers -----
    def _is_bitpack_header( self,
                            meta_path: Path ) -> bool:
        try:
            with self._open_text_metadata(meta_path) as f:
                obj = json.load(f)
            return isinstance(obj, Dict) and "meta" in obj and obj["meta"].get("format") == "ardal.bitpack.v1"
        except Exception:
            return False


    @staticmethod
    def _is_meta_path(path: Path) -> bool:
        suffixes = [suffix.lower() for suffix in path.suffixes]
        if suffixes[-2:] == [".meta", ".gz"]:
            return True
        return bool(suffixes and suffixes[-1] in {".meta", ".json"})


    @staticmethod
    def _open_text_metadata(path: Path):
        suffixes = [suffix.lower() for suffix in path.suffixes]
        if suffixes[-2:] == [".meta", ".gz"]:
            return gzip.open(path, "rt")
        return open(path, "r")
        
        
    def _normalise_missing_masks( self,
                                  missing_obj: Union[Dict, None],
                                  headers: Dict ) -> Dict[str, Dict]:
        guids: list[str] = []
        
        if isinstance(headers, Dict):
            raw_guids = headers.get("guids", [])
            if isinstance(raw_guids, list):
                guids = [g for g in raw_guids if isinstance(g, str)]
        else:
            raise LoadMatrixError(f"headers should be <dict>, not <{type(headers)}>.")

        if isinstance(missing_obj, BinaryMissingMaskStore):
            return missing_obj

        per_guid = {guid: [] for guid in guids}

        if not missing_obj or not isinstance(missing_obj, Dict):
            return {"column_masks": per_guid}

        if isinstance(missing_obj.get("missing"), Dict):
            missing_obj = missing_obj["missing"]

        if isinstance(missing_obj.get("column_masks"), Dict):
            missing_obj = missing_obj["column_masks"]
        elif isinstance(missing_obj.get("missing_masks"), Dict):
            missing_obj = missing_obj["missing_masks"]

        if isinstance(missing_obj, Dict):
            iterable = missing_obj.items()
            for guid, sites in iterable:
                if guid not in per_guid:
                    log.warning(f"Missing-sites data references unknown GUID '{guid}'. Skipping.")
                    continue
                per_guid[guid] = self._coerce_site_list(sites)
        else:
            raise LoadMatrixError(f"missing sites data malformed.") 

        return {"column_masks": per_guid}
    
    
    def _ensure_meta( self,
                      matrix: np.ndarray,
                      headers: Dict,
                      meta_obj: Union[Dict, None],
                      *,
                      is_bitpacked: bool ) -> Dict:
        """
        Ensure metadata is available for downstream consumers.
        If a meta dictionary already exists, return it unchanged.
        Otherwise synthesise a placeholder tailored to the matrix type.
        """
        if isinstance(meta_obj, Dict):
            return meta_obj
        else:
            fmt = "ardal.bitpack.v1" if is_bitpacked else "ardal.dense.v1"
            return make_meta(matrix,
                             headers,
                             generated_by="ardal::ArdalParser",
                             format_name=fmt,
                             matrix_file=None)
    
    
    @staticmethod
    def _matrix_path_format(path: Path) -> str:
        suffixes = [suffix.lower() for suffix in path.suffixes]
        if suffixes[-2:] == [".bin", ".zst"]:
            return ".bin.zst"
        if suffixes:
            return suffixes[-1]
        return ""


    def _load_bitpack_zstd(self,
                           bin_path: Path,
                           n_rows: int,
                           words: int) -> np.ndarray:
        zstandard = require_package(
            "zstandard",
            error_message="The package 'zstandard' is required to load .bin.zst matrices. Install it with `pip install zstandard`."
        )
        expected_bytes = n_rows * words * 8
        dctx = zstandard.ZstdDecompressor()
        data = bytearray()

        try:
            with open(bin_path, "rb") as fin, dctx.stream_reader(fin) as reader:
                while True:
                    chunk = reader.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    data.extend(chunk)
        except Exception as e:
            raise LoadMatrixError(f"Failed to decompress packed matrix '{bin_path.name}': {e}")

        if len(data) != expected_bytes:
            raise LoadMatrixError(
                f"Decompressed binary size mismatch for {bin_path.name}: expected {expected_bytes}, got {len(data)}"
            )

        arr = np.frombuffer(data, dtype=np.dtype("<u8")).reshape((n_rows, words))
        return np.ascontiguousarray(arr)


    def _load_compressed_allele_section( self,
                                         bin_path: Path,
                                         section: Dict ) -> np.memmap:
        """Decompress a compressed allele_matrix section to a temporary memmap file."""
        compression = section.get("compression")
        if compression != "zstd":
            raise UnsupportedFormatError(f"Unsupported allele_matrix compression: {compression}")

        zstandard = require_package(
            "zstandard",
            error_message="The package 'zstandard' is required to load zstd-compressed allele matrix sections. Install it with `pip install zstandard`."
        )
        offset = int(section["offset"])
        nbytes = int(section["nbytes"])
        expected_bytes = int(section["uncompressed_nbytes"])
        dtype = np.dtype(section.get("dtype", "<u8"))
        shape = tuple(section["shape"])

        try:
            with open(bin_path, "rb") as fin:
                fin.seek(offset)
                compressed = fin.read(nbytes)
        except Exception as exc:
            raise LoadMatrixError(f"Failed to read compressed allele_matrix section: {exc}") from exc

        try:
            dctx = zstandard.ZstdDecompressor()
            tmp_path = None
            with tempfile.NamedTemporaryFile(prefix="ardal_allele_matrix_", suffix=".bin", delete=False) as tmp:
                tmp_path = tmp.name
                with dctx.stream_reader(io.BytesIO(compressed)) as reader:
                    total = 0
                    while True:
                        chunk = reader.read(8 * 1024 * 1024)
                        if not chunk:
                            break
                        total += len(chunk)
                        tmp.write(chunk)
                if total != expected_bytes:
                    raise LoadMatrixError(
                        f"Decompressed allele_matrix size mismatch: expected {expected_bytes}, got {total}."
                    )
        except LoadMatrixError:
            try:
                Path(tmp_path).unlink(missing_ok=True)
            except Exception:
                pass
            raise
        except Exception as exc:
            try:
                Path(tmp_path).unlink(missing_ok=True)
            except Exception:
                pass
            raise LoadMatrixError(f"Failed to decompress allele_matrix section: {exc}") from exc

        self._temporary_matrix_files.append(tmp_path)
        try:
            return np.memmap(
                tmp_path,
                mode="r",
                dtype=dtype,
                offset=0,
                shape=shape,
                order="C",
            )
        except Exception as exc:
            raise LoadMatrixError(f"Failed to memmap decompressed allele_matrix section: {exc}") from exc


    @staticmethod
    def _is_le_uint64( arr: np.ndarray ) -> bool:
        """True if dtype is little-endian uint64 (<u8) or native-little uint64."""
        dt = arr.dtype
        if dt == np.dtype("<u8"):
            return True
        if dt == np.uint64 and byteorder == "little":
            return True
        return False


    @staticmethod
    def _is_bitpacked_candidate( arr: np.ndarray ) -> bool:
        return (
            isinstance(arr, np.ndarray)
            and arr.ndim == 2
            and arr.dtype.kind == "u"
            and arr.dtype.itemsize == 8
        )


    @staticmethod
    def _generated_by_string() -> str:
        try:
            return "Ardal v" + version("ardal")
        except PackageNotFoundError:
            return "Ardal"
    
    
    @staticmethod
    def _order_input( a : Union[Dict,np.ndarray],
                      b : Union[Dict,np.ndarray] ) -> Tuple[np.ndarray, Dict]:
        """ Returns the data tuple in a preDictable order
        """
        if isinstance(a, np.ndarray) and isinstance(b, Dict):
            return [a, b]
        elif isinstance(a, Dict) and isinstance(b, np.ndarray):
            return [b, a]
        else:
            raise MalformedInputError("Input list from memory must contain two elements: matrix (np.ndarray) and headers (Dict).")


    @staticmethod
    def _coerce_site_list( sites : Union[None, list, tuple, set] ) -> list:
        if sites is None:
            return []
        if isinstance(sites, list):
            return sites
        if isinstance(sites, (tuple, set)):
            return list(sites)
        return [sites]


    @staticmethod
    def _sha256_file( path: Path,
                      chunk_mb: int = 8 ) -> str:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_mb * 1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()
