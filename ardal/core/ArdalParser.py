# core/ArdalParser.py
import os
import csv
import json
import hashlib
from pathlib import Path
from typing import Union, Tuple

import pandas as pd
import numpy as np

from ..utils.misc import require_package
from ..utils.exceptions import MalformedInputError, UnsupportedFormatError, LoadMatrixError
from ..utils.logger import get_logger

try:
    from _ardal_pack import pack_dense_to_words as _pack_words_cpp
except Exception:
    _pack_words_cpp = None
    
    
log = get_logger()


class ArdalParser:
    """
    Parses Ardal inputs:
      - Dense: CSV, Parquet, NPY/JSON, NPZ/JSON
      - Bitpacked: {headers.json, matrix.bin} with meta.format == "ardal.bitpack.v1"

    Returns:
      self.matrix  -> np.ndarray
        - dense: shape (n_rows, n_cols), dtype arbitrary
        - bitpack: memmap '<u8', shape (n_rows, words_per_row)
      self.headers -> dict with "guids", "alleles"
      self.meta    -> dict (bitpack only), else {}
      self.is_bitpacked -> bool
    """

    def __init__(self,
                 input_data_structure: Union[list, str],
                 verify_hash: bool = False):
        self.input_data = input_data_structure
        self.file_format = None
        self.matrix: np.ndarray | None = None
        self.headers: dict = {}
        self.meta: dict = {}
        self.is_bitpacked: bool = False
        self.verify_hash = verify_hash
        self._parse()



    ## ------------- parsing -------------

    def _parse(self) -> Union[int, None]:
        if self.input_data is None:
            raise MalformedInputError("Input data structure cannot be None.")

        ## in-memory [array, headers] or [headers, array]
        if isinstance(self.input_data, list):
            
            if len(self.input_data) != 2:
                raise MalformedInputError("Input list must contain two elements: matrix and headers (or two file paths).")

            a, b = self.input_data

            if isinstance(a, np.ndarray) and isinstance(b, dict):
                log.debug(f"Parsing matrix data from list in memory of form: {[type(i) for i in self.input_data]}")
                self.matrix, self.headers = np.ascontiguousarray(a), b
                self.file_format = "dense"
                self._validate_dense()
                return 0
            if isinstance(b, np.ndarray) and isinstance(a, dict):
                log.debug(f"Parsing matrix data from list in memory of form: {[type(i) for i in self.input_data]}")
                self.matrix, self.headers = np.ascontiguousarray(b), a
                self.file_format = "dense"
                self._validate_dense()
                return 0

            ## two file paths
            if isinstance(a, str) and isinstance(b, str):
                log.debug(f"Parsing matrix data from two file paths: {self.input_data}")
                
                a, b = Path(a), Path(b)
                if not a.exists() or not b.exists():
                    raise FileNotFoundError(f"One or more file paths do not exist: {a}, {b}")
                exts = {a.suffix.lower(), b.suffix.lower()}
                
                json_path = a if a.suffix.lower() == ".json" else b
                mat_path  = b if json_path == a else a

                ## bitpack pair
                if exts == {".json", ".bin"}:
                    log.debug(f"Detected matrix format: .bit")
                    self.matrix, self.headers, self.meta = self._load_bitpack_pair(json_path, mat_path)
                    self.file_format = "bitpack"
                    self.is_bitpacked = True
                    self._validate_bitpack()
                    return 0

                ## dense pairs
                if exts == {".json", ".npy"}:
                    log.debug(f"Detected matrix format: .npy")
                    dense, headers = self._load_npy_pair(str(json_path), str(mat_path))
                    self._validate_dense_pair(dense, headers)
                    self.matrix = self._pack_dense_to_words(dense)
                    self.headers = headers
                    self.meta = self._mk_bitpack_meta_from_dense(self.matrix, headers)
                    self.file_format = "bitpack"
                    self.is_bitpacked = True
                    self._validate_bitpack()
                    return 0

                if exts == {".json", ".npz"}:
                    log.debug(f"Detected matrix format: .npz")
                    dense, headers = self._load_npz_pair(str(json_path), str(mat_path))
                    self._validate_dense_pair(dense, headers)
                    self.matrix = self._pack_dense_to_words(dense)
                    self.headers = headers
                    self.meta = self._mk_bitpack_meta_from_dense(self.matrix, headers)
                    self.file_format = "bitpack"
                    self.is_bitpacked = True
                    self._validate_bitpack()
                    return 0

                raise UnsupportedFormatError(f"Unrecognized file pair: {a.suffix}, {b.suffix}")

            raise MalformedInputError(
                "If list input, it must be [headers::dict, matrix::np.ndarray] or [headers.json, matrix.bin] or dense pairs."
            )

        ## single path
        if isinstance(self.input_data, str):
            log.debug(f"Parsing matrix data from one file path: {self.input_data}")
            
            p = Path(self.input_data)
            if not p.exists():
                raise FileNotFoundError(f"File does not exist: {p}")

            if p.suffix.lower() == ".json":
                ## might be a bitpack header
                if self._is_bitpack_header(p):
                    self.matrix, self.headers, self.meta = self._load_bitpack_header(p)
                    self.file_format = "bitpack"
                    self.is_bitpacked = True
                    self._validate_bitpack()
                    return 0
                ## else it's a dense header without its matrix -> reject
                raise MalformedInputError("Dense JSON headers provided without matrix path.")

            if p.suffix.lower() == ".csv":
                log.debug(f"Detected matrix format: .csv")
                words, headers, meta = self._load_csv_bitpacked(str(p))
                self.matrix, self.headers, self.meta = words, headers, meta
                self.file_format = "bitpack"
                self.is_bitpacked = True
                self._validate_bitpack()
                return 0

            if p.suffix.lower() == ".parquet":
                log.debug(f"Detected matrix format: .parquet")
                dense, headers = self._load_parquet(str(p))
                self._validate_dense_pair(dense, headers)
                self.matrix = self._pack_dense_to_words(dense)
                self.headers = headers
                self.meta = self._mk_bitpack_meta_from_dense(self.matrix, headers)
                self.file_format = "bitpack"
                self.is_bitpacked = True
                self._validate_bitpack()
                return 0

            if p.suffix.lower() in (".npy", ".npz", ".bin"):
                raise MalformedInputError("Binary file provided without matching headers.json; provide both.")

            raise UnsupportedFormatError(f"Unsupported file extension: {p.suffix}")

        raise MalformedInputError("Input must be list or string.")



    ## ------------- validation -------------

    def _validate_dense(self) -> None:
        if not isinstance(self.matrix, np.ndarray):
            raise LoadMatrixError("Matrix must be a NumPy array.")
        if self.matrix.ndim != 2:
            raise LoadMatrixError("Matrix must be 2-dimensional.")
        if not self.matrix.flags['C_CONTIGUOUS']:
            try:
                self.matrix = np.ascontiguousarray(self.matrix)
            except Exception as e:
                raise LoadMatrixError(f"Failed to make array contiguous in memory : {e}")

        n_rows, n_cols = self.matrix.shape

        if not isinstance(self.headers, dict):
            raise LoadMatrixError("Headers must be a dictionary.")
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
        if len(alleles) != n_cols:
            raise LoadMatrixError(f"Mismatch: {n_cols} matrix columns vs {len(alleles)} alleles.")

        if len(set(guids)) != len(guids):
            raise LoadMatrixError("GUIDs must be unique.")
        if len(set(alleles)) != len(alleles):
            raise LoadMatrixError("Alleles must be unique.")


    def _validate_bitpack(self) -> None:
        if not isinstance(self.matrix, np.ndarray):
            raise LoadMatrixError("Bitpack matrix must be a NumPy memmap.")
        if self.matrix.dtype != np.dtype("<u8"):
            raise LoadMatrixError(f"Bitpack dtype must be little-endian uint64 ('<u8'), got {self.matrix.dtype}.")
        if self.matrix.ndim != 2:
            raise LoadMatrixError("Bitpack matrix must be 2-dimensional.")
        if not self.matrix.flags['C_CONTIGUOUS']:
            raise LoadMatrixError("Bitpack memmap must be C-contiguous.")

        n_rows, words = self.matrix.shape

        if not isinstance(self.meta, dict):
            raise LoadMatrixError("Bitpack JSON missing 'meta' dictionary.")
        if self.meta.get("format") != "ardal.bitpack.v1":
            raise LoadMatrixError(f"Unsupported bitpack format: {self.meta.get('format')}")

        if self.meta.get("dtype") != "<u8" or self.meta.get("endianness") != "little":
            raise LoadMatrixError("dtype/endianness must be '<u8' and 'little'.")
        if not bool(self.meta.get("row_major", True)):
            raise LoadMatrixError("Only row-major bitpack is supported.")
        if int(self.meta.get("bits_per_word", 64)) != 64:
            raise LoadMatrixError("bits_per_word must be 64 for uint64 packing.")

        n_cols = int(self.meta["n_cols"])
        expected_words = (n_cols + 63) // 64
        if words != expected_words:
            raise LoadMatrixError(f"words_per_row mismatch: header {expected_words}, file {words}")

        ## file size check if loaded from a .bin
        bin_resolved = self.meta.get("data_file_resolved")
        if bin_resolved:
            bin_path = Path(bin_resolved)
            expected_bytes = n_rows * expected_words * 8
            size = bin_path.stat().st_size
            if size != expected_bytes:
                raise LoadMatrixError(f"Binary size mismatch: expected {expected_bytes}, got {size}")


        ## headers check
        if not isinstance(self.headers, dict):
            raise LoadMatrixError("Headers must be a dictionary.")
        if "guids" not in self.headers or "alleles" not in self.headers:
            raise LoadMatrixError("Headers must contain 'guids' and 'alleles' keys.")

        guids = self.headers["guids"]
        alleles = self.headers["alleles"]
        if len(guids) != n_rows:
            raise LoadMatrixError(f"Mismatch: {n_rows} rows vs {len(guids)} GUIDs.")
        if len(alleles) != n_cols:
            raise LoadMatrixError(f"Mismatch: {n_cols} allele bits vs {len(alleles)} labels.")
        if len(set(guids)) != len(guids):
            raise LoadMatrixError("GUIDs must be unique.")
        if len(set(alleles)) != len(alleles):
            raise LoadMatrixError("Alleles must be unique.")

        ## optional integrity check
        want_hash = self.verify_hash and self.meta.get("data_sha256")
        if want_hash:
            digest = self._sha256_file(bin_path)
            if digest != self.meta["data_sha256"]:
                raise LoadMatrixError(
                    f"SHA256 mismatch for {bin_path}: expected {self.meta['data_sha256']} got {digest}"
                )

    
    def _validate_dense_pair(self, dense: np.ndarray, headers: dict) -> None:
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

        n_rows, n_cols = dense.shape
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
        if not isinstance(headers, dict):
            raise LoadMatrixError("Headers must be a dictionary.")
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
        if len(alleles) != n_cols:
            raise LoadMatrixError(f"Mismatch: {n_cols} matrix columns vs {len(alleles)} alleles.")

        if len(set(guids)) != len(guids):
            raise LoadMatrixError("GUIDs must be unique.")
        if len(set(alleles)) != len(alleles):
            raise LoadMatrixError("Alleles must be unique.")



    ## ------------- loaders -------------

    def _load_csv_bitpacked(self, csv_path: str, max_chunk_mb: int = 64) -> Tuple[np.memmap, dict, dict]:
        """
        AUTOGENERATED HEADER.
        Stream a wide CSV: first column = GUID, remaining columns = {0,1}.
        Build headers, then process rows in chunks:
        - fill a dense uint8 buffer of shape (chunk_rows, n_cols)
        - pack with self._pack_dense_to_words (your existing function)
        - write into a <u8> memmap on disk

        Returns:
        (memmap<uint64 little>[n_rows, ceil(n_cols/64)], headers, meta)
        """
        ## pandas is horrible for memory, and Ardal handles very large matrices
        ## it is not recommended to store as a csv, but for small matrices this is acceptable
        ## for anything other than very small matrices, pandas is terrible and blooms memory enormously
        ## consequently, this all needs to be manual
        
        path = Path(csv_path)

        ## --- header ---
        with path.open("r", newline="") as f:
            reader = csv.reader(f)
            try:
                header = next(reader)
            except StopIteration:
                raise LoadMatrixError("CSV is empty.")
        if len(header) < 2:
            raise LoadMatrixError("CSV must have an index column + at least one allele column.")
        alleles = [h.strip() for h in header[1:]]
        n_cols = len(alleles)
        if len(set(alleles)) != n_cols:
            raise LoadMatrixError("Allele headers must be unique.")

        ## --- row count ---
        with path.open("r", newline="") as f:
            reader = csv.reader(f)
            _ = next(reader, None)  ## skip header
            n_rows = sum(1 for _ in reader)

        words_per_row = (n_cols + 63) // 64
        bin_path = path.with_suffix(path.suffix + ".bitpack.bin")
        mm = np.memmap(bin_path, mode="w+", dtype=np.dtype("<u8"),
                    shape=(n_rows, words_per_row), order="C")
        mm[:] = 0

        ## --- choose chunk size to bound RAM to ~max_chunk_mb ---
        target_bytes = max(1, int(max_chunk_mb) * 1024 * 1024)
        chunk_rows = max(1, target_bytes // max(1, n_cols))

        guids: list[str] = []
        seen_guids: set[str] = set()

        ## --- pass 2: stream + chunk pack ---
        with path.open("r", newline="") as f:
            reader = csv.reader(f)
            _ = next(reader, None)  ## skip header

            ## allocate once and reuse
            dense_chunk = np.zeros((chunk_rows, n_cols), dtype=np.uint8)
            r_global = 0
            r_in_chunk = 0

            for row in reader:
                if not row:
                    continue
                if len(row) != n_cols + 1:
                    raise LoadMatrixError(f"Row {r_global} has {len(row)-1} allele columns; expected {n_cols}.")

                guid = row[0].strip()
                if guid in seen_guids:
                    raise LoadMatrixError(f"Duplicate GUID encountered: {guid}")
                seen_guids.add(guid)
                guids.append(guid)

                ## fill current dense row
                ## strict: accept "0"/"1" (optionally "0.0"/"1.0"); treat "" as 0
                vals = row[1:]
                dr = dense_chunk[r_in_chunk]
                dr.fill(0)
                for j, tok in enumerate(vals):
                    t = tok.strip()
                    if t == "1" or t == "1.0":
                        dr[j] = 1
                    elif t == "" or t == "0" or t == "0.0":
                        ## will (should) be zero already
                        continue
                    else:
                        raise LoadMatrixError(f"Non-binary token at row {r_global}, col {j+1}: '{tok}'")

                r_in_chunk += 1
                r_global += 1

                ## if chunk full, pack and flush
                if r_in_chunk == chunk_rows:
                    words_chunk = self._pack_dense_to_words(dense_chunk[:r_in_chunk])
                    mm[r_global - r_in_chunk : r_global, :] = words_chunk
                    r_in_chunk = 0  ## reset

            ## flush any tail rows
            if r_in_chunk > 0:
                words_chunk = self._pack_dense_to_words(dense_chunk[:r_in_chunk])
                mm[r_global - r_in_chunk : r_global, :] = words_chunk

        mm.flush()

        headers = {"guids": guids, "alleles": alleles}
        meta = {
            "format": "ardal.bitpack.v1",
            "dtype": "<u8",
            "endianness": "little",
            "row_major": True,
            "n_rows": n_rows,
            "n_cols": n_cols,
            "words_per_row": words_per_row,
            "bits_per_word": 64,
            "row_stride_bytes": words_per_row * 8,
            "data_file": bin_path.name,
            "data_file_resolved": str(bin_path.resolve()),
            "data_nbytes": int(n_rows * words_per_row * 8),
            "data_sha256": None,
        }
        return mm, headers, meta



    def _load_parquet(self, filepath: str) -> Tuple[np.ndarray, dict]:
        df = pd.read_parquet(filepath, engine="fastparquet")
        matrix = df.values
        headers = {"guids": list(df.index), "alleles": list(df.columns)}
        return matrix, headers


    def _load_npy_pair(self, json_path: str, npy_path: str) -> Tuple[np.ndarray, dict]:
        try:
            matrix = np.ascontiguousarray(np.load(npy_path))
            with open(json_path, "r") as f:
                headers = json.load(f)
        except Exception as e:
            raise LoadMatrixError(f"Failed to load npy/json pair: {e}")
        return matrix, headers


    def _load_npz_pair(self, json_path: str, npz_path: str) -> Tuple[np.ndarray, dict]:
        sp_sparse = require_package("scipy", attr="sparse")
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
        try:
            with open(json_path, "r") as f:
                headers = json.load(f)
        except Exception as e:
            raise LoadMatrixError(f"Failed to load JSON headers: {e}")
        return matrix, headers



    ## ----- bitpack helpers -----

    def _is_bitpack_header(self, json_path: Path) -> bool:
        try:
            with open(json_path, "r") as f:
                obj = json.load(f)
            return isinstance(obj, dict) and "meta" in obj and obj["meta"].get("format") == "ardal.bitpack.v1"
        except Exception:
            return False


    def _load_bitpack_header(self, json_path: Path) -> Tuple[np.memmap, dict, dict]:
        with open(json_path, "r") as f:
            obj = json.load(f)
        if "meta" not in obj or "headers" not in obj:
            raise LoadMatrixError("Bitpack JSON must contain 'meta' and 'headers' keys.")
        meta = obj["meta"]
        headers = obj["headers"]

        bin_name = meta.get("data_file")
        if not bin_name:
            raise LoadMatrixError("Bitpack meta missing 'data_file'.")
        bin_path = (json_path.parent / bin_name).resolve()
        if not bin_path.exists():
            raise LoadMatrixError(f"Bitpack binary file not found: {bin_path}")
        meta["data_file_resolved"] = str(bin_path)

        n_rows = int(meta["n_rows"])
        words  = int(meta["words_per_row"])

        arr = np.memmap(bin_path, mode="r", dtype=np.dtype("<u8"), shape=(n_rows, words), order="C")
        return arr, headers, meta


    def _load_bitpack_pair(self, json_path: Path, bin_path: Path) -> Tuple[np.memmap, dict, dict]:
        with open(json_path, "r") as f:
            obj = json.load(f)
        if "meta" not in obj or "headers" not in obj:
            raise LoadMatrixError("Bitpack JSON must contain 'meta' and 'headers' keys.")
        meta = obj["meta"]
        headers = obj["headers"]

        hdr_bin = meta.get("data_file")
        if hdr_bin and Path(hdr_bin).name != bin_path.name:
            log.warning(f"Header data_file='{hdr_bin}' != provided bin '{bin_path.name}'. Using provided bin.")
        meta["data_file_resolved"] = str(bin_path.resolve())

        n_rows = int(meta["n_rows"])
        words  = int(meta["words_per_row"])

        arr = np.memmap(bin_path, mode="r", dtype=np.dtype("<u8"), shape=(n_rows, words), order="C")
        return arr, headers, meta


    @staticmethod
    def _sha256_file(path: Path, chunk_mb: int = 8) -> str:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_mb * 1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()
    
    
    def _pack_dense_to_words(self, dense: np.ndarray) -> np.ndarray:
        log.debug(f"Packing dense array.")
        
        ## expect dense uint8/bool
        ## safely coerces if not
        if dense.dtype == np.bool_:
            arr = dense.view(np.uint8)
        elif dense.dtype == np.uint8:
            arr = dense
        else:
            ## tolerate 0/1 int types
            if np.can_cast(dense.dtype, np.uint8, 'safe'):
                arr = dense.astype(np.uint8, copy=False)
            else:
                raise LoadMatrixError(f"dense dtype {dense.dtype} not packable; use bool/uint8")
        arr = np.ascontiguousarray(arr)
        if _pack_words_cpp is not None:
            log.debug(f"Running C++ packer.")
            return _pack_words_cpp(arr)  ## returns <u8 words>
        
        ## Python fallback if the C++ module is not found
        ## this is slower but safe
        log.debug(f"Running Python packer.")
        
        n_rows, n_cols = arr.shape
        words = (n_cols + 63) // 64
        out = np.zeros((n_rows, words), dtype=np.dtype("<u8"))
        for r in range(n_rows):
            row = arr[r]
            for w in range(words):
                base = w * 64
                limit = min(64, n_cols - base)
                word = 0
                ## vectorize this slice
                if limit > 0:
                    bits = row[base:base+limit].astype(np.uint64)
                    ## accumulate bits -> positions
                    for b in range(limit):
                        if bits[b]:
                            word |= (1 << b)
                out[r, w] = word
        
        log.debug(f"Finished packing matrix.")
        return out

    def _mk_bitpack_meta_from_dense(self, words: np.ndarray, headers: dict) -> dict:
        n_rows, words_per_row = words.shape
        n_cols = len(headers["alleles"])
        return {
            "format": "ardal.bitpack.v1",
            "dtype": "<u8",
            "endianness": "little",
            "row_major": True,
            "n_rows": n_rows,
            "n_cols": n_cols,
            "words_per_row": words_per_row,
            "bits_per_word": 64,
            "row_stride_bytes": words_per_row * 8,
            "data_file": None,            # in-memory; fill if you persist
            "data_nbytes": int(words.nbytes),
            "data_sha256": None
    }