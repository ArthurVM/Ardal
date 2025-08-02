import sys
import os
import pandas as pd
import numpy as np
import scipy as sp
import json

from ..exceptions.exceptions import *


class ArdalParser:
    """ Parses Ardal data from CSV, Parquet, or NPY/JSON file pairs.
    """

    def __init__( self,
                  input_data_structure,
                  file_format=None ):
        """ ArdalParser constructor
        """
        self.input_data = input_data_structure
        self.file_format = file_format.lower() if file_format else None
        self.matrix = None
        self.headers = None

        self._parse()


    def _parse( self ) -> int:
        """ Parses the data based on the specified file format.
        """

        ## check input data is not None
        if self.input_data is None:
            raise MalformedInputError("Input data structure cannot be None.")

        ## handle [np.ndarray, dict] or [dict, np.ndarray]
        if isinstance(self.input_data, list):
            if len(self.input_data) != 2:
                raise MalformedInputError("Input list must contain two elements: matrix and headers.")

            a, b = self.input_data
            if isinstance(a, np.ndarray) and isinstance(b, dict):
                self.matrix, self.headers = np.ascontiguousarray(a), b
                self._validate()
                return 0

            elif isinstance(b, np.ndarray) and isinstance(a, dict):
                self.matrix, self.headers = np.ascontiguousarray(b), a
                self._validate()
                return 0

            ## check for file-based input
            elif isinstance(a, str) and isinstance(b, str):
                if not os.path.exists(a) or not os.path.exists(b):
                    raise FileNotFoundError(f"One or more file paths do not exist: {a}, {b}")
                self.input_data = sorted([a, b], key=lambda x: x.split(".")[-1])
                self.file_format = self.input_data[1].split(".")[-1].lower()
            else:
                raise MalformedInputError("If list input, must contain either [headers, matrix] or two file paths.")

        elif isinstance(self.input_data, str):
            if not os.path.exists(self.input_data):
                raise FileNotFoundError(f"File does not exist: {self.input_data}")
            self.file_format = self.input_data.split(".")[-1].lower()

        if self.file_format not in ["csv", "parquet", "npy", "npz"]:
            raise UnsupportedFormatError("Unsupported file format: must be csv, parquet, npy, or npz")

        if self.file_format == "csv":
            self.matrix, self.headers = self._load_csv(self.input_data)
        elif self.file_format == "parquet":
            self.matrix, self.headers = self._load_parquet(self.input_data)
        elif self.file_format == "npy":
            self.matrix, self.headers = self._load_npy_pair(*self.input_data)
        elif self.file_format == "npz":
            self.matrix, self.headers = self._load_npz_pair(*self.input_data)

        else:
            raise MalformedInputError("Input must be list or string.")

        self._validate()


    def _validate(self):
        ## matrix checks
        if not isinstance(self.matrix, np.ndarray):
            raise LoadMatrixError("Matrix must be a NumPy array.")
        if self.matrix.ndim != 2:
            raise LoadMatrixError("Matrix must be 2-dimensional.")
        
        ## check for matrix memory contiguity
        if not self.matrix.flags['C_CONTIGUOUS']:
            try:
                self.matrix = np.ascontiguousarray(self.matrix)
                ## check again
                if not self.matrix.flags['C_CONTIGUOUS']:
                    raise LoadMatrixError(f"Failed to make array contiguous in memory : {e}")
            except Exception as e:
                raise LoadMatrixError(f"Failed to make array contiguous in memory : {e}")

        n_rows, n_cols = self.matrix.shape

        ## header structure
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


    def _load_csv(self, filepath):
        df = pd.read_csv(filepath, index_col=0)
        matrix = df.values
        headers = {"guids": list(df.index), "alleles": list(df.columns)}
        return matrix, headers


    def _load_parquet(self, filepath):
        df = pd.read_parquet(filepath, engine="fastparquet")
        matrix = df.values
        headers = {"guids": list(df.index), "alleles": list(df.columns)}
        return matrix, headers


    def _load_npy_pair(self, json_path, npy_path):
        try:
            matrix = np.ascontiguousarray(np.load(npy_path))
            with open(json_path, "r") as f:
                headers = json.load(f)
        except Exception as e:
            raise LoadMatrixError(f"Failed to load npy/json pair: {e}")
        return matrix, headers


    def _load_npz_pair(self, json_path, npz_path):
        try:
            # sparse = np.load(npz_path)['matrix']
            sparse = sp.sparse.load_npz(npz_path).toarray()
            matrix = np.ascontiguousarray(sparse)
            with open(json_path, "r") as f:
                headers = json.load(f)
        except Exception as e:
            raise LoadMatrixError(f"Failed to load npz/json pair: {e}")
        return matrix, headers