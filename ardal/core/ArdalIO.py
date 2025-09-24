""" ArdalIO.py
This module provides functionality for reading and writing allele matrices in the Ardal framework.
"""
import os
import pathlib
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Union
from importlib.metadata import version, PackageNotFoundError

from ..utils.misc import require_package
from ..utils.decorators import check_backend_argument
from ..utils.exceptions import MatrixWriteError, ParameterError
from ..utils.logger import get_logger

log = get_logger()


## core/ArdalIO.py
class ArdalIO:
    """ Class for reading and writing allele matrices in the Ardal framework.
    """

    def __init__( self,
                  headerUtils,
                  hybrid_matrix,
                  roaring_enabled : bool ):

        self._headerUtils = headerUtils
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled


    def to_dataframe( self ) -> pd.DataFrame:
        """ Return the allele matrix as a Pandas DataFrame.
        """
        return pd.DataFrame(self._matrix.getBitMatrix(), index=self._headerUtils.headers["guids"], columns=self._headerUtils.headers["alleles"])


    @check_backend_argument
    def to_dict( self,
                 backend : str = "auto"
                 ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        allele_dict = defaultdict(list)
        for guid_idx, guid_name in enumerate(self._headerUtils.headers["guids"]):
            allele_indices = self._matrix.getSetBitIndices(guid_idx, backend=backend)
            for allele_idx in allele_indices:
                allele_dict[guid_name].append(self._headerUtils.decode_allele(allele_idx))
        return dict(allele_dict)
    

    def write( self,
               output_prefix : str,
               out_directory : str = "./",
               format : str = "bin",
               ) -> None:
        """ 
        Write the allele matrix to disk.
        Writes as a numpy/JSON pair.
        The format flag defines what 
        """
        """ Write the allele matrix to disk as a matrix/headers pair in specified format.
        Format must be one of:
        `npy` : outputs a dense uncompressed .npy matrix and .json headers.
        `npz` : outputs a dense compressed .npz matrix and .json headers.
        `bin` : outputs a bitpacked uncompressed .bin matrix with a .json containing headers and packing metadata.

        Args:
            file_path (str): _description_
            headers (str): _description_
            format (str) : _description_

        Returns:
            None
        """
        json = require_package("json", "json")
        
        ALLOWED_FORMATS = ["npy", "npz", "bin"]

        if not os.path.exists(out_directory):
            raise ParameterError(f"Path '{out_directory}' does not exist.")

        json_out_path = os.path.join(out_directory, output_prefix + ".json")
        
        if format not in ALLOWED_FORMATS:
            raise ParameterError(f"Invalid matrix output format : {format}. Must be one of {ALLOWED_FORMATS}.")
        
        matrix_out_path = os.path.join(out_directory, output_prefix + "." + format)

        if os.path.exists(json_out_path):
            raise MatrixWriteError(f"File '{json_out_path}' already exists.")
        if os.path.exists(matrix_out_path):
            raise MatrixWriteError(f"File '{matrix_out_path}' already exists.")

        matrix_to_save = self._matrix.getBitMatrix()
        if format == "npy":
            log.info("Writing dense matrix as .npy")
            np.save(matrix_out_path, matrix_to_save)
            headers_meta = self._headerUtils.headers
        elif format == "npz":
            log.info("Writing dense matrix as .npz")
            np.savez_compressed(matrix_out_path, matrix=matrix_to_save)
            headers_meta = self._headerUtils.headers
        elif format == "bin":
            log.info("Writing packed matrix as .bin")
            headers_meta = self._write_packed(matrix_out_path)
            
        log.info(f"Wrote allele matrix to disk : {matrix_out_path}")

        with open(json_out_path, 'w') as fout:
            json.dump(headers_meta, fout, indent=4)
            
        log.info(f"Wrote headers to disk : {json_out_path}")
        
        
    def _write_packed( self,
                       matrix_out_path : str ):
        arr = self._matrix.getPackedMatrix()
        
        ## ensure little endian
        if arr.dtype.byteorder not in ('<', '='):
            arr = arr.byteswap().newbyteorder('<')
        
        ## write bin
        bin_path = pathlib.Path(matrix_out_path)
        arr.tofile(bin_path)  ## raw little-endian uint64 words

        meta = self._mk_bitpack_meta_from_dense(arr, self._headerUtils.headers, matrix_out_path)
        headers_meta = { "meta" : meta, "headers" : self._headerUtils.headers }
        
        return headers_meta
        
        
    def _mk_bitpack_meta_from_dense( self,
                                     words : np.ndarray,
                                     headers : dict,
                                     data_file : Union[str, None] = None) -> dict:
        """ Compute metadata for writing a bitpacked binary matrix.

        Args:
            words (np.ndarray): _description_
            headers (dict): _description_

        Returns:
            dict: _description_
        """
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
            "data_file": data_file,
            "data_nbytes": int(words.nbytes),
            "data_sha256": None,
            "generated_by" : "Ardal v" + version("ardal"),
    }


    def make_fastas( self,
                     guids : list = [],
                     ref : Union[str, None] = None,
                     allele_id_format : str = "{ref}.{chr}.{start}.{alt}"
                     ) -> None:
        """ Takes a set of guids and a reference fasta and constructs a simulated fasta from each dataset using alleles stored
        in the allele matrix. The position of each allele is determined by the allele_id_format argument, whereby the keywords:
        ref, alt, chr, start, end, are used to outline the allele id naming convention.

        E.g. for the allele `A.chr1.100.101.T` would be decoded using the allele_id_format string `{ref}.{chr}.{start}.{end}.{alt}`.
        """
        raise NotImplementedError("makeFastas function not yet implemented.")
        return None