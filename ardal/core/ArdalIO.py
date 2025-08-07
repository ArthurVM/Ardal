""" ArdalIO.py
This module provides functionality for reading and writing allele matrices in the Ardal framework.
"""
import os
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Union

from ..utils.misc import require_package
from ..utils.decorators import check_backend_argument
from ..utils.exceptions import MatrixWriteError
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
               file_path : str,
               output_prefix : str,
               npz : bool = False
               ) -> None:
        """ Write the allele matrix to disk.
        Writes as a numpy/JSON pair.
        The npz flag writes the numpy matrix as a compressed npz.
        """
        json = require_package("json", "json")

        if not os.path.exists(file_path):
            raise MatrixWriteError(f"Path '{file_path}' does not exist.")

        json_out_path = os.path.join(file_path, output_prefix + "_headers.json")
        matrix_ext = ".npz" if npz else ".npy"
        matrix_out_path = os.path.join(file_path, output_prefix + "_matrix" + matrix_ext)

        if os.path.exists(json_out_path):
            raise MatrixWriteError(f"File '{json_out_path}' already exists.")
        if os.path.exists(matrix_out_path):
            raise MatrixWriteError(f"File '{matrix_out_path}' already exists.")

        matrix_to_save = self._matrix.getBitMatrix()
        if npz:
            np.savez_compressed(matrix_out_path, matrix=matrix_to_save)
        else:
            np.save(matrix_out_path, matrix_to_save)
            
        log.info(f"Wrote allele matrix to disk : {matrix_out_path}")

        with open(json_out_path, 'w') as fout:
            json.dump(self._headerUtils.headers, fout, indent=4)
            
        log.info(f"Wrote headers to disk : {json_out_path}")


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