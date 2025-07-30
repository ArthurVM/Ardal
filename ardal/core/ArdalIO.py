""" ArdalIO.py
This module provides functionality for reading and writing allele matrices in the Ardal framework.
"""

import pandas as pd
import json
import os
import numpy as np
from collections import defaultdict

from .utilities import *
from ..exceptions.exceptions import *


# core/ArdalIO.py
class ArdalIO:
    """ Class for reading and writing allele matrices in the Ardal framework.
    """

    def __init__(self, headers, hybrid_matrix, roaring_enabled):

        self._headers = headers
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled


    def toDataFrame( self ) -> pd.DataFrame:
        """ Return the allele matrix as a Pandas DataFrame.
        """
        return pd.DataFrame(self._matrix.getBitMatrix(), index=self._headers["guids"], columns=self._headers["alleles"])


    def toDict( self,
                force_bit_backend : bool = False ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        allele_dict = defaultdict(list)
        for guid_idx, guid_name in enumerate(self._headers["guids"]):
            snp_indices = self._matrix.getSetBitIndices(guid_idx, force_bit_backend=force_bit_backend)
            for snp_idx in snp_indices:
                allele_dict[guid_name].append(decodeAllele(snp_idx, self._headers))
        return dict(allele_dict)
    

    def write( self,
               file_path : str,
               output_prefix : str,
               npz : bool = False ) -> None:
        """ Write the allele matrix to disk.
        Writes as a numpy/JSON pair.
        The npz flag writes the numpy matrix as a compressed npz.
        """
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

        with open(json_out_path, 'w') as fout:
            json.dump(self._headers, fout, indent=4)


    def makeFastas( self,
                    guids : list = [],
                    ref : str = None,
                    allele_id_format : str = "{ref}.{chr}.{start}.{alt}" ) -> None:
        """ Takes a set of guids and a reference fasta and constructs a simulated fasta from each dataset using SNPs stored
        in the allele matrix. The position of each SNP is determined by the allele_id_format argument, whereby the keywords:
        ref, alt, chr, start, end, are used to outline the allele id naming convention.

        E.g. for the allele `A.chr1.100.101.T` would be decoded using the allele_id_format string `{ref}.{chr}.{start}.{end}.{alt}`.
        """
        None