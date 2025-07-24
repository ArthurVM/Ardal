import pandas as pd
import json
import os
import numpy as np
from collections import defaultdict


# core/IO.py
class ArdalIO:

    def __init__(self, headers, hybrid_matrix, roaring_enabled):

        self._headers = headers
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled


    def toDataFrame( self ) -> pd.DataFrame:
        """ Return the allele matrix as a Pandas DataFrame.
        """
        return pd.DataFrame(self.__hybrid_matrix.getBitMatrix(), index=self.__headers["guids"], columns=self.__headers["alleles"])


    def toDict( self,
                force_bit_backend : bool = False ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        allele_dict = defaultdict(list)
        for guid_idx, guid_name in enumerate(self.__headers["guids"]):
            snp_indices = self.__hybrid_matrix.getSetBitIndices(guid_idx, force_bit_backend=force_bit_backend)
            for snp_idx in snp_indices:
                allele_dict[guid_name].append(self._decodeAllele(snp_idx))
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
            raise ValueError(f"Path '{file_path}' does not exist.")

        json_out_path = os.path.join(file_path, output_prefix + "_headers.json")
        matrix_out_path = os.path.join(file_path, output_prefix + "_matrix.npy")

        if os.path.exists(json_out_path):
            raise ValueError(f"File '{json_out_path}' already exists.")
        if os.path.exists(matrix_out_path):
            raise ValueError(f"File '{matrix_out_path}' already exists.")

        np.save(matrix_out_path, self.__hybrid_matrix.getMatrix())
        with open(os.path.join(file_path, output_prefix + "_headers.json")) as fout:
            json.dump(self.__headers, fout, indent=4)