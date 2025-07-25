import numpy as np
import os
from scipy.spatial.distance import squareform
from humanize import naturalsize
import pandas as pd

from .utilities import *


# core/ArdalCompare.py
class ArdalCompare:

    def __init__(self, headers, hybrid_matrix, roaring_enabled):

        self._headers = headers
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled

    def pairwise( self,
                  metric: str = "hamming",
                  use_simd : bool=True,
                  threads : int=1,
                  force_bit_backend: bool=False ) -> pd.DataFrame:
        """ Calculates a pairwise distance matrix.
        Pairwise distance can be calculated using Hamming, Jaccard, or Inner Product (number of shared SNPs) functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """
        if not isinstance(metric, str):
            raise ValueError("metric must be a string.")
        if not isinstance(use_simd, bool):
            raise ValueError("use_simd must be a boolean.")
        if not isinstance(threads, int):
            raise ValueError("threads must be an integer.")
        if not isinstance(force_bit_backend, bool):
            raise ValueError("force_bit_backend must be a boolean.")
        if threads < 1:
            raise ValueError("threads must be at least 1.")

        ## check there is enough memory to store the pairwise matrix
        mat_size = len(self._headers["guids"])**2
        total_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        if mat_size * 8 > total_memory * 0.8: ## 8 bytes per float64, 80% of total memory
            raise MemoryError(f"Pairwise distance matrix of scale {len(self._headers['guids'])}x{len(self._headers['guids'])} will requre {naturalsize(mat_size * 8, binary=True)} memory and will exceed system memory limits. Please subset your data.")

        ## check the specified distance function is valid
        accepted_dist_functions = ["jaccard", "hamming", "inner_product"]
        if metric not in accepted_dist_functions:
            raise ValueError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")
        
        ## calculate the distance matrix using _ardal
        if metric == "jaccard":
            dist_tri = self._matrix.jaccard(use_simd=use_simd,
                                            threads=threads,
                                            force_bit_backend=force_bit_backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.float32)

        elif metric == "hamming":
            dist_tri = self._matrix.hamming(use_simd=use_simd,
                                            threads=threads,
                                            force_bit_backend=force_bit_backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)

        elif metric == "inner_product":
            dist_tri = self._matrix.innerProduct(use_simd=use_simd,
                                                 threads=threads,
                                                 force_bit_backend=force_bit_backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)
        
        dist_df = pd.DataFrame(dist_matrix, columns=self._headers["guids"], index=self._headers["guids"])
        
        return dist_df
    

    def snvNeighbourhood( self,
                          guid : str,
                          n : int ) -> dict:
        """ find all GUIDs which lie within n SNVs of a given GUID
        WARNING : NOT PRODUCTION READY
        assumes allele ID of form {ref_nucleotide}{pos}{alt_nucleotide} and so the pos can be indexed out with [1:-1]
        """
        checkGUIDs([guid], self._headers)

        snv_neighbourhood = {}

        guid_coord = encodeGuid(guid, self._headers)
        guid_snv_positions = set([decodeAllele(allele_coord, self._headers)[1:-1] for allele_coord in self._matrix.getSetBitIndices(guid_coord)])

        for guid_idx in range(len(self._headers["guids"])):
            if guid_idx == guid_coord:
                continue
            other_snv_positions = set([decodeAllele(allele_coord, self._headers)[1:-1] for allele_coord in self._matrix.getSetBitIndices(guid_idx)])
            snv_dist = guid_snv_positions ^ other_snv_positions
            if len(snv_dist) <= n:
                snv_neighbourhood[decodeGuid(guid_idx, self._headers)] = len(snv_dist)
        
        return snv_neighbourhood
    

    def neighbourhood( self,
                       guid : str,
                       n : int,
                       use_simd : bool = True,
                       threads : int = 1,
                       force_bit_backend: bool=False ) -> dict:
        """ get the allele neighbourhood of a GUID using hamming distance
        """

        if not isinstance(guid, str):
            raise ValueError("guid must be a string.")
        if not isinstance(n, int):
            raise ValueError("n must be an integer.")
        if n < 0:
            raise ValueError("n must be non-negative.")
        if n == 0:
            return {}
        
        checkGUIDs([guid], self._headers)

        guid_coord = encodeGuid(guid, self._headers)
        ncoords = self._matrix.neighbourhood(row_idx=guid_coord,
                                             epsilon=n,
                                             use_simd=use_simd,
                                             threads=threads,
                                             force_bit_backend=force_bit_backend)
        neighbourhood = {decodeGuid(coord, self._headers) : hdist for coord, hdist in ncoords}

        return neighbourhood
