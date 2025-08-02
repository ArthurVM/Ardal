import numpy as np
import os
from collections import defaultdict
import pandas as pd

from .utilities import *
from ..exceptions.exceptions import *


# core/ArdalCompare.py
class ArdalCompare:

    def __init__(self, headers, hybrid_matrix, roaring_enabled):

        self._headers = headers
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled

    def pairwise( self,
                  guid_ids : list = None,
                  allele_ids : list = None,
                  intervals : list = None,
                  coords_bed : str = None,
                  metric: str = "hamming",
                  use_simd : bool=True,
                  threads : int=1,
                  force_bit_backend: bool=False,
                  allele_id_format : str = "{chr}.{start}.{ref}.{alt}") -> pd.DataFrame:
        """ Calculates the distance between pairs of guids. If guids_ids are provided then it will only compute distances between
        specified guids. If allele_ids are provided then it will only compute using specified alleles. If intervals are provided then
        it will find alleles which fall within the given genomic intervals and only compute using these alleles.

        Pairwise distance can be calculated using Hamming, Jaccard, or Inner Product (number of shared SNPs) functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """    
        naturalsize = require_package("humanize", attr="naturalsize")
    
        if not isinstance(metric, str):
            raise TypeError("metric must be a string.")
        if not isinstance(use_simd, bool):
            raise TypeError("use_simd must be a boolean.")
        if not isinstance(threads, int):
            raise TypeError("threads must be an integer.")
        if not isinstance(force_bit_backend, bool):
            raise ParameterError("force_bit_backend must be a boolean.")
        if threads < 1:
            raise ParameterError("threads must be at least 1.")
        
        kwargs = { "metric" : metric,
                   "use_simd" : use_simd,
                   "threads" : threads,
                   "force_bit_backend" : force_bit_backend }
        
        ## specify whether to run local mode
        local_flag = False

        ## local params
        if allele_ids and intervals:
            raise ParameterError("Cannot run pairwise using both allele_ids and intervals kwargs. Please provide only one.")
        
        if coords_bed and not intervals:
            raise ParameterError("Please provide a set of intervals when using bed allele coordinates.")

        local_args = defaultdict(list)
        
        ## check for local mode
        if guid_ids or allele_ids or intervals:
            local_flag = True

        ## get guids to run on
        if guid_ids:
            checkGUIDs(guid_ids, self._headers)
            local_args["guid_ids"] = guid_ids
        else:
            guid_ids = self._headers['guids']
            local_args["guid_ids"] = guid_ids

        ## get alleles to run on
        if allele_ids:
            checkAlleles(allele_ids, self._headers)
            local_args["allele_ids"] = allele_ids
        else:
            local_args["allele_ids"] = self._headers['alleles']

        if intervals:
            allele_ids = getIntervalAlleles(intervals, self._headers, allele_id_format, coords_bed)
            local_args["allele_ids"] = allele_ids

        ## check there is enough memory to store the pairwise matrix
        mat_size = len(guid_ids)**2
        total_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        if mat_size * 8 > total_memory * 0.8: ## 8 bytes per float64, 80% of total memory
            raise MemoryError(f"Pairwise distance matrix of scale {len(guid_ids)}x{len(guid_ids)} will requre {naturalsize(mat_size * 8, binary=True)} memory and will exceed system memory limits. Please subset your data.")

        ## check the specified distance function is valid
        accepted_dist_functions = ["jaccard", "hamming", "inner_product"]
        if metric not in accepted_dist_functions:
            raise ParameterError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")
        
        if not local_flag:
            return self._pairwise_global( **kwargs )

        if local_flag:
            return self._pairwise_local( **local_args,
                                         **kwargs )
    

    def _pairwise_global( self,
                         metric: str = "hamming",
                         use_simd : bool=True,
                         threads : int=1,
                         force_bit_backend: bool=False ) -> pd.DataFrame:
        squareform = require_package("scipy", "scipy.spatial.distance", attr="squareform")
        
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
        
        return pd.DataFrame(dist_matrix, columns=self._headers["guids"], index=self._headers["guids"])


    def _pairwise_local( self,
                         guid_ids : list = None,
                         allele_ids : list = None,
                         metric: str = "hamming",
                         use_simd : bool=True,
                         threads : int=1,
                         force_bit_backend: bool=False ) -> pd.DataFrame:
        squareform = require_package("scipy", "scipy.spatial.distance", attr="squareform")
        
        row_indices = [encodeGuid(g, self._headers) for g in guid_ids]
        col_indices = [encodeAllele(a, self._headers) for a in allele_ids]
        
        if metric == "jaccard":
            dist_tri = self._matrix.jaccard_subset(row_indices=row_indices,
                                                   col_indices=col_indices,
                                                   use_simd=use_simd,
                                                   threads=threads,
                                                   force_bit_backend=force_bit_backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.float32)

        elif metric == "hamming":
            dist_tri = self._matrix.hamming_subset(row_indices=row_indices,
                                                   col_indices=col_indices,
                                                   use_simd=use_simd,
                                                   threads=threads,
                                                   force_bit_backend=force_bit_backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)

        elif metric == "inner_product":
            dist_tri = self._matrix.innerProduct_subset(row_indices=row_indices,
                                                        col_indices=col_indices,
                                                        use_simd=use_simd,
                                                        threads=threads,
                                                        force_bit_backend=force_bit_backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)
        
        return pd.DataFrame(dist_matrix, columns=guid_ids, index=guid_ids)
    

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
            raise TypeError("guid must be a string.")
        if not isinstance(n, int):
            raise TypeError("n must be an integer.")
        if n < 0:
            raise ParameterError("n must be non-negative.")
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
