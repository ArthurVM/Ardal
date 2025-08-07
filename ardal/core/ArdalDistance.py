import os
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Union, List

from ..utils.misc import *
from ..utils.exceptions import *
from ..utils.decorators import check_backend_argument, check_allele_id_format, check_thread_count, check_use_simd, check_guids_list, check_alleles_list
from ..utils.validators import validate_type
from ..utils.logger import get_logger

log = get_logger()


# core/ArdalDistance.py
class ArdalDistance:

    def __init__( self,
                  headerUtils,
                  hybrid_matrix,
                  roaring_enabled : bool ):

        self._headerUtils = headerUtils
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled

    
    @check_backend_argument
    @check_allele_id_format
    @check_thread_count
    @check_use_simd
    @check_guids_list
    @check_alleles_list
    def pairwise( self,
                  guids : Union[list, None] = None,
                  alleles : Union[list, None] = None,
                  intervals : Union[list, None] = None,
                  intervals_bed : Union[str, None] = None,
                  allele_coords_bed : Union[str, None] = None,
                  metric: str = "hamming",
                  use_simd : bool = True,
                  threads : int = 1,
                  backend: str = "auto",
                  allele_id_format : str = "{chr}.{start}.{ref}.{alt}"
                  ) -> Union[pd.DataFrame, None]:
        """ Calculates the distance between pairs of guids. If guids_ids are provided then it will only compute distances between
        specified guids. If alleles are provided then it will only compute using specified alleles. If intervals are provided then
        it will find alleles which fall within the given genomic intervals and only compute using these alleles.

        Pairwise distance can be calculated using Hamming, Jaccard, or Inner Product (number of shared alleles) functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """    
        naturalsize = require_package("humanize", attr="naturalsize")
        time = require_package("time", "time")
        
        log.debug("[PAIRWISE] Starting pairwise")
        
        ## handle function specific parameter checking logic
        self._check_pairwise_args(guids=guids,
                                  alleles=alleles,
                                  intervals=intervals,
                                  intervals_bed=intervals_bed,
                                  allele_coords_bed=allele_coords_bed,
                                  metric=metric)
        
        kwargs = { "metric" : metric,
                   "use_simd" : use_simd,
                   "threads" : threads,
                   "backend" : backend }
        
        ## specify whether to run local mode
        local_flag = False
        local_args = defaultdict(list)
        
        ## check for local mode
        if guids or alleles or intervals:
            local_flag = True

        ## get guids to run on
        if guids:
            self._headerUtils.check_guids(guids)
            local_args["guids"] = guids
        else:
            guids = self._headerUtils.headers['guids']
            local_args["guids"] = guids

        ## get alleles to run on
        if alleles:
            self._headerUtils.check_alleles(alleles)
            local_args["alleles"] = alleles
        else:
            local_args["alleles"] = self._headerUtils.headers['alleles']

        if intervals:
            log.debug(f"[PAIRWISE] Constructing intervals")
            alleles = self._headerUtils.get_interval_alleles(intervals=intervals,
                                                             intervals_bed=intervals_bed,
                                                             allele_coords_bed=allele_coords_bed,
                                                             allele_id_format=allele_id_format)
            local_args["alleles"] = alleles
            log.debug(f"[PAIRWISE] Interval constructed.")
        
        ## check there is enough memory to store the pairwise matrix
        mat_size = len(guids)**2
        natural_mat_size = naturalsize(mat_size * 8, binary=True)
        total_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        log.debug(f"[PAIRWISE] Computing {metric} distance matrix of scale {len(guids)}x{len(guids)} using {natural_mat_size} memory.")
        if mat_size * 8 > total_memory * 0.8: ## 8 bytes per float64, 80% of total memory
            raise MemoryError(f"Pairwise distance matrix of scale {len(guids)}x{len(guids)} will requre {natural_mat_size} memory and will exceed system memory limits. Please subset your data.")
                        
        if not local_flag:
            return self._pairwise_global( **kwargs )

        if local_flag:
            return self._pairwise_local( **local_args,
                                         **kwargs )
            
    
    def _check_pairwise_args( self,
                              guids : Union[list, None] = None,
                              alleles : Union[list, None] = None,
                              intervals : Union[list, None] = None,
                              intervals_bed : Union[str, None] = None,
                              allele_coords_bed : Union[str, None] = None,
                              metric: str = "hamming"
                              ) -> None:
        ACCEPTED_DIST_FUNCTIONS = ["jaccard", "hamming", "inner_product"]
        
        ## check the specified distance function is valid
        validate_type(metric, str, "metric")
        metric_lower = metric.lower()
        if metric_lower not in ACCEPTED_DIST_FUNCTIONS:
            raise ParameterError(f"{metric_lower} not an accepted distance function. Must be one of {ACCEPTED_DIST_FUNCTIONS}.")

        if alleles and intervals:
            raise ParameterError("alleles and intervals arguments are mutually exclusive.")
        
        if intervals and intervals_bed:
            raise ParameterError("intervals and intervals_bed arguments are mutually exclusive.")
        
        if allele_coords_bed and not intervals:
            raise ParameterError("intervals argument cannot be None when allele_coords_bed argument is not None.")


    def _pairwise_global( self,
                          metric: str = "hamming",
                          use_simd : bool=True,
                          threads : int=1,
                          backend : str="auto" ) -> pd.DataFrame:
        squareform = require_package("scipy", "scipy.spatial.distance", attr="squareform")
        
        log.debug(f"[PAIRWISE] Starting global {metric} distance calculations...")
        
        if metric == "jaccard":
            dist_tri = self._matrix.jaccard(use_simd=use_simd,
                                            threads=threads,
                                            backend=backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.float32)

        elif metric == "hamming":
            dist_tri = self._matrix.hamming(use_simd=use_simd,
                                            threads=threads,
                                            backend=backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)

        elif metric == "inner_product":
            dist_tri = self._matrix.innerProduct(use_simd=use_simd,
                                                 threads=threads,
                                                 backend=backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)
        
        log.debug(f"[PAIRWISE] Finished {metric} distance calculations.")
        return pd.DataFrame(dist_matrix, columns=self._headerUtils.headers["guids"], index=self._headerUtils.headers["guids"])


    def _pairwise_local( self,
                         guids : list = None,
                         alleles : list = None,
                         metric: str = "hamming",
                         use_simd : bool = True,
                         threads : int = 1,
                         backend: str = "auto" ) -> pd.DataFrame:
        squareform = require_package("scipy", "scipy.spatial.distance", attr="squareform")
        time = require_package("time", "time")
        
        s = time.time()
                
        row_indices = [self._headerUtils.encode_guid(g) for g in guids]
        col_indices = [self._headerUtils.encode_allele(a) for a in alleles]
        
        log.debug(f"[PAIRWISE] Starting local {metric} distance calculations on {len(col_indices)}x{len(row_indices)} subset matrix...")
                
        if metric == "jaccard":
            dist_tri = self._matrix.jaccard_subset(row_indices=row_indices,
                                                   col_indices=col_indices,
                                                   use_simd=use_simd,
                                                   threads=threads,
                                                   backend=backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.float32)

        elif metric == "hamming":
            dist_tri = self._matrix.hamming_subset(row_indices=row_indices,
                                                   col_indices=col_indices,
                                                   use_simd=use_simd,
                                                   threads=threads,
                                                   backend=backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)

        elif metric == "inner_product":
            dist_tri = self._matrix.innerProduct_subset(row_indices=row_indices,
                                                        col_indices=col_indices,
                                                        use_simd=use_simd,
                                                        threads=threads,
                                                        backend=backend)
            dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)
        
        log.debug(f"[PAIRWISE] Finished {metric} distance calculations.")
        return pd.DataFrame(dist_matrix, columns=guids, index=guids)
    

    def snv_neighbourhood( self,
                           guid : str,
                           n : int ) -> dict:
        """ find all GUIDs which lie within n SNVs of a given GUID
        WARNING : NOT PRODUCTION READY
        assumes allele ID of form {ref_nucleotide}{pos}{alt_nucleotide} and so the pos can be indexed out with [1:-1]
        """
        log.critical("snvNeighbourhood not production ready. May produce unstable results.")
        
        validate_type(guid, str, "guid")
        validate_type(n, int, "n")
        if n < 0:
            raise ParameterError("n must be non-negative.")
        if n == 0:
            return {}
        
        self._headerUtils.check_guids([guid])

        snv_neighbourhood = {}

        guid_coord = self._headerUtils.encode_guid(guid)
        guid_snv_positions = set([self._headerUtils.decode_allele(allele_coord)[1:-1] for allele_coord in self._matrix.getSetBitIndices(guid_coord)])

        for guid_idx in range(len(self._headerUtils.headers["guids"])):
            if guid_idx == guid_coord:
                continue
            other_snv_positions = set([self._headerUtils.decode_allele(allele_coord)[1:-1] for allele_coord in self._matrix.getSetBitIndices(guid_idx)])
            snv_dist = guid_snv_positions ^ other_snv_positions
            if len(snv_dist) <= n:
                snv_neighbourhood[self._headerUtils.decode_guid(guid_idx)] = len(snv_dist)
        
        return snv_neighbourhood
    

    @check_backend_argument
    @check_thread_count
    @check_use_simd
    def neighbourhood( self,
                       guid : str,
                       n : int,
                       use_simd : bool = True,
                       threads : int = 1,
                       backend : str = "auto" ) -> dict:
        """ get the allele neighbourhood of a GUID using hamming distance
        """

        validate_type(guid, str, "guid")
        validate_type(n, int, "n")
        if n < 0:
            raise ParameterError("n must be non-negative.")
        if n == 0:
            return {}
        
        self._headerUtils.check_guids([guid])

        guid_coord = self._headerUtils.encode_guid(guid)
        ncoords = self._matrix.neighbourhood(row_idx=guid_coord,
                                             epsilon=n,
                                             use_simd=use_simd,
                                             threads=threads,
                                             backend=backend)
        neighbourhood = {self._headerUtils.decode_guid(coord) : hdist for coord, hdist in ncoords}

        return neighbourhood
