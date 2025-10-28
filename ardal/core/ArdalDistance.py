import os
import numpy as np
import pandas as pd
import time
from collections import defaultdict
from typing import Union, List, Dict

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
                  guids: Union[list, None] = None,
                  alleles: Union[list, None] = None,
                  intervals: Union[list, None] = None,
                  intervals_bed: Union[str, None] = None,
                  allele_coords_bed: Union[str, None] = None,
                  metric: str = "hamming",
                  use_simd: bool = True,
                  threads: int = 1,
                  backend: str = "auto",
                  allele_id_format: str = "{chr}.{start}.{ref}.{alt}",
                  *,
                  return_square: bool = False,
                  as_dataframe: bool = False,
                  ) -> Union[np.ndarray, pd.DataFrame]:
        """
        Compute pairwise distances.
        Can be Hamming (default), Jaccard, or Inner-Product distance.

        Returns:
            - If return_square=False (default): 1D condensed NumPy array of length n*(n-1)//2
            - If return_square=True:
                * NumPy (n x n) array if as_dataframe=False
                * pandas.DataFrame if as_dataframe=True

        Notes:
            - Large N strongly favors returning condensed. Expanding to square is O(n^2) memory and time.
            - I would refrain from asking a pandas.DataFrame to be built unless you are absolutely sure it is required.
              This can be even more expensive for large N.
        """
        naturalsize = require_package("humanize", attr="naturalsize")

        log.info("[PAIRWISE] Starting pairwise")

        ## function specific checks
        self._check_pairwise_args(
            guids=guids,
            alleles=alleles,
            intervals=intervals,
            intervals_bed=intervals_bed,
            allele_coords_bed=allele_coords_bed,
            metric=metric,
        )
        
        ## raise a warning if the user requests a lower triangle dataframe
        if not return_square and as_dataframe:
            log.warning(f"'return_square' cannot be False when 'as_dataframe' is True. Setting 'return_square' to True.")
            return_square = True

        kwargs = {
            "metric": metric,
            "use_simd": use_simd,
            "threads": threads,
            "backend": backend,
        }

        ## args storage
        local_flag = bool(guids or alleles or intervals)
        local_args = {}

        ## check input guids
        if guids:
            self._headerUtils.check_guids(guids)
            local_args["guids"] = guids
        else:
            guids = self._headerUtils.headers["guids"]
            local_args["guids"] = guids

        ## check input alleles
        if alleles:
            self._headerUtils.check_alleles(alleles)
            local_args["alleles"] = alleles
        else:
            local_args["alleles"] = self._headerUtils.headers["alleles"]

        ## get interval alleles
        if intervals:
            log.info("[PAIRWISE] Constructing intervals")
            alleles = self._headerUtils.get_interval_alleles(
                intervals=intervals,
                intervals_bed=intervals_bed,
                allele_coords_bed=allele_coords_bed,
                allele_id_format=allele_id_format,
            )
            local_args["alleles"] = alleles
            log.info("[PAIRWISE] Interval constructed.")

        ## ------ memory guard (based on requested output shape & real dtype) ------
        n = len(guids)
        out_dtype = self._metric_dtype(metric)
        if return_square:
            count = n * n
        else:
            count = n * (n - 1) // 2
        bytes_per = np.dtype(out_dtype).itemsize
        est_bytes = count * bytes_per
        natural_est = naturalsize(est_bytes, binary=True)

        total_memory = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        log.info(
            f"[PAIRWISE] Planning {n}x{n} {metric} distance matrix ({out_dtype}) with shape "
            f"{'square' if return_square else 'condensed'}, est ~{natural_est}."
        )
        if est_bytes > total_memory * 0.8:
            raise MemoryError(
                f"Requested output ({'nxn' if return_square else 'condensed'}) "
                f"would require ~{natural_est}, exceeding safe memory limits. "
                f"Either set return_square=False or subset your data."
            )
        ## ------------------------------------------------------------------------

        ## run compute
        if not local_flag:
            dist_condensed = self._pairwise_global(**kwargs)
        else:
            dist_condensed = self._pairwise_local(**local_args, **kwargs)

        ## cast once if the backend returns a different dtype
        ## could happen since Jaccard and Hamming use different dtypes
        if dist_condensed.dtype != out_dtype:
            dist_condensed = dist_condensed.astype(out_dtype, copy=False)

        ## return condensed by default
        if not return_square:
            return dist_condensed

        ## expand to square with correct dtype
        mat = self._expand_condensed(dist_condensed, n, out_dtype)

        ## construct dataframe
        ## NOTE: this could be extremely expensive for large matrices
        if as_dataframe:
            log.info("[PAIRWISE] Building pandas dataframe. For large N this could take some time.")
            return pd.DataFrame(mat, index=guids, columns=guids)
        return mat


    @staticmethod
    def _metric_dtype(metric: str):
        m = metric.lower()
        if m == "jaccard":
            return np.float32
        elif m in ("hamming", "inner_product"):
            return np.uint32
        elif m == "cosine":
            return np.float64
        raise ParameterError(f"Unknown metric: {metric}")


    @staticmethod
    def _expand_condensed(condensed: np.ndarray, n: int, dtype) -> np.ndarray:
        """
        Build nxn square array from condensed vector.
        """
        log.info("[PAIRWISE] Expanding the condensed matrix.")
        mat = np.zeros((n, n), dtype=dtype)
        iu, ju = np.triu_indices(n, 1)
        mat[iu, ju] = condensed
        mat[ju, iu] = condensed
        return mat


    def _pairwise_global( self,
                          metric: str = "hamming",
                          use_simd: bool = True,
                          threads: int = 1,
                          backend: str = "auto",
                          ) -> np.ndarray:
        """ Pairwise distance computation for a matrix.

        Args:
            metric (str, optional): the distance metric to use {jaccard/hamming/inner_product}. Defaults to "hamming".
            use_simd (bool, optional): whether to force SIMD. Defaults to True.
            threads (int, optional): number of threads to run computation on. Defaults to 1.
            backend (str, optional): which backend to use {auto/bit/roaring}. Defaults to "auto".

        Raises:
            ParameterError: raised if an unknown metric is passed to this function.

        Returns:
            np.ndarray: a condensed (lower triangle) distance matrix.
        """
        s = time.time()
        log.info(f"[PAIRWISE] Starting global {metric} distance calculations...")

        if metric == "jaccard":
            dist_tri = self._matrix.jaccard(use_simd=use_simd, threads=threads, backend=backend)
        elif metric == "hamming":
            dist_tri = self._matrix.hamming(use_simd=use_simd, threads=threads, backend=backend)
        elif metric == "inner_product":
            dist_tri = self._matrix.innerProduct(use_simd=use_simd, threads=threads, backend=backend)
        elif metric == "cosine":
            dist_tri = self._matrix.cosineDistance(use_simd=use_simd, threads=threads, backend=backend)
        else:
            raise ParameterError(f"Unknown metric: {metric}")

        log.info(f"[PAIRWISE] Finished {metric} distance calculations in {time.time()-s:.3f}s.")
        return dist_tri  ## condensed with correct dtype from backend


    def _pairwise_local( self,
                         guids: List,
                         alleles: List,
                         metric: str = "hamming",
                         use_simd: bool = True,
                         threads: int = 1,
                         backend: str = "auto",
                         ) -> np.ndarray:
        """ Pairwise distance computation for a local (subsetted) matrix.

        Args:
            guids (list): a set of guids to compute over.
            alleles (list): a set of alleles to compute over.
            metric (str, optional): the distance metric to use {jaccard/hamming/inner_product}. Defaults to "hamming".
            use_simd (bool, optional): whether to force SIMD. Defaults to True.
            threads (int, optional): number of threads to run computation on. Defaults to 1.
            backend (str, optional): which backend to use {auto/bit/roaring}. Defaults to "auto".

        Raises:
            ParameterError: raised if an unknown metric is passed to this function.

        Returns:
            np.ndarray: a condensed (lower triangle) distance matrix.
        """
        s = time.time()

        row_indices = sorted([self._headerUtils.encode_guid(g) for g in guids], reverse=False)
        col_indices = sorted([self._headerUtils.encode_allele(a) for a in alleles], reverse=False)

        log.info(
            f"[PAIRWISE] Starting local {metric} on {len(col_indices)}x{len(row_indices)} subset."
        )

        if metric == "jaccard":
            dist_tri = self._matrix.jaccard_subset(
                row_indices=row_indices, col_indices=col_indices, use_simd=use_simd, threads=threads, backend=backend
            )
        elif metric == "hamming":
            dist_tri = self._matrix.hamming_subset(
                row_indices=row_indices, col_indices=col_indices, use_simd=use_simd, threads=threads, backend=backend
            )
        elif metric == "inner_product":
            dist_tri = self._matrix.innerProduct_subset(
                row_indices=row_indices, col_indices=col_indices, use_simd=use_simd, threads=threads, backend=backend
            )
        elif metric == "cosine":
            dist_tri = self._matrix.cosineDistance_subset(
                row_indices=row_indices, col_indices=col_indices, use_simd=use_simd, threads=threads, backend=backend
            )
        else:
            raise ParameterError(f"Unknown metric: {metric}")

        log.info(f"[PAIRWISE] Finished local {metric} distance calculations in {time.time()-s:.3f}s.")
        return dist_tri
    
    
    def _check_pairwise_args( self,
                              guids : Union[list, None] = None,
                              alleles : Union[list, None] = None,
                              intervals : Union[list, None] = None,
                              intervals_bed : Union[str, None] = None,
                              allele_coords_bed : Union[str, None] = None,
                              metric : str = "hamming",
                              return_square : bool = False,
                              as_dataframe : bool = False
                              ) -> None:
        ACCEPTED_DIST_FUNCTIONS = ["jaccard", "hamming", "inner_product", "cosine"]
        
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
        
        

    def snv_neighbourhood( self,
                           guid : str,
                           n : int ) -> dict:
        """ find all GUIDs which lie within n SNVs of a given GUID
        WARNING : NOT PRODUCTION READY
        assumes allele ID of form {ref_nucleotide}{pos}{alt_nucleotide} and so the pos can be indexed out with [1:-1]
        """
        log.critical("snv_neighbourhood not production ready. May produce unstable results.")
        
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
                       backend : str = "auto",
                       metric : str = "hamming" ) -> dict:
        """ get the allele neighbourhood of a GUID using specified metric
        """
        validate_type(guid, str, "guid")
        validate_type(n, int, "n")
        validate_type(metric, str, "metric")
        if n < 0:
            raise ParameterError("n must be non-negative.")
        if n == 0:
            return {}
        metric_lower = metric.lower()
        acceptable_metrics = {"hamming", "inner_product", "innerproduct"}
        if metric_lower not in acceptable_metrics:
            raise ParameterError(f"metric '{metric}' not supported for neighbourhood. Choose from {sorted(acceptable_metrics)}.")
        
        self._headerUtils.check_guids([guid])

        guid_coord = self._headerUtils.encode_guid(guid)
        if metric_lower == "hamming":
            ncoords = self._matrix.neighbourhood(row_idx=guid_coord,
                                                 epsilon=n,
                                                 use_simd=use_simd,
                                                 threads=threads,
                                                 backend=backend)
            neighbourhood = {self._headerUtils.decode_guid(coord) : int(dist) for coord, dist in ncoords}
        else:
            if not hasattr(self._matrix, "innerProductNeighbourhood"):
                raise RuntimeError("inner_product neighbourhood requires a backend with innerProductNeighbourhood support.")
            ncoords = self._matrix.innerProductNeighbourhood(row_idx=guid_coord,
                                                             ip_epsilon=n,
                                                             use_simd=use_simd,
                                                             backend=backend)
            neighbourhood = {self._headerUtils.decode_guid(coord) : int(dist) for coord, dist in ncoords}

        return neighbourhood


    @check_backend_argument
    @check_thread_count
    @check_use_simd
    def knn( self,
             guid : str,
             k : int,
             use_simd : bool = True,
             threads : int = 1,
             backend : str = "auto",
             metric : str = "hamming" ) -> Dict:
        """knn
        Find the k nearest-neighbours

        Args:
            guid (str): _description_
            k (int): _description_
            use_simd (bool, optional): _description_. Defaults to True.
            threads (int, optional): _description_. Defaults to 1.
            backend (str, optional): _description_. Defaults to "auto".
            metric (str, optional): similarity/distance metric ('hamming', 'inner_product',
                'jaccard', 'cosine'). Defaults to "hamming".

        Returns:
            List: _description_
        """
        validate_type(guid, str, "guid")
        validate_type(k, int, "k")
        validate_type(metric, str, "metric")
        if k < 0:
            raise ParameterError("k must be non-negative.")
        if k == 0:
            return {}
        metric_lower = metric.lower()
        acceptable_metrics = {"hamming", "inner_product", "innerproduct", "jaccard", "cosine"}
        if metric_lower not in acceptable_metrics:
            raise ParameterError(f"metric '{metric}' not supported. Choose from {sorted(acceptable_metrics)}.")

        self._headerUtils.check_guids([guid])

        guid_coord = self._headerUtils.encode_guid(guid)

        def _legacy_knn(include_metric: bool):  # pragma: no cover
            try:
                if include_metric:
                    return self._matrix.knn(row_idx=guid_coord,
                                             k=k,
                                             metric=metric_lower,
                                             use_simd=use_simd,
                                             threads=threads,
                                             backend=backend)
                else:
                    return self._matrix.knn(row_idx=guid_coord,
                                             k=k,
                                             use_simd=use_simd,
                                             threads=threads,
                                             backend=backend)
            except TypeError as exc:
                raise RuntimeError(
                    f"knn metric '{metric_lower}' requires a rebuilt backend with metric-aware bindings."
                ) from exc

        if metric_lower == "hamming":
            if hasattr(self._matrix, "knn_hamming"):
                ncoords = self._matrix.knn_hamming(row_idx=guid_coord,
                                                   k=k,
                                                   use_simd=use_simd,
                                                   threads=threads,
                                                   backend=backend)
            else:
                ncoords = _legacy_knn(include_metric=False)
        elif metric_lower in ("inner_product", "innerproduct"):
            if hasattr(self._matrix, "knn_inner_product"):
                ncoords = self._matrix.knn_inner_product(row_idx=guid_coord,
                                                         k=k,
                                                         use_simd=use_simd,
                                                         threads=threads,
                                                         backend=backend)
            else:
                ncoords = _legacy_knn(include_metric=True)
        elif metric_lower == "jaccard":
            if hasattr(self._matrix, "knn_jaccard"):
                ncoords = self._matrix.knn_jaccard(row_idx=guid_coord,
                                                   k=k,
                                                   use_simd=use_simd,
                                                   threads=threads,
                                                   backend=backend)
            else:
                ncoords = _legacy_knn(include_metric=True)
        else:  # cosine
            if hasattr(self._matrix, "knn_cosine"):
                ncoords = self._matrix.knn_cosine(row_idx=guid_coord,
                                                  k=k,
                                                  use_simd=use_simd,
                                                  threads=threads,
                                                  backend=backend)
            else:
                ncoords = _legacy_knn(include_metric=True)
        neighbours = {self._headerUtils.decode_guid(coord) : dist for coord, dist in ncoords}

        return neighbours
