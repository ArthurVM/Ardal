""" ArdalGet.py
This module provides functionality for retrieving and manipulating allele matrices in the Ardal framework.
"""
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import Union, Tuple, List, TYPE_CHECKING

from ..utils.misc import require_package
from ..utils.decorators import check_thread_count, check_alleles_list, check_guids_list
from ..utils.exceptions import ParameterError, RoaringError
from ..utils.logger import get_logger
from ..utils.make_meta import make_meta

log = get_logger()


## core/ArdalGet.py
class ArdalGet:
    """ Class for retrieving and manipulating allele matrices in the Ardal framework.
    """

    def __init__( self,
                  headerUtils,
                  hybrid_matrix,
                  roaring_enabled : bool ):

        self._headerUtils = headerUtils
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled

    
    def _subset_chunk( self,
                       guids: List[str],
                       guid_indices: List[int],
                       alleles: List[str],
                       allele_indices: List[int],
                       data_only: bool,
                       threads: int,
                       child_verbosity: str,
                       child_quiet_init: bool
                       ) -> Union[Tuple[np.ndarray, dict], "Ardal"]:
        """Packed subset helper used for chunked extraction."""
        from ..Ardal import Ardal

        rows_only_flag = len(allele_indices) == 0

        if rows_only_flag:
            sub_matrix = self._matrix.getSubsetPackedMatrix_rows(guid_indices, threads)
        else:
            sub_matrix = self._matrix.getSubsetPackedMatrix(guid_indices, allele_indices, threads)
        if sub_matrix.dtype.byteorder not in ('<', '='):
            sub_matrix = sub_matrix.byteswap().newbyteorder('<')
        sub_matrix = np.asarray(sub_matrix, dtype=np.uint64, copy=False)
        if sub_matrix.ndim == 1:
            sub_matrix = sub_matrix.reshape(1, sub_matrix.shape[0])

        headers_meta: dict = {"headers": {"guids": guids, "alleles": alleles}}

        if self._headerUtils.has_missing_mask():
            per_guid_masks = {}
            if rows_only_flag:
                for guid in guids:
                    per_guid_masks[guid] = list(self._headerUtils.get_guid_missing_mask(guid))
            else:
                idx_map = {orig: new for new, orig in enumerate(allele_indices)}
                for guid in guids:
                    cols = self._headerUtils.get_guid_missing_mask(guid)
                    remapped = []
                    for c in cols:
                        if c in idx_map:
                            remapped.append(idx_map[c])
                    per_guid_masks[guid] = remapped
            headers_meta["column_masks"] = per_guid_masks

        headers_meta["meta"] = make_meta(sub_matrix,
                                         headers_meta["headers"],
                                         generated_by="ardal::subset",
                                         format_name="ardal.bitpack.v1",
                                         matrix_file=None)

        if not data_only:
            return Ardal(data_source=[sub_matrix, headers_meta],
                         allele_id_format=self._headerUtils._allele_id_format,
                         allele_positions=self._headerUtils.get_allele_positions(),
                         verbosity=child_verbosity,
                         quiet_init=child_quiet_init)
        return [sub_matrix, headers_meta]

    
    @check_thread_count
    @check_guids_list
    @check_alleles_list
    def subset( self,
                guids : List[str] = [],
                alleles : List[str] = [],
                drop_zero_cols : bool = True,
                chunk_size : Union[int, None] = None,
                yield_chunks : bool = False,
                data_only : bool = False,
                threads : int = 1,
                child_verbosity : str = "silent",
                child_quiet_init : bool = True
                ) -> Union[Tuple[np.ndarray, dict], "Ardal", List]:
        """ Take a list of GUIDs and subset the allele matrix to include only these GUIDs, allowing for standard operations.
        drop_zero_cols removes alleles with zero frequency in the selected GUIDs.
        chunk_size enables packed column chunking to reduce peak memory.
        yield_chunks returns a generator over chunked outputs.
        Returns a numpy matrix/JSON pair for feeding into Ardal.
        """
        ## check input
        if len(guids) == 0 and len(alleles) == 0:
            raise ParameterError("guids and alleles cannot both be empty.")

        ## check GUIDs
        if guids:
            self._headerUtils.check_guids(guids)
        else:
            guids = list(self._headerUtils.headers["guids"])

        guid_indices = sorted([self._headerUtils.encode_guid(guid) for guid in guids], reverse=False)

        ## check alleles
        if alleles:
            self._headerUtils.check_alleles(alleles)
            rows_only_flag = False
        else:
            ## default to all alleles when no explicit subset requested
            alleles = list(self._headerUtils.headers["alleles"])
            rows_only_flag = True
        
        ## remove columns which have zero mass in the subset matrix
        if drop_zero_cols:
            col_freqs = self._matrix.colFrequency(guid_indices)
            if rows_only_flag:
                alleles = [
                    allele
                    for allele, freq in zip(self._headerUtils.headers["alleles"], col_freqs)
                    if freq > 0.0
                ]
                rows_only_flag = False
            else:
                allele_indices = [self._headerUtils.encode_allele(allele) for allele in alleles]
                alleles = [
                    allele
                    for allele, idx in zip(alleles, allele_indices)
                    if col_freqs[idx] > 0.0
                ]
            
        log.info(f"Subsetting the {self._headerUtils.n_guids}x{self._headerUtils.n_alleles} matrix to {len(guids)}x{len(alleles)}")

        if rows_only_flag:
            allele_indices = []
        else:
            allele_indices = sorted([self._headerUtils.encode_allele(allele) for allele in alleles], reverse=False)
            if alleles:
                index_to_allele = {self._headerUtils.encode_allele(allele): allele for allele in alleles}
                alleles = [index_to_allele[idx] for idx in allele_indices]

        if chunk_size is not None and not rows_only_flag and allele_indices:
            if chunk_size <= 0:
                raise ParameterError("chunk_size must be a positive integer.")

            def _iter_chunks():
                for start in range(0, len(allele_indices), chunk_size):
                    chunk_indices = allele_indices[start:start + chunk_size]
                    chunk_alleles = alleles[start:start + chunk_size]
                    yield self._subset_chunk(guids=guids,
                                             guid_indices=guid_indices,
                                             alleles=chunk_alleles,
                                             allele_indices=chunk_indices,
                                             data_only=data_only,
                                             threads=threads,
                                             child_verbosity=child_verbosity,
                                             child_quiet_init=child_quiet_init)

            if yield_chunks:
                return _iter_chunks()
            n_rows = len(guid_indices)
            n_cols = len(allele_indices)
            wpr = (n_cols + 63) // 64
            out_matrix = np.zeros((n_rows, wpr), dtype=np.uint64)
            col_offset = 0

            for start in range(0, len(allele_indices), chunk_size):
                chunk_indices = allele_indices[start:start + chunk_size]
                if not chunk_indices:
                    continue
                chunk_packed = self._matrix.getSubsetPackedMatrix(guid_indices, chunk_indices, threads)
                if chunk_packed.dtype.byteorder not in ('<', '='):
                    chunk_packed = chunk_packed.byteswap().newbyteorder('<')
                chunk_packed = np.asarray(chunk_packed, dtype=np.uint64, copy=False)
                if chunk_packed.ndim == 1:
                    chunk_packed = chunk_packed.reshape(1, chunk_packed.shape[0])

                word_offset = col_offset // 64
                bit_offset = col_offset % 64
                if bit_offset == 0:
                    out_matrix[:, word_offset:word_offset + chunk_packed.shape[1]] |= chunk_packed
                else:
                    lo = np.left_shift(chunk_packed, np.uint64(bit_offset))
                    out_matrix[:, word_offset:word_offset + lo.shape[1]] |= lo
                    carry = np.right_shift(chunk_packed, np.uint64(64 - bit_offset))
                    dst_start = word_offset + 1
                    dst_end = min(dst_start + carry.shape[1], out_matrix.shape[1])
                    if dst_start < dst_end:
                        out_matrix[:, dst_start:dst_end] |= carry[:, :dst_end - dst_start]

                col_offset += len(chunk_indices)

            if wpr > 0:
                tail_bits = n_cols % 64
                if tail_bits:
                    tail_mask = (np.uint64(1) << np.uint64(tail_bits)) - np.uint64(1)
                    out_matrix[:, -1] &= tail_mask

            sub_matrix = out_matrix
            headers_meta: dict = {"headers": {"guids": guids, "alleles": alleles}}

            if self._headerUtils.has_missing_mask():
                per_guid_masks = {}
                idx_map = {orig: new for new, orig in enumerate(allele_indices)}
                for guid in guids:
                    cols = self._headerUtils.get_guid_missing_mask(guid)
                    remapped = []
                    for c in cols:
                        if c in idx_map:
                            remapped.append(idx_map[c])
                    per_guid_masks[guid] = remapped
                headers_meta["column_masks"] = per_guid_masks

            headers_meta["meta"] = make_meta(sub_matrix,
                                             headers_meta["headers"],
                                             generated_by="ardal::subset",
                                             format_name="ardal.bitpack.v1",
                                             matrix_file=None)

            if not data_only:
                from ..Ardal import Ardal
                return Ardal(data_source=[sub_matrix, headers_meta],
                             allele_id_format=self._headerUtils._allele_id_format,
                             allele_positions=self._headerUtils.get_allele_positions(),
                             verbosity=child_verbosity,
                             quiet_init=child_quiet_init)
            return [sub_matrix, headers_meta]

        sub_matrix = self._matrix.getSubsetPackedMatrix_rows(guid_indices, threads)
        if sub_matrix.dtype.byteorder not in ('<', '='):
            sub_matrix = sub_matrix.byteswap().newbyteorder('<')
        sub_matrix = np.asarray(sub_matrix, dtype=np.uint64, copy=False)
        if sub_matrix.ndim == 1:
            sub_matrix = sub_matrix.reshape(1, sub_matrix.shape[0])

        if not rows_only_flag:
            allele_idx = np.asarray(allele_indices, dtype=np.uint64)
            word_idx = (allele_idx >> np.uint64(6)).astype(np.intp, copy=False)
            bit_offsets = (allele_idx & np.uint64(63)).astype(np.uint64, copy=False)
            selected_words = np.take(sub_matrix, word_idx, axis=1)
            shifted = np.right_shift(selected_words, bit_offsets[np.newaxis, :])
            sub_matrix = (shifted & np.uint64(1)).astype(np.uint8, copy=False)

        ## base headers/meta container
        headers_meta: dict = {"headers": {"guids": guids, "alleles": alleles}}

        ## propagate missing masks for subset
        if self._headerUtils.has_missing_mask():
            per_guid_masks = {}
            if rows_only_flag:
                for guid in guids:
                    per_guid_masks[guid] = list(self._headerUtils.get_guid_missing_mask(guid))
            else:
                ## remap column indices to the subset
                idx_map = {orig: new for new, orig in enumerate(allele_indices)}
                for guid in guids:
                    cols = self._headerUtils.get_guid_missing_mask(guid)
                    remapped = []
                    for c in cols:
                        if c in idx_map:
                            remapped.append(idx_map[c])
                    per_guid_masks[guid] = remapped
            headers_meta["column_masks"] = per_guid_masks
        
        if rows_only_flag and data_only:
            ## materialise dense form for users that requested data only
            dense_full = self._matrix.getBitMatrix()
            sub_matrix = dense_full[guid_indices, :].astype(np.uint8, copy=False)

        ## construct metadata for child instance
        headers_meta["meta"] = make_meta(sub_matrix,
                                         headers_meta["headers"],
                                         generated_by="ardal::subset",
                                         format_name="ardal.bitpack.v1" if (sub_matrix.dtype == np.uint64 and sub_matrix.ndim == 2 and sub_matrix.shape[1] == ((len(alleles) + 63) // 64)) else "ardal.dense.v1",
                                         matrix_file=None)

        if not data_only:
            ## return an ardal object initialised with the subset data
            from ..Ardal import Ardal
            return Ardal(data_source=[sub_matrix, headers_meta],
                         allele_id_format=self._headerUtils._allele_id_format,
                         allele_positions=self._headerUtils.get_allele_positions(),
                         verbosity=child_verbosity,
                         quiet_init=child_quiet_init)
        else:
            ## return the new subset matrix/JSON pair
            return [sub_matrix, headers_meta]
        
    
    def matrix_stats( self,
                      print_table : bool = False ) -> dict:
        """ Return a dictionary containing information about the database and its size in memory.
        """
        naturalsize = require_package("humanize", attr="naturalsize")
        
        n_guids = len(self._headerUtils.headers["guids"])
        n_alleles = len(self._headerUtils.headers["alleles"])
        density = self.density()
        roaring = self.roaring
        bit_matrix_size_bytes = self.packed_matrix().nbytes
        if self.roaring:
            roaring_mat = self.roaring_matrix(decode=False)
            roaring_size_bytes = sum(arr.nbytes for arr in roaring_mat)
            total_size_bytes = bit_matrix_size_bytes + roaring_size_bytes
            roaring_matrix_size = naturalsize(roaring_size_bytes, binary=True)
        else: 
            roaring_matrix_size = None
            total_size_bytes = bit_matrix_size_bytes
        bit_matrix_size = naturalsize(bit_matrix_size_bytes, binary=True)
        total_size = naturalsize(total_size_bytes, binary=True)
        allele_id_format = self._headerUtils._allele_id_format
        
        if self._headerUtils.allele_positions:
            sequences = len(self._headerUtils.allele_positions.keys())
            sites = sum([len(s) for chr, s in self._headerUtils.allele_positions.items()])
        else:
            sequences = None
            sites = None

        stats = {
            "Number of GUIDs"     : n_guids,
            "Number of Alleles"   : n_alleles,
            "Allele id format"    : allele_id_format,
            "Number of Sequences" : sequences,
            "Number of Sites"     : sites,
            "Matrix Density"      : density,
            "Roaring Enabled"     : roaring,
            "Bit Matrix Size"     : bit_matrix_size,
            "Roaring Matrix Size" : roaring_matrix_size,
            "Total Matrix Size"   : total_size
        }
        
        if print_table:
            ## pretty print the stats
            max_key_len = max(len(k) for k in stats.keys())
            print("\n--- Ardal Matrix Statistics ---")
            for k, v in stats.items():
                print(f"{k.ljust(max_key_len)} : {v}")
            print("-----------------------------\n")
        
        return stats
    

    def density( self ) -> float:
        """ Computes the sparsity of the allele matrix.
        """
        return self._matrix.getDensity()
    

    def bit_matrix( self ) -> np.ndarray:
        """ Return the bit allele matrix.
        """
        return self._matrix.getBitMatrix()
    

    def packed_matrix( self ):
        """ Return the 64-bit packed matrix.
        """
        arr = self._matrix.getPackedMatrix()
        
        ## ensure little endian
        if arr.dtype.byteorder not in ('<', '='):
            arr = arr.byteswap().newbyteorder('<')
            
        return arr
        
    
    def roaring_matrix( self,
                       decode : bool=True
                       ) -> dict:
        """ Return the roaring allele matrix.
        """
        roaring_dict = defaultdict(list)
        guids = self._headerUtils.headers["guids"]

        if self.roaring:
            rormat = self._matrix.getRoaringMatrix()

            if decode:
                for i, mat in enumerate(rormat):
                    allele_ids = [self._headerUtils.decode_allele(idx) for idx in mat]
                    roaring_dict[guids[i]]=allele_ids
                return roaring_dict
            else:
                return rormat
                
        else:
            raise RoaringError("Ardal object was instantialised with 'roaring=False'. Cannot retrieve roaring matrix.")
    

    def headers( self ) -> dict:
        """ Return the allele _headers.
        """
        return self._headerUtils.headers
    
    
    def meta( self ) -> dict:
        """ Return the allele matrix metadata.
        """
        return self._headerUtils.meta
    
    
    def row_masses( self ) -> dict:
        """ Return the row masses.
        """
        return self._matrix.getRowMasses()
    
    
    def column_masses( self ) -> dict:
        """ Return the col masses.
        """
        return self._matrix.getColumnMasses()
    
    
    # def missing_sites( self ) -> dict:
    #     """ Returns a dictionary of missing sites.
    #     """
    #     missing_sites = defaultdict(list)
    #     for guid, site_idxs in self._headerUtils._guid_missing_sites.items():
    #         missing_sites[guid] = [self._headerUtils._missing_site_keys.get(str(idx)) for idx in site_idxs]
    #     return missing_sites
