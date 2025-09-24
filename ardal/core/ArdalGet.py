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

    
    @check_thread_count
    @check_guids_list
    @check_alleles_list
    def subset( self,
                guids : List[str] = [],
                alleles : List[str] = [],
                data_only : bool = False,
                threads : int = 1,
                child_verbosity : str = "silent",
                child_quiet_init : bool = True
                ) -> Union[Tuple[np.ndarray, dict], "Ardal"]:
        """ Take a list of GUIDs and subset the allele matrix to include only these GUIDs, allowing for standard operations.
        Returns a numpy matrix/JSON pair for feeding into Ardal.
        """
        from ..Ardal import Ardal

        ## check input
        if len(guids) == 0 and len(alleles) == 0:
            raise ParameterError("guids and alleles cannot both be empty.")

        ## check GUIDs
        if guids:
            self._headerUtils.check_guids(guids)
        else:
            guids = self._headerUtils.headers["guids"]

        ## check alleles
        if alleles:
            self._headerUtils.check_alleles(alleles)
        else:
            alleles = self._headerUtils.headers["alleles"]
            
        log.info(f"Subsetting the {self._headerUtils.n_guids}x{self._headerUtils.n_alleles} matrix to {len(guids)}x{len(alleles)}")
        
        guid_indices = sorted([self._headerUtils.encode_guid(guid) for guid in guids], reverse=False)
        allele_indices = sorted([self._headerUtils.encode_allele(allele) for allele in alleles], reverse=False)
        
        ## subset the DataFrame
        ## TODO: this is pretty grim, could be done better in C++
        # subset_df = pd.DataFrame(self._matrix.getBitMatrix(), index=self._headerUtils.headers["guids"], columns=self._headerUtils.headers["alleles"]).loc[guids, alleles]
        # sub_matrix = subset_df.values.astype(np.uint8)
        
        sub_matrix = self._matrix.getSubsetPackedMatrix(guid_indices,
                                                        allele_indices,
                                                        threads)

        ## create new headers and matrix for the new Ardal object
        sub_headers = {"guids": guids, "alleles": alleles}
        
        if not data_only:
            ## return an ardal object initialised with the subset data
            return Ardal(data_source=[sub_matrix, sub_headers],
                         allele_positions=self._headerUtils.get_allele_positions(),
                         verbosity=child_verbosity,
                         quiet_init=child_quiet_init)
        else:
            ## return the new subset matrix/JSON pair
            return [sub_matrix, sub_headers]
        
    
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

        stats = {
            "Number of GUIDs"     : n_guids,
            "Number of Alleles"   : n_alleles,
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