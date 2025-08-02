""" ArdalGet.py
This module provides functionality for retrieving and manipulating allele matrices in the Ardal framework.
"""
import pandas as pd
import numpy as np
from collections import defaultdict
from humanize import naturalsize

from .utilities import *
from ..exceptions.exceptions import *


# core/ArdalGet.py
class ArdalGet:
    """ Class for retrieving and manipulating allele matrices in the Ardal framework.
    """

    def __init__(self, headers, hybrid_matrix, roaring_enabled):

        self._headers = headers
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled

    
    def subset( self,
                guid_list : list = [],
                allele_list : list = [] ):
        """ Take a list of GUIDs and subset the allele matrix to include only these GUIDs, allowing for standard operations.
        Returns a numpy matrix/JSON pair for feeding into Ardal.
        """
        from ..Ardal import Ardal

        ## check input
        if len(guid_list) == 0 and len(allele_list) == 0:
            raise ParameterError("guid_list and allele_list cannot both be empty.")

        ## check GUIDs
        if guid_list:
            checkGUIDs(guid_list, self._headers)
        else:
            guid_list = self._headers["guids"]

        ## check alleles
        if allele_list:
            checkAlleles(allele_list, self._headers)
        else:
            allele_list = self._headers["alleles"]

        ## subset the DataFrame
        ## TODO: this is pretty grim, could be done better in C++
        subset_df = pd.DataFrame(self._matrix.getBitMatrix(), index=self._headers["guids"], columns=self._headers["alleles"]).loc[guid_list, allele_list]

        ## create new headers and matrix for the new Ardal object
        new_headers = {"guids": subset_df.index.tolist(), "alleles": subset_df.columns.tolist()}
        new_matrix = subset_df.values.astype(np.uint8)
        
        ## return the new subset matrix/JSON pair
        return Ardal([new_matrix, new_headers], quiet=True)
        
    
    def matrixStats( self,
                     print_stats : bool = False ) -> dict:
        """ Return a dictionary containing information about the database and its size in memory.
        """
        n_guids = len(self._headers["guids"])
        n_alleles = len(self._headers["alleles"])
        density = self.density()
        roaring = self.roaring
        bit_matrix_size_bytes = self.bitMatrix().nbytes
        if self.roaring:
            roaring_mat = self.roaringMatrix(decode=False)
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
        
        if print_stats:
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
    


    def bitMatrix( self ) -> np.array:
        """ Return the bit allele matrix.
        """
        return self._matrix.getBitMatrix()
    

    def hybridMatrix( self ):
        """ Return the hybridMatrix.
        """
        return self._matrix

    
    def roaringMatrix( self,
                       decode : bool=True ) -> dict:
        """ Return the roaring allele matrix.
        """
        roaring_dict = defaultdict(list)
        guids = self._headers["guids"]

        if self.roaring:
            rormat = self._matrix.getRoaringMatrix()

            if decode:
                for i, mat in enumerate(rormat):
                    allele_ids = [decodeAllele(idx, self._headers) for idx in mat]
                    roaring_dict[guids[i]]=allele_ids
                return roaring_dict
            else:
                return rormat
                
        else:
            raise RoaringError("Ardal object was instantialised with 'use_roaring_if_sparse=False'. Cannot retrieve roaring matrix.")
    

    def headers( self ) -> dict:
        """ Return the allele _headers.
        """
        return self._headers
    

    def snpCount( self ) -> dict:
        """ Return a dictionary of SNP counts for each GUID.
        """
        guid_mass_vec = self._matrix.getRowMasses()
        return {guid : mass for guid, mass in zip(self._headers["guids"], guid_mass_vec)}
    
    
    def intervalAlleles( self,
                         intervals : list,
                         coords_bed : str = None,
                         allele_id_format : str = "{chr}.{start}.{ref}.{alt}" ) -> list:
        """ Return a list of alleles which fall within the given genomic intervals.
        """
        return getIntervalAlleles(intervals, self._headers, allele_id_format, coords_bed)
    