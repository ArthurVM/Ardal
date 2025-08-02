""" ArdalAllele.py
This module provides functionality for working with alleles in the Ardal framework.
"""
from collections import defaultdict

from .utilities import *
from ..exceptions.exceptions import *


# core/ArdalAllele.py
class ArdalAllele:

    def __init__(self, headers, hybrid_matrix, roaring_enabled):

        self._headers = headers
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled


    def unique( self,
                guid_list : list,
                force_bit_backend: bool=False ) -> set:
        """
        Finds the set of SNPs unique to a given set of GUIDs.

        A SNP is considered unique if it is present in ANY of the specified
        GUIDs and absent in all other GUIDs.

        INPUT:
            guids (list): A list of GUIDs.

        OUTPUT:
            dict: A dictionary of {guid : unique_snps} pairs.

        EXCEPTIONS:
            ValueError: If guids is not a list or set, if guids is empty, or if any GUID is not found.
        """
        ## input checks
        if guid_list:
            checkGUIDs(guid_list, self._headers)
        else:
            raise EmptySelectionError("guids list cannot be empty.")
                    
        unqiue_snps = defaultdict(set)
        for guid in guid_list:
            guid_unique_snps = self._matrix.uniqueSharedBits([encodeGuid(guid, self._headers)], 
                                                             force_bit_backend=force_bit_backend)
            unqiue_snps[guid] = {decodeAllele(idx, self._headers) for idx in guid_unique_snps}
        
        return unqiue_snps
    

    def guidsWithAlleles( self,
                          allele_list : list,
                          force_bit_backend : bool=False ) -> set:
        """ Take a set of alleles and return all GUIDs containing all of those alleles.
        """

        ## check input
        if not allele_list:
            raise EmptySelectionError("allele list cannot be empty.")
        else:
            checkAlleles(allele_list, self._headers)
                    
        ## get the set of all guids which contain all of the specified alleles
        allele_indices = {encodeAllele(allele, self._headers) for allele in allele_list}

        result_guids = set()
        for guid_idx, guid_name in enumerate(self._headers["guids"]):
            present_alleles = set(self._matrix.getSetBitIndices(guid_idx, force_bit_backend=force_bit_backend))
            if allele_indices.issubset(present_alleles):
                result_guids.add(guid_name)
        
        return result_guids


    def matchNames( self,
                    expression: str ) -> list:
        """ Return all allele names that match the given expression with wildcards.
        """
        re = require_package("re", "re")

        if not isinstance(expression, str):
            raise TypeError("Expression must be a string.")

        pattern = re.compile(expression.replace('*', '.*'))
        return set([allele for allele in self._headers["alleles"] if pattern.match(allele)])
    

    def uniqueCore( self, guids : list, force_bit_backend: bool=False ) -> set:
        """
        Finds the set of SNPs unique to a given set of GUIDs and shared by all GUIDs in that set.

        A SNP is considered unique core if it is present in ALL of the specified
        GUIDs and absent in all other GUIDs.

        INPUT:
            guids (list): A list of GUIDs.

        OUTPUT:
            set: A set of unique SNPs.

        EXCEPTIONS:
            ValueError: If guids is not a list or set, if guids is empty, or if any GUID is not found.
        """

        ## input checks
        if not isinstance(guids, list):
            raise ValueError("guids must be a list.")
        if len(guids) == 0:
            raise ValueError("guids set cannot be empty.")
        for guid in guids:
            if guid not in self._headers["guids"]:
                raise ValueError(f"guid '{guid}' not found in allele matrix.")
            
        guid_coords = [encodeGuid(guid, self._headers) for guid in guids]
        unique_snp_indices = self._matrix.uniqueSharedBits(guid_coords, force_bit_backend=force_bit_backend)
        return {decodeAllele(idx, self._headers) for idx in unique_snp_indices}
    

    # def core( self, guids : list, missingness : float = 0.0, return_counts : bool = False ) -> set:
    #     """ Take a set of guids and return alleles common to this subset.
    #     """

    #     ## check input
    #     if not isinstance(guids, list) and not isinstance(missingness, set):
    #         raise ValueError("guids must be a list or set.")
    #     if len(guids) == 0:
    #         raise ValueError("guids set cannot be empty.")
    #     for guid in guids:
    #         if guid not in self._headers["guids"]:
    #             raise ValueError(f"guid '{guid}' not found in allele matrix.")
    #     if missingness < 0 or missingness > 1:
    #         raise ValueError("missingness must be between 0 and 1.")
        
    #     core_alleles, accessory_alleles = self._getCoreAndAccessory(guids, missingness)
        
    #     ## return a dictionary containing counts on the number of guids which exhibit this allele
    #     if return_counts:
    #         return core_alleles
        
    #     ## return snps with counts exceeding the missingness threshold
    #     return {allele for allele, count in core_alleles.items()}


    # def accessory( self, guids : list, missingness : float = 0.0, return_counts : bool = False ) -> set:
    #     """ Take a set of guids and return the accessory alleles (the symmetric set of the core alleles).
    #     """

    #     ## check input
    #     if not isinstance(guids, list) and not isinstance(guids, set):
    #         raise ValueError("guids must be a list or set.")
    #     if len(guids) == 0:
    #         raise ValueError("guids set cannot be empty.")
    #     for guid in guids:
    #         if guid not in self._headers["guids"]:
    #             raise ValueError(f"guid '{guid}' not found in allele matrix.")
    #     # if isinstance(missingness, float) and not isinstance(missingness, int):
    #     #     raise ValueError("missingness must be a float or integer.")
    #     if missingness < 0 or missingness > 1:
    #         raise ValueError("missingness must be between 0 and 1.")
        
    #     core_alleles, accessory_alleles = self._getCoreAndAccessory(guids, missingness)
        
    #     ## return a dictionary containing counts on the number of guids which exhibit this allele
    #     if return_counts:
    #         return accessory_alleles
        
    #     ## return snps with counts exceeding the missingness threshold
    #     return {allele for allele, count in accessory_alleles.items()}
    

    # def _getCoreAndAccessory( self, guids : list, missingness : float = 0.0 ) -> tuple:
    #     """ Take a set of guids and return both the core and accessory allele sets.
    #     """
    #     snp_count_threshold = (1-missingness) * len(guids)

    #     ## preprocess the guid list
    #     encoded_guids = np.array([self._encodeGuid(guid) for guid in guids if guid in self._headers["guids"]])

    #     allele_count_dict = defaultdict(int)
    #     for guid_coord in encoded_guids:
    #         snp_indices = self._matrix.getSetBitIndices(guid_coord)
    #         for snp_idx in snp_indices:
    #             allele_count_dict[self._decodeAllele(snp_idx)] += 1
        
    #     ## initialise core and accessory dictionaries
    #     core_alleles = defaultdict(int)
    #     accessory_alleles = defaultdict(int)

    #     ## populate core and accessoty dicts
    #     for alleles, count in allele_count_dict.items():
    #         if count >= snp_count_threshold:
    #             core_alleles[alleles] = count
    #         elif count < snp_count_threshold:
    #             accessory_alleles[alleles] = count
        
    #     return core_alleles, accessory_alleles