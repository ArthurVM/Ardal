""" ArdalAllele.py
This module provides functionality for working with alleles in the Ardal framework.
"""
from collections import defaultdict
from typing import List, Union, Dict

from ..utils.misc import require_package
from ..utils.decorators import check_backend_argument, check_guids_list, check_alleles_list, check_allele_id_format, check_bed_paths
from ..utils.exceptions import EmptySelectionError, RegexError
from ..utils.validators import validate_type
from ..utils.logger import get_logger

log = get_logger()


## core/ArdalAllele.py
class ArdalAllele:

    def __init__( self,
                  headerUtils,
                  hybrid_matrix,
                  roaring_enabled : bool ):

        self._headerUtils = headerUtils
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled


    @check_backend_argument
    @check_guids_list
    def unique( self,
                guids : List,
                backend : str = "auto" ) -> Dict:
        """
        Finds the set of alleles unique to a given set of GUIDs.

        An allele is considered unique if it is present in ANY of the specified
        GUIDs and absent in all other GUIDs.

        INPUT:
            guids (list): A list of GUIDs.

        OUTPUT:
            dict: A dictionary of {guid : unique_alleles} pairs.

        EXCEPTIONS:
            ValueError: If guids is not a list or set, if guids is empty, or if any GUID is not found.
        """
        ## input checks
        if guids:
            self._headerUtils.check_guids(guids)
        else:
            raise EmptySelectionError("guids list cannot be empty.")
                    
        unique_alleles = defaultdict(set)
        for guid in guids:
            guid_unique_alleles = self._matrix.uniqueSharedBits([self._headerUtils.encode_guid(guid)], 
                                                              backend=backend)
            unique_alleles[guid] = {self._headerUtils.decode_allele(idx) for idx in guid_unique_alleles}
        
        return unique_alleles
    

    @check_backend_argument
    @check_alleles_list
    def guids_with_alleles( self,
                            alleles : List,
                            backend : str = "auto" ) -> List:
        """ Take a set of alleles and return all GUIDs containing all of those alleles.
        """

        ## check input
        if not alleles:
            raise EmptySelectionError("allele list cannot be empty.")
        else:
            self._headerUtils.check_alleles(alleles)
                    
        ## get the set of all guids which contain all of the specified alleles
        allele_indices = {self._headerUtils.encode_allele(allele) for allele in alleles}

        result_guids = set()
        for guid_idx, guid_name in enumerate(self._headerUtils.headers["guids"]):
            present_alleles = set(self._matrix.getSetBitIndices(guid_idx, backend=backend))
            if allele_indices.issubset(present_alleles):
                result_guids.add(guid_name)
        
        return list(result_guids)


    def match_names( self,
                     expression: str ) -> Union[List, None]:
        """ Return all allele names that match the given expression with wildcards.
        """
        re = require_package("re", "re")
        
        validate_type(expression, str, "expression")
        try:
            pattern = re.compile(expression.replace('*', '.*'))
            return list(set([allele for allele in self._headerUtils.headers["alleles"] if pattern.match(allele)]))
        except Exception as e:
            RegexError(f"Regex using expression {expression} failed.")
            

    @check_backend_argument
    @check_guids_list
    def unique_core( self,
                     guids : List,
                     backend : str = "auto" ) -> list:
        """
        Finds the set of alleles unique to a given set of GUIDs and shared by all GUIDs in that set.

        An allele is considered unique core if it is present in ALL of the specified
        GUIDs and absent in all other GUIDs.

        INPUT:
            guids (list): A list of GUIDs.

        OUTPUT:
            list: A list of unique alleles.

        EXCEPTIONS:
            ValueError: If guids is not a list or set, if guids is empty, or if any GUID is not found.
        """
        ## input checks
        self._headerUtils.check_guids(guids)
        
        ## if all guids are selected, just return all alleles
        if guids == self._headerUtils.headers["guids"]:
            return self._headerUtils.headers["alleles"]
            
        guid_coords = [self._headerUtils.encode_guid(guid) for guid in guids]
        unique_allele_indices = self._matrix.uniqueSharedBits(guid_coords, backend=backend)
        results = list({self._headerUtils.decode_allele(idx) for idx in unique_allele_indices})
        
        if len(results) == 0:
            log.warning(f"No unique core alleles found. This is common in large databases. Consider using the stats module if you are looking for alleles which are over-represented in a set of guids.")
            
        return results
    
    
    def count( self ) -> dict:
        """ Return a dictionary of allele counts for each GUID.
        """
        guid_mass_vec = self._matrix.getRowMasses()
        return {guid : mass for guid, mass in zip(self._headerUtils.headers["guids"], guid_mass_vec)}
    
    
    @check_allele_id_format
    @check_bed_paths
    def interval( self,
                  intervals : list,
                  intervals_bed : Union[str, None] = None,
                  allele_coords_bed : Union[str, None] = None,
                  allele_id_format : str = "{chr}.{start}.{ref}.{alt}"
                  ) -> List:
        """ Return a list of alleles which fall within the given genomic intervals.
        """
        return self._headerUtils.get_interval_alleles(intervals=intervals,
                                                      allele_id_format=allele_id_format,
                                                      intervals_bed=intervals_bed,
                                                      allele_coords_bed=allele_coords_bed)
        
        

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
        
    #     ## return alleles with counts exceeding the missingness threshold
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
        
    #     ## return alleles with counts exceeding the missingness threshold
    #     return {allele for allele, count in accessory_alleles.items()}
    

    # def _getCoreAndAccessory( self, guids : list, missingness : float = 0.0 ) -> tuple:
    #     """ Take a set of guids and return both the core and accessory allele sets.
    #     """
    #     allele_count_threshold = (1-missingness) * len(guids)

    #     ## preprocess the guid list
    #     encoded_guids = np.array([self._encode_guid(guid) for guid in guids if guid in self._headers["guids"]])

    #     allele_count_dict = defaultdict(int)
    #     for guid_coord in encoded_guids:
    #         allele_indices = self._matrix.getSetBitIndices(guid_coord)
    #         for allele_idx in allele_indices:
    #             allele_count_dict[self._decode_allele(allele_idx)] += 1
        
    #     ## initialise core and accessory dictionaries
    #     core_alleles = defaultdict(int)
    #     accessory_alleles = defaultdict(int)

    #     ## populate core and accessoty dicts
    #     for alleles, count in allele_count_dict.items():
    #         if count >= allele_count_threshold:
    #             core_alleles[alleles] = count
    #         elif count < allele_count_threshold:
    #             accessory_alleles[alleles] = count
        
    #     return core_alleles, accessory_alleles