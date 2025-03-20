import re
import pandas as pd
import numpy as np
from humanize import naturalsize
from sys import stdout, stderr, getsizeof
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict, namedtuple

import _ardal
from .ArdalParser import ArdalParser


class Ardal(object):
    """ Class for handling Ardal-neighbourhood objects. """


    def __init__( self, data_source : str, __ref : bool=False, file_format : str=None ):
        """ Ardal constructor
        """
        super(Ardal, self).__init__()

        self.__allele_matrix = None
        self.__headers = None
        self.__ref = __ref

        parser = ArdalParser(data_source, file_format)

        if parser.matrix is not None:
            self.__allele_matrix = _ardal.AlleleMatrix(parser.matrix)
            self.__headers = parser.headers
        else:
            ## raise an error if parsing fails to prevent unexpected behaviour down the line
            raise ValueError(f"Failed to parse data from: {data_source}") 
    

    ########## COMPUTE FUNCTIONS ##########
            
    
    def pairwise(self, guids: list = [], metric: str = "hamming", chunk_size : int=100 ) -> pd.DataFrame:
        """ Takes a set of GUIDs and calculates the pairwise distance between them, returning a distance matrix.
        Pairwise distance can be calculated using Jaccard, Euclidean, or absolute distance functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """

        ## check the specified distance function is valid
        accepted_dist_functions = ["jaccard", "hamming"]
        if metric not in accepted_dist_functions:
            raise ValueError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")
        
        ## calculate the distance matrix using _ardal
        if metric == "jaccard":
            dist_tri = self.__allele_matrix.jaccard()
        elif metric == "hamming":
            dist_tri = self.__allele_matrix.hamming()

        dist_matrix = np.array(squareform(dist_tri))
        dist_df = pd.DataFrame(dist_matrix, columns=self.__headers["guids"], index=self.__headers["guids"])
        
        return dist_df

    
    def neighbourhood( self, guid : str, n : int, simd : bool = True ) -> dict:
        """ get the SNP neighbourhood of a GUID
        """

        guid_coord = self._encodeGuid(guid)
        
        ## run using SIMD, requires AVX2 support
        if simd:
            ncoords = self.__allele_matrix.neighbourhoodSIMD(guid_coord, n)

        ## run standard
        else:
            ncoords = self.__allele_matrix.neighbourhood(guid_coord, n)

        neighbourhood = {self._decodeGuid(coord) : hdist for coord, hdist in ncoords}

        return neighbourhood


    def unique( self, guids : list ) -> set:
        """
        Finds the set of SNPs unique to a given set of GUIDs.

        A SNP is considered unique if it is present in all of the specified
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
            if guid not in self.__headers["guids"]:
                raise ValueError(f"guid '{guid}' not found in allele matrix.")
            
        ## get the intersection of SNPs for the guids
        ## core SNPs for input guids
        intersection = self.core(guids)

        ## convert the guids into a set to enable set algebra
        guids_set = set(guids)

        ## construct the set of other guids
        other_guids = set(self.__headers["guids"]) - guids_set
        
        ## access the SNPs present in the other guids
        ## this will be any snp which is not unique to the guid set
        if other_guids:
            other_guid_coords = np.array([self._encodeGuid(guid) for guid in other_guids])
            union_of_others_coords = self.__allele_matrix.gatherSNPs(other_guid_coords)
            union_of_others = {self._decodeAllele(coord) for coord in union_of_others_coords}
        else:
            union_of_others = set() ## return empty set if there are no other guids to compare against

        ## difference the set
        ## subtract the SNPs present in ANY of the other GUIDs from the SNPs present in ALL of the specified GUIDs.
        unique_snps = intersection - union_of_others

        return unique_snps
    

    def core( self, guids : list, missingness : float = 0.0, return_counts : bool = False ) -> set:
        """ Take a set of guids and return alleles common to this subset.
        """

        ## check input
        if not isinstance(guids, list) and not isinstance(missingness, set):
            raise ValueError("guids must be a list or set.")
        if len(guids) == 0:
            raise ValueError("guids set cannot be empty.")
        for guid in guids:
            if guid not in self.__headers["guids"]:
                raise ValueError(f"guid '{guid}' not found in allele matrix.")
        if missingness < 0 or missingness > 1:
            raise ValueError("missingness must be between 0 and 1.")
        
        core_alleles, accessory_alleles = self._getCoreAndAccessory(guids, missingness)
        
        ## return a dictionary containing counts on the number of guids which exhibit this allele
        if return_counts:
            return core_alleles
        
        ## return snps with counts exceeding the missingness threshold
        return {allele for allele, count in core_alleles.items()}


    def accessory( self, guids : list, missingness : float = 0.0, return_counts : bool = False ) -> set:
        """ Take a set of guids and return the accessory alleles (the symmetric set of the core alleles).
        """

        ## check input
        if not isinstance(guids, list) and not isinstance(guids, set):
            raise ValueError("guids must be a list or set.")
        if len(guids) == 0:
            raise ValueError("guids set cannot be empty.")
        for guid in guids:
            if guid not in self.__headers["guids"]:
                raise ValueError(f"guid '{guid}' not found in allele matrix.")
        # if isinstance(missingness, float) and not isinstance(missingness, int):
        #     raise ValueError("missingness must be a float or integer.")
        if missingness < 0 or missingness > 1:
            raise ValueError("missingness must be between 0 and 1.")
        
        core_alleles, accessory_alleles = self._getCoreAndAccessory(guids, missingness)
        
        ## return a dictionary containing counts on the number of guids which exhibit this allele
        if return_counts:
            return accessory_alleles
        
        ## return snps with counts exceeding the missingness threshold
        return {allele for allele, count in accessory_alleles.items()}
            

    def allele( self, alleles : list ) -> set:
        """ Take a set of alleles and return all GUIDs containing all of those alleles.
        """

        ## check input
        if not isinstance(alleles, list) and not isinstance(alleles, set):
            raise ValueError("alleles must be a list or set.")
        if len(alleles) == 0:
            raise ValueError("guids set cannot be empty.")
        for allele in alleles:
            if allele not in self.__headers["alleles"]:
                raise ValueError(f"allele '{allele}' not found in allele matrix.")
            
        ## get the set of all guids which contain all of the specified alleles
        n = len(alleles)
        allele_coords = [self._encodeAllele(allele) for allele in alleles]
        input_coords = np.array([[self._encodeGuid(guid), allele] for guid in self.__headers["guids"] for allele in allele_coords])
        access_result = zip(input_coords, self.__allele_matrix.access(input_coords))
        decoded_results = [self._decodeGuid(guid_c) for (guid_c, allele_c), r in access_result if r==1]

        return {guid for guid in set(decoded_results) if decoded_results.count(guid) == n}
    

    def _getCoreAndAccessory( self, guids : list, missingness : float = 0.0 ) -> tuple:
        """ Take a set of guids and return both the core and accessory allele sets.
        """
        
        snp_count_threshold = (1-missingness) * len(guids)

        ## preprocess the guid list
        encoded_guids = np.array([self._encodeGuid(guid) for guid in guids if guid in self.__headers["guids"]])

        snp_vector = self.__allele_matrix.gatherSNPs(encoded_guids)
        
        ## count allele occurrances and return the set of SNPs which exceed the missingness threshold
        allele_count_dict = defaultdict(int)
        for i in snp_vector:
            allele_count_dict[self._decodeAllele(i)] += 1
        

        ## initialise core and accessory dictionaries
        core_alleles = defaultdict(int)
        accessory_alleles = defaultdict(int)

        ## populate core and accessoty dicts
        for alleles, count in allele_count_dict.items():
            if count >= snp_count_threshold:
                core_alleles[alleles] = count
            elif count < snp_count_threshold:
                accessory_alleles[alleles] = count
        
        return core_alleles, accessory_alleles
        

    ########## PUBLIC UTILITY FUNCTIONS ##########

    def matchAlleleNames(self, expression: str) -> list:
        """ Return all allele names that match the given expression with wildcards.
        """
        if not isinstance(expression, str):
            raise ValueError("Expression must be a string.")

        pattern = re.compile(expression.replace('*', '.*'))
        return set([allele for allele in self.__headers["alleles"] if pattern.match(allele)])
    

    def subsetbyGUID( self, guid_list : list ):
        """ Take a list of GUIDs and subset the allele matrix to include only these GUIDs, allowing for standard operations.
        Returns an Ardal object with the subsetted matrix.
        """
        return None


    def toDict( self ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        return None
    

    def stats( self ) -> dict:
        """ Return a dictionary containing information about the database and its size in memory.
        """
        stats = {
            "n_guids"     : len(self.__headers["guids"]),
            "n_alleles"   : len(self.__headers["alleles"]),
            "matrix_size" : naturalsize(self.getMatrix().nbytes, binary=True)
        }

        return stats
    

    def getMatrix( self ) -> np.array:
        """ Return the allele matrix.
        """
        return self.__allele_matrix.getMatrix()
    

    def getHeaders( self ) -> dict:
        """ Return the allele __headers.
        """
        return self.__headers
    

    def snpCount( self ) -> dict:
        """ Return a dictionary of SNP counts for each GUID.
        """
        guid_mass_vec = self.__allele_matrix.getMass()
        return {guid : mass for guid, mass in zip(self.__headers["guids"], guid_mass_vec)}
    

    def toDataFrame(self) -> pd.DataFrame:
        """ Return the allele matrix as a Pandas DataFrame.
        """
        return pd.DataFrame(self.getMatrix(), index=self.__headers["guids"], columns=self.__headers["alleles"])
    

    def flushCache(self) -> None:
        """ flushes the distance cache.
        """
        self.__allele_matrix.flushCache()
    

    ########## PRIVATE UTILITY FUNCTIONS ##########

    def _checkGUIDs( self, guids : list ) -> list:
        """ Check guids are present within the matrix and construct a list of present guids to proceed with
        """

        present_guids = []
        for id in guids:
            if id in self.__headers["guids"]:
                present_guids.append(id)
            else:
                print(f"{id} not present in allele matrix.")

        return present_guids
    

    def _checkAlleles( self, alleles : list ) -> list:
        """ Check alleles are present within the matrix and construct a list of present alleles to proceed with
        """
        present_alleles = []
        for allele_id in alleles:
            if allele_id in self.__headers["alleles"]:
                present_alleles.append(allele_id)
            else:
                print(f"{allele_id} not present in allele matrix.")

        return present_alleles
    

    def _encodeGuid( self, guid : str ):
        return self.__headers["guids"].index(guid)

    def _decodeGuid( self, row_coord : int ):
        return self.__headers["guids"][row_coord]

    def _encodeAllele( self, allele : str ):
        return self.__headers["alleles"].index(allele)

    def _decodeAllele( self, col_coord : int ):
        return self.__headers["alleles"][col_coord]