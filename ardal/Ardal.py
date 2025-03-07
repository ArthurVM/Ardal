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


    def __init__( self, data_source : str, _ref : bool=False, file_format : str=None ):
        """ Ardal constructor
        """
        super(Ardal, self).__init__()

        self._allele_matrix = None
        self._headers = None
        self._ref = _ref

        parser = ArdalParser(data_source, file_format)

        if parser.matrix is not None:
            self._allele_matrix = _ardal.AlleleMatrix(parser.matrix)
            self._headers = parser.headers
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
            dist_tri = self._allele_matrix.jaccard()
        elif metric == "hamming":
            dist_tri = self._allele_matrix.hamming()

        dist_matrix = np.array(squareform(dist_tri))
        dist_df = pd.DataFrame(dist_matrix, columns=self._headers["guids"], index=self._headers["guids"])
        
        return dist_df

    
    def neighbourhood( self, guid : str, n : int, simd : bool = True ) -> dict:
        """ get the SNP neighbourhood of a GUID
        """

        guid_coord = self._encodeGuid(guid)
        
        ## run using SIMD, requires AVX2 support
        if simd:
            ncoords = self._allele_matrix.neighbourhoodSIMD(guid_coord, n)

        ## run standard
        else:
            ncoords = self._allele_matrix.neighbourhood(guid_coord, n)

        neighbourhood = {self._decodeGuid(coord) : hdist for coord, hdist in ncoords}

        return neighbourhood


    def unique( self, guids : list ) -> set:
        """ Take a set of guids and return alleles unique to this subset.
        """
        intersection, symmetric = self._getIntersectionAndSymmetric(guids)

        return intersection.difference(symmetric)
    

    def core( self, guids : list, missingness : float = 0.0, return_counts : bool = False ) -> set:
        """ Take a set of guids and return alleles common to this subset.
        """

        ## check input
        if not isinstance(guids, list) and not isinstance(missingness, set):
            raise ValueError("guids must be a list or set.")
        if len(guids) == 0:
            raise ValueError("guids set cannot be empty.")
        for guid in guids:
            if guid not in self._headers["guids"]:
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
        if not isinstance(guids, list) and not isinstance(missingness, set):
            raise ValueError("guids must be a list or set.")
        if len(guids) == 0:
            raise ValueError("guids set cannot be empty.")
        for guid in guids:
            if guid not in self._headers["guids"]:
                raise ValueError(f"guid '{guid}' not found in allele matrix.")
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

        ## check guids are present within the matrix and construct a list of present guids to proceed with
        present_alleles = self._checkAlleles(alleles)
        allele_count = len(present_alleles)

        ## initialise a dictionary for storing present alleles
        allele_dict = { guid : [] for guid in self._allele_matrix.index.values }

        for allele_id in present_alleles:
            mask = self._allele_matrix[allele_id].astype(bool)
            for guid in self._allele_matrix.loc[mask].index.values:
                allele_dict[guid].append(allele_id)

        return set([guid for guid, alleles in allele_dict.items() if len(alleles) == allele_count])
    

    def _getCoreAndAccessory( self, guids : list, missingness : float = 0.0 ) -> tuple:
        """ Take a set of guids and return both the core and accessory allele sets.
        """
        
        snp_count_threshold = (1-missingness) * len(guids)

        ## preprocess the guid list
        encoded_guids = np.array([self._encodeGuid(guid) for guid in guids if guid in self._headers["guids"]])

        snp_vector = self._allele_matrix.gatherSNPs(encoded_guids)
        
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
        return set([allele for allele in self._headers["alleles"] if pattern.match(allele)])
    

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
            "n_guids"     : len(self._headers["guids"]),
            "n_alleles"   : len(self._headers["alleles"]),
            "matrix_size" : naturalsize(self.getMatrix().nbytes, binary=True)
        }

        return stats
    

    def getMatrix( self ) -> np.array:
        """ Return the allele matrix.
        """
        return self._allele_matrix.getMatrix()
    

    def getHeaders( self ) -> dict:
        """ Return the allele _headers.
        """
        return self._headers
    

    def snpCount( self ) -> dict:
        """ Return a dictionary of SNP counts for each GUID.
        """
        guid_mass_vec = self._allele_matrix.getMass()
        return {guid : mass for guid, mass in zip(self._headers["guids"], guid_mass_vec)}
    

    def toDataFrame(self) -> pd.DataFrame:
        """ Return the allele matrix as a Pandas DataFrame.
        """
        return pd.DataFrame(self.getMatrix(), index=self._headers["guids"], columns=self._headers["alleles"])
    

    ########## PRIVATE UTILITY FUNCTIONS ##########

    def _checkGUIDs( self, guids : list ) -> list:
        """ Check guids are present within the matrix and construct a list of present guids to proceed with
        """

        present_guids = []
        for id in guids:
            if id in self._headers["guids"]:
                present_guids.append(id)
            else:
                print(f"{id} not present in allele matrix.")

        return present_guids
    

    def _checkAlleles( self, alleles : list ) -> list:
        """ Check alleles are present within the matrix and construct a list of present alleles to proceed with
        """
        present_alleles = []
        for allele_id in alleles:
            if allele_id in self._headers["alleles"]:
                present_alleles.append(allele_id)
            else:
                print(f"{allele_id} not present in allele matrix.")

        return present_alleles
    

    def _encodeGuid( self, guid : str ):
        return self._headers["guids"].index(guid)

    def _decodeGuid( self, row_coord : int ):
        return self._headers["guids"][row_coord]

    def _encodeAllele( self, allele : str ):
        return self._headers["alleles"].index(allele)

    def _decodeAllele( self, col_coord : int ):
        return self._headers["alleles"][col_coord]