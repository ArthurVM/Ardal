import re
import os
import json
import pandas as pd
import numpy as np
from humanize import naturalsize
from sys import stdout, stderr, getsizeof
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict, namedtuple
from typing import List, Union

from . import _ardal
from .ArdalParser import ArdalParser


class Ardal(object):
    """ Class for handling Ardal-neighbourhood objects. """


    def __init__( self,
                  data_source : str,
                  use_roaring_if_sparse : bool=True,
                  density_threshold : float=0.02,
                  force_roaring : bool=False,
                  __ref : bool=False, 
                  file_format : str=None ):
        """ Ardal constructor
        """
        super(Ardal, self).__init__()

        self.__hybrid_matrix = None
        self.__headers = None
        self.__ref = __ref

        parser = ArdalParser(data_source, file_format)

        ## do some parameter fiddling to enforce the generation of a roaring matrix
        if force_roaring:
            use_roaring_if_sparse = True
            density_threshold = 1.0

        if parser.matrix is not None:
            self.__hybrid_matrix = _ardal.HybridMatrix(parser.matrix,
                                                       use_roaring_if_sparse=use_roaring_if_sparse,
                                                       density_threshold=density_threshold)
            self.__headers = parser.headers
        else:
            ## raise an error if parsing fails to prevent unexpected behaviour down the line
            raise ValueError(f"Failed to parse data from: {data_source}") 
        
        self.roaring = self.__hybrid_matrix.roaringEnabled()
        self.stats(print_stats=True)




    ########## COMPUTE FUNCTIONS ##########
            
    
    def pairwise( self,
                  guids: list = [],
                  metric: str = "hamming",
                  use_simd : bool=True,
                  threads : int=1,
                  force_bit_backend: bool=False,
                  chunk_size : int=100 ) -> pd.DataFrame:
        """ Takes a set of GUIDs and calculates the pairwise distance between them, returning a distance matrix.
        Pairwise distance can be calculated using Hamming, Jaccard, or Inner Product (number of shared SNPs) functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """
        ## check there is enough memory to store the pairwise matrix
        mat_size = len(self.__headers["guids"])**2
        total_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        if mat_size * 8 > total_memory * 0.8: ## 8 bytes per float64, 80% of total memory
            raise MemoryError(f"Pairwise distance matrix of scale {len(self.__headers['guids'])}x{len(self.__headers['guids'])} will requre {naturalsize(mat_size * 8, binary=True)} memory and will exceed system memory limits. Please subset your data.")

        ## check the specified distance function is valid
        accepted_dist_functions = ["jaccard", "hamming", "inner_product"]
        if metric not in accepted_dist_functions:
            raise ValueError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")
        
        ## calculate the distance matrix using _ardal
        if metric == "jaccard":
            # dist_tri = self.__hybrid_matrix.jaccard()
            raise NotImplementedError("Jaccard distance is not yet implemented.")
        elif metric == "hamming":
            dist_tri = self.__hybrid_matrix.hamming(use_simd=use_simd,
                                                    threads=threads,
                                                    force_bit_backend=force_bit_backend)
        elif metric == "inner_product":
            dist_tri = self.__hybrid_matrix.innerProduct(use_simd=use_simd,
                                                         threads=threads,
                                                         force_bit_backend=force_bit_backend)
        
        dist_matrix = np.array(squareform(dist_tri), dtype=np.uint32)
        dist_df = pd.DataFrame(dist_matrix, columns=self.__headers["guids"], index=self.__headers["guids"])
        
        return dist_df

    
    def neighbourhood( self,
                       guid : str,
                       n : int,
                       use_simd : bool = True,
                       threads : int = 1,
                       force_bit_backend: bool=False ) -> dict:
        """ get the allele neighbourhood of a GUID using hamming distance
        """

        if not isinstance(guid, str):
            raise ValueError("guid must be a string.")
        if not isinstance(n, int):
            raise ValueError("n must be an integer.")
        if n < 0:
            raise ValueError("n must be non-negative.")
        if n == 0:
            return {}
        if guid not in self.__headers["guids"]:
            raise ValueError(f"guid '{guid}' not found in allele matrix.")


        guid_coord = self._encodeGuid(guid)
        ncoords = self.__hybrid_matrix.neighbourhood(row_idx=guid_coord,
                                                        epsilon=n,
                                                        use_simd=use_simd,
                                                        threads=threads,
                                                        force_bit_backend=force_bit_backend)
        neighbourhood = {self._decodeGuid(coord) : hdist for coord, hdist in ncoords}

        return neighbourhood
    

    def snvDist( self,
                 guid : str,
                 n : int ) -> dict:
        """ find all GUIDs which lie within n SNVs of a given GUID
        WARNING : NOT PRODUCTION READY
        assumes allele ID of form {ref_nucleotide}{pos}{alt_nucleotide} and so the pos can be indexed out with [1:-1]
        """
        snv_neighbourhood = {}

        guid_coord = self._encodeGuid(guid)
        guid_snv_positions = set([self._decodeAllele(allele_coord)[1:-1] for allele_coord in self.__hybrid_matrix.getSetBitIndices(guid_coord)])

        for guid_idx in range(len(self.__headers["guids"])):
            if guid_idx == guid_coord:
                continue
            other_snv_positions = set([self._decodeAllele(allele_coord)[1:-1] for allele_coord in self.__hybrid_matrix.getSetBitIndices(guid_idx)])
            snv_dist = guid_snv_positions ^ other_snv_positions
            if len(snv_dist) <= n:
                snv_neighbourhood[self._decodeGuid(guid_idx)] = len(snv_dist)
        
        return snv_neighbourhood


    def uniqueSNPs( self,
                    guids : list,
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
        if not isinstance(guids, list):
            raise ValueError("guids must be a list.")
        if len(guids) == 0:
            raise ValueError("guids set cannot be empty.")
        for guid in guids:
            if guid not in self.__headers["guids"]:
                raise ValueError(f"guid '{guid}' not found in allele matrix.")
            
        unqiue_snps = defaultdict(set)
        for guid in guids:
            guid_unique_snps = self.__hybrid_matrix.uniqueSharedBits([self._encodeGuid(guid)], 
                                                                     force_bit_backend=force_bit_backend)
            unqiue_snps[guid] = {self._decodeAllele(idx) for idx in guid_unique_snps}
        
        return unqiue_snps
    

    def hasSNPs( self,
                 alleles : list,
                 force_bit_backend : bool=False ) -> set:
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
        allele_indices = {self._encodeAllele(allele) for allele in alleles}

        result_guids = set()
        for guid_idx, guid_name in enumerate(self.__headers["guids"]):
            present_alleles = set(self.__hybrid_matrix.getSetBitIndices(guid_idx, force_bit_backend=force_bit_backend))
            if allele_indices.issubset(present_alleles):
                result_guids.add(guid_name)
        
        return result_guids


    # def uniqueCore( self, guids : list, force_bit_backend: bool=False ) -> set:
    #     """
    #     Finds the set of SNPs unique to a given set of GUIDs and shared by all GUIDs in that set.

    #     A SNP is considered unique core if it is present in ALL of the specified
    #     GUIDs and absent in all other GUIDs.

    #     INPUT:
    #         guids (list): A list of GUIDs.

    #     OUTPUT:
    #         set: A set of unique SNPs.

    #     EXCEPTIONS:
    #         ValueError: If guids is not a list or set, if guids is empty, or if any GUID is not found.
    #     """

    #     ## input checks
    #     if not isinstance(guids, list):
    #         raise ValueError("guids must be a list.")
    #     if len(guids) == 0:
    #         raise ValueError("guids set cannot be empty.")
    #     for guid in guids:
    #         if guid not in self.__headers["guids"]:
    #             raise ValueError(f"guid '{guid}' not found in allele matrix.")
            
    #     guid_coords = [self._encodeGuid(guid) for guid in guids]
    #     unique_snp_indices = self.__hybrid_matrix.uniqueSharedBits(guid_coords, force_bit_backend=force_bit_backend)
    #     return {self._decodeAllele(idx) for idx in unique_snp_indices}
    

    # def core( self, guids : list, missingness : float = 0.0, return_counts : bool = False ) -> set:
    #     """ Take a set of guids and return alleles common to this subset.
    #     """

    #     ## check input
    #     if not isinstance(guids, list) and not isinstance(missingness, set):
    #         raise ValueError("guids must be a list or set.")
    #     if len(guids) == 0:
    #         raise ValueError("guids set cannot be empty.")
    #     for guid in guids:
    #         if guid not in self.__headers["guids"]:
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
    #         if guid not in self.__headers["guids"]:
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
    #     encoded_guids = np.array([self._encodeGuid(guid) for guid in guids if guid in self.__headers["guids"]])

    #     allele_count_dict = defaultdict(int)
    #     for guid_coord in encoded_guids:
    #         snp_indices = self.__hybrid_matrix.getSetBitIndices(guid_coord)
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
    



    ########## STATISTICAL FUNCTIONS ##########


    def af( self,
            guids : list=[] ) -> dict:
        """ Calculates the allele frequency for each allele in the matrix across the group of guids.
        if an empty list is provided then all guids will be used.
        """
        if guids == []:
            guids = self.__headers["guids"]
        else:
            ## check input
            for guid in guids:
                if guid not in self.__headers["guids"]:
                    raise ValueError(f"guid '{guid}' not found in allele matrix.")
        
        guid_coords = [self._encodeGuid(guid) for guid in guids]
        allele_freq = self.__hybrid_matrix.colFrequency(guid_coords)

        return dict(sorted(zip(self.__headers["alleles"], allele_freq), key=lambda x: x[1], reverse=True))
        

    def alleleEntropy( self ) -> dict:
        """ Computes the Shannon entropy for each allele in the matrix.
        Entropy H(X) = -p * log2(p) - (1-p) * log2(1-p), where p is the allele frequency.
        """
        allele_entropy = self.__hybrid_matrix.columnEntropy()
        entropies_dict = dict(sorted(zip(self.__headers["alleles"], allele_entropy), key=lambda x: x[1], reverse=True))
        return entropies_dict

    
    def alleleCooccurrance( self, threshold : float = 0.95, use_simd : bool = True, cache_bytes : int = 1024**3 ) -> dict:
        """ Computes the co-occurance of alleles.
        """
        print("HERE")
        cooc_dict = self.__hybrid_matrix.bitCooccurrence(threshold=threshold,
                                                         use_simd=use_simd,
                                                         cache_bytes=cache_bytes)

        decoded_dict = defaultdict(list)
        for ref, cooc_vec in cooc_dict.items():
            decoded_dict[self._decodeAllele(ref)] = [self._decodeAllele(i) for i in cooc_vec]

        return decoded_dict


    def snpInform( self,
                   guids: list,
                   metric: str = "kullbackleibler" ) -> dict:
        """
        Calculates various scores that measure the association of each SNP
        with a specified group of samples (guids).

        This function serves as a unified interface for metrics that return a single
        float value per SNP, representing the strength or nature of the association.
        It compares an "in-group" (the provided guids) against an "out-group"
        (all other samples).

        Args:
            guids (list): A list of sample identifiers for the in-group.
            metrics (string): Specify which metric to use. Default : 'kullbackleibler'.
                Available metrics:
                - 'kullbackleibler': Kullback-Leibler divergence.
                - 'jensenshannon': Jensen-Shannon divergence.
                - 'informationgain': Information Gain.

        Returns:
            dict: A dictionary of snp : score pairs.
        """
        ## input validation
        available_metrics = {'kullbackleibler', 'jensenshannon', 'informationgain'}
        if metric not in available_metrics:
            raise ValueError(f"Metric '{metric}' is not supported. Available metrics: {list(available_metrics)}")
        
        if not isinstance(guids, list) or not guids:
            raise ValueError("guids must be a non-empty list.")
        for guid in guids:
            if guid not in self.__headers["guids"]:
                raise ValueError(f"guid '{guid}' not found in allele matrix.")

        if metric == 'kullbackleibler':
            return self._klDivergence(guids)
        if metric == 'jensenshannon':
            return self._jsDivergence(guids)
        if metric == 'informationgain':
            return self._informationGain(guids)
    
    
    def _klDivergence( self,
                       guids : list ) -> dict:
        """ Computes the Kullbeck Liebler divergence between the in group (input guids) and out group (all others)
        allele frequency distributions for each allele in the matrix.
        D_{kl}(P||Q) = sum_{x\inX}(P(x) * log2(P(x)/Q(x)))
        """
        guid_coords = [self._encodeGuid(guid) for guid in guids]
        kl_divergence = self.__hybrid_matrix.klDivergence(guid_coords)
        kl_dict = dict(sorted(zip(self.__headers["alleles"], kl_divergence), key=lambda x: x[1], reverse=True))
        return kl_dict

    
    def _jsDivergence( self,
                       guids : list ) -> dict:
        """
        Computes Jensen-Shannon divergence for each SNP between target_guids and others.
        
        Args:
            guids: list of sample identifiers to define the target group
        """
        from scipy.spatial.distance import jensenshannon

        all_guids = self.getHeaders()["guids"]
        X = self.getMatrix()
        
        idx_map = np.array([guid in guids for guid in all_guids])
        guid_coords = [self._encodeGuid(guid) for guid in guids]
        
        ## split into in-group and out-group
        X_target = X[idx_map]
        X_other  = X[~idx_map]

        ## calculate allele frequencies per SNP
        P = X_target.mean(axis=0) + 1e-9   # Add small epsilon to avoid log(0)
        Q = X_other.mean(axis=0) + 1e-9

        ## create 2 row matrix and compute JS over each column
        M = 0.5 * (P + Q)
        js_scores = 0.5 * (P * np.log2(P / M) + Q * np.log2(Q / M))

        return dict(sorted(zip(self.__headers["alleles"], js_scores), key=lambda x: x[1], reverse=True))


    def _informationGain( self,
                          guids : list ) -> dict:
        """
        Compute information gain for each SNP column in binary matrix X,
        with respect to whether a sample is in target_guids.

        Returns:
            np.ndarray of shape (n_snps,) - information gain per SNP
        """
        from scipy.stats import entropy

        all_guids = self.getHeaders()["guids"]
        X = self.getMatrix()

        y = np.array([1 if guid in guids else 0 for guid in all_guids])
        H_C = entropy(np.bincount(y, minlength=2) / len(y), base=2)

        IGs = np.zeros(X.shape[1])
        for j in range(X.shape[1]):
            snp = X[:, j]
            H_C_given_snp = 0
            for val in [0, 1]:
                idx = (snp == val)
                if np.any(idx):
                    y_sub = y[idx]
                    probs = np.bincount(y_sub, minlength=2) / len(y_sub)
                    H_C_given_snp += (len(y_sub) / len(y)) * entropy(probs, base=2)
            IGs[j] = H_C - H_C_given_snp

        return dict(sorted(zip(self.__headers["alleles"], IGs), key=lambda x: x[1], reverse=True))


    def testSnpAssociations( self,
                             guids: list,
                             tests: list = None) -> dict:
        """
        Performs statistical tests to evaluate the significance of the association
        of each SNP with a specified group of samples (guids).

        This function is for metrics that typically return a test statistic and a p-value.

        Args:
            guids (list): A list of sample identifiers for the in-group.
            metrics (string): Specify which test to run. Default : 'chi2'.
                Available tests:
                - 'chi2': Chi-squared test of independence (Not Implemented).
                - 'fisher': Fisher's Exact Test (Not Implemented).

        Returns:
            dict: A dictionary of snp : [score, p-value] pairs.
        """
        raise NotImplementedError("This function is not yet implemented.")




    ########## PUBLIC UTILITY FUNCTIONS ##########


    def matchAlleleNames( self,
                          expression: str ) -> list:
        """ Return all allele names that match the given expression with wildcards.
        """
        if not isinstance(expression, str):
            raise ValueError("Expression must be a string.")

        pattern = re.compile(expression.replace('*', '.*'))
        return set([allele for allele in self.__headers["alleles"] if pattern.match(allele)])
    

    def subset( self,
                guid_list : list = [],
                allele_list : list = [] ) -> list[np.array, dict]:
        """ Take a list of GUIDs and subset the allele matrix to include only these GUIDs, allowing for standard operations.
        Returns a numpy matrix/JSON pair for feeding into Ardal.
        """
        ## check input
        if not isinstance(guid_list, list):
            raise ValueError("guid_list must be a list.")
        if not isinstance(allele_list, list):
            raise ValueError("allele_list must be a list.")
        if len(guid_list) == 0 and len(allele_list) == 0:
            raise ValueError("guid_list and allele_list cannot both be empty.")

        ## check GUIDs
        if guid_list:
            final_guids = self._checkGUIDs(guid_list)
            if not final_guids:
                raise ValueError("None of the provided GUIDs were found in the matrix.")
        else:
            final_guids = self.getHeaders()["guids"]

        ## check alleles
        if allele_list:
            final_alleles = self._checkAlleles(allele_list)
            if not final_alleles:
                raise ValueError("None of the provided alleles were found in the matrix.")
        else:
            final_alleles = self.getHeaders()["alleles"]

        ## subset the DataFrame
        ## TODO: this is pretty grim, could be done better in C++
        subset_df = self.toDataFrame().loc[final_guids, final_alleles]

        ## create new headers and matrix for the new Ardal object
        new_headers = {"guids": subset_df.index.tolist(), "alleles": subset_df.columns.tolist()}
        new_matrix = subset_df.values.astype(np.uint8)
        
        ## return the new subset matrix/JSON pair
        return [new_matrix, new_headers]
        
    
    def write( self,
               file_path : str,
               output_prefix : str,
               npz : bool = False ) -> None:
        """ Write the allele matrix to disk.
        Writes as a numpy/JSON pair.
        The npz flag writes the numpy matrix as a compressed npz.
        """
        if not os.path.exists(file_path):
            raise ValueError(f"Path '{file_path}' does not exist.")

        json_out_path = os.path.join(file_path, output_prefix + "_headers.json")
        matrix_out_path = os.path.join(file_path, output_prefix + "_matrix.npy")

        if os.path.exists(json_out_path):
            raise ValueError(f"File '{json_out_path}' already exists.")
        if os.path.exists(matrix_out_path):
            raise ValueError(f"File '{matrix_out_path}' already exists.")

        np.save(matrix_out_path, self.__hybrid_matrix.getMatrix())
        with open(os.path.join(file_path, output_prefix + "_headers.json")) as fout:
            json.dump(self.__headers, fout, indent=4)


    def toDict( self,
                force_bit_backend : bool = False ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        allele_dict = defaultdict(list)
        for guid_idx, guid_name in enumerate(self.__headers["guids"]):
            snp_indices = self.__hybrid_matrix.getSetBitIndices(guid_idx, force_bit_backend=force_bit_backend)
            for snp_idx in snp_indices:
                allele_dict[guid_name].append(self._decodeAllele(snp_idx))
        return dict(allele_dict)
    

    def stats( self,
               print_stats : bool = False ) -> dict:
        """ Return a dictionary containing information about the database and its size in memory.
        """
        n_guids = len(self.__headers["guids"])
        n_alleles = len(self.__headers["alleles"])
        density = self.getDensity()
        roaring = self.roaring
        bit_matrix_size_bytes = self.getBitMatrix().nbytes
        if self.roaring:
            roaring_mat = self.getRoaringMatrix(decode=False)
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
    

    def getDensity( self ) -> float:
        """ Computes the sparsity of the allele matrix.
        """
        return self.__hybrid_matrix.getDensity()
    

    def getBitMatrix( self ) -> np.array:
        """ Return the bit allele matrix.
        """
        return self.__hybrid_matrix.getBitMatrix()
        
    
    def getRoaringMatrix( self,
                          decode : bool=True ) -> dict:
        """ Return the roaring allele matrix.
        """
        roaring_dict = defaultdict(list)
        guids = self.__headers["guids"]

        if self.roaring:
            rormat = self.__hybrid_matrix.getRoaringMatrix()

            if decode:
                for i, mat in enumerate(rormat):
                    allele_ids = [self._decodeAllele(idx) for idx in mat]
                    roaring_dict[guids[i]]=allele_ids
                return roaring_dict
            else:
                return rormat
                
        else:
            raise ValueError("Ardal object was instantialised with 'use_roaring_if_sparse=False'. Cannot retrieve roaring matrix.")
    

    def getHeaders( self ) -> dict:
        """ Return the allele __headers.
        """
        return self.__headers
    

    def snpCount( self ) -> dict:
        """ Return a dictionary of SNP counts for each GUID.
        """
        guid_mass_vec = self.__hybrid_matrix.getRowMasses()
        return {guid : mass for guid, mass in zip(self.__headers["guids"], guid_mass_vec)}
    

    def toDataFrame( self ) -> pd.DataFrame:
        """ Return the allele matrix as a Pandas DataFrame.
        """
        return pd.DataFrame(self.getMatrix(), index=self.__headers["guids"], columns=self.__headers["alleles"])
    

    def flushCache( self ) -> None:
        """ flushes the distance cache.
        """
        # flushCache is not exposed anymore and the cache is not implemented.
        pass
    



    ########## PRIVATE UTILITY FUNCTIONS ##########


    def _checkGUIDs( self,
                     guids : list ) -> list:
        """ Check guids are present within the matrix and construct a list of present guids to proceed with
        """

        present_guids = []
        for id in guids:
            if id in self.__headers["guids"]:
                present_guids.append(id)
            else:
                print(f"{id} not present in allele matrix.")

        return present_guids
    

    def _checkAlleles( self,
                       alleles : list ) -> list:
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