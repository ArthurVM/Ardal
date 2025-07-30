""" ArdalStats.py
This module provides functionality for calculating statistics on allele matrices in the Ardal framework.
"""
import numpy as np
import os
from collections import defaultdict

from .utilities import *


# core/ArdalStats.py
class ArdalStats:
    """ Class for calculating statistics on allele matrices in the Ardal framework.
    """

    def __init__(self, headers, hybrid_matrix, roaring_enabled):

        self._headers = headers
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled


    def af( self,
            guids : list=[] ) -> dict:
        """ Calculates the allele frequency for each allele in the matrix across the group of guids.
        if an empty list is provided then all guids will be used.
        """
        if guids == []:
            guids = self._headers["guids"]
        else:
            checkGUIDs(guids, self._headers)
        
        guid_coords = [encodeGuid(guid, self._headers) for guid in guids]
        allele_freq = self._matrix.colFrequency(guid_coords)

        return dict(sorted(zip(self._headers["alleles"], allele_freq), key=lambda x: x[1], reverse=True))
        

    def alleleEntropy( self ) -> dict:
        """ Computes the Shannon entropy for each allele in the matrix.
        Entropy H(X) = -p * log2(p) - (1-p) * log2(1-p), where p is the allele frequency.
        """
        allele_entropy = self._matrix.columnEntropy()
        entropies_dict = dict(sorted(zip(self._headers["alleles"], allele_entropy), key=lambda x: x[1], reverse=True))
        return entropies_dict

    
    def alleleCooc( self,
                    allele_indices : list = [],
                    threshold : float = 0.95,
                    threads : int = 1 ) -> dict:
        """ Computes the co-occurrence of pairs of alleles.
        """
        if not isinstance(allele_indices, list):
            raise ValueError("allele_indices must be a list.")
        if not isinstance(threshold, float):
            raise ValueError("threshold must be a float.")
        if not isinstance(threads, int):
            raise ValueError("threads must be an integer.")
        if threshold < 0 or threshold > 1:
            raise ValueError("threshold must be between 0 and 1.")
        if threads < 1:
            raise ValueError("threads must be at least 1.")

        if allele_indices != []:
            checkAlleles(allele_indices, self._headers)

            ## encode the alleles
            allele_indices = [encodeAllele(allele, self._headers) for allele in allele_indices]

            cooc_dict = self._matrix.bitCooccurrence_subset(col_indices=allele_indices,
                                                            threshold=threshold,
                                                            threads=threads)
            
        else:
            cooc_dict = self._matrix.bitCooccurrence_all(threshold=threshold,
                                                         threads=threads)

        decoded_dict = defaultdict(list)
        for ref, cooc_vec in cooc_dict.items():
            decoded_dict[decodeAllele(ref, self._headers)] = [decodeAllele(allele, self._headers) for allele in cooc_vec]

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
        
        checkGUIDs(guids, self._headers)

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
        D_{kl}(P||Q) = sum_{x in X}(P(x) * log2(P(x)/Q(x)))
        """
        guid_coords = [encodeGuid(guid, self._headers) for guid in guids]
        kl_divergence = self._matrix.klDivergence(guid_coords)
        kl_dict = dict(sorted(zip(self._headers["alleles"], kl_divergence), key=lambda x: x[1], reverse=True))
        return kl_dict

    
    def _jsDivergence( self,
                       guids : list ) -> dict:
        """
        Computes Jensen-Shannon divergence for each SNP between target_guids and others.
        
        Args:
            guids: list of sample identifiers to define the target group
        """
        guid_coords = [encodeGuid(guid, self._headers) for guid in guids]
        js_divergence = self._matrix.jsDivergence(guid_coords)
        js_dict = dict(sorted(zip(self._headers["alleles"], js_divergence), key=lambda x: x[1], reverse=True))
        return js_dict


    def _informationGain( self,
                          guids : list ) -> dict:
        """
        Compute information gain for each SNP column in binary matrix X,
        with respect to whether a sample is in target_guids.

        Returns:
            np.ndarray of shape (n_snps,) - information gain per SNP
        """
        guid_coords = [encodeGuid(guid, self._headers) for guid in guids]
        ig = self._matrix.informationGain(guid_coords)
        ig_dict = dict(sorted(zip(self._headers["alleles"], ig), key=lambda x: x[1], reverse=True))
        return ig_dict


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
    