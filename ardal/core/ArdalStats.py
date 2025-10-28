""" ArdalStats.py
This module provides functionality for calculating statistics on allele matrices in the Ardal framework.
"""
import numpy as np
from collections import defaultdict
from typing import Union, List, Dict

from ..utils.exceptions import ParameterError, EmptySelectionError
from ..utils.decorators import check_generic_threshold, check_thread_count, check_guids_list
from ..utils.logger import get_logger
from ..utils.misc import require_package

log = get_logger()


## core/ArdalStats.py
class ArdalStats:
    """ Class for calculating statistics on allele matrices in the Ardal framework.
    """

    def __init__( self,
                  headerUtils,
                  hybrid_matrix,
                  roaring_enabled : bool ):

        self._headerUtils = headerUtils
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled


    def allele_frequency( self,
                          guids : List = []
                          ) -> Dict:
        """ Calculates the allele frequency for each allele in the matrix across the group of guids.
        if an empty list is provided then all guids will be used.
        """
        if guids == []:
            guids = self._headerUtils.headers["guids"]
        else:
            self._headerUtils.check_guids(guids)
        
        guid_coords = [self._headerUtils.encode_guid(guid) for guid in guids]
        allele_freq = self._matrix.colFrequency(guid_coords)

        return dict(sorted(zip(self._headerUtils.headers["alleles"], allele_freq), key=lambda x: x[1], reverse=True))
        

    def allele_entropy( self ) -> dict:
        """ Computes the Shannon entropy for each allele in the matrix.
        Entropy H(X) = -p * log2(p) - (1-p) * log2(1-p), where p is the allele frequency.
        """
        log.info(f"Computing Shannon entropy for all alleles.")
        allele_entropy = self._matrix.columnEntropy()
        entropies_dict = dict(sorted(zip(self._headerUtils.headers["alleles"], allele_entropy), key=lambda x: x[1], reverse=True))
        return entropies_dict

    
    @check_thread_count
    @check_generic_threshold
    def allele_cooc( self,
                     alleles : list = [],
                     threshold : float = 0.95,
                     threads : int = 1
                     ) -> dict:
        """ Computes the co-occurrence of pairs of alleles.
        """
        log.info(f"Computing allele co-occurrence using threshold {threshold}")
        
        if alleles != []:
            ## input validation
            self._headerUtils.check_alleles(alleles)

            ## encode the alleles
            allele_indices = [self._headerUtils.encode_allele(allele) for allele in alleles]

            cooc_dict = self._matrix.bitCooccurrence_subset(col_indices=allele_indices,
                                                            threshold=threshold,
                                                            threads=threads)
            
        else:
            cooc_dict = self._matrix.bitCooccurrence_all(threshold=threshold,
                                                         threads=threads)

        decoded_dict = defaultdict(list)
        for ref, cooc_vec in cooc_dict.items():
            decoded_dict[self._headerUtils.decode_allele(ref)] = [self._headerUtils.decode_allele(allele) for allele in cooc_vec]

        return decoded_dict


    @check_guids_list
    def allele_inform( self,
                       guids: List,
                       method: str = "kullbackleibler"
                       ) -> Union[Dict, None]:
        """
        Calculates scores that measure the association of each allele
        with a specified group of samples (guids).

        This function serves as a unified interface for metrics that return a single
        float value per allele, representing the strength or nature of the association.
        It compares an "in-group" (the provided guids) against an "out-group"
        (all other samples).

        Args:
            guids (list): A list of sample identifiers for the in-group.
            method (string): Specify which metric to use. Default : 'kullbackleibler'.
                Available methods:
                - 'kullbackleibler': Kullback-Leibler divergence.
                - 'jensenshannon': Jensen-Shannon divergence.
                - 'informationgain': Information Gain.

        Returns:
            dict: A dictionary of allele_id : score pairs.
        """
        ## input validation
        available_methods = {'kullbackleibler', 'jensenshannon', 'informationgain'}
        lower_method = method.lower()
        if lower_method not in available_methods:
            raise ParameterError(f"Method '{lower_method}' is not supported. Available method: {list(available_methods)}")
        
        self._headerUtils.check_guids(guids)
        
        log.info(f"Computing allele information.")

        if lower_method == 'kullbackleibler':
            log.debug(f"""method = Kullback-Leibler divergence
                                 n_guids = {len(guids)}""")
            return self._kl_divergence(guids)
        if lower_method == 'jensenshannon':
            log.debug(f"""method = Jensen-Shannon divergence
                                 n_guids = {len(guids)}""")
            return self._js_divergence(guids)
        if lower_method == 'informationgain':
            log.debug(f"""method = Information Gain
                                 n_guids = {len(guids)}""")
            return self._information_gain(guids)
    
    
    def _kl_divergence( self,
                        guids : list
                        ) -> dict:
        """ Computes the Kullbeck Liebler divergence between the in group (input guids) and out group (all others)
        allele frequency distributions for each allele in the matrix.
        D_{kl}(P||Q) = sum_{x in X}(P(x) * log2(P(x)/Q(x)))
        """
        guid_coords = [self._headerUtils.encode_guid(guid) for guid in guids]
        kl_divergence = self._matrix.klDivergence(guid_coords)
        kl_dict = dict(sorted(zip(self._headerUtils.headers["alleles"], kl_divergence), key=lambda x: x[1], reverse=True))
        return kl_dict

    
    def _js_divergence( self,
                        guids : List
                        ) -> Dict:
        """
        Computes Jensen-Shannon divergence for each allele between target_guids and others.
        
        Args:
            guids: list of sample identifiers to define the target group
        """
        guid_coords = [self._headerUtils.encode_guid(guid) for guid in guids]
        js_divergence = self._matrix.jsDivergence(guid_coords)
        js_dict = dict(sorted(zip(self._headerUtils.headers["alleles"], js_divergence), key=lambda x: x[1], reverse=True))
        return js_dict


    def _information_gain( self,
                           guids : List
                           ) -> Dict:
        """
        Compute information gain for each allele column in binary matrix X,
        with respect to whether a sample is in target_guids.

        Returns:
            np.ndarray of shape (n_alleles,) - information gain per allele
        """
        guid_coords = [self._headerUtils.encode_guid(guid) for guid in guids]
        ig = self._matrix.informationGain(guid_coords)
        ig_dict = dict(sorted(zip(self._headerUtils.headers["alleles"], ig), key=lambda x: x[1], reverse=True))
        return ig_dict


    def test_associations( self,
                           guids: List,
                           tests: List = None
                           ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """
        Performs statistical tests to evaluate the association of each allele
        with a specified group of samples (guids). Each test returns both a
        test statistic and a p-value for a 2x2 contingency table comparing the
        in-group against the out-group.

        Args:
            guids (list): Identifiers for the in-group samples.
            tests (list): Iterable of test names to run. Supported values:
                - 'chi2': Pearson's Chi-squared test with Yates correction disabled.
                - 'fisher': Fisher's exact test (two-sided).
                If None, defaults to ['chi2'].

        Returns:
            dict: Mapping allele_id -> {test_name: {"statistic": float, "pvalue": float}}.

        Raises:
            ParameterError: If requested tests are unsupported or the in-group
                covers all samples (i.e. no out-group remains).
        """
        if not guids:
            raise ParameterError("At least one GUID must be provided to perform association tests.")

        available_tests = {"chi2", "fisher"}
        if tests is None:
            tests_to_run = ["chi2"]
        else:
            tests_to_run = []
            for test in tests:
                if not isinstance(test, str):
                    raise ParameterError(f"Test name must be a string, received {type(test)}.")
                lower = test.lower()
                if lower not in available_tests:
                    raise ParameterError(f"Unsupported test '{test}'. Available tests: {sorted(available_tests)}")
                if lower not in tests_to_run:
                    tests_to_run.append(lower)

        total_samples = self._headerUtils.n_guids
        if total_samples == 0:
            raise EmptySelectionError("Association testing requires a non-empty matrix.")

        self._headerUtils.check_guids(guids)
        guid_coords = [self._headerUtils.encode_guid(guid) for guid in guids]
        n_in = len(guid_coords)
        n_out = total_samples - n_in

        if n_out <= 0:
            raise ParameterError("Association tests require at least one out-group sample.")

        allele_ids = self._headerUtils.headers["alleles"]

        freq_in = np.asarray(self._matrix.colFrequency(guid_coords), dtype=np.float64)
        total_masses = np.asarray(self._matrix.getColumnMasses(), dtype=np.int64)

        count_in_pos = np.rint(freq_in * n_in).astype(np.int64, copy=False)
        np.clip(count_in_pos, 0, n_in, out=count_in_pos)

        count_out_pos = total_masses - count_in_pos
        np.clip(count_out_pos, 0, n_out, out=count_out_pos)

        count_in_neg = n_in - count_in_pos
        count_out_neg = n_out - count_out_pos

        results: Dict[str, Dict[str, Dict[str, float]]] = {allele: {} for allele in allele_ids}

        stats_module = require_package("scipy.stats", "scipy.stats")
        chi2_fn = getattr(stats_module, "chi2_contingency", None)
        fisher_fn = getattr(stats_module, "fisher_exact", None)

        for idx, allele in enumerate(allele_ids):
            contingency = np.array([
                [int(count_in_pos[idx]), int(count_in_neg[idx])],
                [int(count_out_pos[idx]), int(count_out_neg[idx])]
            ], dtype=np.int64)

            allele_results: Dict[str, Dict[str, float]] = {}

            if "chi2" in tests_to_run and chi2_fn is not None:
                try:
                    chi2_stat, chi2_pvalue, _, _ = chi2_fn(contingency, correction=False)
                except ValueError:
                    chi2_stat, chi2_pvalue = np.nan, np.nan
                allele_results["chi2"] = {
                    "statistic": float(chi2_stat),
                    "pvalue": float(chi2_pvalue)
                }

            if "fisher" in tests_to_run and fisher_fn is not None:
                try:
                    odds_ratio, fisher_pvalue = fisher_fn(contingency, alternative="two-sided")
                except ValueError:
                    odds_ratio, fisher_pvalue = np.nan, np.nan
                allele_results["fisher"] = {
                    "statistic": float(odds_ratio),
                    "pvalue": float(fisher_pvalue)
                }

            results[allele] = allele_results

        return results
    
