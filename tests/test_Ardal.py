import sys
import pytest
import numpy as np
import pandas as pd
from ardal import Ardal, ArdalParser

from .conftest import *


def test_unique(ardal_object):
    """
    Tests the unique method of the ardal class.
    """

    ## test case 1: SNPs unique to GUID1 and GUID2
    unique_snps1 = ardal_object.unique(["GUID1", "GUID2"])
    assert unique_snps1 == {"SNP1"}

    ## test case 2: SNPs unique to GUID3, GUID4, and GUID5
    unique_snps2 = ardal_object.unique(["GUID3", "GUID4", "GUID5"])
    assert unique_snps2 == {"SNP2"}

    ## test case 3: SNPs unique to GUID6 and GUID7
    unique_snps3 = ardal_object.unique(["GUID6", "GUID7"])
    assert unique_snps3 == {"SNP3"}

    ### test case 4: SNPs unique to GUID8, GUID9, and GUID10
    unique_snps4 = ardal_object.unique(["GUID8", "GUID9", "GUID10"])
    assert unique_snps4 == {"SNP4"}

    ## test case 5: SNPs unique to GUID1, GUID2, GUID3, GUID4, GUID5, GUID6, GUID7
    unique_snps5 = ardal_object.unique(["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7"])
    assert unique_snps5 == {"SNP6"}

    ## test case 6: SNPs unique to GUID6, GUID7, GUID8, GUID9, GUID10
    unique_snps6 = ardal_object.unique(["GUID6", "GUID7", "GUID8", "GUID9", "GUID10"])
    assert unique_snps6 == {"SNP9"}

    ## test case 7: SNPs unique to GUID1, GUID2, GUID3, GUID4, GUID5, GUID6, GUID7, GUID8, GUID9, GUID10
    unique_snps7 = ardal_object.unique(["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9", "GUID10"])
    assert unique_snps7 == {"SNP5", "SNP7"}

    ## test case 8: No unique SNPs
    unique_snps8 = ardal_object.unique(["GUID1", "GUID5"])
    assert unique_snps8 == set()


def test_pairwise(ardal_object):
    """
    Tests the pairwise method of the ardal class.
    """

    ## test case 1: Hamming distance
    dist_matrix = ardal_object.pairwise(metric="hamming")
    assert isinstance(dist_matrix, pd.DataFrame)
    assert dist_matrix.shape == (10, 10)
    assert dist_matrix.to_dict() == {'GUID1': {'GUID1': 0, 'GUID2': 0, 'GUID3': 2, 'GUID4': 2, 'GUID5': 2, 'GUID6': 4, 'GUID7': 4, 'GUID8': 5, 'GUID9': 5, 'GUID10': 6},
                                     'GUID2': {'GUID1': 0, 'GUID2': 0, 'GUID3': 2, 'GUID4': 2, 'GUID5': 2, 'GUID6': 4, 'GUID7': 4, 'GUID8': 5, 'GUID9': 5, 'GUID10': 6},
                                     'GUID3': {'GUID1': 2, 'GUID2': 2, 'GUID3': 0, 'GUID4': 0, 'GUID5': 0, 'GUID6': 4, 'GUID7': 4, 'GUID8': 5, 'GUID9': 5, 'GUID10': 6},
                                     'GUID4': {'GUID1': 2, 'GUID2': 2, 'GUID3': 0, 'GUID4': 0, 'GUID5': 0, 'GUID6': 4, 'GUID7': 4, 'GUID8': 5, 'GUID9': 5, 'GUID10': 6},
                                     'GUID5': {'GUID1': 2, 'GUID2': 2, 'GUID3': 0, 'GUID4': 0, 'GUID5': 0, 'GUID6': 4, 'GUID7': 4, 'GUID8': 5, 'GUID9': 5, 'GUID10': 6},
                                     'GUID6': {'GUID1': 4, 'GUID2': 4, 'GUID3': 4, 'GUID4': 4, 'GUID5': 4, 'GUID6': 0, 'GUID7': 0, 'GUID8': 3, 'GUID9': 3, 'GUID10': 4},
                                     'GUID7': {'GUID1': 4, 'GUID2': 4, 'GUID3': 4, 'GUID4': 4, 'GUID5': 4, 'GUID6': 0, 'GUID7': 0, 'GUID8': 3, 'GUID9': 3, 'GUID10': 4},
                                     'GUID8': {'GUID1': 5, 'GUID2': 5, 'GUID3': 5, 'GUID4': 5, 'GUID5': 5, 'GUID6': 3, 'GUID7': 3, 'GUID8': 0, 'GUID9': 0, 'GUID10': 1},
                                     'GUID9': {'GUID1': 5, 'GUID2': 5, 'GUID3': 5, 'GUID4': 5, 'GUID5': 5, 'GUID6': 3, 'GUID7': 3, 'GUID8': 0, 'GUID9': 0, 'GUID10': 1},
                                     'GUID10': {'GUID1': 6, 'GUID2': 6, 'GUID3': 6, 'GUID4': 6, 'GUID5': 6, 'GUID6': 4, 'GUID7': 4, 'GUID8': 1, 'GUID9': 1, 'GUID10': 0}}

    ## test case 2: Jaccard distance
    dist_matrix = ardal_object.pairwise(metric="jaccard")
    assert isinstance(dist_matrix, pd.DataFrame)
    assert dist_matrix.shape == (10, 10)

    ## test case 3: Invalid metric
    with pytest.raises(ValueError):
        ardal_object.pairwise(metric="invalid")


def test_neighbourhood(ardal_object):
    """
    Tests the neighbourhood method of the ardal class.
    """

    ## test case 1: SIMD
    neighbourhood = ardal_object.neighbourhood("GUID1", n=2, simd=True)
    assert isinstance(neighbourhood, dict)
    assert neighbourhood == {'GUID2': 0, 'GUID3': 2, 'GUID4': 2, 'GUID5': 2}

    ## test case 2: Non-SIMD
    neighbourhood = ardal_object.neighbourhood("GUID1", n=2, simd=False)
    assert isinstance(neighbourhood, dict)
    assert neighbourhood == {'GUID2': 0, 'GUID3': 2, 'GUID4': 2, 'GUID5': 2}

    ## test case 3: Invalid GUID
    with pytest.raises(ValueError):
        ardal_object.neighbourhood("INVALID", n=2)


def test_core(ardal_object):
    """
    Tests the core method of the ardal class.
    """

    ## test case 1: Core alleles
    core_alleles = ardal_object.core(["GUID1", "GUID2"])
    assert isinstance(core_alleles, set)

    ## test case 2: Invalid GUID
    with pytest.raises(ValueError):
        ardal_object.core(["INVALID"])


def test_accessory(ardal_object):
    """
    Tests the accessory method of the ardal class.
    """

    ## test case 1: Accessory alleles
    accessory_alleles = ardal_object.accessory(["GUID1", "GUID2"])
    assert isinstance(accessory_alleles, set)

    ## test case 2: Invalid GUID
    with pytest.raises(ValueError):
        ardal_object.accessory(["INVALID"])


def test_allele(ardal_object):
    """
    Tests the allele method of the ardal class.
    """

    ## test case 1: SNP1
    alleles = ardal_object.allele(["SNP1"])
    assert alleles == {"GUID1", "GUID2"}

    ## test case 2: SNP2
    alleles = ardal_object.allele(["SNP2"])
    assert alleles == {"GUID3", "GUID4", "GUID5"}

    ## test case 3: SNP3
    alleles = ardal_object.allele(["SNP3"])
    assert alleles == {"GUID6", "GUID7"}

    ## test case 4: SNP4
    alleles = ardal_object.allele(["SNP4"])
    assert alleles == {"GUID8", "GUID9", "GUID10"}

    ## test case 5: SNP5
    alleles = ardal_object.allele(["SNP5"])
    assert alleles == {"GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9", "GUID10"}

    ## test case 6: SNP6
    alleles = ardal_object.allele(["SNP6"])
    assert alleles == {"GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7"}

    ## test case 7: SNP7
    alleles = ardal_object.allele(["SNP7"])
    assert alleles == {"GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9", "GUID10"}

    ## test case 8: SNP8
    alleles = ardal_object.allele(["SNP8"])
    assert alleles == {"GUID1", "GUID2", "GUID3", "GUID4", "GUID5"}

    ## test case 9: SNP9
    alleles = ardal_object.allele(["SNP9"])
    assert alleles == {"GUID6", "GUID7", "GUID8", "GUID9", "GUID10"}

    ## test case 10: SNP10
    alleles = ardal_object.allele(["SNP10"])
    assert alleles == {"GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9"}

    ## test case 10: SNP10
    alleles = ardal_object.allele(["SNP5", "SNP10"])
    assert alleles == {'GUID6', 'GUID9', 'GUID2', 'GUID8', 'GUID1', 'GUID3', 'GUID5', 'GUID7', 'GUID4'}

    ## test case 11: Invalid Allele
    with pytest.raises(ValueError):
        ardal_object.allele(["INVALID"])


def test_matchAlleleNames(ardal_object):
    """
    Tests the matchAlleleNames method of the ardal class.
    """

    ## test case 1: Match Allele Names
    matched_alleles = ardal_object.matchAlleleNames("SNP*")
    assert isinstance(matched_alleles, set)

    ## test case 2: Invalid Input
    with pytest.raises(ValueError):
        ardal_object.matchAlleleNames(1)


def test_stats(ardal_object):
    """
    Tests the stats method of the ardal class.
    """

    ## test case 1: Stats
    stats = ardal_object.stats()
    assert isinstance(stats, dict)


def test_getMatrix(ardal_object):
    """
    Tests the getMatrix method of the ardal class.
    """

    ## test case 1: Get Matrix
    matrix = ardal_object.getMatrix()
    assert isinstance(matrix, np.ndarray)


def test_getHeaders(ardal_object):
    """
    Tests the getHeaders method of the ardal class.
    """

    ## test case 1: Get Headers
    headers = ardal_object.getHeaders()
    assert isinstance(headers, dict)


def test_snpCount(ardal_object):
    """
    Tests the snpCount method of the ardal class.
    """

    ## test case 1: Snp Count
    snp_count = ardal_object.snpCount()
    assert snp_count == {'GUID1': 6,
                         'GUID2': 6,
                         'GUID3': 6,
                         'GUID4': 6,
                         'GUID5': 6,
                         'GUID6': 6,
                         'GUID7': 6,
                         'GUID8': 5,
                         'GUID9': 5,
                         'GUID10': 4}


def test_toDataFrame(ardal_object):
    """
    Tests the toDataFrame method of the ardal class.
    """

    ## test case 1: To Data Frame
    df = ardal_object.toDataFrame()
    assert isinstance(df, pd.DataFrame)


def test_flushCache(ardal_object):
    """
    Tests the flushCache method of the ardal class.
    """

    ## test case 1: Flush Cache
    ardal_object.flushCache()
