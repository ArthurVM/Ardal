import sys
import pytest
import numpy as np
from ardal import Ardal, ArdalParser


def test_access(ardal_object):
    """
    Tests the access method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Access
    coords = np.array([[0, 0], [1, 1]], dtype=np.uint64)
    result = matrix.access(coords)
    assert isinstance(result, np.ndarray)

    ## test case 2: Invalid Coordinates
    coords = np.array([[100, 100]], dtype=np.uint64)
    with pytest.raises(RuntimeError):
        matrix.access(coords)


def test_accessGUID(ardal_object):
    """
    Tests the accessGUID method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Access GUID
    result = matrix.accessGUID(0)
    assert isinstance(result, set)


def test_getMatrix(ardal_object):
    """
    Tests the getMatrix method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Get Matrix
    result = matrix.getMatrix()
    assert isinstance(result, np.ndarray)


def test_getMass(ardal_object):
    """
    Tests the getMass method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Get Mass
    result = matrix.getMass()
    assert isinstance(result, list)


def test_hamming(ardal_object):
    """
    Tests the hamming method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Hamming
    result = matrix.hamming()
    assert isinstance(result, np.ndarray)


def test_jaccard(ardal_object):
    """
    Tests the jaccard method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Jaccard
    result = matrix.jaccard()
    assert isinstance(result, np.ndarray)


def test_neighbourhood(ardal_object):
    """
    Tests the neighbourhood method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Neighbourhood
    result = matrix.neighbourhood(0, 2)
    assert isinstance(result, np.ndarray)

    ## test case 2: Invalid Row Coord
    with pytest.raises(RuntimeError):
        matrix.neighbourhood(100, 2)

    ## test case 3: Invalid Epsilon
    with pytest.raises(RuntimeError):
        matrix.neighbourhood(0, -1)


def test_neighbourhoodSIMD(ardal_object):
    """
    Tests the neighbourhoodSIMD method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Neighbourhood SIMD
    result = matrix.neighbourhoodSIMD(0, 2)
    assert isinstance(result, list)

    ## test case 2: Invalid Row Coord
    with pytest.raises(RuntimeError):
        matrix.neighbourhoodSIMD(100, 2)

    ## test case 3: Invalid Epsilon
    with pytest.raises(RuntimeError):
        matrix.neighbourhoodSIMD(0, -1)


def test_gatherSNPs(ardal_object):
    """
    Tests the gatherSNPs method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Gather SNPs
    guid_indices = np.array([0, 1], dtype=np.int32)
    result = matrix.gatherSNPs(guid_indices)
    assert isinstance(result, list)


def test_flushCache(ardal_object):
    """
    Tests the flushCache method of the AlleleMatrix class.
    """
    obj = ardal_object
    matrix = obj._Ardal__allele_matrix

    ## test case 1: Flush Cache
    matrix.flushCache()
